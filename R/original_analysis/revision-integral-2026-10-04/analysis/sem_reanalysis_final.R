# ==============================================================================
# REANALISIS CANONICO CFA / SEM - REVISION 2026-10-04
# ==============================================================================
# Objetivos
#   1. Auditar la base sin modificarla.
#   2. Separar el modelo de medicion del modelo estructural.
#   3. Evaluar por pais el CFA de SC, SBL y PSB y compararlo con modelos rivales.
#   4. Usar MLR como analisis principal y WLSMV como sensibilidad ordinal.
#   5. Evaluar fiabilidad, AVE, HTMT, correlaciones latentes e invarianza.
#   6. Estimar un modelo asociativo pooled con valores -> actitudes -> prosocialidad,
#      reportando por separado admisibilidad, ajuste y limites sustantivos.
#   7. Registrar cada corrida, hashes, manifiesto y diagnosticos bootstrap validos.
#
# Decisiones importantes
#   - No se fijan varianzas residuales en 0.001.
#   - Benevolencia, Universalismo, Poder y Logro se usan como variables observadas
#     en el modelo estructural principal. No se fuerzan factores de dos indicadores.
#   - Los compuestos de valores se estiman solo como sensibilidad y se etiquetan
#     expresamente como exploratorios.
#   - Las reespecificaciones no se realizan automaticamente a partir de indices de
#     modificacion. El script exporta diagnosticos para una decision teorica.
#   - La base .sav se lee, pero nunca se sobrescribe.
#
# Ejecucion sugerida desde la raiz del proyecto:
#   Rscript analysis/sem_reanalysis_final.R
#
# Variables de entorno opcionales:
#   SEM_SAV_PATH            Ruta explicita al archivo .sav.
#   SEM_BOOT_HTMT           Replicas bootstrap para IC de HTMT/alfa (1000).
#   SEM_BOOT_INDIRECT       Replicas bootstrap ML para indirectos pooled (1000).
#   SEM_RUN_ORDINAL         1 = ejecutar WLSMV/policoricas; 0 = omitir esa
#                           sensibilidad en una corrida MLR rapida.
#   SEM_RUN_RIVALS          1 = estimar tres modelos de medicion rivales; 0 =
#                           omitirlos en una prueba rapida del flujo.
#   SEM_OVERRIDE_GATE       1 = ejecutar modelos posteriores solo como diagnostico
#                           aunque una compuerta de ajuste no se cumpla.
# ==============================================================================

options(stringsAsFactors = FALSE, scipen = 999)

required_packages <- c(
  "haven", "lavaan", "psych", "dplyr", "tidyr", "purrr", "readr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Faltan paquetes requeridos: ", paste(missing_packages, collapse = ", "),
    ". Instálelos antes de ejecutar el script.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(haven)
  library(lavaan)
  library(psych)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
})

# ------------------------------------------------------------------------------
# 0. Rutas, configuracion y utilidades generales
# ------------------------------------------------------------------------------

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1L) {
    encoded_path <- sub("^--file=", "", file_arg)
    decoded_path <- gsub("~+~", " ", encoded_path, fixed = TRUE)
    return(normalizePath(decoded_path, mustWork = TRUE))
  }
  normalizePath(file.path(getwd(), "analysis", "sem_reanalysis_final.R"),
                mustWork = FALSE)
}

script_path <- get_script_path()
script_dir <- dirname(script_path)
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

# Una corrida nunca reutiliza archivos de otra corrida.
run_started_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
run_id <- Sys.getenv("SEM_RUN_ID", "")
if (!nzchar(run_id)) run_id <- paste0(format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "-", Sys.getpid())
if (!grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", run_id)) stop("SEM_RUN_ID contiene caracteres no permitidos.", call. = FALSE)
output_root <- Sys.getenv("SEM_OUTPUT_ROOT", file.path(script_dir, "outputs"))
output_root <- normalizePath(output_root, mustWork = FALSE)
results_dir <- file.path(output_root, run_id)
if (dir.exists(results_dir)) stop("La corrida ya existe; use un SEM_RUN_ID nuevo: ", results_dir, call. = FALSE)
audit_dir <- file.path(results_dir, "01_audit")
cfa_dir <- file.path(results_dir, "02_cfa")
invariance_dir <- file.path(results_dir, "03_invariance")
structural_dir <- file.path(results_dir, "04_structural")

invisible(lapply(
  c(results_dir, audit_dir, cfa_dir, invariance_dir, structural_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

boot_htmt <- suppressWarnings(as.integer(Sys.getenv("SEM_BOOT_HTMT", "1000")))
boot_indirect <- suppressWarnings(as.integer(Sys.getenv("SEM_BOOT_INDIRECT", "1000")))
if (is.na(boot_htmt) || boot_htmt < 0L) boot_htmt <- 1000L
if (is.na(boot_indirect) || boot_indirect < 0L) boot_indirect <- 1000L
run_ordinal <- !identical(Sys.getenv("SEM_RUN_ORDINAL", "1"), "0")
run_rivals <- !identical(Sys.getenv("SEM_RUN_RIVALS", "1"), "0")
override_gate <- identical(Sys.getenv("SEM_OVERRIDE_GATE", "0"), "1")
ncpus <- suppressWarnings(as.integer(Sys.getenv("SEM_NCPUS", "2")))
if (is.na(ncpus) || ncpus < 1L) ncpus <- 1L
set.seed(20261004)

ascii_text <- function(x) {
  out <- iconv(as.character(x), from = "", to = "ASCII//TRANSLIT")
  out
}

normalize_token <- function(x) {
  x <- tolower(trimws(ascii_text(x)))
  gsub("[^a-z0-9]+", "", x)
}

as_scalar <- function(x, default = NA_real_) {
  if (length(x) == 0L || all(is.na(x))) return(default)
  unname(as.numeric(x[[1L]]))
}

safe_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  x <- as.data.frame(x)
  if (ncol(x) == 0L) x <- data.frame(note = character(), stringsAsFactors = FALSE)
  readr::write_csv(x, path, na = "")
}

safe_text <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(enc2utf8(as.character(x)), con = path, useBytes = TRUE)
}

write_manifest <- function() {
  files <- list.files(results_dir, recursive = TRUE, full.names = TRUE)
  files <- files[!dir.exists(files) & basename(files) != "manifest.csv"]
  manifest <- tibble(path = substring(files, nchar(results_dir) + 2L),
                     bytes = unname(file.info(files)$size),
                     md5 = unname(tools::md5sum(files)))
  safe_csv(manifest, file.path(results_dir, "manifest.csv"))
}
write_run_status <- function(state, detail = "") {
  safe_csv(tibble(run_id = run_id, state = state, started_utc = run_started_utc,
                  updated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), detail = detail),
           file.path(results_dir, "run_status.csv"))
}
invisible(file.copy(script_path, file.path(results_dir, "source_script.R"), overwrite = FALSE))
invisible(file.copy(script_path, file.path(results_dir, "executed_script_snapshot.R"), overwrite = FALSE))
write_run_status("STARTED")
options(error = function() {
  detail <- geterrmessage()
  try(write_run_status("FAILED", detail), silent = TRUE)
  try(write_manifest(), silent = TRUE)
  quit(save = "no", status = 1L, runLast = FALSE)
})
message("RUN_ID: ", run_id, " | OUTPUT: ", results_dir)


matrix_to_long <- function(mat, value_name = "value", diagonal = TRUE) {
  if (is.null(mat) || length(mat) == 0L) return(tibble())
  mat <- as.matrix(mat)
  out <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(out) <- c("row", "column", value_name)
  if (!diagonal) out <- out[out$row != out$column, , drop = FALSE]
  tibble::as_tibble(out)
}

capture_conditions <- function(expr) {
  warnings_seen <- character()
  value <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        warnings_seen <<- c(warnings_seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) structure(list(message = conditionMessage(e)), class = "safe_error")
  )
  list(
    value = value,
    error = if (inherits(value, "safe_error")) value$message else NA_character_,
    warnings = unique(warnings_seen)
  )
}

is_lavaan_fit <- function(x) inherits(x, "lavaan")

fit_measure_or_na <- function(fm, name) {
  if (name %in% names(fm)) as_scalar(fm[[name]]) else NA_real_
}

fit_measure_cache <- new.env(parent = emptyenv())
extract_fit_row <- function(fit, model, estimator, country = "Pooled") {
  if (!is_lavaan_fit(fit)) {
    return(tibble(
      country = country, model = model, estimator = estimator,
      n = NA_real_, npar = NA_real_, chisq = NA_real_, df = NA_real_,
      pvalue = NA_real_, chisq_scaled = NA_real_, pvalue_scaled = NA_real_,
      cfi = NA_real_, tli = NA_real_, rmsea = NA_real_, rmsea_lo90 = NA_real_,
      rmsea_hi90 = NA_real_, cfi_robust = NA_real_, tli_robust = NA_real_,
      rmsea_robust = NA_real_, rmsea_robust_lo90 = NA_real_,
      rmsea_robust_hi90 = NA_real_, srmr = NA_real_, aic = NA_real_,
      bic = NA_real_
    ))
  }
  cache_key <- paste(country, model, estimator, sep = "|")
  if (exists(cache_key, fit_measure_cache, inherits = FALSE)) {
    fm <- get(cache_key, fit_measure_cache, inherits = FALSE)
  } else {
    # Evita calcular indices que no se exportan y repetir eigendescomposiciones.
    wanted <- c("ntotal", "npar", "chisq", "df", "pvalue", "chisq.scaled", "pvalue.scaled",
      "cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "cfi.robust", "tli.robust",
      "rmsea.robust", "rmsea.ci.lower.robust", "rmsea.ci.upper.robust", "srmr", "aic", "bic")
    fm <- lavaan::fitMeasures(fit, fit.measures = wanted)
    assign(cache_key, fm, fit_measure_cache)
  }
  tibble(
    country = country,
    model = model,
    estimator = estimator,
    n = fit_measure_or_na(fm, "ntotal"),
    npar = fit_measure_or_na(fm, "npar"),
    chisq = fit_measure_or_na(fm, "chisq"),
    df = fit_measure_or_na(fm, "df"),
    pvalue = fit_measure_or_na(fm, "pvalue"),
    chisq_scaled = fit_measure_or_na(fm, "chisq.scaled"),
    pvalue_scaled = fit_measure_or_na(fm, "pvalue.scaled"),
    cfi = fit_measure_or_na(fm, "cfi"),
    tli = fit_measure_or_na(fm, "tli"),
    rmsea = fit_measure_or_na(fm, "rmsea"),
    rmsea_lo90 = fit_measure_or_na(fm, "rmsea.ci.lower"),
    rmsea_hi90 = fit_measure_or_na(fm, "rmsea.ci.upper"),
    cfi_robust = fit_measure_or_na(fm, "cfi.robust"),
    tli_robust = fit_measure_or_na(fm, "tli.robust"),
    rmsea_robust = fit_measure_or_na(fm, "rmsea.robust"),
    rmsea_robust_lo90 = fit_measure_or_na(fm, "rmsea.ci.lower.robust"),
    rmsea_robust_hi90 = fit_measure_or_na(fm, "rmsea.ci.upper.robust"),
    srmr = fit_measure_or_na(fm, "srmr"),
    aic = fit_measure_or_na(fm, "aic"),
    bic = fit_measure_or_na(fm, "bic")
  )
}

preferred_fit_value <- function(row, robust_name, conventional_name) {
  z <- row[[robust_name]][[1L]]
  if (is.na(z)) row[[conventional_name]][[1L]] else z
}

extract_admissibility <- function(fit, country, model, estimator) {
  if (!is_lavaan_fit(fit)) {
    return(tibble(
      country = country, model = model, estimator = estimator,
      converged = FALSE, post_check = FALSE, negative_residuals = NA_integer_,
      max_abs_latent_correlation = NA_real_, min_eigen_cov_lv = NA_real_,
      admissible = FALSE, fit_acceptable = FALSE
    ))
  }

  converged <- isTRUE(tryCatch(lavInspect(fit, "converged"), error = function(e) FALSE))
  post_check <- isTRUE(tryCatch(lavInspect(fit, "post.check"), error = function(e) FALSE))
  pe <- parameterEstimates(fit)
  negative_residuals <- sum(
    pe$op == "~~" & pe$lhs == pe$rhs & pe$est < -1e-6,
    na.rm = TRUE
  )

  cor_lv <- tryCatch(lavInspect(fit, "cor.lv"), error = function(e) NULL)
  max_abs_cor <- NA_real_
  cor_mats <- if (is.matrix(cor_lv)) list(cor_lv) else cor_lv
  if (is.list(cor_mats) && length(cor_mats)) {
    cor_values <- unlist(lapply(cor_mats, function(m) {
      if (!is.matrix(m) || nrow(m) <= 1L) return(numeric())
      abs(m[row(m) != col(m)])
    }))
    if (length(cor_values)) max_abs_cor <- max(cor_values, na.rm = TRUE)
  }

  cov_lv <- tryCatch(lavInspect(fit, "cov.lv"), error = function(e) NULL)
  min_eigen <- NA_real_
  cov_mats <- if (is.matrix(cov_lv)) list(cov_lv) else cov_lv
  if (is.list(cov_mats) && length(cov_mats)) {
    eig_values <- unlist(lapply(cov_mats, function(m) {
      if (!is.matrix(m)) return(NA_real_)
      tryCatch(
        min(eigen(m, symmetric = TRUE, only.values = TRUE)$values),
        error = function(e) NA_real_
      )
    }))
    if (any(is.finite(eig_values))) min_eigen <- min(eig_values, na.rm = TRUE)
  }

  fit_row <- extract_fit_row(fit, model, estimator, country)
  cfi_eval <- preferred_fit_value(fit_row, "cfi_robust", "cfi")
  rmsea_eval <- preferred_fit_value(fit_row, "rmsea_robust", "rmsea")
  srmr_eval <- fit_row$srmr[[1L]]

  admissible <- converged && post_check && negative_residuals == 0L &&
    (is.na(max_abs_cor) || max_abs_cor < 1) &&
    (is.na(min_eigen) || min_eigen > 1e-8)

  # Son heuristicas transparentes, no leyes universales. Se exportan por separado
  # de la admisibilidad matematica.
  fit_acceptable <- admissible && !is.na(cfi_eval) && !is.na(rmsea_eval) &&
    !is.na(srmr_eval) && cfi_eval >= .90 && rmsea_eval <= .08 && srmr_eval <= .08

  tibble(
    country = country,
    model = model,
    estimator = estimator,
    converged = converged,
    post_check = post_check,
    negative_residuals = negative_residuals,
    max_abs_latent_correlation = max_abs_cor,
    min_eigen_cov_lv = min_eigen,
    admissible = admissible,
    fit_acceptable = fit_acceptable
  )
}

extract_parameters <- function(fit, country, model) {
  if (!is_lavaan_fit(fit)) return(tibble())
  unstd_raw <- parameterEstimates(
    fit, standardized = FALSE, ci = TRUE, level = .95
  ) %>% as_tibble()
  if (!"group" %in% names(unstd_raw)) unstd_raw$group <- 1L
  unstd <- unstd_raw %>%
    transmute(
      lhs, op, rhs, group,
      est_unstd = est, se_unstd = se, z_unstd = z, p_unstd = pvalue,
      ci_unstd_low = ci.lower, ci_unstd_high = ci.upper
    )

  std_raw <- standardizedSolution(
    fit, type = "std.all", se = TRUE, zstat = TRUE, pvalue = TRUE,
    ci = TRUE, level = .95
  ) %>% as_tibble()
  if (!"group" %in% names(std_raw)) std_raw$group <- 1L
  std <- std_raw %>%
    transmute(
      lhs, op, rhs, group,
      est_std = est.std, se_std = se, z_std = z, p_std = pvalue,
      ci_std_low = ci.lower, ci_std_high = ci.upper
    )

  full_join(unstd, std, by = c("lhs", "op", "rhs", "group")) %>%
    mutate(country = country, model = model, .before = 1)
}

extract_r2 <- function(fit, country, model) {
  if (!is_lavaan_fit(fit)) return(tibble())
  r2 <- tryCatch(lavInspect(fit, "r2"), error = function(e) NULL)
  if (is.null(r2)) return(tibble())
  if (is.list(r2) && !is.atomic(r2)) {
    return(map_dfr(seq_along(r2), function(g) {
      tibble(
        country = names(r2)[g] %||% paste0("group_", g),
        model = model,
        outcome = names(r2[[g]]),
        r2 = as.numeric(r2[[g]])
      )
    }))
  }
  tibble(country = country, model = model, outcome = names(r2), r2 = as.numeric(r2))
}

# Operador auxiliar local para nombres posiblemente nulos.
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) return(y)
  if (length(x) == 1L && is.na(x)) return(y)
  x
}

extract_latent_correlations <- function(fit, country, model) {
  if (!is_lavaan_fit(fit)) return(tibble())
  mats <- tryCatch(lavInspect(fit, "cor.lv"), error = function(e) NULL)
  if (is.null(mats)) return(tibble())
  if (is.matrix(mats)) mats <- list(mats)
  group_names <- tryCatch(lavInspect(fit, "group.label"), error = function(e) NULL)
  if (is.null(group_names) || length(group_names) != length(mats)) {
    group_names <- if (length(mats) == 1L) country else paste0("group_", seq_along(mats))
  }
  map_dfr(seq_along(mats), function(i) {
    matrix_to_long(mats[[i]], "latent_correlation", diagonal = FALSE) %>%
      filter(row < column) %>%
      mutate(country = group_names[i], model = model, .before = 1)
  })
}

extract_top_residual_correlations <- function(fit, country, model, n_top = 50L) {
  if (!is_lavaan_fit(fit)) return(tibble())
  rr <- tryCatch(lavResiduals(fit, type = "cor"), error = function(e) NULL)
  if (is.null(rr)) return(tibble())
  mat <- rr$cov
  if (is.list(mat) && !is.matrix(mat)) mat <- mat[[1L]]
  if (!is.matrix(mat)) return(tibble())
  matrix_to_long(mat, "residual_correlation", diagonal = FALSE) %>%
    filter(row < column) %>%
    mutate(abs_residual = abs(residual_correlation)) %>%
    arrange(desc(abs_residual)) %>%
    slice_head(n = n_top) %>%
    mutate(country = country, model = model, .before = 1)
}

extract_modification_indices <- function(fit, country, model, n_top = 100L) {
  if (!is_lavaan_fit(fit)) return(tibble())
  out <- tryCatch(
    modificationIndices(fit, standardized = TRUE, sort. = TRUE),
    error = function(e) NULL
  )
  if (is.null(out)) return(tibble())
  as_tibble(out) %>%
    filter(mi >= 3.84) %>%
    slice_head(n = n_top) %>%
    mutate(country = country, model = model, .before = 1)
}

# ------------------------------------------------------------------------------
# 1. Localizar y leer la base sin modificarla
# ------------------------------------------------------------------------------

explicit_sav <- Sys.getenv("SEM_SAV_PATH", "")
if (nzchar(explicit_sav) && !file.exists(explicit_sav)) {
  stop("SEM_SAV_PATH no existe: ", explicit_sav, call. = FALSE)
}
exact_candidate <- file.path(project_dir, "upload", "Base_Final(20260806-032528).sav")
search_candidates <- unique(c(
  exact_candidate,
  file.path(project_dir, "Base_Final.sav"),
  file.path(project_dir, "upload", "Base_Final.sav"),
  file.path(getwd(), "Base_Final.sav"),
  file.path(getwd(), "upload", "Base_Final.sav")
))
search_candidates <- search_candidates[file.exists(search_candidates)]
if (nzchar(explicit_sav)) {
  sav_path <- normalizePath(explicit_sav, mustWork = TRUE)
} else if (length(search_candidates) == 1L) {
  sav_path <- normalizePath(search_candidates[[1L]], mustWork = TRUE)
} else if (length(search_candidates) == 0L) {
  stop("No se encontro la base .sav. Defina SEM_SAV_PATH de manera explicita.", call. = FALSE)
} else {
  stop("Varias bases candidatas: ", paste(search_candidates, collapse = " | "),
       ". Defina SEM_SAV_PATH para no elegir una version por prioridad arbitraria.", call. = FALSE)
}
# Manifiesto inicial: aun si una etapa posterior falla, conserva fuente exacta.
safe_csv(tibble(kind = c("script", "data"), path = c(script_path, sav_path),
                md5 = unname(tools::md5sum(c(script_path, sav_path)))),
         file.path(results_dir, "source_hashes.csv"))
data_raw <- haven::read_sav(sav_path)
original_names <- names(data_raw)
# En algunas plataformas iconv translitera las vocales acentuadas con un
# apostrofo (p. ej., "País" -> "Pa'is"). Para los nombres analiticos se
# conservan solo letras, numeros y guion bajo; las etiquetas originales quedan
# intactas en el diccionario exportado.
safe_names <- gsub("[^A-Za-z0-9_]+", "", ascii_text(original_names))
if (anyNA(safe_names) || anyDuplicated(safe_names)) {
  stop("La transliteracion de nombres produjo NA o duplicados; revise el diccionario.",
       call. = FALSE)
}
names(data_raw) <- safe_names

value_labels_text <- function(x) {
  labs <- attr(x, "labels", exact = TRUE)
  if (is.null(labs)) return(NA_character_)
  paste0(names(labs), "=", unname(labs), collapse = "; ")
}

dictionary <- map_dfr(seq_along(data_raw), function(i) {
  x <- data_raw[[i]]
  tibble(
    position = i,
    original_name = original_names[i],
    analysis_name = names(data_raw)[i],
    variable_label = attr(x, "label", exact = TRUE) %||% NA_character_,
    storage_class = paste(class(x), collapse = ";"),
    value_labels = value_labels_text(x)
  )
})
safe_csv(dictionary, file.path(audit_dir, "variable_dictionary.csv"))

cat("Base leida desde: ", sav_path, "\n", sep = "")
cat("Dimensiones: ", nrow(data_raw), " filas x ", ncol(data_raw), " variables\n", sep = "")

# ------------------------------------------------------------------------------
# 2. Definicion de variables, recodificacion y auditoria de calidad
# ------------------------------------------------------------------------------

pp_items <- c("PP1", "PP2", "PP3")
sc_items <- c(
  "APS1", "APS2", "APS3", "APS4", "APS5", "APS6", "APS7", "APS9",
  "APS10", "APS14"
)
sbl_items <- c(
  "APS11", "APS12", "APS13", "APS15", "APS16", "APS18", "APS19",
  "APS20", "APS21"
)
aps_items <- c(sc_items, sbl_items)
value_items <- c("Benevolencia", "Universalismo", "Poder", "Logro")
indicator_items <- c(sc_items, sbl_items, pp_items)
all_model_items <- c(indicator_items, value_items)
control_vars <- c("Sexo", "Edad", "Estudios", "NSE")

required_vars <- unique(c(
  "Pais", "CSR1_rec", all_model_items, control_vars
))
missing_required <- setdiff(required_vars, names(data_raw))
if (length(missing_required) > 0L) {
  # CSR1_rec puede sustituirse por CSR1; las demas variables son obligatorias.
  missing_without_csr <- setdiff(missing_required, "CSR1_rec")
  if (length(missing_without_csr) > 0L || !"CSR1" %in% names(data_raw)) {
    stop(
      "Faltan variables requeridas: ", paste(missing_required, collapse = ", "),
      call. = FALSE
    )
  }
}

# Mapeo semantico a partir de las etiquetas reales, antes de eliminar atributos.
# Si las etiquetas no permiten un mapeo inequívoco, solo se acepta un 0/1 sin etiquetas.
binary_recode_rows <- list()
recode_binary <- function(x, variable_name, positive, negative) {
  raw <- if (is.factor(x) || is.character(x)) as.character(x) else as.numeric(x)
  labels <- attr(x, "labels", exact = TRUE)
  if (!is.null(labels)) {
    label_tokens <- normalize_token(names(labels))
    label_map <- ifelse(label_tokens %in% positive, 1,
                        ifelse(label_tokens %in% negative, 0, NA_real_))
    pos <- match(as.numeric(x), as.numeric(labels))
    mapped <- as.numeric(label_map[pos])
    if (any(!is.na(x) & is.na(mapped))) {
      stop(variable_name, ": etiquetas no reconocidas o codigos sin etiqueta; no se asume 1/2 ni 0/1.", call. = FALSE)
    }
    binary_recode_rows[[variable_name]] <<- tibble(variable = variable_name,
      original_code = as.character(labels), original_label = names(labels),
      normalized_label = label_tokens, recoded = label_map)
    return(mapped)
  }
  if (is.factor(x) || is.character(x)) {
    tokens <- normalize_token(raw)
    mapped <- ifelse(tokens %in% positive, 1,
                     ifelse(tokens %in% negative, 0, NA_real_))
    if (any(!is.na(x) & is.na(mapped))) stop(variable_name, ": categorias textuales ambiguas.", call. = FALSE)
    return(as.numeric(mapped))
  }
  observed <- sort(unique(raw[!is.na(raw)]))
  if (length(observed) && all(observed %in% c(0, 1))) {
    binary_recode_rows[[variable_name]] <<- tibble(variable = variable_name,
      original_code = as.character(observed), original_label = NA_character_,
      normalized_label = "unlabelled_0_1", recoded = observed)
    return(as.numeric(raw))
  }
  stop(variable_name, ": no puede recodificarse sin etiquetas validas; codigos: ",
       paste(observed, collapse = ", "), call. = FALSE)
}
recode_yes_no <- function(x, variable_name) recode_binary(x, variable_name,
  positive = c("si", "yes", "true", "afirmativo", "1"),
  negative = c("no", "false", "negativo", "0"))

recode_country <- function(x) {
  if (inherits(x, "haven_labelled")) {
    out <- as.character(haven::as_factor(x, levels = "labels"))
  } else if (is.factor(x) || is.character(x)) {
    out <- as.character(x)
  } else {
    x_num <- suppressWarnings(as.numeric(x))
    out <- dplyr::recode(
      x_num, `1` = "Argentina", `2` = "Colombia", `3` = "Peru", `4` = "Mexico",
      .default = NA_character_
    )
  }
  out <- ascii_text(out)
  token <- normalize_token(out)
  out <- dplyr::case_when(
    token == "argentina" ~ "Argentina",
    token == "colombia" ~ "Colombia",
    token == "peru" ~ "Peru",
    token == "mexico" ~ "Mexico",
    TRUE ~ NA_character_
  )
  factor(out, levels = c("Argentina", "Colombia", "Peru", "Mexico"))
}

data <- data_raw
data$Pais <- recode_country(data$Pais)
csr_source <- if ("CSR1_rec" %in% names(data)) "CSR1_rec" else "CSR1"
data$CSR1_bin <- recode_yes_no(data[[csr_source]], csr_source)
if (all(c("CSR1", "CSR1_rec") %in% names(data))) {
  source_csr <- recode_yes_no(data$CSR1, "CSR1")
  recoded_csr <- recode_yes_no(data$CSR1_rec, "CSR1_rec")
  disagreements <- sum(source_csr != recoded_csr, na.rm = TRUE)
  safe_csv(tibble(n = nrow(data), disagreements = disagreements,
                  missing_disagreements = sum(xor(is.na(source_csr), is.na(recoded_csr)))),
           file.path(audit_dir, "csr1_source_recode_consistency.csv"))
  if (disagreements || any(xor(is.na(source_csr), is.na(recoded_csr)))) {
    stop("CSR1 y CSR1_rec no coinciden despues de mapear etiquetas; revise el origen.", call. = FALSE)
  }
}
data$Sexo <- recode_binary(data$Sexo, "Sexo",
  c("mujer", "femenino", "female", "woman", "1"),
  c("hombre", "masculino", "male", "man", "0"))
data$Estudios <- recode_binary(data$Estudios, "Estudios",
  c("postgrado", "posgrado", "postgraduate", "graduate", "1"),
  c("pregrado", "undergraduate", "0"))
safe_csv(bind_rows(binary_recode_rows), file.path(audit_dir, "binary_recode_dictionary.csv"))

numeric_vars <- unique(c(all_model_items, control_vars))
for (v in numeric_vars) data[[v]] <- suppressWarnings(as.numeric(data[[v]]))

if (anyNA(data$Pais)) {
  stop("Pais contiene valores perdidos o no reconocidos despues de recodificar.",
       call. = FALSE)
}
csr_observed_levels <- sort(unique(na.omit(data$CSR1_bin)))
if (length(csr_observed_levels) == 0L ||
    !all(csr_observed_levels %in% c(0, 1))) {
  stop("CSR1_bin contiene codigos diferentes de 0/1.", call. = FALSE)
}
if (length(csr_observed_levels) < 2L) {
  warning("CSR1_bin solo presenta un nivel; sus regresiones no seran identificables.")
}

# Las etiquetas originales se preservan en data_raw y en el diccionario.
# Sexo y Estudios se mapean semánticamente a 0/1; edad/NSE se estandarizan.
# No se imputan valores.
data <- data %>%
  mutate(
    Sexo_num = as.numeric(Sexo),
    Estudios_num = as.numeric(Estudios),
    Edad_z = as.numeric(scale(Edad)),
    NSE_z = as.numeric(scale(NSE))
  )

expected_ranges <- c(
  setNames(rep("1:5", length(indicator_items)), indicator_items),
  setNames(rep("1:6", length(value_items)), value_items),
  CSR1_bin = "0:1",
  Sexo = "0:1",
  Estudios = "0:1",
  NSE = "1:5"
)

range_audit <- imap_dfr(expected_ranges, function(spec, variable) {
  expected <- eval(parse(text = spec), envir = baseenv())
  x <- data[[variable]]
  observed <- sort(unique(x[!is.na(x)]))
  invalid <- observed[!observed %in% expected]
  tibble(
    variable = variable,
    expected_codes = paste(expected, collapse = ","),
    n_nonmissing = sum(!is.na(x)),
    n_missing = sum(is.na(x)),
    observed_min = if (length(observed)) min(observed) else NA_real_,
    observed_max = if (length(observed)) max(observed) else NA_real_,
    invalid_codes = if (length(invalid)) paste(invalid, collapse = ",") else NA_character_,
    n_invalid = sum(!is.na(x) & !x %in% expected)
  )
})
safe_csv(range_audit, file.path(audit_dir, "range_and_missing_audit.csv"))

critical_range_failure <- any(range_audit$n_invalid > 0L)
if (critical_range_failure) {
  stop(
    "Se detectaron codigos fuera de rango. Revise results_sem/01_audit/",
    "range_and_missing_audit.csv antes de estimar modelos.",
    call. = FALSE
  )
}

sample_by_country <- data %>%
  count(Pais, name = "n") %>%
  mutate(proportion = n / sum(n))
safe_csv(sample_by_country, file.path(audit_dir, "sample_by_country.csv"))

csr_by_country <- data %>%
  count(Pais, CSR1_bin, name = "n") %>%
  group_by(Pais) %>%
  mutate(proportion_within_country = n / sum(n)) %>%
  ungroup()
safe_csv(csr_by_country, file.path(audit_dir, "csr1_by_country.csv"))

csr_variation_by_country <- data %>%
  group_by(Pais) %>%
  summarise(
    observed_levels = n_distinct(CSR1_bin, na.rm = TRUE),
    both_levels_present = observed_levels == 2L,
    .groups = "drop"
  )
safe_csv(
  csr_variation_by_country,
  file.path(audit_dir, "csr1_variation_by_country.csv")
)
if (!all(csr_variation_by_country$both_levels_present)) {
  warning("CSR1 no presenta ambos niveles en todos los paises; revise rutas por grupo.")
}

category_frequencies <- data %>%
  select(Pais, all_of(all_model_items)) %>%
  pivot_longer(-Pais, names_to = "item", values_to = "response") %>%
  count(Pais, item, response, .drop = FALSE, name = "n") %>%
  group_by(Pais, item) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()
safe_csv(category_frequencies, file.path(audit_dir, "item_category_frequencies.csv"))

missing_by_country <- data %>%
  select(Pais, all_of(c(all_model_items, "CSR1_bin", control_vars))) %>%
  pivot_longer(-Pais, names_to = "variable", values_to = "value") %>%
  group_by(Pais, variable) %>%
  summarise(
    n = n(), n_missing = sum(is.na(value)), proportion_missing = mean(is.na(value)),
    .groups = "drop"
  )
safe_csv(missing_by_country, file.path(audit_dir, "missing_by_country.csv"))

# Casos analiticos: se exige informacion completa para que MLR, WLSMV y las
# comparaciones anidadas utilicen exactamente las mismas observaciones.
analysis_vars <- c(
  "Pais", "CSR1_bin", all_model_items, "Sexo_num", "Estudios_num", "Edad_z", "NSE_z"
)
complete_flag <- complete.cases(data[, analysis_vars])
flow_table <- tibble(
  stage = c("Base original", "Casos completos en variables analiticas"),
  n = c(nrow(data), sum(complete_flag)),
  excluded = c(0L, sum(!complete_flag))
)
safe_csv(flow_table, file.path(audit_dir, "analytic_sample_flow.csv"))
if (mean(!complete_flag) > .05) {
  warning(
    "Mas de 5% de casos quedaria fuera por datos faltantes. Revise el mecanismo ",
    "de perdida antes de continuar."
  )
}
dat <- droplevels(data[complete_flag, , drop = FALSE])
saveRDS(dat, file.path(audit_dir, "analytic_data.rds"))

# Auditoria etaria. La presencia de menores no es un criterio estadistico de
# exclusion; exige verificar la aprobacion etica, el asentimiento y, cuando
# corresponda, el consentimiento parental o la exencion aplicable.
minor_audit <- data %>%
  mutate(is_minor = Edad < 18) %>%
  group_by(Pais) %>%
  summarise(
    n = n(),
    n_minors = sum(is_minor, na.rm = TRUE),
    proportion_minors = mean(is_minor, na.rm = TRUE),
    minimum_age = min(Edad, na.rm = TRUE),
    .groups = "drop"
  )
safe_csv(minor_audit, file.path(audit_dir, "participants_under_18_by_country.csv"))
safe_csv(
  data %>%
    filter(Edad < 18) %>%
    transmute(
      participant_id = if ("Codigo" %in% names(data)) as.character(Codigo) else NA_character_,
      Pais, Edad, Sexo = Sexo_num, Estudios = Estudios_num
    ),
  file.path(audit_dir, "participants_under_18_records.csv")
)

# Auditoria de duplicados. Se compara el perfil completo de las 99 variables
# originales restantes, no solo los indicadores del SEM. Se marcan tanto todos
# los miembros de cada conjunto como las copias posteriores; nada se elimina en
# el analisis principal.
id_var <- if ("Codigo" %in% names(dat)) "Codigo" else NA_character_
duplicate_id <- rep(FALSE, nrow(dat))
if (!is.na(id_var)) {
  duplicate_id <- duplicated(dat[[id_var]]) | duplicated(dat[[id_var]], fromLast = TRUE)
}
profile_vars <- setdiff(safe_names, id_var)
profile_key <- apply(
  dat[, profile_vars, drop = FALSE],
  1,
  paste,
  collapse = "|"
)
duplicate_profile_member <- duplicated(profile_key) |
  duplicated(profile_key, fromLast = TRUE)
duplicate_profile_copy <- duplicated(profile_key)

row_sd <- function(df) apply(df, 1, stats::sd, na.rm = TRUE)
max_response_share <- function(df) {
  apply(df, 1, function(z) {
    z <- z[!is.na(z)]
    if (!length(z)) return(NA_real_)
    max(table(z)) / length(z)
  })
}

quality_flags <- tibble(
  row_number = seq_len(nrow(dat)),
  participant_id = if (!is.na(id_var)) as.character(dat[[id_var]]) else NA_character_,
  Pais = dat$Pais,
  duplicate_id = duplicate_id,
  duplicate_full_profile_member = duplicate_profile_member,
  duplicate_full_profile_copy = duplicate_profile_copy,
  aps_zero_within_person_sd = row_sd(dat[, aps_items, drop = FALSE]) == 0,
  psb_zero_within_person_sd = row_sd(dat[, pp_items, drop = FALSE]) == 0,
  sc_sbl_psb_zero_within_person_sd =
    row_sd(dat[, indicator_items, drop = FALSE]) == 0,
  values_zero_within_person_sd = row_sd(dat[, value_items, drop = FALSE]) == 0,
  aps_same_response_share = max_response_share(dat[, aps_items, drop = FALSE]),
  aps_extreme_straightlining = aps_same_response_share >= .95
)
safe_csv(quality_flags, file.path(audit_dir, "quality_flags_no_auto_exclusion.csv"))

quality_summary <- quality_flags %>%
  summarise(
    n = n(),
    duplicate_ids = sum(duplicate_id),
    duplicate_full_profile_members = sum(duplicate_full_profile_member),
    duplicate_full_profile_copies = sum(duplicate_full_profile_copy),
    aps_zero_sd = sum(aps_zero_within_person_sd),
    psb_zero_sd = sum(psb_zero_within_person_sd),
    sc_sbl_psb_zero_sd = sum(sc_sbl_psb_zero_within_person_sd),
    values_zero_sd = sum(values_zero_within_person_sd),
    aps_extreme_straightlining = sum(aps_extreme_straightlining)
  )
safe_csv(quality_summary, file.path(audit_dir, "quality_summary.csv"))

sensitivity_scenarios <- bind_rows(
  tibble(scenario = "all_cases", keep = rep(TRUE, nrow(dat))),
  tibble(scenario = "exclude_duplicate_copies", keep = !duplicate_profile_copy),
  tibble(
    scenario = "exclude_invariant_SC_SBL_PSB_block",
    keep = !quality_flags$sc_sbl_psb_zero_within_person_sd
  ),
  tibble(
    scenario = "exclude_duplicates_and_invariant_block",
    keep = !duplicate_profile_copy &
      !quality_flags$sc_sbl_psb_zero_within_person_sd
  )
) %>%
  mutate(Pais = rep(dat$Pais, times = 4L)) %>%
  group_by(scenario, Pais) %>%
  summarise(n_retained = sum(keep), n_excluded = sum(!keep), .groups = "drop")
safe_csv(
  sensitivity_scenarios,
  file.path(audit_dir, "predefined_sensitivity_scenario_counts.csv")
)

# Verificacion de las dummies de CSR2, aunque no entren en este SEM. Esta
# auditoria evita que una variable derivada incorrecta se reutilice en analisis
# posteriores.
if (all(c("CSR2", "CR2_Si", "CR2_No") %in% names(data))) {
  csr2_num <- suppressWarnings(as.numeric(data$CSR2))
  cr2_si <- suppressWarnings(as.numeric(data$CR2_Si))
  cr2_no <- suppressWarnings(as.numeric(data$CR2_No))
  csr2_dummy_check <- tibble(
    variable = c("CR2_Si", "CR2_No"),
    discrepancies_vs_CSR2 = c(
      sum(cr2_si != as.numeric(csr2_num == 1), na.rm = TRUE),
      sum(cr2_no != as.numeric(csr2_num == 2), na.rm = TRUE)
    )
  )
  safe_csv(csr2_dummy_check, file.path(audit_dir, "csr2_dummy_integrity_check.csv"))
}

# Verificacion de las puntuaciones preexistentes. Estas comprobaciones no alteran
# la especificacion; sirven para detectar discrepancias entre versiones.
score_checks <- list()
if ("Prosocialidad" %in% names(data)) {
  old_psb <- suppressWarnings(as.numeric(data$Prosocialidad))
  score_checks[["Prosocialidad_equals_mean_PP1_PP2"]] <- isTRUE(all.equal(
    old_psb,
    rowMeans(data[, c("PP1", "PP2")], na.rm = FALSE),
    check.attributes = FALSE
  ))
}
if ("VAL_AutoTrascendencia" %in% names(data)) {
  st_old <- suppressWarnings(as.numeric(data$VAL_AutoTrascendencia))
  score_checks[["VAL_ST_equals_mean_Benevolencia_Universalismo"]] <- isTRUE(all.equal(
    st_old,
    rowMeans(data[, c("Benevolencia", "Universalismo")], na.rm = FALSE),
    check.attributes = FALSE
  ))
  if ("Quijotismo" %in% names(data)) {
    score_checks[["VAL_ST_equals_mean_Benevolencia_Universalismo_Quijotismo"]] <-
      isTRUE(all.equal(
        st_old,
        rowMeans(data[, c("Benevolencia", "Universalismo", "Quijotismo")],
                 na.rm = FALSE),
        check.attributes = FALSE
      ))
  }
}
if ("VAL_AutoPromocion" %in% names(data)) {
  se_old <- suppressWarnings(as.numeric(data$VAL_AutoPromocion))
  score_checks[["VAL_SE_equals_mean_Poder_Logro"]] <- isTRUE(all.equal(
    se_old,
    rowMeans(data[, c("Poder", "Logro")], na.rm = FALSE),
    check.attributes = FALSE
  ))
}
if (length(score_checks)) {
  safe_csv(
    tibble::enframe(unlist(score_checks), name = "check", value = "result"),
    file.path(audit_dir, "preexisting_score_checks.csv")
  )
}

# ------------------------------------------------------------------------------
# 3. Modelos de medicion: modelo teorico y rivales
# ------------------------------------------------------------------------------

measurement_model_3f <- paste0(
  "SC =~ ", paste(sc_items, collapse = " + "), "\n",
  "SBL =~ ", paste(sbl_items, collapse = " + "), "\n",
  "PSB =~ ", paste(pp_items, collapse = " + ")
)

# Rival 1: SC y SBL colapsados en una sola dimension ASP; PSB permanece separado.
measurement_model_asp1_psb <- paste0(
  "ASP =~ ", paste(aps_items, collapse = " + "), "\n",
  "PSB =~ ", paste(pp_items, collapse = " + ")
)

# Rival 2: SBL y PSB colapsados; SC permanece separado. Diagnostica el posible
# solapamiento de la disposicion de liderazgo con la prosocialidad autorreportada.
measurement_model_sc_sblpsb <- paste0(
  "SC =~ ", paste(sc_items, collapse = " + "), "\n",
  "SBL_PSB =~ ", paste(c(sbl_items, pp_items), collapse = " + ")
)

# Rival 3: un unico factor comun para los 22 indicadores.
measurement_model_1f <- paste0(
  "GENERAL =~ ", paste(indicator_items, collapse = " + ")
)

measurement_models <- c(
  list(theoretical_3factor = measurement_model_3f),
  if (run_rivals) {
    list(
      rival_ASP1_plus_PSB = measurement_model_asp1_psb,
      rival_SC_plus_SBLPSB = measurement_model_sc_sblpsb,
      rival_all_1factor = measurement_model_1f
    )
  } else {
    list()
  }
)

safe_text(
  c(
    "MODELO TEORICO DE TRES FACTORES", measurement_model_3f, "",
    "RIVAL: ASP UNIDIMENSIONAL + PSB", measurement_model_asp1_psb, "",
    "RIVAL: SC + SBL/PSB COMBINADOS", measurement_model_sc_sblpsb, "",
    "RIVAL: FACTOR GENERAL", measurement_model_1f
  ),
  file.path(cfa_dir, "measurement_model_syntax.txt")
)

fit_cfa_safe <- function(model, data, estimator = "MLR", ordered = NULL) {
  capture_conditions({
    args <- list(
      model = model,
      data = data,
      estimator = estimator,
      std.lv = TRUE,
      meanstructure = TRUE
    )
    if (estimator == "MLR") {
      args$missing <- "fiml"
    } else if (estimator == "WLSMV") {
      args$ordered <- ordered
      args$parameterization <- "theta"
      args$missing <- "pairwise"
    }
    do.call(lavaan::cfa, args)
  })
}

safe_lrt <- function(restricted_fit, full_fit, comparison, country) {
  if (!is_lavaan_fit(restricted_fit) || !is_lavaan_fit(full_fit)) return(tibble())
  ans <- capture_conditions(
    lavaan::lavTestLRT(
      restricted_fit, full_fit, method = "satorra.bentler.2010"
    )
  )
  if (!is.na(ans$error)) {
    return(tibble(
      country = country, comparison = comparison, error = ans$error
    ))
  }
  as.data.frame(ans$value) %>%
    tibble::rownames_to_column("model_row") %>%
    as_tibble() %>%
    mutate(
      country = country,
      comparison = comparison,
      warnings = paste(ans$warnings, collapse = " | "),
      .before = 1
    )
}

construct_blocks <- list(SC = sc_items, SBL = sbl_items, PSB = pp_items)

compute_cr_ave <- function(fit, country, model = "theoretical_3factor") {
  if (!is_lavaan_fit(fit)) return(tibble())
  ss <- standardizedSolution(fit, type = "std.all", se = FALSE, zstat = FALSE, pvalue = FALSE, ci = FALSE) %>% as_tibble()
  map_dfr(names(construct_blocks), function(latent) {
    items <- construct_blocks[[latent]]
    lambda <- ss %>%
      filter(op == "=~", lhs == latent, rhs %in% items) %>%
      pull(est.std)
    theta <- ss %>%
      filter(op == "~~", lhs == rhs, lhs %in% items) %>%
      arrange(match(lhs, items)) %>%
      pull(est.std)
    if (length(lambda) != length(items) || length(theta) != length(items)) {
      return(tibble(
        country = country, model = model, construct = latent,
        n_indicators = length(items), composite_reliability = NA_real_,
        ave = NA_real_, admissible_components = FALSE
      ))
    }
    denom_cr <- sum(lambda)^2 + sum(theta)
    denom_ave <- sum(lambda^2) + sum(theta)
    tibble(
      country = country,
      model = model,
      construct = latent,
      n_indicators = length(items),
      composite_reliability = sum(lambda)^2 / denom_cr,
      ave = sum(lambda^2) / denom_ave,
      admissible_components = all(is.finite(c(lambda, theta))) && all(theta >= 0)
    )
  })
}

cronbach_alpha <- function(df) {
  # Formula cruda directa. Evita reconstruir el objeto completo de psych::alpha
  # en cada replica bootstrap y produce el mismo alfa sin inversion automatica.
  x <- as.matrix(df)
  storage.mode(x) <- "double"
  k <- ncol(x)
  if (k < 2L) return(NA_real_)
  item_variance <- sum(apply(x, 2, stats::var, na.rm = TRUE))
  total_variance <- stats::var(rowSums(x, na.rm = FALSE), na.rm = TRUE)
  if (!is.finite(total_variance) || total_variance <= 0) return(NA_real_)
  k / (k - 1) * (1 - item_variance / total_variance)
}

bootstrap_alpha <- function(df, b = 1000L, seed = 20260805L) {
  point <- cronbach_alpha(df)
  if (b <= 0L) {
    return(c(alpha = point, ci_low = NA_real_, ci_high = NA_real_, valid_replicates = 0))
  }
  set.seed(seed)
  vals <- replicate(b, {
    rows <- sample.int(nrow(df), nrow(df), replace = TRUE)
    tryCatch(cronbach_alpha(df[rows, , drop = FALSE]), error = function(e) NA_real_)
  })
  valid <- is.finite(vals)
  enough <- sum(valid) >= .95 * b
  ci <- if (enough) stats::quantile(vals[valid], probs = c(.025, .975), names = FALSE) else rep(NA_real_, 2L)
  c(alpha = point, ci_low = ci[1], ci_high = ci[2], valid_replicates = sum(valid))
}

htmt_from_cor <- function(r, blocks = construct_blocks) {
  r <- as.matrix(r)
  pairs <- combn(names(blocks), 2, simplify = FALSE)
  map_dfr(pairs, function(pair) {
    a <- blocks[[pair[1]]]
    b <- blocks[[pair[2]]]
    cross <- abs(r[a, b, drop = FALSE])
    mono_a <- abs(r[a, a, drop = FALSE][upper.tri(r[a, a, drop = FALSE])])
    mono_b <- abs(r[b, b, drop = FALSE][upper.tri(r[b, b, drop = FALSE])])
    denom <- sqrt(mean(mono_a, na.rm = TRUE) * mean(mono_b, na.rm = TRUE))
    tibble(
      construct_1 = pair[1], construct_2 = pair[2],
      htmt = mean(cross, na.rm = TRUE) / denom
    )
  })
}

bootstrap_htmt_pearson <- function(df, country, b = 1000L, seed = 20260805L) {
  point <- htmt_from_cor(stats::cor(df, use = "pairwise.complete.obs"))
  if (b <= 0L) {
    return(point %>% mutate(ci_low = NA_real_, ci_high = NA_real_, country = country,
                           requested_replicates = b, valid_replicates = 0L, valid_fraction = NA_real_))
  }
  set.seed(seed)
  boot_list <- vector("list", b)
  for (i in seq_len(b)) {
    rows <- sample.int(nrow(df), nrow(df), replace = TRUE)
    boot_list[[i]] <- tryCatch(
      htmt_from_cor(stats::cor(df[rows, , drop = FALSE], use = "pairwise.complete.obs")) %>%
        mutate(rep = i),
      error = function(e) NULL
    )
  }
  boot <- bind_rows(boot_list)
  cis <- boot %>%
    group_by(construct_1, construct_2) %>%
    summarise(
      requested_replicates = b,
      valid_replicates = sum(is.finite(htmt)),
      valid_fraction = valid_replicates / b,
      ci_low = if (valid_fraction >= .95) quantile(htmt[is.finite(htmt)], .025) else NA_real_,
      ci_high = if (valid_fraction >= .95) quantile(htmt[is.finite(htmt)], .975) else NA_real_,
      .groups = "drop"
    )
  point %>%
    left_join(cis, by = c("construct_1", "construct_2")) %>%
    mutate(country = country, method = "Pearson bootstrap", .before = 1)
}

ordinal_htmt <- function(df, country) {
  poly <- capture_conditions(psych::polychoric(df, correct = .5, progress = FALSE))
  if (!is.na(poly$error)) {
    return(tibble(
      country = country, method = "Polychoric point estimate",
      construct_1 = NA_character_, construct_2 = NA_character_, htmt = NA_real_,
      error = poly$error, warnings = paste(poly$warnings, collapse = " | ")
    ))
  }
  htmt_from_cor(poly$value$rho) %>%
    mutate(
      country = country,
      method = "Polychoric point estimate",
      error = NA_character_,
      warnings = paste(poly$warnings, collapse = " | "),
      .before = 1
    )
}

alpha_table_for_data <- function(df, country, b = boot_htmt) {
  imap_dfr(construct_blocks, function(items, construct) {
    out <- bootstrap_alpha(df[, items, drop = FALSE], b = b,
                           seed = 20260805L + match(construct, names(construct_blocks)))
    tibble(
      country = country,
      construct = construct,
      alpha = unname(out["alpha"]),
      ci_low = unname(out["ci_low"]),
      ci_high = unname(out["ci_high"]),
      bootstrap_replicates = b,
      valid_replicates = as.integer(out["valid_replicates"]),
      valid_fraction = if (b > 0) as.numeric(out["valid_replicates"]) / b else NA_real_
    )
  })
}

# Ajuste pooled del modelo de medicion y de los rivales.
pooled_cfa_results <- imap(measurement_models, function(model_syntax, model_name) {
  message("CFA MLR pooled: ", model_name)
  fit_cfa_safe(model_syntax, dat, estimator = "MLR")
})
pooled_cfa_fits <- map(pooled_cfa_results, "value")

country_levels <- levels(dat$Pais)
country_cfa_results <- setNames(vector("list", length(country_levels)), country_levels)
country_wlsmv_results <- setNames(vector("list", length(country_levels)), country_levels)

for (country in country_levels) {
  d_country <- droplevels(dat[dat$Pais == country, , drop = FALSE])
  country_cfa_results[[country]] <- imap(
    measurement_models,
    function(model_syntax, model_name) {
      message("CFA MLR ", country, ": ", model_name)
      fit_cfa_safe(model_syntax, d_country, estimator = "MLR")
    }
  )

  # WLSMV se usa como sensibilidad y nunca colapsa categorias automaticamente.
  # Las celdas vacias/escasas quedan registradas como advertencias o errores.
  if (run_ordinal) {
    message("CFA WLSMV ", country, ": theoretical_3factor")
    country_wlsmv_results[[country]] <- fit_cfa_safe(
      measurement_model_3f,
      d_country,
      estimator = "WLSMV",
      ordered = indicator_items
    )
  }
}

# ------------------------------------------------------------------------------
# 4. Exportar ajuste, diagnosticos y validez de los CFA
# ------------------------------------------------------------------------------

cfa_fit_rows <- list()
cfa_status_rows <- list()
cfa_warning_rows <- list()
cfa_parameter_rows <- list()
cfa_r2_rows <- list()
cfa_latent_cor_rows <- list()
cfa_residual_rows <- list()
cfa_mi_rows <- list()
cfa_cr_ave_rows <- list()
cfa_alpha_rows <- list()
cfa_htmt_rows <- list()
cfa_htmt_ordinal_rows <- list()
cfa_lrt_rows <- list()
checkpoint_cfa <- function(stage) {
  checkpoint_dir <- file.path(cfa_dir, "checkpoints")
  safe_csv(bind_rows(cfa_fit_rows), file.path(checkpoint_dir, "cfa_fit_progress.csv"))
  safe_csv(bind_rows(cfa_status_rows), file.path(checkpoint_dir, "cfa_admissibility_progress.csv"))
  safe_csv(bind_rows(cfa_warning_rows), file.path(checkpoint_dir, "cfa_errors_progress.csv"))
  safe_text(paste0("PARTIAL CHECKPOINT: ", stage, "; not final until run_status=COMPLETED"),
            file.path(checkpoint_dir, "PARTIAL_STATUS.txt"))
}

for (model_name in names(pooled_cfa_results)) {
  message("Exportando diagnosticos CFA pooled: ", model_name)
  res <- pooled_cfa_results[[model_name]]
  fit <- res$value
  cfa_fit_rows[[length(cfa_fit_rows) + 1L]] <-
    extract_fit_row(fit, model_name, "MLR", "Pooled")
  cfa_status_rows[[length(cfa_status_rows) + 1L]] <-
    extract_admissibility(fit, "Pooled", model_name, "MLR")
  if (length(res$warnings) || !is.na(res$error)) {
    cfa_warning_rows[[length(cfa_warning_rows) + 1L]] <- tibble(
      country = "Pooled", model = model_name, estimator = "MLR",
      error = res$error, warnings = paste(res$warnings, collapse = " | ")
    )
  }
}
checkpoint_cfa("pooled fit/status exported")

pooled_theoretical <- pooled_cfa_fits$theoretical_3factor
if (is_lavaan_fit(pooled_theoretical)) {
  cfa_parameter_rows[[length(cfa_parameter_rows) + 1L]] <-
    extract_parameters(pooled_theoretical, "Pooled", "theoretical_3factor") %>%
    filter(op %in% c("=~", "~~"))
  cfa_r2_rows[[length(cfa_r2_rows) + 1L]] <-
    extract_r2(pooled_theoretical, "Pooled", "theoretical_3factor")
  cfa_latent_cor_rows[[length(cfa_latent_cor_rows) + 1L]] <-
    extract_latent_correlations(pooled_theoretical, "Pooled", "theoretical_3factor")
  cfa_residual_rows[[length(cfa_residual_rows) + 1L]] <-
    extract_top_residual_correlations(
      pooled_theoretical, "Pooled", "theoretical_3factor"
    )
  cfa_mi_rows[[length(cfa_mi_rows) + 1L]] <-
    extract_modification_indices(pooled_theoretical, "Pooled", "theoretical_3factor")
  cfa_cr_ave_rows[[length(cfa_cr_ave_rows) + 1L]] <-
    compute_cr_ave(pooled_theoretical, "Pooled")
  for (rival_name in setdiff(names(pooled_cfa_results), "theoretical_3factor")) {
    cfa_lrt_rows[[length(cfa_lrt_rows) + 1L]] <- safe_lrt(
      pooled_cfa_results[[rival_name]]$value,
      pooled_theoretical,
      paste0(rival_name, " vs theoretical_3factor"),
      "Pooled"
    )
  }
}
cfa_alpha_rows[[1L]] <- alpha_table_for_data(dat, "Pooled")
message("Bootstrap HTMT pooled: ", boot_htmt, " replicas")
cfa_htmt_rows[[1L]] <- bootstrap_htmt_pearson(
  dat[, indicator_items, drop = FALSE], "Pooled", b = boot_htmt
)
if (run_ordinal) {
  cfa_htmt_ordinal_rows[[1L]] <- ordinal_htmt(
    dat[, indicator_items, drop = FALSE], "Pooled"
  )
}

for (country in country_levels) {
  message("Diagnosticos/HTMT/policoricas CFA: ", country)
  d_country <- droplevels(dat[dat$Pais == country, , drop = FALSE])
  country_results <- country_cfa_results[[country]]

  for (model_name in names(country_results)) {
    res <- country_results[[model_name]]
    fit <- res$value
    cfa_fit_rows[[length(cfa_fit_rows) + 1L]] <-
      extract_fit_row(fit, model_name, "MLR", country)
    cfa_status_rows[[length(cfa_status_rows) + 1L]] <-
      extract_admissibility(fit, country, model_name, "MLR")
    if (length(res$warnings) || !is.na(res$error)) {
      cfa_warning_rows[[length(cfa_warning_rows) + 1L]] <- tibble(
        country = country, model = model_name, estimator = "MLR",
        error = res$error, warnings = paste(res$warnings, collapse = " | ")
      )
    }
  }

  checkpoint_cfa(paste0(country, " fit/status exported"))
  main_fit <- country_results$theoretical_3factor$value
  if (is_lavaan_fit(main_fit)) {
    cfa_parameter_rows[[length(cfa_parameter_rows) + 1L]] <-
      extract_parameters(main_fit, country, "theoretical_3factor") %>%
      filter(op %in% c("=~", "~~"))
    cfa_r2_rows[[length(cfa_r2_rows) + 1L]] <-
      extract_r2(main_fit, country, "theoretical_3factor")
    cfa_latent_cor_rows[[length(cfa_latent_cor_rows) + 1L]] <-
      extract_latent_correlations(main_fit, country, "theoretical_3factor")
    cfa_residual_rows[[length(cfa_residual_rows) + 1L]] <-
      extract_top_residual_correlations(main_fit, country, "theoretical_3factor")
    cfa_mi_rows[[length(cfa_mi_rows) + 1L]] <-
      extract_modification_indices(main_fit, country, "theoretical_3factor")
    cfa_cr_ave_rows[[length(cfa_cr_ave_rows) + 1L]] <-
      compute_cr_ave(main_fit, country)
  }

  cfa_alpha_rows[[length(cfa_alpha_rows) + 1L]] <-
    alpha_table_for_data(d_country, country)
  cfa_htmt_rows[[length(cfa_htmt_rows) + 1L]] <-
    bootstrap_htmt_pearson(
      d_country[, indicator_items, drop = FALSE], country, b = boot_htmt,
      seed = 20260805L + match(country, country_levels)
    )
  if (run_ordinal) {
    cfa_htmt_ordinal_rows[[length(cfa_htmt_ordinal_rows) + 1L]] <-
      ordinal_htmt(d_country[, indicator_items, drop = FALSE], country)
  }

  # Comparaciones teorico vs rivales bajo MLR. Las pruebas se registran solo si
  # lavaan reconoce correctamente la anidacion; AIC/BIC siempre quedan en la
  # tabla general de ajuste.
  for (rival_name in setdiff(names(country_results), "theoretical_3factor")) {
    cfa_lrt_rows[[length(cfa_lrt_rows) + 1L]] <- safe_lrt(
      country_results[[rival_name]]$value,
      main_fit,
      paste0(rival_name, " vs theoretical_3factor"),
      country
    )
  }

  # Sensibilidad WLSMV.
  if (run_ordinal) {
    wres <- country_wlsmv_results[[country]]
    wfit <- wres$value
    cfa_fit_rows[[length(cfa_fit_rows) + 1L]] <-
      extract_fit_row(wfit, "theoretical_3factor", "WLSMV", country)
    cfa_status_rows[[length(cfa_status_rows) + 1L]] <-
      extract_admissibility(wfit, country, "theoretical_3factor", "WLSMV")
    if (length(wres$warnings) || !is.na(wres$error)) {
      cfa_warning_rows[[length(cfa_warning_rows) + 1L]] <- tibble(
        country = country, model = "theoretical_3factor", estimator = "WLSMV",
        error = wres$error, warnings = paste(wres$warnings, collapse = " | ")
      )
    }
    if (is_lavaan_fit(wfit)) {
      cfa_parameter_rows[[length(cfa_parameter_rows) + 1L]] <-
        extract_parameters(wfit, country, "theoretical_3factor_WLSMV") %>%
        filter(op %in% c("=~", "~~", "|"))
    }
  }
}

cfa_fit_table <- bind_rows(cfa_fit_rows)
cfa_status_table <- bind_rows(cfa_status_rows)
cfa_parameter_table <- bind_rows(cfa_parameter_rows)
cfa_r2_table <- bind_rows(cfa_r2_rows)
cfa_latent_cor_table <- bind_rows(cfa_latent_cor_rows)
cfa_residual_table <- bind_rows(cfa_residual_rows)
cfa_mi_table <- bind_rows(cfa_mi_rows)
cfa_cr_ave_table <- bind_rows(cfa_cr_ave_rows)
cfa_alpha_table <- bind_rows(cfa_alpha_rows)
cfa_htmt_table <- bind_rows(cfa_htmt_rows)
cfa_htmt_ordinal_table <- bind_rows(cfa_htmt_ordinal_rows)
cfa_lrt_table <- bind_rows(cfa_lrt_rows)

safe_csv(cfa_fit_table, file.path(cfa_dir, "cfa_fit_complete.csv"))
safe_csv(cfa_status_table, file.path(cfa_dir, "cfa_admissibility_and_fit_gate.csv"))
safe_csv(bind_rows(cfa_warning_rows), file.path(cfa_dir, "cfa_errors_and_warnings.csv"))
safe_csv(cfa_parameter_table, file.path(cfa_dir, "cfa_parameters_se_ci.csv"))
safe_csv(cfa_r2_table, file.path(cfa_dir, "cfa_indicator_r2.csv"))
safe_csv(cfa_latent_cor_table, file.path(cfa_dir, "cfa_latent_correlations.csv"))
safe_csv(cfa_residual_table, file.path(cfa_dir, "cfa_largest_residual_correlations.csv"))
safe_csv(cfa_mi_table, file.path(cfa_dir, "cfa_modification_indices_diagnostic_only.csv"))
safe_csv(cfa_cr_ave_table, file.path(cfa_dir, "cfa_composite_reliability_ave.csv"))
safe_csv(cfa_alpha_table, file.path(cfa_dir, "cfa_alpha_bootstrap_ci.csv"))
safe_csv(cfa_htmt_table, file.path(cfa_dir, "cfa_htmt_pearson_bootstrap_ci.csv"))
safe_csv(cfa_htmt_ordinal_table, file.path(cfa_dir, "cfa_htmt_polychoric_sensitivity.csv"))
safe_csv(cfa_lrt_table, file.path(cfa_dir, "cfa_rival_model_robust_lrt.csv"))

# Fornell-Larcker y banderas de validez. Se reportan junto con HTMT; ninguna
# regla se utiliza de forma aislada para eliminar constructos o items.
fornell_larcker <- cfa_latent_cor_table %>%
  transmute(
    country,
    construct_1 = pmin(row, column),
    construct_2 = pmax(row, column),
    latent_correlation
  ) %>%
  distinct(country, construct_1, construct_2, .keep_all = TRUE) %>%
  left_join(
    cfa_cr_ave_table %>%
      transmute(country, construct_1 = construct, sqrt_ave_1 = sqrt(ave)),
    by = c("country", "construct_1")
  ) %>%
  left_join(
    cfa_cr_ave_table %>%
      transmute(country, construct_2 = construct, sqrt_ave_2 = sqrt(ave)),
    by = c("country", "construct_2")
  ) %>%
  mutate(
    abs_latent_correlation = abs(latent_correlation),
    fornell_larcker_passed = abs_latent_correlation < sqrt_ave_1 &
      abs_latent_correlation < sqrt_ave_2
  )
safe_csv(fornell_larcker, file.path(cfa_dir, "cfa_fornell_larcker.csv"))

validity_flags <- cfa_htmt_table %>%
  mutate(
    pair_1 = pmin(construct_1, construct_2),
    pair_2 = pmax(construct_1, construct_2),
    construct_1 = pair_1,
    construct_2 = pair_2
  ) %>%
  select(-pair_1, -pair_2) %>%
  mutate(
    htmt_point_below_085 = htmt < .85,
    htmt_point_below_090 = htmt < .90,
    htmt_ci_upper_below_085 = ci_high < .85,
    htmt_ci_upper_below_090 = ci_high < .90
  ) %>%
  left_join(
    fornell_larcker %>%
      select(country, construct_1, construct_2, fornell_larcker_passed),
    by = c("country", "construct_1", "construct_2")
  )
safe_csv(validity_flags, file.path(cfa_dir, "cfa_discriminant_validity_flags.csv"))

saveRDS(pooled_cfa_fits, file.path(cfa_dir, "pooled_cfa_models_mlr.rds"))
saveRDS(
  map(country_cfa_results, ~ map(.x, "value")),
  file.path(cfa_dir, "country_cfa_models_mlr.rds")
)
if (run_ordinal) {
  saveRDS(
    map(country_wlsmv_results, "value"),
    file.path(cfa_dir, "country_cfa_models_wlsmv_sensitivity.rds")
  )
} else {
  safe_text(
    "Sensibilidad ordinal omitida en esta corrida (SEM_RUN_ORDINAL=0).",
    file.path(cfa_dir, "ORDINAL_SENSITIVITY_SKIPPED.txt")
  )
}

message("Tablas CFA/HTMT exportadas; evaluando compuertas e invariancia.")
# Compuertas de medicion.
country_main_status <- cfa_status_table %>%
  filter(
    country %in% country_levels,
    model == "theoretical_3factor",
    estimator == "MLR"
  )
country_cfa_admissible <- nrow(country_main_status) == length(country_levels) &&
  all(country_main_status$admissible)
country_cfa_fit_acceptable <- country_cfa_admissible &&
  all(country_main_status$fit_acceptable)

pooled_main_status <- cfa_status_table %>%
  filter(country == "Pooled", model == "theoretical_3factor", estimator == "MLR")
pooled_cfa_admissible <- nrow(pooled_main_status) == 1L && pooled_main_status$admissible[[1L]]
pooled_cfa_fit_acceptable <- pooled_cfa_admissible &&
  pooled_main_status$fit_acceptable[[1L]]

country_convergent_validity_screen <- cfa_cr_ave_table %>%
  filter(country %in% country_levels) %>%
  summarise(
    passed = n() == length(country_levels) * length(construct_blocks) &&
      all(composite_reliability >= .70, na.rm = FALSE) &&
      all(ave >= .50, na.rm = FALSE)
  ) %>%
  pull(passed)
country_discriminant_htmt_screen <- validity_flags %>%
  filter(country %in% country_levels) %>%
  summarise(
    passed = n() == length(country_levels) * choose(length(construct_blocks), 2) &&
      all(htmt_ci_upper_below_090, na.rm = FALSE)
  ) %>%
  pull(passed)
country_fornell_larcker_screen <- fornell_larcker %>%
  filter(country %in% country_levels) %>%
  summarise(
    passed = n() == length(country_levels) * choose(length(construct_blocks), 2) &&
      all(fornell_larcker_passed, na.rm = FALSE)
  ) %>%
  pull(passed)

measurement_gate <- tibble(
  gate = c(
    "country_CFA_mathematically_admissible",
    "country_CFA_fit_heuristics_acceptable",
    "country_convergent_validity_CR_AVE_screen",
    "country_discriminant_validity_HTMT_CI90_screen",
    "country_discriminant_validity_Fornell_Larcker_screen",
    "pooled_CFA_mathematically_admissible",
    "pooled_CFA_fit_heuristics_acceptable",
    "override_requested"
  ),
  passed = c(
    country_cfa_admissible,
    country_cfa_fit_acceptable,
    country_convergent_validity_screen,
    country_discriminant_htmt_screen,
    country_fornell_larcker_screen,
    pooled_cfa_admissible,
    pooled_cfa_fit_acceptable,
    override_gate
  ),
  implication = c(
    "Debe cumplirse para intentar invarianza multigrupo.",
    "Si falla, invarianza/SEM se consideran diagnosticos, no confirmatorios.",
    "Bandera descriptiva: CR >= .70 y AVE >= .50 en todos los grupos/constructos.",
    "Bandera descriptiva: limite superior IC95% de HTMT < .90.",
    "Bandera descriptiva complementaria; no decide por si sola.",
    "Debe cumplirse para estimar SEM pooled.",
    "Si falla, SEM pooled se etiqueta como diagnostico.",
    "El override nunca convierte un modelo deficiente en evidencia confirmatoria."
  )
)
safe_csv(measurement_gate, file.path(results_dir, "gatekeeping_measurement.csv"))

# ------------------------------------------------------------------------------
# 5. Invarianza multigrupo, condicionada a CFA admisibles por pais
# ------------------------------------------------------------------------------

invariance_results <- list()
metric_invariance_supported <- FALSE

if (country_cfa_admissible || override_gate) {
  inv_configural <- capture_conditions(
    cfa(
      measurement_model_3f,
      data = dat,
      group = "Pais",
      estimator = "MLR",
      std.lv = TRUE,
      meanstructure = TRUE,
      missing = "fiml"
    )
  )
  inv_metric <- capture_conditions(
    cfa(
      measurement_model_3f,
      data = dat,
      group = "Pais",
      estimator = "MLR",
      std.lv = TRUE,
      meanstructure = TRUE,
      missing = "fiml",
      group.equal = "loadings"
    )
  )
  inv_scalar <- capture_conditions(
    cfa(
      measurement_model_3f,
      data = dat,
      group = "Pais",
      estimator = "MLR",
      std.lv = TRUE,
      meanstructure = TRUE,
      missing = "fiml",
      group.equal = c("loadings", "intercepts")
    )
  )

  invariance_results <- list(
    configural = inv_configural$value,
    metric = inv_metric$value,
    scalar = inv_scalar$value
  )

  inv_fit <- bind_rows(
    extract_fit_row(inv_configural$value, "Configural", "MLR", "Multigroup"),
    extract_fit_row(inv_metric$value, "Metric", "MLR", "Multigroup"),
    extract_fit_row(inv_scalar$value, "Scalar", "MLR", "Multigroup")
  )
  inv_status <- bind_rows(
    extract_admissibility(inv_configural$value, "Multigroup", "Configural", "MLR"),
    extract_admissibility(inv_metric$value, "Multigroup", "Metric", "MLR"),
    extract_admissibility(inv_scalar$value, "Multigroup", "Scalar", "MLR")
  )

  robust_or_plain <- function(tbl, row_id, robust_col, plain_col) {
    val <- tbl[[robust_col]][row_id]
    if (is.na(val)) tbl[[plain_col]][row_id] else val
  }
  cfi_cfg <- robust_or_plain(inv_fit, 1L, "cfi_robust", "cfi")
  cfi_met <- robust_or_plain(inv_fit, 2L, "cfi_robust", "cfi")
  cfi_sca <- robust_or_plain(inv_fit, 3L, "cfi_robust", "cfi")
  rmsea_cfg <- robust_or_plain(inv_fit, 1L, "rmsea_robust", "rmsea")
  rmsea_met <- robust_or_plain(inv_fit, 2L, "rmsea_robust", "rmsea")
  rmsea_sca <- robust_or_plain(inv_fit, 3L, "rmsea_robust", "rmsea")

  inv_deltas <- tibble(
    comparison = c("Metric - Configural", "Scalar - Metric"),
    delta_cfi = c(cfi_met - cfi_cfg, cfi_sca - cfi_met),
    delta_rmsea = c(rmsea_met - rmsea_cfg, rmsea_sca - rmsea_met),
    delta_srmr = c(
      inv_fit$srmr[2L] - inv_fit$srmr[1L],
      inv_fit$srmr[3L] - inv_fit$srmr[2L]
    ),
    cfi_rule_passed = delta_cfi >= -.010,
    rmsea_rule_passed = delta_rmsea <= .015,
    srmr_rule_passed = c(delta_srmr[1L] <= .030, delta_srmr[2L] <= .010)
  )

  lrt_cfg_metric <- safe_lrt(
    inv_metric$value, inv_configural$value,
    "Configural vs Metric", "Multigroup"
  )
  lrt_metric_scalar <- safe_lrt(
    inv_scalar$value, inv_metric$value,
    "Metric vs Scalar", "Multigroup"
  )

  inv_conditions <- bind_rows(
    tibble(
      model = "Configural", error = inv_configural$error,
      warnings = paste(inv_configural$warnings, collapse = " | ")
    ),
    tibble(
      model = "Metric", error = inv_metric$error,
      warnings = paste(inv_metric$warnings, collapse = " | ")
    ),
    tibble(
      model = "Scalar", error = inv_scalar$error,
      warnings = paste(inv_scalar$warnings, collapse = " | ")
    )
  )

  safe_csv(inv_fit, file.path(invariance_dir, "invariance_fit_complete.csv"))
  safe_csv(inv_status, file.path(invariance_dir, "invariance_admissibility.csv"))
  safe_csv(inv_deltas, file.path(invariance_dir, "invariance_delta_indices.csv"))
  safe_csv(
    bind_rows(lrt_cfg_metric, lrt_metric_scalar),
    file.path(invariance_dir, "invariance_robust_lrt.csv")
  )
  safe_csv(inv_conditions, file.path(invariance_dir, "invariance_errors_warnings.csv"))

  if (is_lavaan_fit(inv_configural$value)) {
    safe_csv(
      extract_parameters(inv_configural$value, "Multigroup", "Configural") %>%
        filter(op == "=~"),
      file.path(invariance_dir, "configural_loadings_se_ci.csv")
    )
    safe_csv(
      extract_latent_correlations(inv_configural$value, "Multigroup", "Configural"),
      file.path(invariance_dir, "configural_latent_correlations.csv")
    )
  }
  if (is_lavaan_fit(inv_metric$value)) {
    safe_csv(
      extract_parameters(inv_metric$value, "Multigroup", "Metric") %>%
        filter(op == "=~"),
      file.path(invariance_dir, "metric_loadings_se_ci.csv")
    )
    score_test <- capture_conditions(lavTestScore(inv_metric$value, epc = TRUE))
    if (!is.na(score_test$error)) {
      safe_text(
        c("No se pudo calcular lavTestScore:", score_test$error),
        file.path(invariance_dir, "metric_equality_score_test_error.txt")
      )
    } else {
      if (!is.null(score_test$value$test)) {
        safe_csv(
          as.data.frame(score_test$value$test),
          file.path(invariance_dir, "metric_equality_score_test_global.csv")
        )
      }
      if (!is.null(score_test$value$uni)) {
        safe_csv(
          as.data.frame(score_test$value$uni),
          file.path(invariance_dir, "metric_equality_score_test_univariate.csv")
        )
      }
    }
  }

  saveRDS(invariance_results, file.path(invariance_dir, "measurement_invariance_models.rds"))

  configural_admissible <- inv_status$admissible[inv_status$model == "Configural"]
  metric_admissible <- inv_status$admissible[inv_status$model == "Metric"]
  configural_fit_ok <- inv_status$fit_acceptable[inv_status$model == "Configural"]
  delta_metric_ok <- nrow(inv_deltas) >= 1L && all(
    unlist(inv_deltas[1L, c(
      "cfi_rule_passed", "rmsea_rule_passed", "srmr_rule_passed"
    )])
  )
  metric_invariance_supported <-
    isTRUE(configural_admissible) && isTRUE(metric_admissible) &&
    isTRUE(configural_fit_ok) && isTRUE(delta_metric_ok)
} else {
  safe_text(
    c(
      "INVARIANZA OMITIDA",
      "Al menos un CFA por pais no fue matematicamente admisible.",
      "No se estimaron modelos multigrupo de invarianza ni rutas comparativas."
    ),
    file.path(invariance_dir, "INVARIANCE_SKIPPED.txt")
  )
}

invariance_gate <- tibble(
  gate = c(
    "country_CFA_admissible",
    "metric_invariance_supported_including_configural_fit",
    "multigroup_structural_allowed",
    "override_requested"
  ),
  passed = c(
    country_cfa_admissible,
    metric_invariance_supported,
    metric_invariance_supported || override_gate,
    override_gate
  ),
  note = c(
    "Condicion minima para estimar invarianza.",
    "Exige admisibilidad, ajuste configural heuristico y deltas aceptables.",
    "Sin esta condicion, cualquier SEM multigrupo es solo diagnostico.",
    "El override no convierte resultados diagnosticos en confirmatorios."
  )
)
safe_csv(invariance_gate, file.path(results_dir, "gatekeeping_invariance.csv"))

# ------------------------------------------------------------------------------
# 6. Modelos estructurales con valores observados
# ------------------------------------------------------------------------------

country_codes <- c("arg", "col", "per", "mex")
names(country_codes) <- c("Argentina", "Colombia", "Peru", "Mexico")

label_modifier <- function(base, group_codes = NULL) {
  if (is.null(group_codes)) return(base)
  paste0("c(", paste0(base, "_", group_codes, collapse = ","), ")")
}

label_name <- function(base, group_code = NULL) {
  if (is.null(group_code)) base else paste0(base, "_", group_code)
}

make_observed_value_sem <- function(group_codes = NULL, controls = FALSE) {
  lm <- function(base) label_modifier(base, group_codes)
  control_terms <- if (controls) {
    " + Sexo_num + Estudios_num + Edad_z + NSE_z"
  } else {
    ""
  }

  regression_lines <- c(
    paste0("Benevolencia ~ ", lm("a_ben"), "*CSR1_bin", control_terms),
    paste0("Universalismo ~ ", lm("a_uni"), "*CSR1_bin", control_terms),
    paste0("Poder ~ ", lm("a_pow"), "*CSR1_bin", control_terms),
    paste0("Logro ~ ", lm("a_log"), "*CSR1_bin", control_terms),
    paste0(
      "SC ~ ", lm("b_ben"), "*Benevolencia + ",
      lm("b_uni"), "*Universalismo + ", lm("b_pow"), "*Poder + ",
      lm("b_log"), "*Logro", control_terms
    ),
    paste0(
      "SBL ~ ", lm("c_ben"), "*Benevolencia + ",
      lm("c_uni"), "*Universalismo + ", lm("c_pow"), "*Poder + ",
      lm("c_log"), "*Logro", control_terms
    ),
    paste0(
      "PSB ~ ", lm("d_sc"), "*SC + ", lm("d_sbl"), "*SBL + ",
      lm("p_ben"), "*Benevolencia + ", lm("p_uni"), "*Universalismo + ",
      lm("p_pow"), "*Poder + ", lm("p_log"), "*Logro", control_terms
    )
  )

  residual_covariances <- c(
    "Benevolencia ~~ Universalismo + Poder + Logro",
    "Universalismo ~~ Poder + Logro",
    "Poder ~~ Logro",
    "SC ~~ SBL"
  )

  values <- c(ben = "Benevolencia", uni = "Universalismo", pow = "Poder", log = "Logro")
  groups_for_defs <- if (is.null(group_codes)) NA_character_ else group_codes
  defined <- character()

  for (g in groups_for_defs) {
    suffix <- if (is.na(g)) NULL else g
    for (v in names(values)) {
      b <- label_name(paste0("b_", v), suffix)
      cc <- label_name(paste0("c_", v), suffix)
      p <- label_name(paste0("p_", v), suffix)
      dsc <- label_name("d_sc", suffix)
      dsbl <- label_name("d_sbl", suffix)
      a <- label_name(paste0("a_", v), suffix)
      out_suffix <- if (is.null(suffix)) "" else paste0("_", suffix)

      defined <- c(
        defined,
        paste0("ind_", v, "_SC", out_suffix, " := ", b, "*", dsc),
        paste0("ind_", v, "_SBL", out_suffix, " := ", cc, "*", dsbl),
        paste0(
          "ind_", v, "_total", out_suffix, " := ",
          b, "*", dsc, " + ", cc, "*", dsbl
        ),
        paste0(
          "ind_", v, "_SBL_minus_SC", out_suffix, " := ",
          cc, "*", dsbl, " - ", b, "*", dsc
        ),
        paste0(
          "total_", v, "_to_PSB", out_suffix, " := ", p, " + ",
          b, "*", dsc, " + ", cc, "*", dsbl
        ),
        paste0(
          "CSR_via_", v, "_to_PSB", out_suffix, " := ", a, "*(", p,
          " + ", b, "*", dsc, " + ", cc, "*", dsbl, ")"
        )
      )
    }

    out_suffix <- if (is.null(suffix)) "" else paste0("_", suffix)
    dsc <- label_name("d_sc", suffix)
    dsbl <- label_name("d_sbl", suffix)
    defined <- c(
      defined,
      paste0("direct_SBL_minus_SC", out_suffix, " := ", dsbl, " - ", dsc),
      paste0(
        "CSR_total_indirect_to_PSB", out_suffix, " := ",
        paste0(
          vapply(names(values), function(v) {
            a <- label_name(paste0("a_", v), suffix)
            p <- label_name(paste0("p_", v), suffix)
            b <- label_name(paste0("b_", v), suffix)
            cc <- label_name(paste0("c_", v), suffix)
            paste0(a, "*(", p, " + ", b, "*", dsc, " + ", cc, "*", dsbl, ")")
          }, character(1)),
          collapse = " + "
        )
      )
    )
  }

  # Contrastes planificados frente a Argentina. Se prueban diferencias de
  # coeficientes; no se infieren diferencias comparando significancias separadas.
  if (!is.null(group_codes) && length(group_codes) >= 2L) {
    comparison_bases <- c(
      "a_ben", "a_uni", "a_pow", "a_log",
      "b_ben", "b_uni", "b_pow", "b_log",
      "c_ben", "c_uni", "c_pow", "c_log",
      "d_sc", "d_sbl", "p_ben", "p_uni", "p_pow", "p_log"
    )
    ref <- group_codes[1L]
    for (other in group_codes[-1L]) {
      for (base in comparison_bases) {
        defined <- c(
          defined,
          paste0(
            "diff_", base, "_", other, "_vs_", ref, " := ",
            base, "_", other, " - ", base, "_", ref
          )
        )
      }
    }
  }

  paste(
    measurement_model_3f,
    paste(regression_lines, collapse = "\n"),
    paste(residual_covariances, collapse = "\n"),
    paste(defined, collapse = "\n"),
    sep = "\n\n"
  )
}

make_composite_value_sem <- function(group_codes = NULL) {
  lm <- function(base) label_modifier(base, group_codes)
  regressions <- c(
    paste0("ST_comp ~ ", lm("a_st"), "*CSR1_bin"),
    paste0("SE_comp ~ ", lm("a_se"), "*CSR1_bin"),
    paste0(
      "SC ~ ", lm("b_st"), "*ST_comp + ", lm("b_se"), "*SE_comp"
    ),
    paste0(
      "SBL ~ ", lm("c_st"), "*ST_comp + ", lm("c_se"), "*SE_comp"
    ),
    paste0(
      "PSB ~ ", lm("d_sc"), "*SC + ", lm("d_sbl"), "*SBL + ",
      lm("p_st"), "*ST_comp + ", lm("p_se"), "*SE_comp"
    ),
    "ST_comp ~~ SE_comp",
    "SC ~~ SBL"
  )

  groups_for_defs <- if (is.null(group_codes)) NA_character_ else group_codes
  defined <- character()
  for (g in groups_for_defs) {
    suffix <- if (is.na(g)) NULL else g
    os <- if (is.null(suffix)) "" else paste0("_", suffix)
    for (v in c("st", "se")) {
      b <- label_name(paste0("b_", v), suffix)
      cc <- label_name(paste0("c_", v), suffix)
      p <- label_name(paste0("p_", v), suffix)
      a <- label_name(paste0("a_", v), suffix)
      dsc <- label_name("d_sc", suffix)
      dsbl <- label_name("d_sbl", suffix)
      defined <- c(
        defined,
        paste0("ind_", v, "_SC", os, " := ", b, "*", dsc),
        paste0("ind_", v, "_SBL", os, " := ", cc, "*", dsbl),
        paste0("ind_", v, "_total", os, " := ", b, "*", dsc, " + ", cc, "*", dsbl),
        paste0("total_", v, "_to_PSB", os, " := ", p, " + ", b, "*", dsc,
               " + ", cc, "*", dsbl),
        paste0("CSR_via_", v, "_to_PSB", os, " := ", a, "*(", p, " + ",
               b, "*", dsc, " + ", cc, "*", dsbl, ")")
      )
    }
    defined <- c(
      defined,
      paste0("direct_SBL_minus_SC", os, " := ", label_name("d_sbl", suffix),
             " - ", label_name("d_sc", suffix)),
      paste0(
        "CSR_total_indirect_to_PSB", os, " := ",
        label_name("a_st", suffix), "*(", label_name("p_st", suffix), " + ",
        label_name("b_st", suffix), "*", label_name("d_sc", suffix), " + ",
        label_name("c_st", suffix), "*", label_name("d_sbl", suffix), ") + ",
        label_name("a_se", suffix), "*(", label_name("p_se", suffix), " + ",
        label_name("b_se", suffix), "*", label_name("d_sc", suffix), " + ",
        label_name("c_se", suffix), "*", label_name("d_sbl", suffix), ")"
      )
    )
  }

  if (!is.null(group_codes) && length(group_codes) >= 2L) {
    bases <- c(
      "a_st", "a_se", "b_st", "b_se", "c_st", "c_se",
      "d_sc", "d_sbl", "p_st", "p_se"
    )
    ref <- group_codes[1L]
    for (other in group_codes[-1L]) {
      for (base in bases) {
        defined <- c(
          defined,
          paste0("diff_", base, "_", other, "_vs_", ref, " := ",
                 base, "_", other, " - ", base, "_", ref)
        )
      }
    }
  }

  paste(
    measurement_model_3f,
    paste(regressions, collapse = "\n"),
    paste(defined, collapse = "\n"),
    sep = "\n\n"
  )
}

fit_sem_safe <- function(
    model, data, estimator = "MLR", group = NULL, group_equal = NULL,
    bootstrap = NULL) {
  capture_conditions({
    args <- list(
      model = model,
      data = data,
      estimator = estimator,
      std.lv = TRUE,
      meanstructure = TRUE,
      fixed.x = TRUE
    )
    if (!is.null(group)) args$group <- group
    if (!is.null(group_equal)) args$group.equal <- group_equal
    if (identical(estimator, "MLR")) {
      args$missing <- "listwise"
    }
    if (identical(estimator, "ML") && !is.null(bootstrap) && bootstrap > 0L) {
      args$se <- "bootstrap"
      args$bootstrap <- bootstrap
      args$test <- "standard"
      args$missing <- "listwise"
    }
    do.call(lavaan::sem, args)
  })
}

export_sem_bundle <- function(result, prefix, interpretation_status) {
  fit <- result$value
  current_status <- extract_admissibility(fit, prefix, prefix, "MLR")
  exploratory <- grepl("EXPLORATORY", interpretation_status)
  measurement_passed <- grepl("MEASUREMENT_GATE_PASSED", interpretation_status)
  own_admissible <- isTRUE(current_status$admissible[[1L]])
  own_fit_ok <- isTRUE(current_status$fit_acceptable[[1L]])
  interpretation_status <- if (!own_admissible) {
    "UNINTERPRETABLE_SEM_NOT_ADMISSIBLE"
  } else if (exploratory) {
    if (own_fit_ok) "EXPLORATORY_COMPOSITE_SEM_FIT_PASSED" else "DIAGNOSTIC_EXPLORATORY_COMPOSITE_SEM_FIT_FAILED"
  } else if (measurement_passed && own_fit_ok) {
    "ASSOCIATIVE_MODEL_MEASUREMENT_AND_SEM_FIT_PASSED"
  } else {
    "DIAGNOSTIC_ONLY_MEASUREMENT_OR_SEM_FIT_FAILED"
  }
  fit_row <- extract_fit_row(fit, prefix, "MLR", prefix) %>%
    mutate(interpretation_status = interpretation_status)
  status <- extract_admissibility(fit, prefix, prefix, "MLR") %>%
    mutate(interpretation_status = interpretation_status)
  safe_csv(fit_row, file.path(structural_dir, paste0(prefix, "_fit.csv")))
  safe_csv(status, file.path(structural_dir, paste0(prefix, "_admissibility.csv")))
  safe_csv(
    tibble(
      model = prefix,
      interpretation_status = interpretation_status,
      error = result$error,
      warnings = paste(result$warnings, collapse = " | ")
    ),
    file.path(structural_dir, paste0(prefix, "_errors_warnings.csv"))
  )
  if (!is_lavaan_fit(fit)) return(invisible(NULL))

  params <- extract_parameters(fit, prefix, prefix)
  safe_csv(params, file.path(structural_dir, paste0(prefix, "_all_parameters_se_ci.csv")))
  safe_csv(
    params %>% filter(op == "~"),
    file.path(structural_dir, paste0(prefix, "_structural_paths_se_ci.csv"))
  )
  safe_csv(
    parameterEstimates(fit, standardized = TRUE, ci = TRUE, level = .95) %>%
      as_tibble() %>%
      filter(op == ":=") %>%
      mutate(
        multiplicity_family = case_when(
          grepl("^diff_", lhs) ~ "planned_cross_country_contrasts",
          grepl("SBL_minus_SC", lhs) ~ "planned_mediator_contrasts",
          TRUE ~ "indirect_and_total_effects"
        )
      ) %>%
      group_by(multiplicity_family) %>%
      mutate(p_holm_within_family = p.adjust(pvalue, method = "holm")) %>%
      ungroup(),
    file.path(structural_dir, paste0(prefix, "_indirects_contrasts_ci.csv"))
  )
  safe_csv(
    extract_r2(fit, prefix, prefix),
    file.path(structural_dir, paste0(prefix, "_r2.csv"))
  )
  safe_csv(
    extract_latent_correlations(fit, prefix, prefix),
    file.path(structural_dir, paste0(prefix, "_latent_correlations.csv"))
  )
  safe_csv(
    extract_modification_indices(fit, prefix, prefix, n_top = 100L),
    file.path(structural_dir, paste0(prefix, "_modification_indices_diagnostic_only.csv"))
  )
  saveRDS(fit, file.path(structural_dir, paste0(prefix, ".rds")))
  invisible(NULL)
}

# Construccion explicita de compuestos para la sensibilidad. No se sustituyen las
# variables originales ni se escriben en el .sav.
dat <- dat %>%
  mutate(
    ST_comp = rowMeans(across(all_of(c("Benevolencia", "Universalismo"))),
                       na.rm = FALSE),
    SE_comp = rowMeans(across(all_of(c("Poder", "Logro"))), na.rm = FALSE)
  )

# Se calculan directamente las fiabilidades de los pares para evitar presentarlos
# como si fueran las escalas SC/SBL/PSB.
value_pair_reliability <- bind_rows(
  tibble(
    country = "Pooled", composite = "ST_comp",
    alpha = cronbach_alpha(dat[, c("Benevolencia", "Universalismo")])
  ),
  tibble(
    country = "Pooled", composite = "SE_comp",
    alpha = cronbach_alpha(dat[, c("Poder", "Logro")])
  ),
  map_dfr(country_levels, function(country) {
    dc <- dat[dat$Pais == country, , drop = FALSE]
    tibble(
      country = country,
      composite = c("ST_comp", "SE_comp"),
      alpha = c(
        cronbach_alpha(dc[, c("Benevolencia", "Universalismo")]),
        cronbach_alpha(dc[, c("Poder", "Logro")])
      )
    )
  })
)
safe_csv(
  value_pair_reliability,
  file.path(structural_dir, "value_composite_reliability_warning.csv")
)
safe_text(
  c(
    "ADVERTENCIA SOBRE LOS COMPUESTOS DE VALORES",
    "ST_comp = promedio de Benevolencia y Universalismo.",
    "SE_comp = promedio de Poder y Logro.",
    "Estas sensibilidades no corrigen el error de medicion.",
    "La baja fiabilidad de ST, especialmente por pais, impide usar el compuesto",
    "como sustituto confirmatorio de un factor latente validado.",
    "Quijotismo no se incorpora automaticamente: su inclusion debe decidirse con",
    "el manual y la teoria del instrumento, no por mejora oportunista del ajuste."
  ),
  file.path(structural_dir, "IMPORTANT_composite_warning.txt")
)

# Sintaxis reproducible de los modelos estructurales.
pooled_primary_syntax <- make_observed_value_sem()
pooled_controlled_syntax <- make_observed_value_sem(controls = TRUE)
pooled_composite_syntax <- make_composite_value_sem()
multigroup_primary_syntax <- make_observed_value_sem(group_codes = unname(country_codes))
multigroup_composite_syntax <- make_composite_value_sem(
  group_codes = unname(country_codes)
)

safe_text(pooled_primary_syntax,
          file.path(structural_dir, "syntax_pooled_observed_values.txt"))
safe_text(pooled_controlled_syntax,
          file.path(structural_dir, "syntax_pooled_observed_values_controls.txt"))
safe_text(pooled_composite_syntax,
          file.path(structural_dir, "syntax_pooled_composite_sensitivity.txt"))
safe_text(multigroup_primary_syntax,
          file.path(structural_dir, "syntax_multigroup_observed_values.txt"))
safe_text(multigroup_composite_syntax,
          file.path(structural_dir, "syntax_multigroup_composite_sensitivity.txt"))

pooled_interpretation_status <- if (pooled_cfa_fit_acceptable) {
  "MEASUREMENT_GATE_PASSED"
} else {
  "DIAGNOSTIC_ONLY_MEASUREMENT_FIT_GATE_FAILED"
}

pooled_sem_primary <- NULL
pooled_sem_controlled <- NULL
pooled_sem_composite <- NULL
pooled_sem_bootstrap <- NULL

if (pooled_cfa_admissible || override_gate) {
  pooled_sem_primary <- fit_sem_safe(
    pooled_primary_syntax, dat, estimator = "MLR"
  )
  export_sem_bundle(
    pooled_sem_primary,
    "pooled_observed_values_MLR",
    pooled_interpretation_status
  )

  # Sensibilidad con controles demograficos. No se interpreta como ajuste causal
  # definitivo: permite verificar estabilidad de signos y magnitudes.
  pooled_sem_controlled <- fit_sem_safe(
    pooled_controlled_syntax, dat, estimator = "MLR"
  )
  export_sem_bundle(
    pooled_sem_controlled,
    "pooled_observed_values_controls_MLR",
    paste0(pooled_interpretation_status, "_COVARIATE_SENSITIVITY")
  )

  # Sensibilidad exploratoria con compuestos de valores.
  pooled_sem_composite <- fit_sem_safe(
    pooled_composite_syntax, dat, estimator = "MLR"
  )
  export_sem_bundle(
    pooled_sem_composite,
    "pooled_value_composites_MLR_sensitivity",
    "EXPLORATORY_LOW_RELIABILITY_COMPOSITE_SENSITIVITY"
  )

  # IC percentiles de efectos definidos bajo ML. El modelo principal sigue MLR.
  # Guarda todas las replicas, incluso las invalidas, para auditar qué entra al IC.
  bootstrap_status <- "SKIPPED_BY_SETTING"
  if (boot_indirect > 0L && is_lavaan_fit(pooled_sem_primary$value)) {
    message("Bootstrap pooled indirectos ML: ", boot_indirect, " replicas, ", ncpus, " CPU(s)")
    ml_base <- capture_conditions(lavaan::sem(pooled_primary_syntax, data = dat,
      estimator = "ML", std.lv = TRUE, meanstructure = TRUE, fixed.x = TRUE,
      missing = "listwise", se = "none", test = "none", baseline = FALSE, h1 = FALSE))
    if (is_lavaan_fit(ml_base$value) && isTRUE(lavInspect(ml_base$value, "converged"))) {
      boot_fun <- function(fit) {
        pe <- lavaan::parameterEstimates(fit, se = FALSE, ci = FALSE, standardized = FALSE)
        defs <- pe[pe$op == ":=", , drop = FALSE]
        post_ok <- isTRUE(suppressWarnings(lavaan::lavInspect(fit, "post.check")))
        neg <- any(pe$op == "~~" & pe$lhs == pe$rhs & pe$est < -1e-6, na.rm = TRUE)
        cl <- lavaan::lavInspect(fit, "cor.lv")
        cv <- lavaan::lavInspect(fit, "cov.lv")
        corr_ok <- all(abs(cl[row(cl) != col(cl)]) < 1)
        eig_ok <- min(eigen(cv, symmetric = TRUE, only.values = TRUE)$values) > 1e-8
        c(converged = 1, admissible = as.numeric(post_ok && !neg && corr_ok && eig_ok),
          setNames(defs$est, defs$lhs))
      }
      bootstrap_result <- capture_conditions(lavaan::lavBootstrap(ml_base$value,
        r = boot_indirect, type = "ordinary", fun = boot_fun,
        parallel = if (ncpus > 1L) "snow" else "no",
        ncpus = ncpus, iseed = 20261004L))
      if (is.matrix(bootstrap_result$value)) {
        boot_mat <- bootstrap_result$value
        failures <- attr(boot_mat, "error.idx")
        nonadmissible <- attr(boot_mat, "nonadmissible")
        finite_rows <- apply(boot_mat[, -(1:2), drop = FALSE], 1L, function(x) all(is.finite(x)))
        valid <- is.finite(boot_mat[, "converged"]) & boot_mat[, "converged"] == 1 &
          is.finite(boot_mat[, "admissible"]) & boot_mat[, "admissible"] == 1 & finite_rows
        if (length(c(failures, nonadmissible))) valid[unique(c(failures, nonadmissible))] <- FALSE
        valid_n_boot <- sum(valid)
        invalid_n_boot <- sum(!valid)
        valid_fraction_boot <- mean(valid)
        sufficient <- valid_fraction_boot >= .95
        bootstrap_status <- if (sufficient) "COMPLETED_SUFFICIENT_VALID_REPLICATES" else "INSUFFICIENT_VALID_REPLICATES_CI_WITHHELD"
        diag <- tibble(requested = boot_indirect, returned = nrow(boot_mat),
          failed_or_nonconverged = length(failures), lavaan_nonadmissible = length(nonadmissible),
          valid = valid_n_boot, invalid = invalid_n_boot, valid_fraction = valid_fraction_boot,
          minimum_fraction_required = .95, sufficient = sufficient, seed = 20261004L,
          resampling = "ordinary individual-level iid bootstrap", ncpus = ncpus,
          warnings = paste(bootstrap_result$warnings, collapse = " | "))
        safe_csv(diag, file.path(structural_dir, "pooled_bootstrap_diagnostics.csv"))
        safe_csv(as_tibble(boot_mat) %>% mutate(rep = seq_len(nrow(boot_mat)), valid = valid, .before = 1),
                 file.path(structural_dir, "pooled_bootstrap_replicates.csv"))
        saveRDS(boot_mat, file.path(structural_dir, "pooled_bootstrap_replicates.rds"))
        point <- boot_fun(ml_base$value)
        effect_names <- names(point)[-(1:2)]
        boot_pe <- map_dfr(effect_names, function(name) {
          vals <- boot_mat[valid, name]
          ci <- if (sufficient && length(vals) > 1L) quantile(vals, c(.025, .975), names = FALSE) else c(NA_real_, NA_real_)
          tibble(lhs = name, op = ":=", estimator = "ML", est = unname(point[name]),
                 bootstrap_se = if (length(vals) > 1L) sd(vals) else NA_real_,
                 ci.lower = ci[1], ci.upper = ci[2], requested_replicates = boot_indirect,
                 valid_replicates = sum(valid), valid_fraction = mean(valid),
                 ci_method = "95% ordinary bootstrap percentile; unstandardized effect",
                 status = bootstrap_status,
                 zero_excluded = if (sufficient) ci[1] > 0 || ci[2] < 0 else NA,
                 interpretation = "Associative sensitivity; no temporal or causal identification; iid does not address university clustering")
        })
        safe_csv(boot_pe, file.path(structural_dir, "pooled_observed_values_bootstrap_indirects_ci.csv"))
        saveRDS(ml_base$value, file.path(structural_dir, "pooled_observed_values_ML_bootstrap_base.rds"))
      } else {
        bootstrap_status <- "FAILED"
        safe_text(c(bootstrap_result$error, bootstrap_result$warnings),
                  file.path(structural_dir, "pooled_bootstrap_error.txt"))
      }
    } else {
      bootstrap_status <- "ML_BASE_NOT_CONVERGED"
      safe_text(c(ml_base$error, ml_base$warnings), file.path(structural_dir, "pooled_bootstrap_error.txt"))
    }
  }
  safe_csv(tibble(stage = "indirect_bootstrap", status = bootstrap_status,
                  requested_replicates = boot_indirect), file.path(structural_dir, "bootstrap_stage_status.csv"))

} else {
  safe_text(
    c(
      "SEM POOLED OMITIDO",
      "El CFA pooled no fue matematicamente admisible.",
      "No se deben interpretar rutas hasta corregir el modelo de medicion."
    ),
    file.path(structural_dir, "POOLED_SEM_SKIPPED.txt")
  )
}

multigroup_sem_primary <- NULL
multigroup_sem_composite <- NULL
multigroup_interpretation_status <- if (metric_invariance_supported) {
  "MEASUREMENT_GATE_PASSED_METRIC"
} else {
  "DIAGNOSTIC_ONLY_METRIC_GATE_FAILED_OVERRIDE"
}

if (metric_invariance_supported || override_gate) {
  multigroup_sem_primary <- fit_sem_safe(
    multigroup_primary_syntax,
    dat,
    estimator = "MLR",
    group = "Pais",
    group_equal = "loadings"
  )
  export_sem_bundle(
    multigroup_sem_primary,
    "multigroup_observed_values_metric_MLR",
    multigroup_interpretation_status
  )

  multigroup_sem_composite <- fit_sem_safe(
    multigroup_composite_syntax,
    dat,
    estimator = "MLR",
    group = "Pais",
    group_equal = "loadings"
  )
  export_sem_bundle(
    multigroup_sem_composite,
    "multigroup_value_composites_metric_MLR_sensitivity",
    "EXPLORATORY_LOW_RELIABILITY_COMPOSITE_SENSITIVITY"
  )
} else {
  safe_text(
    c(
      "SEM MULTIGRUPO OMITIDO",
      "No se obtuvo soporte suficiente para invarianza metrica bajo las compuertas",
      "predefinidas. No se compararon rutas entre paises.",
      "Para ejecutar el modelo exclusivamente como diagnostico, defina",
      "SEM_OVERRIDE_GATE=1. Ese override no legitima inferencias comparativas."
    ),
    file.path(structural_dir, "MULTIGROUP_SEM_SKIPPED.txt")
  )
}

saveRDS(
  list(
    pooled_primary = if (!is.null(pooled_sem_primary)) pooled_sem_primary$value else NULL,
    pooled_controls = if (!is.null(pooled_sem_controlled)) pooled_sem_controlled$value else NULL,
    pooled_composite = if (!is.null(pooled_sem_composite)) pooled_sem_composite$value else NULL,
    multigroup_primary = if (!is.null(multigroup_sem_primary)) multigroup_sem_primary$value else NULL,
    multigroup_composite = if (!is.null(multigroup_sem_composite)) multigroup_sem_composite$value else NULL
  ),
  file.path(structural_dir, "all_structural_models.rds")
)

# ------------------------------------------------------------------------------
# 7. Registro final, criterios y guia de interpretacion
# ------------------------------------------------------------------------------

analysis_settings <- tibble(
  setting = c(
    "run_id", "script_path", "script_md5", "sav_path", "sav_md5", "analysis_timestamp", "R_version", "primary_estimator",
    "ordinal_sensitivity_estimator", "bootstrap_htmt_alpha",
    "bootstrap_indirects", "run_ordinal_sensitivity", "run_rival_models",
    "override_gate", "country_order",
    "measurement_factors", "value_model_primary", "value_model_sensitivity"
  ),
  value = c(
    run_id, script_path, unname(tools::md5sum(script_path)),
    sav_path,
    unname(tools::md5sum(sav_path)),
    format(Sys.time(), tz = "UTC", usetz = TRUE),
    R.version.string,
    "MLR",
    "WLSMV by country; no automatic category collapsing",
    as.character(boot_htmt),
    as.character(boot_indirect),
    as.character(run_ordinal),
    as.character(run_rivals),
    as.character(override_gate),
    paste(country_levels, collapse = " | "),
    "SC, SBL, PSB",
    "Benevolencia, Universalismo, Poder, Logro as observed variables",
    "ST_comp and SE_comp; exploratory because reliability is limited"
  )
)
safe_csv(analysis_settings, file.path(results_dir, "analysis_settings.csv"))

gate_summary <- bind_rows(
  measurement_gate %>% transmute(stage = "Measurement", gate, passed, note = implication),
  invariance_gate %>% transmute(stage = "Invariance", gate, passed, note)
)
safe_csv(gate_summary, file.path(results_dir, "gatekeeping_summary.csv"))

session_lines <- capture.output(sessionInfo())
safe_text(session_lines, file.path(results_dir, "sessionInfo.txt"))

package_versions <- tibble(
  package = required_packages,
  version = vapply(
    required_packages,
    function(pkg) as.character(utils::packageVersion(pkg)),
    character(1)
  )
)
safe_csv(package_versions, file.path(results_dir, "package_versions.csv"))

readme_lines <- c(
  "RESULTADOS DEL REANALISIS CFA / SEM",
  "===================================",
  "",
  "Orden de lectura recomendado:",
  "1. 01_audit/range_and_missing_audit.csv y quality_summary.csv.",
  "2. 02_cfa/cfa_admissibility_and_fit_gate.csv.",
  "3. 02_cfa/cfa_fit_complete.csv, correlaciones latentes, AVE y HTMT.",
  "4. 02_cfa/cfa_rival_model_robust_lrt.csv y diagnosticos locales.",
  "5. gatekeeping_summary.csv.",
  "6. Solo si las compuertas lo permiten: 03_invariance y 04_structural.",
  "",
  "Reglas de gobierno analitico:",
  "- Admisibilidad significa convergencia, post-check correcto, ausencia de",
  "  varianzas residuales negativas y matriz latente admisible.",
  "- Los umbrales CFI >= .90, RMSEA <= .08 y SRMR <= .08 se usan como",
  "  heuristicas transparentes, no como reglas mecanicas universales.",
  "- La invarianza metrica requiere, ademas, deterioro limitado en CFI, RMSEA",
  "  y SRMR. El chi-cuadrado robusto se informa, pero no decide por si solo.",
  "- Los indices de modificacion no autorizan cambios automaticos. Cualquier",
  "  error correlacionado, carga cruzada o eliminacion de item exige teoria y",
  "  replicacion entre paises.",
  "- WLSMV es una sensibilidad. La escasez y ausencia de categorias inferiores",
  "  por pais puede volver inestables los umbrales; el script no fusiona",
  "  categorias de forma oportunista.",
  "- Los contrastes entre paises se estiman directamente. No se concluyen",
  "  diferencias porque un p sea significativo en un pais y no en otro.",
  "- Los valores se modelan como items observados. Sus contrastes entre paises",
  "  presuponen comparabilidad de respuesta y deben interpretarse con cautela,",
  "  pues no existe aqui una escala latente de valores con invarianza demostrada.",
  "- Los coeficientes no estandarizados son la base de los contrastes formales;",
  "  los estandarizados se reportan como magnitudes descriptivas.",
  "- Los indirectos MLR tienen IC delta robustos. La disponibilidad del bootstrap",
  paste0("  ML depende de la etapa: replicas solicitadas = ", boot_indirect, "; revise bootstrap_stage_status.csv y diagnosticos."),
  paste0("- HTMT/alfa bootstrap: ", boot_htmt, " replicas solicitadas; ver valid_replicates y valid_fraction."),
  paste0("- Sensibilidad WLSMV solicitada: ", run_ordinal, "; modelos rivales solicitados: ", run_rivals, "."),
  "- Los IC bootstrap son percentiles, individuales iid, condicionados a >=95%",
  "  de replicas convergentes y matematicamente admisibles; no corrigen agrupamiento universitario.",
  "- Con 1000 replicas, cada cola del IC95% usa aproximadamente 25 replicas;",
  "  los extremos tienen error Monte Carlo. Efectos cercanos a cero requieren mas replicas.",
  "- Las etiquetas de ajuste aprobadas no establecen causalidad ni validacion externa.",
  "- El modelo curricular principal impone ausencia de rutas directas CSR1->SC/SBL/PSB;",
  "  se requieren sensibilidades de mediacion parcial y ajuste por pais/institucion.",
  "- Los compuestos ST/SE son exploratorios y no resuelven error de medicion.",
  "- Un override ejecuta modelos diagnosticos, pero nunca convierte un mal",
  "  modelo de medicion en evidencia confirmatoria.",
  "",
  paste0("Estado CFA por pais admisible: ", country_cfa_admissible),
  paste0("Estado ajuste CFA por pais heuristico: ", country_cfa_fit_acceptable),
  paste0("Estado invarianza metrica: ", metric_invariance_supported),
  paste0("Override solicitado: ", override_gate),
  "",
  "La base original no fue modificada."
)
safe_text(readme_lines, file.path(results_dir, "README_results.txt"))

safe_text(c("LRT de modelos rivales: diagnostico exploratorio.",
  "Colapsar factores correlacionados puede imponer una correlacion de frontera;",
  "no interpretar el p de chi-cuadrado como prueba regular unica. Comparar AIC/BIC, ajuste y teoria.",
  "Un factor general con peor ajuste no descarta por si solo sesgo de metodo comun."),
  file.path(cfa_dir, "rival_model_interpretation_note.txt"))
write_run_status("COMPLETED", "Etapas opcionales y fallos de modelos documentados; COMPLETED significa flujo concluido, no modelo cientificamente validado.")
write_manifest()
cat("\nReanalisis terminado. Resultados escritos en:\n", results_dir, "\n", sep = "")
cat("Revise primero gatekeeping_summary.csv y README_results.txt.\n")
