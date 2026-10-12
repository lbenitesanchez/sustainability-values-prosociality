# Postproceso reproducible: indices ordinales escalados, etapas reales,
# sensibilidad Monte Carlo de los extremos bootstrap y manifiesto actualizado.
suppressPackageStartupMessages({library(lavaan); library(readr); library(dplyr); library(tibble); library(purrr)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Uso: Rscript finalize_canonical.R /ruta/outputs/run_id")
results_dir <- normalizePath(args[[1]], mustWork = TRUE)
status <- read_csv(file.path(results_dir, "run_status.csv"), show_col_types = FALSE)
if (!identical(status$state[[1]], "COMPLETED")) stop("Postprocesar solo un flujo concluido.")
write_csv_local <- function(x, name) write_csv(as.data.frame(x), file.path(results_dir, name), na = "")
models_path <- file.path(results_dir, "02_cfa", "country_cfa_models_wlsmv_sensitivity.rds")
if (file.exists(models_path)) {
  models <- readRDS(models_path)
  wanted <- c("chisq.scaled", "df.scaled", "pvalue.scaled", "cfi.scaled", "tli.scaled", "rmsea.scaled",
              "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "srmr")
  rows <- imap_dfr(models, function(fit, country) {
    if (!inherits(fit, "lavaan")) return(tibble(country = country, status = "FIT_OBJECT_UNAVAILABLE"))
    fm <- tryCatch(fitMeasures(fit, fit.measures = wanted), error = function(e) NULL)
    if (is.null(fm)) return(tibble(country = country, status = "FIT_MEASURES_UNAVAILABLE"))
    out <- as_tibble_row(as.list(fm))
    mutate(out, country = country, status = "SCALED_INDICES", .before = 1)
  })
  write_csv_local(rows, "02_cfa/ordinal_scaled_fit_supplement.csv")
  writeLines(c("WLSMV: separar indices convencionales, escalados y robustos.",
     "El chi-cuadrado WLSMV es escalado/desplazado; el CFI convencional alto no reemplaza el CFI escalado/robusto.",
     "ordinal_scaled_fit_supplement.csv conserva indices escalados; cfa_fit_complete.csv conserva convencionales/robustos.",
     "No seleccionar la variante del indice por ajuste favorable.",
     "La sensibilidad ordinal es CFA por país; no se afirma un SEM ordinal multigrupo válido.",
     "Revisar escasez/ausencia de categorias, vcov y estimaciones locales junto con indices globales."),
     file.path(results_dir, "02_cfa", "ordinal_fit_interpretation.txt"))
  dat <- readRDS(file.path(results_dir, "01_audit", "analytic_data.rds"))
  blocks <- c("APS1", "APS2", "APS3", "APS4", "APS5", "APS6", "APS7", "APS9", "APS10", "APS14",
              "APS11", "APS12", "APS13", "APS15", "APS16", "APS18", "APS19", "APS20", "APS21", "PP1", "PP2", "PP3")
  sparsity <- map_dfr(levels(dat$Pais), function(country) map_dfr(blocks, function(item) {
    z <- dat[[item]][dat$Pais == country]
    tab <- table(factor(z, levels = 1:5))
    tibble(country = country, item = item, n = length(z), observed_categories = sum(tab > 0),
           absent_categories = sum(tab == 0), categories_n_below_5 = sum(tab < 5),
           minimum_count_including_zero = min(tab))
  }))
  write_csv_local(sparsity, "02_cfa/ordinal_category_sparsity.csv")
}
boot_path <- file.path(results_dir, "04_structural", "pooled_bootstrap_replicates.csv")
if (file.exists(boot_path)) {
  reps <- read_csv(boot_path, show_col_types = FALSE)
  # Reconstruye exclusivamente contadores metadata desde las mascaras guardadas.
  # La primera corrida completa tuvo sombreado secuencial dentro de tibble();
  # nunca se recalculan parametros, replicas o IC durante esta correccion.
  diag_path <- file.path(results_dir, "04_structural", "pooled_bootstrap_diagnostics.csv")
  if (file.exists(diag_path)) {
    boot_diag <- read_csv(diag_path, show_col_types = FALSE)
    valid_n_boot <- sum(reps$valid)
    invalid_n_boot <- sum(!reps$valid)
    valid_fraction_boot <- mean(reps$valid)
    metadata_changed <- !isTRUE(all.equal(as.numeric(boot_diag$valid[[1]]), as.numeric(valid_n_boot))) ||
      !isTRUE(all.equal(as.numeric(boot_diag$invalid[[1]]), as.numeric(invalid_n_boot))) ||
      !isTRUE(all.equal(as.numeric(boot_diag$valid_fraction[[1]]), as.numeric(valid_fraction_boot)))
    if (metadata_changed) {
      raw_path <- file.path(results_dir, "04_structural", "pooled_bootstrap_diagnostics_raw.csv")
      if (!file.exists(raw_path)) invisible(file.copy(diag_path, raw_path))
      boot_diag$valid <- valid_n_boot
      boot_diag$invalid <- invalid_n_boot
      boot_diag$valid_fraction <- valid_fraction_boot
      write_csv_local(boot_diag, "04_structural/pooled_bootstrap_diagnostics.csv")
      write_csv_local(tibble(stage = "postprocessing_metadata_only", affected_file = "04_structural/pooled_bootstrap_diagnostics.csv",
        reason = "Sequential tibble data masking of column valid; rebuilt valid/invalid/fraction from saved logical replicate flags",
        parameters_reestimated = FALSE, bootstrap_reestimated = FALSE, intervals_changed = FALSE,
        valid = valid_n_boot, invalid = invalid_n_boot, valid_fraction = valid_fraction_boot), "patch_metadata.csv")
    }
  }
  effects <- setdiff(names(reps), c("rep", "valid", "converged", "admissible"))
  half <- floor(nrow(reps) / 2)
  mc <- map_dfr(effects, function(name) {
    all_values <- reps[[name]][reps$valid]
    first <- reps[[name]][reps$valid & reps$rep <= half]
    second <- reps[[name]][reps$valid & reps$rep > half]
    qa <- quantile(all_values, c(.025,.975), names = FALSE)
    q1 <- quantile(first, c(.025,.975), names = FALSE)
    q2 <- quantile(second, c(.025,.975), names = FALSE)
    tibble(effect = name, valid_replicates = length(all_values), first_half_valid = length(first), second_half_valid = length(second),
       pooled_low = qa[1], pooled_high = qa[2], first_half_low = q1[1], first_half_high = q1[2],
       second_half_low = q2[1], second_half_high = q2[2],
       maximum_endpoint_half_difference = max(abs(q1-q2)),
       endpoint_difference_fraction_of_interval_width = max(abs(q1-q2))/diff(qa),
       first_half_zero_excluded = q1[1] > 0 || q1[2] < 0,
       second_half_zero_excluded = q2[1] > 0 || q2[2] < 0,
       note = "Descriptive Monte Carlo stability diagnostic, not an inferential test or a rule for choosing favorable intervals")
  })
  write_csv_local(mc, "04_structural/bootstrap_montecarlo_stability.csv")
}
cfa <- read_csv(file.path(results_dir, "02_cfa", "cfa_admissibility_and_fit_gate.csv"), show_col_types = FALSE)
htmt <- read_csv(file.path(results_dir, "02_cfa", "cfa_htmt_pearson_bootstrap_ci.csv"), show_col_types = FALSE)
settings <- read_csv(file.path(results_dir, "analysis_settings.csv"), show_col_types = FALSE)
get_setting <- function(n) settings$value[match(n, settings$setting)]
stage_rows <- bind_rows(
  cfa %>% group_by(estimator) %>% summarise(stage = paste0("CFA_", first(estimator)), attempted_models = n(),
     converged_models = sum(converged), admissible_models = sum(admissible), fit_passed_models = sum(fit_acceptable),
     stage_status = "MODELS_ASSESSED_SEE_MODEL_LEVEL_RESULTS", .groups = "drop") %>% select(-estimator),
  tibble(stage = "HTMT_BOOTSTRAP", attempted_models = nrow(htmt), converged_models = NA_integer_,
     admissible_models = NA_integer_, fit_passed_models = NA_integer_,
     stage_status = if (as.integer(get_setting("bootstrap_htmt_alpha")) > 0L) "COMPLETED_SEE_VALID_REPLICATE_COUNTS" else "SKIPPED_BY_SETTING"),
  tibble(stage = "MULTIGROUP_SEM", attempted_models = if (file.exists(file.path(results_dir, "04_structural", "MULTIGROUP_SEM_SKIPPED.txt"))) 0L else 2L, converged_models = NA_integer_,
     admissible_models = NA_integer_, fit_passed_models = NA_integer_,
     stage_status = if (file.exists(file.path(results_dir, "04_structural", "MULTIGROUP_SEM_SKIPPED.txt"))) "SKIPPED_METRIC_MEASUREMENT_GATE_FAILED" else "SEE_MULTIGROUP_MODEL_FILES")
)
for (prefix in c("pooled_observed_values_MLR", "pooled_observed_values_controls_MLR", "pooled_value_composites_MLR_sensitivity")) {
  fit_file <- file.path(results_dir, "04_structural", paste0(prefix, "_admissibility.csv"))
  if (file.exists(fit_file)) {
    x <- read_csv(fit_file, show_col_types = FALSE)
    stage_rows <- bind_rows(stage_rows, tibble(stage = paste0("SEM_", prefix), attempted_models = 1L,
      converged_models = as.integer(x$converged[[1]]), admissible_models = as.integer(x$admissible[[1]]),
      fit_passed_models = as.integer(x$fit_acceptable[[1]]), stage_status = x$interpretation_status[[1]]))
  }
}
inv_file <- file.path(results_dir, "03_invariance", "invariance_admissibility.csv")
if (file.exists(inv_file)) {
  inv <- read_csv(inv_file, show_col_types = FALSE)
  stage_rows <- bind_rows(stage_rows, inv %>% transmute(stage = paste0("INVARIANCE_", model),
    attempted_models = 1L, converged_models = as.integer(converged), admissible_models = as.integer(admissible),
    fit_passed_models = as.integer(fit_acceptable), stage_status = "SEE_FIT_AND_DELTA_GATES_NO_INVARIANCE_CLAIM_FROM_CONVERGENCE"))
}
bs <- file.path(results_dir, "04_structural", "bootstrap_stage_status.csv")
if (file.exists(bs)) {
  b <- read_csv(bs, show_col_types = FALSE)
  stage_rows <- bind_rows(stage_rows, tibble(stage = "INDIRECT_BOOTSTRAP", attempted_models = 1L,
    converged_models = NA_integer_, admissible_models = NA_integer_, fit_passed_models = NA_integer_, stage_status = b$status[[1]]))
}
write_csv_local(stage_rows, "stage_status.csv")
# Hash y copia del postproceso: tambien forma parte de los artefactos verificables.
file_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
this_path <- normalizePath(gsub("~+~", " ", sub("^--file=", "", file_args[[1]]), fixed = TRUE), mustWork = TRUE)
invisible(file.copy(this_path, file.path(results_dir, "source_finalizer.R"), overwrite = TRUE))
write_csv_local(tibble(kind = "finalizer", path = this_path, md5 = unname(tools::md5sum(this_path))), "finalizer_hash.csv")
executed_snapshot <- file.path(results_dir, "executed_script_snapshot.R")
if (!file.exists(executed_snapshot) && file.exists(file.path(results_dir, "source_script.R"))) {
  invisible(file.copy(file.path(results_dir, "source_script.R"), executed_snapshot))
}
canonical_path <- file.path(dirname(this_path), "sem_reanalysis_final.R")
source_history <- tibble(kind = c("executed_model_script_snapshot", "current_canonical_for_future_runs", "executed_finalizer"),
  path = c(executed_snapshot, canonical_path, file.path(results_dir, "source_finalizer.R")),
  md5 = unname(tools::md5sum(c(executed_snapshot, canonical_path, file.path(results_dir, "source_finalizer.R")))))
write_csv_local(source_history, "source_version_history.csv")
files <- list.files(results_dir, recursive = TRUE, full.names = TRUE)
files <- files[!dir.exists(files) & basename(files) != "manifest.csv"]
write_csv_local(tibble(path = substring(files, nchar(results_dir)+2L), bytes = unname(file.info(files)$size),
                      md5 = unname(tools::md5sum(files))), "manifest.csv")
cat("CANONICAL_POSTPROCESS_COMPLETED\n")
