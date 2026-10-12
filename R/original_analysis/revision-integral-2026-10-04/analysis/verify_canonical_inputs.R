# Pruebas de regresion de etiquetas binarias y codificacion macOS de rutas.
# Ejecucion: R_LIBS_USER="biblioteca_privada:biblioteca_usuario" Rscript verify_canonical_inputs.R
suppressPackageStartupMessages({ library(haven); library(dplyr); library(tibble) })
script_file <- file.path(getwd(), "sem_reanalysis_final.R")
if (!file.exists(script_file)) stop("Ejecutar desde analysis/.")
expr <- parse(script_file)
selected <- c("ascii_text", "normalize_token", "recode_binary", "recode_yes_no", "binary_recode_rows", "get_script_path")
for (e in expr) if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=")) {
  lhs <- as.character(e[[2]])
  if (length(lhs) == 1L && lhs %in% selected) eval(e)
}
data_path <- Sys.getenv("SEM_SAV_PATH", "/Users/luisbenites/Dropbox/Paper con MA/Base_Final.sav")
x <- haven::read_sav(data_path)
stopifnot(identical(recode_yes_no(x$CSR1, "CSR1"), recode_yes_no(x$CSR1_rec, "CSR1_rec")))
stopifnot(sum(recode_yes_no(x$CSR1, "CSR1")) == 402)
inverted <- haven::labelled(c(1, 2, NA), c("No" = 1, "Sí" = 2))
stopifnot(identical(recode_yes_no(inverted, "inverted_labels"), c(0, 1, NA_real_)))
sex <- haven::labelled(c(1, 2), c("Hombre" = 1, "Mujer" = 2))
stopifnot(identical(recode_binary(sex, "Sexo", c("mujer"), c("hombre")), c(0, 1)))
err <- try(recode_yes_no(haven::labelled(c(1, 2), c("quizá" = 1, "ninguna" = 2)), "unknown"), silent = TRUE)
stopifnot(inherits(err, "try-error"))
stopifnot(identical(gsub("~+~", " ", "/Users/Test/Paper~+~con~+~MA/test.R", fixed = TRUE),
                    "/Users/Test/Paper con MA/test.R"))
stopifnot(identical(get_script_path(), normalizePath(file.path(getwd(), "verify_canonical_inputs.R"))))
# Ejecuta el bloque real de calculo metadata de la fuente, sin copiar la formula.
source_lines <- readLines(script_file, warn = FALSE)
start <- grep("valid_n_boot <- sum(valid)", source_lines, fixed = TRUE)[[1]]
end <- grep("safe_csv(diag, file.path", source_lines, fixed = TRUE)[[1]] - 1L
regression_env <- new.env(parent = globalenv())
regression_env$valid <- c(TRUE, TRUE, TRUE, FALSE)
regression_env$boot_mat <- matrix(0, nrow = 4, ncol = 2)
regression_env$boot_indirect <- 4L
regression_env$failures <- 4L
regression_env$nonadmissible <- integer()
regression_env$ncpus <- 1L
regression_env$bootstrap_result <- list(warnings = character())
eval(parse(text = source_lines[start:end]), regression_env)
stopifnot(regression_env$diag$valid[[1]] == 3,
          regression_env$diag$invalid[[1]] == 1,
          regression_env$diag$valid_fraction[[1]] == .75,
          !regression_env$diag$sufficient[[1]])
cat("BINARY_LABELS_PATH_AND_BOOTSTRAP_METADATA_REGRESSIONS_PASS\n")
