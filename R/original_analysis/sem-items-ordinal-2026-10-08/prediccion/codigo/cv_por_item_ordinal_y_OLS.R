# Sensibilidad exploratoria: predicción de PP1, PP2 y PP3 por separado.
# Conserva los 653 participantes y los pliegues 5 x 5 ya utilizados.
# No se recopilan datos, seleccionan predictores ni cambian especificaciones.
# El CFA ordinal y las inferencias SEM se realizan en otro script del paquete.
options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
project <- "/Users/luisbenites/Documents/Codex/Paper con MA - 1/reproduccion-integral-V3-2026-10-10/proyecto"
root <- file.path(project, "sem-items-ordinal-2026-10-08/prediccion")
output <- file.path(root, "resultados")
private <- file.path(output, "privados")
.libPaths(c(file.path(project, "revision-revista-2026-10-04/auditoria-codigos/R-library"),
           "/Users/luisbenites/Library/R/arm64/4.6/library", .libPaths()))
suppressPackageStartupMessages({library(lavaan); library(MASS); library(jsonlite); library(digest)})
dir.create(private, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(private, "fits"), recursive = TRUE, showWarnings = FALSE)
writeLines("RUNNING", file.path(output, "estado.txt"))
protected <- c(
  analytic_data = file.path(project, "revision-integral-2026-10-04/analysis/outputs/canonical-full-20261004-02/01_audit/analytic_data.rds"),
  P3 = file.path(project, "revision-integral-2026-10-04/sensitivities/P3_country_partial_demographics.rds"),
  prior_benchmark = file.path(project, "sem-alternativas-2026-10-08/codigo/benchmark_incremental_observado.R"),
  prior_cv = file.path(project, "sem-alternativas-2026-10-08/codigo/cross_validate_models.R"))
file_hash <- function(p) digest::digest(p, algo = "sha256", file = TRUE)
hash_before <- vapply(protected, file_hash, character(1))
old_p3 <- readRDS(protected[["P3"]])
audited <- as.data.frame(readRDS(protected[["analytic_data"]]))
# La matriz P3 contiene las recodificaciones originales y siete controles.
data <- as.data.frame(lavaan::lavInspect(old_p3, "data"))
common <- intersect(names(data), names(audited))
stopifnot(nrow(data) == 653L, nrow(audited) == 653L,
          max(abs(as.matrix(data[, common, drop = FALSE]) - as.matrix(audited[, common, drop = FALSE]))) == 0)
data$Pais <- audited$Pais
items <- c("PP1", "PP2", "PP3")
categories <- 1:5
controls <- c("C_Colombia", "C_Peru", "C_Mexico", "Sexo_num", "Estudios_num", "Edad_z", "NSE_z")
values <- c("Benevolencia", "Universalismo", "Poder", "Logro")
sc <- c("APS1", "APS2", "APS3", "APS4", "APS5", "APS6", "APS7", "APS9", "APS10", "APS14")
sbl <- c("APS11", "APS12", "APS13", "APS15", "APS16", "APS18", "APS19", "APS20", "APS21")
aps <- c(sc, sbl)
blocks <- list(B0_CSR_controls = c("CSR1_bin", controls),
               B1_add_values = c("CSR1_bin", controls, values),
               B2_add_APS = c("CSR1_bin", controls, values, aps),
               B3_APS_without_values = c("CSR1_bin", controls, aps))
all_variables <- unique(c(items, unlist(blocks, use.names = FALSE)))
stopifnot(!any(items %in% unlist(blocks, use.names = FALSE)),
          all(vapply(data[all_variables], is.numeric, logical(1))),
          !anyNA(data[all_variables]), all(is.finite(as.matrix(data[all_variables]))),
          all(as.matrix(data[items]) %in% categories))
source(protected[["prior_cv"]], local = TRUE)
strata <- .sem_cv_strata(data)
folds <- 5L; repeats <- 5L; seed <- 20261008L
assignment <- vector("list", repeats)
for (r in seq_len(repeats)) {
  set.seed(seed + r); a <- integer(nrow(data))
  for (s in unique(strata)) {
    rows <- which(strata == s)
    a[rows] <- sample(rep(seq_len(folds), length.out = length(rows)))
  }
  assignment[[r]] <- a
}
# Cotejo exacto contra las expresiones del generador SEM original, sin ajustar SEM.
expressions <- as.list(body(run_sem_cv))
ai <- which(vapply(expressions, function(e) grepl("^assignment <- vector", paste(deparse(e), collapse = " ")), logical(1)))
stopifnot(length(ai) == 1L)
reference <- new.env(parent = environment(run_sem_cv))
reference$data <- data; reference$strata <- strata; reference$repeats <- repeats
reference$folds <- folds; reference$seed <- seed
eval(expressions[[ai]], reference); eval(expressions[[ai + 1L]], reference)
stopifnot(identical(assignment, reference$assignment))
category_counts <- do.call(rbind, lapply(items, function(item) data.frame(item = item, category = categories,
                   n = as.integer(table(factor(data[[item]], levels = categories))))))
training_counts <- do.call(rbind, lapply(seq_len(repeats), function(r) do.call(rbind,
  lapply(seq_len(folds), function(k) do.call(rbind, lapply(items, function(item)
    data.frame(repetition = r, fold = k, item = item, category = categories,
               n_train = as.integer(table(factor(data[assignment[[r]] != k, item], levels = categories))),
               n_test = as.integer(table(factor(data[assignment[[r]] == k, item], levels = categories))))))))))
logloss_floor <- 1e-12
# Se documenta ANTES de los ajustes. Es un análisis posterior a explorar la base,
# no un prerregistro ni una validación externa. Las categorías raras se mantienen.
protocol <- list(status = "SPECIFIED_BEFORE_NEW_FITS", recorded_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  design = "Exploratory internal sensitivity after data and previous results examined",
  N = nrow(data), items = items, categories = categories, blocks = blocks,
  estimators = c("MASS polr probit", "OLS", "training empirical category baseline"),
  folds = folds, repetitions = repeats, seed = seed, exact_prior_assignment_verified = TRUE,
  outcome_predictor_leakage_prohibited = TRUE, predictors_standardized_using_training_only = TRUE,
  no_selection_regularization_or_tuning = TRUE, optimization = list(method = "BFGS", maxit = 5000L, reltol = 1e-10),
  initialization_rescue = "If default preliminary binary GLM fails to initialize polr, use zero predictor coefficients and qnorm of empirical cumulative training category proportions as starts, same ordinal model and optimization.",
  empty_category_rule = "Fit ordinal levels present in training; map probabilities to global levels 1-5 with zero for absent levels. Never alter folds or exclude observations.",
  n_empty_training_categories = sum(training_counts$n_train == 0),
  scoring = list(expected_score = "sum actual category value * probability", RMSE = "sqrt mean squared error", MAE = "mean absolute error",
    RPS = "mean sum over 4 cumulative cutpoints (Fpred - Fobserved)^2 divided by 4; smaller is better",
    logloss_raw = "mean negative log probability of observed category; Inf when exact zero",
    logloss_numerical = "mean negative log max(p observed,1e-12); numerical evaluation floor only, not a model change",
    prioritization = "RPS and common RMSE/MAE; clipped logloss has sensitivity to rare empty categories"),
  admissibility = "all folds required; convergence code 0, finite coefficients/ordered thresholds/probabilities and positive definite finite Hessian; OLS full rank",
  protected_paths = protected, original_sha256_before = hash_before)
protocol_path <- file.path(output, "protocolo_antes_ajustes.json")
if (file.exists(protocol_path)) {
  # Se conserva el protocolo original tras una incidencia de inicialización.
  # El rescate determinista se documenta antes de nuevos ajustes de reanudación.
  protocol_path <- file.path(output, "protocolo_antes_reanudacion.json")
}
jsonlite::write_json(protocol, protocol_path, pretty = TRUE, auto_unbox = TRUE, digits = NA)
write.csv(category_counts, file.path(output, "frecuencias_categorias.csv"), row.names = FALSE)
write.csv(training_counts, file.path(output, "frecuencias_por_pliegue.csv"), row.names = FALSE)
saveRDS(assignment, file.path(private, "asignaciones_cv.rds"))

metrics <- function(y, expected, p = NULL) {
  errors <- expected - y
  n <- length(y)
  out <- list(n_test = n, sum_squared_error = sum(errors^2), sum_absolute_error = sum(abs(errors)),
              RMSE = sqrt(mean(errors^2)), MAE = mean(abs(errors)),
              sum_RPS = NA_real_, RPS = NA_real_, sum_logloss_raw = NA_real_, logloss_raw = NA_real_,
              sum_logloss_floor = NA_real_, logloss_floor = NA_real_, n_probability_zero = NA_integer_, n_probability_below_floor = NA_integer_)
  if (!is.null(p)) {
    stopifnot(nrow(p) == n, ncol(p) == length(categories), all(is.finite(p)), all(p >= 0),
              all(p <= 1 + 1e-12), max(abs(rowSums(p) - 1)) < 1e-10)
    observed_p <- p[cbind(seq_len(n), match(y, categories))]
    cumulative <- t(apply(p, 1L, cumsum))[, seq_len(length(categories) - 1L), drop = FALSE]
    true_cumulative <- sapply(head(categories, -1L), function(cutpoint) as.numeric(y <= cutpoint))
    rps <- rowSums((cumulative - true_cumulative)^2) / (length(categories) - 1L)
    out$sum_RPS <- sum(rps); out$RPS <- mean(rps)
    out$sum_logloss_raw <- sum(-log(observed_p)); out$logloss_raw <- mean(-log(observed_p))
    out$sum_logloss_floor <- sum(-log(pmax(observed_p, logloss_floor)))
    out$logloss_floor <- out$sum_logloss_floor / n
    out$n_probability_zero <- sum(observed_p == 0)
    out$n_probability_below_floor <- sum(observed_p < logloss_floor)
  }
  out
}
fold_records <- list(); prediction_records <- list(); index <- 0L
record <- function(r, k, item, id, estimator, test, expected, p, diag) {
  index <<- index + 1L
  scored <- metrics(data[test, item], expected, p)
  fold_records[[index]] <<- cbind(data.frame(repetition = r, fold = k, item = item, model = id, estimator = estimator,
                               n_train = nrow(data) - length(test)), as.data.frame(scored), as.data.frame(diag))
  predicted <- data.frame(row_id = test, repetition = r, fold = k, item = item, model = id, estimator = estimator,
                         observed = data[test, item], predicted_expected = expected)
  for (c in categories) predicted[[paste0("probability_", c)]] <- if (is.null(p)) NA_real_ else p[, c]
  prediction_records[[index]] <<- predicted
}
for (r in seq_len(repeats)) {
  for (k in seq_len(folds)) {
    test <- which(assignment[[r]] == k); train <- which(assignment[[r]] != k)
    for (item in items) {
      counts <- table(factor(data[train, item], levels = categories))
      present <- categories[counts > 0]
      baseline_p <- matrix(rep(as.numeric(counts) / length(train), each = length(test)), nrow = length(test),
                           dimnames = list(NULL, as.character(categories)))
      baseline_expected <- as.vector(baseline_p %*% categories)
      empty_diag <- list(converged = TRUE, coefficient_rank = 0L, min_hessian_eigen = NA_real_, hessian_condition_number = NA_real_,
                         coefficients_finite = TRUE, threshold_ordered = TRUE, missing_training_categories = paste(setdiff(categories, present), collapse = ";"),
                         warnings = "", optimization_function_calls = NA_integer_, initialization_rescued = FALSE)
      record(r, k, item, "B_empirical", "empirical", test, baseline_expected, baseline_p, empty_diag)
      for (id in names(blocks)) {
        predictors <- blocks[[id]]
        x_train <- data[train, predictors, drop = FALSE]; x_test <- data[test, predictors, drop = FALSE]
        stopifnot(!any(items %in% names(x_test)))
        # Se aprende la escala exclusivamente en entrenamiento, para ambos estimadores.
        center <- vapply(x_train, mean, numeric(1)); scale <- vapply(x_train, stats::sd, numeric(1))
        stopifnot(all(is.finite(center)), all(is.finite(scale)), all(scale > 0))
        z_train <- as.data.frame(sweep(sweep(as.matrix(x_train), 2L, center, "-"), 2L, scale, "/"))
        z_test <- as.data.frame(sweep(sweep(as.matrix(x_test), 2L, center, "-"), 2L, scale, "/"))
        training <- data.frame(response = ordered(data[train, item], levels = present), z_train, check.names = FALSE)
        warnings <- character()
        ordinal_path <- file.path(private, "fits", sprintf("%s_%s_ordinal_r%02d_f%02d.rds", item, id, r, k))
        rescued <- FALSE
        if (file.exists(ordinal_path)) {
          cached <- readRDS(ordinal_path)
          stopifnot(identical(cached$training_rows, train), identical(cached$testing_rows, test),
                    identical(cached$predictors, predictors), identical(cached$training_center, center),
                    identical(cached$training_scale, scale), identical(cached$category_support, present))
          fit <- cached$fit
          warnings <- cached$diagnostics$warnings
          rescued <- isTRUE(cached$diagnostics$initialization_rescued)
        } else {
          fit <- tryCatch(withCallingHandlers(MASS::polr(response ~ ., data = training, method = "probit", Hess = TRUE,
                         model = TRUE, control = list(maxit = 5000L, reltol = 1e-10)),
                      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
                     error = function(e) e)
          if (inherits(fit, "error")) {
            if (!grepl("attempt to find suitable starting values failed", conditionMessage(fit), fixed = TRUE)) stop(fit)
            # Rescate del inicializador, sin cambiar predictor, categorías ni likelihood.
            # Los umbrales iniciales provienen sólo de frecuencias del entrenamiento.
            starts <- c(rep(0, length(predictors)), qnorm(head(cumsum(as.numeric(counts[counts > 0])) / length(train), -1L)))
            stopifnot(all(is.finite(starts)))
            rescued <- TRUE
            fit <- withCallingHandlers(MASS::polr(response ~ ., data = training, method = "probit", Hess = TRUE,
                         model = TRUE, start = starts, control = list(maxit = 5000L, reltol = 1e-10)),
                      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
          }
        }
        hessian_eigen <- eigen((fit$Hessian + t(fit$Hessian)) / 2, symmetric = TRUE, only.values = TRUE)$values
        stopifnot(fit$convergence == 0L, all(is.finite(fit$coefficients)), all(is.finite(fit$zeta)),
                  all(diff(fit$zeta) > 0), all(is.finite(fit$Hessian)), all(is.finite(hessian_eigen)), min(hessian_eigen) > 0)
        fitted_p <- as.matrix(predict(fit, newdata = z_test, type = "probs"))
        p <- matrix(0, nrow = length(test), ncol = length(categories), dimnames = list(NULL, as.character(categories)))
        p[, colnames(fitted_p)] <- fitted_p
        expected <- as.vector(p %*% categories)
        diagnostics <- list(converged = fit$convergence == 0L, coefficient_rank = length(fit$coefficients),
          min_hessian_eigen = min(hessian_eigen), hessian_condition_number = max(hessian_eigen) / min(hessian_eigen),
          coefficients_finite = all(is.finite(c(fit$coefficients, fit$zeta))), threshold_ordered = all(diff(fit$zeta) > 0),
          missing_training_categories = paste(setdiff(categories, present), collapse = ";"),
          warnings = paste(unique(warnings[nzchar(warnings)]), collapse = " | "), optimization_function_calls = as.integer(fit$niter[1L]),
          initialization_rescued = rescued)
        record(r, k, item, id, "ordinal_probit", test, expected, p, diagnostics)
        wrapper <- list(fit = fit, training_rows = train, testing_rows = test, predictors = predictors,
                        training_center = center, training_scale = scale, category_support = present, diagnostics = diagnostics)
        saveRDS(wrapper, ordinal_path)
        # OLS conserva los mismos predictores y puntajes originales; no recorta predicciones.
        ols_training <- data.frame(response = data[train, item], z_train, check.names = FALSE)
        ols_path <- file.path(private, "fits", sprintf("%s_%s_OLS_r%02d_f%02d.rds", item, id, r, k))
        if (file.exists(ols_path)) {
          cached <- readRDS(ols_path)
          stopifnot(identical(cached$training_rows, train), identical(cached$testing_rows, test),
                    identical(cached$predictors, predictors), identical(cached$training_center, center),
                    identical(cached$training_scale, scale))
          ols <- cached$fit
        } else ols <- stats::lm(response ~ ., data = ols_training, na.action = stats::na.fail)
        stopifnot(ols$rank == length(predictors) + 1L, all(is.finite(coef(ols))))
        ols_expected <- as.numeric(predict(ols, newdata = z_test))
        stopifnot(all(is.finite(ols_expected)))
        ols_diag <- empty_diag; ols_diag$coefficient_rank <- ols$rank
        record(r, k, item, id, "OLS", test, ols_expected, NULL, ols_diag)
        wrapper$fit <- ols; wrapper$diagnostics <- ols_diag
        saveRDS(wrapper, ols_path)
      }
    }
    cat(sprintf("Completed repetition %d/%d fold %d/%d: all 3 items and 4 blocks\n", r, repeats, k, folds))
    flush.console()
  }
}
fold_table <- do.call(rbind, fold_records)
predictions <- do.call(rbind, prediction_records)
stopifnot(nrow(fold_table) == repeats * folds * length(items) * (2L * length(blocks) + 1L),
          nrow(predictions) == repeats * nrow(data) * length(items) * (2L * length(blocks) + 1L),
          all(fold_table$converged), all(fold_table$coefficients_finite), all(fold_table$threshold_ordered))
write.csv(fold_table, file.path(output, "metricas_y_diagnosticos_por_pliegue.csv"), row.names = FALSE)
saveRDS(predictions, file.path(private, "predicciones_fuera_muestra.rds"))
# Agregar primero las sumas de error: evita promediar incorrectamente RMSE de folds.
keys <- unique(fold_table[, c("repetition", "item", "model", "estimator")])
repetition_table <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
  key <- keys[i, ]; x <- fold_table[fold_table$repetition == key$repetition & fold_table$item == key$item &
                                    fold_table$model == key$model & fold_table$estimator == key$estimator, ]
  n <- sum(x$n_test); stopifnot(nrow(x) == folds, n == nrow(data))
  cbind(key, data.frame(folds_valid = nrow(x), n_test = n, sum_squared_error = sum(x$sum_squared_error),
     sum_absolute_error = sum(x$sum_absolute_error), RMSE = sqrt(sum(x$sum_squared_error) / n),
     MAE = sum(x$sum_absolute_error) / n, RPS = sum(x$sum_RPS) / n,
     logloss_raw = sum(x$sum_logloss_raw) / n, logloss_floor = sum(x$sum_logloss_floor) / n,
     n_probability_zero = sum(x$n_probability_zero), n_probability_below_floor = sum(x$n_probability_below_floor)))
}))
overall_keys <- unique(repetition_table[, c("item", "model", "estimator")])
overall <- do.call(rbind, lapply(seq_len(nrow(overall_keys)), function(i) {
  key <- overall_keys[i, ]; x <- repetition_table[repetition_table$item == key$item &
                                                 repetition_table$model == key$model & repetition_table$estimator == key$estimator, ]
  stopifnot(nrow(x) == repeats)
  cbind(key, data.frame(repetitions_valid = nrow(x), total_test_evaluations = sum(x$n_test),
    RMSE = mean(x$RMSE), RMSE_partition_SD = sd(x$RMSE), RMSE_pooled = sqrt(sum(x$sum_squared_error) / sum(x$n_test)),
    MAE = mean(x$MAE), MAE_partition_SD = sd(x$MAE), RPS = mean(x$RPS), RPS_partition_SD = sd(x$RPS),
    logloss_raw = mean(x$logloss_raw), logloss_floor = mean(x$logloss_floor), logloss_floor_partition_SD = sd(x$logloss_floor),
    n_probability_zero = sum(x$n_probability_zero), n_probability_below_floor = sum(x$n_probability_below_floor)))
}))
pairs <- list(values_added = c("B1_add_values", "B0_CSR_controls"), APS_added_after_values = c("B2_add_APS", "B1_add_values"),
              values_and_APS_added = c("B2_add_APS", "B0_CSR_controls"), values_added_after_APS = c("B2_add_APS", "B3_APS_without_values"),
              all_information_vs_empirical = c("B2_add_APS", "B_empirical"))
paired <- do.call(rbind, lapply(items, function(item) do.call(rbind, lapply(c("ordinal_probit", "OLS"), function(estimator)
  do.call(rbind, lapply(names(pairs), function(comparison) {
    pair <- pairs[[comparison]]
    a <- repetition_table[repetition_table$item == item & repetition_table$estimator == estimator & repetition_table$model == pair[1], ]
    # Para OLS, el baseline de medias es el mismo baseline empírico calculado antes.
    reference_estimator <- if (pair[2] == "B_empirical") "empirical" else estimator
    b <- repetition_table[repetition_table$item == item & repetition_table$estimator == reference_estimator & repetition_table$model == pair[2], ]
    stopifnot(nrow(a) == repeats, identical(a$repetition, b$repetition))
    data.frame(item = item, estimator = estimator, comparison = comparison, model = pair[1], reference_model = pair[2],
      repetition = a$repetition, difference_RMSE = a$RMSE - b$RMSE, relative_RMSE_change_percent = 100 * (a$RMSE / b$RMSE - 1),
      difference_MAE = a$MAE - b$MAE, difference_RPS = if (estimator == "OLS") NA_real_ else a$RPS - b$RPS,
      difference_logloss_floor = if (estimator == "OLS") NA_real_ else a$logloss_floor - b$logloss_floor)
  }))))))
paired_keys <- unique(paired[, c("item", "estimator", "comparison", "model", "reference_model")])
paired_summary <- do.call(rbind, lapply(seq_len(nrow(paired_keys)), function(i) {
  key <- paired_keys[i, ]; x <- paired[paired$item == key$item & paired$estimator == key$estimator & paired$comparison == key$comparison, ]
  cbind(key, data.frame(repetitions = nrow(x), mean_difference_RMSE = mean(x$difference_RMSE),
     SD_difference_RMSE = sd(x$difference_RMSE), min_difference_RMSE = min(x$difference_RMSE), max_difference_RMSE = max(x$difference_RMSE),
     mean_relative_RMSE_change_percent = mean(x$relative_RMSE_change_percent), mean_difference_MAE = mean(x$difference_MAE),
     mean_difference_RPS = mean(x$difference_RPS), mean_difference_logloss_floor = mean(x$difference_logloss_floor),
     repetitions_better_RMSE = sum(x$difference_RMSE < 0), repetitions_better_MAE = sum(x$difference_MAE < 0),
     repetitions_better_RPS = if (key$estimator == "OLS") NA_integer_ else sum(x$difference_RPS < 0)))
}))
write.csv(repetition_table, file.path(output, "metricas_por_repeticion.csv"), row.names = FALSE)
write.csv(overall, file.path(output, "resumen_modelos_por_item.csv"), row.names = FALSE)
write.csv(paired, file.path(output, "diferencias_pareadas_por_repeticion.csv"), row.names = FALSE)
write.csv(paired_summary, file.path(output, "resumen_diferencias_pareadas.csv"), row.names = FALSE)
# Se comparan todos los tamaños de fold con los 100 registros del benchmark previo.
prior <- read.csv(file.path(project, "sem-alternativas-2026-10-08/resultados/benchmark_incremental/benchmark_fold_metrics.csv"))
if (!all(c("repetition", "fold", "n_train", "n_test") %in% names(prior))) stop("Unexpected prior benchmark columns")
new_fold <- unique(fold_table[, c("repetition", "fold", "n_train", "n_test")])
comparison <- merge(prior[, c("repetition", "fold", "n_train", "n_test")], new_fold,
                    by = c("repetition", "fold"), suffixes = c("_prior", "_new"), all.x = TRUE)
stopifnot(nrow(comparison) == nrow(prior), !anyNA(comparison),
          all(comparison$n_train_prior == comparison$n_train_new), all(comparison$n_test_prior == comparison$n_test_new))
hash_after <- vapply(protected, file_hash, character(1))
stopifnot(identical(hash_before, hash_after))
audit <- list(status = "COMPLETED", N = nrow(data), protected_paths = protected, sha256_before = hash_before,
  sha256_after = hash_after, originals_unchanged = identical(hash_before, hash_after),
  exact_prior_assignments_verified = TRUE, prior_benchmark_fold_sizes_verified = TRUE,
  n_ordinal_fits = sum(fold_table$estimator == "ordinal_probit"), n_OLS_fits = sum(fold_table$estimator == "OLS"),
  n_empirical_folds = sum(fold_table$estimator == "empirical"), n_private_predictions = nrow(predictions),
  all_training_converged = all(fold_table$converged), all_coefficients_finite = all(fold_table$coefficients_finite),
  all_thresholds_ordered = all(fold_table$threshold_ordered),
  all_ordinal_Hessians_positive_definite = all(fold_table$min_hessian_eigen[fold_table$estimator == "ordinal_probit"] > 0),
  minimum_ordinal_Hessian_eigenvalue = min(fold_table$min_hessian_eigen, na.rm = TRUE),
  maximum_ordinal_Hessian_condition_number = max(fold_table$hessian_condition_number, na.rm = TRUE),
  fits_with_warnings = sum(nzchar(fold_table$warnings)), complete_folds_no_omission = TRUE,
  fits_using_deterministic_starting_value_rescue = sum(fold_table$initialization_rescued),
  no_inference_tests_on_overlapping_CV = TRUE, n_empty_training_categories = sum(training_counts$n_train == 0),
  caveat = "Two training folds lack category 1. Zero mass yields infinite raw logloss. Finite logloss uses predetermined 1e-12 scoring floor only. RPS and RMSE/MAE retained.")
jsonlite::write_json(audit, file.path(output, "auditoria_final.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA)
capture.output(sessionInfo(), file = file.path(output, "sessionInfo.txt"))
writeLines("COMPLETED", file.path(output, "estado.txt"))
cat("Completed all ordinal, OLS and empirical folds; protected sources unchanged.\n")
