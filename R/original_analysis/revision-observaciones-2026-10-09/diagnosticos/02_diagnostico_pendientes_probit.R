# Observación 4: diagnóstico PARCIAL de pendientes comunes del probit acumulativo.
# Tablas: frecuencias_categorias.csv, bases_probit_comunes.csv,
# pruebas_nominales_valores.csv y parametros_probit_diagnosticos.csv.
# Se ejecutan 3 bases B2 y 12 alternativas (4 valores x 3 PP), fijadas antes
# de estos ajustes en protocolo_antes_diagnosticos.json. Resultados previos conocidos:
# esto no es prerregistro, validación externa ni un diagnóstico global de 31 efectos.
options(stringsAsFactors = FALSE)
read_only_existing <- '--diagnose-existing' %in% commandArgs(trailingOnly = TRUE)
Sys.setenv(VECLIB_MAXIMUM_THREADS = '1', OPENBLAS_NUM_THREADS = '1', OMP_NUM_THREADS = '1')
project <- '/Users/luisbenites/Documents/Codex/Paper con MA - 1/reproduccion-integral-V3-2026-10-10/proyecto'
dest <- file.path(project, 'revision-observaciones-2026-10-09/diagnosticos')
.libPaths(c(file.path(dest, 'R-library'),
            file.path(project, 'revision-revista-2026-10-04/auditoria-codigos/R-library'),
            '/Users/luisbenites/Library/R/arm64/4.6/library', .libPaths()))
suppressPackageStartupMessages({library(ordinal); library(lavaan); library(jsonlite); library(digest)})
stopifnot(file.exists(file.path(dest, 'protocolo_antes_diagnosticos.json')))
out <- file.path(dest, 'resultados'); private <- file.path(dest, 'privados')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(private, recursive = TRUE, showWarnings = FALSE)
write_out <- function(x, n) write.csv(x, file.path(out, n), row.names = FALSE, na = '')

# P3 mantiene las recodificaciones exactas del artículo. La segunda fuente auditada
# se usa únicamente para verificar igualdad de datos. Ninguna fuente se modifica.
p3file <- file.path(project, 'revision-integral-2026-10-04/sensitivities/P3_country_partial_demographics.rds')
auditfile <- file.path(project, 'revision-integral-2026-10-04/analysis/outputs/canonical-full-20261004-02/01_audit/analytic_data.rds')
raw <- as.data.frame(lavInspect(readRDS(p3file), 'data'))
audited <- as.data.frame(readRDS(auditfile))
shared <- intersect(names(raw), names(audited))
stopifnot(max(abs(as.matrix(raw[, shared]) - as.matrix(audited[, shared]))) == 0,
          nrow(raw) == 653L, !anyNA(raw), all(is.finite(as.matrix(raw))))
controls <- c('C_Colombia','C_Peru','C_Mexico','Sexo_num','Estudios_num','Edad_z','NSE_z')
values <- c('Benevolencia','Universalismo','Poder','Logro')
aps <- c('APS1','APS2','APS3','APS4','APS5','APS6','APS7','APS9','APS10','APS14',
         'APS11','APS12','APS13','APS15','APS16','APS18','APS19','APS20','APS21')
predictors <- c('CSR1_bin', controls, values, aps)
outcomes <- c('PP1', 'PP2', 'PP3')
stopifnot(length(predictors) == 31L, all(c(predictors, outcomes) %in% names(raw)))

# Escalado en la muestra completa SOLO para este diagnóstico de ajuste.
# No produce predicción de casos retenidos ni sustituye los pliegues CV del artículo.
dat <- raw
centers <- vapply(raw[, predictors], mean, numeric(1))
scales <- vapply(raw[, predictors], sd, numeric(1))
stopifnot(all(scales > 0))
dat[, predictors] <- as.data.frame(scale(raw[, predictors], center = centers, scale = scales))
for (y in outcomes) {
  stopifnot(all(raw[[y]] %in% 1:5))
  dat[[y]] <- ordered(raw[[y]], levels = 1:5)
}
write_out(data.frame(predictor = predictors, center = centers, scale = scales), 'escalado_solo_diagnostico.csv')
freq <- do.call(rbind, lapply(outcomes, function(y) {
  tt <- table(dat[[y]]); data.frame(outcome = y, category = 1:5, n = as.integer(tt), proportion = as.integer(tt)/653)
}))
write_out(freq, 'frecuencias_categorias.csv')
stopifnot(identical(freq$n[freq$category == 1], c(3L, 4L, 2L)))

control <- clm.control(maxIter = 200, maxLineIter = 100, gradTol = 1e-7,
                       relTol = 1e-7, convergence = 'warn')
safe <- function(expr) {
  warnings <- character()
  fit <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning')
  }), error = function(e) e)
  list(fit = fit, warnings = unique(warnings), error = if (inherits(fit, 'error')) conditionMessage(fit) else '')
}
conv_code <- function(f) {
  x <- f$convergence
  if (is.list(x)) {
    if (!is.null(x$code)) return(as.integer(x$code))
    if (!is.null(x$convergence)) return(as.integer(x$convergence))
    return(NA_integer_)
  }
  if (is.numeric(x)) return(as.integer(x)[1L])
  NA_integer_
}

# Aplicar los criterios definidos ANTES de las alternativas. Una prueba fallida
# queda no evaluable (p=NA). No se colapsan categorías ni se suprimen predictores.
# La columna de ubicación del valor nominal puede quedar estructuralmente aliased:
# sus 4 pendientes de umbral la reemplazan, con +3 parámetros netos. Se comprueba
# que no existan alias adicionales y se revisan sólo parámetros estimables.
diagnose <- function(z, expected_alias) {
  empty <- data.frame(N = NA_integer_, edf = NA_integer_, convergence_code = NA_integer_,
    finite_coef = FALSE, finite_SE = FALSE, aliased_count = NA_integer_, expected_alias = expected_alias,
    max_gradient = NA_real_, Hessian_rank = NA_integer_, Hessian_dimension = NA_integer_,
    Hessian_min_eigen = NA_real_, Hessian_max_eigen = NA_real_, Hessian_tolerance = NA_real_,
    condition_number = NA_real_, probability_min = NA_real_, probability_max = NA_real_,
    simplex_max_error = NA_real_, ordered_cumulative_probabilities = FALSE,
    ordered_thresholds_observed = FALSE, logLik = NA_real_, admissible = FALSE,
    diagnostic_reason = z$error, warnings = paste(z$warnings, collapse = ' | '))
  if (!inherits(z$fit, 'clm')) return(empty)
  f <- z$fit; r <- empty
  co <- coef(f, na.rm = TRUE)
  H <- as.matrix(f$Hessian)
  eig <- if (all(is.finite(H))) eigen((H+t(H))/2, symmetric = TRUE, only.values = TRUE)$values else NA_real_
  tol <- if (all(is.finite(eig))) 1e-8 * max(1, max(abs(eig))) else NA_real_
  vv <- safe(vcov(f)); se <- if (is.matrix(vv$fit)) sqrt(diag(vv$fit)) else NA_real_
  pr <- safe(predict(f, newdata = dat[, predictors, drop = FALSE], type = 'prob')$fit)
  pp <- pr$fit
  valid_shape <- is.matrix(pp) && identical(dim(pp), c(653L, 5L))
  finite_pp <- valid_shape && all(is.finite(pp))
  # Las probabilidades acumuladas deben aumentar; equivalen al orden de los
  # umbrales probit para cada perfil observado. Se revisa además el predictor
  # lineal de los 4 umbrales (eta1: límite superior; eta2: límite inferior).
  cumulative_ok <- FALSE; threshold_ok <- FALSE
  if (finite_pp) {
    cumulative <- t(apply(pp, 1L, cumsum))
    cumulative_ok <- all(pp >= -1e-10) && all(apply(cumulative, 1L, diff) >= -1e-10)
  }
  lp <- safe(predict(f, newdata = dat[, predictors, drop = FALSE], type = 'linear.predictor'))
  if (is.list(lp$fit) && is.matrix(lp$fit$eta1) && is.matrix(lp$fit$eta2)) {
    e1 <- lp$fit$eta1; e2 <- lp$fit$eta2
    # Los extremos se representan con ±1e5 en predict.clm. Los 4 umbrales
    # interiores eta1[,1:4] (equivalentes a eta2[,2:5]) deben ser finitos
    # y estrictamente crecientes; no incluimos los extremos artificiales.
    if (nrow(e1) == 653L && ncol(e1) == 5L) {
      interior <- e1[, 1:4, drop = FALSE]
      threshold_ok <- all(is.finite(interior)) && all(apply(interior, 1L, diff) > 0)
      threshold_ok <- threshold_ok && max(abs(interior - e2[, 2:5, drop = FALSE])) < 1e-8
    }
  }
  r$N <- nobs(f); r$edf <- f$edf; r$convergence_code <- conv_code(f)
  r$finite_coef <- all(is.finite(co)); r$finite_SE <- all(is.finite(se)) && length(se) == length(co)
  r$aliased_count <- sum(unlist(f$aliased)); r$max_gradient <- max(abs(f$gradient))
  r$Hessian_rank <- if (is.finite(tol)) sum(abs(eig) > tol) else NA_integer_
  r$Hessian_dimension <- nrow(H); r$Hessian_min_eigen <- min(eig); r$Hessian_max_eigen <- max(eig)
  r$Hessian_tolerance <- tol
  r$condition_number <- if (all(is.finite(eig)) && min(eig) > 0) max(eig)/min(eig) else NA_real_
  r$probability_min <- if (finite_pp) min(pp) else NA_real_
  r$probability_max <- if (finite_pp) max(pp) else NA_real_
  r$simplex_max_error <- if (finite_pp) max(abs(rowSums(pp)-1)) else NA_real_
  r$ordered_cumulative_probabilities <- cumulative_ok; r$ordered_thresholds_observed <- threshold_ok
  r$logLik <- as.numeric(logLik(f))
  checks <- c(convergence = isTRUE(r$convergence_code == 0L), finite_coefficients = r$finite_coef,
    finite_SE = r$finite_SE, no_unexpected_alias = isTRUE(r$aliased_count == expected_alias),
    gradient = isTRUE(r$max_gradient <= 1e-5), Hessian_full_rank = isTRUE(r$Hessian_rank == nrow(H)),
    Hessian_positive = isTRUE(r$Hessian_min_eigen > tol), probability_simplex = finite_pp &&
      min(pp) >= -1e-10 && max(pp) <= 1+1e-10 && max(abs(rowSums(pp)-1)) <= 1e-8,
    cumulative_order = cumulative_ok, threshold_order = threshold_ok, sample = r$N == 653L)
  r$admissible <- all(checks)
  r$diagnostic_reason <- if (r$admissible) 'All preplanned numerical criteria pass' else
    paste(names(checks)[!checks], collapse = '; ')
  r$warnings <- paste(unique(c(z$warnings, vv$warnings, pr$warnings, lp$warnings,
      vv$error, pr$error, lp$error))[nzchar(unique(c(z$warnings, vv$warnings, pr$warnings, lp$warnings,
      vv$error, pr$error, lp$error)))], collapse = ' | ')
  r
}
parameter_rows <- list()
save_fit <- function(z, key) {
  saveRDS(z, file.path(private, paste0(key, '.rds')))
  if (inherits(z$fit, 'clm')) {
    f <- z$fit; co <- coef(f, na.rm = TRUE)
    v <- safe(vcov(f)); se <- if (is.matrix(v$fit)) sqrt(diag(v$fit)) else rep(NA_real_, length(co))
    parameter_rows[[key]] <<- data.frame(fit = key, parameter = names(co), estimate = unname(co),
      SE = unname(se), lower_95_Wald = unname(co - qnorm(.975)*se), upper_95_Wald = unname(co + qnorm(.975)*se))
  }
}

bases <- list(); altrows <- list()
for (y in outcomes) {
  formula <- reformulate(predictors, response = y)
  cat('BASELINE', y, '\n'); flush.console()
  z0 <- if (read_only_existing) readRDS(file.path(private, paste0(y, '_B2_comun.rds'))) else
    safe(clm(formula = formula, data = dat, link = 'probit', threshold = 'flexible',
             control = control, na.action = na.fail))
  save_fit(z0, paste0(y, '_B2_comun'))
  d0 <- diagnose(z0, expected_alias = 0L)
  bases[[y]] <- cbind(data.frame(outcome = y, model = 'B2 common slopes; 31 predictors'), d0)
  write_out(do.call(rbind, bases), 'bases_probit_comunes.csv')
  for (value in values) {
    cat('NOMINAL', y, value, '\n'); flush.console()
    za <- if (read_only_existing) readRDS(file.path(private, paste0(y, '_', value, '.rds'))) else
      safe(clm(formula = formula, nominal = reformulate(value), data = dat,
               link = 'probit', threshold = 'flexible', control = control, na.action = na.fail))
    key <- paste(y, value, sep = '_'); save_fit(za, key)
    da <- diagnose(za, expected_alias = 1L)
    df_add <- da$edf - d0$edf
    ll_improvement <- da$logLik - d0$logLik
    same_sample <- inherits(z0$fit, 'clm') && inherits(za$fit, 'clm') &&
      identical(rownames(z0$fit$model), rownames(za$fit$model))
    lr_valid <- isTRUE(d0$admissible) && isTRUE(da$admissible) && isTRUE(df_add == 3L) &&
      isTRUE(same_sample) && isTRUE(ll_improvement >= -1e-7)
    LR <- if (lr_valid) 2*max(0, ll_improvement) else NA_real_
    p <- if (lr_valid) pchisq(LR, df = df_add, lower.tail = FALSE) else NA_real_
    altrows[[key]] <- cbind(data.frame(outcome = y, nominal_value = value, planned_family_size = 12L,
      baseline_admissible = d0$admissible, baseline_logLik = d0$logLik, baseline_edf = d0$edf,
      df_added = df_add, same_sample = same_sample, logLik_gain = ll_improvement,
      LR_inference_valid = lr_valid, LR = LR, p_raw = p, p_Holm_12 = NA_real_), da)
    write_out(do.call(rbind, altrows), 'pruebas_nominales_valores.csv')
  }
}
tbl <- do.call(rbind, altrows)
valid <- is.finite(tbl$p_raw)
# Familia PLANIFICADA de 12, incluyendo fallos como hipótesis no evaluadas.
# p.adjust con n=12 impide reducir la penalización por haber fallado alternativas.
tbl$p_Holm_12[valid] <- p.adjust(tbl$p_raw[valid], method = 'holm', n = 12L)
write_out(tbl, 'pruebas_nominales_valores.csv')
write_out(do.call(rbind, parameter_rows), 'parametros_probit_diagnosticos.csv')
writeLines(capture.output(sessionInfo()), file.path(out, 'sessionInfo_probit.txt'))
write_json(list(completed_at = format(Sys.time(), '%Y-%m-%d %H:%M:%S %Z'),
  ordinal_version = as.character(packageVersion('ordinal')), attempts = 12L,
  valid_partial_tests = sum(valid), failed_partial_tests = sum(!valid),
  Holm_family = 12L, no_global_assumption_test = TRUE,
  source_P3_sha256 = digest(p3file, algo = 'sha256', file = TRUE),
  source_data_sha256 = digest(auditfile, algo = 'sha256', file = TRUE),
  limitation = 'Only four value effects were examined conditional on unchanged B2; rare category counts 3/4/2 weaken asymptotic tests. Failure is not evidence of common slopes.'),
  file.path(out, 'estado_probit.json'), pretty = TRUE, auto_unbox = TRUE)
cat('COMPLETED 12 partial diagnostics;', sum(valid), 'valid tests;', sum(!valid), 'not evaluable\n')
