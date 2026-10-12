# Sensibilidad exploratoria: PP1, PP2 y PP3 como respuestas separadas.
# Conserva exactamente la muestra, los controles P3 y el bloque APS de M2.
# El modelo principal anterior NO se modifica: este script estima dos modelos
# adicionales fijados antes del ajuste, sin seleccionar rutas ni borrar ítems.
options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           OMP_NUM_THREADS = "1")
project <- "/Users/luisbenites/Documents/Codex/Paper con MA - 1/reproduccion-integral-V3-2026-10-10/proyecto"
root <- file.path(project, "sem-items-ordinal-2026-10-08/items")
out <- file.path(root, "resultados")
private <- file.path(root, "privados")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(private, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(file.path(project, "revision-revista-2026-10-04/auditoria-codigos/R-library"),
            "/Users/luisbenites/Library/R/arm64/4.6/library", .libPaths()))
suppressPackageStartupMessages({library(lavaan); library(jsonlite); library(digest)})
state <- file.path(out, "estado.txt")
resume_exports <- "--reanudar-exportacion" %in% commandArgs(trailingOnly = TRUE)
writeLines("STARTED", state)
options(error = function() {
  writeLines(c("FAILED", geterrmessage()), state)
  quit(status = 1, save = "no")
})
source_data <- file.path(project, "revision-integral-2026-10-04/analysis/outputs/canonical-full-20261004-02/01_audit/analytic_data.rds")
source_p3 <- file.path(project, "revision-integral-2026-10-04/sensitivities/P3_country_partial_demographics.rds")
source_syntax <- file.path(project, "sem-alternativas-2026-10-08/resultados/M2_puentes2_syntax.txt")
source_reference <- file.path(project, "sem-alternativas-2026-10-08/resultados/M2_puentes2_SEM_MLR.rds")
protected <- c(source_data, source_p3, source_syntax, source_reference)
sha <- function(p) digest::digest(p, algo = "sha256", file = TRUE)
before <- vapply(protected, sha, character(1))

# Se lee la matriz YA recodificada de P3, incluidos z de edad/NSE y dummies.
# La base auditada se utiliza sólo para comprobar identidad y recuperar país.
old_p3 <- readRDS(source_p3)
audited <- as.data.frame(readRDS(source_data))
data <- as.data.frame(lavInspect(old_p3, "data"))
common <- intersect(names(data), names(audited))
stopifnot(nrow(data) == 653L, nrow(audited) == nrow(data),
          max(abs(as.matrix(data[, common, drop = FALSE]) -
                    as.matrix(audited[, common, drop = FALSE]))) == 0)
data$Pais <- audited$Pais
pp <- c("PP1", "PP2", "PP3")
sc <- c("APS1", "APS2", "APS3", "APS4", "APS5", "APS6", "APS7", "APS9", "APS10", "APS14")
sbl <- c("APS11", "APS12", "APS13", "APS15", "APS16", "APS18", "APS19", "APS20", "APS21")
aps <- c(sc, sbl)
values <- c("Benevolencia", "Universalismo", "Poder", "Logro")
focal <- c("SC", "SBL", values)
controls <- c("CSR1_bin", "C_Colombia", "C_Peru", "C_Mexico", "Sexo_num", "Estudios_num", "Edad_z", "NSE_z")
stopifnot(!anyNA(data[, c(aps, pp, values, controls)]),
          all(vapply(data[, c(aps, pp, values, controls)], is.numeric, logical(1))),
          all(is.finite(as.matrix(data[, c(aps, pp, values, controls)]))),
          all(as.matrix(data[, c(aps, pp)]) %in% 1:5))
saveRDS(data, file.path(private, "datos_P3_para_items_y_CV.rds"))

# Se conservan las seis regresiones upstream y las covarianzas originales.
# Se reemplaza PSB por tres regresiones con todos los mismos predictores.
# Los residuos PP se correlacionan libremente para retener su dependencia.
original <- readLines(source_syntax)
keep <- !grepl("^PSB ", original) & !grepl(":=", original, fixed = TRUE)
syntax_lines <- original[keep]
label <- function(item, predictor) paste0("item", match(item, pp), "_", predictor)
for (item in pp) {
  slopes <- c(paste0(vapply(focal, function(x) label(item, x), character(1)), "*", focal), controls)
  syntax_lines <- c(syntax_lines, paste(item, "~", paste(slopes, collapse = " + ")))
}
syntax_lines <- c(syntax_lines, "PP1 ~~ PP2 + PP3", "PP2 ~~ PP3")
syntax <- paste(syntax_lines, collapse = "\n")
stopifnot(!grepl("PSB", syntax, fixed = TRUE),
          sum(grepl(" ~ ", syntax_lines, fixed = TRUE)) == 9L)
writeLines(syntax_lines, file.path(out, "modelo_PP_items_syntax.txt"))

# Registro ANTES de los dos ajustes nuevos. La base ya se había analizado:
# esto fija esta sensibilidad, pero no constituye un prerregistro del estudio.
protocol <- list(
  timestamp_UTC = format(Sys.time(), tz = "UTC"), N = nrow(data),
  status = "Exploratory sensitivity fixed before the two new full-sample fits; previous study results known",
  new_full_fits = 2L, estimators = c("MLR", "WLSMV"),
  ordinal_parameterization = "theta; underlying response residual variances fixed to 1 for ordinal PP outcomes",
  ordinal_variables = c(aps, pp), observed_levels = 1:5,
  same_P3_recoding = TRUE, controls_in_each_PP_equation = controls,
  measurement = "M2 SC/SBL APS block, including APS2 crossloading on SBL and APS13 on SC, unchanged; no PSB factor in item sensitivity",
  structural = "Six unchanged upstream equations; all six focal predictors plus CSR1 and seven P3 controls in each PP regression; three free PP residual covariances",
  primary_inference = "Six 3-df robust Wald tests of three item slopes jointly zero per predictor; Holm across the six predictors separately within estimator",
  secondary_inference = "18 pointwise robust item slopes and 95% CIs; Holm across the 18 focal item paths separately within estimator",
  heterogeneity = "MLR only: six 2-df Wald tests of equality of the three raw slopes, Holm across six; common observed 1-5 PP units; exploratory diagnostic",
  ordinal_heterogeneity_limit = "No equality test of WLSMV slopes as a substantive cross-item difference: latent ordinal outcome scales are arbitrary and can differ in residual dispersion",
  interpretation = "No causal effects; no item deletion; no model selection; no conventional AVE for observed PP outcomes; this does not repair convergent validity",
  prediction = "Delegated separate observed-composite ordinal-probit benchmark using exactly previous 5x5 country-stratified folds; not prediction from latent SEM",
  no_new_data = TRUE, no_preregistration_or_external_validation = TRUE,
  original_sha256 = as.list(before),
  code_sha256 = sha(file.path(root, "codigo/analizar_PP_por_item.R")))
if (!resume_exports) {
  write_json(protocol, file.path(out, "protocolo_antes_de_ajustes.json"),
             pretty = TRUE, auto_unbox = TRUE, digits = NA)
  writeLines("PROTOCOL_FROZEN", file.path(out, "protocolo_estado.txt"))
} else {
  stopifnot(file.exists(file.path(out, "protocolo_antes_de_ajustes.json")))
  frozen <- read_json(file.path(out, "protocolo_antes_de_ajustes.json"), simplifyVector = TRUE)
  stopifnot(frozen$N == 653L, frozen$new_full_fits == 2L,
            identical(frozen$original_sha256, as.list(before)))
  write_json(list(timestamp_UTC = format(Sys.time(), tz = "UTC"),
    reason = "Export checks repaired only: fixed zero theta entries for dummy endogenous observed variables are structural representation; residual variances are parameter-table diagonal covariances, often in psi. Exact list lookup avoids implied$cov partially matching cov.x for conditional WLSMV; its 26-response covariance is reconstructed with delta and joined to eight exogenous controls.",
    reuse_two_cached_fits = TRUE, additional_fits = 0L,
    frozen_protocol_unchanged = TRUE, current_code_sha256 = sha(file.path(root, "codigo/analizar_PP_por_item.R"))),
    file.path(out, "reparacion_exportacion_sin_nuevos_ajustes.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA)
}
csv <- function(x, name) write.csv(x, file.path(out, name), row.names = FALSE, na = "")
counts <- do.call(rbind, lapply(c(aps, pp), function(x) {
  data.frame(variable = x, category = 1:5,
             frequency = as.numeric(table(factor(data[[x]], levels = 1:5))))
}))
csv(counts, "frecuencias_categorias_agregadas.csv")

fits <- if (resume_exports) readRDS(file.path(private, "SEM_items_ambos_ajustes.rds")) else list()
if (!resume_exports) for (id in c("MLR", "WLSMV")) {
  warnings <- character()
  args <- list(model = syntax, data = data, estimator = id, meanstructure = TRUE,
               fixed.x = TRUE, control = list(iter.max = 2000))
  if (id == "WLSMV") {
    args$ordered <- c(aps, pp)
    args$parameterization <- "theta"
  }
  cat("FULL FIT", id, "\n"); flush.console()
  fit <- withCallingHandlers(do.call(lavaan::sem, args), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  attr(fit, "captured_warnings") <- unique(warnings)
  fits[[id]] <- fit
  saveRDS(fit, file.path(private, paste0("SEM_items_", id, ".rds")))
}
saveRDS(fits, file.path(private, "SEM_items_ambos_ajustes.rds"))

# Exportación auditable: IC pointwise y pruebas conjuntas con covarianza robusta.
# Las cifras estandarizadas ordinales describen una escala latente diferente.
all_paths <- all_joints <- all_heterogeneity <- fit_rows <- diag_rows <- list()
loads <- residuals <- all_r2 <- list()
for (id in names(fits)) {
  f <- fits[[id]]
  pe <- parameterEstimates(f, ci = TRUE)
  ss <- standardizedSolution(f, type = "std.all", ci = TRUE)
  csv(pe, paste0(id, "_todos_parametros_IC.csv"))
  csv(ss, paste0(id, "_estandarizados_IC.csv"))
  mat <- lavInspect(f, "est")
  imp <- lavInspect(f, "implied")
  conditional_min <- delta_identity_difference <- NA_real_
  # conditional.x=TRUE en WLSMV devuelve res.cov y cov.x, NO un campo cov.
  # El acceso exacto impide que $cov seleccione parcialmente sólo cov.x.
  implied <- imp[["cov", exact = TRUE]]
  if (is.null(implied)) {
    residual_cov <- imp[["res.cov", exact = TRUE]]
    slopes <- imp[["res.slopes", exact = TRUE]]
    xcov <- imp[["cov.x", exact = TRUE]]
    stopifnot(!is.null(residual_cov), !is.null(slopes), !is.null(xcov))
    U <- solve(diag(nrow(mat$beta)) - mat$beta)
    conditional_unscaled <- mat$lambda %*% U %*% mat$psi %*% t(U) %*% t(mat$lambda) + mat$theta
    D <- diag(as.numeric(mat$delta))
    delta_identity_difference <- max(abs(D %*% conditional_unscaled %*% D - residual_cov))
    stopifnot(delta_identity_difference < 1e-8)
    marginal_y <- residual_cov + slopes %*% xcov %*% t(slopes)
    cross_yx <- slopes %*% xcov
    implied <- rbind(cbind(marginal_y, cross_yx), cbind(t(cross_yx), xcov))
    conditional_min <- min(eigen(residual_cov, symmetric = TRUE, only.values = TRUE)$values)
  }
  vc <- lavInspect(f, "vcov")
  ev <- function(x) min(eigen((x + t(x))/2, symmetric = TRUE, only.values = TRUE)$values)
  # Las variables endógenas observadas se representan como dummy eta: su
  # theta puede estar fijada en cero y su varianza residual vive en psi.
  # Se comprueban las varianzas reales desde la tabla y theta semidefinida.
  variances <- pe$est[pe$op == "~~" & pe$lhs == pe$rhs]
  free <- parTable(f)$free > 0L
  all_finite <- all(is.finite(pe$est)) && all(is.finite(pe$se)) &&
    all(is.finite(pe$ci.lower)) && all(is.finite(pe$ci.upper))
  diag_rows[[id]] <- data.frame(estimator = id, N = lavInspect(f, "nobs"),
    converged = lavInspect(f, "converged"), post_check = lavInspect(f, "post.check"),
    min_variance = min(variances), min_theta_eigen = ev(mat$theta),
    min_psi_eigen = ev(mat$psi), min_implied_cov_eigen = ev(implied),
    implied_cov_dimension = nrow(implied),
    implied_cov_scope = "Marginal joint covariance of responses and exogenous controls, ordinal response scale where relevant",
    min_conditional_response_cov_eigen = conditional_min,
    max_delta_identity_difference = delta_identity_difference,
    min_vcov_eigen = ev(vc), all_parameter_SE_CI_finite = all_finite,
    warnings = paste(attr(f, "captured_warnings"), collapse = " | "))
  stopifnot(lavInspect(f, "converged"), lavInspect(f, "post.check"),
            all_finite, min(variances) > 0, ev(mat$theta) >= -1e-10,
            ev(mat$psi) > 0, ev(implied) > 0, ev(vc) > 0)
  fm <- fitMeasures(f, c("chisq", "df", "pvalue", "chisq.scaled", "df.scaled", "pvalue.scaled",
                         "cfi", "tli", "rmsea", "cfi.scaled", "tli.scaled", "rmsea.scaled",
                         "cfi.robust", "tli.robust", "rmsea.robust", "srmr", "npar"))
  fit_rows[[id]] <- data.frame(estimator = id, N = nrow(data), as.list(fm),
    note = "Describe each fixed representation; do not select estimator from fit or compare MLR/WLSMV indices as same-scale improvement", check.names = FALSE)
  foc <- pe[pe$op == "~" & pe$lhs %in% pp & pe$rhs %in% focal, ]
  stopifnot(nrow(foc) == 18L)
  matched <- ss[match(paste(foc$lhs, foc$op, foc$rhs), paste(ss$lhs, ss$op, ss$rhs)), ]
  foc$estimator <- id
  foc$beta_std_all <- matched$est.std
  foc$beta_CI_lower <- matched$ci.lower
  foc$beta_CI_upper <- matched$ci.upper
  foc$p_Holm_18 <- p.adjust(foc$pvalue, "holm")
  foc$CI_note <- "Pointwise 95% CI, not multiplicity-adjusted; p_Holm_18 is separate"
  foc$scale_note <- if (id == "MLR") "Observed PP 1-5 units" else "Probit underlying PP response; theta residual variance fixed to 1; cross-item slopes not substantively common scale"
  all_paths[[id]] <- foc
  joint <- lapply(focal, function(predictor) {
    con <- paste(vapply(pp, function(item) paste0(label(item, predictor), " == 0"), character(1)), collapse = "\n")
    w <- lavTestWald(f, constraints = con)
    data.frame(estimator = id, predictor = predictor, statistic = w$stat, df = w$df,
               pvalue = w$p.value, constraints = gsub("\n", "; ", con, fixed = TRUE),
               covariance = "Estimator-specific robust parameter covariance")
  })
  joint <- do.call(rbind, joint)
  joint$p_Holm_6 <- p.adjust(joint$pvalue, "holm")
  all_joints[[id]] <- joint
  if (id == "MLR") {
    heterogeneity <- lapply(focal, function(predictor) {
      con <- paste(paste0(label(pp[1], predictor), " == ", label(pp[2:3], predictor)), collapse = "\n")
      w <- lavTestWald(f, constraints = con)
      data.frame(estimator = id, predictor = predictor, statistic = w$stat, df = w$df,
                 pvalue = w$p.value, constraints = gsub("\n", "; ", con, fixed = TRUE),
                 scale_note = "Equality on common observed PP 1-5 units; exploratory conditional association heterogeneity")
    })
    heterogeneity <- do.call(rbind, heterogeneity)
    heterogeneity$p_Holm_6 <- p.adjust(heterogeneity$pvalue, "holm")
    all_heterogeneity[[id]] <- heterogeneity
  }
  load <- ss[ss$op == "=~", ]; load$estimator <- id; loads[[id]] <- load
  resid <- pe[pe$op == "~~" & pe$lhs %in% c(pp, "SC", "SBL") &
                pe$rhs %in% c(pp, "SC", "SBL"), ]
  resid$estimator <- id
  resid$std_all <- ss$est.std[match(paste(resid$lhs, resid$op, resid$rhs), paste(ss$lhs, ss$op, ss$rhs))]
  residuals[[id]] <- resid
  rs <- lavInspect(f, "rsquare")
  all_r2[[id]] <- data.frame(estimator = id, variable = names(rs), R_squared = as.numeric(rs),
                              scale_note = if (id == "MLR") "Observed outcomes/indicators" else "Underlying ordinal response where relevant")
}
bind <- function(x) do.call(rbind, x)
csv(bind(all_paths), "rutas_focales_18_por_estimador.csv")
csv(bind(all_joints), "Wald_conjuntos_6_Holm.csv")
csv(bind(all_heterogeneity), "Wald_heterogeneidad_MLR_6_Holm.csv")
csv(bind(fit_rows), "indices_ajuste_descriptivos.csv")
csv(bind(diag_rows), "diagnosticos_ajustes.csv")
csv(bind(loads), "cargas_APS_ambos_estimadores.csv")
csv(bind(residuals), "residuales_PP_SC_SBL.csv")
csv(bind(all_r2), "R2_ambos_estimadores.csv")
after <- vapply(protected, sha, character(1))
stopifnot(identical(before, after))
write_json(list(status = "COMPLETED", timestamp_UTC = format(Sys.time(), tz = "UTC"),
  N = nrow(data), full_fits = 2L, all_admissible = TRUE,
  all_SE_CI_finite = TRUE, original_files_unchanged = TRUE,
  protected_sha256 = as.list(after),
  note = "Item outcomes are observed responses and have no AVE; CFA convergence is assessed separately in original factor structure"),
  file.path(out, "auditoria_y_finalizacion.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
writeLines("COMPLETED", state)
cat("COMPLETED: two fixed full-sample SEM item sensitivities\n")
