# Verificación algebraica de resultados ya guardados: NO estima modelos.
project <- "/Users/luisbenites/Documents/Codex/Paper con MA - 1/reproduccion-integral-V3-2026-10-10/proyecto"
root <- file.path(project, "sem-items-ordinal-2026-10-08/items")
.libPaths(c(file.path(project, "revision-revista-2026-10-04/auditoria-codigos/R-library"),
            "/Users/luisbenites/Library/R/arm64/4.6/library", .libPaths()))
suppressPackageStartupMessages({library(lavaan); library(jsonlite); library(digest)})
out <- file.path(root, "resultados")
fits <- readRDS(file.path(root, "privados/SEM_items_ambos_ajustes.rds"))
data <- readRDS(file.path(root, "privados/datos_P3_para_items_y_CV.rds"))
joints <- read.csv(file.path(out, "Wald_conjuntos_6_Holm.csv"))
heterogeneity <- read.csv(file.path(out, "Wald_heterogeneidad_MLR_6_Holm.csv"))
paths <- read.csv(file.path(out, "rutas_focales_18_por_estimador.csv"))
records <- list(); i <- 0L
for (id in names(fits)) {
  fit <- fits[[id]]; table <- parTable(fit); coefficients <- coef(fit)
  covariance <- lavInspect(fit, "vcov")
  rows <- joints[joints$estimator == id, ]
  for (predictor in rows$predictor) {
    labels <- paste0("item", 1:3, "_", predictor)
    free_indices <- table$free[match(labels, table$label)]
    stopifnot(length(free_indices) == 3L, all(free_indices > 0L))
    b <- coefficients[free_indices]; V <- covariance[free_indices, free_indices]
    statistic <- as.numeric(t(b) %*% solve(V, b))
    expected <- rows$statistic[rows$predictor == predictor]
    i <- i + 1L
    records[[i]] <- data.frame(estimator = id, predictor = predictor,
      test = "Joint association", manually_computed = statistic,
      lavaan_export = expected, absolute_difference = abs(statistic - expected))
    if (id == "MLR") {
      contrast <- matrix(c(1, -1, 0, 1, 0, -1), nrow = 2, byrow = TRUE)
      cb <- contrast %*% b; CV <- contrast %*% V %*% t(contrast)
      statistic <- as.numeric(t(cb) %*% solve(CV, cb))
      expected <- heterogeneity$statistic[heterogeneity$predictor == predictor]
      i <- i + 1L
      records[[i]] <- data.frame(estimator = id, predictor = predictor,
        test = "MLR slope heterogeneity", manually_computed = statistic,
        lavaan_export = expected, absolute_difference = abs(statistic - expected))
    }
  }
  stopifnot(max(abs(p.adjust(rows$pvalue, "holm") - rows$p_Holm_6)) < 1e-12)
  pp <- paths[paths$estimator == id, ]
  stopifnot(nrow(pp) == 18L,
            max(abs(p.adjust(pp$pvalue, "holm") - pp$p_Holm_18)) < 1e-12)
}
audit <- do.call(rbind, records)
stopifnot(max(audit$absolute_difference) < 1e-8)
write.csv(audit, file.path(out, "auditoria_Wald_calculo_algebraico.csv"), row.names = FALSE)
source <- readRDS(file.path(project, "revision-integral-2026-10-04/sensitivities/P3_country_partial_demographics.rds"))
source_matrix <- lavInspect(source, "data")
identical_data <- max(abs(as.matrix(data[, colnames(source_matrix)]) - source_matrix)) == 0
stopifnot(identical_data, nrow(data) == 653L)
write_json(list(all_pass = TRUE, no_new_model_fits = TRUE, N = nrow(data),
  exact_P3_data = identical_data, manual_Wald_tests = nrow(audit),
  maximum_manual_Wald_difference = max(audit$absolute_difference),
  Holm_families_verified = TRUE,
  note = "This is an algebraic audit by the analysis author; independent review should inspect the same cached fits"),
  file.path(out, "auditoria_Wald_y_P3.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA)
cat("AUDIT_OK; 18 algebraic Wald verifications; no new fits\n")
