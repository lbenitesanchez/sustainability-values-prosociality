# Observación 1: identificación y calidad numérica de los 28 ajustes YA archivados.
# No estima modelos, no modifica datos ni archivos fuente. Produce las tablas
# identificacion_28_ajustes.csv y restricciones_7_especificaciones.csv.
options(stringsAsFactors = FALSE)
Sys.setenv(VECLIB_MAXIMUM_THREADS = '1', OPENBLAS_NUM_THREADS = '1', OMP_NUM_THREADS = '1')
project <- '/Users/luisbenites/Documents/Codex/Paper con MA - 1/reproduccion-integral-V3-2026-10-10/proyecto'
dest <- file.path(project, 'revision-observaciones-2026-10-09/diagnosticos')
source_dir <- file.path(project, 'sem-alternativas-2026-10-08/resultados')
.libPaths(c(file.path(dest, 'R-library'),
            file.path(project, 'revision-revista-2026-10-04/auditoria-codigos/R-library'),
            '/Users/luisbenites/Library/R/arm64/4.6/library', .libPaths()))
suppressPackageStartupMessages({library(lavaan); library(jsonlite); library(digest)})
dir.create(file.path(dest, 'resultados'), recursive = TRUE, showWarnings = FALSE)
write_out <- function(x, name) write.csv(x, file.path(dest, 'resultados', name), row.names = FALSE, na = '')
ids <- c('M0_original', 'M1_ASP_unico', 'M2_puentes2', 'M3_puentes4',
         'M4_S1_refSC', 'M5_S1_refSBL', 'M6_ESEM_target')
sc <- c('APS1','APS2','APS3','APS4','APS5','APS6','APS7','APS9','APS10','APS14')
sbl <- c('APS11','APS12','APS13','APS15','APS16','APS18','APS19','APS20','APS21')
items <- c(sc, sbl, 'PP1', 'PP2', 'PP3')

# Rango numérico: tolerancia relativa 1e-8, documentada en cada resultado.
# Un rango completo es evidencia LOCAL en el punto estimado, no una demostración
# de identificación global, ni de significado sustantivo, ni de buen ajuste.
matrix_diag <- function(x, symmetric = FALSE) {
  x <- as.matrix(x)
  if (!all(is.finite(x)) || !length(x)) return(list(rank = NA_integer_, dim = ncol(x),
      min = NA_real_, max = NA_real_, tol = NA_real_, negative = NA_integer_, cond = NA_real_))
  if (symmetric) {
    ev <- eigen((x + t(x))/2, symmetric = TRUE, only.values = TRUE)$values
    tol <- 1e-8 * max(1, max(abs(ev)))
    list(rank = sum(abs(ev) > tol), dim = ncol(x), min = min(ev), max = max(ev),
         tol = tol, negative = sum(ev < -tol), cond = if (all(ev > tol)) max(ev)/min(ev) else NA_real_)
  } else {
    sv <- svd(x, nu = 0, nv = 0)$d
    tol <- 1e-8 * max(1, max(sv))
    list(rank = sum(sv > tol), dim = ncol(x), min = min(sv), max = max(sv),
         tol = tol, negative = NA_integer_, cond = if (all(sv > tol)) max(sv)/min(sv) else NA_real_)
  }
}
one_matrix <- function(x) if (is.list(x)) x[[1L]] else x

# Las restricciones activas forman J. Z es una base ortonormal de null(J).
# Las dos restricciones rotacionales ESEM hacen que el rango bruto pueda ser
# menor que el número de coordenadas sin que eso pruebe no identificación.
tangent_basis <- function(fit, p) {
  J <- fit@Model@con.jac
  if (!length(J) || !nrow(J)) return(list(J = matrix(0, 0, p), Z = diag(p), rank = 0L))
  stopifnot(ncol(J) == p, all(is.finite(J)))
  q <- qr(t(J), tol = 1e-10)
  Q <- qr.Q(q, complete = TRUE)
  Z <- if (q$rank < p) Q[, seq.int(q$rank + 1L, p), drop = FALSE] else matrix(0, p, 0)
  stopifnot(max(abs(J %*% Z), 0) < 1e-6)
  list(J = J, Z = Z, rank = q$rank)
}

rows <- list(); constraint_rows <- list()
for (id in ids) for (kind in c('CFA', 'SEM')) for (estimator in c('MLR', 'WLSMV')) {
  key <- paste(id, kind, estimator, sep = '_')
  file <- file.path(source_dir, paste0(key, '.rds'))
  f <- readRDS(file)
  stopifnot(inherits(f, 'lavaan'), sum(lavInspect(f, 'nobs')) == 653L)
  pt <- parameterTable(f)
  I <- as.matrix(lavInspect(f, 'information'))
  H <- as.matrix(lavInspect(f, 'hessian'))
  D <- as.matrix(one_matrix(lavInspect(f, 'delta')))
  V <- as.matrix(lavInspect(f, 'vcov'))
  basis <- tangent_basis(f, nrow(I)); Z <- basis$Z
  stopifnot(ncol(D) == nrow(I), identical(dim(H), dim(I)))
  rawI <- matrix_diag(I, TRUE); rawH <- matrix_diag(H, TRUE); rawD <- matrix_diag(D)
  tI <- matrix_diag(crossprod(Z, I %*% Z), TRUE)
  tH <- matrix_diag(crossprod(Z, H %*% Z), TRUE)
  tD <- matrix_diag(D %*% Z)
  # Covarianza sándwich/robusta: su finitud no sustituye el diagnóstico Hessiano.
  free_rows <- pt$free > 0 & !pt$op %in% c(':=', '==', '>', '<')
  freeSE <- pt$se[free_rows]
  vr <- pt$op == '~~' & pt$lhs == pt$rhs
  freevar <- vr & pt$free > 0
  fixed0 <- vr & pt$free == 0 & abs(pt$est) < 1e-12
  # Sólo los 22 indicadores tienen errores estocásticos de medición aquí.
  # Los ceros adicionales en theta del SEM corresponden a predictores observados
  # representados por lavaan como variables auxiliares; no son varianzas libres negativas.
  T <- as.matrix(one_matrix(lavInspect(f, 'theta')))
  idx <- match(items, rownames(T)); stopifnot(!anyNA(idx))
  td <- matrix_diag(T[idx, idx, drop = FALSE], TRUE)
  mats <- lavInspect(f, 'est')
  psi <- if (!is.null(mats$psi)) mats$psi else mats[[1L]]$psi
  pd <- matrix_diag(psi, TRUE)
  fm <- suppressWarnings(tryCatch(fitMeasures(f), error = function(e) numeric()))
  getfm <- function(n) if (n %in% names(fm)) unname(fm[n]) else NA_real_
  tests <- f@test
  robust <- if (length(tests) > 1L) tests[[2L]] else NULL
  robust_stat <- if (!is.null(robust$stat)) robust$stat else NA_real_
  robust_p <- if (!is.null(robust$pvalue)) robust$pvalue else NA_real_
  proper_variance <- !any(pt$est[freevar] < -1e-8) && td$negative == 0 && pd$negative == 0
  finite_se <- all(is.finite(freeSE))
  local_rank <- tD$rank == ncol(Z) && tI$rank == ncol(Z)
  tangent_H_positive <- tH$rank == ncol(Z) && tH$negative == 0 && tH$min > tH$tol
  opts <- f@Options
  scale <- if (kind == 'CFA') 'std.lv=TRUE; factor variances fixed to 1' else
    paste(paste(pt$lhs[pt$op == '=~' & pt$free == 0 & abs(pt$est - 1) < 1e-10],
                pt$rhs[pt$op == '=~' & pt$free == 0 & abs(pt$est - 1) < 1e-10], sep = ':'), collapse = '; ')
  if (id == 'M6_ESEM_target') scale <- if (kind == 'CFA')
    'EFA and PSB unit variances; two active rotation restrictions' else
    'EFA unit disturbance variances before rotation; PSB marker PP1; two active rotation restrictions'
  summary_status <- if (!proper_variance) 'Improper variance/covariance diagnostic' else
    if (!finite_se) 'Non-finite standard errors' else if (!local_rank) 'Local projected rank deficient' else
      if (!tangent_H_positive) 'Projected observed Hessian not positive definite' else
        if (!is.finite(robust_stat) || !is.finite(robust_p)) 'Robust global test unavailable' else
          'Converged; local numerical diagnostics pass'
  rows[[key]] <- data.frame(model = id, kind = kind, estimator = estimator, N = 653,
    scale = scale, crossloading_policy = if (id == 'M6_ESEM_target') 'All 38 APS loadings estimated; soft target zeros' else
      'Unspecified crossloadings fixed to zero; see specification table',
    rotation = if (id == 'M6_ESEM_target') 'Oblique target (lavaan pst); 20 starts in archived fit' else 'None',
    df = tests[[1L]]$df, coordinates = nrow(I), active_constraints = basis$rank, free_tangent_dimension = ncol(Z),
    converged = lavInspect(f, 'converged'), native_post_check = suppressWarnings(lavInspect(f, 'post.check')),
    finite_free_SE = finite_se, zero_free_SE_count = sum(freeSE == 0, na.rm = TRUE), vcov_finite = all(is.finite(V)),
    raw_information_rank = rawI$rank, raw_hessian_rank = rawH$rank, raw_Delta_rank = rawD$rank,
    tangent_information_rank = tI$rank, tangent_information_min_eigen = tI$min, tangent_information_tol = tI$tol,
    tangent_Hessian_rank = tH$rank, tangent_Hessian_min_eigen = tH$min, tangent_Hessian_negative_eigen_count = tH$negative,
    tangent_Hessian_tol = tH$tol, tangent_Delta_rank = tD$rank, tangent_Delta_min_singular = tD$min,
    tangent_Delta_tol = tD$tol, local_rank_pass = local_rank, tangent_Hessian_positive = tangent_H_positive,
    negative_free_variances = sum(pt$est[freevar] < -1e-8), fixed_zero_variance_count = sum(fixed0),
    theta_indicator_min_eigen = td$min, psi_min_eigen = pd$min, proper_variance = proper_variance,
    robust_test = if (!is.null(robust$test)) robust$test else '', robust_chisq = robust_stat, robust_p = robust_p,
    cfi_robust = getfm('cfi.robust'), rmsea_robust = getfm('rmsea.robust'), srmr = getfm('srmr'),
    warnings_archived = paste(attr(f, 'captured_warnings'), collapse = ' | '), status = summary_status,
    source_sha256 = digest(file, algo = 'sha256', file = TRUE))
  cat('AUDITED', key, summary_status, '\n'); flush.console()
}
write_out(do.call(rbind, rows), 'identificacion_28_ajustes.csv')

# Tabla descriptiva de las restricciones. Deriva los marcadores del objeto SEM MLR
# y detalla restricciones explícitas del constructor original sin reestimarlo.
cross <- c('SC/SBL primary loading only; PSB PP1-PP3; all other crossloadings zero',
  'Single SA factor for all 19 APS; separate PSB; other crossloadings zero',
  'M0 plus APS2 on SBL and APS13 on SC; all other crossloadings zero',
  'M0 plus APS2 and APS14 on SBL; APS13 and APS18 on SC; other crossloadings zero',
  'SCref on all 19 APS; SBLspec on 9 SBL items; PSB separate; others zero',
  'SBLref on all 19 APS; SCspec on 10 SC items; PSB separate; others zero',
  'Both APS factors load on all 19 APS; PSB CFA PP1-PP3; APS target zeros are soft, not fixed')
for (i in seq_along(ids)) {
  id <- ids[i]; f <- readRDS(file.path(source_dir, paste0(id, '_SEM_MLR.rds'))); p <- parameterTable(f)
  marker <- p[p$op == '=~' & p$free == 0 & abs(p$est - 1) < 1e-10, ]
  constraint_rows[[id]] <- data.frame(model = id, indicators = '19 APS plus 3 PP; all retained',
    CFA_scaling = 'std.lv=TRUE: factor variances fixed to 1',
    SEM_scaling = if (i == 7) 'EFA unit disturbance variances before rotation; PSB loading PP1 fixed to 1' else
      paste(paste0(marker$lhs, ': loading ', marker$rhs, '=1'), collapse = '; '),
    crossloadings = cross[i],
    covariance_constraint = if (i %in% 5:6) 'Reference-specific covariance zero in CFA; disturbance covariance zero in SEM' else
      if (i == 2) 'SA and PSB freely covary in CFA; PSB regressed on SA in SEM' else
        'APS factors freely covary in CFA; APS disturbance covariance estimated in SEM',
    rotation_constraints = if (i == 7) 'Oblique target; 2 active constraints; factor sign/label orientation is rotational' else 'None',
    WLSMV_identification = 'Theta parameterization: ordinal response residual variances fixed to 1; four thresholds/item; fixed.x TRUE',
    MLR_identification = 'Continuous working-response scale; meanstructure TRUE; fixed.x TRUE',
    scope = 'Local numerical evidence only; not proof of global identification or invariance')
}
write_out(do.call(rbind, constraint_rows), 'restricciones_7_especificaciones.csv')
writeLines(capture.output(sessionInfo()), file.path(dest, 'resultados', 'sessionInfo_identificacion.txt'))
cat('COMPLETED identification audit: 28 archived objects; no model fits\n')
