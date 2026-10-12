# Evidencia para la revisión neutral, observaciones 2, 3, 6 y 7.
# EJECUCIÓN: /usr/local/bin/Rscript 01_exportar_evidencia_existente.R
# Este archivo SOLO lee ajustes lavaan YA ejecutados y agregados archivados.
# No llama a cfa(), sem(), update(), ni realiza selección o estimación nueva.
# No exporta datos individuales, identificadores, particiones ni archivos RDS.
# Los nombres de cada CSV explicitan la tabla/resultado que pueden alimentar.

proyecto <- '/Users/luisbenites/Documents/Codex/Paper con MA - 1/reproduccion-integral-V3-2026-10-10/proyecto'
salida <- file.path(proyecto, 'revision-observaciones-2026-10-09/evidencia/observaciones_2_3_6_7')
dir.create(salida, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(file.path(proyecto, 'revision-revista-2026-10-04/auditoria-codigos/R-library'), .libPaths()))
suppressPackageStartupMessages(library(lavaan))
suppressPackageStartupMessages(library(semTools))
library(jsonlite)
library(digest)
base <- file.path(proyecto, 'sem-alternativas-2026-10-08/resultados')
sens <- file.path(proyecto, 'revision-3-revista-colombiana-psicologia-2026-10-08/analisis_sensibilidad')
inv <- file.path(proyecto, 'sem-multigrupo-2026-10-08/invariancia/resultados')
fuentes <- character()
hashes_antes <- character()
leer_rds <- function(ruta) {
  fuentes <<- c(fuentes, ruta)
  hashes_antes <<- c(hashes_antes, digest(file = ruta, algo = 'sha256'))
  readRDS(ruta)
}
leer_csv <- function(ruta) {
  fuentes <<- c(fuentes, ruta)
  hashes_antes <<- c(hashes_antes, digest(file = ruta, algo = 'sha256'))
  read.csv(ruta, check.names = FALSE, stringsAsFactors = FALSE)
}
guardar <- function(x, nombre) write.csv(x, file.path(salida, nombre), row.names = FALSE, na = '')
pegar <- function(x) do.call(rbind, x)

# TABLAS 1 / S1 / S4: ajuste, tamaño, grados de libertad y complejidad.
# RMSEA: intervalo del 90%, configurado explícitamente para reproducir los
# intervalos archivados. No confundirlo con los IC puntuales del 95% de parámetros.
campos <- c('chisq', 'chisq.scaled', 'df', 'npar', 'cfi.robust', 'tli.robust',
            'rmsea.robust', 'rmsea.ci.lower.robust', 'rmsea.ci.upper.robust',
            'cfi.scaled', 'tli.scaled', 'rmsea.scaled',
            'rmsea.ci.lower.scaled', 'rmsea.ci.upper.scaled', 'srmr', 'aic', 'bic')
ajuste <- function(f, modelo, tipo, estimador, escenario = 'all') {
  fm <- suppressWarnings(tryCatch(fitMeasures(f, fm.args = list(rmsea.ci.level = .90)),
                                 error = function(e) numeric()))
  valores <- vapply(campos, function(z) if (z %in% names(fm)) unname(fm[z]) else NA_real_, numeric(1))
  data.frame(model = modelo, kind = tipo, estimator = estimador, scenario = escenario,
             N = sum(lavInspect(f, 'nobs')), converged = lavInspect(f, 'converged'),
             rmsea_ci_level = .90, as.list(valores), check.names = FALSE)
}
modelos <- c('M0_original', 'M1_ASP_unico', 'M2_puentes2', 'M3_puentes4',
             'M4_S1_refSC', 'M5_S1_refSBL', 'M6_ESEM_target')
ajustes <- cargas <- correlaciones <- list()
objetos <- list()
for (m in modelos) for (tipo in c('CFA', 'SEM')) for (e in c('MLR', 'WLSMV')) {
  key <- paste(m, tipo, e, sep = '_')
  f <- leer_rds(file.path(base, paste0(key, '.rds')))
  objetos[[key]] <- f
  ajustes[[key]] <- ajuste(f, m, tipo, e)
  # Cargas de CFA, con IC puntuales 95% en la escala estandarizada std.all.
  # MLR usa categorías numéricas; WLSMV usa respuestas continuas subyacentes.
  if (tipo == 'CFA') {
    st <- suppressWarnings(standardizedSolution(f, type = 'std.all', ci = TRUE, level = .95))
    x <- st[st$op == '=~', c('lhs', 'rhs', 'est.std', 'se', 'ci.lower', 'ci.upper')]
    x$model <- m; x$estimator <- e; x$N <- sum(lavInspect(f, 'nobs'))
    x$ci_level <- .95; x$ci_scope <- 'Pointwise; no model-selection adjustment'
    x$scale <- if (e == 'MLR') 'Numeric category scores' else 'Underlying ordinal responses'
    cargas[[key]] <- x
    lv <- lavNames(f, 'lv')
    x <- st[st$op == '~~' & st$lhs != st$rhs & st$lhs %in% lv & st$rhs %in% lv,
            c('lhs', 'rhs', 'est.std', 'se', 'ci.lower', 'ci.upper')]
    if (nrow(x)) {
      x$model <- m; x$kind <- tipo; x$estimator <- e; x$N <- sum(lavInspect(f, 'nobs'))
      x$ci_level <- .95; x$definition <- 'CFA latent factor correlation; pointwise CI'
      # Una correlación fijada por definición no tiene incertidumbre estimada:
      # no presentar [0,0] como intervalo que prueba independencia empírica.
      pt <- parameterTable(f); ix <- match(paste(x$lhs,'~~',x$rhs),paste(pt$lhs,pt$op,pt$rhs))
      fijo <- !is.na(ix) & pt$free[ix] == 0
      x$se[fijo] <- NA_real_; x$ci.lower[fijo] <- NA_real_; x$ci.upper[fijo] <- NA_real_
      x$ci_level[fijo] <- NA_real_; x$definition[fijo] <- 'Fixed correlation constraint; no sampling CI or empirical independence claim'
      correlaciones[[key]] <- x
    }
  } else {
    # SEM: cor.lv es correlación TOTAL implicada, no la covarianza residual ~~.
    # No trasladar el IC de una covarianza residual a una correlación total.
    cc <- lavInspect(f, 'cor.lv'); inds <- which(upper.tri(cc), arr.ind = TRUE)
    if (nrow(inds)) {
      x <- data.frame(lhs = rownames(cc)[inds[,1]], rhs = colnames(cc)[inds[,2]],
                      est.std = cc[inds], se = NA_real_, ci.lower = NA_real_, ci.upper = NA_real_,
                      model = m, kind = tipo, estimator = e, N = sum(lavInspect(f, 'nobs')),
                      ci_level = NA_real_, definition = 'SEM total model-implied correlation; CI not archived')
      correlaciones[[key]] <- x
    }
  }
}
guardar(pegar(ajustes), 'T1_S1_S4_ajuste_N_df_npar_IC90.csv')
guardar(pegar(cargas), 'S_medicion_cargas_CFA_IC95_por_estimador.csv')
guardar(pegar(correlaciones), 'S_medicion_correlaciones_CFA_IC95_SEM_total_sin_IC.csv')

# TABLAS S3 / S5: AVE y fiabilidad con las definiciones de la versión neutral.
# AVE no es la media de las cargas cuadradas cuando los indicadores difieren
# en varianza. semTools usa la razón de varianza común/varianza pertinente.
# omega ordinal en categorías observadas y omega latente son magnitudes distintas.
# Reutilizamos los cálculos ya regenerados y auditados por el código base de
# esta misma revisión, evitando repetir integración ordinal costosa. Este
# archivo regenera la presentación y trazabilidad, no la estimación de modelos.
audit_base <- file.path(proyecto,'revision-observaciones-2026-10-09/resultados_base')
r3 <- leer_csv(file.path(audit_base,'S3_AVE_omega_recalculadas.csv'))
r3$N <- 653L; r3$ci_lower <- NA_real_; r3$ci_upper <- NA_real_
r3$ci_note <- 'Point intervals for AVE/omega not archived; no new bootstrap'
guardar(r3,'S3_AVE_fiabilidad_fuente_regenerada_auditada.csv')
r5 <- leer_csv(file.path(audit_base,'02_S5_medicion_numerica.csv'))
r5$AVE_ci_lower <- NA_real_; r5$AVE_ci_upper <- NA_real_
r5$omega_ci_lower <- NA_real_; r5$omega_ci_upper <- NA_real_
r5$uncertainty_note <- 'Reliability/AVE CIs not archived; no new bootstrap was performed'
r5$definition <- 'PSB composite reliability; SC/SBL conventional AVE unavailable with cross-loadings'
guardar(r5,'S3_S5_AVE_fiabilidad_escalas_sin_IC_inventados.csv')
especifica <- list()
for (m in c('M4_S1_refSC','M5_S1_refSBL')) for(e in c('MLR','WLSMV')) {
  f <- objetos[[paste(m,'CFA',e,sep='_')]]; est <- lavInspect(f,'est')
  dominio <- if(m=='M4_S1_refSC') 'SBLspec' else 'SCspec'
  items <- if(dominio=='SBLspec') paste0('APS',c(11,12,13,15,16,18,19,20,21)) else paste0('APS',c(1,2,3,4,5,6,7,9,10,14))
  Sigma <- est$lambda %*% est$psi %*% t(est$lambda) + est$theta
  valor <- sum(est$lambda[items,dominio])^2 * est$psi[dominio,dominio] / sum(Sigma[items,items])
  especifica[[paste(m,e)]] <- data.frame(model=m,estimator=e,specific_domain=dominio,N=sum(lavInspect(f,'nobs')),
    omega_HS=valor, scale=if(e=='MLR') 'Numeric category scores' else 'Underlying ordinal responses',
    definition='Reliability of orthogonal domain-specific component; not full-domain reliability',
    ci_lower=NA_real_,ci_upper=NA_real_)
}
guardar(pegar(especifica), 'S3_Smenos1_fiabilidad_componente_especifico.csv')

# OBSERVACIÓN 2 / nuevas tablas S10+: sensibilidad por banderas de calidad.
# Se reabren ajustes guardados; NO se aplican exclusiones ni se vuelven a ajustar.
escenarios <- c('all','without_duplicate_copies','without_constant22','without_both')
sens_ajuste <- list()
for(s in escenarios) for(k in c('CFA_MLR','CFA_WLSMV','SEM_MLR','SEM_WLSMV','ITEMS_MLR','ITEMS_WLSMV')) {
  f <- leer_rds(file.path(sens,'privados',paste0(s,'__',k,'.rds')))
  a <- strsplit(k,'_')[[1]]
  sens_ajuste[[paste(s,k)]] <- ajuste(f,'M2_puentes2',a[1],a[2],s)
  # Prueba de identidad exacta para referencias full sample CFA/SEM.
  if(s=='all' && a[1] %in% c('CFA','SEM')) {
    ref <- objetos[[paste('M2_puentes2',a[1],a[2],sep='_')]]
    stopifnot(identical(coef(f),coef(ref)), identical(lavInspect(f,'nobs'),lavInspect(ref,'nobs')))
  }
}
guardar(pegar(sens_ajuste), 'S10_sensibilidad_ajuste_N_df_npar_IC90.csv')
copias <- c('escenarios_muestra.csv','modelos_centrales/medicion_CFA_sensibilidad.csv',
 'modelos_centrales/rutas_focales_sensibilidad.csv','modelos_centrales/Wald_conjuntos_6_Holm.csv',
 'modelos_centrales/Wald_heterogeneidad_MLR_6_Holm.csv','auditoria_fuga_particiones_historicas.csv',
 'grouped_cv/resumen_modelos_por_item.csv','grouped_cv/resumen_diferencias_pareadas.csv',
 'grouped_cv/metricas_por_repeticion.csv','grouped_cv/tamanos_y_agrupacion_folds.csv',
 'tablas_para_manuscrito/S10_medicion_y_SEM_sensibilidad.csv',
 'tablas_para_manuscrito/S11_asociaciones_PP_sensibilidad.csv',
 'tablas_para_manuscrito/S12_incremento_predictivo_sensibilidad.csv')
for(p in copias) guardar(leer_csv(file.path(sens,p)), paste0('archivado_',basename(p)))

# OBSERVACIÓN 6 / invariancia M6: solo los resultados ya calculados con MLR.
# País e institución están confundidos; no sustenta efectos nacionales.
# Métodos robustos explícitos, nunca restar dos chi-cuadrados escalados.
x <- leer_csv(file.path(inv,'comparacion_invariancia.csv'))
x <- x[x$model %in% c('configural','metric','scalar'),]
x$N <- 653L; x$rmsea_ci_level <- .90
guardar(x,'S_invariancia_M6_MLR_ajuste_N_df_npar_IC90.csv')
x <- leer_csv(file.path(inv,'robust_LRT_explicit_and_delta.csv'))
guardar(x[!is.na(x[['Chisq diff']]),], 'S_invariancia_M6_MLR_LRT_metodos_explicitos.csv')
flujo <- data.frame(stage=c('Invited students','Initial survey responses','Archived records',
 'Complete analytic records','Full sample analysis','Without subsequent matching copies',
 'Without constant 22-item profiles','Without both flagged sets'),
 N=c(NA,NA,653,653,653,650,618,615),
 status=c('Unknown; not supplied','Unknown; not supplied',rep('Observed',6)),
 note=c('Response rate cannot be calculated','Prior exclusions/duplicates before archival unknown',
 'Convenience sample from four institutions','Complete analysis variables in archived dataset',
 'Primary reference; flags do not diagnose inattention','3 subsequent copies; identity unverified',
 '35 constant profiles; do not imply invalid responding','38 flagged records; sets do not overlap'))
guardar(flujo,'S_flujo_reclutamiento_conocido_y_desconocido.csv')

# Verificación final de fuentes: se incluye únicamente ruta+SHA, nunca contenido.
despues <- vapply(fuentes, function(p) digest(file=p,algo='sha256'),character(1))
stopifnot(identical(unname(hashes_antes),unname(despues)))
guardar(data.frame(path=fuentes,sha256_before=hashes_antes,sha256_after=despues),
        'manifest_fuentes_solo_lectura.csv')
write_json(list(status='EXPORTED_FROM_EXISTING_FITS_NO_NEW_ESTIMATION',
 lavaan_version=as.character(packageVersion('lavaan')),semTools_version=as.character(packageVersion('semTools')),
 original_sources_unchanged=TRUE,full_sample_M2_coefficients_identical=TRUE,
 RMSEA_ci_level=.90,parameter_ci_level=.95,
 files=sort(list.files(salida,pattern='\\.csv$')),
 limits=c('Intervals are pointwise, not model-selection adjusted',
 'Reliability/AVE CIs absent, not manufactured',
 'SEM total latent correlations have no archived CI',
 'No grouped B0/B1 benchmark is available',
 'Metric evidence mixed; scalar invariance unsupported',
 'Recruitment denominators unknown; country and institution confounded')),
 file.path(salida,'reporte_exportacion_R.json'),pretty=TRUE,auto_unbox=TRUE,na='null')
cat('PASS: aggregate evidence exported; original sources unchanged; no new fits.\n')
