# [LISO_PREVIO] / tabla provisional S20: reutilizar sensibilidad ya ejecutada.
# SOLO copia cifras agregadas auditadas; no cambia modelos/datos ni adopta
# orientacion editorial de la carpeta donde se archivaron esos calculos.
# Cuatro splits: una institucion fuera, sin indicadores de pais en predictores.
# El probit/OLS previo no es el modelo ridge de cinco categorias del script 01.
project <- Sys.getenv("PAPER_PROYECTO","/Users/luisbenites/Documents/Codex/Paper con MA - 1/reproduccion-integral-V3-2026-10-10/proyecto")
.libPaths(c(file.path(project,"revision-revista-2026-10-04/auditoria-codigos/R-library"),.libPaths()))
stopifnot(requireNamespace("jsonlite",quietly=TRUE))
src <- file.path(project,"revision-4-ARLA-2026-10-09/transportabilidad_institucional")
dst <- file.path(project,"revision-observaciones-2026-10-09/prediccion/LISO_previo_agregados")
dir.create(dst,recursive=TRUE,showWarnings=FALSE)
sha <- function(path)strsplit(system2("/usr/bin/shasum",c("-a","256",shQuote(path)),stdout=TRUE)," +")[[1]][1]
aud <- jsonlite::read_json(file.path(src,"auditoria_final.json"),simplifyVector=TRUE)
qpath <- file.path(project,"revision-3-revista-colombiana-psicologia-2026-10-08/analisis_sensibilidad/privados/datos_y_banderas_para_sensibilidades.rds")
rpath <- file.path(project,"revision-integral-2026-10-04/analysis/outputs/canonical-full-20261004-02/01_audit/analytic_data.rds")
stopifnot(aud$N==653,aud$n_sites==4,isTRUE(aud$all_models_scored),isTRUE(aud$all_numerically_valid),
          aud$identical_groups_crossing_splits==0,sha(qpath)==aud$sources_sha256$quality,sha(rpath)==aud$sources_sha256$archived)
names_to_copy <- c("contrastes_pareados_por_institucion.csv","tablas_para_manuscrito/S14_LISO_resumen_agrupado_numerico.csv",
  "metricas_y_diagnosticos_por_institucion.csv","resumen_errores_agrupados.csv","resumen_incremento_predictivo.csv",
  "particiones_y_duplicados.csv","frecuencias_train_test.csv","auditoria_final.json","auditoria_independiente.json")
provenance <- lapply(names_to_copy,function(name) {
  from <- file.path(src,name);to<-file.path(dst,basename(name))
  stopifnot(file.exists(from));file.copy(from,to,overwrite=TRUE);stopifnot(sha(from)==sha(to))
  list(file=basename(name),source=from,sha256=sha(to))
})
jsonlite::write_json(list(status="Identical audited aggregate reuse; no new estimation",N=653,sites=4,
  files=provenance,quality_sha256=sha(qpath),raw_sha256=sha(rpath),
  limitation="Descriptive four-institution sensitivity of previous probit/OLS; not external validation of new ridge model"),
  file.path(dst,"procedencia_reutilizacion_LISO.json"),auto_unbox=TRUE,pretty=TRUE)
cat("LISO: nueve archivos agregados auditados copiados; datos exactos verificados.\n")
