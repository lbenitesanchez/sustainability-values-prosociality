# MÓDULO 06 — TABLAS V3 DEFINITIVAS T3 Y S9.
# Función pública: generar_prediccion_agrupada().
# Lee resúmenes ya calculados con las mismas particiones agrupadas; no entrena.
# Se elige scenario==all: las exclusiones se presentan por separado en S14.
# B3=APS sin bloque de valores; B2=APS más cuatro valores (mismos controles).

# [T3] y [S9] Prediccion con grupos de registros identicos ----
# El archivo numerico procede de la evaluacion agrupada ya ejecutada.
# Cada media compara B2 y B3 dentro de las mismas particiones.
# No se mezclan las antiguas particiones con estas cifras.
PROYECTO <- Sys.getenv("MA_PROJECT_DIR",unset=getwd())
REVISION <- Sys.getenv("MA_REVISION_DIR",unset=file.path(PROYECTO,"revision-observaciones-2026-10-09"))
generar_prediccion_agrupada <- function() {
  src <- file.path(PROYECTO,"revision-3-revista-colombiana-psicologia-2026-10-08/analisis_sensibilidad/grouped_cv")
  out <- file.path(REVISION,"tablas_finales")
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  agregados <- file.path(REVISION,"evidencia/observaciones_2_3_6_7")
# [FUENTES / T3, S9] Primero se buscan copias auditadas de resúmenes.
# Sólo si no existen, se usa la salida original grouped_cv del mismo proyecto.
  ruta_a <- file.path(agregados,"archivado_resumen_modelos_por_item.csv")
  ruta_b <- file.path(agregados,"archivado_resumen_diferencias_pareadas.csv")
  if(!file.exists(ruta_a))ruta_a <- file.path(src,"resumen_modelos_por_item.csv")
  if(!file.exists(ruta_b))ruta_b <- file.path(src,"resumen_diferencias_pareadas.csv")
  a <- read.csv(ruta_a,stringsAsFactors=FALSE)
  b <- read.csv(ruta_b,stringsAsFactors=FALSE)
  a <- a[a$scenario=="all",];b <- b[b$scenario=="all",]
  f <- function(x,d=4)formatC(x,format="f",digits=d,decimal.mark=".")
# [AUXILIAR / T3, S9] metric() extrae una métrica por ítem/estimador/bloque;
# delta() extrae la diferencia pareada. Ambas exigen exactamente una fila.
# Se conserva el promedio original de diferencias, no un nuevo cociente.
  metric <- function(item,est,model,col){
    r <- a[a$item==item & a$estimator==est & a$model==model,,drop=FALSE]
    stopifnot(nrow(r)==1L);r[[col]]
  }
  delta <- function(item,est,col){
    r <- b[b$item==item & b$estimator==est,,drop=FALSE]
    stopifnot(nrow(r)==1L);r[[col]]
  }
  # T3: errores en la escala original de cada respuesta.
  rows <- list();k<-0L
  for(item in paste0("PP",1:3)) for(est in c("ordinal_probit","OLS")){
    k<-k+1L
    rows[[k]] <- c(item,if(est=="OLS")"OLS" else "Probit",
      f(metric(item,est,"B3_APS_without_values","RMSE")),
      f(metric(item,est,"B2_add_APS","RMSE")),
      f(delta(item,est,"mean_difference_RMSE")),
      f(delta(item,est,"mean_relative_RMSE_change_percent"),2))
  }
  T3 <- rbind(c("Item","Model","B3 RMSE","B2 RMSE","Delta RMSE","Change %"),do.call(rbind,rows))
  # S9: probabilidades ordinales y comparador numerico; RPS de OLS no se
  # calcula porque OLS no proporciona distribucion por categoria.
  S9 <- rbind(c("Item","RPS B3","RPS B2","Delta RPS","MAE B3","MAE B2"),
    do.call(rbind,lapply(paste0("PP",1:3),function(item)c(item,
      f(metric(item,"ordinal_probit","B3_APS_without_values","RPS"),5),
      f(metric(item,"ordinal_probit","B2_add_APS","RPS"),5),
      f(delta(item,"ordinal_probit","mean_difference_RPS"),5),
      f(metric(item,"ordinal_probit","B3_APS_without_values","MAE")),
      f(metric(item,"ordinal_probit","B2_add_APS","MAE"))))))
# [SALIDAS] T3.csv y S9.csv tienen encabezado dentro de la primera fila.
# RMSE y MAE:4 decimales; RPS:5; cambio porcentual:2; sin reestimar métricas.
  for(id in c("T3","S9"))write.table(get(id),file.path(out,paste0(id,".csv")),
    sep=",",row.names=FALSE,col.names=FALSE,fileEncoding="UTF-8")
  invisible(list(T3=T3,S9=S9))
}
if(sys.nframe()==0L)generar_prediccion_agrupada()
