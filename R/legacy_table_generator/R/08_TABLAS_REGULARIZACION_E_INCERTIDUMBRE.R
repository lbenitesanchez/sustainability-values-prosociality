# MÓDULO 08 — TABLAS S17, S18, S19, S20.
# Función pública: formatear_extension_predictiva().
# Sólo formatea CSV guardados de ridge, OOB, calibración y LISO.
# No vuelve a entrenar, ajustar lambda ni calcular incertidumbre con datos nuevos.
# Una tabla sintética sería una demostración aparte y no reproduce estas cifras.

# Presentacion de los resultados nuevos sin copiar numeros desde Word.
# Se ejecuta DESPUES de prediccion/01_regularizacion_e_incertidumbre_agrupada.R.
REVISION <- Sys.getenv("MA_REVISION_DIR", unset = file.path(getwd(),"revision-observaciones-2026-10-09"))
formatear_extension_predictiva <- function() {
  indir <- file.path(REVISION,"prediccion");out <- file.path(REVISION,"tablas_finales")
  dir.create(out,showWarnings=FALSE)
# [AUXILIARES / S17–S20] f conserva NA y redondea sólo texto final.
# rango() ensambla límites ya calculados; no crea intervalos nuevos.
  f <- function(x,n=4L)ifelse(is.finite(x),formatC(x,format="f",digits=n,decimal.mark="."),"NA")
  rango <- function(a,b,n=4L)paste0("[",f(a,n),", ",f(b,n),"]")
  guardar <- function(x,id)write.table(x,file.path(out,paste0(id,".csv")),sep=",",
    row.names=FALSE,col.names=FALSE,fileEncoding="UTF-8")
# [S17] CV_regularizado_comparacion.csv: B3 y B2, diferencias absolutas
# y relativas pareadas por PP y métrica. RMSE/MAE/RPS a 5 decimales; cambio a2.
  # [REGULARIZED_CV] Media de cinco repeticiones, cada una con 653 evaluados.
  a <- read.csv(file.path(indir,"CV_regularizado_comparacion.csv"),stringsAsFactors=FALSE)
  stopifnot(nrow(a)==3L,all(a$Complete_repetitions==5L))
  filas <- list();i <- 0L
  for(j in seq_len(nrow(a)))for(metrica in c("RMSE","MAE","RPS")) {
    i <- i+1L
    filas[[i]] <- c(a$Outcome[j],metrica,f(a[[paste0("B3_",metrica)]][j],5),
      f(a[[paste0("B2_",metrica)]][j],5),f(a[[paste0("Delta_",metrica)]][j],5),
      f(a[[paste0("Delta_",metrica,"_percent")]][j],2))
  }
  x <- rbind(c("Item","Metric","B3","B2","Difference","Change %"),do.call(rbind,filas))
  guardar(x,"S17")
# [S18] OOB_regularizado_resumen.csv:200 réplicas solicitadas,9 filas.
# El control de completitud de la auditoría acredita las réplicas existentes;
# este formateador conserva puntos y percentiles del mismo diseño OOB.
  # [OOB_UNCERTAINTY] El punto y el rango pertenecen SOLO al diseno OOB.
  # Los percentiles son variacion de remuestreo con reentrenamiento, sin
  # cobertura nominal validada y sin inferencia poblacional de instituciones.
  a <- read.csv(file.path(indir,"OOB_regularizado_resumen.csv"),stringsAsFactors=FALSE)
  stopifnot(nrow(a)==9L,all(a$Replicates_requested==200L))
  x <- rbind(c("Item","Metric","Mean delta","2.5–97.5 range","Mean change %","2.5–97.5 range %"),
     cbind(a$Outcome,a$Metric,f(a$Delta_mean),rango(a$Delta_lower,a$Delta_upper),
       f(a$Relative_percent_mean,2),rango(a$Relative_percent_lower,a$Relative_percent_upper,2)))
  guardar(x,"S18")
# [S19] Se unen calibraciones Probit y Ridge. Para cada PP y corte 1–4,
# compara porcentaje observado/predicho B3/B2 y brecha en puntos porcentuales.
# Se exige misma frecuencia observada entre bloques antes de mostrarla.
  # [CALIBRATION] Cuatro puntos acumulados, con frecuencia observada como
  # referencia. Los CSV de bandas incluyen sus denominadores; no se tratan
  # las cinco repeticiones como 3265 personas independientes.
  a <- read.csv(file.path(indir,"CV_regularizado_calibracion_cortes.csv"),stringsAsFactors=FALSE)
  a$Model <- "Ridge"
  pr <- read.csv(file.path(indir,"CV_probit_calibracion_cortes.csv"),stringsAsFactors=FALSE)
  pr$Model <- "Probit"
  a <- rbind(pr[,names(a)],a)
  rows <- list();i<-0L
  for(model in c("Probit","Ridge"))for(y in paste0("PP",1:3))for(k in 1:4){
    z <- a[a$Model==model&a$Outcome==y&a$Cut==k,,drop=FALSE]
    b3 <- z[z$Block=="B3",];b2<-z[z$Block=="B2",]
    stopifnot(nrow(b3)==1L,nrow(b2)==1L,abs(b3$Observed-b2$Observed)<1e-10)
    i<-i+1L
    rows[[i]] <- c(model,y,paste0("Y <= ",k),f(100*b3$Observed,2),f(100*b3$Predicted,2),
                  f(100*b2$Predicted,2),f(100*b3$Gap,2),f(100*b2$Gap,2))
  }
  guardar(rbind(c("Model","Item","Cut","Observed %","B3 predicted %","B2 predicted %","B3 gap pp","B2 gap pp"),do.call(rbind,rows)),"S19")
  # [S20] Se reutiliza la sensibilidad por institución ejecutada con estos datos.
  # No se presentan las cuatro instituciones como muestra aleatoria de países.
# [S20 / FUENTES] Contrastes por institución y agregado LISO guardados.
# Selección explícita ordinal_probit:12 filas institucionales + 3 agrupadas.
# El N agrupado 653 no convierte instituciones en muestra nacional.
  lisodir <- file.path(indir,"LISO_previo_agregados")
  a <- read.csv(file.path(lisodir,"contrastes_pareados_por_institucion.csv"),stringsAsFactors=FALSE)
  a <- a[a$estimator=="ordinal_probit",]
  b <- read.csv(file.path(lisodir,"S14_LISO_resumen_agrupado_numerico.csv"),stringsAsFactors=FALSE)
  b <- b[b$Estimator=="ordinal_probit",]
  stopifnot(nrow(a)==12L,nrow(b)==3L)
  z <- rbind(cbind(paste(a$held_out_site,a$item,sep=" / "),a$n_test,
      f(a$B3_RMSE),f(a$B2_RMSE),f(a$delta_RMSE_percent,2),f(a$delta_RPS,5)),
    cbind(paste("Pooled",b$Outcome,sep=" / "),653,
      f(b$B3_RMSE),f(b$B2_RMSE),f(b$Delta_RMSE_percent,2),f(b$B2_RPS-b$B3_RPS,5)))
  guardar(rbind(c("Institution / item","N","B3 RMSE","B2 RMSE","Change %","Delta RPS"),z),"S20")
  invisible(TRUE)
}
if(sys.nframe()==0L)formatear_extension_predictiva()
