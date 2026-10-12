# [CALIBRACION_CV_PROBIT] / S19: calibracion de probabilidades probit YA
# archivadas, con exactamente las mismas 5x5 particiones agrupadas del ridge.
# NO ajusta modelos, NO reconstruye probabilidades aproximadas, NO usa Word.
# Cinco columnas de probabilidades originales, incluidos ceros por categorias
# ausentes. Calibracion marginal descriptiva; no prueba calibracion condicional.
options(stringsAsFactors=FALSE)
project <- Sys.getenv("PAPER_PROYECTO","/Users/luisbenites/Documents/Codex/Paper con MA - 1/reproduccion-integral-V3-2026-10-10/proyecto")
.libPaths(c(file.path(project,"revision-revista-2026-10-04/auditoria-codigos/R-library"),.libPaths()))
stopifnot(requireNamespace("jsonlite",quietly=TRUE))
out <- file.path(project,"revision-observaciones-2026-10-09/prediccion")
source_dir <- file.path(project,"revision-3-revista-colombiana-psicologia-2026-10-08/analisis_sensibilidad/grouped_cv/all/privados")
predpath <- file.path(source_dir,"predicciones_fuera_muestra.rds")
qpath <- file.path(project,"revision-3-revista-colombiana-psicologia-2026-10-08/analisis_sensibilidad/privados/datos_y_banderas_para_sensibilidades.rds")
sha <- function(path)strsplit(system2("/usr/bin/shasum",c("-a","256",shQuote(path)),stdout=TRUE)," +")[[1]][1]
before <- c(predictions=sha(predpath),quality=sha(qpath),assignments=sha(file.path(source_dir,"asignaciones_por_grupo.rds")))
q <- readRDS(qpath);assign <- readRDS(file.path(source_dir,"asignaciones_por_grupo.rds"));pred<-readRDS(predpath)
stopifnot(before["quality"]=="283d622d631faabfce1da2e5c947188d6b006d40df0f6b2690945d24f99be76e",
          identical(as.integer(assign$group),as.integer(q$group)))
selected<-pred[pred$estimator=="ordinal_probit"&pred$model%in%c("B2_add_APS","B3_APS_without_values"),]
if(!nrow(selected)) stop("No hay filas ordinal_probit; revisar etiqueta exacta sin aproximar probabilidades")
checkcols<-paste0("probability_",1:5)
stopifnot(all(checkcols%in%names(selected)))
rows<-list();metrics<-list();max_expected_diff<-0;max_probability_error<-0;zero_observed<-0L
for(r in 1:5)for(y in c("PP1","PP2","PP3"))for(b in c("B3","B2")) {
  model<-if(b=="B2")"B2_add_APS" else "B3_APS_without_values"
  z<-selected[selected$repetition==r&selected$item==y&selected$model==model,]
  stopifnot(nrow(z)==653L,identical(sort(z$original_row),seq_len(653L)),!anyDuplicated(z$original_row),
    all(z$fold==assign$assignment[[r]][z$original_row]),
    all(z$observed==q$data[[y]][z$original_row]),
    all(vapply(split(z$fold,q$group[z$original_row]),function(v)length(unique(v))==1L,logical(1))))
  P<-as.matrix(z[,checkcols]);stopifnot(all(is.finite(P)),all(P>=0))
  pe<-max(abs(rowSums(P)-1));stopifnot(pe<1e-10);max_probability_error<-max(max_probability_error,pe)
  expectation<-as.vector(P%*%(1:5));ed<-max(abs(expectation-z$predicted_expected));stopifnot(ed<1e-10);max_expected_diff<-max(max_expected_diff,ed)
  zero_observed<-zero_observed+sum(P[cbind(seq_len(653L),as.integer(z$observed))]==0)
  for(k in 1:4) {
    Fp<-rowSums(P[,seq_len(k),drop=FALSE]);observed<-mean(z$observed<=k)
    rows[[length(rows)+1L]]<-data.frame(Replication=r,Outcome=y,Block=b,Cut=k,N=653,
      Predicted=mean(Fp),Observed=observed,Gap=mean(Fp)-observed)
  }
  rps<-mean(vapply(seq_len(653L),function(i)sum((cumsum(P[i,])[1:4]-as.numeric(z$observed[i]<=1:4))^2)/4,numeric(1)))
  metrics[[length(metrics)+1L]]<-data.frame(Replication=r,Outcome=y,Block=b,RMSE=sqrt(mean((expectation-z$observed)^2)),MAE=mean(abs(expectation-z$observed)),RPS=rps)
}
detail<-do.call(rbind,rows)
aggregate_rows<-aggregate(cbind(Predicted,Observed,Gap)~Outcome+Block+Cut,data=detail,FUN=mean)
write.csv(aggregate_rows,file.path(out,"CV_probit_calibracion_cortes.csv"),row.names=FALSE,fileEncoding="UTF-8")
write.csv(detail,file.path(out,"CV_probit_calibracion_cortes_por_repeticion.csv"),row.names=FALSE,fileEncoding="UTF-8")
write.csv(do.call(rbind,metrics),file.path(out,"CV_probit_metricas_reconstruidas_por_repeticion.csv"),row.names=FALSE,fileEncoding="UTF-8")
# Segunda formula: sumar cada columna probit original directamente en los
# 3265 tests por item/bloque. Como cada repeticion contiene 653, coincide con
# la media de las cinco calibraciones, sin darles una interpretacion de SE.
errors<-numeric(nrow(aggregate_rows))
for(i in seq_len(nrow(aggregate_rows))) {
  zz<-aggregate_rows[i,];model<-if(zz$Block=="B2")"B2_add_APS" else "B3_APS_without_values"
  z<-selected[selected$item==zz$Outcome&selected$model==model,]
  P<-as.matrix(z[,checkcols]);direct<-sum(P[,seq_len(zz$Cut),drop=FALSE])/nrow(z)
  errors[i]<-abs(direct-zz$Predicted)
  stopifnot(abs(mean(z$observed<=zz$Cut)-zz$Observed)<1e-10,abs(direct-mean(z$observed<=zz$Cut)-zz$Gap)<1e-10)
}
after<-c(predictions=sha(predpath),quality=sha(qpath),assignments=sha(file.path(source_dir,"asignaciones_por_grupo.rds")))
stopifnot(identical(before,after),max(errors)<1e-10)
jsonlite::write_json(list(status="PASS",new_fits=0,source_sha256_before=as.list(before),source_sha256_after=as.list(after),
  expected_rows=19590,actual_rows=nrow(selected),repetitions=5,N_each_repetition=653,aggregate_rows=nrow(aggregate_rows),
  probability_sum_max_error=max_probability_error,expected_score_max_error=max_expected_diff,
  independent_calibration_formula_max_error=max(errors),zero_observed_probabilities=zero_observed,
  data_order_and_group_fold_checks=TRUE,
  description="Exact archived probabilities only; same grouped CV5x5; original observed proportion unchanged across repetitions; no approximations or new fits"),
  file.path(out,"auditoria_calibracion_probit_CV.json"),auto_unbox=TRUE,pretty=TRUE,digits=16)
cat("PASS: calibracion exacta probit",nrow(aggregate_rows),"filas agregadas;",nrow(selected),"tests archivados.\n")
