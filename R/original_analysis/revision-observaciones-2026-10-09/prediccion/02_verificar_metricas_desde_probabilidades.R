# AUDITORIA INDEPENDIENTE DE CIFRAS -- no entrena ni selecciona modelos.
# Reconstruye RMSE, MAE y RPS desde probabilidades privadas usando otra
# implementacion; coteja las matrices CSV usadas por el manuscrito.
# Revisa separacion de grupos, folds interiores y pesos de clones bootstrap.
# Ejecutar despues de finalizar 01_regularizacion_e_incertidumbre_agrupada.R.
options(stringsAsFactors=FALSE)
project <- Sys.getenv("PAPER_PROYECTO","/Users/luisbenites/Documents/Codex/Paper con MA - 1/reproduccion-integral-V3-2026-10-10/proyecto")
out <- file.path(project,"revision-observaciones-2026-10-09/prediccion")
.libPaths(c(file.path(project,"revision-revista-2026-10-04/auditoria-codigos/R-library"),.libPaths()))
stopifnot(requireNamespace("jsonlite",quietly=TRUE))
read_csv <- function(name) read.csv(file.path(out,name),check.names=FALSE)
sha <- function(path)strsplit(system2("/usr/bin/shasum",c("-a","256",shQuote(path)),stdout=TRUE)," +")[[1]][1]
report <- jsonlite::read_json(file.path(out,"reporte_regularizacion_incertidumbre.json"),simplifyVector=TRUE)
qpath <- file.path(project,"revision-3-revista-colombiana-psicologia-2026-10-08/analisis_sensibilidad/privados/datos_y_banderas_para_sensibilidades.rds")
q <- readRDS(qpath)
rawpath <- file.path(project,"revision-integral-2026-10-04/analysis/outputs/canonical-full-20261004-02/01_audit/analytic_data.rds")
raw <- readRDS(rawpath)
dat <- q$data;dat$Edad_original<-raw$Edad;dat$NSE_original<-raw$NSE
stopifnot(sha(qpath)==report$source_sha256_before$q,
          sha(rawpath)==report$source_sha256_before$raw,
          identical(report$source_sha256_before,report$source_sha256_after))
errors <- character(); diffs <- numeric(); reconstruction <- list();cal_reconstruction<-list()
add_error <- function(x)errors <<- c(errors,x)
check <- function(x,description)if(!isTRUE(x))add_error(description)
new_metrics <- function(y,P) {
  check(ncol(P)==5L,"Cinco categorias no retenidas")
  check(all(is.finite(P))&&all(P>0),"Probabilidad no positiva o no finita")
  check(max(abs(rowSums(P)-1))<=1e-10,"Probabilidades no suman uno")
  expected <- rowSums(sweep(P,2,1:5,"*"))
  rps_total <- 0
  for(k in 1:4)rps_total <- rps_total+sum((rowSums(P[,seq_len(k),drop=FALSE])-as.numeric(y<=k))^2)
  c(RMSE=sqrt(sum((y-expected)^2)/length(y)),MAE=sum(abs(y-expected))/length(y),RPS=rps_total/(4*length(y)))
}
compare_row <- function(z,y,P,description) {
  m <- new_metrics(y,P)
  for(v in names(m)) {
    delta <- abs(m[[v]]-z[[v]])
    diffs <<- c(diffs,delta)
    if(!is.finite(delta)||delta>1e-10)add_error(paste(description,v,"mismatch"))
  }
  for(k in 1:4) {
    cumulative<-rowSums(P[,seq_len(k),drop=FALSE])
    cal_reconstruction[[length(cal_reconstruction)+1L]] <<- data.frame(
      Design=if(grepl("^CV",z$Design))"CV" else "OOB",Replication=z$Replication,Outcome=z$Outcome,Block=z$Block,Cut=k,
      Predicted=mean(cumulative),Observed=mean(y<=k),Gap=mean(cumulative)-mean(y<=k))
  }
  m
}
cv <- read_csv("CV_regularizado_metricas_por_repeticion.csv")
cvRDS <- readRDS(file.path(out,"privados","CV_predicciones_agrupadas.rds"))
check(length(cvRDS)==5L,"No hay cinco repeticiones CV")
for(r in seq_along(cvRDS)) {
  zz <- cvRDS[[r]]
  check(identical(as.integer(zz$original_row),as.integer(q$original_row)),"Orden CV no coincide")
  for(y in c("PP1","PP2","PP3"))for(b in c("B3","B2")) {
    z <- cv[cv$Replication==r&cv$Outcome==y&cv$Block==b,,drop=FALSE]
    check(nrow(z)==1L&&z$Valid,"CV incompleta")
    if(nrow(z)==1L&&z$Valid) {
      m<-compare_row(z,zz$y[[y]],zz$p[[y]][[b]],paste("CV",r,y,b))
      reconstruction[[length(reconstruction)+1L]]<-data.frame(Design="CV",Replication=r,Outcome=y,Block=b,t(m))
    }
  }
}
boot <- read_csv("OOB_regularizado_metricas_por_replica.csv")
files <- list.files(file.path(out,"privados","OOB_predicciones"),pattern="^replica_[0-9]+[.]rds$",full.names=TRUE)
check(length(files)==report$bootstrap_repetitions,"Numero de replicas OOB inconsistente")
partition_checks <- list(); cal_bands <- list()
for(file in files) {
  r <- as.integer(sub("replica_([0-9]+)[.]rds$","\\1",basename(file)))
  zz <- readRDS(file)
  check(length(intersect(q$group[zz$train],q$group[zz$test]))==0L,paste("Grupo entre train/test",r))
  check(identical(sort(c(zz$train,zz$test)),seq_len(653L)),paste("Particion incompleta",r))
  check(identical(as.integer(zz$original_row),as.integer(q$original_row[zz$test])),paste("Orden OOB",r))
  check(all(vapply(split(zz$weights,q$group[zz$train]),function(v)length(unique(v))==1L,logical(1))),paste("Pesos por grupo",r))
  check(all(vapply(split(zz$inner_folds,q$group[zz$train]),function(v)length(unique(v))==1L,logical(1))),paste("Clones cruzan folds interiores",r))
  for(y in c("PP1","PP2","PP3"))for(b in c("B3","B2")) {
    pred <- zz$predictions[[paste(y,b,sep="_")]]
    z <- boot[boot$Replication==r&boot$Outcome==y&boot$Block==b,,drop=FALSE]
    check(nrow(z)==1L,paste("Registro OOB ausente",r,y,b))
    if(nrow(z)==1L&&pred$valid) {
      check(z$Valid,paste("Estado OOB discordante",r,y,b))
      m<-compare_row(z,zz$y[[y]],pred$p,paste("OOB",r,y,b))
      reconstruction[[length(reconstruction)+1L]]<-data.frame(Design="OOB",Replication=r,Outcome=y,Block=b,t(m))
      # Auditar que la escala final usa solo los pesos y filas de entrenamiento.
      vars <- names(pred$scaler$mean)
      X <- as.matrix(dat[zz$train,vars,drop=FALSE])
      expected_mean <- colSums(X*zz$weights)/sum(zz$weights)
      expected_sd <- sqrt(colSums(sweep(X,2,expected_mean,"-")^2*zz$weights)/(sum(zz$weights)-1))
      expected_sd[!is.finite(expected_sd)|expected_sd<1e-12]<-1
      check(max(abs(pred$scaler$mean-expected_mean))<1e-10,paste("Media fuera train",r,y,b))
      check(max(abs(pred$scaler$sd-expected_sd))<1e-10,paste("SD fuera train",r,y,b))
      for(k in 1:4) {
        cumulative <- rowSums(pred$p[,seq_len(k),drop=FALSE]);bin<-pmin(5L,floor(cumulative*5)+1L)
        for(band in 1:5) {
          id <- which(bin==band)
          cal_bands[[length(cal_bands)+1L]]<-data.frame(Outcome=y,Block=b,Cut=k,Band_index=band,N=length(id),
             Predicted_sum=sum(cumulative[id]),Observed_sum=sum(zz$y[[y]][id]<=k))
        }
      }
    }
  }
  partition_checks[[length(partition_checks)+1L]]<-data.frame(Replication=r,Distinct_train_rows=length(zz$train),Train_weight=sum(zz$weights),N_OOB=length(zz$test),Groups_OOB=length(unique(q$group[zz$test])))
}
rec<-do.call(rbind,reconstruction)
csvsummary <- read_csv("OOB_regularizado_resumen.csv")
for(i in seq_len(nrow(csvsummary))) {
  z<-csvsummary[i,];metric<-z$Metric
  z3<-rec[rec$Design=="OOB"&rec$Outcome==z$Outcome&rec$Block=="B3",];z2<-rec[rec$Design=="OOB"&rec$Outcome==z$Outcome&rec$Block=="B2",]
  pair<-merge(z3[,c("Replication",metric)],z2[,c("Replication",metric)],by="Replication",suffixes=c("_B3","_B2"))
  d<-pair[[paste0(metric,"_B2")]]-pair[[paste0(metric,"_B3")]];relative<-100*d/pair[[paste0(metric,"_B3")]]
  qq<-as.numeric(quantile(d,c(.025,.975),type=7));rq<-as.numeric(quantile(relative,c(.025,.975),type=7))
  values<-c(B3_mean=mean(pair[[paste0(metric,"_B3")]]),B2_mean=mean(pair[[paste0(metric,"_B2")]]),Delta_mean=mean(d),
       Delta_lower=qq[1],Delta_upper=qq[2],Relative_percent_mean=mean(relative),Relative_percent_lower=rq[1],Relative_percent_upper=rq[2])
  for(v in names(values)) {
    dd<-abs(z[[v]]-values[[v]]);diffs<-c(diffs,dd);check(is.finite(dd)&&dd<1e-10,paste("Resumen OOB",z$Outcome,metric,v))
  }
}
cvsum<-read_csv("CV_regularizado_comparacion.csv")
for(i in seq_len(nrow(cvsum)))for(metric in c("RMSE","MAE","RPS")) {
  z<-cvsum[i,];z3<-rec[rec$Design=="CV"&rec$Outcome==z$Outcome&rec$Block=="B3",];z2<-rec[rec$Design=="CV"&rec$Outcome==z$Outcome&rec$Block=="B2",]
  pair<-merge(z3[,c("Replication",metric)],z2[,c("Replication",metric)],by="Replication",suffixes=c("_B3","_B2"))
  values<-c(mean(pair[[paste0(metric,"_B3")]]),mean(pair[[paste0(metric,"_B2")]]),
      mean(pair[[paste0(metric,"_B2")]]-pair[[paste0(metric,"_B3")]]),mean(100*(pair[[paste0(metric,"_B2")]]-pair[[paste0(metric,"_B3")]])/pair[[paste0(metric,"_B3")]]))
  names(values)<-c(paste0("B3_",metric),paste0("B2_",metric),paste0("Delta_",metric),paste0("Delta_",metric,"_percent"))
  for(v in names(values)) {
    dd<-abs(z[[v]]-values[[v]]);diffs<-c(diffs,dd);check(is.finite(dd)&&dd<1e-10,paste("Resumen CV",z$Outcome,metric,v))
  }
}
if(length(cal_bands)) {
  cb<-aggregate(cbind(N,Predicted_sum,Observed_sum)~Outcome+Block+Cut+Band_index,data=do.call(rbind,cal_bands),FUN=sum)
  cb$Predicted<-ifelse(cb$N>0,cb$Predicted_sum/cb$N,NA_real_);cb$Observed<-ifelse(cb$N>0,cb$Observed_sum/cb$N,NA_real_);cb$Gap<-cb$Predicted-cb$Observed
  write.csv(cb,file.path(out,"OOB_regularizado_calibracion_bandas.csv"),row.names=FALSE,fileEncoding="UTF-8",na="")
}
cal_all<-do.call(rbind,cal_reconstruction)
cv_cal_csv<-read_csv("CV_regularizado_calibracion_cortes.csv")
oob_cal_csv<-read_csv("OOB_regularizado_calibracion_cortes.csv")
for(i in seq_len(nrow(cv_cal_csv))) {
  z<-cv_cal_csv[i,];xx<-cal_all[cal_all$Design=="CV"&cal_all$Outcome==z$Outcome&cal_all$Block==z$Block&cal_all$Cut==z$Cut,]
  for(v in c("Predicted","Observed","Gap")) {
    dd<-abs(z[[v]]-mean(xx[[v]]));diffs<-c(diffs,dd);check(dd<1e-10,paste("CV calibracion",z$Outcome,z$Block,z$Cut,v))
  }
}
for(i in seq_len(nrow(oob_cal_csv))) {
  z<-oob_cal_csv[i,];xx<-cal_all[cal_all$Design=="OOB"&cal_all$Outcome==z$Outcome&cal_all$Block==z$Block&cal_all$Cut==z$Cut,]
  qq<-as.numeric(quantile(xx$Gap,c(.025,.975),type=7))
  v<-c(Predicted_mean=mean(xx$Predicted),Observed_mean=mean(xx$Observed),Gap_mean=mean(xx$Gap),Gap_lower=qq[1],Gap_upper=qq[2])
  for(nm in names(v)) {
    dd<-abs(z[[nm]]-v[[nm]]);diffs<-c(diffs,dd);check(dd<1e-10,paste("OOB calibracion",z$Outcome,z$Block,z$Cut,nm))
  }
}
# Frecuencia de lambda finalmente elegida: agregado, sin coeficientes/filas.
d<-read_csv("diagnosticos_ajustes.csv");d<-d[grepl("/final$",d$Tag),]
tags<-strsplit(d$Tag,"/",fixed=TRUE)
selection<-data.frame(Design=vapply(tags,`[`,character(1),1),Outcome=vapply(tags,function(x)x[length(x)-2L],character(1)),
   Block=vapply(tags,function(x)x[length(x)-1L],character(1)),Lambda=d$Lambda,Count=1L)
selection<-aggregate(Count~Design+Outcome+Block+Lambda,data=selection,FUN=sum)
write.csv(selection,file.path(out,"frecuencia_lambda_seleccionada.csv"),row.names=FALSE,fileEncoding="UTF-8")
frequencies <- do.call(rbind,lapply(c("PP1","PP2","PP3"),function(y)data.frame(Outcome=y,Category=1:5,N=as.integer(table(factor(q$data[[y]],levels=1:5))))))
write.csv(frequencies,file.path(out,"frecuencias_categorias_PP.csv"),row.names=FALSE,fileEncoding="UTF-8")
write.csv(do.call(rbind,partition_checks),file.path(out,"auditoria_particiones_OOB.csv"),row.names=FALSE,fileEncoding="UTF-8")
jsonlite::write_json(list(status=if(length(errors))"FAIL" else "PASS",errors=errors,checks_metric_values=length(diffs),
  max_absolute_difference=if(length(diffs))max(diffs) else NA_real_,
  reconstructed_CV_rows=sum(rec$Design=="CV"),reconstructed_OOB_rows=sum(rec$Design=="OOB"),
  checked_OOB_replicates=length(files),all_categories_positive=!any(grepl("Probabilidad no positiva|Cinco categorias no retenidas",errors)),
  note="Independent metric formulas and reconstruction, conditional grouped partitions, final train-only weighted scalers. Coverage not validated; no new fits.",
  sources_unchanged=sha(qpath)==report$source_sha256_after$q&&sha(rawpath)==report$source_sha256_after$raw,
  sessionInfo=capture.output(sessionInfo())),file.path(out,"auditoria_independiente_probabilidades.json"),auto_unbox=TRUE,pretty=TRUE,digits=16)
if(length(errors))stop(paste(errors,collapse="\n"))
report$status<-"Completed; independent probability reconstruction PASS"
report$independent_audit_status<-"PASS"
report$independent_audit_file<-"auditoria_independiente_probabilidades.json"
jsonlite::write_json(report,file.path(out,"reporte_regularizacion_incertidumbre.json"),auto_unbox=TRUE,pretty=TRUE,digits=16,null="null",na="null")
cat("PASS: reconstruccion independiente",nrow(rec),"filas; max diferencia",max(diffs),"\n")
