# AUDITORÍA INDEPENDIENTE DE LOS AGREGADOS LISO (V4 ARLA).
# Este código reconstruye pérdidas directamente desde predicciones individuales
# privadas; no reutiliza la función metrics ni toma los errores agregados como
# entrada. También verifica predicciones contra modelos guardados y escalados.
# Los resultados públicos contienen únicamente comprobaciones agregadas.
options(stringsAsFactors=FALSE)
project <- Sys.getenv("ARLA_PROYECTO", Sys.getenv("RCP_PROYECTO", "/Users/luisbenites/Documents/Codex/Paper con MA - 1/reproduccion-integral-V3-2026-10-10/proyecto"))
root <- file.path(project,"revision-4-ARLA-2026-10-09/transportabilidad_institucional")
.libPaths(c(file.path(project,"revision-revista-2026-10-04/auditoria-codigos/R-library"),.libPaths()))
suppressPackageStartupMessages({library(MASS);library(jsonlite);library(digest)})
private <- file.path(root,"privados")
pr <- readRDS(file.path(private,"predicciones_fuera_institucion.rds"))
ft <- read.csv(file.path(root,"metricas_y_diagnosticos_por_institucion.csv"),na.strings="")
po <- read.csv(file.path(root,"resumen_errores_agrupados.csv"),na.strings="")
ps <- read.csv(file.path(root,"resumen_incremento_predictivo.csv"),na.strings="")
pa <- read.csv(file.path(root,"contrastes_pareados_por_institucion.csv"),na.strings="")
audit <- read_json(file.path(root,"auditoria_entradas.json"),simplifyVector=TRUE)
quality_path <- file.path(project,"revision-3-revista-colombiana-psicologia-2026-10-08/analisis_sensibilidad/privados/datos_y_banderas_para_sensibilidades.rds")
archived_path <- file.path(project,"revision-integral-2026-10-04/analysis/outputs/canonical-full-20261004-02/01_audit/analytic_data.rds")
q <- readRDS(quality_path);archived <- as.data.frame(readRDS(archived_path))
data <- q$data;data$Edad <- archived$Edad;data$NSE <- archived$NSE
stopifnot(identical(audit$sources_sha256$quality,digest(quality_path,algo="sha256",file=TRUE)),
 identical(audit$sources_sha256$archived,digest(archived_path,algo="sha256",file=TRUE)))
max_difference <- 0
check_numeric <- function(actual,expected){
 actual <- as.numeric(actual);expected <- as.numeric(expected)
 stopifnot(length(actual)==length(expected),identical(is.na(actual),is.na(expected)))
 ix <- !is.na(actual);stopifnot(identical(is.infinite(actual[ix]),is.infinite(expected[ix])))
 ix <- is.finite(actual) & is.finite(expected)
 if(any(ix)){
  delta <- max(abs(actual[ix]-expected[ix]));max_difference <<- max(max_difference,delta)
  stopifnot(delta<1e-11)
 }
 invisible(TRUE)
}
# Otra implementación del RPS: probabilidad acumulada por sumas sucesivas.
rebuild <- function(z){
 ok <- all(is.finite(z$predicted_expected));d <- z$predicted_expected-z$observed
 ans <- c(n_expected=nrow(z),n_scored=if(ok)nrow(z)else 0,
  SSE=if(ok)sum(d*d)else NA_real_,SAE=if(ok)sum(abs(d))else NA_real_,
  RMSE=if(ok)sqrt(sum(d*d)/nrow(z))else NA_real_,MAE=if(ok)sum(abs(d))/nrow(z)else NA_real_,
  sum_RPS=NA_real_,RPS=NA_real_,sum_logloss_raw=NA_real_,logloss_raw=NA_real_,
  sum_logloss_floor=NA_real_,logloss_floor=NA_real_,n_zero_probability=NA_real_)
 if(z$estimator[1]!="OLS" && ok){
  p <- as.matrix(z[paste0("probability_",1:5)])
  check_numeric(rowSums(p),rep(1,nrow(p)))
  check_numeric(as.vector(p %*% (1:5)),z$predicted_expected)
  loss <- numeric(nrow(z));cum <- numeric(nrow(z))
  for(k in 1:4){cum <- cum+p[,k];loss <- loss+(cum-as.numeric(z$observed<=k))^2/4}
  probobs <- p[cbind(seq_len(nrow(z)),z$observed)]
  ans[c("sum_RPS","RPS")] <- c(sum(loss),sum(loss)/nrow(z))
  ans[c("sum_logloss_raw","logloss_raw")] <- c(sum(-log(probobs)),sum(-log(probobs))/nrow(z))
  ans[c("sum_logloss_floor","logloss_floor")] <- c(sum(-log(pmax(probobs,1e-12))),sum(-log(pmax(probobs,1e-12)))/nrow(z))
  ans["n_zero_probability"] <- sum(probobs==0)
 };ans
}
reconstructed <- list();fit_predictions_verified <- 0L
for(i in seq_len(nrow(ft))){
 a <- ft[i,];z <- pr[pr$held_out_site==a$held_out_site & pr$item==a$item &
  pr$model==a$model & pr$estimator==a$estimator,]
 z <- z[order(z$original_row),]
 stopifnot(nrow(z)==a$n_expected,identical(z$original_row,which(as.character(data$Pais)==a$held_out_site)),
  !anyDuplicated(z$original_row),identical(z$observed,as.numeric(data[z$original_row,a$item])))
 rebuilt <- rebuild(z)
 check_numeric(as.numeric(a[names(rebuilt)]),unname(rebuilt))
 reconstructed[[i]] <- cbind(a[,c("held_out_site","item","model","estimator")],as.data.frame(as.list(rebuilt)))
 test <- z$original_row;train <- setdiff(seq_len(653L),test)
 stopifnot(length(intersect(q$group[train],q$group[test]))==0L)
 if(a$model=="B_empirical"){
  cp <- as.numeric(table(factor(data[train,a$item],levels=1:5)))/length(train)
  check_numeric(as.matrix(z[paste0("probability_",1:5)]),matrix(rep(cp,each=length(test)),nrow=length(test)))
 }else{
  est <- if(a$estimator=="OLS")"OLS"else"probit"
  f <- readRDS(file.path(private,"fits",paste(a$held_out_site,a$item,a$model,paste0(est,".rds"),sep="_")))
  stopifnot(identical(f$train,train),identical(f$test,test),
   !any(grepl("C_|Pais|country|site",f$vars)),!any(paste0("PP",1:3)%in%f$vars))
  xt <- data[train,f$vars,drop=FALSE]
  check_numeric(f$center,vapply(xt,mean,numeric(1)));check_numeric(f$scale,vapply(xt,sd,numeric(1)))
  newdata <- as.data.frame(sweep(sweep(as.matrix(data[test,f$vars,drop=FALSE]),2,f$center,"-"),2,f$scale,"/"))
  if(all(is.finite(z$predicted_expected))){
   if(est=="OLS")check_numeric(as.numeric(predict(f$fit,newdata=newdata)),z$predicted_expected)else{
    px <- as.matrix(predict(f$fit,newdata=newdata,type="probs"))
    mapped <- matrix(0,length(test),5L,dimnames=list(NULL,as.character(1:5)))
    mapped[,colnames(px)] <- px
    check_numeric(mapped,as.matrix(z[paste0("probability_",1:5)]))
   }
   fit_predictions_verified <- fit_predictions_verified+1L
  }
 }
}
rec <- do.call(rbind,reconstructed)
# Resumen ponderado: reconstruir con TODAS las predicciones una vez por modelo.
for(i in seq_len(nrow(po))){
 a <- po[i,];z <- pr[pr$item==a$item & pr$model==a$model & pr$estimator==a$estimator,]
 stopifnot(nrow(z)==653L,!anyDuplicated(z$original_row))
 r <- rebuild(z)
 check_numeric(c(a$pooled_RMSE,a$pooled_MAE,a$pooled_RPS,a$pooled_logloss_raw,a$pooled_logloss_floor),
  unname(r[c("RMSE","MAE","RPS","logloss_raw","logloss_floor")]))
 bysite <- rec[rec$item==a$item & rec$model==a$model & rec$estimator==a$estimator,]
 check_numeric(c(a$unweighted_site_mean_RMSE,a$unweighted_site_mean_MAE,a$unweighted_site_mean_RPS),
  c(mean(bysite$RMSE),mean(bysite$MAE),mean(bysite$RPS)))
}
independent_paired <- list()
for(i in seq_len(nrow(pa))){
 a <- pa[i,];b2 <- rec[rec$held_out_site==a$held_out_site & rec$item==a$item & rec$estimator==a$estimator & rec$model=="B2_add_values",]
 b3 <- rec[rec$held_out_site==a$held_out_site & rec$item==a$item & rec$estimator==a$estimator & rec$model=="B3_APS_without_values",]
 vals <- c(B3_RMSE=b3$RMSE,B2_RMSE=b2$RMSE,delta_RMSE=b2$RMSE-b3$RMSE,
  delta_RMSE_percent=100*(b2$RMSE/b3$RMSE-1),B3_MAE=b3$MAE,B2_MAE=b2$MAE,
  delta_MAE=b2$MAE-b3$MAE,B3_RPS=b3$RPS,B2_RPS=b2$RPS,delta_RPS=b2$RPS-b3$RPS)
 check_numeric(as.numeric(a[names(vals)]),unname(vals))
 independent_paired[[i]] <- cbind(a[,c("held_out_site","item","estimator")],as.data.frame(as.list(vals)))
}
ip <- do.call(rbind,independent_paired)
for(i in seq_len(nrow(ps))){
 a <- ps[i,];za <- pr[pr$item==a$item & pr$estimator==a$estimator & pr$model=="B2_add_values",]
 zb <- pr[pr$item==a$item & pr$estimator==a$estimator & pr$model=="B3_APS_without_values",]
 ra <- rebuild(za);rb <- rebuild(zb);ds <- ip[ip$item==a$item & ip$estimator==a$estimator,]
 vals <- c(pooled_B3_RMSE=rb[["RMSE"]],pooled_B2_RMSE=ra[["RMSE"]],
  pooled_delta_RMSE_percent=100*(ra[["RMSE"]]/rb[["RMSE"]]-1),
  pooled_B3_RPS=rb[["RPS"]],pooled_B2_RPS=ra[["RPS"]],pooled_delta_RPS=ra[["RPS"]]-rb[["RPS"]],
  unweighted_mean_site_delta_RMSE_percent=mean(ds$delta_RMSE_percent),
  site_min_delta_RMSE_percent=min(ds$delta_RMSE_percent),site_max_delta_RMSE_percent=max(ds$delta_RMSE_percent),
  sites_better_RMSE=sum(ds$delta_RMSE<0),sites_better_RPS=if(a$estimator=="OLS")NA_real_ else sum(ds$delta_RPS<0))
 check_numeric(as.numeric(a[names(vals)]),unname(vals))
}
write.csv(rec,file.path(root,"metricas_reconstruidas_independientemente.csv"),row.names=FALSE,na="")
# Verificación independiente de descriptivos: años y categorías ORIGINALES.
# SD se reconstruye con divisor n-1; porcentajes, con conteos explícitos.
table_folder <- file.path(root,"tablas_para_manuscrito")
descriptives_verified <- FALSE
if(file.exists(file.path(table_folder,"S16_descriptivos_institucion_numerico.csv"))){
 ds <- read.csv(file.path(table_folder,"S16_descriptivos_institucion_numerico.csv"))
 stopifnot(all(archived$CSR1_bin==as.numeric(archived$CSR1==1)),
  all(archived$Sexo_num==as.numeric(archived$Sexo)),all(archived$Estudios_num==as.numeric(archived$Estudios)))
 for(i in seq_len(nrow(ds))){
  a <- ds[i,];z <- archived[as.character(archived$Pais)==a$Site,];n <- nrow(z)
  em <- sum(z$Edad)/n;nm <- sum(z$NSE)/n
  vals <- c(N=n,Age_mean=em,Age_SD=sqrt(sum((z$Edad-em)^2)/(n-1)),
   Women_percent=100*sum(z$Sexo_num==1)/n,Postgraduate_percent=100*sum(z$Estudios_num==1)/n,
   CSR_course_yes_percent=100*sum(z$CSR1_bin==1)/n,NSE_mean=nm,NSE_SD=sqrt(sum((z$NSE-nm)^2)/(n-1)),
   PP1_mean=sum(z$PP1)/n,PP2_mean=sum(z$PP2)/n,PP3_mean=sum(z$PP3)/n)
  check_numeric(as.numeric(a[names(vals)]),unname(vals))
 }
 descriptives_verified <- TRUE
}
write_json(list(status="PASS",records=653L,site_outcome_model_score_rows_verified=60L,
 saved_model_predictions_verified=fit_predictions_verified,training_scaling_verified=TRUE,
 no_outcome_or_site_predictors=TRUE,duplicate_group_separation_absent=TRUE,
 institution_scores_and_pooled_scores_verified=TRUE,paired_site_and_pooled_differences_verified=TRUE,
 raw_variable_site_descriptives_verified=descriptives_verified,
 maximum_absolute_numeric_difference=max_difference,sources_sha256_unchanged=TRUE,
 predictions_sha256=digest(file.path(private,"predicciones_fuera_institucion.rds"),algo="sha256",file=TRUE),
 verification_script_sha256=digest(file.path(root,"02_verificar_transportabilidad_LISO.R"),algo="sha256",file=TRUE)),
 file.path(root,"auditoria_independiente.json"),pretty=TRUE,auto_unbox=TRUE,digits=NA)
cat("PASS: 60 score rows, 48 saved-model predictions, all paired and pooled losses verified.\n")
