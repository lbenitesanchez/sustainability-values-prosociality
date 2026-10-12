# VERSIÓN 4 PARA ARLA: SENSIBILIDAD DE TRANSPORTABILIDAD ENTRE INSTITUCIONES.
# Entrada recomendada: Rscript 01_transportabilidad_LISO.R
# Se fija el protocolo en protocolo_antes_ajustes.json ANTES de este script.
# No cambia datos, modelos SEM ni resultados de versiones anteriores.
# Pais identifica cuatro sitios muestrales (una institución por país): no permite
# separar diferencias institucionales de diferencias nacionales.
# Las cuatro particiones entrenan en tres sitios y evalúan el sitio restante.
# No se calculan pruebas ni intervalos con cuatro instituciones como clústeres.
options(stringsAsFactors=FALSE, warn=1)
Sys.setenv(VECLIB_MAXIMUM_THREADS="1", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
project <- Sys.getenv("ARLA_PROYECTO", Sys.getenv("RCP_PROYECTO", "/Users/luisbenites/Documents/Codex/Paper con MA - 1/reproduccion-integral-V3-2026-10-10/proyecto"))
root <- file.path(project,"revision-4-ARLA-2026-10-09/transportabilidad_institucional")
private <- file.path(root,"privados");dir.create(file.path(private,"fits"),recursive=TRUE,showWarnings=FALSE)
.libPaths(c(file.path(project,"revision-revista-2026-10-04/auditoria-codigos/R-library"),.libPaths()))
suppressPackageStartupMessages({library(MASS);library(jsonlite);library(digest)})
state <- file.path(root,"estado.txt");writeLines("RUNNING",state)
options(error=function(){writeLines(c("FAILED",geterrmessage()),state);quit(status=1,save="no")})
stopifnot(file.exists(file.path(root,"protocolo_antes_ajustes.json")))
paths <- c(quality=file.path(project,"revision-3-revista-colombiana-psicologia-2026-10-08/analisis_sensibilidad/privados/datos_y_banderas_para_sensibilidades.rds"),
 archived=file.path(project,"revision-integral-2026-10-04/analysis/outputs/canonical-full-20261004-02/01_audit/analytic_data.rds"))
sha <- function(p) digest(p,algo="sha256",file=TRUE)
hash_before <- vapply(paths,sha,character(1))
q <- readRDS(paths[["quality"]]);archived <- as.data.frame(readRDS(paths[["archived"]]))
data <- q$data
stopifnot(nrow(data)==653L,nrow(archived)==653L,identical(q$original_row,seq_len(653L)))
# Edad y NSE originales evitan usar constantes de estandarización del conjunto
# completo. Ambas se centran y escalan NUEVAMENTE dentro de cada entrenamiento.
data$Edad <- archived$Edad;data$NSE <- archived$NSE
items <- paste0("PP",1:3);categories <- 1:5
aps <- c("APS1","APS2","APS3","APS4","APS5","APS6","APS7","APS9","APS10","APS14",
 "APS11","APS12","APS13","APS15","APS16","APS18","APS19","APS20","APS21")
controls <- c("Sexo_num","Estudios_num","Edad","NSE")
values <- c("Benevolencia","Universalismo","Poder","Logro")
blocks <- list(B3_APS_without_values=c("CSR1_bin",controls,aps),
 B2_add_values=c("CSR1_bin",controls,aps,values))
common <- unique(c(items,unlist(blocks)))
stopifnot(all(common %in% names(archived)),
 max(abs(as.matrix(data[,common])-as.matrix(archived[,common])))==0,
 !anyNA(data[,common]),!any(items %in% unlist(blocks)),
 all(as.matrix(data[,items]) %in% categories))
sites <- c("Argentina","Colombia","Peru","Mexico")
site <- as.character(data$Pais)
stopifnot(setequal(unique(site),sites),all(vapply(split(site,q$group),function(z)length(unique(z))==1L,logical(1))))
sample <- do.call(rbind,lapply(sites,function(s){
 train <- which(site!=s);test <- which(site==s)
 data.frame(held_out_site=s,N=653L,n_train=length(train),n_test=length(test),
  identical_groups_crossing_split=length(intersect(q$group[train],q$group[test])),
  n_test_duplicate_members=sum(q$dup_member[test]),n_test_subsequent_copies=sum(q$dup_copy[test]),
  n_train_duplicate_members=sum(q$dup_member[train]),n_train_subsequent_copies=sum(q$dup_copy[train]))
}))
stopifnot(all(sample$identical_groups_crossing_split==0L),sum(sample$n_test)==653L)
write.csv(sample,file.path(root,"particiones_y_duplicados.csv"),row.names=FALSE)
write_json(list(status="INPUTS_VERIFIED_BEFORE_FITS",recorded_at=format(Sys.time(),tz="UTC",usetz=TRUE),
 sources_sha256=as.list(hash_before),script_sha256=sha(file.path(root,"01_transportabilidad_LISO.R")),
 protocol_sha256=sha(file.path(root,"protocolo_antes_ajustes.json")),blocks=blocks,
 N=653L,site_sizes=as.list(setNames(sample$n_test,sample$held_out_site)),
 raw_age_and_NSE_verified=TRUE,no_identical_groups_crossing=TRUE),
 file.path(root,"auditoria_entradas.json"),pretty=TRUE,auto_unbox=TRUE,digits=NA)

# Esta función corresponde a las columnas de error de la tabla LISO/S14.
# RPS se normaliza por cuatro umbrales y se interpreta como pérdida, menor mejor.
metrics <- function(y,expected,p=NULL){
 n <- length(y);err <- expected-y;ok <- length(expected)==n && all(is.finite(expected))
 x <- list(n_expected=n,n_scored=if(ok)n else 0L,
  SSE=if(ok)sum(err^2)else NA_real_,SAE=if(ok)sum(abs(err))else NA_real_,
  RMSE=if(ok)sqrt(mean(err^2))else NA_real_,MAE=if(ok)mean(abs(err))else NA_real_,
  sum_RPS=NA_real_,RPS=NA_real_,sum_logloss_raw=NA_real_,logloss_raw=NA_real_,
  sum_logloss_floor=NA_real_,logloss_floor=NA_real_,n_zero_probability=NA_integer_)
 if(!is.null(p) && all(is.finite(p))){
  stopifnot(nrow(p)==n,ncol(p)==5L,all(p>=0),max(abs(rowSums(p)-1))<1e-10)
  cp <- t(apply(p,1,cumsum))[,1:4,drop=FALSE]
  truth <- sapply(1:4,function(cut)as.numeric(y<=cut))
  rps <- rowSums((cp-truth)^2)/4;obs <- p[cbind(seq_len(n),y)]
  x$sum_RPS <- sum(rps);x$RPS <- mean(rps)
  x$sum_logloss_raw <- sum(-log(obs));x$logloss_raw <- mean(-log(obs))
  x$sum_logloss_floor <- sum(-log(pmax(obs,1e-12)));x$logloss_floor <- mean(-log(pmax(obs,1e-12)))
  x$n_zero_probability <- sum(obs==0)
 };x
}
rows <- predictions <- frequency <- list();i <- 0L
record <- function(s,item,id,estimator,train,test,expected,p,diagnostics){
 i <<- i+1L
 rows[[i]] <<- cbind(data.frame(held_out_site=s,item=item,model=id,estimator=estimator,n_train=length(train)),
  as.data.frame(metrics(data[test,item],expected,p)),as.data.frame(diagnostics))
 z <- data.frame(original_row=q$original_row[test],held_out_site=s,item=item,model=id,
  estimator=estimator,observed=data[test,item],predicted_expected=expected)
 for(c in categories)z[[paste0("probability_",c)]] <- if(is.null(p))NA_real_ else p[,c]
 predictions[[i]] <<- z
}
for(s in sites){
 train <- which(site!=s);test <- which(site==s)
 stopifnot(length(intersect(q$group[train],q$group[test]))==0L)
 for(item in items){
  counts <- table(factor(data[train,item],levels=categories));present <- categories[counts>0]
  frequency[[paste(s,item,sep="_")]] <- data.frame(held_out_site=s,item=item,category=categories,
   n_train=as.integer(counts),n_test=as.integer(table(factor(data[test,item],levels=categories))))
  bp <- matrix(rep(as.numeric(counts)/length(train),each=length(test)),nrow=length(test))
  empty <- list(converged=TRUE,coefficients_finite=TRUE,thresholds_ordered=TRUE,
   min_hessian_eigen=NA_real_,hessian_condition=NA_real_,initialization_rescued=FALSE,
   missing_training_categories=paste(setdiff(categories,present),collapse=";"),warnings="",error="")
  record(s,item,"B_empirical","empirical",train,test,as.vector(bp %*% categories),bp,empty)
  for(id in names(blocks)){
   vars <- blocks[[id]];xt <- data[train,vars,drop=FALSE];xv <- data[test,vars,drop=FALSE]
   center <- vapply(xt,mean,numeric(1));scale <- vapply(xt,sd,numeric(1))
   stopifnot(all(is.finite(scale)),all(scale>0))
   zt <- as.data.frame(sweep(sweep(as.matrix(xt),2,center,"-"),2,scale,"/"))
   zv <- as.data.frame(sweep(sweep(as.matrix(xv),2,center,"-"),2,scale,"/"))
   training <- data.frame(response=ordered(data[train,item],levels=present),zt)
   warns <- character();rescued <- FALSE
   fit <- tryCatch(withCallingHandlers(polr(response ~ .,data=training,method="probit",Hess=TRUE,
    model=TRUE,control=list(maxit=5000L,reltol=1e-10)),warning=function(w){
     warns <<- c(warns,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
   if(inherits(fit,"error")){
    permitted <- grepl("attempt to find suitable starting values failed",conditionMessage(fit),fixed=TRUE) ||
     grepl("initial value in 'vmmin' is not finite",conditionMessage(fit),fixed=TRUE)
    if(permitted && length(present)>=3L){
     starts <- c(rep(0,length(vars)),qnorm(head(cumsum(as.numeric(counts[counts>0]))/length(train),-1L)))
     rescued <- TRUE
     fit <- tryCatch(withCallingHandlers(polr(response ~ .,data=training,method="probit",Hess=TRUE,
      model=TRUE,start=starts,control=list(maxit=5000L,reltol=1e-10)),warning=function(w){
       warns <<- c(warns,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
    }
   }
   di <- empty;di$warnings <- paste(unique(warns),collapse=" | ");di$initialization_rescued <- rescued
   expected <- rep(NA_real_,length(test));p <- matrix(NA_real_,length(test),5L)
   if(inherits(fit,"error")){
    di$error <- conditionMessage(fit);di$converged <- FALSE;di$coefficients_finite <- FALSE;di$thresholds_ordered <- FALSE
   }else{
    he <- tryCatch(eigen((fit$Hessian+t(fit$Hessian))/2,symmetric=TRUE,only.values=TRUE)$values,error=function(e)NA_real_)
    di$converged <- fit$convergence==0L;di$coefficients_finite <- all(is.finite(c(fit$coefficients,fit$zeta)))
    di$thresholds_ordered <- all(diff(fit$zeta)>0);di$min_hessian_eigen <- min(he)
    di$hessian_condition <- max(he)/min(he)
    valid <- di$converged && di$coefficients_finite && di$thresholds_ordered && all(is.finite(he)) && min(he)>0
    if(valid){
     fp <- tryCatch(as.matrix(predict(fit,newdata=zv,type="probs")),error=function(e)e)
     if(inherits(fp,"error")){di$error <- conditionMessage(fp)}else{
      p <- matrix(0,nrow=length(test),ncol=5L,dimnames=list(NULL,as.character(categories)))
      p[,colnames(fp)] <- fp;expected <- as.vector(p %*% categories)
     }
    }else di$error <- "Numerical validity criterion failed; predictions retained as unavailable"
   }
   record(s,item,id,"ordinal_probit",train,test,expected,p,di)
   saveRDS(list(fit=fit,train=q$original_row[train],test=q$original_row[test],center=center,scale=scale,
    vars=vars,warnings=warns,rescued=rescued,diagnostics=di),
    file.path(private,"fits",paste(s,item,id,"probit.rds",sep="_")))
   # OLS es otra sensibilidad: no se recortan valores fuera de 1-5 y no produce
   # probabilidades ordinales. Por ello RPS/entropía no se atribuyen a OLS.
   ols <- tryCatch(lm(response ~ .,data=data.frame(response=data[train,item],zt),na.action=na.fail),error=function(e)e)
   od <- empty;oe <- rep(NA_real_,length(test))
   if(inherits(ols,"error")){
    od$error <- conditionMessage(ols);od$converged <- FALSE;od$coefficients_finite <- FALSE
   }else{
    od$coefficients_finite <- all(is.finite(coef(ols)))
    od$converged <- ols$rank==length(vars)+1L && od$coefficients_finite
    if(od$converged)oe <- as.numeric(predict(ols,newdata=zv))else od$error <- "OLS rank or finite-coefficient criterion failed"
   }
   record(s,item,id,"OLS",train,test,oe,NULL,od)
   saveRDS(list(fit=ols,train=q$original_row[train],test=q$original_row[test],center=center,scale=scale,
    vars=vars,diagnostics=od),file.path(private,"fits",paste(s,item,id,"OLS.rds",sep="_")))
  }
 }
 cat("LISO",s,"complete\n");flush.console()
}
ft <- do.call(rbind,rows);pr <- do.call(rbind,predictions)
stopifnot(nrow(ft)==60L,nrow(pr)==653L*3L*5L)
write.csv(ft,file.path(root,"metricas_y_diagnosticos_por_institucion.csv"),row.names=FALSE,na="")
write.csv(do.call(rbind,frequency),file.path(root,"frecuencias_train_test.csv"),row.names=FALSE)
saveRDS(pr,file.path(private,"predicciones_fuera_institucion.rds"))
# CONTRASTES PAREADOS POR SITIO: bloques B2 y B3 evaluados sobre los mismos casos.
paired <- do.call(rbind,lapply(sites,function(s)do.call(rbind,lapply(items,function(it)
 do.call(rbind,lapply(c("ordinal_probit","OLS"),function(e){
  a <- ft[ft$held_out_site==s & ft$item==it & ft$estimator==e & ft$model=="B2_add_values",]
  b <- ft[ft$held_out_site==s & ft$item==it & ft$estimator==e & ft$model=="B3_APS_without_values",]
  stopifnot(nrow(a)==1L,nrow(b)==1L,a$n_expected==b$n_expected)
  data.frame(held_out_site=s,item=it,estimator=e,n_test=a$n_expected,
   B3_RMSE=b$RMSE,B2_RMSE=a$RMSE,delta_RMSE=a$RMSE-b$RMSE,
   delta_RMSE_percent=100*(a$RMSE/b$RMSE-1),B3_MAE=b$MAE,B2_MAE=a$MAE,
   delta_MAE=a$MAE-b$MAE,B3_RPS=b$RPS,B2_RPS=a$RPS,delta_RPS=a$RPS-b$RPS,
   both_models_scored=all(c(a$n_scored,b$n_scored)==a$n_expected))
 }))))))
write.csv(paired,file.path(root,"contrastes_pareados_por_institucion.csv"),row.names=FALSE,na="")
# RESUMEN: N pondera las pérdidas de cada registro una sola vez. La media no
# ponderada entre cuatro sitios se registra aparte, sin SE, intervalos ni p.
keys <- unique(ft[,c("item","model","estimator")])
pooled <- do.call(rbind,lapply(seq_len(nrow(keys)),function(i){
 key <- keys[i,];z <- ft[ft$item==key$item & ft$model==key$model & ft$estimator==key$estimator,]
 ok <- nrow(z)==4L && sum(z$n_scored)==653L
 cbind(key,data.frame(n_sites=4L,N_expected=653L,N_scored=sum(z$n_scored),
  pooled_RMSE=if(ok)sqrt(sum(z$SSE)/653)else NA_real_,
  pooled_MAE=if(ok)sum(z$SAE)/653 else NA_real_,
  pooled_RPS=if(ok && key$estimator!="OLS")sum(z$sum_RPS)/653 else NA_real_,
  pooled_logloss_raw=if(ok && key$estimator!="OLS")sum(z$sum_logloss_raw)/653 else NA_real_,
  pooled_logloss_floor=if(ok && key$estimator!="OLS")sum(z$sum_logloss_floor)/653 else NA_real_,
  unweighted_site_mean_RMSE=if(ok)mean(z$RMSE)else NA_real_,
  unweighted_site_mean_MAE=if(ok)mean(z$MAE)else NA_real_,
  unweighted_site_mean_RPS=if(ok && key$estimator!="OLS")mean(z$RPS)else NA_real_))
}))
write.csv(pooled,file.path(root,"resumen_errores_agrupados.csv"),row.names=FALSE,na="")
summary_pairs <- do.call(rbind,lapply(items,function(it)do.call(rbind,lapply(c("ordinal_probit","OLS"),function(e){
 a <- pooled[pooled$item==it & pooled$estimator==e & pooled$model=="B2_add_values",]
 b <- pooled[pooled$item==it & pooled$estimator==e & pooled$model=="B3_APS_without_values",]
 z <- paired[paired$item==it & paired$estimator==e,];ok <- all(z$both_models_scored)
 data.frame(item=it,estimator=e,N=653L,sites_scored=sum(z$both_models_scored),
  pooled_B3_RMSE=b$pooled_RMSE,pooled_B2_RMSE=a$pooled_RMSE,
  pooled_delta_RMSE_percent=100*(a$pooled_RMSE/b$pooled_RMSE-1),
  pooled_B3_RPS=b$pooled_RPS,pooled_B2_RPS=a$pooled_RPS,pooled_delta_RPS=a$pooled_RPS-b$pooled_RPS,
  unweighted_mean_site_delta_RMSE_percent=if(ok)mean(z$delta_RMSE_percent)else NA_real_,
  site_min_delta_RMSE_percent=if(ok)min(z$delta_RMSE_percent)else NA_real_,
  site_max_delta_RMSE_percent=if(ok)max(z$delta_RMSE_percent)else NA_real_,
  sites_better_RMSE=if(ok)sum(z$delta_RMSE<0)else NA_integer_,
  sites_better_RPS=if(ok && e!="OLS")sum(z$delta_RPS<0)else NA_integer_)
}))))
write.csv(summary_pairs,file.path(root,"resumen_incremento_predictivo.csv"),row.names=FALSE,na="")
stopifnot(identical(hash_before,vapply(paths,sha,character(1))))
ordinal <- ft[ft$estimator=="ordinal_probit",]
write_json(list(status=if(all(ft$n_scored==ft$n_expected))"COMPLETED_ALL_REQUESTED_RECORDS_SCORED"else"COMPLETED_WITH_EXPLICIT_FAILURES",
 N=653L,n_sites=4L,n_probit_attempts=24L,n_OLS_attempts=24L,n_empirical_references=12L,
 site_country_confounded=TRUE,country_indicators_excluded=TRUE,training_only_scaling=TRUE,
 identical_groups_crossing_splits=0L,no_omitted_sites=TRUE,no_omitted_cases=TRUE,
 all_models_scored=all(ft$n_scored==ft$n_expected),all_numerically_valid=all(ft$converged & ft$coefficients_finite & ft$thresholds_ordered),
 ordinal_minimum_hessian_eigen=min(ordinal$min_hessian_eigen,na.rm=TRUE),
 ordinal_initializer_rescues=sum(ordinal$initialization_rescued),ordinal_fits_with_warnings=sum(nzchar(ordinal$warnings)),
 absent_training_categories=sum(do.call(rbind,frequency)$n_train==0),
 total_zero_probability_evaluations=sum(ft$n_zero_probability,na.rm=TRUE),
 sources_unchanged=TRUE,sources_sha256=as.list(hash_before),
 limitations=c("Four sampled sites only; institutional and country effects inseparable", "Common probit slopes not established by these analyses", "Predictive sensitivity does not establish measurement invariance", "No four-cluster uncertainty or independent-fold tests")),
 file.path(root,"auditoria_final.json"),pretty=TRUE,auto_unbox=TRUE,digits=NA)
writeLines(capture.output(sessionInfo()),file.path(root,"sessionInfo.txt"))
writeLines("COMPLETED",state);cat("ALL LISO FITS AND AGGREGATES COMPLETED\n")
