# Verificación independiente desde objetos y predicciones, sin estimación.
# Recalcula RMSE, MAE y RPS; coteja cada caso y las asignaciones por grupos.
# Verifica que las sensibilidades conservaron las opciones y estructura M2.
options(stringsAsFactors=FALSE)
project <- Sys.getenv("RCP_PROYECTO", "/Users/luisbenites/Documents/Codex/Paper con MA - 1/reproduccion-integral-V3-2026-10-10/proyecto")
root <- file.path(project,"revision-3-revista-colombiana-psicologia-2026-10-08/analisis_sensibilidad")
.libPaths(c(file.path(project,"revision-revista-2026-10-04/auditoria-codigos/R-library"),.libPaths()))
suppressPackageStartupMessages({library(lavaan);library(jsonlite);library(digest)})
q <- readRDS(file.path(root,"privados/datos_y_banderas_para_sensibilidades.rds"))
reported <- read.csv(file.path(root,"grouped_cv/metricas_por_repeticion.csv"))
cv_checks <- list();reconstructed <- list();index <- 0L
for(s in names(q$scenarios)){
 base <- file.path(root,"grouped_cv",s,"privados")
 a <- readRDS(file.path(base,"asignaciones_por_grupo.rds"))
 p <- readRDS(file.path(base,"predicciones_fuera_muestra.rds"))
 stopifnot(identical(a$original_row,which(q$scenarios[[s]])),
  identical(a$group,q$group[q$scenarios[[s]]]),
  nrow(p)==5L*length(a$original_row)*3L*5L)
 observed <- vapply(seq_len(nrow(p)),function(i)q$data[p$original_row[i],p$item[i]],numeric(1))
 stopifnot(identical(as.numeric(observed),as.numeric(p$observed)))
 for(r in 1:5){
  assignment <- a$assignment[[r]]
  stopifnot(all(vapply(split(assignment,a$group),function(x)length(unique(x))==1L,logical(1))))
  for(k in 1:5){
   train <- which(assignment!=k);test <- which(assignment==k)
   stopifnot(length(intersect(a$group[train],a$group[test]))==0L)
   z <- p[p$repetition==r & p$fold==k,]
   stopifnot(setequal(unique(z$original_row),a$original_row[test]))
  }
 }
 ks <- unique(p[,c("repetition","item","model","estimator")])
 for(i in seq_len(nrow(ks))){
  key <- ks[i,];z <- p[p$repetition==key$repetition & p$item==key$item &
   p$model==key$model & p$estimator==key$estimator,]
  stopifnot(nrow(z)==length(a$original_row),!anyDuplicated(z$original_row))
  err <- z$predicted_expected-z$observed;rmse <- sqrt(mean(err^2));mae <- mean(abs(err));rps <- NA_real_
  if(key$estimator!="OLS"){
   prob <- as.matrix(z[,paste0("probability_",1:5)])
   stopifnot(all(is.finite(prob)),all(prob>=0),max(abs(rowSums(prob)-1))<1e-10,
    max(abs(as.vector(prob %*% (1:5))-z$predicted_expected))<1e-12)
   cumulative <- t(apply(prob,1,cumsum))[,1:4,drop=FALSE]
   truth <- sapply(1:4,function(cut)as.numeric(z$observed<=cut))
   rps <- mean(rowSums((cumulative-truth)^2)/4)
  }
  index <- index+1L
  reconstructed[[index]] <- cbind(data.frame(scenario=s,N=length(a$original_row)),key,
   data.frame(RMSE=rmse,MAE=mae,RPS=rps))
 }
 # Las envolturas polr prueban que la escala y los casos proceden exclusivamente
 # del entrenamiento, sin utilizar PP como predictor. No vuelve a ajustar.
 cached <- list.files(file.path(base,"fits"),pattern="\\.rds$",full.names=TRUE)
 stopifnot(length(cached)==150L)
 for(path in cached){
  fit <- readRDS(path)
  stopifnot(!any(paste0("PP",1:3) %in% fit$vars),length(intersect(q$group[fit$train],q$group[fit$test]))==0L)
  xt <- q$data[fit$train,fit$vars,drop=FALSE]
  stopifnot(max(abs(vapply(xt,mean,numeric(1))-fit$center))<1e-12,
   max(abs(vapply(xt,sd,numeric(1))-fit$scale))<1e-12,
   fit$fit$convergence==0L,all(is.finite(c(fit$fit$coefficients,fit$fit$zeta))))
 }
 cv_checks[[s]] <- list(N=length(a$original_row),groups=length(unique(a$group)),
  private_prediction_rows=nrow(p),probit_fits_inspected=length(cached),
  groups_split=0L,training_only_scaling_verified=TRUE,outcome_leakage=FALSE)
}
computed <- do.call(rbind,reconstructed)
key <- function(x)paste(x$scenario,x$repetition,x$item,x$model,x$estimator)
ix <- match(key(reported),key(computed));stopifnot(!anyNA(ix),nrow(computed)==nrow(reported))
diff <- vapply(c("RMSE","MAE","RPS"),function(v){
 a <- reported[[v]];b <- computed[[v]][ix];stopifnot(identical(is.na(a),is.na(b)))
 max(abs(a-b),na.rm=TRUE)
},numeric(1))
stopifnot(all(diff<1e-12))
write.csv(computed,file.path(root,"grouped_cv/agregados_reconstruidos_independientemente.csv"),row.names=FALSE,na="")
model_checks <- list()
for(s in names(q$scenarios))for(id in c("CFA_MLR","CFA_WLSMV","SEM_MLR","SEM_WLSMV","ITEMS_MLR","ITEMS_WLSMV")){
 f <- readRDS(file.path(root,"privados",paste0(s,"__",id,".rds")))
 ref <- readRDS(file.path(root,"privados",paste0("all__",id,".rds")))
 opts <- c("estimator.orig","std.lv","fixed.x","conditional.x","parameterization","meanstructure","missing")
 stopifnot(identical(lavInspect(f,"options")[opts],lavInspect(ref,"options")[opts]),
  identical(as.character(f@call$model),as.character(ref@call$model)),
  lavInspect(f,"nobs")==sum(q$scenarios[[s]]),lavInspect(f,"converged"),lavInspect(f,"post.check"))
 model_checks[[paste(s,id)]] <- list(N=lavInspect(f,"nobs"),unchanged_model_and_options=TRUE)
}
pa <- read.csv(file.path(root,"modelos_centrales/rutas_focales_sensibilidad.csv"))
for(s in unique(pa$scenario))for(id in unique(pa$model)){
 z <- pa[pa$scenario==s & pa$model==id,]
 expected <- if(grepl("^SEM",id))21L else 18L
 stopifnot(nrow(z)==expected,all(z$Holm_family_size==expected),
  max(abs(p.adjust(z$pvalue,"holm")-z$p_Holm_family))<1e-12)
}
report <- list(status="PASS",new_fits=0L,maximum_reconstructed_metric_differences=as.list(diff),
 CV_scenarios=cv_checks,model_specification_checks=model_checks,Holm_families_verified=c(21L,18L),
 rows_in_reconstructed_repetition_table=nrow(computed),
 limitation="Internal retrospective exploratory validation; grouped identical records do not establish participant identity; no external validation or tests treating repeated folds as independent")
write_json(report,file.path(root,"auditoria_independiente_agregados_y_modelos.json"),pretty=TRUE,auto_unbox=TRUE,digits=NA)
cat("INDEPENDENT CHECKS PASS\n");print(diff)
