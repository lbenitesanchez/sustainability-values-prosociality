# REVISIÓN 3: B2/B3 predictivos con conjuntos de registros idénticos agrupados.
# Ejecutar DESPUÉS de 01_modelos_centrales_sensibilidad.R (o cuando ya haya
# creado privados/datos_y_banderas_para_sensibilidades.rds).
# B2: CSR1 + siete controles + cuatro valores + 19 APS.
# B3: mismos predictores salvo los cuatro valores.
# Cada escenario utiliza 5 folds x 5 repeticiones propios. Dentro del escenario
# los folds son idénticos para B2, B3, probit y OLS. No se interpretan diferencias
# entre escenarios como efectos de eliminar casos, ni se prueban folds como iid.
# Se conservan categorías raras y todos los casos definidos por cada escenario.
options(stringsAsFactors=FALSE,warn=1)
Sys.setenv(VECLIB_MAXIMUM_THREADS="1",OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1")
project <- Sys.getenv("RCP_PROYECTO", "/Users/luisbenites/Documents/Codex/Paper con MA - 1/reproduccion-integral-V3-2026-10-10/proyecto")
root <- file.path(project,"revision-3-revista-colombiana-psicologia-2026-10-08/analisis_sensibilidad")
out <- file.path(root,"grouped_cv")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
.libPaths(c(file.path(project,"revision-revista-2026-10-04/auditoria-codigos/R-library"),.libPaths()))
suppressPackageStartupMessages({library(MASS);library(jsonlite);library(digest)})
state <- file.path(out,"estado.txt");writeLines("RUNNING",state)
options(error=function(){writeLines(c("FAILED",geterrmessage()),state);quit(status=1,save="no")})
quality_path <- file.path(root,"privados/datos_y_banderas_para_sensibilidades.rds")
q <- readRDS(quality_path); original_hash <- digest(quality_path,algo="sha256",file=TRUE)
items <- paste0("PP",1:3);categories <- 1:5
sc <- c("APS1","APS2","APS3","APS4","APS5","APS6","APS7","APS9","APS10","APS14")
sbl <- c("APS11","APS12","APS13","APS15","APS16","APS18","APS19","APS20","APS21")
controls <- c("C_Colombia","C_Peru","C_Mexico","Sexo_num","Estudios_num","Edad_z","NSE_z")
values <- c("Benevolencia","Universalismo","Poder","Logro")
blocks <- list(B2_add_APS=c("CSR1_bin",controls,values,sc,sbl),
 B3_APS_without_values=c("CSR1_bin",controls,sc,sbl))
stopifnot(!any(items %in% unlist(blocks)),nrow(q$data)==653L)
folds <- 5L;repeats <- 5L;seed <- 20261008L;logloss_floor <- 1e-12

# Se aleatorizan los grupos completos dentro de cada institución/país; se
# asigna cada grupo al fold con menos registros, con empate aleatorizado.
# La posible identidad de participantes NO se conoce: es agrupación precautoria.
make_assignment <- function(data,group,scenario_index){
 assignments <- vector("list",repeats)
 for(r in seq_len(repeats)){
  set.seed(seed+1000L*scenario_index+r);a <- integer(nrow(data))
  for(s in unique(data$Pais)){
   ix <- which(data$Pais==s); ids <- unique(group[ix]); ids <- sample(ids,length(ids))
   sizes <- integer(folds)
   for(g in ids){
    rows <- ix[group[ix]==g];candidates <- which(sizes==min(sizes))
    k <- candidates[sample.int(length(candidates),1L)]
    a[rows] <- k;sizes[k] <- sizes[k]+length(rows)
   }
  }
  stopifnot(all(a %in% seq_len(folds)),all(vapply(split(a,group),function(z)length(unique(z))==1L,logical(1))))
  assignments[[r]] <- a
 }
 assignments
}
all_folds <- all_repeats <- all_summary <- all_pairs <- all_pair_summary <- list()
diagnostics <- list(); assignments_all <- list(); fold_counts <- list()
for(si in seq_along(q$scenarios)){
 scenario <- names(q$scenarios)[si];keep <- q$scenarios[[si]]
 data <- q$data[keep,,drop=FALSE];group <- q$group[keep];original_row <- q$original_row[keep]
 assignment <- make_assignment(data,group,si);assignments_all[[scenario]] <- assignment
 output <- file.path(out,scenario);private <- file.path(output,"privados")
 dir.create(file.path(private,"fits"),recursive=TRUE,showWarnings=FALSE)
 saveRDS(list(assignment=assignment,group=group,original_row=original_row),file.path(private,"asignaciones_por_grupo.rds"))
 training_counts <- do.call(rbind,lapply(seq_len(repeats),function(r)do.call(rbind,
  lapply(seq_len(folds),function(k)do.call(rbind,lapply(items,function(item)
   data.frame(scenario=scenario,repetition=r,fold=k,item=item,category=categories,
    n_train=as.integer(table(factor(data[assignment[[r]]!=k,item],levels=categories))),
    n_test=as.integer(table(factor(data[assignment[[r]]==k,item],levels=categories))))))))))
 write.csv(training_counts,file.path(output,"frecuencias_por_pliegue.csv"),row.names=FALSE)
 fold_counts[[scenario]] <- do.call(rbind,lapply(seq_len(repeats),function(r)do.call(rbind,
  lapply(seq_len(folds),function(k)data.frame(scenario=scenario,repetition=r,fold=k,
   N=nrow(data),n_test=sum(assignment[[r]]==k),n_train=sum(assignment[[r]]!=k),
   n_test_clusters=length(unique(group[assignment[[r]]==k])),identical_sets_split=0L)))))
 protocol <- list(status="SPECIFIED_BEFORE_NEW_FITS",recorded_at=format(Sys.time(),tz="UTC",usetz=TRUE),
  scenario=scenario,N=nrow(data),number_of_clusters=length(unique(group)),blocks=blocks,
  estimators=c("ordinal cumulative probit parallel slopes","OLS","empirical category reference"),
  folds=folds,repetitions=repeats,seed_base=seed,scenario_seed_offset=1000L*si,
  grouping="Identical on all 99 original fields except Codigo; unknown participant identity; keep whole group in one fold",
  country_stratification=TRUE,training_only_scaling=TRUE,no_outcome_predictors=TRUE,
  no_selection_tuning_or_regularization=TRUE,
  initialization_rescue="Default polr initializer failures (suitable starting values failed or initial vmmin value not finite): zero slopes and training empirical CDF probit thresholds; same likelihood and optimization",
  categories="Fit support present in training, map to original 1-5 with zero for absent categories; no case/category deletion",
  scoring="RMSE/MAE on expected original categories; normalized RPS over four cutpoints; raw logloss Inf for zero probability; 1e-12 floor diagnostic only",
  parallel_slopes="Fixed cumulative-probit benchmark assumes common slopes across cumulative category thresholds; not checked as a confirmed substantive measurement model",
  no_tests_on_dependent_folds=TRUE,scenario_comparisons="Different retained cases and folds; do not interpret cross-scenario absolute errors as exclusion effects",
  original_quality_sha256=original_hash)
 protocol_path <- file.path(output,"protocolo_antes_ajustes.json")
 if(file.exists(protocol_path)){
  protocol$resumption_note <- "Preserve original protocol and completed fits; add deterministic zero-slope training-CDF start for specific nonfinite vmmin initialization error only, without changing likelihood, data or partitions"
  protocol_path <- file.path(output,"protocolo_reanudacion_inicializacion.json")
 }
 write_json(protocol,protocol_path,pretty=TRUE,auto_unbox=TRUE,digits=NA)
 metrics <- function(y,expected,p=NULL){
  err <- expected-y;n <- length(y)
  x <- list(n_test=n,sum_squared_error=sum(err^2),sum_absolute_error=sum(abs(err)),
   RMSE=sqrt(mean(err^2)),MAE=mean(abs(err)),sum_RPS=NA_real_,RPS=NA_real_,
   sum_logloss_raw=NA_real_,logloss_raw=NA_real_,sum_logloss_floor=NA_real_,logloss_floor=NA_real_,
   n_probability_zero=NA_integer_)
  if(!is.null(p)){
   stopifnot(nrow(p)==n,ncol(p)==5L,all(is.finite(p)),all(p>=0),max(abs(rowSums(p)-1))<1e-10)
   obs <- p[cbind(seq_len(n),y)]
   cp <- t(apply(p,1,cumsum))[,1:4,drop=FALSE]
   truth <- sapply(1:4,function(cut)as.numeric(y<=cut));rps <- rowSums((cp-truth)^2)/4
   x$sum_RPS <- sum(rps);x$RPS <- mean(rps)
   x$sum_logloss_raw <- sum(-log(obs));x$logloss_raw <- x$sum_logloss_raw/n
   x$sum_logloss_floor <- sum(-log(pmax(obs,logloss_floor)));x$logloss_floor <- x$sum_logloss_floor/n
   x$n_probability_zero <- sum(obs==0)
  };x
 }
 fold_records <- list();predictions <- list();index <- 0L
 record <- function(r,k,item,id,estimator,test,expected,p,diag){
  index <<- index+1L
  fold_records[[index]] <<- cbind(data.frame(scenario=scenario,repetition=r,fold=k,item=item,
   model=id,estimator=estimator,N=nrow(data),n_train=nrow(data)-length(test)),
   as.data.frame(metrics(data[test,item],expected,p)),as.data.frame(diag))
  z <- data.frame(original_row=original_row[test],repetition=r,fold=k,item=item,model=id,
   estimator=estimator,observed=data[test,item],predicted_expected=expected)
  for(c in categories)z[[paste0("probability_",c)]] <- if(is.null(p))NA_real_ else p[,c]
  predictions[[index]] <<- z
 }
 for(r in seq_len(repeats))for(k in seq_len(folds)){
  test <- which(assignment[[r]]==k);train <- which(assignment[[r]]!=k)
  stopifnot(length(intersect(group[train],group[test]))==0L)
  for(item in items){
   counts <- table(factor(data[train,item],levels=categories));present <- categories[counts>0]
   bp <- matrix(rep(as.numeric(counts)/length(train),each=length(test)),nrow=length(test))
   empty <- list(converged=TRUE,coefficients_finite=TRUE,threshold_ordered=TRUE,
    min_hessian_eigen=NA_real_,hessian_condition_number=NA_real_,initialization_rescued=FALSE,
    missing_training_categories=paste(setdiff(categories,present),collapse=";"),warnings="")
   record(r,k,item,"B_empirical","empirical",test,as.vector(bp %*% categories),bp,empty)
   for(id in names(blocks)){
    vars <- blocks[[id]];xt <- data[train,vars,drop=FALSE];xv <- data[test,vars,drop=FALSE]
    center <- vapply(xt,mean,numeric(1));scale <- vapply(xt,sd,numeric(1))
    stopifnot(all(is.finite(scale)),all(scale>0))
    zt <- as.data.frame(sweep(sweep(as.matrix(xt),2,center,"-"),2,scale,"/"))
    zv <- as.data.frame(sweep(sweep(as.matrix(xv),2,center,"-"),2,scale,"/"))
    training <- data.frame(response=ordered(data[train,item],levels=present),zt)
    path <- file.path(private,"fits",sprintf("%s_%s_probit_r%02d_f%02d.rds",item,id,r,k))
    warns <- character();rescued <- FALSE
    if(file.exists(path)){
     cached <- readRDS(path)
     stopifnot(identical(cached$train,original_row[train]),identical(cached$test,original_row[test]),
      identical(cached$center,center),identical(cached$scale,scale),identical(cached$vars,vars))
     fit <- cached$fit;warns <- cached$warnings;rescued <- cached$rescued
    }else{
     fit <- tryCatch(withCallingHandlers(polr(response ~ .,data=training,method="probit",Hess=TRUE,
      model=TRUE,control=list(maxit=5000L,reltol=1e-10)),warning=function(w){
       warns <<- c(warns,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
     if(inherits(fit,"error")){
      allowed_initializer_error <- grepl("attempt to find suitable starting values failed",conditionMessage(fit),fixed=TRUE) ||
       grepl("initial value in 'vmmin' is not finite",conditionMessage(fit),fixed=TRUE)
      if(!allowed_initializer_error)stop(fit)
      starts <- c(rep(0,length(vars)),qnorm(head(cumsum(as.numeric(counts[counts>0]))/length(train),-1L)))
      stopifnot(all(is.finite(starts)));rescued <- TRUE
      fit <- withCallingHandlers(polr(response ~ .,data=training,method="probit",Hess=TRUE,
       model=TRUE,start=starts,control=list(maxit=5000L,reltol=1e-10)),warning=function(w){
        warns <<- c(warns,conditionMessage(w));invokeRestart("muffleWarning")})
     }
    }
    he <- eigen((fit$Hessian+t(fit$Hessian))/2,symmetric=TRUE,only.values=TRUE)$values
    stopifnot(fit$convergence==0L,all(is.finite(c(fit$coefficients,fit$zeta))),
     all(diff(fit$zeta)>0),all(is.finite(he)),min(he)>0)
    fp <- as.matrix(predict(fit,newdata=zv,type="probs"))
    p <- matrix(0,nrow=length(test),ncol=5L,dimnames=list(NULL,as.character(categories)))
    p[,colnames(fp)] <- fp
    di <- empty;di$min_hessian_eigen <- min(he);di$hessian_condition_number <- max(he)/min(he)
    di$warnings <- paste(unique(warns),collapse=" | ");di$initialization_rescued <- rescued
    record(r,k,item,id,"ordinal_probit",test,as.vector(p %*% categories),p,di)
    saveRDS(list(fit=fit,train=original_row[train],test=original_row[test],center=center,
     scale=scale,vars=vars,warnings=warns,rescued=rescued),path)
    ols <- lm(response ~ .,data=data.frame(response=data[train,item],zt),na.action=na.fail)
    stopifnot(ols$rank==length(vars)+1L,all(is.finite(coef(ols))))
    expected <- as.numeric(predict(ols,newdata=zv));stopifnot(all(is.finite(expected)))
    record(r,k,item,id,"OLS",test,expected,NULL,empty)
   }
  }
  cat("GROUPED CV",scenario,"repeat",r,"fold",k,"complete\n");flush.console()
 }
 ft <- do.call(rbind,fold_records);pr <- do.call(rbind,predictions)
 stopifnot(nrow(ft)==375L,nrow(pr)==repeats*nrow(data)*3L*5L)
 write.csv(ft,file.path(output,"metricas_y_diagnosticos_por_pliegue.csv"),row.names=FALSE)
 saveRDS(pr,file.path(private,"predicciones_fuera_muestra.rds"))
 keys <- unique(ft[,c("scenario","N","repetition","item","model","estimator")])
 rt <- do.call(rbind,lapply(seq_len(nrow(keys)),function(i){
  key <- keys[i,];z <- ft[ft$repetition==key$repetition & ft$item==key$item &
   ft$model==key$model & ft$estimator==key$estimator,];n <- sum(z$n_test)
  stopifnot(nrow(z)==folds,n==nrow(data))
  cbind(key,data.frame(RMSE=sqrt(sum(z$sum_squared_error)/n),MAE=sum(z$sum_absolute_error)/n,
   RPS=sum(z$sum_RPS)/n,logloss_raw=sum(z$sum_logloss_raw)/n,
   logloss_floor=sum(z$sum_logloss_floor)/n,n_probability_zero=sum(z$n_probability_zero)))
 }))
 keys <- unique(rt[,c("scenario","N","item","model","estimator")])
 su <- do.call(rbind,lapply(seq_len(nrow(keys)),function(i){
  key <- keys[i,];z <- rt[rt$item==key$item & rt$model==key$model & rt$estimator==key$estimator,]
  stopifnot(nrow(z)==repeats)
  cbind(key,data.frame(repetitions=repeats,RMSE=mean(z$RMSE),RMSE_partition_SD=sd(z$RMSE),
   MAE=mean(z$MAE),RPS=mean(z$RPS),RPS_partition_SD=sd(z$RPS),
   logloss_raw=mean(z$logloss_raw),logloss_floor=mean(z$logloss_floor)))
 }))
 pairs <- do.call(rbind,lapply(items,function(it)do.call(rbind,lapply(c("ordinal_probit","OLS"),function(e){
  a <- rt[rt$item==it & rt$estimator==e & rt$model=="B2_add_APS",]
  b <- rt[rt$item==it & rt$estimator==e & rt$model=="B3_APS_without_values",]
  stopifnot(identical(a$repetition,b$repetition))
  data.frame(scenario=scenario,N=nrow(data),item=it,estimator=e,repetition=a$repetition,
   difference_RMSE=a$RMSE-b$RMSE,relative_RMSE_change_percent=100*(a$RMSE/b$RMSE-1),
   difference_MAE=a$MAE-b$MAE,difference_RPS=a$RPS-b$RPS)
 }))))
 ps <- do.call(rbind,lapply(items,function(it)do.call(rbind,lapply(c("ordinal_probit","OLS"),function(e){
  z <- pairs[pairs$item==it & pairs$estimator==e,]
  data.frame(scenario=scenario,N=nrow(data),item=it,estimator=e,repetitions=repeats,
   mean_difference_RMSE=mean(z$difference_RMSE),SD_difference_RMSE=sd(z$difference_RMSE),
   min_difference_RMSE=min(z$difference_RMSE),max_difference_RMSE=max(z$difference_RMSE),
   mean_relative_RMSE_change_percent=mean(z$relative_RMSE_change_percent),
   mean_difference_MAE=mean(z$difference_MAE),mean_difference_RPS=mean(z$difference_RPS),
   repetitions_better_RMSE=sum(z$difference_RMSE<0),
   repetitions_better_RPS=if(e=="OLS")NA_integer_ else sum(z$difference_RPS<0))
 }))))
 all_folds[[scenario]] <- ft;all_repeats[[scenario]] <- rt;all_summary[[scenario]] <- su
 all_pairs[[scenario]] <- pairs;all_pair_summary[[scenario]] <- ps
 for(nm in c("all_folds","all_repeats","all_summary","all_pairs","all_pair_summary")){
  dest <- c(all_folds="metricas_por_pliegue_todos_escenarios.csv",all_repeats="metricas_por_repeticion.csv",
   all_summary="resumen_modelos_por_item.csv",all_pairs="diferencias_pareadas_por_repeticion.csv",
   all_pair_summary="resumen_diferencias_pareadas.csv")[[nm]]
  write.csv(do.call(rbind,get(nm)),file.path(out,dest),row.names=FALSE,na="")
 }
 diagnostics[[scenario]] <- list(N=nrow(data),clusters=length(unique(group)),n_probit_fits=150L,
  n_OLS_fits=150L,n_empirical_folds=75L,identical_sets_split=0L,
  all_converged=all(ft$converged),all_coefficients_finite=all(ft$coefficients_finite),
  all_thresholds_ordered=all(ft$threshold_ordered),minimum_Hessian_eigen=min(ft$min_hessian_eigen,na.rm=TRUE),
  initialization_rescues=sum(ft$initialization_rescued),fits_with_warnings=sum(nzchar(ft$warnings)),
  absent_training_categories=sum(training_counts$n_train==0),zero_scoring_probabilities=sum(ft$n_probability_zero,na.rm=TRUE),
  no_omitted_folds=TRUE,no_independence_tests=TRUE)
 write_json(diagnostics,file.path(out,"diagnosticos_escenarios.json"),pretty=TRUE,auto_unbox=TRUE,digits=NA)
}
stopifnot(identical(original_hash,digest(quality_path,algo="sha256",file=TRUE)))
write.csv(do.call(rbind,fold_counts),file.path(out,"tamanos_y_agrupacion_folds.csv"),row.names=FALSE)
write_json(list(status="COMPLETED",N_primary=653L,N_sensitivities=c(650L,618L,615L),
 n_probit_fits=600L,n_OLS_fits=600L,identical_sets_crossing_folds=0L,original_quality_unchanged=TRUE,
 all_requested_cases_scored=TRUE,source_quality_sha256=original_hash,
 no_parallel_slope_test="Probit benchmark assumes parallel cumulative slopes; not an established substantive outcome model",
 no_between_scenario_exclusion_effect="Scenario-specific folds and samples; paired incremental contrasts only within scenario"),
 file.path(out,"auditoria_final.json"),pretty=TRUE,auto_unbox=TRUE,digits=NA)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
writeLines("COMPLETED",state);cat("ALL GROUPED BENCHMARKS COMPLETED\n")
