# REVISION NEUTRAL ORIGINAL -- PREDICCION, OBSERVACIONES 4 Y 5
# No cambia manuscritos, bases ni versiones anteriores.
# Ejecutar: Rscript 01_regularizacion_e_incertidumbre_agrupada.R
# PAPER_PROYECTO permite cambiar la raiz; PRED_BOOT_REPS solo sirve para
# cronometrajes de desarrollo. La ejecucion definitiva preespecificada usa 200.
#
# MAPA DE TABLAS PROVISIONALES (numeracion final a cargo del manuscrito):
# [CV_REGULARIZADO] CV_regularizado_comparacion.csv: B2/B3, RMSE/MAE/RPS,
#    media de cinco repeticiones exteriores agrupadas de cinco folds.
# [OOB_PAREADO] OOB_regularizado_resumen.csv: punto OOB e intervalo de
#    remuestreo 2.5-97.5, SEPARADOS del estimando CV; ningun IC poblacional.
# [CALIBRACION_CV] CV_regularizado_calibracion_cortes.csv y _bandas.csv:
#    cuatro cortes de la distribucion ordinal y bandas de probabilidad fijas.
# [CALIBRACION_OOB] OOB_regularizado_calibracion_cortes.csv: calibracion
#    marginal descriptiva en las 200 replicas OOB con reentrenamiento.
# [LISO_PREVIO] LISO_previo_agregados/: solo agregados del analisis ya auditado.
# RDS en privados contienen predicciones/asignaciones individuales: NO compartir.

options(stringsAsFactors=FALSE, warn=1)
project <- Sys.getenv("PAPER_PROYECTO", "/Users/luisbenites/Documents/Codex/Paper con MA - 1/reproduccion-integral-V3-2026-10-10/proyecto")
out <- file.path(project, "revision-observaciones-2026-10-09", "prediccion")
.libPaths(c(file.path(project,"revision-revista-2026-10-04/auditoria-codigos/R-library"), .libPaths()))
stopifnot(requireNamespace("nnet",quietly=TRUE), requireNamespace("jsonlite",quietly=TRUE))
dir.create(file.path(out,"privados","OOB_predicciones"),recursive=TRUE,showWarnings=FALSE)
write_csv <- function(x,name) write.csv(x,file.path(out,name),row.names=FALSE,fileEncoding="UTF-8",na="")
write_json <- function(x,name) jsonlite::write_json(x,file.path(out,name),auto_unbox=TRUE,pretty=TRUE,null="null",na="null",digits=16)
sha <- function(path) strsplit(system2("/usr/bin/shasum",c("-a","256",shQuote(path)),stdout=TRUE)," +")[[1]][1]
clock_start <- Sys.time()
protocol_path <- file.path(out,"protocolo_antes_ajustes.json")
stopifnot(file.exists(protocol_path))
protocol_sha <- sha(protocol_path)
src_q <- file.path(project,"revision-3-revista-colombiana-psicologia-2026-10-08/analisis_sensibilidad/privados/datos_y_banderas_para_sensibilidades.rds")
src_raw <- file.path(project,"revision-integral-2026-10-04/analysis/outputs/canonical-full-20261004-02/01_audit/analytic_data.rds")
src_assign <- file.path(project,"revision-3-revista-colombiana-psicologia-2026-10-08/analisis_sensibilidad/grouped_cv/all/privados/asignaciones_por_grupo.rds")
hash_before <- c(q=sha(src_q),raw=sha(src_raw),assignments=sha(src_assign))
stopifnot(hash_before["q"]=="283d622d631faabfce1da2e5c947188d6b006d40df0f6b2690945d24f99be76e",
          hash_before["raw"]=="64220098e861a191aea484e1a470f2d28133e08d7472eb3fc51fd91ffbe8c501")
q <- readRDS(src_q); raw <- readRDS(src_raw); a <- readRDS(src_assign)
dat <- q$data
dat$Edad_original <- raw$Edad
dat$NSE_original <- raw$NSE
group <- q$group
site <- as.character(dat$Pais)
aps <- c("APS1","APS2","APS3","APS4","APS5","APS6","APS7","APS9","APS10","APS14",
         "APS11","APS12","APS13","APS15","APS16","APS18","APS19","APS20","APS21")
controls <- c("CSR1_bin","C_Colombia","C_Peru","C_Mexico","Sexo_num","Estudios_num","Edad_original","NSE_original")
blocks <- list(B3=c(controls,aps),B2=c(controls,aps,"Benevolencia","Universalismo","Poder","Logro"))
outcomes <- c("PP1","PP2","PP3")
lambda_grid <- c(.001,.01,.1)
nboot <- as.integer(Sys.getenv("PRED_BOOT_REPS","200"))
stopifnot(nrow(dat)==653L,length(unique(group))==650L,
          identical(as.integer(a$group),as.integer(group)),identical(as.integer(a$original_row),as.integer(q$original_row)),
          identical(as.character(raw$Pais),site),all(vapply(split(site,group),function(z)length(unique(z))==1L,logical(1))),
          !anyNA(dat[,unique(c(unlist(blocks),outcomes))]),all(dat[,outcomes]>=1&dat[,outcomes]<=5),
          all(vapply(a$assignment,function(z)all(vapply(split(z,group),function(v)length(unique(v))==1L,logical(1))),logical(1))))
shared_vars <- intersect(c(aps,outcomes,"Benevolencia","Universalismo","Poder","Logro","CSR1_bin","Sexo_num","Estudios_num"),names(raw))
stopifnot(all(vapply(shared_vars,function(v)isTRUE(all.equal(as.numeric(dat[[v]]),as.numeric(raw[[v]]))),logical(1))))

# La red NO tiene capa oculta: es regresion multinomial softmax lineal.
# K se fija a 5 usando una matriz one-hot aunque alguna columna este vacia.
# nnet penaliza TODOS los pesos, incluidos los cinco sesgos/interceptos.
# Todos los logits se estiman simetricamente; el ridge identifica su traslacion.
# Esto evita que una categoria ausente fuerce interceptos infinitos. No impone
# orden proporcional y por ello es comparador distinto del SEM/probit ordinal.
diag_env <- new.env(parent=emptyenv()); diag_env$rows <- list(); diag_env$n <- 0L
fit_ridge <- function(x,y,w,lambda,tag) {
  Y <- matrix(0,nrow(x),5,dimnames=list(NULL,as.character(1:5)))
  Y[cbind(seq_along(y),as.integer(y))] <- 1
  warning_text <- character()
  fit <- tryCatch(withCallingHandlers(nnet::nnet.default(x,Y,weights=w,size=0,skip=TRUE,softmax=TRUE,
                 rang=0,decay=lambda*sum(w),maxit=300,trace=FALSE,MaxNWts=10000,
                 abstol=1e-8,reltol=1e-8),warning=function(e){warning_text <<- c(warning_text,conditionMessage(e));invokeRestart("muffleWarning")}),error=function(e)e)
  retry <- FALSE
  if(!inherits(fit,"error") && fit$convergence!=0L) {
    retry <- TRUE
    fit <- tryCatch(withCallingHandlers(nnet::nnet.default(x,Y,weights=w,size=0,skip=TRUE,softmax=TRUE,
                     Wts=fit$wts,decay=lambda*sum(w),maxit=1000,trace=FALSE,MaxNWts=10000,
                     abstol=1e-8,reltol=1e-8),warning=function(e){warning_text <<- c(warning_text,conditionMessage(e));invokeRestart("muffleWarning")}),error=function(e)e)
  }
  valid <- !inherits(fit,"error") && fit$convergence==0L && all(is.finite(fit$wts)) && is.finite(fit$value)
  diag_env$n <- diag_env$n+1L
  diag_env$rows[[diag_env$n]] <- data.frame(Tag=tag,Lambda=lambda,Train_rows=nrow(x),Train_weight=sum(w),
      Min_category_weight=min(colSums(Y*w)),Valid=valid,Retry=retry,
      Convergence=if(inherits(fit,"error"))NA_integer_ else fit$convergence,
      Objective=if(inherits(fit,"error"))NA_real_ else fit$value,
      Warnings=paste(unique(warning_text),collapse=" | "),Error=if(inherits(fit,"error"))conditionMessage(fit) else "")
  if(!valid) return(NULL)
  fit
}
scale_fit <- function(x,w) {
  mu <- colSums(x*w)/sum(w)
  xc <- sweep(x,2,mu,"-")
  sdv <- sqrt(colSums(xc^2*w)/max(1,sum(w)-1))
  sdv[!is.finite(sdv)|sdv<1e-12] <- 1
  list(mean=mu,sd=sdv)
}
scale_apply <- function(x,s) sweep(sweep(x,2,s$mean,"-"),2,s$sd,"/")
pred_check <- function(fit,x) {
  p <- as.matrix(predict(fit,x,type="raw"))
  if(ncol(p)!=5L || any(!is.finite(p)) || any(p<=0) || max(abs(rowSums(p)-1))>1e-10) stop("Probabilidades invalidas")
  colnames(p) <- as.character(1:5)
  p
}
rps_loss <- function(y,p) {
  Fp <- t(apply(p,1,cumsum))[,1:4,drop=FALSE]
  indicator <- outer(as.integer(y),1:4,"<=")*1
  rowMeans((Fp-indicator)^2)
}
scores <- function(y,p) {
  m <- as.vector(p%*%(1:5)); cumulative <- t(apply(p,1,cumsum))[,1:4,drop=FALSE]
  obs <- outer(as.integer(y),1:4,"<=")*1
  c(RMSE=sqrt(mean((m-y)^2)),MAE=mean(abs(m-y)),RPS=mean(rps_loss(y,p)),
    Mean_error=mean(m-y),Cum_gap_max=max(abs(colMeans(cumulative-obs))),
    Min_probability=min(p),Zero_observed_probability=sum(p[cbind(seq_along(y),as.integer(y))]<=0))
}
# Cada grupo original conserva su fold interior, incluso si fue muestreado
# varias veces. Ningun escalado/tuning consulta la particion exterior de test.
inner_folds <- function(indices,seed) {
  set.seed(seed)
  first <- indices[!duplicated(group[indices])]
  assignment <- integer(length(first))
  for(s in unique(site)) {
    zz <- which(site[first]==s)
    assignment[sample(zz)] <- rep(1:2,length.out=length(zz))
  }
  assignment[match(group[indices],group[first])]
}
train_tune_predict <- function(train,test,weights,vars,yvar,folds,tag) {
  X <- as.matrix(dat[,vars,drop=FALSE]); storage.mode(X) <- "double"
  y <- as.integer(dat[[yvar]])
  interior <- matrix(Inf,length(lambda_grid),2)
  interior_n <- numeric(2)
  for(k in 1:2) {
    ii <- which(folds!=k); vv <- which(folds==k)
    sc <- scale_fit(X[train[ii],,drop=FALSE],weights[ii])
    xtrain <- scale_apply(X[train[ii],,drop=FALSE],sc)
    xval <- scale_apply(X[train[vv],,drop=FALSE],sc)
    interior_n[k] <- sum(weights[vv])
    for(l in seq_along(lambda_grid)) {
      fit <- fit_ridge(xtrain,y[train[ii]],weights[ii],lambda_grid[l],paste(tag,"inner",k,l,sep="/"))
      if(!is.null(fit)) {
        pp <- tryCatch(pred_check(fit,xval),error=function(e)NULL)
        if(!is.null(pp)) interior[l,k] <- weighted.mean(rps_loss(y[train[vv]],pp),weights[vv])
      }
    }
  }
  tune_score <- as.vector(interior%*%interior_n/sum(interior_n))
  if(!any(is.finite(tune_score))) return(list(valid=FALSE,error="Todos los candidatos de tuning fallaron"))
  best <- order(tune_score,-lambda_grid)[1L]
  sc <- scale_fit(X[train,,drop=FALSE],weights)
  fit <- fit_ridge(scale_apply(X[train,,drop=FALSE],sc),y[train],weights,lambda_grid[best],paste0(tag,"/final"))
  if(is.null(fit)) return(list(valid=FALSE,error="Ajuste final no convergente"))
  pp <- tryCatch(pred_check(fit,scale_apply(X[test,,drop=FALSE],sc)),error=function(e)NULL)
  if(is.null(pp)) return(list(valid=FALSE,error="Probabilidades finales invalidas"))
  list(valid=TRUE,p=pp,lambda=lambda_grid[best],tune_score=tune_score,scaler=sc)
}
cal_cuts <- function(y,p,item,block,replication,design) {
  Fp <- t(apply(p,1,cumsum))[,1:4,drop=FALSE]
  data.frame(Design=design,Replication=replication,Outcome=item,Block=block,Cut=1:4,N=length(y),
      Predicted=colMeans(Fp),Observed=vapply(1:4,function(k)mean(y<=k),numeric(1)),
      Gap=colMeans(Fp)-vapply(1:4,function(k)mean(y<=k),numeric(1)))
}
cal_bands <- function(y,p,item,block,replication,design) {
  Fp <- t(apply(p,1,cumsum))[,1:4,drop=FALSE]
  ans <- list(); at <- 0L
  for(k in 1:4) {
    bin <- pmin(5L,floor(Fp[,k]*5)+1L)
    for(b in 1:5) {
      z <- which(bin==b); at <- at+1L
      ans[[at]] <- data.frame(Design=design,Replication=replication,Outcome=item,Block=block,Cut=k,
         Band=paste0(if(b==1)"[" else "[",sprintf("%.1f",(b-1)/5),",",sprintf("%.1f",b/5),if(b==5)"]" else ")"),
         N=length(z),Predicted=if(length(z))mean(Fp[z,k]) else NA_real_,
         Observed=if(length(z))mean(y[z]<=k) else NA_real_)
    }
  }
  do.call(rbind,ans)
}
to_metric_row <- function(obj,y,item,block,replication,design,ntrain,nweight,ntest,ngtest) {
  z <- if(obj$valid)scores(y,obj$p) else rep(NA_real_,7)
  names(z) <- c("RMSE","MAE","RPS","Mean_error","Cum_gap_max","Min_probability","Zero_observed_probability")
  data.frame(Design=design,Replication=replication,Outcome=item,Block=block,
      N_train_distinct_rows=ntrain,N_train_weight=nweight,N_test=ntest,N_test_groups=ngtest,
      Valid=obj$valid,Lambda=if(obj$valid)obj$lambda else NA_real_,Error=if(obj$valid)"" else obj$error,t(z),check.names=FALSE)
}
save_progress <- function(stage,done,target,metric_rows) {
  write_json(list(stage=stage,completed=done,target=target,elapsed_seconds=as.numeric(difftime(Sys.time(),clock_start,units="secs")),
      fits=diag_env$n,valid_fits=sum(vapply(diag_env$rows,function(z)isTRUE(z$Valid),logical(1))),
      retries=sum(vapply(diag_env$rows,function(z)isTRUE(z$Retry),logical(1)))),"estado_ejecucion.json")
  if(length(metric_rows)) write_csv(do.call(rbind,metric_rows),paste0(stage,"_metricas_checkpoint.csv"))
  write_csv(do.call(rbind,diag_env$rows),"diagnosticos_ajustes.csv")
  cat(stage,done,"/",target,"elapsed",round(as.numeric(difftime(Sys.time(),clock_start,units="secs")),1),"seconds; fits",diag_env$n,"\n");flush.console()
}

# [CV_REGULARIZADO] Reutilizar exactamente las 25 particiones agrupadas V3.
# Se agregan predicciones de cinco folds antes de calcular la metrica de cada
# repeticion. Las 5 repeticiones NO producen un SE artificial entre folds.
cv_rows <- list(); cv_cuts <- list(); cv_bands <- list(); cv_preds <- list()
for(r in seq_along(a$assignment)) {
  pred_rep <- lapply(outcomes,function(y)lapply(blocks,function(v)matrix(NA_real_,nrow(dat),5)))
  names(pred_rep) <- outcomes
  lambda_rep <- lapply(outcomes,function(y)lapply(blocks,function(v)numeric(5)));names(lambda_rep)<-outcomes
  valid_rep <- matrix(TRUE,3,2,dimnames=list(outcomes,names(blocks)))
  for(f in 1:5) {
    train <- which(a$assignment[[r]]!=f); test <- which(a$assignment[[r]]==f)
    folds <- inner_folds(train,2026100904L+1000L*r+f)
    for(y in outcomes) for(b in names(blocks)) {
      obj <- train_tune_predict(train,test,rep(1,length(train)),blocks[[b]],y,folds,paste("CV",r,f,y,b,sep="/"))
      if(obj$valid) {pred_rep[[y]][[b]][test,] <- obj$p;lambda_rep[[y]][[b]][f]<-obj$lambda} else valid_rep[y,b]<-FALSE
    }
    if(r==1L&&f==1L) save_progress("CV",1,25,cv_rows)
  }
  for(y in outcomes) for(b in names(blocks)) {
    valid <- valid_rep[y,b]&&all(is.finite(pred_rep[[y]][[b]]))
    obj <- list(valid=valid,p=pred_rep[[y]][[b]],lambda=if(valid)mean(lambda_rep[[y]][[b]]) else NA_real_,error="Alguna particion exterior fallo")
    cv_rows[[length(cv_rows)+1L]] <- to_metric_row(obj,dat[[y]],y,b,r,"CV5x5_grouped",NA_real_,NA_real_,nrow(dat),length(unique(group)))
    if(valid) {
      cv_cuts[[length(cv_cuts)+1L]] <- cal_cuts(dat[[y]],obj$p,y,b,r,"CV5x5_grouped")
      cv_bands[[length(cv_bands)+1L]] <- cal_bands(dat[[y]],obj$p,y,b,r,"CV5x5_grouped")
    }
  }
  cv_preds[[r]] <- list(original_row=q$original_row,y=dat[,outcomes],p=pred_rep,lambda=lambda_rep,valid=valid_rep)
  saveRDS(cv_preds,file.path(out,"privados","CV_predicciones_agrupadas.rds"))
  save_progress("CV",r*5,25,cv_rows)
}
cv <- do.call(rbind,cv_rows)
write_csv(cv,"CV_regularizado_metricas_por_repeticion.csv")
write_csv(do.call(rbind,cv_cuts),"CV_regularizado_calibracion_cortes_por_repeticion.csv")
write_csv(do.call(rbind,cv_bands),"CV_regularizado_calibracion_bandas_por_repeticion.csv")
cv_comparison <- lapply(outcomes,function(y) {
  z3 <- cv[cv$Outcome==y&cv$Block=="B3",];z2 <- cv[cv$Outcome==y&cv$Block=="B2",]
  ans <- data.frame(Outcome=y,Design="CV5x5_grouped",Complete_repetitions=sum(z3$Valid&z2$Valid),N=653)
  for(m in c("RMSE","MAE","RPS")) {
    ok <- z3$Valid&z2$Valid
    ans[[paste0("B3_",m)]] <- mean(z3[[m]][ok]);ans[[paste0("B2_",m)]] <- mean(z2[[m]][ok])
    ans[[paste0("Delta_",m)]] <- mean(z2[[m]][ok]-z3[[m]][ok])
    ans[[paste0("Delta_",m,"_percent")]] <- mean(100*(z2[[m]][ok]-z3[[m]][ok])/z3[[m]][ok])
  }
  ans
})
write_csv(do.call(rbind,cv_comparison),"CV_regularizado_comparacion.csv")
cv_cal <- do.call(rbind,cv_cuts)
write_csv(aggregate(cbind(Predicted,Observed,Gap)~Outcome+Block+Cut,data=cv_cal,FUN=mean),"CV_regularizado_calibracion_cortes.csv")
band <- do.call(rbind,cv_bands)
band$Predicted_sum <- ifelse(band$N>0,band$N*band$Predicted,0)
band$Observed_sum <- ifelse(band$N>0,band$N*band$Observed,0)
bandagg <- aggregate(cbind(N,Predicted_sum,Observed_sum)~Outcome+Block+Cut+Band,data=band,FUN=sum)
bandagg$Predicted <- with(bandagg,ifelse(N>0,Predicted_sum/N,NA_real_));bandagg$Observed<-with(bandagg,ifelse(N>0,Observed_sum/N,NA_real_))
bandagg$Gap<-bandagg$Predicted-bandagg$Observed
write_csv(bandagg[,c("Outcome","Block","Cut","Band","N","Predicted","Observed","Gap")],"CV_regularizado_calibracion_bandas.csv")

# [OOB_PAREADO] NO es el intervalo del punto CV anterior. En cada replica se
# vuelve a muestrear, escalar, tunear y entrenar. Los grupos se muestrean dentro
# de las cuatro instituciones FIJAS. Los test OOB solo evalua grupos no usados.
# Las multiplicidades se aplican como pesos; nunca repartir clones entre folds.
boot_rows <- list(); boot_cuts <- list(); partitions <- list()
group_first <- which(!duplicated(group)); gids <- group[group_first]
gsite <- site[group_first]
for(r in seq_len(nboot)) {
  set.seed(2026100904L+r)
  draws <- unlist(lapply(unique(site),function(s)sample(gids[gsite==s],sum(gsite==s),replace=TRUE)),use.names=FALSE)
  counts <- table(factor(draws,levels=gids))
  row_weight <- as.numeric(counts[match(group,gids)])
  train <- which(row_weight>0);test <- which(row_weight==0)
  stopifnot(length(intersect(group[train],group[test]))==0L,length(test)>0L)
  weights <- row_weight[train]
  folds <- inner_folds(train,2026100904L+100000L+r)
  partitions[[r]] <- data.frame(Replication=r,N_train_distinct_rows=length(train),N_train_weight=sum(weights),N_train_groups=length(unique(group[train])),
      N_test=length(test),N_test_groups=length(unique(group[test])),Crossing_groups=length(intersect(unique(group[train]),unique(group[test]))))
  preds <- list()
  for(y in outcomes) for(b in names(blocks)) {
    obj <- train_tune_predict(train,test,weights,blocks[[b]],y,folds,paste("OOB",r,y,b,sep="/"))
    boot_rows[[length(boot_rows)+1L]] <- to_metric_row(obj,dat[[y]][test],y,b,r,"OOB_group_bootstrap_conditional_4sites",length(train),sum(weights),length(test),length(unique(group[test])))
    preds[[paste(y,b,sep="_")]] <- list(valid=obj$valid,p=if(obj$valid)obj$p else NULL,lambda=if(obj$valid)obj$lambda else NULL,
        scaler=if(obj$valid)obj$scaler else NULL,error=if(obj$valid)"" else obj$error)
    if(obj$valid)boot_cuts[[length(boot_cuts)+1L]] <- cal_cuts(dat[[y]][test],obj$p,y,b,r,"OOB_group_bootstrap_conditional_4sites")
  }
  saveRDS(list(original_row=q$original_row[test],y=dat[test,outcomes],train=train,test=test,group=group,
      weights=weights,inner_folds=folds,predictions=preds),file.path(out,"privados","OOB_predicciones",sprintf("replica_%03d.rds",r)))
  if(r%%20L==0L||r==1L||r==nboot) {
    saveRDS(list(metrics=boot_rows,cuts=boot_cuts,partitions=partitions),file.path(out,"privados","OOB_checkpoint.rds"))
    save_progress("OOB",r,nboot,boot_rows)
  }
}
boot <- do.call(rbind,boot_rows)
write_csv(boot,"OOB_regularizado_metricas_por_replica.csv")
write_csv(do.call(rbind,partitions),"OOB_particiones_agregadas.csv")
write_csv(do.call(rbind,boot_cuts),"OOB_regularizado_calibracion_cortes_por_replica.csv")
summary_rows <- list(); paired_rows <- list()
for(y in outcomes) {
  z3 <- boot[boot$Outcome==y&boot$Block=="B3",];z2 <- boot[boot$Outcome==y&boot$Block=="B2",]
  stopifnot(identical(z3$Replication,z2$Replication))
  for(m in c("RMSE","MAE","RPS")) {
    ok <- z3$Valid&z2$Valid
    delta <- z2[[m]][ok]-z3[[m]][ok]
    relative <- 100*delta/z3[[m]][ok]
    qq <- function(x)if(length(x))as.numeric(quantile(x,c(.025,.975),type=7)) else c(NA_real_,NA_real_)
    dq <- qq(delta);rq <- qq(relative)
    summary_rows[[length(summary_rows)+1L]] <- data.frame(Outcome=y,Metric=m,Design="OOB paired group resampling conditional on 4 institutions",
       Point_definition="Mean of OOB paired replicate metrics; separate from CV5x5",
       Interval_label="OOB paired resampling 2.5-97.5 percentile range; coverage not validated",
       Replicates_requested=nboot,Replicates_complete=sum(ok),Replicates_failed=sum(!ok),
       B3_mean=mean(z3[[m]][ok]),B2_mean=mean(z2[[m]][ok]),Delta_mean=mean(delta),Delta_median=median(delta),
       Delta_lower=dq[1],Delta_upper=dq[2],Relative_percent_mean=mean(relative),Relative_percent_lower=rq[1],Relative_percent_upper=rq[2],
       Better_count=sum(delta<0),Equal_count=sum(delta==0),Worse_count=sum(delta>0),
       Min_N_OOB=min(z3$N_test[ok]),Max_N_OOB=max(z3$N_test[ok]),Mean_N_OOB=mean(z3$N_test[ok]))
    paired_rows[[length(paired_rows)+1L]] <- data.frame(Outcome=y,Metric=m,Replication=z3$Replication,
         Valid_pair=ok,B3=z3[[m]],B2=z2[[m]],Delta=z2[[m]]-z3[[m]],Relative_percent=100*(z2[[m]]-z3[[m]])/z3[[m]])
  }
}
write_csv(do.call(rbind,summary_rows),"OOB_regularizado_resumen.csv")
write_csv(do.call(rbind,paired_rows),"OOB_diferencias_emparejadas_por_replica.csv")
bc <- do.call(rbind,boot_cuts)
cal_summary <- list()
for(y in outcomes)for(b in names(blocks))for(k in 1:4) {
  z <- bc[bc$Outcome==y&bc$Block==b&bc$Cut==k,]
  qq <- as.numeric(quantile(z$Gap,c(.025,.975),type=7))
  cal_summary[[length(cal_summary)+1L]] <- data.frame(Outcome=y,Block=b,Cut=k,Complete_replicates=nrow(z),
      Predicted_mean=mean(z$Predicted),Observed_mean=mean(z$Observed),Gap_mean=mean(z$Gap),Gap_lower=qq[1],Gap_upper=qq[2],
      Interval_label="OOB resampling 2.5-97.5 percentile range; coverage not validated")
}
write_csv(do.call(rbind,cal_summary),"OOB_regularizado_calibracion_cortes.csv")

# [LISO_PREVIO] Reutilizacion documentada, sin recalcular ni cambiar resultados.
# Los coeficientes/dummies son distintos del comparador ridge y no se mezclan.
liso_dir <- file.path(project,"revision-4-ARLA-2026-10-09/transportabilidad_institucional")
liso_dest <- file.path(out,"LISO_previo_agregados")
dir.create(liso_dest,showWarnings=FALSE)
liso_files <- c("metricas_y_diagnosticos_por_institucion.csv","resumen_errores_agrupados.csv","resumen_incremento_predictivo.csv",
                "contrastes_pareados_por_institucion.csv","particiones_y_duplicados.csv","frecuencias_train_test.csv")
liso_provenance <- lapply(liso_files,function(name) {
  src <- file.path(liso_dir,name);dst<-file.path(liso_dest,name)
  stopifnot(file.exists(src));file.copy(src,dst,overwrite=TRUE);stopifnot(sha(src)==sha(dst))
  list(file=name,source=src,sha256=sha(dst))
})
hash_after <- c(q=sha(src_q),raw=sha(src_raw),assignments=sha(src_assign))
stopifnot(identical(hash_before,hash_after),sha(protocol_path)==protocol_sha)
diagnostics <- do.call(rbind,diag_env$rows)
write_csv(diagnostics,"diagnosticos_ajustes.csv")
write_json(list(status="Completed; pending independent reconstruction audit",protocol_sha256=protocol_sha,source_sha256_before=as.list(hash_before),source_sha256_after=as.list(hash_after),
   n=653,groups=650,sites=as.list(table(site)),lambda_grid=lambda_grid,inner_folds=2,outer_folds=5,outer_repetitions=5,bootstrap_repetitions=nboot,
   model="Symmetric five-category linear softmax with ridge on all slopes and all intercepts; different from SEM and cumulative probit",
   decay_definition="lambda times weighted training N",fits=nrow(diagnostics),failed_fits=sum(!diagnostics$Valid),retried_fits=sum(diagnostics$Retry),
   warning_fits=sum(nzchar(diagnostics$Warnings)),missing_training_categories_fits=sum(diagnostics$Min_category_weight==0),
   final_fit_count=sum(grepl("/final$",diagnostics$Tag)),failed_final_fits=sum(!diagnostics$Valid&grepl("/final$",diagnostics$Tag)),
   max_probability_sum_error="Checked <= 1e-10 in every prediction",zero_probability_policy="All five predicted probabilities must be finite and strictly positive",
   cv_estimand="Mean of five grouped CV repeat metrics, each evaluating all 653 original records once",
   oob_estimand="Mean of paired OOB replicate metrics after bootstrap group-weighted fitting and full internal tuning, conditional on four fixed institutions",
   interval="OOB paired resampling 2.5-97.5 percentile range, descriptive and exploratory; no validated nominal coverage",
   limitations=c("Not the interval for the CV point estimate","Groups exchangeable within institutions assumed, not proof of all residual independence",
     "Four fixed institutions, no institution-population inference","Country/institution confounded","Rare category support is model-based, not validated by few category-1 records",
     "Grouped identical records do not establish unique participants"),LISO_reuse=liso_provenance,
   elapsed_seconds=as.numeric(difftime(Sys.time(),clock_start,units="secs")),sessionInfo=capture.output(sessionInfo())),"reporte_regularizacion_incertidumbre.json")
save_progress("Completed",nboot,nboot,boot_rows)
cat("Completado. Verificar de forma independiente las metricas desde RDS privados antes de integrar cifras.\n")
