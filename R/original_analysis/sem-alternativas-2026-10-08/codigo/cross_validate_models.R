# Comparación predictiva interna de un conjunto FIJO de modelos SEM.
# Ningún modelo se elige ni modifica a partir de los pliegues de evaluación.
# La función usa PP1–PP3 del test únicamente como respuestas para calcular errores.
# Las predicciones condicionales usan TODOS los demás indicadores observados.

.sem_cv_numeric_matrix <- function(data, variables) {
  if (!all(variables %in% names(data)))
    stop("Missing observed variables: ", paste(setdiff(variables, names(data)), collapse=", "))
  bad <- variables[!vapply(data[variables], is.numeric, logical(1))]
  if(length(bad)) stop("Observed variables must be numeric: ",paste(bad,collapse=", "))
  x <- as.matrix(data[,variables,drop=FALSE])
  if(any(!is.finite(x))) stop("Non-finite values in observed data")
  storage.mode(x) <- "double"
  x
}

.sem_cv_spd <- function(x, tolerance=1e-8) {
  if(!is.matrix(x) || nrow(x)!=ncol(x) || any(!is.finite(x)))
    return(list(ok=FALSE, min_eigen=NA_real_, max_eigen=NA_real_, chol=NULL))
  if(max(abs(x-t(x))) > 1e-7*max(1,max(abs(x))))
    return(list(ok=FALSE, min_eigen=NA_real_, max_eigen=NA_real_, chol=NULL))
  values <- eigen((x+t(x))/2, symmetric=TRUE, only.values=TRUE)$values
  min_eig <- min(values); max_eig <- max(values)
  ok <- min_eig > tolerance*max(1,max_eig)
  upper <- if(ok) tryCatch(chol(x),error=function(e)NULL) else NULL
  list(ok=ok && !is.null(upper),min_eigen=min_eig,max_eigen=max_eig,chol=upper)
}

# Media condicional Gaussiana: E[Y|X] = muY + (X-muX) SigmaXX^-1 SigmaXY.
# Esta función no recibe ningún valor de Y del conjunto de prueba.
.sem_cv_conditional_pp <- function(x_test, mean, covariance, pp=c("PP1","PP2","PP3")) {
  ov <- names(mean)
  if(is.null(ov) || !all(pp %in% ov)) stop("Missing PP names in implied means")
  predictors <- setdiff(ov,pp)
  if(!identical(colnames(x_test),predictors)) stop("Predictor columns do not match implied variable order")
  sigma_xx <- covariance[predictors,predictors,drop=FALSE]
  checked <- .sem_cv_spd(sigma_xx)
  if(!checked$ok) stop("Conditional predictor covariance is not positive definite")
  # No ridge regularisation: an improper covariance invalidates the fold.
  coefficients <- backsolve(checked$chol,
                            forwardsolve(t(checked$chol),covariance[predictors,pp,drop=FALSE]))
  predicted <- sweep(sweep(x_test,2,mean[predictors],"-") %*% coefficients,
                     2,mean[pp],"+")
  colnames(predicted) <- pp
  predicted
}

# Puntuación logarítmica conjunta por observación para las mismas variables OV.
# Las covariables fixed.x contribuyen la misma referencia de entrenamiento a
# todos los modelos; se exige la matriz completa, incluidos sus momentos.
.sem_cv_gaussian_logscore <- function(y_test, mean, covariance) {
  if(!identical(colnames(y_test),names(mean))) stop("Log-score variable ordering mismatch")
  checked <- .sem_cv_spd(covariance)
  if(!checked$ok) stop("Joint implied covariance is not positive definite")
  centered <- sweep(y_test,2,mean,"-")
  whitened <- backsolve(checked$chol,t(centered),transpose=TRUE)
  quadratic <- colSums(whitened^2)
  log_det <- 2*sum(log(diag(checked$chol)))
  -0.5*(ncol(y_test)*log(2*pi)+log_det+quadratic)
}

.sem_cv_fit_check <- function(fit) {
  result <- list(converged=FALSE,admissible=FALSE,reason="",
                 min_residual_variance=NA_real_,min_theta_eigen=NA_real_,
                 min_psi_eigen=NA_real_,min_implied_cov_eigen=NA_real_,
                 max_abs_latent_correlation=NA_real_,ov=NULL,mean=NULL,covariance=NULL)
  if(inherits(fit,"error")) { result$reason <- conditionMessage(fit); return(result) }
  if(!inherits(fit,"lavaan")) { result$reason <- "fit_fun did not return a lavaan fit"; return(result) }
  result$converged <- isTRUE(lavaan::lavInspect(fit,"converged"))
  if(!result$converged) { result$reason <- "Model did not converge"; return(result) }
  tryCatch({
    pe <- lavaan::parameterEstimates(fit,se=FALSE)
    if(any(!is.finite(pe$est))) stop("Non-finite parameter estimates")
    variances <- pe$est[pe$op=="~~" & pe$lhs==pe$rhs]
    result$min_residual_variance <- min(variances)
    if(!length(variances) || any(variances<=0)) stop("Non-positive estimated or fixed variance")
    theta_check <- .sem_cv_spd(lavaan::lavInspect(fit,"theta"))
    result$min_theta_eigen <- theta_check$min_eigen
    if(!theta_check$ok) stop("Residual observed covariance is not positive definite")
    matrices <- lavaan::lavInspect(fit,"est")
    if(is.list(matrices) && !is.null(matrices$psi)) {
      psi_check <- .sem_cv_spd(matrices$psi)
      result$min_psi_eigen <- psi_check$min_eigen
      if(!psi_check$ok) stop("Residual latent covariance is not positive definite")
    }
    latent_cor <- lavaan::lavInspect(fit,"cor.lv")
    if(nrow(latent_cor)>1L) {
      result$max_abs_latent_correlation <- max(abs(latent_cor[row(latent_cor)!=col(latent_cor)]))
      if(!is.finite(result$max_abs_latent_correlation) || result$max_abs_latent_correlation>=1)
        stop("Improper latent correlation")
    } else result$max_abs_latent_correlation <- 0
    implied <- lavaan::lavInspect(fit,"implied")
    if(is.null(implied$cov) || is.null(implied$mean)) stop("Single-group implied means and covariance are required")
    result$ov <- lavaan::lavNames(fit,"ov")
    if(!setequal(colnames(implied$cov),result$ov) || !setequal(rownames(implied$cov),result$ov) ||
       !setequal(names(implied$mean),result$ov)) stop("Full named implied moments do not include all observed variables")
    result$covariance <- implied$cov[result$ov,result$ov,drop=FALSE]
    result$mean <- implied$mean[result$ov]
    if(any(!is.finite(result$mean))) stop("Non-finite implied means")
    covariance_check <- .sem_cv_spd(result$covariance)
    result$min_implied_cov_eigen <- covariance_check$min_eigen
    if(!covariance_check$ok) stop("Full implied covariance is not positive definite")
    result$admissible <- TRUE
    result
  },error=function(e) { result$admissible <- FALSE; result$reason <- conditionMessage(e); result })
}

.sem_cv_strata <- function(data) {
  if("Pais" %in% names(data)) {
    strata <- as.character(data$Pais)
  } else {
    sites <- c("C_Colombia","C_Peru","C_Mexico")
    if(!all(sites %in% names(data))) stop("Pais or all three site dummies are required")
    z <- .sem_cv_numeric_matrix(data,sites)
    if(any(!z %in% c(0,1)) || any(rowSums(z)>1)) stop("Invalid site dummy coding")
    strata <- ifelse(z[,1]==1,"Colombia",ifelse(z[,2]==1,"Peru",ifelse(z[,3]==1,"Mexico","Argentina")))
  }
  if(anyNA(strata) || any(!nzchar(strata))) stop("Missing site strata")
  strata
}

run_sem_cv <- function(data, model_specs, fit_fun, out_dir, folds=5, repeats=5, seed=20261008) {
  if(!requireNamespace("lavaan",quietly=TRUE)) stop("lavaan is required")
  if(!is.data.frame(data)) data <- as.data.frame(data)
  if(!is.list(model_specs) || !length(model_specs) || !is.function(fit_fun)) stop("Invalid fixed model set or fit function")
  ids <- names(model_specs)
  if(is.null(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop("model_specs must have unique nonempty model names")
  folds <- as.integer(folds); repeats <- as.integer(repeats); seed <- as.integer(seed)
  if(length(folds)!=1L || is.na(folds) || folds<2L || length(repeats)!=1L || is.na(repeats) || repeats<1L ||
     length(seed)!=1L || is.na(seed) || seed<0L || seed+repeats>=.Machine$integer.max) stop("Invalid folds, repeats or seed")
  pp <- c("PP1","PP2","PP3")
  .sem_cv_numeric_matrix(data,pp)
  strata <- .sem_cv_strata(data)
  if(any(table(strata)<folds)) stop("Every site must have at least one test observation in every fold")
  if(!dir.exists(out_dir) && !dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)) stop("Cannot create output directory")
  # Fijar previamente todos los pliegues evita que las rotaciones aleatorias
  # de un modelo alteren los conjuntos de evaluación de otro modelo.
  had_seed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if(had_seed) old_seed <- get(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  on.exit(if(had_seed) assign(".Random.seed",old_seed,envir=.GlobalEnv) else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
  assignment <- vector("list",repeats)
  for(r in seq_len(repeats)) {
    set.seed(seed+r)
    a <- integer(nrow(data))
    for(s in unique(strata)) {
      rows <- which(strata==s)
      a[rows] <- sample(rep(seq_len(folds),length.out=length(rows)))
    }
    assignment[[r]] <- a
  }
  diagnostics <- list(); summaries <- list(); shared_ov <- NULL
  d_index <- 0L; s_index <- 0L
  for(r in seq_len(repeats)) {
    for(id in ids) {
      total_n <- 0L; sum_squared_error <- 0; sum_absolute_error <- 0; sum_logscore <- 0
      valid_folds <- 0L
      for(k in seq_len(folds)) {
        test_rows <- which(assignment[[r]]==k)
        train_rows <- which(assignment[[r]]!=k)
        warning_messages <- character()
        fit <- tryCatch(withCallingHandlers(
          fit_fun(model_specs[[id]],data[train_rows,,drop=FALSE],estimator="ML",fast=TRUE),
          warning=function(w){warning_messages<<-c(warning_messages,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
        checked <- .sem_cv_fit_check(fit)
        record <- data.frame(repetition=r,model=id,fold=k,n_train=length(train_rows),n_test=length(test_rows),
                             converged=checked$converged,admissible=checked$admissible,
                             min_residual_variance=checked$min_residual_variance,
                             min_theta_eigen=checked$min_theta_eigen,min_psi_eigen=checked$min_psi_eigen,
                             min_implied_cov_eigen=checked$min_implied_cov_eigen,
                             max_abs_latent_correlation=checked$max_abs_latent_correlation,
                             RMSE_PP_mean=NA_real_,MAE_PP_mean=NA_real_,mean_joint_logscore=NA_real_,
                             reason=checked$reason,warnings=paste(unique(warning_messages),collapse=" | "),stringsAsFactors=FALSE)
        if(checked$admissible) {
          score <- tryCatch({
            if(is.null(shared_ov)) shared_ov <- checked$ov
            if(!setequal(checked$ov,shared_ov)) stop("Models do not use the same observed variable set")
            # Un único orden compartido para todas las puntuaciones y modelos.
            mu <- checked$mean[shared_ov]
            sigma <- checked$covariance[shared_ov,shared_ov,drop=FALSE]
            observed <- .sem_cv_numeric_matrix(data[test_rows,,drop=FALSE],shared_ov)
            predictor_names <- setdiff(shared_ov,pp)
            predicted_pp <- .sem_cv_conditional_pp(observed[,predictor_names,drop=FALSE],mu,sigma,pp)
            target <- rowMeans(observed[,pp,drop=FALSE])
            prediction <- rowMeans(predicted_pp)
            error <- prediction-target
            logscore <- .sem_cv_gaussian_logscore(observed,mu,sigma)
            if(any(!is.finite(c(error,logscore)))) stop("Non-finite held-out score")
            list(n=length(error),squared=sum(error^2),absolute=sum(abs(error)),logscore=sum(logscore))
          },error=function(e)e)
          if(inherits(score,"error")) {
            record$admissible <- FALSE; record$reason <- conditionMessage(score)
          } else {
            record$RMSE_PP_mean <- sqrt(score$squared/score$n)
            record$MAE_PP_mean <- score$absolute/score$n
            record$mean_joint_logscore <- score$logscore/score$n
            total_n <- total_n+score$n; sum_squared_error <- sum_squared_error+score$squared
            sum_absolute_error <- sum_absolute_error+score$absolute; sum_logscore <- sum_logscore+score$logscore
            valid_folds <- valid_folds+1L
          }
        }
        d_index <- d_index+1L; diagnostics[[d_index]] <- record
      }
      # Si un solo pliegue falla, no se resume únicamente el subconjunto fácil.
      complete <- valid_folds==folds && total_n==nrow(data)
      s_index <- s_index+1L
      summaries[[s_index]] <- data.frame(repetition=r,model=id,folds_requested=folds,
        folds_valid=valid_folds,n_test_scored=total_n,complete_repetition=complete,
        RMSE_PP_mean=if(complete)sqrt(sum_squared_error/total_n)else NA_real_,
        MAE_PP_mean=if(complete)sum_absolute_error/total_n else NA_real_,
        mean_joint_logscore=if(complete)sum_logscore/total_n else NA_real_,stringsAsFactors=FALSE)
      cat("CV repetition ",r," model ",id,": ",valid_folds,"/",folds," valid folds\n",sep="")
      utils::write.csv(do.call(rbind,diagnostics),file.path(out_dir,"cv_fold_diagnostics.csv"),row.names=FALSE,na="")
      utils::write.csv(do.call(rbind,summaries),file.path(out_dir,"cv_repeat_summary.csv"),row.names=FALSE,na="")
    }
  }
  diagnostic_table <- do.call(rbind,diagnostics)
  summary_table <- do.call(rbind,summaries)
  common_complete <- vapply(seq_len(repeats),function(r)all(summary_table$complete_repetition[summary_table$repetition==r]),logical(1))
  overall <- do.call(rbind,lapply(ids,function(id) {
    rows <- summary_table[summary_table$model==id,,drop=FALSE]
    all_complete <- all(rows$complete_repetition)
    # Estos SD describen variación de particiones reutilizadas; no son SE ni CI.
    data.frame(model=id,repetitions_requested=repeats,repetitions_valid=sum(rows$complete_repetition),
      all_repetitions_complete=all_complete,
      RMSE_PP_mean=if(all_complete)mean(rows$RMSE_PP_mean)else NA_real_,
      RMSE_partition_SD=if(all_complete && repeats>1L)sd(rows$RMSE_PP_mean)else NA_real_,
      MAE_PP_mean=if(all_complete)mean(rows$MAE_PP_mean)else NA_real_,
      MAE_partition_SD=if(all_complete && repeats>1L)sd(rows$MAE_PP_mean)else NA_real_,
      mean_joint_logscore=if(all_complete)mean(rows$mean_joint_logscore)else NA_real_,
      joint_logscore_partition_SD=if(all_complete && repeats>1L)sd(rows$mean_joint_logscore)else NA_real_,
      common_complete_repetitions=sum(common_complete),stringsAsFactors=FALSE)
  }))
  metadata <- list(N=nrow(data),models=ids,folds=folds,repeats=repeats,seed=seed,
    strata_counts=as.list(table(strata)),observed_variables=shared_ov,
    target="Observed row mean of PP1 PP2 PP3",predictors="All observed variables other than PP1 PP2 PP3",
    estimator="ML training fit without SE or tests",score="Full joint Gaussian mean log density; larger is better",
    RMSE_MAE="Unclipped conditional predicted PP mean against observed PP mean; smaller is better",
    scaling="Input covariates retained as supplied; no transformation learned from test responses",
    incomplete_policy="Any failed fold invalidates its repetition; any failed repetition prevents the overall metric for that model",
    uncertainty="Partition SD is descriptive only; no independent-fold tests, confidence intervals or SE are computed",
    interpretation="Repeated internal validation of questionnaire responses; not external validation or observed behavioural prediction",
    privacy="Only aggregate summaries exported; no participant-level predictions, scores or fold assignments exported")
  utils::write.csv(overall,file.path(out_dir,"cv_overall_summary.csv"),row.names=FALSE,na="")
  if(requireNamespace("jsonlite",quietly=TRUE)) jsonlite::write_json(metadata,file.path(out_dir,"cv_metadata.json"),pretty=TRUE,auto_unbox=TRUE)
  saveRDS(metadata,file.path(out_dir,"cv_metadata.rds"))
  invisible(list(overall=overall,repetitions=summary_table,fold_diagnostics=diagnostic_table,metadata=metadata))
}
