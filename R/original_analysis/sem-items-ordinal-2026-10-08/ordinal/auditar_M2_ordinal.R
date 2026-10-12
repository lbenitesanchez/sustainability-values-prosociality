# Auditoría de objetos M2 ya estimados: no ajusta ni modifica ningún modelo.
# Se distinguen la escala numérica observada, la respuesta latente ordinal,
# y la fiabilidad de la suma de los tres ítems en sus categorías originales.
.libPaths(c('/Users/luisbenites/Documents/Codex/Paper con MA - 1/reproduccion-integral-V3-2026-10-10/proyecto/revision-revista-2026-10-04/auditoria-codigos/R-library',
            '/Users/luisbenites/Library/R/arm64/4.6/library', .libPaths()))
suppressPackageStartupMessages(library(lavaan))
suppressPackageStartupMessages(library(semTools))
suppressPackageStartupMessages(library(jsonlite))
base <- '/Users/luisbenites/Documents/Codex/Paper con MA - 1/reproduccion-integral-V3-2026-10-10/proyecto/sem-alternativas-2026-10-08/resultados'
out <- '/Users/luisbenites/Documents/Codex/Paper con MA - 1/reproduccion-integral-V3-2026-10-10/proyecto/sem-items-ordinal-2026-10-08/ordinal'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
ids <- c('CFA_MLR','CFA_WLSMV','SEM_MLR','SEM_WLSMV')
items <- paste0('PP',1:3)
sha <- function(x) unname(system2('/usr/bin/shasum', c('-a','256',shQuote(x)),stdout=TRUE))
safe <- function(expr) tryCatch(expr,error=function(e)structure(list(error=conditionMessage(e)),class='audit_error'))
eigen_diag <- function(x) {
  if(inherits(x,'audit_error'))return(list(error=x$error))
  v <- eigen((x+t(x))/2,symmetric=TRUE,only.values=TRUE)$values
  list(min=min(v),max=max(v),condition=max(v)/min(v),positive=min(v)>0)
}
flatten_rel <- function(x) {
  if(inherits(x,'audit_error')) return(NA_real_)
  if(is.list(x))return(as.numeric(x$PSB)[1])
  as.numeric(x['PSB'])[1]
}
all_diagnostics <- list(); measure <- list(); fits <- list(); loadings <- list(); paths <- list()
freq <- list(); fit_rows <- list(); checks <- list(); expected <- list()
for(id in ids){
  path <- file.path(base,paste0('M2_puentes2_',id,'.rds'))
  original_hash <- sha(path)
  f <- readRDS(path)
  ordinal <- grepl('WLSMV',id)
  est <- lavInspect(f,'est'); options <- lavInspect(f,'options')
  par <- parameterEstimates(f,standardized=TRUE,ci=TRUE)
  std <- standardizedSolution(f,ci=TRUE)
  if(!'label' %in% names(std)) std$label <- ''
  std <- std[,c('lhs','op','rhs','label','est.std','se','z','pvalue','ci.lower','ci.upper')]
  pt <- parTable(f); free <- pt$free>0
  # Las estimaciones y SE exportadas incluyen parámetros definidos, pero
  # admisibilidad/SE se evalúan sólo para parámetros libres estimados.
  free_ok <- all(is.finite(pt$est[free])) && all(is.finite(pt$se[free])) && all(pt$se[free]>0)
  lambda <- est$lambda; psi <- est$psi
  beta <- if(is.null(est$beta)) matrix(0,nrow(psi),ncol(psi),dimnames=dimnames(psi)) else est$beta
  U <- solve(diag(nrow(beta))-beta)
  phi_cond <- U %*% psi %*% t(U)
  phi_marg <- phi_cond
  if(!is.null(est$gamma)) phi_marg <- U %*% (psi + est$gamma %*% est$cov.x %*% t(est$gamma)) %*% t(U)
  dimnames(phi_marg) <- dimnames(phi_cond) <- dimnames(psi)
  sigma_marg <- lambda %*% phi_marg %*% t(lambda)+est$theta
  sigma_cond <- lambda %*% phi_cond %*% t(lambda)+est$theta
  covlv <- lavInspect(f,'cov.lv'); lv <- colnames(covlv)
  # Si lavaan entrega una covarianza marginal o condicional, se documenta
  # la igualdad; nunca se confunde la varianza residual PSB con su total.
  diff_lv_marg <- max(abs(covlv-phi_marg[lv,lv,drop=FALSE]))
  diff_lv_cond <- max(abs(covlv-phi_cond[lv,lv,drop=FALSE]))
  idx <- match(items,rownames(lambda)); k <- match('PSB',colnames(lambda)); lam <- lambda[idx,k]
  ave_model_marg <- sum(lam^2*phi_marg[k,k])/sum(diag(sigma_marg)[idx])
  rel_model_marg <- sum(lam)^2*phi_marg[k,k]/sum(sigma_marg[idx,idx])
  load_std_marg <- lam*sqrt(phi_marg[k,k])/sqrt(diag(sigma_marg)[idx])
  # semTools usa varianzas del modelo saturado de referencia por defecto.
  # Sólo se utiliza su AVE si la función soporta el objeto; en SEM WLSMV
  # conditional.x, el bloque sampstat es res.cov, no cov, y falla su API.
  ave_obs <- safe(semTools::AVE(f,obs.var=TRUE))
  ave_imp <- safe(semTools::AVE(f,obs.var=FALSE))
  rel_obs_ord <- safe(semTools::compRelSEM(f,obs.var=TRUE,ord.scale=TRUE))
  rel_imp_ord <- safe(semTools::compRelSEM(f,obs.var=FALSE,ord.scale=TRUE))
  rel_obs_lat <- safe(semTools::compRelSEM(f,obs.var=TRUE,ord.scale=FALSE))
  rel_imp_lat <- safe(semTools::compRelSEM(f,obs.var=FALSE,ord.scale=FALSE))
  ave_obs_psb <- if(inherits(ave_obs,'audit_error'))NA_real_ else as.numeric(ave_obs['PSB'])
  ave_imp_psb <- if(inherits(ave_imp,'audit_error'))NA_real_ else as.numeric(ave_imp['PSB'])
  info <- safe(lavInspect(f,'information'))
  hessian <- safe(lavInspect(f,'hessian'))
  vcov_est <- safe(lavInspect(f,'vcov'))
  post <- safe(lavInspect(f,'post.check'))
  gradients <- safe(lavInspect(f,'gradient'))
  information <- eigen_diag(info); hessdiag <- eigen_diag(hessian); vcdiag <- eigen_diag(vcov_est)
  psidiag <- eigen_diag(psi); margdiag <- eigen_diag(phi_marg); sigdiag <- eigen_diag(sigma_marg)
  # Escalas ordinales: los umbrales deben estar en orden estricto por ítem.
  thresholds <- par[par$op=='|',c('lhs','rhs','est','se','ci.lower','ci.upper')]
  threshold_order <- if(nrow(thresholds)) all(vapply(split(thresholds$est,thresholds$lhs),function(x)all(diff(x)>0),logical(1))) else TRUE
  stdpsb <- std[std$lhs=='PSB'&std$op=='=~',]
  rownames(stdpsb)<-NULL
  loadings[[id]] <- cbind(model=id, stdpsb, matrix_marginal_std=load_std_marg,
                         matrix_conditional_std=lam*sqrt(phi_cond[k,k])/sqrt(diag(sigma_cond)[idx]))
  focal_rhs <- c('SC','SBL','Benevolencia','Universalismo','Poder','Logro','CSR1_bin')
  pf <- std[std$lhs=='PSB'&std$op=='~'&std$rhs %in% focal_rhs,]
  # Familia reproducible: los mismos 21 caminos sustantivos del modelo M2.
  all_focal <- std[(std$op=='~' & ((std$lhs %in% c('Benevolencia','Universalismo','Poder','Logro') & std$rhs=='CSR1_bin') |
                   (std$lhs %in% c('SC','SBL') & std$rhs %in% c('Benevolencia','Universalismo','Poder','Logro','CSR1_bin')) |
                   (std$lhs=='PSB' & std$rhs %in% focal_rhs))),]
  stopifnot(nrow(all_focal)==if(grepl('SEM',id))21L else 0L)
  # Contrastes/Holm usan Wald robusto del parámetro no estandarizado;
  # los IC de beta estandarizada se conservan como IC puntuales de descripción.
  parkeys <- paste(par$lhs,par$op,par$rhs)
  all_focal$pvalue_unstandardized <- par$pvalue[match(paste(all_focal$lhs,all_focal$op,all_focal$rhs),parkeys)]
  all_focal$p_Holm_21 <- p.adjust(all_focal$pvalue_unstandardized,method='holm')
  write.csv(cbind(model=rep(id,nrow(all_focal)),all_focal),file.path(out,paste0(id,'_rutas_21_Holm.csv')),row.names=FALSE)
  pf$p_Holm_21 <- all_focal$p_Holm_21[match(paste(pf$lhs,pf$op,pf$rhs),paste(all_focal$lhs,all_focal$op,all_focal$rhs))]
  pf$pvalue_unstandardized <- par$pvalue[match(paste(pf$lhs,pf$op,pf$rhs),parkeys)]
  paths[[id]] <- cbind(model=rep(id,nrow(pf)),pf)
  fit_stats <- fitMeasures(f)
  fit_rows[[id]] <- data.frame(model=id,metric=names(fit_stats),value=as.numeric(fit_stats))
  measure[[id]] <- data.frame(model=id,ordinal=ordinal,
      response_scale=if(ordinal)'underlying ordinal response' else 'observed numeric score',
      AVE_semTools_observed_reference=ave_obs_psb,AVE_semTools_model_implied=ave_imp_psb,
      AVE_matrix_model_marginal=ave_model_marg,
      reliability_matrix_model_marginal=rel_model_marg,
      reliability_semTools_observed_score_reference=flatten_rel(rel_obs_ord),
      reliability_semTools_observed_score_model=flatten_rel(rel_imp_ord),
      reliability_semTools_latent_response_reference=flatten_rel(rel_obs_lat),
      reliability_semTools_latent_response_model=flatten_rel(rel_imp_lat),
      corr_SC_SBL_lavaan=lavInspect(f,'cor.lv')['SC','SBL'],
      corr_SC_SBL_marginal=cov2cor(phi_marg)['SC','SBL'],
      corr_SC_SBL_conditional=cov2cor(phi_cond)['SC','SBL'],
      covlv_max_difference_marginal=diff_lv_marg,covlv_max_difference_conditional=diff_lv_cond)
  all_diagnostics[[id]] <- list(source_file=path,source_hash_before=original_hash,
     estimator_original=options$estimator.orig,estimator_internal=options$estimator,
     parameterization=options$parameterization,conditional_x=options$conditional.x,fixed_x=options$fixed.x,
     missing=options$missing,se=options$se,test=options$test,nobs=lavInspect(f,'nobs'),
     converged=lavInspect(f,'converged'),post_check=post,optimizer_warning=f@optim$warn.txt,
     all_free_parameters_estimate_SE_finite_positive=free_ok,
     free_count=sum(free),thresholds_strictly_ordered=threshold_order,
     max_absolute_gradient=if(inherits(gradients,'audit_error'))NA_real_ else max(abs(gradients)),
     information=information,hessian=hessdiag,robust_parameter_covariance=vcdiag,
     residual_factor_covariance=psidiag,marginal_factor_covariance=margdiag,
     implied_response_covariance=sigdiag,indicator_residual_variance_min=min(diag(est$theta)),
     semTools_AVE_obsvar_error=if(inherits(ave_obs,'audit_error'))ave_obs$error else NULL,
     semTools_AVE_implied_error=if(inherits(ave_imp,'audit_error'))ave_imp$error else NULL,
     semTools_reliability_error=if(inherits(rel_obs_ord,'audit_error'))rel_obs_ord$error else NULL,
     source_hash_after=sha(path))
  data <- lavInspect(f,'data')
  if(ordinal){
     freq[[id]] <- do.call(rbind,lapply(lavNames(f,'ov.ord'),function(n){
       tb <- table(data[,n]); data.frame(model=id,item=n,category=names(tb),n=as.integer(tb),pct=as.integer(tb)/nrow(data)*100)
     }))
     write.csv(thresholds,file.path(out,paste0(id,'umbrales.csv')),row.names=FALSE)
  }
  checks[[id]] <- list(converged=lavInspect(f,'converged'),post_check=isTRUE(post),
      finite_SE=free_ok,information_positive=isTRUE(information$positive),hessian_positive=isTRUE(hessdiag$positive),
      vcov_positive=isTRUE(vcdiag$positive),phi_positive=isTRUE(margdiag$positive),sigma_positive=isTRUE(sigdiag$positive),
      residuals_nonnegative=min(diag(est$theta))>=0,
      ordered_thresholds=threshold_order,unchanged_source=identical(original_hash,sha(path)))
  write.csv(par,file.path(out,paste0(id,'_parametros_CI.csv')),row.names=FALSE)
  write.csv(std,file.path(out,paste0(id,'_estandarizados_CI.csv')),row.names=FALSE)
}
write.csv(do.call(rbind,measure),file.path(out,'PSB_AVE_fiabilidad_escalas.csv'),row.names=FALSE)
write.csv(do.call(rbind,loadings),file.path(out,'PSB_cargas_CI.csv'),row.names=FALSE)
write.csv(do.call(rbind,paths),file.path(out,'PSB_rutas_CI_Holm.csv'),row.names=FALSE)
write.csv(do.call(rbind,freq),file.path(out,'frecuencias_ordinales_agregadas.csv'),row.names=FALSE)
write.csv(do.call(rbind,fit_rows),file.path(out,'indices_ajuste_por_estimador.csv'),row.names=FALSE)
write_json(all_diagnostics,file.path(out,'diagnosticos_objetos_existentes.json'),pretty=TRUE,auto_unbox=TRUE,digits=16,null='null',na='null')
write_json(list(checks=checks,all_pass=all(unlist(checks))),file.path(out,'auditoria_estado.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
print(do.call(rbind,measure),digits=12)
print(do.call(rbind,paths)[,c('model','rhs','est.std','ci.lower','ci.upper','p_Holm_21')],digits=9)
print(list(all_pass=all(unlist(checks))))
