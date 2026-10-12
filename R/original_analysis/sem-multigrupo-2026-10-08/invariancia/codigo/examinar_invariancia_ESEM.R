# Invariancia del ESEM de medición entre Argentina, Colombia, Perú y México.
# Extensión exploratoria del análisis ya realizado; no es un prerregistro.
# Se conservan todos los casos/ítems y no se liberan restricciones por índices MI.
options(stringsAsFactors=FALSE,warn=1)
Sys.setenv(VECLIB_MAXIMUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
project <- '/Users/luisbenites/Documents/Codex/Paper con MA - 1/reproduccion-integral-V3-2026-10-10/proyecto'
parent <- file.path(project,'sem-multigrupo-2026-10-08')
root <- file.path(parent,'invariancia')
out <- file.path(root,'resultados')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
.libPaths(c(file.path(project,'revision-revista-2026-10-04/auditoria-codigos/R-library'),
            '/Users/luisbenites/Library/R/arm64/4.6/library',.libPaths()))
suppressPackageStartupMessages({library(lavaan);library(jsonlite);library(digest)})
state <- file.path(out,'estado.txt')
writeLines('STARTED',state)
options(error=function(){writeLines(c('FAILED',geterrmessage()),state);quit(status=1,save='no')})
csv <- function(x,name)write.csv(x,file.path(out,name),row.names=FALSE,na='')
json <- function(x,name)write_json(x,file.path(out,name),pretty=TRUE,auto_unbox=TRUE,null='null')
data_path <- file.path(project,'revision-integral-2026-10-04/analysis/outputs/canonical-full-20261004-02/01_audit/analytic_data.rds')
p3_path <- file.path(project,'revision-integral-2026-10-04/sensitivities/P3_country_partial_demographics.rds')
reference_path <- file.path(parent,'resultados/M6_configural_private.rds')
protected <- c(data_path,p3_path,reference_path,
               '/Users/luisbenites/Dropbox/Paper con MA/Base_Final.sav',
               file.path(project,'sem-alternativas-2026-10-08/documentos/Manuscript_SEM_comparison_JHV.docx'),
               file.path(project,'sem-alternativas-2026-10-08/documentos/Supplementary_SEM_comparison_JHV.docx'))
before <- vapply(protected,digest,character(1),algo='sha256',file=TRUE)
dat_full <- as.data.frame(readRDS(data_path))
p3 <- readRDS(p3_path)
dat <- as.data.frame(lavInspect(p3,'data'))
shared <- intersect(names(dat),names(dat_full))
stopifnot(max(abs(as.matrix(dat[,shared])-as.matrix(dat_full[,shared])))==0)
dat$Pais <- dat_full$Pais
sc <- c('APS1','APS2','APS3','APS4','APS5','APS6','APS7','APS9','APS10','APS14')
sbl <- c('APS11','APS12','APS13','APS15','APS16','APS18','APS19','APS20','APS21')
aps <- c(sc,sbl);pp <- c('PP1','PP2','PP3');indicators <- c(aps,pp)
dat <- dat[,c(indicators,'Pais')]
dat[,indicators] <- lapply(dat[,indicators],as.numeric)
stopifnot(nrow(dat)==653L,!anyNA(dat))
syntax <- paste(readLines(file.path(parent,'resultados/M6_configural_syntax.txt'),warn=FALSE),collapse='\n')
target <- matrix(0,19,2,dimnames=list(aps,c('SC','SBL')))
target[sc,'SC'] <- 1;target[sbl,'SBL'] <- 1
reference <- readRDS(reference_path)

# Congelar antes de nuevos ajustes dos niveles completos de igualdad y dos
# réplicas numéricas del configural. La invariancia parcial nativa se excluye:
# su interacción con la rotación en esta versión requiere una formulación propia.
candidates <- list(metric=c('loadings'),scalar=c('loadings','intercepts'))
seeds <- c(metric=2026100812L,scalar=2026100813L,replica1=2026100810L,replica2=2026100811L)
protocol <- list(timestamp_UTC=format(Sys.time(),tz='UTC'),
  status='Prospective extension after observed configural results; exploratory, not preregistered',
  N=653,group='Pais',indicators=indicators,lavaan=as.character(packageVersion('lavaan')),
  estimator='MLR',missing='fiml',std_lv=TRUE,rotation='target',rstarts=20,
  constraints=candidates,seeds=as.list(seeds),
  interpretation='Joint assessment of robust LRT, delta fit, residuals, identification, numerical equality and content',
  exact_invariance='Evaluate robust nested likelihood-ratio tests, report failures and disagreements',
  approximate_invariance='Delta CFI/RMSEA/SRMR are descriptive and do not by themselves prove invariance; no automatic cutoff-based selection',
  partial_native_excluded=TRUE,no_MI_search=TRUE,no_item_deletion=TRUE,
  no_structural_or_mean_country_contrasts_without_measurement_support=TRUE,
  original_sha256=as.list(before),script_sha256=digest(file.path(root,'codigo/examinar_invariancia_ESEM.R'),algo='sha256',file=TRUE))
json(protocol,'protocol_before_fits.json')

safe_fit <- function(equal,seed) {
  set.seed(seed);w <- character()
  z <- tryCatch(withCallingHandlers(
    sem(syntax,data=dat,group='Pais',estimator='MLR',meanstructure=TRUE,
        std.lv=TRUE,fixed.x=TRUE,missing='fiml',group.equal=equal,
        rotation='target',rotation.args=list(target=target,rstarts=20,orthogonal=FALSE),
        control=list(iter.max=2000)),
    warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart('muffleWarning')}),
    error=function(e)e)
  list(fit=z,warnings=unique(w),error=if(inherits(z,'error'))conditionMessage(z)else NULL)
}
measures <- c('chisq','chisq.scaled','df','cfi.robust','tli.robust','rmsea.robust',
              'rmsea.ci.lower.robust','rmsea.ci.upper.robust','srmr','aic','bic','npar','logl')
fit_stats <- function(f,name,warnings=character()) {
  if(!inherits(f,'lavaan'))return(data.frame(model=name,error=conditionMessage(f)))
  fm <- fitMeasures(f)
  values <- setNames(vapply(measures,function(m)if(m%in%names(fm))unname(fm[m])else NA_real_,numeric(1)),measures)
  cbind(data.frame(model=name,converged=lavInspect(f,'converged'),post_check=lavInspect(f,'post.check'),
                   warnings=paste(warnings,collapse=' | '),error=''),as.data.frame(as.list(values)))
}
bind <- function(rows) {
  cols <- unique(unlist(lapply(rows,names)))
  do.call(rbind,lapply(rows,function(r){r[setdiff(cols,names(r))]<-NA;r[,cols,drop=FALSE]}))
}
fits <- list(configural=reference);messages <- list();rows <- list(configural=fit_stats(reference,'configural'))
for(name in names(candidates)) {
  cat('FITTING',name,'\n');flush.console()
  z <- safe_fit(candidates[[name]],seeds[name])
  fits[[name]] <- z$fit;messages[[name]]<-z[c('warnings','error')]
  if(inherits(z$fit,'lavaan'))saveRDS(z$fit,file.path(out,paste0(name,'_private.rds')))
  rows[[name]] <- fit_stats(z$fit,name,z$warnings)
  csv(bind(rows),'comparacion_invariancia.csv');json(messages,'fit_messages.json')
  print(rows[[name]],row.names=FALSE);flush.console()
}

# Las réplicas contrastan estabilidad del ajuste y de la orientación target.
for(name in c('replica1','replica2')) {
  cat('FITTING',name,'\n');flush.console()
  z <- safe_fit(character(),seeds[name])
  fits[[name]] <- z$fit;messages[[name]] <- z[c('warnings','error')]
  if(inherits(z$fit,'lavaan'))saveRDS(z$fit,file.path(out,paste0(name,'_private.rds')))
  rows[[name]] <- fit_stats(z$fit,name,z$warnings)
  csv(bind(rows),'comparacion_invariancia.csv');json(messages,'fit_messages.json')
}
saveRDS(fits,file.path(out,'all_fits_private.rds'))

# La igualdad debe comprobarse en cargas NO estandarizadas. Sus versiones
# estandarizadas pueden diferir al variar la dispersión del factor o del ítem.
equalities <- list();diag <- list();stability <- list()
eig <- function(m)min(eigen(m,symmetric=TRUE,only.values=TRUE)$values)
for(name in names(fits)) {
  f <- fits[[name]];if(!inherits(f,'lavaan'))next
  pe <- parameterEstimates(f,standardized=TRUE);ss <- standardizedSolution(f)
  labels <- lavInspect(f,'group.label');pe$Pais<-labels[pe$group];ss$Pais<-labels[ss$group]
  csv(pe,paste0(name,'_parameters.csv'));csv(ss,paste0(name,'_standardized.csv'))
  writeLines(capture.output(summary(f,fit.measures=TRUE,standardized=TRUE)),file.path(out,paste0(name,'_summary.txt')))
  mm <- lavInspect(f,'est');implied <- lavInspect(f,'implied');thetas <- lavInspect(f,'theta')
  lmax <- max(vapply(mm[-1],function(m)max(abs(m$lambda-mm[[1]]$lambda)),numeric(1)))
  imax <- max(vapply(mm[-1],function(m)max(abs(m$nu-mm[[1]]$nu)),numeric(1)))
  equalities[[name]] <- data.frame(model=name,max_unstandardized_loading_difference=lmax,
                                 max_item_intercept_difference=imax,
                                 effective_npar=fitMeasures(f,'npar'))
  diag[[name]] <- do.call(rbind,lapply(seq_along(labels),function(g){
    m<-mm[[g]];v<-pe$est[pe$group==g & pe$op=='~~' & pe$lhs==pe$rhs]
    data.frame(model=name,Pais=labels[g],N=lavInspect(f,'nobs')[g],min_variance=min(v),
               min_theta_eigen=eig(thetas[[g]]),min_psi_eigen=eig(m$psi),min_implied_eigen=eig(implied[[g]]$cov),
               all_SE_finite=all(is.finite(pe$se[pe$group==g])),
               all_CI_finite=all(is.finite(pe$ci.lower[pe$group==g])&is.finite(pe$ci.upper[pe$group==g])),
               SC_variance=m$psi['SC','SC'],SBL_variance=m$psi['SBL','SBL'],PSB_variance=m$psi['PSB','PSB'],
               SC_mean=m$alpha['SC',1],SBL_mean=m$alpha['SBL',1],PSB_mean=m$alpha['PSB',1])
  }))
  if(name%in%c('replica1','replica2')) {
    mr<-lavInspect(reference,'est');ir<-lavInspect(reference,'implied')
    stability[[name]] <- data.frame(model=name,loglik_difference=fitMeasures(f,'logl')-fitMeasures(reference,'logl'),
      max_lambda_difference=max(vapply(seq_along(mm),function(g)max(abs(mm[[g]]$lambda-mr[[g]]$lambda)),numeric(1))),
      max_covariance_difference=max(vapply(seq_along(mm),function(g)max(abs(implied[[g]]$cov-ir[[g]]$cov)),numeric(1))))
  }
}
csv(bind(equalities),'igualdades_numericas.csv')
csv(bind(diag),'diagnosticos_por_pais.csv')
csv(bind(stability),'estabilidad_configural.csv')

# MLR requiere una diferencia robusta calculada conjuntamente: no se restan
# directamente los chi-cuadrados escalados de los dos modelos.
lrt_results <- list();lrt_messages<-list()
for(pair in list(c('configural','metric'),c('metric','scalar'))) {
  if(!all(vapply(fits[pair],inherits,logical(1),'lavaan')))next
  for(method in c('satorra.bentler.2001','satorra.bentler.2010')) {
    key<-paste(c(pair,method),collapse='_');w<-character()
    test<-tryCatch(withCallingHandlers(
      lavTestLRT(fits[[pair[1]]],fits[[pair[2]]],method=method),
      warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart('muffleWarning')}),error=function(e)e)
    lrt_messages[[key]]<-list(warnings=unique(w),error=if(inherits(test,'error'))conditionMessage(test)else NULL)
    if(!inherits(test,'error')) {
      z<-as.data.frame(test);z$model<-rownames(z);z$comparison<-paste(pair,collapse='_to_');z$method<-method
      rownames(z)<-NULL;lrt_results[[key]]<-z
      writeLines(capture.output(test),file.path(out,paste0('LRT_',key,'.txt')))
    }
  }
}
if(length(lrt_results))csv(bind(lrt_results),'robust_LRT.csv')
json(lrt_messages,'LRT_messages.json')
deltas <- lapply(list(c('configural','metric'),c('metric','scalar')),function(pair){
  if(!all(vapply(fits[pair],inherits,logical(1),'lavaan')))return(NULL)
  a<-fitMeasures(fits[[pair[1]]]);b<-fitMeasures(fits[[pair[2]]])
  data.frame(comparison=paste(pair,collapse='_to_'),delta_CFI_robust=b['cfi.robust']-a['cfi.robust'],
             delta_TLI_robust=b['tli.robust']-a['tli.robust'],delta_RMSEA_robust=b['rmsea.robust']-a['rmsea.robust'],
             delta_SRMR=b['srmr']-a['srmr'],delta_df=b['df']-a['df'])
})
csv(do.call(rbind,deltas),'cambios_ajuste.csv')
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
after<-vapply(protected,digest,character(1),algo='sha256',file=TRUE)
stopifnot(identical(before,after))
json(list(paths=protected,sha256_before=unname(before),sha256_after=unname(after),all_originals_unchanged=TRUE),'originals_hash_audit.json')
writeLines('COMPLETED',state)
cat('COMPLETED: invariance comparisons and stability checks saved; interpretation requires independent audit.\n')
