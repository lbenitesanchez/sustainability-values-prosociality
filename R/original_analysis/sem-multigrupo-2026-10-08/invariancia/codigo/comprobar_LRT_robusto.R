# Comprobación adicional del contraste robusto de invariancia.
# No se ajustan modelos nuevos: se usan los tres objetos ya estimados.
# ESEM usa identificaciones rotacionales diferentes entre niveles; por eso
# se comprueba también el anidamiento en medias/covarianzas con método delta.
options(stringsAsFactors=FALSE,warn=1)
Sys.setenv(VECLIB_MAXIMUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
project <- '/Users/luisbenites/Documents/Codex/Paper con MA - 1/reproduccion-integral-V3-2026-10-10/proyecto'
root <- file.path(project,'sem-multigrupo-2026-10-08/invariancia')
out <- file.path(root,'resultados')
.libPaths(c(file.path(project,'revision-revista-2026-10-04/auditoria-codigos/R-library'),
            '/Users/luisbenites/Library/R/arm64/4.6/library',.libPaths()))
suppressPackageStartupMessages({library(lavaan);library(jsonlite);library(digest)})
stopifnot(readLines(file.path(out,'estado.txt'))=='COMPLETED')
fits <- readRDS(file.path(out,'all_fits_private.rds'))
stopifnot(all(vapply(fits[c('configural','metric','scalar')],inherits,logical(1),'lavaan')))
protocol <- list(timestamp_UTC=format(Sys.time(),tz='UTC'),
  status='Supplementary check after initial SB2001/SB2010 outputs, fixed before new delta tests',
  motivation='ESEM rotational parameterisation changes; initial metric SB2001 was borderline and SB2010 highly different',
  methods=c('satorra.bentler.2001 explicit Yuan Bentler Mplus','satorra.2000 delta unshifted'),
  test='yuan.bentler.mplus',a_method='delta',scaled_shifted=FALSE,
  selection='Retain every result and warning, do not choose method by significance',
  inputs_sha256=digest(file.path(out,'all_fits_private.rds'),algo='sha256',file=TRUE),
  no_new_model_estimation=TRUE,script_sha256=digest(file.path(root,'codigo/comprobar_LRT_robusto.R'),algo='sha256',file=TRUE))
write_json(protocol,file.path(out,'protocol_before_delta_LRT.json'),pretty=TRUE,auto_unbox=TRUE,null='null')
rows <- list();messages <- list()
for(pair in list(c('configural','metric'),c('metric','scalar'))) {
  stopifnot(identical(lavInspect(fits[[pair[1]]],'nobs'),lavInspect(fits[[pair[2]]],'nobs')),
            identical(lavNames(fits[[pair[1]]],'ov'),lavNames(fits[[pair[2]]],'ov')))
  for(method in c('satorra.bentler.2001','satorra.2000')) {
    key <- paste(c(pair,method),collapse='_');w <- character()
    # Usar nombres de objetos cortos evita que lavaan deparse los objetos S4
    # completos (que contienen respuestas) como etiquetas de las filas.
    fit_a <- fits[[pair[1]]];fit_b <- fits[[pair[2]]]
    z <- tryCatch(withCallingHandlers(
      if(method=='satorra.2000')lavTestLRT(fit_a,fit_b,method=method,
        test='yuan.bentler.mplus',a_method='delta',scaled_shifted=FALSE)
      else lavTestLRT(fit_a,fit_b,method=method,test='yuan.bentler.mplus'),
      warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart('muffleWarning')}),error=function(e)e)
    messages[[key]]<-list(warnings=unique(w),error=if(inherits(z,'error'))conditionMessage(z)else NULL)
    if(!inherits(z,'error')) {
      r<-as.data.frame(z);r$model<-pair;r$comparison<-paste(pair,collapse='_to_');r$method<-method
      rownames(r)<-NULL;rows[[key]]<-r
      writeLines(capture.output(z),file.path(out,paste0('LRT_explicit_',key,'.txt')))
      cat('LRT',key,'\n');print(z);flush.console()
    }
  }
}
cols<-unique(unlist(lapply(rows,names)))
all<-do.call(rbind,lapply(rows,function(r){r[setdiff(cols,names(r))]<-NA;r[,cols,drop=FALSE]}))
write.csv(all,file.path(out,'robust_LRT_explicit_and_delta.csv'),row.names=FALSE,na='')
write_json(messages,file.path(out,'LRT_explicit_and_delta_messages.json'),pretty=TRUE,auto_unbox=TRUE,null='null')
cat('COMPLETED: no new models fitted.\n')
