# Extensión exploratoria: medición ESEM multigrupo por país.
# Mantiene los 653 casos y los 22 indicadores de la comparación anterior.
# La matriz target guía la rotación; no impone cargas iguales entre países.
# No prueba diferencias de rutas ni de medias latentes entre países.
options(stringsAsFactors=FALSE, warn=1)
Sys.setenv(VECLIB_MAXIMUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
project <- '/Users/luisbenites/Documents/Codex/Paper con MA - 1/reproduccion-integral-V3-2026-10-10/proyecto'
root <- file.path(project,'sem-multigrupo-2026-10-08')
prior <- file.path(project,'sem-alternativas-2026-10-08')
.libPaths(c(file.path(project,'revision-revista-2026-10-04/auditoria-codigos/R-library'),
            '/Users/luisbenites/Library/R/arm64/4.6/library',.libPaths()))
suppressPackageStartupMessages({library(lavaan);library(jsonlite);library(digest)})
out <- file.path(root,'resultados')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
state <- file.path(out,'estado.txt')
writeLines('STARTED',state)
options(error=function(){writeLines(c('FAILED',geterrmessage()),state);quit(status=1,save='no')})
write_table <- function(x,name)write.csv(x,file.path(out,name),row.names=FALSE,na='')

# Recuperar exactamente las recodificaciones del modelo anterior.
data_path <- file.path(project,'revision-integral-2026-10-04/analysis/outputs/canonical-full-20261004-02/01_audit/analytic_data.rds')
p3_path <- file.path(project,'revision-integral-2026-10-04/sensitivities/P3_country_partial_demographics.rds')
raw_path <- '/Users/luisbenites/Dropbox/Paper con MA/Base_Final.sav'
prior_fit_path <- file.path(prior,'resultados/model_fit.csv')
protected_paths <- c(data_path,p3_path,raw_path,prior_fit_path,
                    file.path(prior,'documentos/Manuscript_SEM_comparison_JHV.docx'),
                    file.path(prior,'documentos/Supplementary_SEM_comparison_JHV.docx'))
hash_before <- vapply(protected_paths,digest,character(1),algo='sha256',file=TRUE)
audited <- as.data.frame(readRDS(data_path))
p3 <- readRDS(p3_path)
dat <- as.data.frame(lavInspect(p3,'data'))
shared <- intersect(names(dat),names(audited))
stopifnot(max(abs(as.matrix(dat[,shared])-as.matrix(audited[,shared])))==0)
dat$Pais <- audited$Pais
sc <- c('APS1','APS2','APS3','APS4','APS5','APS6','APS7','APS9','APS10','APS14')
sbl <- c('APS11','APS12','APS13','APS15','APS16','APS18','APS19','APS20','APS21')
aps <- c(sc,sbl)
pp <- c('PP1','PP2','PP3')
stopifnot(nrow(dat)==653L,!anyNA(dat[,c(aps,pp,'Pais')]))
dat[,c(aps,pp)] <- lapply(dat[,c(aps,pp)],as.numeric)

# Usar la misma sintaxis de medición M6 que antes, sin las regresiones P3.
syntax <- paste(readLines(file.path(prior,'resultados/M6_ESEM_target_syntax.txt'),warn=FALSE)[1:2],collapse='\n')
stopifnot(grepl('efa("asp")',syntax,fixed=TRUE),!grepl(' ~ ',syntax,fixed=TRUE))
writeLines(syntax,file.path(out,'M6_configural_syntax.txt'))
target <- matrix(0,nrow=length(aps),ncol=2,dimnames=list(aps,c('SC','SBL')))
target[sc,'SC'] <- 1
target[sbl,'SBL'] <- 1
seed <- 2026100809L
protocol <- list(timestamp_UTC=format(Sys.time(),tz='UTC'),
  status='Exploratory extension fixed before this new fit, after prior full-sample analyses',
  question='Can the existing M6 target ESEM provide a more adequate configural model across countries?',
  model='M6 ESEM ASP two factors, plus separate PSB CFA',N=653,
  group_variable='Pais',indicators_retained=c(aps,pp),estimator='MLR',
  group_equal=character(0),group_partial=character(0),rotation='target',
  random_starts=20,seed=seed,missing='fiml',std_lv=TRUE,
  no_item_deletion=TRUE,no_MI_search=TRUE,no_residual_correlation_additions=TRUE,
  no_structural_country_comparisons=TRUE,no_metric_or_scalar_claim=TRUE,
  interpretation='Configuration and residual diagnostics; target rotation is not measurement invariance',
  original_paths=protected_paths,original_sha256=unname(hash_before))
write_json(protocol,file.path(out,'protocol_before_fit.json'),pretty=TRUE,auto_unbox=TRUE)

set.seed(seed)
warnings <- character()
cat('FITTING M6 ESEM CONFIGURAL BY Pais, MLR\n');flush.console()
fit <- tryCatch(withCallingHandlers(
  sem(syntax,data=dat,group='Pais',estimator='MLR',meanstructure=TRUE,
      std.lv=TRUE,fixed.x=TRUE,missing='fiml',control=list(iter.max=2000),
      rotation='target',rotation.args=list(target=target,rstarts=20,orthogonal=FALSE)),
  warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),
  error=function(e)e)
write_json(list(warnings=unique(warnings),error=if(inherits(fit,'error'))conditionMessage(fit)else NULL),
           file.path(out,'fit_messages.json'),pretty=TRUE,auto_unbox=TRUE)
if(!inherits(fit,'lavaan'))stop(conditionMessage(fit))

# El objeto lavaan contiene respuestas; permanece privado en la carpeta local.
saveRDS(fit,file.path(out,'M6_configural_private.rds'))
fm <- fitMeasures(fit)
measure_names <- c('chisq.scaled','df','cfi.robust','tli.robust','rmsea.robust',
                  'rmsea.ci.lower.robust','rmsea.ci.upper.robust','srmr','aic','bic','npar')
new_fit <- cbind(data.frame(model='M6_ESEM_target',kind='configural',estimator='MLR',
                           N=sum(lavInspect(fit,'nobs')),converged=lavInspect(fit,'converged'),
                           post_check=lavInspect(fit,'post.check'),warnings=paste(unique(warnings),collapse=' | ')),
                 as.data.frame(as.list(fm[measure_names])))
old_fit <- read.csv(prior_fit_path,check.names=FALSE)
old_fit <- old_fit[old_fit$kind=='configural',names(new_fit),drop=FALSE]
comparison <- rbind(old_fit,new_fit)
write_table(comparison,'comparacion_configural.csv')
cat('FIT INDICES\n');print(new_fit,row.names=FALSE);flush.console()

# Guardar solo parámetros y diagnósticos agregados, sin respuestas individuales.
ss <- standardizedSolution(fit)
pe <- parameterEstimates(fit,standardized=TRUE)
group_labels <- lavInspect(fit,'group.label')
ss$Pais <- group_labels[ss$group]
pe$Pais <- group_labels[pe$group]
write_table(ss,'M6_configural_standardized.csv')
write_table(pe,'M6_configural_parameters.csv')
group_ns <- lavInspect(fit,'nobs')
cors <- lavInspect(fit,'cor.lv')
thetas <- lavInspect(fit,'theta')
matrices <- lavInspect(fit,'est')
eig <- function(x)min(eigen(x,symmetric=TRUE,only.values=TRUE)$values)
diagnostics <- lapply(seq_along(group_labels),function(g){
  loads <- ss[ss$group==g & ss$op=='=~',]
  primary <- loads[(loads$lhs=='SC' & loads$rhs%in%sc)|(loads$lhs=='SBL' & loads$rhs%in%sbl),]
  primary_dominant <- vapply(seq_len(nrow(primary)),function(i){
    r <- primary[i,]
    other <- loads[loads$rhs==r$rhs & loads$lhs%in%c('SC','SBL') & loads$lhs!=r$lhs,]
    nrow(other)==1L && abs(r$est.std)>abs(other$est.std)
  },logical(1))
  variances <- pe$est[pe$group==g & pe$op=='~~' & pe$lhs==pe$rhs]
  data.frame(Pais=group_labels[g],N=group_ns[g],
             SC_SBL_correlation=cors[[g]]['SC','SBL'],
             min_variance=min(variances),min_theta_eigen=eig(thetas[[g]]),
             min_psi_eigen=eig(matrices[[g]]$psi),
             primary_positive=sum(primary$est.std>0),primary_total=nrow(primary),
             primary_absolute_dominant=sum(primary_dominant),
             all_standardized_SE_finite=all(is.finite(ss$se[ss$group==g])))
})
diag <- do.call(rbind,diagnostics)
write_table(diag,'diagnostico_por_pais.csv')
cat('COUNTRY DIAGNOSTICS\n');print(diag,row.names=FALSE);flush.console()

# Congruencia descriptiva entre patrones SC/SBL; no es prueba de invariancia.
get_pattern <- function(g,f) {
  z <- ss[ss$group==g & ss$op=='=~' & ss$lhs==f,]
  setNames(z$est.std,z$rhs)[aps]
}
pairs <- combn(seq_along(group_labels),2,simplify=FALSE)
congruence <- do.call(rbind,lapply(pairs,function(pair)do.call(rbind,lapply(c('SC','SBL'),function(f){
  a <- get_pattern(pair[1],f);b <- get_pattern(pair[2],f)
  data.frame(Pais_1=group_labels[pair[1]],Pais_2=group_labels[pair[2]],factor=f,
             descriptive_cosine=sum(a*b)/sqrt(sum(a^2)*sum(b^2)))
}))))
write_table(congruence,'congruencia_descriptiva.csv')
writeLines(capture.output(summary(fit,fit.measures=TRUE,standardized=TRUE)),file.path(out,'M6_configural_summary.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
hash_after <- vapply(protected_paths,digest,character(1),algo='sha256',file=TRUE)
stopifnot(identical(hash_before,hash_after))
write_json(list(paths=protected_paths,sha256_before=unname(hash_before),
                sha256_after=unname(hash_after),all_originals_unchanged=TRUE),
           file.path(out,'originals_hash_audit.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines('COMPLETED',state)
cat('COMPLETED. Metric/scalar invariance and country path comparisons were not tested.\n')
