# Comparación exploratoria de SEM para valores actitudes sostenibles y prosocialidad
# Los datos originales y los Word anteriores se conservan sin sobrescribir.
# La lista de alternativas se fija por contenido antes de estimar resultados nuevos.
# Este protocolo es prospectivo para esta comparación, no un prerregistro del estudio.
options(stringsAsFactors=FALSE, warn=1)
Sys.setenv(VECLIB_MAXIMUM_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
project <- '/Users/luisbenites/Documents/Codex/Paper con MA - 1/reproduccion-integral-V3-2026-10-10/proyecto'
root <- file.path(project,'sem-alternativas-2026-10-08')
.libPaths(c(file.path(project,'revision-revista-2026-10-04/auditoria-codigos/R-library'),
            '/Users/luisbenites/Library/R/arm64/4.6/library',.libPaths()))
suppressPackageStartupMessages({library(lavaan);library(semTools);library(jsonlite);library(digest)})
out <- file.path(root,'resultados')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
state <- file.path(out,'estado.txt')
writeLines('STARTED',state)
options(error=function(){writeLines(c('FAILED',geterrmessage()),state);quit(status=1,save='no')})
set.seed(20261008)

# Los registros analíticos ya fueron cotejados con la base SPSS por un auditor independiente.
data_path <- file.path(project,'revision-integral-2026-10-04/analysis/outputs/canonical-full-20261004-02/01_audit/analytic_data.rds')
p3_path <- file.path(project,'revision-integral-2026-10-04/sensitivities/P3_country_partial_demographics.rds')
audited_data <- as.data.frame(readRDS(data_path))
old_p3 <- readRDS(p3_path)
dat <- as.data.frame(lavInspect(old_p3,'data'))
# P3 conserva los 34 campos y recodificaciones exactos usados en los resultados.
# La base analítica completa conserva además el factor de sede para estratificar.
shared <- intersect(names(dat),names(audited_data))
stopifnot(max(abs(as.matrix(dat[,shared])-as.matrix(audited_data[,shared])))==0)
dat$Pais <- audited_data$Pais
sc <- c('APS1','APS2','APS3','APS4','APS5','APS6','APS7','APS9','APS10','APS14')
sbl <- c('APS11','APS12','APS13','APS15','APS16','APS18','APS19','APS20','APS21')
aps <- c(sc,sbl)
pp <- c('PP1','PP2','PP3')
values <- c('Benevolencia','Universalismo','Poder','Logro')
controls <- c('C_Colombia','C_Peru','C_Mexico','Sexo_num','Estudios_num','Edad_z','NSE_z')
ov <- c(aps,pp,values,'CSR1_bin',controls)
stopifnot(nrow(dat)==653L,all(ov%in%names(dat)),!anyNA(dat[,ov]))
dat[,ov] <- lapply(dat[,ov],as.numeric)
stopifnot(max(abs(as.matrix(dat[,colnames(lavInspect(old_p3,'data'))])-lavInspect(old_p3,'data')))==0)
join <- function(x)paste(x,collapse=' + ')
measurement <- function(a,b=NULL,extra='') {
 paste(c(paste0(a[[1]],' =~ ',join(a[[2]])),
         if(!is.null(b))paste0(b[[1]],' =~ ',join(b[[2]])),
         paste0('PSB =~ ',join(pp)),extra),collapse='\n')
}

# Mantener todos los indicadores evita comprar un AVE mayor a costa de cobertura.
# M1 prueba colapsar SC/SBL; no se presupone que sea una solución adecuada.
specs <- list(
 M0_original=list(factors=c('SC','SBL'),measurement=measurement(list('SC',sc),list('SBL',sbl)),covariance='SC ~~ SBL',description='Modelo publicado de tres factores y SEM P3 de referencia'),
 M1_ASP_unico=list(factors='SA',measurement=measurement(list('SA',aps)),covariance='',description='Una actitud sostenible común y PSB separado'),
 M2_puentes2=list(factors=c('SC','SBL'),measurement=measurement(list('SC',sc),list('SBL',sbl),'SBL =~ APS2\nSC =~ APS13'),covariance='SC ~~ SBL',description='Dos cargas cruzadas justificadas por contenidos empresa sociedad y comunidad'),
 M3_puentes4=list(factors=c('SC','SBL'),measurement=measurement(list('SC',sc),list('SBL',sbl),'SBL =~ APS2 + APS14\nSC =~ APS13 + APS18'),covariance='SC ~~ SBL',description='Sensibilidad con cuatro ítems puente fijados por contenido'),
 M4_S1_refSC=list(factors=c('SCref','SBLspec'),measurement=measurement(list('SCref',aps),list('SBLspec',sbl),'SCref ~~ 0*SBLspec'),covariance='SCref ~~ 0*SBLspec',description='S1: orientación definida por SC y componente residual de SBL'),
 M5_S1_refSBL=list(factors=c('SBLref','SCspec'),measurement=measurement(list('SBLref',aps),list('SCspec',sc),'SBLref ~~ 0*SCspec'),covariance='SBLref ~~ 0*SCspec',description='Sensibilidad S1 con dominio de referencia invertido'),
 M6_ESEM_target=list(factors=c('SC','SBL'),measurement=paste0('efa("asp")*SC + efa("asp")*SBL =~ ',join(aps),'\nPSB =~ ',join(pp)),covariance='SC ~~ SBL',description='Dos factores ASP con rotación target oblicua y PSB CFA separado',esem=TRUE)
)

# Repetir la misma parte estructural permite atribuir diferencias a la medición.
# S1 impone ortogonalidad de perturbaciones después de ajustar predictores,
# no necesariamente ortogonalidad marginal entre factores.
structural <- function(factors,covariance) {
 lines <- vapply(seq_along(values),function(j)paste0(values[j],' ~ a',j,'*CSR1_bin + ',join(controls)),character(1))
 for(k in seq_along(factors)) {
   rhs <- c(paste0('b',k,'_',seq_along(values),'*',values),paste0('q',k,'*CSR1_bin'),controls)
   lines <- c(lines,paste0(factors[k],' ~ ',join(rhs)))
 }
 lines <- c(lines,paste0('PSB ~ ',join(c(paste0('d',seq_along(factors),'*',factors),paste0('p',seq_along(values),'*',values),'qpsb*CSR1_bin',controls))),
            'Benevolencia ~~ Universalismo + Poder + Logro','Universalismo ~~ Poder + Logro','Poder ~~ Logro',covariance)
 for(j in seq_along(values)) {
   products <- paste0('b',seq_along(factors),'_',j,'*d',seq_along(factors))
   lines <- c(lines,paste0('ind_',j,'_',seq_along(factors),' := ',products),paste0('ind_',j,'_total := ',paste(products,collapse=' + ')))
 }
 via_values <- paste(vapply(seq_along(values),function(j)paste0('a',j,'*(p',j,' + ',paste0('b',seq_along(factors),'_',j,'*d',seq_along(factors),collapse=' + '),')'),character(1)),collapse=' + ')
 lines <- c(lines,paste0('CSR_via_values := ',via_values),paste0('CSR_via_att',seq_along(factors),' := q',seq_along(factors),'*d',seq_along(factors)),paste0('CSR_total_indirect := ',via_values,' + ',paste0('q',seq_along(factors),'*d',seq_along(factors),collapse=' + ')))
 paste(lines[nzchar(lines)],collapse='\n')
}
target <- matrix(0,nrow=length(aps),ncol=2,dimnames=list(aps,c('SC','SBL')))
target[sc,'SC'] <- 1;target[sbl,'SBL'] <- 1
for(id in names(specs)) {
 specs[[id]]$id <- id
 sem_measurement <- specs[[id]]$measurement
 if(grepl('_S1_',id))sem_measurement <- sub(paste0('\n',specs[[id]]$covariance),'',sem_measurement,fixed=TRUE)
 specs[[id]]$syntax <- paste(sem_measurement,structural(specs[[id]]$factors,specs[[id]]$covariance),sep='\n')
 writeLines(specs[[id]]$syntax,file.path(out,paste0(id,'_syntax.txt')))
}
protocol <- list(date='2026-10-08',status='Frozen before new candidate estimates, after seeing original study results',N=653,indicators_retained=c(aps,pp),candidate_ids=names(specs),descriptions=lapply(specs,`[[`,'description'),crossloadings_by_content=c('APS2 on SBL','APS13 on SC','APS14 on SBL in four-bridge sensitivity','APS18 on SC in four-bridge sensitivity'),no_automatic_MI_search=TRUE,no_item_deletion=TRUE,no_constraints_to_force_positive_residuals=TRUE,comparison='MLR same observations, plus ordinal WLSMV sensitivity, admissibility, interpretability, internal repeated cross-validation',selection='No selection based on p values or AVE threshold; ESEM/S1 AVE not directly interchangeable with CFA',validation='Internal only; full dataset was previously analysed. No independent or external validation.',data_sha256=digest(data_path,algo='sha256',file=TRUE),p3_sha256=digest(p3_path,algo='sha256',file=TRUE),script_sha256=digest(file.path(root,'codigo/comparar_modelos_SEM.R'),algo='sha256',file=TRUE))
write_json(protocol,file.path(out,'protocol_before_new_fits.json'),pretty=TRUE,auto_unbox=TRUE)

safe <- function(expr) {
 warns <- character()
 fit <- tryCatch(withCallingHandlers(expr,warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
 list(fit=fit,warnings=unique(warns),error=if(inherits(fit,'error'))conditionMessage(fit)else'')
}
fit_fun <- function(spec,training,estimator='MLR',fast=FALSE,kind='SEM',group=NULL) {
 args <- list(model=if(kind=='CFA')spec$measurement else spec$syntax,data=training,estimator=estimator,meanstructure=TRUE,fixed.x=TRUE,control=list(iter.max=2000))
 if(kind=='CFA') {
   args$std.lv <- TRUE
   # Mantener FIML en CFA reproduce los índices robustos de la revisión anterior.
   # No hay faltantes; cambia el cálculo de los índices robustos, no los casos ni las cargas.
   if(estimator%in%c('MLR','ML'))args$missing <- 'fiml'
 }
 if(fast) {args$se <- 'none';args$test <- 'none'}
 if(estimator=='WLSMV') {args$ordered <- c(aps,pp);args$parameterization <- 'theta'}
 # Esta interfaz explícita es compatible con lavaan 0.7.2 instalada.
 if(isTRUE(spec$esem)) {
   args$rotation <- 'target'
   args$rotation.args <- list(target=target,rstarts=20,orthogonal=FALSE)
 }
 if(!is.null(group))args$group <- group
 z <- safe(do.call(lavaan::sem,args))
 if(inherits(z$fit,'lavaan'))attr(z$fit,'captured_warnings') <- z$warnings
 z$fit
}
isfit <- function(x)inherits(x,'lavaan')
fit_stats <- function(fit,id,kind,estimator) {
 if(!isfit(fit))return(data.frame(model=id,kind=kind,estimator=estimator,error=conditionMessage(fit)))
 fm <- tryCatch(fitMeasures(fit),error=function(e)numeric())
 # Guardar por separado robust y scaled evita sustituir etiquetas cuando
 # lavaan no calcula índices robustos para el SEM de respuestas mixtas WLSMV.
 measures <- c('chisq.scaled','df','cfi.robust','tli.robust','rmsea.robust','rmsea.ci.lower.robust','rmsea.ci.upper.robust','cfi.scaled','tli.scaled','rmsea.scaled','rmsea.ci.lower.scaled','rmsea.ci.upper.scaled','srmr','aic','bic','npar')
 nums <- setNames(vapply(measures,function(n)if(n%in%names(fm))unname(fm[n])else NA_real_,numeric(1)),measures)
 pe <- parameterEstimates(fit,se=FALSE)
 resid <- pe$est[pe$op=='~~'&pe$lhs==pe$rhs]
 cor <- tryCatch(lavInspect(fit,'cor.lv'),error=function(e)matrix(NA_real_,1,1))
 cors <- if(is.list(cor))cor else list(cor)
 eig <- function(x)tryCatch(min(eigen(x,symmetric=TRUE,only.values=TRUE)$values),error=function(e)NA_real_)
 theta <- lavInspect(fit,'theta');thetas<-if(is.list(theta))theta else list(theta)
 matrices <- lavInspect(fit,'est');psis <- if(!is.null(matrices$psi))list(matrices$psi)else lapply(matrices,`[[`,'psi')
 cbind(data.frame(model=id,kind=kind,estimator=estimator,N=sum(lavInspect(fit,'nobs')),converged=lavInspect(fit,'converged'),post_check=lavInspect(fit,'post.check'),min_residual=min(resid),min_theta_eigen=min(vapply(thetas,eig,numeric(1))),min_psi_eigen=min(vapply(psis,eig,numeric(1))),max_abs_cor=max(vapply(cors,function(m)max(abs(m[lower.tri(m)])),numeric(1)),na.rm=TRUE),error='',warnings=paste(attr(fit,'captured_warnings'),collapse=' | ')),as.data.frame(as.list(nums)))
}
bind <- function(rows) {
 cols <- unique(unlist(lapply(rows,names)))
 do.call(rbind,lapply(rows,function(r){r[setdiff(cols,names(r))]<-NA;r[,cols,drop=FALSE]}))
}
write <- function(x,name)write.csv(x,file.path(out,name),row.names=FALSE,na='')
fits <- list();fitrows <- list();reliability <- list();paths <- list();loads <- list()
for(id in names(specs)) {
 spec <- specs[[id]]
 for(kind in c('CFA','SEM'))for(estimator in c('MLR','WLSMV')) {
  key <- paste(id,kind,estimator,sep='_')
  cat('FITTING',key,'\n');flush.console()
  f <- fit_fun(spec,dat,estimator=estimator,kind=kind)
  fits[[key]] <- f
  fitrows[[key]] <- fit_stats(f,id,kind,estimator)
  saveRDS(f,file.path(out,paste0(key,'.rds')))
  write(bind(fitrows),'model_fit.csv')
  if(!isfit(f)||!lavInspect(f,'converged'))next
  ss <- standardizedSolution(f,type='std.all',ci=TRUE)
  write(parameterEstimates(f,ci=TRUE),paste0(key,'_parameters.csv'))
  write(ss,paste0(key,'_standardized.csv'))
  lf <- ss[ss$op=='=~',];lf$model<-id;lf$kind<-kind;lf$estimator<-estimator;loads[[key]]<-lf
  if(kind=='CFA') {
   av <- safe(semTools::AVE(f))
   if(is.numeric(av$fit))reliability[[paste0(key,'AVE')]]<-data.frame(model=id,estimator=estimator,measure='AVE',factor=names(av$fit),value=as.numeric(av$fit),note=paste(if(id%in%c('M0_original','M1_ASP_unico'))'Simple structure AVE based on observed variances'else'Crossloaded factors return NA; never interpreted as repaired validity',if(estimator=='WLSMV')'Ordinal latent response scale'else'Continuous working response scale'))
   # Reliability with cross-loadings is evaluated explicitly using the model matrices below.
   mats <- lavInspect(f,'est');L <- mats$lambda;V <- mats$psi;T <- mats$theta
   observed <- L%*%V%*%t(L)+T
   for(fac in spec$factors) {
    items <- if(grepl('spec$',fac))if(fac=='SBLspec')sbl else sc else if(fac%in%c('SA','SCref','SBLref'))aps else if(fac=='SC')sc else sbl
    idx <- match(items,rownames(L));k <- match(fac,colnames(L))
    # With correlated factors this is a component ratio, NOT a causal or uniquely reliable share.
    numerator <- sum(L[idx,k])^2*V[k,k]
    den <- sum(observed[idx,idx])
    reliability[[paste(key,fac,'component',sep='_')]] <- data.frame(model=id,estimator=estimator,measure='component_variance_ratio',factor=fac,value=numerator/den,note=if(grepl('spec$',fac))'Omega hierarchical subscale for orthogonal S1 CFA; domain residual only'else'Factor component ratio; correlated-factor shared covariance prevents unique interpretation')
   }
  }
  if(kind=='SEM') {
   foc <- ss[ss$op=='~' & ss$lhs%in%c(values,spec$factors,'PSB') & ss$rhs%in%c(values,spec$factors,'CSR1_bin'),]
   pe <- parameterEstimates(f,ci=TRUE)
   un <- pe[pe$op=='~' & pe$lhs%in%c(values,spec$factors,'PSB') & pe$rhs%in%c(values,spec$factors,'CSR1_bin'),]
   foc$p_holm_B <- p.adjust(un$pvalue[match(paste(foc$lhs,foc$rhs),paste(un$lhs,un$rhs))],method='holm')
   foc$family_size <- nrow(foc);foc$model<-id;foc$estimator<-estimator;paths[[key]]<-foc
   ind <- ss[ss$op==':=',];ind$p_holm<-p.adjust(ind$pvalue,'holm');ind$family_size<-nrow(ind)
   write(ind,paste0(key,'_indirects.csv'))
  }
 }
}
write(bind(reliability),'reliability_and_AVE.csv')
write(bind(paths),'focal_paths_all_models.csv')
write(bind(loads),'loadings_all_models.csv')
# Reproduction of P3 is a numerical provenance check, not a statistical test.
new <- standardizedSolution(fits$M0_original_SEM_MLR)
old <- standardizedSolution(old_p3)
new <- new[new$op=='~',];old <- old[old$op=='~',]
key <- function(x)paste(x$lhs,x$rhs)
delta <- max(abs(new$est.std-old$est.std[match(key(new),key(old))]))
stopifnot(delta<1e-5)
write(data.frame(check='P3 reproduced all regression betas',max_difference=delta),'baseline_reproduction.csv')

# Configural comparisons remain a measurement gate, not a test of national differences.
for(id in c('M0_original','M2_puentes2','M3_puentes4','M4_S1_refSC','M5_S1_refSBL')) {
 cat('CONFIGURAL',id,'\n');flush.console()
 f <- fit_fun(specs[[id]],dat,kind='CFA',group='Pais')
 fits[[paste0(id,'_configural')]]<-f
 saveRDS(f,file.path(out,paste0(id,'_configural.rds')))
 fitrows[[paste0(id,'_configural')]]<-fit_stats(f,id,'configural','MLR')
 write(bind(fitrows),'model_fit.csv')
}
saveRDS(fits,file.path(out,'all_fits_private.rds'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
source(file.path(root,'codigo/cross_validate_models.R'))
cat('CROSS VALIDATION 5 folds x 5 repetitions\n');flush.console()
run_sem_cv(dat,specs,fit_fun,file.path(out,'cross_validation'),folds=5,repeats=5,seed=20261008)
stopifnot(identical(digest(data_path,algo='sha256',file=TRUE),protocol$data_sha256),identical(digest(p3_path,algo='sha256',file=TRUE),protocol$p3_sha256))
writeLines('COMPLETED',state)
cat('COMPLETED',out,'\n')
