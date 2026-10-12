# Complete M2 evidence: extraction from saved fits and aggregate descriptions.
# This script NEVER calls cfa(), sem(), lavaan(), lm(), or any fitting routine.
# Usage: Rscript --vanilla 01_extract_complete_M2.R /path/to/project /path/to/output [old_S7_S8_csv]
# Optional R library: Sys.setenv(M2_R_LIBRARY='/path/to/existing/R-library').
# Private inputs stay outside the public package; no participant-level output.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2L) stop('Supply private project directory and output directory.')
project <- normalizePath(args[1L],mustWork=TRUE)
out <- normalizePath(args[2L],mustWork=FALSE)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
lib <- Sys.getenv('M2_R_LIBRARY','')
if(nzchar(lib)) .libPaths(c(lib,.libPaths()))
if(dir.exists(file.path(project,'revision-revista-2026-10-04/auditoria-codigos/R-library')))
  .libPaths(c(file.path(project,'revision-revista-2026-10-04/auditoria-codigos/R-library'),.libPaths()))
needed <- c('lavaan','haven','jsonlite','digest')
stopifnot(all(vapply(needed,requireNamespace,logical(1),quietly=TRUE)))
options(stringsAsFactors=FALSE,warn=1)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
jsonlite::write_json(list(operation='existing-fit extraction and aggregate descriptions; no refitting',
  RNGkind=RNGkind(),seed='None set: no stochastic operation',
  packages=as.list(setNames(vapply(needed,function(x)as.character(packageVersion(x)),character(1)),needed))),
  file.path(out,'environment.json'),pretty=TRUE,auto_unbox=TRUE)
csv <- function(x,name) utils::write.csv(x,file.path(out,paste0(name,'.csv')),row.names=FALSE,na='NA',fileEncoding='UTF-8')
fmt <- function(x) ifelse(is.na(x),'NA',formatC(x,format='f',digits=3))
fp <- function(x) ifelse(is.na(x),'—',ifelse(x<.001,'<.001',fmt(x)))
cell <- function(x,lo,hi) paste0(fmt(x),' [',fmt(lo),', ',fmt(hi),']')
manifest <- list()
add_table <- function(id,caption,notes,data,file,model,estimator,source,n=653L){
  csv(data,file)
  manifest[[length(manifest)+1L]] <<- list(id=id,caption=caption,notes=notes,
    note=notes,path=paste0(file,'.csv'),formatted_csv=paste0(file,'.csv'),
    full_csv=if(id=='T3')'T3_S8_item_paths_complete.csv' else
      if(id=='S22')'S22_descriptives_numeric.csv' else
      if(id%in%c('S23','S24','S25','S26'))paste0(id,'_loadings_complete.csv') else
      if(id=='S27')'S27_correlations_complete.csv' else
      if(startsWith(id,'S28'))paste0('S28_',if(grepl('A',id))'MLR' else 'WLSMV','_all_regressions.csv') else
      paste0('S29_',if(grepl('A',id))'MLR' else 'WLSMV','_indirects_complete.csv'),
    landscape=ncol(data)>4L,rows=nrow(data),columns=ncol(data),model=model,
    estimator=estimator,N=n,source=source,script='01_extract_complete_M2.R')
}
sources <- c(raw='datos_originales/Base_Final.sav',
  analytic='revision-integral-2026-10-04/analysis/outputs/canonical-full-20261004-02/01_audit/analytic_data.rds',
  item_MLR='sem-items-ordinal-2026-10-08/items/privados/SEM_items_MLR.rds',
  item_WLSMV='sem-items-ordinal-2026-10-08/items/privados/SEM_items_WLSMV.rds')
for(kind in c('CFA','SEM')) for(est in c('MLR','WLSMV'))
  sources[paste(kind,est,sep='_')] <- paste0('sem-alternativas-2026-10-08/resultados/M2_puentes2_',kind,'_',est,'.rds')
stopifnot(all(file.exists(file.path(project,sources))))
hashes <- vapply(file.path(project,sources),digest::digest,character(1),algo='sha256',file=TRUE)
csv(data.frame(input=names(sources),relative_path=unname(sources),sha256=hashes),'input_manifest')
fits <- lapply(sources[grepl('^(item|CFA|SEM)_',names(sources))],function(p)readRDS(file.path(project,p)))
status <- do.call(rbind,lapply(names(fits),function(id){f<-fits[[id]];op<-lavaan::lavInspect(f,'options');
  data.frame(model=id,N=sum(lavaan::lavInspect(f,'nobs')),converged=all(lavaan::lavInspect(f,'converged')),
    admissible=all(lavaan::lavInspect(f,'post.check')),estimator=op$estimator,se=op$se,
    conditional_x=op$conditional.x,fixed_x=op$fixed.x,parameterization=op$parameterization)}))
csv(status,'model_status');stopifnot(all(status$N==653),all(status$converged),all(status$admissible))

# Join b with its own b confidence interval and beta with its own delta-method CI.
# Key includes group and operator; a key mismatch is a hard failure.
extract <- function(f){
  p <- lavaan::parameterEstimates(f,ci=TRUE)
  s <- lavaan::standardizedSolution(f,type='std.all',ci=TRUE,cov.std=TRUE)
  key <- function(x)paste(x$lhs,x$op,x$rhs,if('group'%in%names(x))x$group else 1L,sep='|')
  stopifnot(!anyDuplicated(key(p)),!anyDuplicated(key(s)))
  j <- match(key(p),key(s));stopifnot(!anyNA(j))
  data.frame(lhs=p$lhs,op=p$op,rhs=p$rhs,label=if('label'%in%names(p))p$label else '',b=p$est,b_se=p$se,
    b_ci_lower=p$ci.lower,b_ci_upper=p$ci.upper,p_b=p$pvalue,
    beta=s$est.std[j],beta_se=s$se[j],beta_ci_lower=s$ci.lower[j],beta_ci_upper=s$ci.upper[j],
    p_beta=s$pvalue[j])
}
params <- lapply(fits,extract)
for(id in names(params))csv(params[[id]],paste0(id,'_all_parameters'))
controls <- c('C_Colombia','C_Peru','C_Mexico','Sexo_num','Estudios_num','Edad_z','NSE_z')
values <- c('Benevolencia','Universalismo','Poder','Logro')
focus <- c('SC','SBL','Universalismo','Benevolencia','Poder','Logro')
english <- c(SC='SC',SBL='SBL',PSB='PSB',Universalismo='Universalism',Benevolencia='Benevolence',Poder='Power',Logro='Achievement',CSR1_bin='CSR1 (yes=1)')
en <- function(x)ifelse(x%in%names(english),english[x],x)

# NEW MAIN TABLE 3: exactly the existing 18-path MLR family used in old S7.
# Reordering rows does not change the family: 6 predictors x 3 outcomes.
item_all <- list()
for(est in c('MLR','WLSMV')){
  z <- params[[paste0('item_',est)]]
  z <- z[z$op=='~' & z$lhs%in%paste0('PP',1:3) & z$rhs%in%focus,]
  stopifnot(nrow(z)==18L)
  z$p_Holm_18 <- p.adjust(z$p_b,'holm');z$estimator<-est
  item_all[[est]]<-z
  mat <- data.frame(Predictor=en(focus),check.names=FALSE)
  for(pp in paste0('PP',1:3)){
    a<-z[match(paste(pp,focus),paste(z$lhs,z$rhs)),];stopifnot(!anyNA(a$lhs))
    outcome_labels<-c(PP1='helping',PP2='sharing',PP3='empathic support')
    mat[[paste(pp,outcome_labels[[pp]],'β [95% CI]')]]<-paste0(fmt(a$beta),ifelse(a$p_Holm_18<.05,'*',''),' [',fmt(a$beta_ci_lower),', ',fmt(a$beta_ci_upper),']')
  }
  if(est=='MLR')add_table('T3','Adjusted associations with the three observed prosocial self-descriptions (M2 item-response SEM, MLR).',
    'N=653. PP1=helping; PP2=sharing; PP3=empathic support. Cells show fully standardized beta and its own pointwise 95% CI. * p_Holm<.05 across all 18 unstandardized focal path tests. All six predictors are shown. Models include CSR1 and the seven controls. This model has three observed outcomes and is distinct from the latent-PSB SEM. No new model was estimated.',mat,'T3_items_MLR_formatted','M2 observed PP1-PP3','MLR',sources[['item_MLR']])
  else csv(mat,'S8_items_WLSMV_formatted_check')
}
csv(do.call(rbind,item_all),'T3_S8_item_paths_complete')
# Cotejo with the previously distributed 18-path family, including both scales.
old <- if(length(args)>=3L)normalizePath(args[3L],mustWork=TRUE)else ''
check <- data.frame(check='previous_S7_S8_numeric_export_available',passed=if(nzchar(old))TRUE else NA,
  detail=if(nzchar(old))old else 'NOT_CHECKED: optional reference CSV was not supplied')
if(nzchar(old)){
  o<-read.csv(old);z<-do.call(rbind,item_all);j<-match(paste(z$estimator,z$lhs,z$rhs),paste(o$estimator,o$lhs,o$rhs));stopifnot(!anyNA(j))
  for(pair in list(c('b','est'),c('beta','beta_std_all'),c('beta_ci_lower','beta_CI_lower'),c('beta_ci_upper','beta_CI_upper'),c('p_Holm_18','p_Holm_18'))){
    delta<-max(abs(z[[pair[1]]]-o[[pair[2]]][j]));check<-rbind(check,data.frame(check=paste('S7/S8',pair[1]),passed=delta<1e-12,detail=format(delta,digits=16)))
  }
}
csv(check,'S7_S8_comparison');stopifnot(all(check$passed,na.rm=TRUE))

# S22: aggregate sample description from SPSS; audit against analytic rows.
# user_na=TRUE preserves SPSS missing definitions before zapping them explicitly.
raw_labeled <- haven::read_sav(file.path(project,sources[['raw']]),user_na=TRUE)
dictionary <- do.call(rbind,lapply(seq_along(raw_labeled),function(j){x<-raw_labeled[[j]];lab<-attr(x,'labels');
  data.frame(original_name=names(raw_labeled)[j],variable_label=if(is.null(attr(x,'label')))'' else attr(x,'label'),
    value_labels=if(is.null(lab))'' else paste(paste(names(lab),unname(lab),sep='='),collapse='; '),
    user_missing_values=paste(attr(x,'na_values'),collapse='; '),user_missing_range=paste(attr(x,'na_range'),collapse='; '))}))
csv(dictionary,'dictionary_original_SPSs_labels')
raw<-as.data.frame(lapply(raw_labeled,haven::zap_missing));names(raw)<-names(raw_labeled)
dat<-as.data.frame(readRDS(file.path(project,sources[['analytic']])))
stopifnot(nrow(raw)==653L,nrow(dat)==653L,all(as.numeric(raw[['Código']])==dat$Codigo),
  all(as.numeric(raw[['CSR1']])%in%c(1,2)),all(as.numeric(raw[['CSR1_rec']])%in%c(0,1)),
  all(as.numeric(raw[['CSR1_rec']])==as.numeric(as.numeric(raw[['CSR1']])==1)),
  all(dat$CSR1_bin==as.numeric(raw[['CSR1_rec']])) )
# Last comparison explicitly checks code 1 (Yes) -> 1; code 2 (No) -> 0.
country_names<-c('Argentina','Colombia','Peru','Mexico')
country<-country_names[as.numeric(raw[['País']])]
stopifnot(all(country==as.character(dat$Pais)))
for(v in c('Sexo','Estudios','Edad','NSE'))stopifnot(all(as.numeric(raw[[v]])==as.numeric(dat[[v]])))
stopifnot(all(dat$Sexo_num==as.numeric(raw$Sexo)),all(dat$Estudios_num==as.numeric(raw$Estudios)),
  max(abs(dat$Edad_z-as.numeric(scale(raw$Edad))))<1e-12,
  max(abs(dat$NSE_z-as.numeric(scale(raw$NSE))))<1e-12)
expected<-dat
for(ct in c('Colombia','Peru','Mexico'))expected[[paste0('C_',ct)]]<-as.numeric(country==ct)
data_checks<-do.call(rbind,lapply(names(fits),function(id){
  observed<-lavaan::lavInspect(fits[[id]],'data');stopifnot(nrow(observed)==653L,all(colnames(observed)%in%names(expected)))
  reference<-as.matrix(expected[,colnames(observed),drop=FALSE]);storage.mode(reference)<-'double'
  delta<-max(abs(observed-reference));stopifnot(is.finite(delta),delta<1e-12)
  data.frame(model=id,N=nrow(observed),fields=ncol(observed),max_absolute_difference=delta,
    comparison='Stored fit data vs row-aligned analytic data and original-label country dummy coding')
}))
csv(data_checks,'fit_data_source_checks')
csv(data.frame(variable=c('CSR1','CSR1_rec','CSR1_bin','Sexo_num','Estudios_num','Edad_z','NSE_z','C_Colombia','C_Peru','C_Mexico'),
  definition=c('SPSS original: Yes=1; No=2','SPSS recode: Yes=1; No=0','Analytic: identical to CSR1_rec; Yes=1; No=0',
  'Woman=1; Man=0','Postgraduate=1; Undergraduate=0','Global sample mean/SD z score for age',
  'Global sample mean/SD z score for 1-5 socioeconomic code','Institution in Colombia=1; Argentina reference',
  'Institution in Peru=1; Argentina reference','Institution in Mexico=1; Argentina reference')),'control_dictionary')
groups<-c('Overall',country_names)
rows<-list()
add_desc<-function(variable,category,values) {rows[[length(rows)+1L]]<<-data.frame(Variable=variable,Category=category,as.list(setNames(values,c('Overall',paste0('Institution in ',country_names)))),check.names=FALSE)}
subset_group<-function(g)if(g=='Overall')rep(TRUE,nrow(raw))else country==g
add_desc('Sample','N',vapply(groups,function(g)as.character(sum(subset_group(g))),character(1)))
age<-as.numeric(raw$Edad)
add_desc('Age','Mean (SD)',vapply(groups,function(g){x<-age[subset_group(g)];sprintf('%.2f (%.2f)',mean(x,na.rm=TRUE),sd(x,na.rm=TRUE))},character(1)))
add_desc('Age','Observed / missing',vapply(groups,function(g){x<-age[subset_group(g)];paste0(sum(!is.na(x)),' / ',sum(is.na(x)))},character(1)))
cats<-list(Sex=list(v='Sexo',codes=c(0,1),labs=c('Man','Woman')),
  Education=list(v='Estudios',codes=c(0,1),labs=c('Undergraduate','Postgraduate')),
  `Socioeconomic level`=list(v='NSE',codes=1:5,labs=c('Low','Lower-middle','Middle','Upper-middle','High')),
  `CSR1: attended ethics/CSR course at university`=list(v='CSR1_rec',codes=c(0,1),labs=c('No','Yes')))
desc_numeric<-list()
for(v in names(cats)){s<-cats[[v]];x<-as.numeric(raw[[s$v]])
  for(k in seq_along(s$codes)){
    vals<-vapply(groups,function(g){ix<-subset_group(g);den<-sum(!is.na(x[ix]));n<-sum(x[ix]==s$codes[k],na.rm=TRUE);
      desc_numeric[[length(desc_numeric)+1L]]<<-data.frame(variable=v,category=s$labs[k],group=g,n=n,denominator=den,missing=sum(is.na(x[ix])),percent=100*n/den)
      sprintf('%d (%.1f%%)',n,100*n/den)},character(1));add_desc(v,s$labs[k],vals)
  }
  add_desc(v,'Denominator / missing',vapply(groups,function(g){a<-x[subset_group(g)];paste0(sum(!is.na(a)),' / ',sum(is.na(a)))},character(1)))
}
add_table('S22','Sample characteristics overall and by participating institution.',
  'N=653. Values are n (%) except age, mean (SD). Percent denominators exclude missing values within each institution; denominator/missing rows are explicit. There is one institution per country, so institution and country effects cannot be separated. Institution names were not verified in the analytic file; neutral geographic labels are used. CSR1 asks whether a course on ethics and/or corporate social responsibility was taken at the university. The analytic source contains all 653 rows but does not establish that no selection occurred before the retained file.',
  do.call(rbind,rows),'S22_descriptives_formatted','Retained sample','Descriptive',sources[['raw']])
csv(do.call(rbind,desc_numeric),'S22_descriptives_numeric')
csv(do.call(rbind,lapply(groups,function(g){x<-age[subset_group(g)];data.frame(group=g,n=sum(!is.na(x)),missing=sum(is.na(x)),mean=mean(x,na.rm=TRUE),sd=sd(x,na.rm=TRUE),min=min(x,na.rm=TRUE),max=max(x,na.rm=TRUE))})),'S22_age_numeric')
csv(data.frame(check=c('raw_n','analytic_n','participant_order','CSR1_semantic_recode','sex_education_age_NSE_match','global_z_scores'),
  result=c('653','653','Exact match','Original Yes=1/No=2 -> analytic Yes=1/No=0','Exact match','Exact within 1e-12')),'S22_audit_checks')

# S23-S26: every M2 loading. No mixing of estimators or CFA and SEM.
loadtables<-c(CFA_MLR='S23',CFA_WLSMV='S24',SEM_MLR='S25',SEM_WLSMV='S26')
for(id in names(loadtables)){
  z<-params[[id]];z<-z[z$op=='=~',];stopifnot(nrow(z)==24L)
  pt<-lavaan::parTable(fits[[id]])
  idx<-match(paste(z$lhs,z$op,z$rhs),paste(pt$lhs,pt$op,pt$rhs))
  stopifnot(!anyNA(idx));z$b_fixed<-pt$free[idx]==0L
  z$role<-ifelse((z$lhs=='SBL'&z$rhs=='APS2')|(z$lhs=='SC'&z$rhs=='APS13'),'Cross-loading','Primary')
  csv(z,paste0(loadtables[id],'_loadings_complete'))
  m<-data.frame(Factor=z$lhs,Item=z$rhs,Role=z$role,
    `b [95% CI]`=ifelse(z$b_fixed,paste0(fmt(z$b),' (fixed)'),cell(z$b,z$b_ci_lower,z$b_ci_upper)),
    `beta [95% CI]`=cell(z$beta,z$beta_ci_lower,z$beta_ci_upper),`p (b)`=fp(z$p_b),check.names=FALSE)
  add_table(loadtables[id],paste('Complete M2 loadings:',gsub('_','; ',id),'.'),
    paste('N=653. b and beta (std.all) each have their own pointwise 95% interval; intervals are not multiplicity-adjusted. A fixed unstandardized marker loading is marked fixed, not presented with a confidence interval; its standardized loading remains an estimated quantity and retains its own interval. Two prespecified cross-loadings are APS2 on SBL and APS13 on SC.',
      if(grepl('WLSMV',id))'WLSMV is DWLS with robust mean/variance-adjusted tests and theta parameterization; loadings concern underlying ordinal responses.' else 'MLR treats the indicators as numeric scores and uses robust inference.'),
    m,paste0(loadtables[id],'_loadings_formatted'),id,sub('.*_','',id),sources[[id]])
}

# S27: total factor correlations and disturbance correlations are separate panels.
# lavInspect(cor.lv) is the model-implied total factor correlation (including
# effects propagated through the structural regressions and observed X covariance).
# The standardized ~~ parameters with cov.std=TRUE describe residual correlation.
# CFA intervals are extracted from the matching standardized ~~ parameter.
# For SEM total correlations cor.lv provides point estimates only; no CI is made up.
corrows<-list()
for(id in names(loadtables)){
 f<-fits[[id]];m<-lavaan::lavInspect(f,'est');lv<-c('SC','SBL','PSB')
 total<-lavaan::lavInspect(f,'cor.lv');psi<-m$psi[lv,lv,drop=FALSE];res<-cov2cor(psi)
 # Independent algebraic check of total covariance: include observed X covariance
 # for conditional.x fits. This only transforms saved fitted matrices.
 A<-if(is.null(m$beta))diag(nrow(m$psi))else solve(diag(nrow(m$beta))-m$beta)
 total_cov<-A%*%(m$psi+if(!is.null(m$gamma))m$gamma%*%m$cov.x%*%t(m$gamma)else 0)%*%t(A)
 dimnames(total_cov)<-dimnames(m$psi)
 stopifnot(max(abs(cov2cor(total_cov)[lv,lv]-total[lv,lv]))<1e-9)
 for(pair in list(c('SC','SBL'),c('SC','PSB'),c('SBL','PSB'))){
   total_lower<-total_upper<-NA_real_;total_status<-'Point estimate; SEM total-correlation CI not estimated'
   if(grepl('CFA_',id)){
     z<-params[[id]];a<-z[z$op=='~~'&((z$lhs==pair[1]&z$rhs==pair[2])|(z$lhs==pair[2]&z$rhs==pair[1])),]
     stopifnot(nrow(a)==1L,abs(a$beta-total[pair[1],pair[2]])<1e-9)
     total_lower<-a$beta_ci_lower;total_upper<-a$beta_ci_upper;total_status<-'Pointwise 95% CI from standardized CFA covariance'
   }
   corrows[[length(corrows)+1L]]<-data.frame(Model=id,Pair=paste(pair,collapse=' / '),Type='Total model-implied',r=total[pair[1],pair[2]],lower=total_lower,upper=total_upper,status=total_status)
   if(grepl('SEM_',id)){
     z<-params[[id]];a<-z[z$op=='~~'&((z$lhs==pair[1]&z$rhs==pair[2])|(z$lhs==pair[2]&z$rhs==pair[1])),]
     if(nrow(a)==1L){lo<-a$beta_ci_lower;hi<-a$beta_ci_upper;st<-if(a$b_se==0)'Fixed to zero' else 'Pointwise 95% CI'}else{lo<-hi<-NA_real_;st<-'Fixed to zero by omitted disturbance covariance'}
     corrows[[length(corrows)+1L]]<-data.frame(Model=id,Pair=paste(pair,collapse=' / '),Type='Disturbance (conditional)',r=res[pair[1],pair[2]],lower=lo,upper=hi,status=st)
   }
 }
}
cort<-do.call(rbind,corrows);csv(cort,'S27_correlations_complete')
corm<-data.frame(Model=gsub('_','; ',cort$Model),Pair=cort$Pair,Type=cort$Type,
  `r [95% CI where available]`=ifelse(is.na(cort$lower),fmt(cort$r),cell(cort$r,cort$lower,cort$upper)),Status=cort$status,check.names=FALSE)
add_table('S27','Total factor correlations and conditional disturbance correlations in the four M2 fits.',
  'N=653. Total correlations are from lavInspect(cor.lv), which propagates the fitted structural regressions and observed exogenous covariance. CFA total correlations and their pointwise 95% intervals match the standardized CFA covariances. SEM total correlations are distinct from conditional disturbance correlations computed from psi; their intervals were not estimated. SEM SC-PSB and SBL-PSB disturbance covariances are fixed to zero, not empirical evidence of independence. Intervals for free residual correlations are pointwise std.all intervals.',
  corm,'S27_correlations_formatted','M2 CFA/SEM','MLR/WLSMV',paste(sources[names(loadtables)],collapse='; '))

# S28 and S29: preserve the originally specified families, separately per fit.
# Regression family is 21 focal ~ paths; 49 control paths are descriptive.
# Indirect family is all 16 := effects; preserve the ORIGINAL standardized
# p-values used by comparar_modelos_SEM.R. A separate audit column computes
# Holm from unstandardized p-values without changing reported decisions.
for(est in c('MLR','WLSMV')){
 id<-paste0('SEM_',est);z<-params[[id]];r<-z[z$op=='~',]
 focal<-r$lhs%in%c(values,'SC','SBL','PSB') & r$rhs%in%c(values,'SC','SBL','CSR1_bin')
 stopifnot(sum(focal)==21L,nrow(r)==70L)
 r$family<-ifelse(focal,'21 focal regression paths','Control path; not in focal Holm family')
 r$p_Holm<-NA_real_;r$p_Holm[focal]<-p.adjust(r$p_b[focal],'holm')
 csv(r,paste0('S28_',est,'_all_regressions'))
 for(panel in c('focal','controls')){
  a<-r[if(panel=='focal')focal else !focal,]
  m<-data.frame(Outcome=en(a$lhs),Predictor=en(a$rhs),`b [95% CI]`=cell(a$b,a$b_ci_lower,a$b_ci_upper),
    `beta [95% CI]`=cell(a$beta,a$beta_ci_lower,a$beta_ci_upper),`p (b)`=fp(a$p_b),`p Holm21`=fp(a$p_Holm),check.names=FALSE)
  add_table(paste0('S28',if(est=='MLR')'A' else 'B',if(panel=='focal')'1' else '2'),
    paste('M2 SEM',est,if(panel=='focal')'complete focal paths.' else 'complete control paths.'),
    'N=653. Unstandardized b and fully standardized beta have their own pointwise 95% CIs. Holm uses unstandardized p values across 21 focal regressions within each estimator; it does not adjust control paths, model selection or analyses across estimators. Controls: institution indicators (Argentina reference), woman=1, postgraduate=1, globally standardized age and socioeconomic level.',
    m,paste0('S28_',est,'_',panel,'_formatted'),id,est,sources[[id]])
 }
 ind<-z[z$op==':=',];stopifnot(nrow(ind)==16L)
 ind$p_Holm16_b_audit<-p.adjust(ind$p_b,'holm')
 ind$p_Holm16_originalStd<-p.adjust(ind$p_beta,'holm')
 historical<-read.csv(file.path(project,paste0('sem-alternativas-2026-10-08/resultados/M2_puentes2_SEM_',est,'_indirects.csv')))
 j<-match(ind$lhs,historical$lhs);stopifnot(!anyNA(j),max(abs(ind$p_Holm16_originalStd-historical$p_holm[j]))<1e-10)
 csv(data.frame(estimator=est,effect=ind$lhs,p_b=ind$p_b,p_beta=ind$p_beta,
  Holm16_b=ind$p_Holm16_b_audit,Holm16_standardized=ind$p_Holm16_originalStd,
  support_b=ind$p_Holm16_b_audit<.05,support_standardized=ind$p_Holm16_originalStd<.05),paste0('S29_',est,'_historical_rule_comparison'))
 ind$meaning<-vapply(ind$lhs,function(name){
   if(grepl('^ind_',name)){v<-as.integer(strsplit(name,'_')[[1]][2]);suffix<-strsplit(name,'_')[[1]][3];paste(en(values[v]),'via',switch(suffix,'1'='SC','2'='SBL','total'='SC and SBL jointly'),'to PSB')}
   else switch(name,CSR_via_values='CSR1 through value portraits (including their SC/SBL paths) to PSB',CSR_via_att1='CSR1 through SC to PSB',CSR_via_att2='CSR1 through SBL to PSB',CSR_total_indirect='Total defined indirect association of CSR1 with PSB',name)},character(1))
 csv(ind,paste0('S29_',est,'_indirects_complete'))
 im<-data.frame(Effect=ind$lhs,Definition=ind$meaning,`b [95% CI]`=cell(ind$b,ind$b_ci_lower,ind$b_ci_upper),
   `beta [95% CI]`=cell(ind$beta,ind$beta_ci_lower,ind$beta_ci_upper),`p (beta)`=fp(ind$p_beta),`p Holm16 (original)`=fp(ind$p_Holm16_originalStd),check.names=FALSE)
 add_table(paste0('S29',if(est=='MLR')'A' else 'B'),paste('All 16 defined M2 indirect associations,',est,'.'),
   'N=653. These are model-defined products/sums of cross-sectional regression coefficients, not causal mediation effects. Each b and beta interval belongs to its own scale. The ORIGINAL procedure applies Holm to the standardized-solution p values for all 16 indirect effects within each estimator; this rule is retained here. This differs from the 21-path and 18-path regression families, which use unstandardized p values. The complete CSV also provides b-based Holm as an audit-only comparison; it is not used to change decisions.',
   im,paste0('S29_',est,'_indirects_formatted'),id,est,sources[[id]])
}

# Machine-readable map for the parent document builder and public README.
jsonlite::write_json(manifest,file.path(out,'tables_manifest.json'),pretty=TRUE,auto_unbox=TRUE,na='null')
stopifnot(identical(hashes,vapply(file.path(project,sources),digest::digest,character(1),algo='sha256',file=TRUE)))
writeLines(c('PASS: all private inputs unchanged by SHA-256.',
  'PASS: 653 cases in raw, analytic and each saved fit.',
  if(nzchar(old))'PASS: new main T3 matches old S7 numeric paths and Holm18.' else 'NOT_CHECKED: optional previous S7 reference CSV not supplied.',
  'PASS: every loading/regression/indirect estimate paired with its own-scale interval.',
  'No model re-estimation, no new exclusions, no altered specification, no participant-level exports.'),file.path(out,'EXTRACTION_STATUS.txt'))
cat('SUCCESS: ',length(manifest),' document-ready table panels; saved fits extracted without re-estimation.\n',sep='')
