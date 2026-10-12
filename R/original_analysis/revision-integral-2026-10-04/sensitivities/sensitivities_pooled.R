# Reproducible pooled SEM sensitivities; originals are read only.
# Run: Rscript --vanilla sensitivities_pooled.R
# Optional environment: SENS_BOOT_B=1000, SENS_BOOT_CORES=2.
# Limit BLAS threads before R starts: VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1.
options(stringsAsFactors=FALSE, warn=1)
# Portable configuration: required SENS_INPUT, optional SENS_AUDIT_LIB / SENS_USER_LIB.
argv<-commandArgs(trailingOnly=FALSE)
script_file<-sub("^--file=","",argv[grepl("^--file=",argv)])
script_dir<-if(length(script_file))dirname(normalizePath(script_file))else getwd()
source(file.path(script_dir,"checkpoint_guard.R"))
input<-Sys.getenv("SENS_INPUT",file.path(script_dir,"Base_Final.sav"))
if(!file.exists(input))stop("Set SENS_INPUT to an existing Base_Final.sav")
input<-normalizePath(input)
audit_lib<-Sys.getenv("SENS_AUDIT_LIB","")
user_lib<-Sys.getenv("SENS_USER_LIB","")
library_paths<-c(audit_lib,user_lib);library_paths<-library_paths[nzchar(library_paths)]
if(any(!dir.exists(library_paths)))stop("Configured R library directory does not exist")
.libPaths(c(library_paths,.libPaths()))
suppressPackageStartupMessages({library(haven);library(lavaan)})
B<-as.integer(Sys.getenv("SENS_BOOT_B","1000"))
seed<-as.integer(Sys.getenv("SENS_SEED","20261004"))
cores_requested<-as.integer(Sys.getenv("SENS_BOOT_CORES","2"))
if(is.na(B)||B<1L||is.na(seed)||seed<0L||seed>=.Machine$integer.max||is.na(cores_requested)||cores_requested<1L)
 stop("Invalid B, seed or core configuration")
cores<-min(2L,cores_requested)
output_root<-Sys.getenv("SENS_OUTPUT_ROOT",file.path(script_dir,"runs"))
run_id<-Sys.getenv("SENS_RUN_ID",paste0(format(Sys.time(),"%Y%m%dT%H%M%OS6",tz="UTC"),"Z-",Sys.getpid()))
if(!grepl("^[A-Za-z0-9_.-]+$",run_id))stop("SENS_RUN_ID must contain only letters, digits, _, . or -")
out<-file.path(output_root,run_id)
if(dir.exists(out)||file.exists(out))stop("Run directory already exists; choose a fresh SENS_RUN_ID: ",out)
if(!dir.create(out,recursive=TRUE,showWarnings=FALSE))stop("Cannot create fresh run directory")
out<-normalizePath(out)
# Recovery reads a validated checkpoint from an older run and writes into this new run.
resume_checkpoint<-Sys.getenv("SENS_RESUME_CHECKPOINT","")
if(nzchar(resume_checkpoint)&&!file.exists(resume_checkpoint))stop("Requested checkpoint does not exist")
Sys.setenv(SENS_RUN_DIR=out,SENS_INPUT=input)
writeLines("STARTED",file.path(out,"run_state.txt"))
options(error=function(){writeLines("FAILED; see log. Prior runs remain unchanged.",file.path(out,"run_state.txt"));traceback(3);quit(status=1L)})
if(!all(file.copy(c(file.path(script_dir,"sensitivities_pooled.R"),file.path(script_dir,"checkpoint_guard.R"),
            file.path(script_dir,"finalize_sensitivity_tables.R"),file.path(script_dir,"export_focal_B_tables.R")),out)))
 stop("Cannot preserve run code snapshot")
set.seed(seed)
write_out <- function(x,name) write.csv(x,file.path(out,name),row.names=FALSE,na="")
raw <- read_sav(input)
names(raw) <- stringi::stri_trans_general(names(raw),"Latin-ASCII")
country_labels <- as.character(haven::as_factor(raw$Pais,levels="labels"))
country_labels <- stringi::stri_trans_general(country_labels,"Latin-ASCII")
stopifnot(all(country_labels %in% c("Argentina","Colombia","Peru","Mexico")))
raw$Pais <- factor(country_labels,levels=c("Argentina","Colombia","Peru","Mexico"))
csr_labels <- tolower(stringi::stri_trans_general(as.character(haven::as_factor(raw$CSR1_rec,levels="labels")),"Latin-ASCII"))
stopifnot(all(csr_labels %in% c("no","si")))
sc <- c("APS1","APS2","APS3","APS4","APS5","APS6","APS7","APS9","APS10","APS14")
sbl <- c("APS11","APS12","APS13","APS15","APS16","APS18","APS19","APS20","APS21")
psb <- c("PP1","PP2","PP3")
indicators <- c(sc,sbl,psb)
values <- c(ben="Benevolencia",uni="Universalismo",pow="Poder",log="Logro")
dat <- as.data.frame(raw[,c("Codigo","Pais",unname(values),indicators,"Sexo","Estudios","Edad","NSE")])
for(v in c(unname(values),indicators,"Sexo","Estudios","Edad","NSE")) dat[[v]] <- as.numeric(dat[[v]])
dat$CSR1_bin <- as.integer(csr_labels=="si")
dat$C_Colombia <- as.integer(dat$Pais=="Colombia")
dat$C_Peru <- as.integer(dat$Pais=="Peru")
dat$C_Mexico <- as.integer(dat$Pais=="Mexico")
dat$Sexo_num <- dat$Sexo
dat$Estudios_num <- dat$Estudios
dat$Edad_z <- as.numeric(scale(dat$Edad))
dat$NSE_z <- as.numeric(scale(dat$NSE))
write_out(data.frame(variable=c("Edad","NSE"),mean=c(mean(dat$Edad),mean(dat$NSE)),
                    SD=c(sd(dat$Edad),sd(dat$NSE)),reference_N=nrow(dat)),"control_centering.csv")
complete <- complete.cases(dat)
stopifnot(all(complete),nrow(dat)==653L,length(indicators)==22L)
stopifnot(all(unlist(dat[,indicators]) %in% 1:5),all(unlist(dat[,unname(values)]) %in% 1:6))
profile <- as.data.frame(raw[,setdiff(names(raw),"Codigo")])
dup_copy <- duplicated(profile)
dup_member <- dup_copy | duplicated(profile,fromLast=TRUE)
constant22 <- apply(dat[,indicators],1,function(z) length(unique(z))==1L)
quality <- data.frame(row_number=seq_len(nrow(dat)),participant_id=as.character(dat$Codigo),
                      Pais=dat$Pais,constant_22=constant22,duplicate_profile_copy=dup_copy,
                      duplicate_profile_member=dup_member)
write_out(quality,"quality_flags.csv")
scenarios <- list(all=rep(TRUE,nrow(dat)),without_duplicate_copies=!dup_copy,
                  without_constant22=!constant22,
                  without_duplicate_copies_and_constant22=!dup_copy&!constant22)
write_out(do.call(rbind,lapply(names(scenarios),function(z)
  data.frame(scenario=z,N=sum(scenarios[[z]]),excluded=sum(!scenarios[[z]])))),
  "sample_scenarios.csv")
write_out(data.frame(file=input,md5=unname(tools::md5sum(input)),N=nrow(raw),
                     columns=ncol(raw),R=R.version.string,lavaan=as.character(packageVersion("lavaan")),
                     haven=as.character(packageVersion("haven")),seed=seed,B=B,cores=cores),
          "provenance.csv")
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))

measure <- paste0("SC =~ ",paste(sc,collapse=" + "),"\nSBL =~ ",paste(sbl,collapse=" + "),
                  "\nPSB =~ ",paste(psb,collapse=" + "))
dummies <- c("C_Colombia","C_Peru","C_Mexico")
demog <- c("Sexo_num","Estudios_num","Edad_z","NSE_z")
make_model <- function(country=TRUE,partial=TRUE,demographics=FALSE) {
  country_terms <- paste(if(country) dummies else paste0("0*",dummies),collapse=" + ")
  controls <- paste0(" + ",country_terms,if(demographics) paste0(" + ",paste(demog,collapse=" + ")) else "")
  rows <- c(measure)
  for(k in names(values)) rows <- c(rows,paste0(values[[k]]," ~ a_",k,"*CSR1_bin",controls))
  sc_vals <- paste(paste0("b_",names(values),"*",unname(values)),collapse=" + ")
  sbl_vals <- paste(paste0("c_",names(values),"*",unname(values)),collapse=" + ")
  psb_vals <- paste(paste0("p_",names(values),"*",unname(values)),collapse=" + ")
  rows <- c(rows,paste0("SC ~ ",sc_vals,if(partial) " + q_sc*CSR1_bin" else " + 0*CSR1_bin",controls),
            paste0("SBL ~ ",sbl_vals,if(partial) " + q_sbl*CSR1_bin" else " + 0*CSR1_bin",controls),
            paste0("PSB ~ d_sc*SC + d_sbl*SBL + ",psb_vals,
                   if(partial) " + q_psb*CSR1_bin" else " + 0*CSR1_bin",controls),
            "Benevolencia ~~ Universalismo + Poder + Logro",
            "Universalismo ~~ Poder + Logro","Poder ~~ Logro","SC ~~ SBL")
  csr_values <- character()
  for(k in names(values)) {
    rows <- c(rows,paste0("ind_",k,"_SC := b_",k,"*d_sc"),
              paste0("ind_",k,"_SBL := c_",k,"*d_sbl"),
              paste0("ind_",k,"_total := b_",k,"*d_sc + c_",k,"*d_sbl"),
              paste0("total_",k,"_PSB := p_",k," + b_",k,"*d_sc + c_",k,"*d_sbl"))
    csr_values <- c(csr_values,paste0("a_",k,"*(p_",k," + b_",k,"*d_sc + c_",k,"*d_sbl)"))
  }
  rows <- c(rows,paste0("CSR_via_values := ",paste(csr_values,collapse=" + ")))
  if(partial) rows <- c(rows,
     "CSR_via_SC_direct := q_sc*d_sc","CSR_via_SBL_direct := q_sbl*d_sbl",
     paste0("CSR_total_indirect := ",paste(c(csr_values,"q_sc*d_sc","q_sbl*d_sbl"),collapse=" + ")),
     "CSR_total_PSB := q_psb + CSR_total_indirect")
  else rows <- c(rows,"CSR_total_indirect := CSR_via_values","CSR_total_PSB := CSR_via_values")
  paste(rows,collapse="\n")
}
# Validate recovery before fitting; every run writes into a new directory.
chosen_id<-"P3_country_partial_demographics"
chosen_syn<-make_model(country=TRUE,partial=TRUE,demographics=TRUE)
country_rows<-split(seq_len(nrow(dat)),dat$Pais)
set.seed(seed+1L)
resample_rows<-replicate(B,unlist(lapply(country_rows,function(ix)sample(ix,length(ix),replace=TRUE))),simplify=FALSE)
point_names<-c("SBL_to_PSB","SC_to_PSB",unlist(lapply(names(values),function(k)
 c(paste0("ind_",k,"_",c("SC","SBL","total")),paste0("CSR_via_",k)))),
 "CSR_via_values","CSR_via_SC_direct","CSR_via_SBL_direct","CSR_total_indirect","CSR_direct_PSB","CSR_total_PSB")
checkpoint_meta<-sens_checkpoint_meta(input,dat,chosen_syn,seed,B,resample_rows,point_names)
boot_results<-if(nzchar(resume_checkpoint))sens_read_checkpoint(resume_checkpoint,checkpoint_meta)else list()
if(length(boot_results))cat("Verified recovery:",length(boot_results),"completed replicas; output",out,"\n")
saveRDS(checkpoint_meta,file.path(out,"bootstrap_metadata.rds"))
saveRDS(resample_rows,file.path(out,"bootstrap_resample_rows_private.rds"))
write.csv(data.frame(field=names(checkpoint_meta),value=vapply(checkpoint_meta,function(z)paste(z,collapse=";"),character(1))),
 file.path(out,"bootstrap_metadata.csv"),row.names=FALSE)
write.csv(data.frame(run_id=run_id,input=input,output=out,resume_checkpoint=resume_checkpoint,
 seed=seed,B=B,cores=cores,audit_lib=audit_lib,user_lib=user_lib),file.path(out,"run_configuration.csv"),row.names=FALSE)
if(Sys.getenv("SENS_PREFLIGHT_ONLY","")=="1") {
 writeLines("PREFLIGHT_VERIFIED; no SEM estimated",file.path(out,"run_state.txt"))
 cat("Preflight metadata verified; no SEM estimated. Output:",out,"\n")
 quit(status=0L)
}
capture_fit <- function(syntax,df,estimator="MLR",fast=FALSE) {
  warns <- character()
  fit <- tryCatch(withCallingHandlers(
    sem(syntax,data=df,estimator=estimator,meanstructure=TRUE,fixed.x=TRUE,
        start=if(fast) coef(chosen) else "default",
        se=if(fast) "none" else "default",test=if(fast) "none" else "default"),
    warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart("muffleWarning")}),
    error=function(e)e)
  list(fit=if(inherits(fit,"lavaan")) fit else NULL,
       error=if(inherits(fit,"error")) conditionMessage(fit) else "",
       warnings=paste(unique(warns),collapse=" | "))
}
admissible <- function(fit) {
  if(is.null(fit)) return(data.frame(converged=FALSE,finite=FALSE,min_residual_variance=NA,
                                    min_theta_eigen=NA,min_psi_eigen=NA,max_abs_latent_correlation=NA,
                                    admissible=FALSE))
  pe <- parameterEstimates(fit,se=FALSE)
  resid <- pe$est[pe$op=="~~" & pe$lhs==pe$rhs]
  eig <- function(m) min(eigen(m,symmetric=TRUE,only.values=TRUE)$values)
  theta <- lavInspect(fit,"theta")
  psi <- lavInspect(fit,"est")$psi
  corlv <- lavInspect(fit,"cor.lv")
  maxcor <- if(nrow(corlv)>1) max(abs(corlv[lower.tri(corlv)])) else 0
  conv <- isTRUE(lavInspect(fit,"converged"))
  fin <- all(is.finite(pe$est))
  mintheta <- eig(theta);minpsi<-eig(psi)
  ok <- conv && fin && all(resid>0) && mintheta>0 && minpsi>0 && maxcor<1
  data.frame(converged=conv,finite=fin,min_residual_variance=min(resid),
             min_theta_eigen=mintheta,min_psi_eigen=minpsi,
             max_abs_latent_correlation=maxcor,admissible=ok)
}
fitrow <- function(fit,id,scenario,N) {
  measures <- c("chisq.scaled","df","pvalue.scaled","cfi.robust","tli.robust","rmsea.robust",
                "rmsea.ci.lower.robust","rmsea.ci.upper.robust","srmr","aic","bic","npar")
  fm <- if(is.null(fit)) setNames(rep(NA_real_,length(measures)),measures) else
      vapply(measures,function(m) tryCatch(unname(fitMeasures(fit,m)),error=function(e) NA_real_),numeric(1))
  cbind(data.frame(model=id,scenario=scenario,N=N),as.data.frame(as.list(fm)),admissible(fit))
}
export_fit <- function(result,id,scenario,df) {
  fit <- result$fit
  write_out(data.frame(model=id,scenario=scenario,error=result$error,warnings=result$warnings),
            paste0(id,"_warnings.csv"))
  row <- fitrow(fit,id,scenario,nrow(df))
  write_out(row,paste0(id,"_fit.csv"))
  if(is.null(fit)) return(row)
  saveRDS(fit,file.path(out,paste0(id,".rds")))
  pe <- parameterEstimates(fit,ci=TRUE,standardized=FALSE)
  ss <- standardizedSolution(fit,type="std.all",ci=TRUE)
  ssn <- standardizedSolution(fit,type="std.nox",ci=TRUE)
  write_out(pe,paste0(id,"_unstandardized_parameters_robust_ci.csv"))
  write_out(ss,paste0(id,"_standardized_parameters_robust_ci.csv"))
  write_out(ssn,paste0(id,"_standardized_nox_parameters_robust_ci.csv"))
  write_out(ss[ss$op=="~",],paste0(id,"_paths_standardized_robust_ci.csv"))
  write_out(ss[ss$op==":=",],paste0(id,"_indirects_standardized_robust_ci.csv"))
  write_out(data.frame(outcome=names(lavInspect(fit,"r2")),R2=as.numeric(lavInspect(fit,"r2"))),
            paste0(id,"_r2.csv"))
  row
}
models <- list(
 P0_no_country_full_chain=list(country=FALSE,partial=FALSE,demographics=FALSE),
 P1_country_full_chain=list(country=TRUE,partial=FALSE,demographics=FALSE),
 P2_country_partial_curricular=list(country=TRUE,partial=TRUE,demographics=FALSE),
 P3_country_partial_demographics=list(country=TRUE,partial=TRUE,demographics=TRUE),
 P4_country_full_chain_demographics=list(country=TRUE,partial=FALSE,demographics=TRUE),
 P5_no_country_partial_demographics=list(country=FALSE,partial=TRUE,demographics=TRUE))
fits <- list(); rows <- list()
for(id in names(models)) {
  cat("Fitting",id,"N=",nrow(dat),"\n")
  syn <- do.call(make_model,models[[id]])
  writeLines(syn,file.path(out,paste0(id,"_syntax.txt")))
  result <- capture_fit(syn,dat)
  fits[[id]]<-result$fit
  rows[[id]]<-export_fit(result,id,"all",dat)
}
chosen_id <- "P3_country_partial_demographics"
chosen_syn <- do.call(make_model,models[[chosen_id]])
for(s in setdiff(names(scenarios),"all")) {
  id <- paste0("S_",s)
  df <- dat[scenarios[[s]],,drop=FALSE]
  cat("Fitting",id,"N=",nrow(df),"\n")
  result <- capture_fit(chosen_syn,df)
  fits[[id]] <- result$fit
  rows[[id]] <- export_fit(result,id,s,df)
}
legacy_sbl <- c("APS11","APS12","APS13","APS16","APS18","APS19","APS21")
legacy_syn <- sub(paste0("SBL =~ ",paste(sbl,collapse=" + ")),
                  paste0("SBL =~ ",paste(legacy_sbl,collapse=" + ")),chosen_syn,fixed=TRUE)
writeLines(legacy_syn,file.path(out,"S_legacy_SBL7_syntax.txt"))
legacy_result <- capture_fit(legacy_syn,dat)
fits$S_legacy_SBL7 <- legacy_result$fit
rows$S_legacy_SBL7 <- export_fit(legacy_result,"S_legacy_SBL7","legacy_SBL7_all",dat)
summary_fit <- do.call(rbind,rows);rownames(summary_fit)<-NULL
write_out(summary_fit,"models_fit_and_admissibility.csv")
lrt <- function(a,b,name) {
  fit1<-fits[[a]];fit2<-fits[[b]]
  ans<-tryCatch(lavTestLRT(fit1,fit2,method="satorra.bentler.2001"),error=function(e)e)
  if(inherits(ans,"error")) return(data.frame(comparison=name,error=conditionMessage(ans)))
  z<-as.data.frame(ans);z$model<-rownames(z);z$comparison<-name;rownames(z)<-NULL
  z
}
write_out(rbind(lrt("P0_no_country_full_chain","P1_country_full_chain","country_fixed_zero_vs_free"),
                lrt("P1_country_full_chain","P2_country_partial_curricular","curricular_direct_fixed_zero_vs_free"),
                lrt("P4_country_full_chain_demographics","P3_country_partial_demographics","curricular_direct_fixed_zero_vs_free_same_demographics"),
                lrt("P5_no_country_partial_demographics","P3_country_partial_demographics","country_fixed_zero_vs_free_same_partial_demographics")),
          "nested_robust_MLR_comparisons.csv")
stability <- do.call(rbind,lapply(names(fits),function(id){
  f<-fits[[id]]
  if(is.null(f)) return(NULL)
  z<-standardizedSolution(f,type="std.all",ci=TRUE)
  z<-z[z$op=="~" & z$lhs %in% c("SC","SBL","PSB",unname(values)) &
             z$rhs %in% c("SC","SBL",unname(values),"CSR1_bin"),]
  z$model<-id;z$N<-rows[[id]]$N;z
}))
write_out(stability,"stability_key_paths_standardized.csv")

# Correlation-based HTMT diagnostic and bootstrap, independent of SEM fitting.
blocks <- list(SC=sc,SBL=sbl,PSB=psb)
htmt <- function(df) {
  R<-cor(df[,indicators])
  inner<-sapply(blocks,function(v)mean(abs(R[v,v][lower.tri(R[v,v])])))
  pairs<-combn(names(blocks),2,simplify=FALSE)
  setNames(sapply(pairs,function(pair)mean(abs(R[blocks[[pair[1]]],blocks[[pair[2]]]]))/
           sqrt(inner[[pair[1]]]*inner[[pair[2]]])),sapply(pairs,paste,collapse="_"))
}
htmtp<-htmt(dat)
htmtboot<-vapply(resample_rows,function(ix) htmt(dat[ix,,drop=FALSE]),htmtp)
htmtci<-t(apply(htmtboot,1,quantile,probs=c(.025,.975),na.rm=TRUE))
write_out(data.frame(pair=names(htmtp),HTMT=htmtp,CI_low=htmtci[,1],CI_high=htmtci[,2],
                     bootstrap_requested=B,bootstrap_finite=rowSums(is.finite(htmtboot)),
                     method="Pearson absolute-correlation HTMT; percentile 95%; stratified within country"),
          "htmt_pooled_bootstrap_95ci.csv")
# Save conditional latent correlations of the chosen structural model.
corlv<-lavInspect(fits[[chosen_id]],"cor.lv")
write_out(data.frame(construct=rownames(corlv),corlv,check.names=FALSE),"chosen_latent_correlations.csv")

effect_vector <- function(fit,type="std.all") {
  ss<-standardizedSolution(fit,type=type,se=FALSE,zstat=FALSE,pvalue=FALSE,ci=FALSE)
  get<-function(lhs,rhs) {
    z<-ss$est.std[ss$op=="~"&ss$lhs==lhs&ss$rhs==rhs]
    if(length(z)!=1L) stop(paste("Missing path",lhs,rhs))
    z
  }
  ds<-get("PSB","SC");dl<-get("PSB","SBL")
  effect<-c(SBL_to_PSB=dl,SC_to_PSB=ds)
  vals<-numeric()
  for(k in names(values)) {
    a<-get(values[[k]],"CSR1_bin");b<-get("SC",values[[k]])
    cc<-get("SBL",values[[k]]);p<-get("PSB",values[[k]])
    effect[paste0("ind_",k,"_SC")]<-b*ds
    effect[paste0("ind_",k,"_SBL")]<-cc*dl
    effect[paste0("ind_",k,"_total")]<-b*ds+cc*dl
    vals[k]<-a*(p+b*ds+cc*dl)
    effect[paste0("CSR_via_",k)]<-vals[k]
  }
  effect["CSR_via_values"]<-sum(vals)
  effect["CSR_via_SC_direct"]<-get("SC","CSR1_bin")*ds
  effect["CSR_via_SBL_direct"]<-get("SBL","CSR1_bin")*dl
  effect["CSR_total_indirect"]<-sum(vals)+effect["CSR_via_SC_direct"]+effect["CSR_via_SBL_direct"]
  effect["CSR_direct_PSB"]<-get("PSB","CSR1_bin")
  effect["CSR_total_PSB"]<-effect["CSR_direct_PSB"]+effect["CSR_total_indirect"]
  effect
}
chosen<-fits[[chosen_id]]
stopifnot(admissible(chosen)$admissible)
point<-effect_vector(chosen)
point_nox<-effect_vector(chosen,"std.nox")
# Check the standardized defined effects against manual standardized products.
ss<-standardizedSolution(chosen,type="std.all",se=FALSE,zstat=FALSE,pvalue=FALSE,ci=FALSE)
shared<-intersect(names(point),ss$lhs[ss$op==":="])
check<-data.frame(effect=shared,manual=point[shared],
                  lavaan=ss$est.std[match(shared,ss$lhs)],
                  difference=point[shared]-ss$est.std[match(shared,ss$lhs)])
write_out(check,"standardized_indirect_consistency_check.csv")
stopifnot(max(abs(check$difference))<1e-7)
cat("Starting",B,"stratified bootstrap fits; cores",cores,"\n")
# ML and MLR have the same likelihood estimates; use ML without SE/test per replicate.
# Accept only converged, finite, positive-residual, positive-definite solutions.
boot_one<-function(i) {
  df<-dat[resample_rows[[i]],,drop=FALSE]
  res<-capture_fit(chosen_syn,df,estimator="ML",fast=TRUE)
  adm<-tryCatch(admissible(res$fit),error=function(e)data.frame(admissible=FALSE))
  valid<-isTRUE(adm$admissible)
  vals<-if(valid) tryCatch(effect_vector(res$fit),error=function(e) rep(NA_real_,length(point))) else rep(NA_real_,length(point))
  valsnox<-if(valid) tryCatch(effect_vector(res$fit,"std.nox"),error=function(e) rep(NA_real_,length(point_nox))) else rep(NA_real_,length(point_nox))
  valid<-valid && all(is.finite(vals)) && all(is.finite(valsnox))
  list(i=i,valid=valid,error=res$error,warnings=res$warnings,
       values=if(valid) vals else rep(NA_real_,length(point)),
       values_nox=if(valid) valsnox else rep(NA_real_,length(point_nox)))
}
checkpoint<-file.path(out,"bootstrap_checkpoint.rds")
stopifnot(identical(names(point),checkpoint_meta$point_names),identical(names(point_nox),checkpoint_meta$point_names))
sens_write_checkpoint(boot_results,checkpoint_meta,checkpoint)
cl<-parallel::makePSOCKcluster(cores)
parallel::clusterExport(cl,c("library_paths","capture_fit","admissible","effect_vector","dat",
                            "resample_rows","chosen_syn","point","point_nox","chosen","values","boot_one"))
parallel::clusterEvalQ(cl,{
 .libPaths(c(library_paths,.libPaths()))
 suppressPackageStartupMessages(library(lavaan))
 NULL
})
batch<-50L
remaining_starts<-if(length(boot_results)<B) seq.int(length(boot_results)+1L,B,by=batch) else integer()
for(start in remaining_starts) {
  ids<-start:min(B,start+batch-1L)
  got<-parallel::parLapply(cl,ids,boot_one)
  boot_results<-c(boot_results,got)
  cat("Bootstrap complete",length(boot_results),"of",B,
      "; valid",sum(vapply(boot_results,function(x)isTRUE(x$valid),logical(1))),"\n")
  flush.console()
  sens_write_checkpoint(boot_results,checkpoint_meta,checkpoint)
}
parallel::stopCluster(cl)
valid<-vapply(boot_results,function(x)isTRUE(x$valid),logical(1))
write_out(data.frame(replica=vapply(boot_results,function(x)x$i,integer(1)),
                     valid=valid,error=vapply(boot_results,function(x)x$error,character(1)),
                     warnings=vapply(boot_results,function(x)x$warnings,character(1))),
          "bootstrap_replicate_status.csv")
mat<-t(vapply(boot_results,function(x)x$values,point));colnames(mat)<-names(point)
matnox<-t(vapply(boot_results,function(x)x$values_nox,point_nox));colnames(matnox)<-names(point_nox)
write_out(data.frame(replica=seq_len(B),valid=valid,mat,check.names=FALSE),"bootstrap_standardized_replicates.csv")
write_out(data.frame(replica=seq_len(B),valid=valid,matnox,check.names=FALSE),"bootstrap_std_nox_replicates.csv")
boot_summary<-function(m,p,type) {
  qs<-t(apply(m[valid,,drop=FALSE],2,quantile,probs=c(.025,.975),na.rm=TRUE))
  data.frame(effect=names(p),standardization=type,point=as.numeric(p),bootstrap_SE=apply(m[valid,,drop=FALSE],2,sd),
             CI_low=qs[,1],CI_high=qs[,2],requested=B,valid=sum(valid),invalid=sum(!valid),
             method="Within-country stratified participant bootstrap; 95% percentile; ML point-estimation equals MLR")
}
write_out(rbind(boot_summary(mat,point,"std.all"),boot_summary(matnox,point_nox,"std.nox")),
          "bootstrap_indirects_and_key_paths_95ci.csv")
saveRDS(list(fits=fits,dat=dat,quality=quality,bootstrap=boot_results,syntax=chosen_syn,
             bootstrap_metadata=checkpoint_meta,bootstrap_resample_rows=resample_rows,run_id=run_id),
        file.path(out,"sensitivities_models_and_bootstrap.rds"))
cat("Done. Valid bootstrap replicas",sum(valid),"/",B,"\n")
source(file.path(script_dir,"finalize_sensitivity_tables.R"),local=FALSE)
source(file.path(script_dir,"export_focal_B_tables.R"),local=FALSE)
writeLines("COMPLETED",file.path(out,"run_state.txt"))
