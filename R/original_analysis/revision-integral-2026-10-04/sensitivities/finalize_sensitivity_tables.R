# Run after sensitivities_pooled.R; safe to run before bootstrap finishes for table preparation.
options(stringsAsFactors=FALSE,warn=1)
# Called by the main runner, or with explicit SENS_RUN_DIR and SENS_INPUT.
configured_libs<-c(Sys.getenv("SENS_AUDIT_LIB",""),Sys.getenv("SENS_USER_LIB",""))
.libPaths(c(configured_libs[nzchar(configured_libs)],.libPaths()))
suppressPackageStartupMessages(library(lavaan))
out<-Sys.getenv("SENS_RUN_DIR","")
input<-Sys.getenv("SENS_INPUT","")
if(!nzchar(out)||!dir.exists(out)||!nzchar(input)||!file.exists(input))
 stop("Finalization requires explicit existing SENS_RUN_DIR and SENS_INPUT")
write_out<-function(x,name)write.csv(x,file.path(out,name),row.names=FALSE,na="")
values<-c("Benevolencia","Universalismo","Poder","Logro")
primary_id<-"P3_country_partial_demographics"
primary<-readRDS(file.path(out,paste0(primary_id,".rds")))
syntax<-paste(readLines(file.path(out,paste0(primary_id,"_syntax.txt"))),collapse="\n")
full_sbl<-c("APS11","APS12","APS13","APS15","APS16","APS18","APS19","APS20","APS21")
verified_legacy<-c("APS11","APS12","APS13","APS16","APS18","APS19","APS21")
original<-haven::read_sav(input)
names(original)<-stringi::stri_trans_general(names(original),"Latin-ASCII")
score_check<-lapply(c("LES","Descartado_LES"),function(v) {
 items<-if(v=="LES") verified_legacy else full_sbl
 if(!v%in%names(original))return(data.frame(score=v,items=paste(items,collapse=";"),max_absolute_difference=NA,exact_within_1e_12=NA))
 rebuilt<-rowMeans(as.data.frame(lapply(original[,items],as.numeric)))
 delta<-as.numeric(original[[v]])-rebuilt
 data.frame(score=v,items=paste(items,collapse=";"),max_absolute_difference=max(abs(delta)),
            exact_within_1e_12=all(abs(delta)<1e-12))
})
write_out(do.call(rbind,score_check),"legacy_score_verification.csv")
corrected<-sub(paste0("SBL =~ ",paste(full_sbl,collapse=" + ")),
               paste0("SBL =~ ",paste(verified_legacy,collapse=" + ")),syntax,fixed=TRUE)
writeLines(corrected,file.path(out,"S_legacy_SBL7_syntax.txt"))
d<-as.data.frame(lavInspect(primary,"data"))
warns<-character()
legacy_path<-file.path(out,"S_legacy_SBL7.rds")
legacy<-if(file.exists(legacy_path))readRDS(legacy_path)else
 withCallingHandlers(sem(corrected,data=d,estimator="MLR",meanstructure=TRUE,fixed.x=TRUE),
  warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart("muffleWarning")})
stopifnot(lavInspect(legacy,"converged"),lavInspect(legacy,"post.check"))
saveRDS(legacy,file.path(out,"S_legacy_SBL7.rds"))
write_out(parameterEstimates(legacy,ci=TRUE),"S_legacy_SBL7_unstandardized_parameters_robust_ci.csv")
ss<-standardizedSolution(legacy,type="std.all",ci=TRUE)
write_out(ss,"S_legacy_SBL7_standardized_parameters_robust_ci.csv")
write_out(standardizedSolution(legacy,type="std.nox",ci=TRUE),"S_legacy_SBL7_standardized_nox_parameters_robust_ci.csv")
write_out(ss[ss$op=="~",],"S_legacy_SBL7_paths_standardized_robust_ci.csv")
write_out(ss[ss$op==":=",],"S_legacy_SBL7_indirects_standardized_robust_ci.csv")
r2<-lavInspect(legacy,"r2")
write_out(data.frame(outcome=names(r2),R2=as.numeric(r2)),"S_legacy_SBL7_r2.csv")
write_out(data.frame(model="S_legacy_SBL7",scenario="legacy_SBL7_all",error="",warnings=paste(unique(warns),collapse=" | ")),"S_legacy_SBL7_warnings.csv")
fitrow<-function(fit,id,scenario,N) {
 m<-c("chisq.scaled","df","pvalue.scaled","cfi.robust","tli.robust","rmsea.robust",
      "rmsea.ci.lower.robust","rmsea.ci.upper.robust","srmr","aic","bic","npar")
 fm<-vapply(m,function(z)unname(fitMeasures(fit,z)),numeric(1))
 pe<-parameterEstimates(fit,se=FALSE)
 resid<-pe$est[pe$op=="~~"&pe$lhs==pe$rhs]
 theta<-lavInspect(fit,"theta");psi<-lavInspect(fit,"est")$psi
 eig<-function(x)min(eigen(x,symmetric=TRUE,only.values=TRUE)$values)
 corlv<-lavInspect(fit,"cor.lv")
 maxcor<-max(abs(corlv[lower.tri(corlv)]))
 cbind(data.frame(model=id,scenario=scenario,N=N),as.data.frame(as.list(fm)),
       data.frame(converged=lavInspect(fit,"converged"),finite=all(is.finite(pe$est)),
                  min_residual_variance=min(resid),min_theta_eigen=eig(theta),min_psi_eigen=eig(psi),
                  max_abs_latent_correlation=maxcor,admissible=all(resid>0)&&eig(theta)>0&&eig(psi)>0&&maxcor<1))
}
write_out(fitrow(legacy,"S_legacy_SBL7","legacy_SBL7_all",nrow(d)),"S_legacy_SBL7_fit.csv")
id<-c("P0_no_country_full_chain","P1_country_full_chain","P2_country_partial_curricular",
      "P3_country_partial_demographics","P4_country_full_chain_demographics",
      "P5_no_country_partial_demographics","S_without_duplicate_copies",
      "S_without_constant22","S_without_duplicate_copies_and_constant22","S_legacy_SBL7")
fits<-setNames(lapply(id,function(x)readRDS(file.path(out,paste0(x,".rds")))),id)
scenarios<-setNames(c(rep("all",6),"without_duplicate_copies","without_constant22",
                      "without_duplicate_copies_and_constant22","legacy_SBL7_all"),id)
fitrows<-lapply(id,function(z)fitrow(fits[[z]],z,scenarios[[z]],lavInspect(fits[[z]],"nobs")))
for(i in seq_along(id))write_out(fitrows[[i]],paste0(id[i],"_fit.csv"))
write_out(do.call(rbind,fitrows),"models_fit_and_admissibility.csv")
focal<-function(fit) {
 z<-standardizedSolution(fit,type="std.all",ci=TRUE)
 z<-z[z$op=="~"&z$lhs%in%c(values,"SC","SBL","PSB") &
               z$rhs%in%c(values,"SC","SBL","CSR1_bin"),]
 stopifnot(nrow(z)==21L)
 z$p_holm_family21<-p.adjust(z$pvalue,method="holm",n=21L)
 z$family_size<-21L
 z
}
stability<-do.call(rbind,lapply(id,function(z){
 x<-focal(fits[[z]]);x$model<-z;x$N<-lavInspect(fits[[z]],"nobs");x}))
write_out(stability,"stability_key_paths_standardized.csv")
write_out(focal(primary),"primary_focal21_holm.csv")
ss<-standardizedSolution(primary,type="std.all",ci=TRUE)
indnames<-c(as.vector(outer(c("ben","uni","pow","log"),c("SC","SBL","total"),function(a,b)paste0("ind_",a,"_",b))),
           "CSR_via_values","CSR_via_SC_direct","CSR_via_SBL_direct","CSR_total_indirect")
ind<-ss[ss$op==":=" & ss$lhs%in%indnames,]
stopifnot(nrow(ind)==16L)
ind$p_holm_family16<-p.adjust(ind$pvalue,method="holm",n=16L)
ind$family_size<-16L
write_out(ind,"primary_indirects_family16_holm.csv")
write_out(data.frame(outcome=names(lavInspect(primary,"r2")),R2=as.numeric(lavInspect(primary,"r2"))),"primary_R2.csv")
compact<-stability[stability$lhs=="PSB"&stability$rhs%in%c("SBL","SC"),]
compact$fit_CFI<-vapply(compact$model,function(z)unname(fitMeasures(fits[[z]],"cfi.robust")),numeric(1))
compact$fit_TLI<-vapply(compact$model,function(z)unname(fitMeasures(fits[[z]],"tli.robust")),numeric(1))
compact$fit_RMSEA<-vapply(compact$model,function(z)unname(fitMeasures(fits[[z]],"rmsea.robust")),numeric(1))
compact$fit_SRMR<-vapply(compact$model,function(z)unname(fitMeasures(fits[[z]],"srmr")),numeric(1))
write_out(compact,"stability_SC_SBL_compact.csv")
writeLines(c(
 "Holm family 21: 4 CSR1->values + 8 values->SC/SBL + 4 values->PSB + 2 SC/SBL->PSB + 3 CSR1->SC/SBL/PSB.",
 "Country and demographic control coefficients are excluded from that substantive family.",
 "Holm family 16: 12 value indirects (SC, SBL, total for each of four values) + CSR via values + CSR via SC directly + CSR via SBL directly + total CSR indirect.",
 "These overlapping indirects are reported together; Holm is valid under arbitrary dependence.",
 "p_holm refers to robust delta tests of standardized parameters; bootstrap percentile IC95 remain pointwise, not multiplicity adjusted.",
 "No model or item set is selected because a p value passes Holm.",
 "Verified legacy LES7 uses APS11,12,13,16,18,19,21 and omits APS15 and APS20.",
 "A preliminary seven-item set supplied before LES verification was replaced and is not included in final summary tables."
),file.path(out,"multiplicity_and_legacy_notes.txt"))
# If the main bootstrap bundle already exists, make its legacy fit consistent.
bundle<-file.path(out,"sensitivities_models_and_bootstrap.rds")
if(file.exists(bundle)) {
 z<-readRDS(bundle);z$fits$S_legacy_SBL7<-legacy;saveRDS(z,bundle)
}
cat("Postprocessing complete; correct LES7; Holm21 and Holm16.\n")
