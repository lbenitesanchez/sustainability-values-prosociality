# Read-only statistical postprocessing of persisted P3; no model or bootstrap is fitted.
configured_libs<-c(Sys.getenv("SENS_AUDIT_LIB",""),Sys.getenv("SENS_USER_LIB",""))
.libPaths(c(configured_libs[nzchar(configured_libs)],.libPaths()))
suppressPackageStartupMessages(library(lavaan))
out<-Sys.getenv("SENS_RUN_DIR","")
if(!nzchar(out)||!dir.exists(out))stop("Set SENS_RUN_DIR explicitly")
primary<-readRDS(file.path(out,"P3_country_partial_demographics.rds"))
pe<-parameterEstimates(primary,ci=TRUE,standardized=FALSE)
ss<-standardizedSolution(primary,type="std.all",ci=TRUE)
key<-function(z)paste(z$lhs,z$op,z$rhs,sep="|")
vals<-c("Benevolencia","Universalismo","Poder","Logro")
sel<-pe$op=="~"&pe$lhs%in%c(vals,"SC","SBL","PSB")&pe$rhs%in%c(vals,"SC","SBL","CSR1_bin")
b<-pe[sel,];s<-ss[match(key(b),key(ss)),];stopifnot(nrow(b)==21L,!anyNA(s$lhs))
join<-function(b,s,n) data.frame(lhs=b$lhs,op=b$op,rhs=b$rhs,
 B=b$est,SE_B=b$se,z_B=b$z,p_B=b$pvalue,CI_B_low=b$ci.lower,CI_B_high=b$ci.upper,
 beta=s$est.std,SE_beta=s$se,z_beta=s$z,p_beta=s$pvalue,
 CI_beta_low=s$ci.lower,CI_beta_high=s$ci.upper,
 p_B_Holm=p.adjust(b$pvalue,"holm",n=n),p_beta_Holm=p.adjust(s$pvalue,"holm",n=n),
 family_size=n,primary_test="Robust Wald z of unstandardized B; Holm")
z<-join(b,s,21L)
write.csv(z,file.path(out,"primary_focal21_holm_B.csv"),row.names=FALSE)
indnames<-c(as.vector(outer(c("ben","uni","pow","log"),c("SC","SBL","total"),function(a,b)paste0("ind_",a,"_",b))),
 "CSR_via_values","CSR_via_SC_direct","CSR_via_SBL_direct","CSR_total_indirect")
b<-pe[pe$op==":="&pe$lhs%in%indnames,];s<-ss[match(key(b),key(ss)),]
stopifnot(nrow(b)==16L,!anyNA(s$lhs))
write.csv(join(b,s,16L),file.path(out,"primary_indirects_family16_holm_B.csv"),row.names=FALSE)
writeLines(c(
 "primary_focal21_holm_B.csv: 21 focal coefficients; primary contrast p_B and p_B_Holm from robust Wald z of B.",
 "B, SE_B and CI_B_* describe the unstandardized parameter. beta, SE_beta and CI_beta_* describe std.all.",
 "p_beta and p_beta_Holm are separate delta-method tests of the standardized estimand. They need not equal p_B.",
 "primary_focal21_holm.csv preserves the earlier std.all contrast table; do not relabel its pvalue as p_B.",
 "primary_indirects_family16_holm_B.csv applies the same clear separation to the 16 defined indirects.",
 "All confidence intervals in these tables are pointwise robust 95% intervals; Holm adjusts p values only.",
 "Bootstrap pointwise percentile intervals are in bootstrap_indirects_and_key_paths_95ci.csv, with the correct standardization field.",
 "No data, specification, sample or model-selection rule changed during this export."
 ),file.path(out,"B_and_standardized_inference_notes.txt"))
cat("Exported B and standardized estimands separately; Holm21 significant:",
 paste(paste(z$rhs[z$p_B_Holm<.05],"->",z$lhs[z$p_B_Holm<.05]),collapse="; "),"\n")
