# REVISIÓN 3: sensibilidad de calidad de registros, sin nuevos datos.
# Ejecutar: Rscript 01_modelos_centrales_sensibilidad.R
# El análisis principal conserva N=653. No se diagnostica inatención ni identidad.
# Las tres exclusiones son sensibilidades fijadas antes de sus nuevos ajustes.
# No se eligen modelos, ítems, rutas o casos según significación estadística.
options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
project <- Sys.getenv("RCP_PROYECTO", "/Users/luisbenites/Documents/Codex/Paper con MA - 1/reproduccion-integral-V3-2026-10-10/proyecto")
root <- file.path(project, "revision-3-revista-colombiana-psicologia-2026-10-08/analisis_sensibilidad")
out <- file.path(root, "modelos_centrales")
private <- file.path(root, "privados")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(private, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(file.path(project, "revision-revista-2026-10-04/auditoria-codigos/R-library"), .libPaths()))
suppressPackageStartupMessages({library(lavaan); library(semTools); library(jsonlite); library(digest)})
state <- file.path(out, "estado.txt")
writeLines("RUNNING", state)
options(error = function() {writeLines(c("FAILED", geterrmessage()), state); quit(status=1,save="no")})
paths <- c(
 data = file.path(project, "revision-integral-2026-10-04/analysis/outputs/canonical-full-20261004-02/01_audit/analytic_data.rds"),
 flags = file.path(project, "revision-integral-2026-10-04/analysis/outputs/canonical-full-20261004-02/01_audit/quality_flags_no_auto_exclusion.csv"),
 p3 = file.path(project, "revision-integral-2026-10-04/sensitivities/P3_country_partial_demographics.rds"),
 folds = file.path(project, "sem-items-ordinal-2026-10-08/prediccion/resultados/privados/asignaciones_cv.rds"))
sha <- function(p) digest(p, algo="sha256", file=TRUE)
hash_before <- vapply(paths, sha, character(1))
audited <- as.data.frame(readRDS(paths[["data"]]))
data <- as.data.frame(lavInspect(readRDS(paths[["p3"]]), "data"))
common <- intersect(names(data), names(audited))
stopifnot(nrow(data)==653L,nrow(audited)==653L,
 max(abs(as.matrix(data[,common])-as.matrix(audited[,common])))==0)
data$Pais <- audited$Pais
pp <- paste0("PP",1:3)
sc <- c("APS1","APS2","APS3","APS4","APS5","APS6","APS7","APS9","APS10","APS14")
sbl <- c("APS11","APS12","APS13","APS15","APS16","APS18","APS19","APS20","APS21")
aps <- c(sc,sbl); indicators <- c(aps,pp)
values <- c("Benevolencia","Universalismo","Poder","Logro")
focal <- c("SC","SBL",values)
# Los primeros 100 campos son las variables originales; los otros cinco
# son recodificaciones de auditoría. Codigo NO forma parte de la comparación.
profile_vars <- setdiff(names(audited)[seq_len(100L)], "Codigo")
stopifnot(length(profile_vars)==99L)
profile <- audited[,profile_vars,drop=FALSE]
dup_copy <- duplicated(profile)
dup_member <- dup_copy | duplicated(profile,fromLast=TRUE)
constant22 <- apply(data[,indicators],1,function(z)length(unique(z))==1L)
archived <- read.csv(paths[["flags"]])
stopifnot(identical(dup_copy,archived$duplicate_full_profile_copy),
 identical(dup_member,archived$duplicate_full_profile_member),
 identical(constant22,archived$sc_sbl_psb_zero_within_person_sd),
 sum(dup_copy)==3L,sum(dup_member)==5L,sum(constant22)==35L)
# Clave con los mismos 99 campos que la auditoría anterior, comprobando que
# genera exactamente las mismas banderas que duplicated(data.frame).
keys <- apply(profile,1,paste,collapse="|")
group <- match(keys,unique(keys))
stopifnot(identical(duplicated(group),dup_copy))
groups <- split(seq_len(nrow(data)),group)
duplicate_groups <- groups[lengths(groups)>1L]
stopifnot(length(duplicate_groups)==2L,all(vapply(groups,function(ix)
 length(unique(data$Pais[ix]))==1L,logical(1))))
scenarios <- list(all=rep(TRUE,653L),without_duplicate_copies=!dup_copy,
 without_constant22=!constant22,without_both=!dup_copy & !constant22)
sample_counts <- do.call(rbind,lapply(names(scenarios),function(s)
 data.frame(scenario=s,N=sum(scenarios[[s]]),excluded=sum(!scenarios[[s]]),
  duplicate_copies_remaining=sum(dup_copy & scenarios[[s]]),
  constant22_remaining=sum(constant22 & scenarios[[s]]))))
write.csv(sample_counts,file.path(root,"escenarios_muestra.csv"),row.names=FALSE)
historical <- readRDS(paths[["folds"]])
leak <- do.call(rbind,lapply(seq_along(historical),function(r){
 a <- historical[[r]]
 spread <- vapply(duplicate_groups,function(ix)length(unique(a[ix]))>1L,logical(1))
 data.frame(repetition=r,n_identical_sets=length(duplicate_groups),
  sets_across_folds=sum(spread),members_exposed=sum(lengths(duplicate_groups)[spread]),
  note="Possible dependence if identical records belong to one person; identity unverified")
}))
write.csv(leak,file.path(root,"auditoria_fuga_particiones_historicas.csv"),row.names=FALSE)
quality <- list(data=data,group=group,dup_copy=dup_copy,dup_member=dup_member,
 constant22=constant22,scenarios=scenarios,original_row=seq_len(653L))
saveRDS(quality,file.path(private,"datos_y_banderas_para_sensibilidades.rds"))
protocol <- list(recorded_at=format(Sys.time(),tz="UTC",usetz=TRUE),
 status="Exploratory specifications fixed before these new fits; earlier results known",
 no_new_data=TRUE,primary_N=653L,sensitivity_N=sample_counts,
 duplicate_rule="Identical on all 99 original fields excluding Codigo; keep first in existing row order; participant identity unverified",
 constant_rule="One identical response across 19 retained APS and PP1-PP3 (22 indicators); no inference of inattention",
 duplicate_sets=2L,duplicate_members=5L,duplicate_copies=3L,constant_profiles=35L,
 adjustment="Original P3 affine recodings and controls preserved; no redefinition by scenario",
 models=c("M2 CFA","M2 latent SEM","M2 with three correlated observed PP outcomes"),
 estimators=c("MLR","WLSMV theta"),new_fits=18L,
 multiplicity="Holm separately within each estimator/scenario: 21 latent-SEM focal paths (before filtering), 18 item slopes, 6 joint predictors, 6 MLR heterogeneity tests; selection family not covered",
 no_selection="No deletion based on coefficients, p-values, fit or category rarity; exclusions only these sensitivities",
 sources_sha256=as.list(hash_before),script_sha256=sha(file.path(root,"01_modelos_centrales_sensibilidad.R")))
write_json(protocol,file.path(root,"protocolo_sensibilidad_antes_ajustes.json"),pretty=TRUE,auto_unbox=TRUE,digits=NA)
references <- list()
for(id in c("CFA_MLR","CFA_WLSMV","SEM_MLR","SEM_WLSMV"))
 references[[id]] <- readRDS(file.path(project,"sem-alternativas-2026-10-08/resultados",paste0("M2_puentes2_",id,".rds")))
for(e in c("MLR","WLSMV"))references[[paste0("ITEMS_",e)]] <- readRDS(file.path(project,
 "sem-items-ordinal-2026-10-08/items/privados",paste0("SEM_items_",e,".rds")))
safe <- function(expr)tryCatch(expr,error=function(e)NULL)
ev <- function(x)if(is.null(x)||!is.matrix(x)||!all(is.finite(x)))NA_real_ else
 min(eigen((x+t(x))/2,symmetric=TRUE,only.values=TRUE)$values)
fit_records <- paths_records <- joint_records <- hetero_records <- measure_records <- list()
fit_status <- list(); counter <- 0L
csv <- function(x,name)write.csv(x,file.path(out,name),row.names=FALSE,na="")
for(s in names(scenarios))for(id in names(references)){
 counter <- counter+1L
 e <- if(grepl("WLSMV",id))"WLSMV" else "MLR"
 dd <- data[scenarios[[s]],,drop=FALSE]
 cache <- file.path(private,paste0(s,"__",id,".rds"))
 warns <- character(); error <- ""; seconds <- 0
 if(s=="all")f <- references[[id]] else {
  cat("FIT",counter,s,id,"N",nrow(dd),"\n");flush.console()
  tm <- system.time(f <- tryCatch(withCallingHandlers(
   stats::update(references[[id]],data=dd),warning=function(w){
    warns <<- c(warns,conditionMessage(w));invokeRestart("muffleWarning")}),
   error=function(err){error <<- conditionMessage(err);NULL}))
  seconds <- unname(tm[["elapsed"]])
 }
 if(is.null(f)){
  fit_status[[paste(s,id)]] <- list(scenario=s,model=id,status="FAILED",error=error,warnings=unique(warns))
  write_json(fit_status,file.path(out,"estado_ajustes.json"),pretty=TRUE,auto_unbox=TRUE)
  next
 }
 saveRDS(f,cache)
 pe <- parameterEstimates(f,ci=TRUE); ss <- standardizedSolution(f,ci=TRUE)
 pe$scenario <- s;pe$model <- id;pe$estimator <- e
 csv(pe,paste0(s,"__",id,"__parametros.csv"))
 ss$scenario <- s;ss$model <- id;ss$estimator <- e
 csv(ss,paste0(s,"__",id,"__estandarizados.csv"))
 fm <- safe(fitMeasures(f,c("chisq","df","chisq.scaled","df.scaled","cfi","tli","rmsea", "cfi.scaled","tli.scaled","rmsea.scaled","cfi.robust","tli.robust","rmsea.robust","srmr","npar")))
 if(is.null(fm))fm <- setNames(rep(NA_real_,15L),c("chisq","df","chisq.scaled","df.scaled","cfi","tli","rmsea","cfi.scaled","tli.scaled","rmsea.scaled","cfi.robust","tli.robust","rmsea.robust","srmr","npar"))
 pt <- parTable(f);free <- pt$free>0
 est <- lavInspect(f,"est");vc <- safe(lavInspect(f,"vcov"))
 allfinite <- all(is.finite(pt$est[free])) && all(is.finite(pt$se[free])) && all(pt$se[free]>0)
 rec <- data.frame(scenario=s,model=id,estimator=e,N=nrow(dd),as.list(fm),
  converged=lavInspect(f,"converged"),post_check=lavInspect(f,"post.check"),
  free_parameters_and_SE_finite=allfinite,min_vcov_eigen=ev(vc),
  min_psi_eigen=ev(est$psi),min_theta_eigen=ev(est$theta),
  elapsed_seconds=seconds,warnings=paste(unique(warns),collapse=" | "),check.names=FALSE)
 fit_records[[paste(s,id)]] <- rec
 fit_status[[paste(s,id)]] <- list(scenario=s,model=id,status="COMPLETED",converged=rec$converged,
  post_check=rec$post_check,free_parameters_and_SE_finite=allfinite,warnings=unique(warns))
 if(grepl("^CFA",id)){
  lv <- lavInspect(f,"cor.lv")
  rel <- safe(semTools::compRelSEM(f,tau.eq=FALSE,ord.scale=(e=="WLSMV")))
  ave <- safe(semTools::AVE(f,obs.var=TRUE))
  measure_records[[paste(s,id)]] <- data.frame(scenario=s,model=id,estimator=e,N=nrow(dd),
   SC_SBL_correlation=lv["SC","SBL"],PSB_AVE=if(is.null(ave))NA_real_ else as.numeric(ave["PSB"]),
   PSB_omega_observed_categories=if(is.null(rel))NA_real_ else as.numeric(unlist(rel)["PSB"]),
   AVE_scale=if(e=="WLSMV")"Underlying ordinal responses" else "Numeric category scores")
 }
 if(grepl("^SEM",id)){
  # Reproducir la familia original COMPLETA de 21 rutas antes de filtrar:
  # CSR1 -> cuatro valores, SC, SBL y PSB (7); valores -> SC/SBL (8),
  # valores -> PSB (4); SC/SBL -> PSB (2). No se ajusta sólo una subtabla.
  keep <- pe$op=="~" & ((pe$lhs %in% values & pe$rhs=="CSR1_bin") |
    (pe$lhs %in% c("SC","SBL") & pe$rhs %in% c(values,"CSR1_bin")) |
    (pe$lhs=="PSB" & pe$rhs %in% c(focal,"CSR1_bin")))
  pa <- pe[keep,,drop=FALSE]
  stopifnot(nrow(pa)==21L)
  ii <- match(paste(pa$lhs,pa$op,pa$rhs),paste(ss$lhs,ss$op,ss$rhs))
  pa$beta <- ss$est.std[ii];pa$p_Holm_family <- p.adjust(pa$pvalue,"holm")
  pa$Holm_family_size <- 21L
  paths_records[[paste(s,id)]] <- pa
 }
 if(grepl("^ITEMS",id)){
  pa <- pe[pe$op=="~" & pe$lhs %in% pp & pe$rhs %in% focal,,drop=FALSE]
  stopifnot(nrow(pa)==18L)
  ii <- match(paste(pa$lhs,pa$op,pa$rhs),paste(ss$lhs,ss$op,ss$rhs))
  pa$beta <- ss$est.std[ii];pa$p_Holm_family <- p.adjust(pa$pvalue,"holm")
  pa$Holm_family_size <- 18L
  paths_records[[paste(s,id)]] <- pa
  joint <- do.call(rbind,lapply(focal,function(pr){
   labs <- pe$label[pe$op=="~" & pe$lhs %in% pp & pe$rhs==pr]
   w <- safe(lavTestWald(f,paste(paste0(labs," == 0"),collapse="\n")))
   data.frame(scenario=s,model=id,estimator=e,predictor=pr,
    statistic=if(is.null(w))NA_real_ else w$stat,df=if(is.null(w))NA_real_ else w$df,
    pvalue=if(is.null(w))NA_real_ else w$p.value)
  }))
  joint$p_Holm_6 <- p.adjust(joint$pvalue,"holm");joint_records[[paste(s,id)]] <- joint
  if(e=="MLR"){
   h <- do.call(rbind,lapply(focal,function(pr){
    z <- pe[pe$op=="~" & pe$lhs %in% pp & pe$rhs==pr,]
    labs <- z$label[match(pp,z$lhs)]
    w <- safe(lavTestWald(f,paste(paste0(labs[1]," == ",labs[2:3]),collapse="\n")))
    data.frame(scenario=s,model=id,estimator=e,predictor=pr,
     statistic=if(is.null(w))NA_real_ else w$stat,df=if(is.null(w))NA_real_ else w$df,
     pvalue=if(is.null(w))NA_real_ else w$p.value)
   }))
   h$p_Holm_6 <- p.adjust(h$pvalue,"holm");hetero_records[[paste(s,id)]] <- h
  }
 }
 csv(do.call(rbind,fit_records),"indices_y_diagnosticos_modelos.csv")
 write_json(fit_status,file.path(out,"estado_ajustes.json"),pretty=TRUE,auto_unbox=TRUE)
}
csv(do.call(rbind,paths_records),"rutas_focales_sensibilidad.csv")
csv(do.call(rbind,joint_records),"Wald_conjuntos_6_Holm.csv")
csv(do.call(rbind,hetero_records),"Wald_heterogeneidad_MLR_6_Holm.csv")
csv(do.call(rbind,measure_records),"medicion_CFA_sensibilidad.csv")
after <- vapply(paths,sha,character(1));stopifnot(identical(hash_before,after))
write_json(list(status="COMPLETED",reference_objects_reused=6L,new_sensitivity_fits=18L,
 primary_N_retained=653L,originals_unchanged=TRUE,source_sha256_before=as.list(hash_before),
 source_sha256_after=as.list(after),all_requested_model_records=length(fit_status),
 caveat="Sensitivity flags do not verify identity or inattention; causal, convergence-validity and cross-site limits remain"),
 file.path(out,"auditoria_final.json"),pretty=TRUE,auto_unbox=TRUE,digits=NA)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
writeLines("COMPLETED",state)
cat("COMPLETED CENTRAL SENSITIVITIES\n")
