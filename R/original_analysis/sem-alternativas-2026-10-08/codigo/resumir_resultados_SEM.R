# Extrae evidencia agregada para las tablas del artículo y el informe.
# No cambia, selecciona ni vuelve a ajustar ningún modelo.
options(stringsAsFactors=FALSE)
project <- '/Users/luisbenites/Documents/Codex/Paper con MA - 1/reproduccion-integral-V3-2026-10-10/proyecto'
.libPaths(c(file.path(project,'revision-revista-2026-10-04/auditoria-codigos/R-library'),
            '/Users/luisbenites/Library/R/arm64/4.6/library',.libPaths()))
suppressPackageStartupMessages(library(lavaan))
r <- file.path(project,'sem-alternativas-2026-10-08/resultados')
ids <- c('M0_original','M1_ASP_unico','M2_puentes2','M3_puentes4','M4_S1_refSC','M5_S1_refSBL','M6_ESEM_target')
rows <- list();ri <- list();locals <- list();j <- 0L
for(id in ids)for(kind in c('CFA','SEM'))for(est in c('MLR','WLSMV')) {
  f <- readRDS(file.path(r,paste(id,kind,est,sep='_') |> paste0('.rds')))
  # Guardar cada correlación con sus dos nombres evita confundir la mayor
  # correlación del modelo con la correlación concreta entre SC y SBL.
  z <- lavInspect(f,'cor.lv')
  for(a in seq_len(nrow(z)-1L))for(b in seq.int(a+1L,nrow(z))) {
    j<-j+1L;rows[[j]]<-data.frame(model=id,kind=kind,estimator=est,
      factor1=rownames(z)[a],factor2=colnames(z)[b],correlation=z[a,b])
  }
  if(kind=='SEM') {
    rr<-lavInspect(f,'rsquare')
    ri[[paste(id,est)]]<-data.frame(model=id,estimator=est,variable=names(rr),R_squared=as.numeric(rr))
  }
  if(est=='MLR') {
    # Sólo 22 indicadores APS/PP, para que el número de pares sea comparable.
    pp <- c('PP1','PP2','PP3')
    aps <- c('APS1','APS2','APS3','APS4','APS5','APS6','APS7','APS9','APS10','APS14',
             'APS11','APS12','APS13','APS15','APS16','APS18','APS19','APS20','APS21')
    observed<-lavInspect(f,'sampstat')$cov
    implied<-lavInspect(f,'implied')$cov
    variables<-c(aps,pp)
    # Residuals en unidades de correlación utilizando la desviación observada.
    res<-(observed[variables,variables]-implied[variables,variables])/
      sqrt(outer(diag(observed[variables,variables]),diag(observed[variables,variables])))
    v<-res[lower.tri(res)]
    locals[[paste(id,kind)]]<-data.frame(model=id,kind=kind,estimator=est,pairs=length(v),
      residual_pairs_abs_above_0_10=sum(abs(v)>.1),max_abs_residual=max(abs(v)))
  }
}
write.csv(do.call(rbind,rows),file.path(r,'latent_correlations_named.csv'),row.names=FALSE)
write.csv(do.call(rbind,ri),file.path(r,'R_squared_all_models.csv'),row.names=FALSE)
write.csv(do.call(rbind,locals),file.path(r,'local_residual_summary.csv'),row.names=FALSE)
# La incertidumbre de diferencias CV se describe con dispersión de las cinco
# particiones. No son cinco muestras independientes ni intervalos de confianza.
cv <- read.csv(file.path(r,'cross_validation/cv_repeat_summary.csv'))
base <- cv[cv$model=='M0_original',]
paired <- do.call(rbind,lapply(ids,function(id){
 z<-cv[cv$model==id,];b<-base[match(z$repetition,base$repetition),]
 data.frame(model=id,repetition=z$repetition,
   delta_RMSE_vs_M0=z$RMSE_PP_mean-b$RMSE_PP_mean,
   delta_joint_logscore_vs_M0=z$mean_joint_logscore-b$mean_joint_logscore)
}))
write.csv(paired,file.path(r,'cross_validation/cv_paired_differences_vs_M0.csv'),row.names=FALSE)
cat('Aggregate named correlations R squared local residuals and paired CV exported\n')
