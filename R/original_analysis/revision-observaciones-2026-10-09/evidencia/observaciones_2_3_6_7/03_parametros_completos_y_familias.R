# Exportaciones completas de las siete especificaciones, sin nuevos ajustes.
# Este bloque permite localizar todos los coeficientes de tablas 2/S1-S8 y
# distinguir IC puntuales estandarizados de p robustos no estandarizados.
# Copia originales agregados/sintaxis y regenera familias focales/indirectas.
proyecto <- '/Users/luisbenites/Documents/Codex/Paper con MA - 1/reproduccion-integral-V3-2026-10-10/proyecto'
.libPaths(c(file.path(proyecto,'revision-revista-2026-10-04/auditoria-codigos/R-library'),.libPaths()))
suppressPackageStartupMessages(library(lavaan));library(digest);library(jsonlite)
fuente <- file.path(proyecto,'sem-alternativas-2026-10-08/resultados')
salida <- file.path(proyecto,'revision-observaciones-2026-10-09/evidencia/observaciones_2_3_6_7/parametros_completos')
dir.create(salida,recursive=TRUE,showWarnings=FALSE)
archivos <- list.files(fuente,pattern='(_parameters\\.csv|_standardized\\.csv|_indirects\\.csv|_syntax\\.txt|^focal_paths_all_models\\.csv)$',full.names=TRUE)
manifest <- do.call(rbind,lapply(archivos,function(p){
  h<-digest(file=p,algo='sha256'); destino<-file.path(salida,basename(p))
  stopifnot(file.copy(p,destino,overwrite=TRUE),identical(h,digest(file=destino,algo='sha256')))
  data.frame(source=p,copy=basename(p),sha256=h,source_unchanged=identical(h,digest(file=p,algo='sha256')))
}))
write.csv(manifest,file.path(salida,'manifest_parametros.csv'),row.names=FALSE)
modelos<-c('M0_original','M1_ASP_unico','M2_puentes2','M3_puentes4','M4_S1_refSC','M5_S1_refSBL','M6_ESEM_target')
todos <- focales <- indirectos <- list(); hashes<-list()
for(m in modelos) for(e in c('MLR','WLSMV')){
  p<-file.path(fuente,paste0(m,'_SEM_',e,'.rds'));h<-digest(file=p,algo='sha256');f<-readRDS(p)
  pe<-parameterEstimates(f,ci=TRUE,level=.95)
  st<-suppressWarnings(standardizedSolution(f,type='std.all',ci=TRUE,level=.95))
  clave<-function(x)paste(x$lhs,x$op,x$rhs,sep='|')
  j<-match(clave(pe),clave(st))
  pe$beta<-st$est.std[j];pe$beta_ci_lower<-st$ci.lower[j];pe$beta_ci_upper<-st$ci.upper[j]
  pe$model<-m;pe$estimator<-e;pe$N<-sum(lavInspect(f,'nobs'));pe$ci_level<-.95
  pe$inference_scope<-if(m=='M6_ESEM_target'&&e=='WLSMV') 'Archived solution; failed robust global test and negative curvature prevent ordinal structural confirmation' else 'Pointwise conditional inference; no model-selection adjustment'
  # Todas las regresiones, incluidos controles, están disponibles en esta tabla.
  x<-pe[pe$op=='~',];todos[[paste(m,e)]]<-x
  # Los caminos focales tienen etiquetas; los controles no. Familias 21/15.
  x<-x[nzchar(x$label),];stopifnot(nrow(x)==if(m=='M1_ASP_unico')15 else 21)
  x$p_Holm_focal<-p.adjust(x$pvalue,method='holm');x$Holm_family_size<-nrow(x)
  focales[[paste(m,e)]]<-x
  # Indirectos definidos por :=: familias 16/11, p de parameterEstimates.
  x<-pe[pe$op==':=',];stopifnot(nrow(x)==if(m=='M1_ASP_unico')11 else 16)
  x$p_Holm_indirect<-p.adjust(x$pvalue,method='holm');x$Holm_family_size<-nrow(x)
  indirectos[[paste(m,e)]]<-x
  stopifnot(identical(h,digest(file=p,algo='sha256')))
  hashes[[paste(m,e)]]<-list(path=p,sha256=h)
}
# ESEM incluye columnas de bloque extra: preservarlas y completar NA en otros
# modelos evita confundir parámetros sin modificar ninguna estimación.
unir <- function(lista){
  columnas<-unique(unlist(lapply(lista,names)))
  lista<-lapply(lista,function(x){for(k in setdiff(columnas,names(x)))x[[k]]<-NA;x[,columnas,drop=FALSE]})
  do.call(rbind,lista)
}
write.csv(unir(todos),file.path(salida,'regenerado_todas_regresiones_b_beta_IC95.csv'),row.names=FALSE,na='')
write.csv(unir(focales),file.path(salida,'regenerado_focales_21_15_Holm_b_beta_IC95.csv'),row.names=FALSE,na='')
write.csv(unir(indirectos),file.path(salida,'regenerado_indirectos_16_11_Holm_b_beta_IC95.csv'),row.names=FALSE,na='')
write_json(list(status='PASS_NO_NEW_FITS',aggregate_source_files=nrow(manifest),all7models_both_estimators=TRUE,
 originals_unchanged=TRUE,source_fits=hashes,ci_level=.95,
 adjustment='Holm of unstandardized robust p, separately per model/estimator; no selection protection',
 M6_ordinal_SEM='Parameters retained for transparency, not ordinal structural confirmation'),
 file.path(salida,'reporte_parametros_completos.json'),pretty=TRUE,auto_unbox=TRUE)
cat('PASS: 78 aggregate/syntax exports; full regression, focal and indirect families regenerated.\n')
