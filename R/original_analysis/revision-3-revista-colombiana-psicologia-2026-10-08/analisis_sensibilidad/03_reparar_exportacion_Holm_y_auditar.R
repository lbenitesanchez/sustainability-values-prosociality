# Reparación de exportación: conserva todos los ajustes y coeficientes.
# El SEM latente utiliza la familia original de 21 rutas; los PP por ítem
# utilizan sus 18 pendientes. La corrección no vuelve a estimar modelos.
options(stringsAsFactors=FALSE)
project <- Sys.getenv("RCP_PROYECTO", "/Users/luisbenites/Documents/Codex/Paper con MA - 1/reproduccion-integral-V3-2026-10-10/proyecto")
root <- file.path(project,"revision-3-revista-colombiana-psicologia-2026-10-08/analisis_sensibilidad")
out <- file.path(root,"modelos_centrales");private <- file.path(root,"privados")
.libPaths(c(file.path(project,"revision-revista-2026-10-04/auditoria-codigos/R-library"),.libPaths()))
suppressPackageStartupMessages({library(lavaan);library(jsonlite);library(digest)})
q <- readRDS(file.path(private,"datos_y_banderas_para_sensibilidades.rds"))
values <- c("Benevolencia","Universalismo","Poder","Logro");focal <- c("SC","SBL",values)
pp <- paste0("PP",1:3)
paths <- list();hashes <- list();checks <- list()
script <- file.path(root,"03_reparar_exportacion_Holm_y_auditar.R")
write_json(list(recorded_at=format(Sys.time(),tz="UTC",usetz=TRUE),
 reason="Latent-SEM Holm exported on 14 selected paths; restore original complete 21-path family before filtering. Item-family18 unchanged. No new fits or coefficient changes.",
 new_fits=0L,script_sha256=digest(script,algo="sha256",file=TRUE)),
 file.path(out,"protocolo_reparacion_Holm_antes_exports.json"),pretty=TRUE,auto_unbox=TRUE)
for(s in names(q$scenarios))for(id in c("SEM_MLR","SEM_WLSMV","ITEMS_MLR","ITEMS_WLSMV")){
 path <- file.path(private,paste0(s,"__",id,".rds"))
 before <- digest(path,algo="sha256",file=TRUE);f <- readRDS(path)
 pe <- parameterEstimates(f,ci=TRUE);ss <- standardizedSolution(f,ci=TRUE)
 pe$scenario <- s;pe$model <- id;pe$estimator <- if(grepl("WLSMV",id))"WLSMV" else "MLR"
 if(grepl("^SEM",id)){
  keep <- pe$op=="~" & ((pe$lhs %in% values & pe$rhs=="CSR1_bin") |
   (pe$lhs %in% c("SC","SBL") & pe$rhs %in% c(values,"CSR1_bin")) |
   (pe$lhs=="PSB" & pe$rhs %in% c(focal,"CSR1_bin")))
  z <- pe[keep,];family <- 21L
 }else{z <- pe[pe$op=="~" & pe$lhs %in% pp & pe$rhs %in% focal,];family <- 18L}
 stopifnot(nrow(z)==family)
 ix <- match(paste(z$lhs,z$op,z$rhs),paste(ss$lhs,ss$op,ss$rhs))
 stopifnot(!anyNA(ix));z$beta <- ss$est.std[ix]
 z$p_Holm_family <- p.adjust(z$pvalue,"holm");z$Holm_family_size <- family
 paths[[paste(s,id)]] <- z
 after <- digest(path,algo="sha256",file=TRUE);stopifnot(identical(before,after))
 hashes[[paste(s,id)]] <- before
 checks[[paste(s,id)]] <- list(N=lavInspect(f,"nobs"),Holm_family=family,
  converged=lavInspect(f,"converged"),post_check=lavInspect(f,"post.check"))
}
old <- read.csv(file.path(out,"rutas_focales_sensibilidad.csv"))
new <- do.call(rbind,paths)
key <- function(x)paste(x$scenario,x$model,x$lhs,x$op,x$rhs)
ix <- match(key(old),key(new));stopifnot(!anyNA(ix))
max_diff <- max(abs(old$est-new$est[ix]),abs(old$se-new$se[ix]),abs(old$pvalue-new$pvalue[ix]))
stopifnot(max_diff<1e-12)
write.csv(new,file.path(out,"rutas_focales_sensibilidad.csv"),row.names=FALSE,na="")
# La ruta Power -> SBL del M2 principal recupera el ajuste correcto sin
# asignarle significación por usar una familia truncada en una tabla.
power <- new[new$scenario=="all" & new$model %in% c("SEM_MLR","SEM_WLSMV") &
 new$lhs=="SBL" & new$rhs=="Poder",c("model","pvalue","p_Holm_family","Holm_family_size")]
write_json(list(status="CORRECTED_AND_VERIFIED",new_fits=0L,
 maximum_raw_coefficient_SE_p_difference=max_diff,saved_fits_unchanged=TRUE,
 fit_sha256=hashes,checks=checks,primary_Power_SBL=power),
 file.path(out,"auditoria_reparacion_Holm.json"),pretty=TRUE,auto_unbox=TRUE,digits=NA)
print(power,row.names=FALSE)
cat("HOLM EXPORT CORRECTED WITHOUT NEW FITS\n")
