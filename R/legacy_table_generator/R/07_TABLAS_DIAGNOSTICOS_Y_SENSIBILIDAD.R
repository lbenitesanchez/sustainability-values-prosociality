# MÓDULO 07 — TABLAS S10, S11, S12, S13, S14, S15, S16, S21.
# Funciones: formatear_evidencia(); formatear_identificacion(); formatear_probit().
# Todas leen CSV de la reproducción auditada. Ninguna ajusta o prueba modelos nuevos.
# Los nombres S10/S11/S12 dentro de algunos CSV son etiquetas históricas:
# la tabla final se identifica únicamente en guardar(...,"Sxx").

# Tablas S10-S16 y S21: formatear resultados analíticos ya obtenidos.
# La evidencia procede de ajustes guardados y de pruebas nuevas; nunca del Word.
# Las cifras sin redondear, intervalos y advertencias permanecen en los CSV fuente.
REVISION <- Sys.getenv("MA_REVISION_DIR", unset = file.path(getwd(),"revision-observaciones-2026-10-09"))
evidencia <- file.path(REVISION,"evidencia/observaciones_2_3_6_7")
diagnosticos <- file.path(REVISION,"diagnosticos/resultados")
destino <- file.path(REVISION,"tablas_finales")
dir.create(destino,recursive=TRUE,showWarnings=FALSE)
# [AUXILIARES / S10–S16, S21] Leer conserva columnas; guardar evita índices.
# f3/f2 redondean3/2 decimales, fp usa<.001. NA sigue siendo dato no disponible.
# escenario() sólo traduce etiquetas: no vuelve a excluir registros.
leer <- function(p) read.csv(p,check.names=FALSE,stringsAsFactors=FALSE)
guardar <- function(m,id) write.table(m,file.path(destino,paste0(id,".csv")),
  sep=",",quote=TRUE,row.names=FALSE,col.names=FALSE,na="NA",fileEncoding="UTF-8")
f3 <- function(x) ifelse(is.finite(x),sprintf("%.3f",x),"NA")
f2 <- function(x) ifelse(is.finite(x),sprintf("%.2f",x),"NA")
fp <- function(x) ifelse(!is.finite(x),"NA",ifelse(x<.001,"<.001",sprintf("%.3f",x)))
escenario <- function(x) unname(c(all="All",without_duplicate_copies="No matching copies",
  without_constant22="No constant profiles",without_both="Neither flag")[x])

formatear_evidencia <- function() {
  # S12: cuatro escenarios x dos estimadores; AVE y omega son escalas distintas.
  x <- leer(file.path(evidencia,"archivado_S10_medicion_y_SEM_sensibilidad.csv"))
  stopifnot(nrow(x)==8L)
  guardar(rbind(c("Sample","N","Est.","SC–SBL r","PSB AVE","PSB score omega","SEM CFI"),
    cbind(escenario(x$scenario),x$N,x$estimator,f3(x$CFA_SC_SBL_correlation),f3(x$PSB_AVE),
      f3(x$PSB_omega_observed_categories),paste0(f3(ifelse(x$estimator=="MLR",x$SEM_CFI_robust,x$SEM_CFI_scaled)),
        ifelse(x$estimator=="MLR"," (R)"," (S)")))),"S12")
  # S13: intervalos puntuales del 95% y ajuste Holm DENTRO de cada modelo/familia.
  x <- leer(file.path(evidencia,"archivado_S11_asociaciones_PP_sensibilidad.csv"))
  guardar(rbind(c("Sample / N","Est.","Benevolence–PP3 b [95% CI]","PP3 p Holm18","Joint p Holm6","Max SBL p Holm18"),
    cbind(paste0(escenario(x$scenario)," / ",x$N),x$estimator,
      paste0(f3(x$Benevolence_PP3_b)," [",f3(x$Benevolence_PP3_CI_lower),", ",f3(x$Benevolence_PP3_CI_upper),"]"),
      fp(x$Benevolence_PP3_p_Holm18),fp(x$Benevolence_joint_p_Holm6),fp(x$maximum_SBL_PP_p_Holm18))),"S13")
  # S14: diferencias pareadas B2-B3, sin tratar las cinco repeticiones como independientes.
  x <- leer(file.path(evidencia,"archivado_S12_incremento_predictivo_sensibilidad.csv"))
  p <- x[x$estimator=="ordinal_probit",];o <- x[x$estimator=="OLS",]
  o <- o[match(paste(p$scenario,p$item),paste(o$scenario,o$item)),]
  stopifnot(nrow(p)==12L,all(p$N==o$N))
  guardar(rbind(c("Sample","N","Item","Probit RMSE change %","OLS RMSE change %","Probit RPS difference"),
    cbind(escenario(p$scenario),p$N,p$item,f2(p$relative_RMSE_change_percent),
      f2(o$relative_RMSE_change_percent),sprintf("%.5f",p$mean_difference_RPS))),"S14")
# [S21 / FUENTES] Ajuste M6 multigrupo y dos métodos robustos LRT ya
# exportados. ptest no recalcula el contraste; una comparación inexistente
# se muestra con raya. No infiere invariancia M2 ni igualdad país-institución.
  # S21: no inferir igualdad de medición M2 a partir de la sensibilidad M6.
  x <- leer(file.path(evidencia,"S_invariancia_M6_MLR_ajuste_N_df_npar_IC90.csv"))
  tests <- leer(file.path(evidencia,"S_invariancia_M6_MLR_LRT_metodos_explicitos.csv"))
  ptest <- function(model,method) {
    z <- tests[tests$model==model & tests$method==method,,drop=FALSE]
    if(nrow(z)==0L)return("—")
    stopifnot(nrow(z)==1L);fp(z[["Pr(>Chisq)"]])
  }
  guardar(rbind(c("Restriction","df","CFI robust","RMSEA robust [90% CI]","SRMR","p SB2001","p S2000"),
    cbind(x$model,x$df,f3(x$cfi.robust),paste0(f3(x$rmsea.robust)," [",f3(x$rmsea.ci.lower.robust),", ",f3(x$rmsea.ci.upper.robust),"]"),
      f3(x$srmr),vapply(x$model,ptest,"",method="satorra.bentler.2001"),
      vapply(x$model,ptest,"",method="satorra.2000"))),"S21")
  invisible(TRUE)
}

formatear_identificacion <- function() {
  # S10: las restricciones y marcadores proceden de los objetos originales.
  x <- leer(file.path(diagnosticos,"restricciones_7_especificaciones.csv"))
  guardar(rbind(c("Model","Measurement / cross-loadings","SEM scale","Factor relation / rotation"),
    cbind(sub("_.*","",x$model),x$crossloadings,x$SEM_scaling,
      paste(x$covariance_constraint,x$rotation_constraints,sep="; "))),"S10")
  # S11: rango en el espacio tangente de TODAS las restricciones activas.
  # d = número efectivo de direcciones libres; H = Hessiano observado proyectado.
  x <- leer(file.path(diagnosticos,"identificacion_28_ajustes.csv"))
  stopifnot(nrow(x)==28L)
# [S11 / REGLA DE REPORTE] Pass exige rango local, curvatura positiva,
# varianzas propias y errores estándar finitos. Convergencia por sí sola
# no basta; el M6 ordinal señalado permanece visible como Flag.
  estado <- ifelse(x$local_rank_pass & x$tangent_Hessian_positive & x$proper_variance &
    x$finite_free_SE,"Pass", "Flag")
  global <- ifelse(is.finite(x$robust_chisq)&is.finite(x$robust_p),"Available","Unavailable")
  guardar(rbind(c("Model / type","Est.","df / d","Delta / I rank","Min H eigenvalue","Local / variance / SE","Robust test"),
    cbind(paste(sub("_.*","",x$model),x$kind),x$estimator,paste(x$df,x$free_tangent_dimension,sep=" / "),
      paste(x$tangent_Delta_rank,x$tangent_information_rank,sep=" / "),
      formatC(x$tangent_Hessian_min_eigen,digits=2,format="e"),estado,global)),"S11")
  invisible(TRUE)
}

formatear_probit <- function() {
  # S15: soporte observado completo; no se fusiona ninguna categoría.
  x <- leer(file.path(diagnosticos,"frecuencias_categorias.csv"))
  filas <- lapply(paste0("PP",1:3),function(item){
    z <- x[x$outcome==item,];z <- z[order(z$category),]
    stopifnot(identical(z$category,1:5),sum(z$n)==653L)
    c(item,z$n,f2(100*sum(z$n[z$category>=4])/sum(z$n)),sum(z$n))
  })
  guardar(rbind(c("Item","1","2","3","4","5","4–5 %","N"),do.call(rbind,filas)),"S15")
  # S16: pruebas parciales preespecificadas, incluidos los intentos no evaluables.
  x <- leer(file.path(diagnosticos,"pruebas_nominales_valores.csv"))
  stopifnot(nrow(x)==12L,all(x$planned_family_size==12L))
# [S16 / REGLA DE REPORTE] Holm12 usa la familia planificada de12 pruebas.
# Un diagnóstico no evaluable se conserva como tal, sin convertirlo en p=1.
  valores <- c(Benevolencia="Benevolence",Universalismo="Universalism",Poder="Power",Logro="Achievement")
  estado <- ifelse(x$LR_inference_valid,"Evaluable",paste0("Not evaluable: ",x$diagnostic_reason))
  guardar(rbind(c("Item","Value","LR","Δdf","p raw","p Holm12","Diagnostic"),
    cbind(x$outcome,unname(valores[x$nominal_value]),f3(x$LR),x$df_added,fp(x$p_raw),fp(x$p_Holm_12),estado)),"S16")
}

if(sys.nframe()==0L) {
  formatear_evidencia()
  formatear_identificacion()
  formatear_probit()
}
