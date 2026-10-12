# AVISO V3: las tablas finales de este módulo son S5, S6, S7, S8.
# T3/S9 de este archivo pertenecen al contrato histórico y se calculan sólo
# para verificar consistencia. El principal las descarta explícitamente y
# usa las tablas actuales de validación AGRUPADA del módulo06.
# No hay entrenamiento aquí: se extraen modelos y predicciones ya guardados.
# Los RDS privados son requisitos locales; nunca se incluyen en la entrega pública.

# Módulo 02: reproducción de las tablas T3 y S5–S9 --------------------------
# Este archivo define una función; Source no ajusta modelos ni cambia la base.
# El controlador llama generar_tablas_items(config), donde config contiene:
#   proyecto: carpeta original Paper con MA - 1;
#   salida: carpeta nueva para CSV agregados;
#   registro: carpeta del contrato de tablas y de comprobaciones.
# Se releen objetos ya estimados y predicciones ya auditadas. Se recalculan
# estadísticos, métricas y tablas; esto NO constituye una nueva estimación.
# Las etiquetas [T3], [S5] ... aparecen como secciones en el esquema de RStudio.

generar_tablas_items <- function(config) {
  stopifnot(is.list(config), all(c("proyecto", "salida", "registro") %in% names(config)))
  for (key in c("proyecto", "salida", "registro")) {
    stopifnot(is.character(config[[key]]), length(config[[key]]) == 1L,
              !is.na(config[[key]]), nzchar(config[[key]]))
  }
  proyecto <- normalizePath(config$proyecto, mustWork = TRUE)
  salida <- config$salida; registro <- config$registro
  dir.create(salida, recursive = TRUE, showWarnings = FALSE)
  dir.create(registro, recursive = TRUE, showWarnings = FALSE)
  # Las bibliotecas del proyecto se priorizan, sin instalar paquetes.
  librerias_previas <- .libPaths()
  on.exit(.libPaths(librerias_previas), add = TRUE)
  .libPaths(c(file.path(proyecto, "revision-revista-2026-10-04/auditoria-codigos/R-library"),
              path.expand("~/Library/R/arm64/4.6/library"), librerias_previas))
  paquetes <- c("lavaan", "semTools", "jsonlite", "digest")
  faltan <- paquetes[!vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)]
  if (length(faltan)) stop("Faltan bibliotecas locales: ", paste(faltan, collapse = ", "))

  raiz <- file.path(proyecto, "sem-items-ordinal-2026-10-08")
# [FUENTES / S5–S8 Y CONTRATO HISTÓRICO T3/S9] CFA M2, SEM por ítem,
# predicciones y asignaciones existentes. Todos se protegen con SHA256.
  fuentes <- c(
    CFA_MLR = file.path(proyecto, "sem-alternativas-2026-10-08/resultados/M2_puentes2_CFA_MLR.rds"),
    CFA_WLSMV = file.path(proyecto, "sem-alternativas-2026-10-08/resultados/M2_puentes2_CFA_WLSMV.rds"),
    items_MLR = file.path(raiz, "items/privados/SEM_items_MLR.rds"),
    items_WLSMV = file.path(raiz, "items/privados/SEM_items_WLSMV.rds"),
    datos_P3 = file.path(raiz, "items/privados/datos_P3_para_items_y_CV.rds"),
    predicciones = file.path(raiz, "prediccion/resultados/privados/predicciones_fuera_muestra.rds"),
    particiones = file.path(raiz, "prediccion/resultados/privados/asignaciones_cv.rds"))
  stopifnot(all(file.exists(fuentes)))
  hash <- function(p) digest::digest(p, algo = "sha256", file = TRUE)
  hashes_antes <- vapply(fuentes, hash, character(1))
  pp <- c("PP1", "PP2", "PP3")
  predictores <- c("SC", "SBL", "Benevolencia", "Universalismo", "Poder", "Logro")
  etiquetas <- c(SC = "SC", SBL = "SBL", Benevolencia = "Benevolence",
                 Universalismo = "Universalism", Poder = "Power", Logro = "Achievement")
  tablas <- list()
  # Los números se formatean al FINAL, para que los CSV mantengan precisión.
# [FORMATO / S5–S8] Decimales y p<.001 se aplican al texto final.
# una_fila() exige clave única; evita elegir silenciosamente una fila repetida.
  fnum <- function(x, decimales = 3L) formatC(as.numeric(x), format = "f", digits = decimales, decimal.mark = ".")
  fp <- function(x) ifelse(x < .001, "<.001", fnum(x))
  exportar <- function(x, nombre) utils::write.csv(x, file.path(salida, nombre), row.names = FALSE, na = "NA")
  una_fila <- function(x, ...) {
    claves <- list(...); mantener <- rep(TRUE, nrow(x))
    for (key in names(claves)) mantener <- mantener & x[[key]] == claves[[key]]
    respuesta <- x[mantener, , drop = FALSE]
    stopifnot(nrow(respuesta) == 1L)
    respuesta
  }

  # [S5] CFA M2 ordinal, AVE y fiabilidad en sus respectivas escalas ----------
  # Sólo se utilizan los CFA sin ecuaciones estructurales. En SEM ordinal
  # conditional.x, semTools puede encontrar res.cov en lugar de cov; ese
  # objeto NO se introduce aquí ni se interpreta su varianza residual como total.
  # AVE(obs.var=TRUE) usa varianzas de referencia. ord.scale=TRUE proporciona
  # fiabilidad de la suma de categorías observadas; FALSE usa respuesta latente.
  extraer_PSB <- function(x) {
    y <- if (is.list(x)) x[["PSB"]] else x["PSB"]
    stopifnot(length(y) == 1L, is.finite(y))
    as.numeric(y)
  }
  medicion <- do.call(rbind, lapply(c("CFA_MLR", "CFA_WLSMV"), function(id) {
    ajuste <- readRDS(fuentes[[id]])
    opciones <- lavaan::lavInspect(ajuste, "options")
    stopifnot(lavaan::lavInspect(ajuste, "converged"),
              isTRUE(lavaan::lavInspect(ajuste, "post.check")), !isTRUE(opciones$conditional.x))
    indices <- lavaan::fitMeasures(ajuste, c("cfi.robust", "cfi.scaled", "rmsea.robust", "rmsea.scaled"))
    ave <- extraer_PSB(semTools::AVE(ajuste, obs.var = TRUE))
    omega_puntaje <- extraer_PSB(semTools::compRelSEM(ajuste, obs.var = TRUE, ord.scale = TRUE))
    omega_latente <- extraer_PSB(semTools::compRelSEM(ajuste, obs.var = TRUE, ord.scale = FALSE))
    # Comprobación adicional del AVE: varianza explicada / varianza de referencia.
    matrices <- lavaan::lavInspect(ajuste, "est")
    lambda <- matrices$lambda[pp, "PSB"]
    referencia <- lavaan::lavInspect(ajuste, "sampstat")[["cov", exact = TRUE]]
    if (isTRUE(lavaan::lavInspect(ajuste, "categorical"))) {
      # sampstat contiene correlaciones de respuestas ordinales. Las cargas
      # theta usan su escala sin normalizar: se restaura con inversos de delta,
      # igual que AVE de semTools, antes de comparar varianzas explicadas.
      stopifnot(identical(rownames(matrices$delta), colnames(referencia)))
      inversos_delta <- diag(1 / as.numeric(matrices$delta))
      dimnames(inversos_delta) <- dimnames(referencia)
      referencia <- inversos_delta %*% referencia %*% inversos_delta
    }
    ave_matricial <- sum(lambda^2 * matrices$psi["PSB", "PSB"]) / sum(diag(referencia)[match(pp, colnames(referencia))])
    stopifnot(abs(ave - ave_matricial) < 1e-8)
    data.frame(model = id, N = lavaan::lavInspect(ajuste, "nobs"),
      response_scale = if (id == "CFA_WLSMV") "underlying ordinal response" else "observed numeric score",
      cfi_robust = as.numeric(indices["cfi.robust"]), cfi_scaled = as.numeric(indices["cfi.scaled"]),
      rmsea_robust = as.numeric(indices["rmsea.robust"]), rmsea_scaled = as.numeric(indices["rmsea.scaled"]),
      cor_SC_SBL = lavaan::lavInspect(ajuste, "cor.lv")["SC", "SBL"],
      AVE_PSB_reference = ave, AVE_PSB_matrix_reference = ave_matricial,
      omega_PSB_observed_score_reference = omega_puntaje,
      omega_PSB_latent_response_reference = omega_latente, stringsAsFactors = FALSE)
  }))
  exportar(medicion, "02_S5_medicion_numerica.csv")
  filas_s5 <- lapply(list(c("CFA_MLR", "robust", "MLR robust"),
                         c("CFA_WLSMV", "scaled", "WLSMV scaled"),
                         c("CFA_WLSMV", "robust", "WLSMV robust")), function(especificacion) {
    x <- una_fila(medicion, model = especificacion[1L])
    c(especificacion[3L], fnum(x[[paste0("cfi_", especificacion[2L])]]),
      fnum(x[[paste0("rmsea_", especificacion[2L])]]), fnum(x$cor_SC_SBL),
      fnum(x$AVE_PSB_reference), fnum(x$omega_PSB_observed_score_reference))
  })
  tablas$S5 <- rbind(c("Analysis", "CFI", "RMSEA", "SC SBL r", "AVE PSB", "Score ω"), do.call(rbind, filas_s5))

  # [S6] Wald conjunto y comparación MLR entre respuestas --------------------
  # Para cada predictor se prueba conjuntamente b_PP1=b_PP2=b_PP3=0 (3 df).
  # La prueba MLR de igualdad usa dos diferencias: b_PP1-b_PP2 y b_PP1-b_PP3.
  # Se aplica Holm a seis predictores dentro de cada familia/estimador.
  # V contiene la covarianza ROBUSTA del ajuste ya estimado. No se ajusta otro SEM.
  wald <- function(b, V, H) {
    diferencias <- as.vector(H %*% b)
    var_diferencias <- H %*% V %*% t(H)
    stopifnot(all(is.finite(var_diferencias)), min(eigen(var_diferencias, symmetric = TRUE, only.values = TRUE)$values) > 0)
    estadistico <- as.numeric(crossprod(diferencias, solve(var_diferencias, diferencias)))
    data.frame(statistic = estadistico, df = nrow(H), pvalue = stats::pchisq(estadistico, nrow(H), lower.tail = FALSE))
  }
# [CÁLCULO COMPARTIDO / S6, S7, S8] Se extraen los18 caminos (6 predictores
# × 3 PP). La covarianza robusta guardada produce Wald; no se estima otro SEM.
  caminos_lista <- conjuntos_lista <- igualdad_lista <- list()
  for (id in c("MLR", "WLSMV")) {
    ajuste <- readRDS(fuentes[[paste0("items_", id)]])
    stopifnot(lavaan::lavInspect(ajuste, "converged"), isTRUE(lavaan::lavInspect(ajuste, "post.check")))
    parametros <- lavaan::parameterEstimates(ajuste, ci = TRUE)
    estandarizados <- lavaan::standardizedSolution(ajuste, type = "std.all", ci = TRUE)
    tabla_parametros <- lavaan::parTable(ajuste)
    V <- lavaan::lavInspect(ajuste, "vcov")
    caminos <- parametros[parametros$op == "~" & parametros$lhs %in% pp & parametros$rhs %in% predictores, ]
    stopifnot(nrow(caminos) == 18L)
    claves <- paste(caminos$lhs, caminos$op, caminos$rhs)
    estandar <- estandarizados[match(claves, paste(estandarizados$lhs, estandarizados$op, estandarizados$rhs)), ]
    caminos$estimator <- id
    caminos$beta_std_all <- estandar$est.std
    caminos$beta_CI_lower <- estandar$ci.lower
    caminos$beta_CI_upper <- estandar$ci.upper
    # Holm utiliza p de coeficientes NO estandarizados; los IC β son puntuales.
    caminos$p_Holm_18 <- stats::p.adjust(caminos$pvalue, method = "holm")
    stopifnot(all(is.finite(as.matrix(caminos[c("est", "se", "pvalue", "beta_std_all", "beta_CI_lower", "beta_CI_upper")]))))
    caminos_lista[[id]] <- caminos
    conjunto <- igualdad <- list()
    for (predictor in predictores) {
      seleccion <- tabla_parametros[tabla_parametros$op == "~" & tabla_parametros$lhs %in% pp & tabla_parametros$rhs == predictor, ]
      seleccion <- seleccion[match(pp, seleccion$lhs), ]
      libres <- seleccion$free
      stopifnot(length(libres) == 3L, !anyNA(libres), all(libres > 0), !anyDuplicated(libres))
      b <- seleccion$est; covarianza <- V[libres, libres, drop = FALSE]
      conjunto[[predictor]] <- cbind(data.frame(estimator = id, predictor = predictor), wald(b, covarianza, diag(3L)))
      if (id == "MLR") {
        H <- rbind(c(1, -1, 0), c(1, 0, -1))
        igualdad[[predictor]] <- cbind(data.frame(estimator = id, predictor = predictor), wald(b, covarianza, H))
      }
    }
    conjunto <- do.call(rbind, conjunto); rownames(conjunto) <- NULL
    conjunto$p_Holm_6 <- stats::p.adjust(conjunto$pvalue, method = "holm")
    conjuntos_lista[[id]] <- conjunto
    if (id == "MLR") {
      igualdad <- do.call(rbind, igualdad); rownames(igualdad) <- NULL
      igualdad$p_Holm_6 <- stats::p.adjust(igualdad$pvalue, method = "holm")
      igualdad_lista[[id]] <- igualdad
    }
  }
  caminos <- do.call(rbind, caminos_lista); conjuntos <- do.call(rbind, conjuntos_lista)
  igualdad <- do.call(rbind, igualdad_lista)
  exportar(caminos, "02_S7_S8_caminos_numericos_IC_Holm18.csv")
  exportar(conjuntos, "02_S6_Wald_conjunto_numerico_Holm6.csv")
  exportar(igualdad, "02_S6_Wald_igualdad_MLR_numerico_Holm6.csv")
  tablas$S6 <- rbind(c("Predictor", "MLR p Holm6", "Ordinal p Holm6", "MLR equality p Holm6"),
    do.call(rbind, lapply(predictores, function(predictor) c(etiquetas[[predictor]],
      fp(una_fila(conjuntos, estimator = "MLR", predictor = predictor)$p_Holm_6),
      fp(una_fila(conjuntos, estimator = "WLSMV", predictor = predictor)$p_Holm_6),
      fp(una_fila(igualdad, estimator = "MLR", predictor = predictor)$p_Holm_6)))))
# [FORMATO / S7, S8] La misma función forma las dos tablas cambiando sólo
# el estimador. Estrella=Holm18; intervalo=IC95% puntual estandarizado.
  crear_tabla_caminos <- function(estimador) {
    rbind(c("Predictor", "PP1 β [95% CI]", "PP2 β [95% CI]", "PP3 β [95% CI]"),
      do.call(rbind, lapply(predictores, function(predictor) {
        c(etiquetas[[predictor]], vapply(pp, function(item) {
          x <- una_fila(caminos, estimator = estimador, lhs = item, rhs = predictor)
          paste0(fnum(x$beta_std_all), if (x$p_Holm_18 < .05) "*" else "", "\n[",
                 fnum(x$beta_CI_lower), ", ", fnum(x$beta_CI_upper), "]")
        }, character(1)))
      })))
  }
  # [S7] Asociaciones por respuesta MLR -------------------------------------
  # β describe la asociación estandarizada; * indica p_Holm18<.05.
  # Los intervalos robustos del 95% NO están ajustados por multiplicidad.
  tablas$S7 <- crear_tabla_caminos("MLR")
  # [S8] Asociaciones por respuesta ordinal ---------------------------------
  # La forma es idéntica, pero la respuesta probit usa una escala latente:
  # no se confunde una diferencia de significación con diferencia entre ítems.
  tablas$S8 <- crear_tabla_caminos("WLSMV")

  # Predicciones privadas existentes: recalcular métricas, no ajustar --------
# [LEGADO INACTIVO PARA TABLAS FINALES T3/S9] Esta extracción conserva
# la auditoría del contrato histórico. No determina T3/S9 actuales.
# La ejecución actual pasa después al módulo06 y usa particiones agrupadas.
  predicciones <- readRDS(fuentes[["predicciones"]])
  datos <- as.data.frame(readRDS(fuentes[["datos_P3"]]))
  asignaciones <- readRDS(fuentes[["particiones"]])
  categorias <- seq_len(5L); repeticiones <- seq_len(5L); pliegues <- seq_len(5L)
  stopifnot(nrow(datos) == 653L, nrow(predicciones) == 88155L, length(asignaciones) == 5L,
            all(predicciones$item %in% pp), all(predicciones$observed %in% categorias))
  respuestas <- as.matrix(datos[pp])
  stopifnot(all(predicciones$observed == respuestas[cbind(predicciones$row_id, match(predicciones$item, pp))]))
  for (r in repeticiones) {
    x <- predicciones[predicciones$repetition == r, ]
    stopifnot(all(x$fold == asignaciones[[r]][x$row_id]))
  }
  # Cada grupo es una sola prueba fuera de muestra del mismo ítem/modelo/fold.
  grupos <- split(predicciones, paste(predicciones$repetition, predicciones$fold,
                                    predicciones$item, predicciones$model, predicciones$estimator, sep = "|"))
  metricas_fold <- do.call(rbind, lapply(grupos, function(x) {
    stopifnot(!anyDuplicated(x$row_id), all(is.finite(x$predicted_expected)))
    claves <- x[1L, c("repetition", "fold", "item", "model", "estimator")]
    esperado <- x$predicted_expected
    suma_RPS <- suma_logloss_bruto <- suma_logloss_piso <- NA_real_
    n_ceros <- NA_integer_
    if (claves$estimator != "OLS") {
      P <- as.matrix(x[paste0("probability_", categorias)])
      stopifnot(all(is.finite(P)), all(P >= 0), all(P <= 1), max(abs(rowSums(P) - 1)) < 1e-10)
      esperado_prob <- as.vector(P %*% categorias)
      stopifnot(max(abs(esperado_prob - esperado)) < 1e-12)
      # RPS compara las cuatro probabilidades acumuladas con lo ocurrido.
      acumulada <- cbind(P[, 1L], rowSums(P[, 1:2, drop = FALSE]),
                        rowSums(P[, 1:3, drop = FALSE]), rowSums(P[, 1:4, drop = FALSE]))
      real <- outer(x$observed, 1:4, "<=") * 1
      suma_RPS <- sum((acumulada - real)^2) / 4
      p_real <- P[cbind(seq_len(nrow(x)), match(x$observed, categorias))]
      # El logloss bruto puede ser Inf por categorías ausentes. El piso fijo
      # es sólo evaluación numérica; no transforma ni regulariza los modelos.
      suma_logloss_bruto <- sum(-log(p_real)); suma_logloss_piso <- sum(-log(pmax(p_real, 1e-12)))
      n_ceros <- sum(p_real == 0)
    } else stopifnot(all(is.na(as.matrix(x[paste0("probability_", categorias)]))))
    error <- esperado - x$observed
    cbind(claves, data.frame(n_test = nrow(x), sum_squared_error = sum(error^2),
      sum_absolute_error = sum(abs(error)), sum_RPS = suma_RPS,
      sum_logloss_raw = suma_logloss_bruto, sum_logloss_floor = suma_logloss_piso,
      n_probability_zero = n_ceros))
  }))
  rownames(metricas_fold) <- NULL
  stopifnot(nrow(metricas_fold) == 675L)
  grupos_rep <- split(metricas_fold, paste(metricas_fold$repetition, metricas_fold$item,
                                         metricas_fold$model, metricas_fold$estimator, sep = "|"))
  metricas_rep <- do.call(rbind, lapply(grupos_rep, function(x) {
    stopifnot(nrow(x) == 5L, setequal(x$fold, pliegues), sum(x$n_test) == nrow(datos))
    clave <- x[1L, c("repetition", "item", "model", "estimator")]
    n <- sum(x$n_test)
    # Se suman los errores antes de calcular RMSE; NO se promedian RMSE de folds.
    cbind(clave, data.frame(n_test = n, RMSE = sqrt(sum(x$sum_squared_error) / n),
      MAE = sum(x$sum_absolute_error) / n, RPS = sum(x$sum_RPS) / n,
      logloss_raw = sum(x$sum_logloss_raw) / n, logloss_floor = sum(x$sum_logloss_floor) / n,
      n_probability_zero = sum(x$n_probability_zero)))
  }))
  rownames(metricas_rep) <- NULL
  grupos_modelo <- split(metricas_rep, paste(metricas_rep$item, metricas_rep$model, metricas_rep$estimator, sep = "|"))
  metricas <- do.call(rbind, lapply(grupos_modelo, function(x) {
    stopifnot(nrow(x) == 5L, setequal(x$repetition, repeticiones))
    cbind(x[1L, c("item", "model", "estimator")], data.frame(RMSE = mean(x$RMSE),
      MAE = mean(x$MAE), RPS = mean(x$RPS), logloss_raw = mean(x$logloss_raw),
      logloss_floor = mean(x$logloss_floor), RMSE_partition_SD = stats::sd(x$RMSE),
      RPS_partition_SD = stats::sd(x$RPS)))
  }))
  rownames(metricas) <- NULL
  exportar(metricas_fold, "02_T3_S9_metricas_fold_recalculadas.csv")
  exportar(metricas_rep, "02_T3_S9_metricas_repeticion_recalculadas.csv")
  exportar(metricas, "02_T3_S9_metricas_modelo_recalculadas.csv")

  # [T3] Predicción por ítem: cuatro bloques ordinales -----------------------
  # Cada celda es la media de RMSE de CINCO repeticiones, con idénticos folds.
  # B0=exposición/controles; B1 añade valores; B3 añade APS sin valores;
  # B2 contiene exposición, controles, cuatro valores y los 19 APS.
  bloques_t3 <- c("B0_CSR_controls", "B1_add_values", "B3_APS_without_values", "B2_add_APS")
  tablas$T3 <- rbind(c("Response", "B0 RMSE", "B1 RMSE", "B3 RMSE", "B2 RMSE"),
    do.call(rbind, lapply(pp, function(item) c(item, vapply(bloques_t3, function(bloque)
      fnum(una_fila(metricas, item = item, model = bloque, estimator = "ordinal_probit")$RMSE, 4L), character(1))))))

  # [S9] Aporte del BLOQUE de cuatro valores después de APS ------------------
  # El porcentaje se calcula dentro de cada repetición: 100*(RMSE_B2/RMSE_B3-1).
  # Después se promedian esos cinco porcentajes. No se usa el cociente de
  # promedios ni se atribuye el cambio a Benevolencia por separado.
  # Valores negativos son mejora. Los folds se superponen entre repeticiones:
  # no se calculan p ni IC tratando estas particiones como muestras independientes.
  cambios <- do.call(rbind, lapply(pp, function(item) do.call(rbind,
    lapply(c("ordinal_probit", "OLS"), function(estimador) do.call(rbind,
      lapply(repeticiones, function(r) {
        b2 <- una_fila(metricas_rep, repetition = r, item = item, model = "B2_add_APS", estimator = estimador)
        b3 <- una_fila(metricas_rep, repetition = r, item = item, model = "B3_APS_without_values", estimator = estimador)
        data.frame(item = item, estimator = estimador, repetition = r,
          relative_RMSE_change_percent = 100 * (b2$RMSE / b3$RMSE - 1),
          difference_RMSE = b2$RMSE - b3$RMSE, difference_RPS = b2$RPS - b3$RPS)
      }))))))
  cambio_medio <- do.call(rbind, lapply(split(cambios, paste(cambios$item, cambios$estimator)), function(x)
    data.frame(item = x$item[1L], estimator = x$estimator[1L],
      mean_relative_RMSE_change_percent = mean(x$relative_RMSE_change_percent),
      mean_difference_RMSE = mean(x$difference_RMSE), mean_difference_RPS = mean(x$difference_RPS),
      repetitions_better_RMSE = sum(x$difference_RMSE < 0))))
  rownames(cambio_medio) <- NULL
  exportar(cambios, "02_S9_cambios_pareados_recalculados.csv")
  exportar(cambio_medio, "02_S9_cambios_medios_recalculados.csv")
  tablas$S9 <- rbind(c("Item", "RPS B3", "RPS B2", "Ordinal RMSE change %", "OLS RMSE change %"),
    do.call(rbind, lapply(pp, function(item) c(item,
      fnum(una_fila(metricas, item = item, model = "B3_APS_without_values", estimator = "ordinal_probit")$RPS, 5L),
      fnum(una_fila(metricas, item = item, model = "B2_add_APS", estimator = "ordinal_probit")$RPS, 5L),
      fnum(una_fila(cambio_medio, item = item, estimator = "ordinal_probit")$mean_relative_RMSE_change_percent, 2L),
      fnum(una_fila(cambio_medio, item = item, estimator = "OLS")$mean_relative_RMSE_change_percent, 2L)))))

  # Contrato Word: comparar resultados calculados, nunca copiar cifras -------
  # El JSON sólo se lee aquí, cuando TODAS las tablas ya fueron recalculadas.
  contrato_path <- file.path(registro, "tablas_word_esperadas.json")
  stopifnot(file.exists(contrato_path))
  contrato <- jsonlite::read_json(contrato_path, simplifyVector = FALSE)
  claves_tablas <- c("T3", "S5", "S6", "S7", "S8", "S9")
  tablas <- tablas[claves_tablas]
  cotejo <- lapply(claves_tablas, function(key) {
    entrada <- Filter(function(x) identical(x$key, key), contrato)
    stopifnot(length(entrada) == 1L)
    esperado <- do.call(rbind, lapply(entrada[[1L]]$cells, function(x) unlist(x, use.names = FALSE)))
    igual <- identical(unname(tablas[[key]]), unname(esperado))
    utils::write.table(tablas[[key]], file.path(salida, paste0("02_", key, "_celdas_para_Word.csv")),
                       sep = ",", quote = TRUE, qmethod = "double", row.names = FALSE, col.names = FALSE, fileEncoding = "UTF-8")
    list(table = key, matches_Word_contract = igual, dimensions = dim(tablas[[key]]),
         source_document = entrada[[1L]]$document, source_docx_sha256 = entrada[[1L]]$docx_sha256)
  })
  hashes_despues <- vapply(fuentes, hash, character(1))
  sin_cambios <- identical(hashes_antes, hashes_despues)
  mapa <- data.frame(table = claves_tablas, RStudio_section = paste0("[", claves_tablas, "]"),
    function_name = "generar_tablas_items(config)",
    numeric_output = c("02_T3_S9_metricas_modelo_recalculadas.csv", "02_S5_medicion_numerica.csv",
      "02_S6_Wald_conjunto_numerico_Holm6.csv; 02_S6_Wald_igualdad_MLR_numerico_Holm6.csv",
      "02_S7_S8_caminos_numericos_IC_Holm18.csv", "02_S7_S8_caminos_numericos_IC_Holm18.csv", "02_S9_cambios_medios_recalculados.csv"),
    calculation = c("RMSE por repetición desde predicciones privadas; media de cinco repeticiones",
      "CFA M2 existente: fitMeasures, AVE y compRelSEM; escalas observada y latente",
      "Wald con covarianza robusta: conjuntos por estimador e igualdad sólo MLR; Holm6",
      "standardizedSolution MLR e IC robusto; estrella por Holm18 del p no estandarizado",
      "standardizedSolution WLSMV e IC robusto; estrella por Holm18 del p no estandarizado",
      "RPS ordinal y cambio relativo de RMSE B2/B3 por repetición; ambos estimadores"), stringsAsFactors = FALSE)
  utils::write.csv(mapa, file.path(registro, "02_mapa_tablas_y_codigo.csv"), row.names = FALSE)
  comprobacion <- list(status = "REPRODUCED_FROM_VERIFIED_EXISTING_OBJECTS", new_model_fits = 0L,
    original_sources_unchanged = sin_cambios, source_paths = fuentes,
    sha256_before = hashes_antes, sha256_after = hashes_despues,
    tables = cotejo, all_Word_cells_match = all(vapply(cotejo, function(x) x$matches_Word_contract, logical(1))),
    N = nrow(datos), existing_private_predictions = nrow(predicciones), recomputed_fold_groups = nrow(metricas_fold),
    rare_category_caveat = "Raw logloss remains infinite for PP2/PP3; fixed 1e-12 scoring floor is diagnostic only. RPS/RMSE/MAE finite.",
    scope = "Secondary exploratory internal reproduction; no new data, estimates, variable selection or p-values from overlapping folds.")
  jsonlite::write_json(comprobacion, file.path(registro, "02_comprobacion_reproduccion.json"),
                       pretty = TRUE, auto_unbox = TRUE, digits = NA)
  stopifnot(sin_cambios, comprobacion$all_Word_cells_match)
  message("Módulo 02: T3 y S5–S9 reproducidas desde objetos existentes; ningún ajuste nuevo.")
  tablas
}
