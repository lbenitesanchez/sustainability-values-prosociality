# =============================================================================
# DATOS INTEGRAMENTE ARTIFICIALES PARA DEMOSTRAR EL ESQUEMA ANALITICO DE V3
# =============================================================================
# Autoría del artefacto: preparado para revisión por los autores, 2026-10-10.
# Uso: Rscript --vanilla 02_SIMULAR_DATOS_ARTIFICIALES.R [carpeta_de_salida]
# Dependencias obligatorias: solamente R base y los paquetes distribuidos con R.
# Dependencia opcional: haven, exclusivamente para exportar una copia .sav.
#
# ADVERTENCIA CIENTIFICA:
# Este script NO reproduce las cifras del manuscrito o del suplemento. Genera
# un ejemplo artificial, sin leer datos de participantes, modelos estimados,
# correlaciones, frecuencias ni exportaciones de resultados del estudio.
# Los parámetros de abajo son HIPOTETICOS. No han sido calibrados a microdatos.
# No es una base anonimizada, una imputación, una réplica de casos reales ni un
# mecanismo con garantía formal de privacidad diferencial. Los nombres/códigos
# de variables se conservan para enseñar el diseño y las recodificaciones.
#
# RELACION CON LAS TABLAS:
# - Las variables permiten ilustrar los bloques de predictores y resultados
#   utilizados por las tablas del manuscrito y del suplemento.
# - Ninguna tabla empírica del artículo se obtiene a partir de esta base.
# - Para las cifras publicadas use 01_GENERAR_TODAS_LAS_TABLAS_V3.R y sus
#   fuentes verificadas. Los modelos reales exigen datos autorizados.
# - Los scripts originales fijan N, hashes y países reales. Esta base NO debe
#   introducirse en ese flujo como si fuese la muestra original.
#
# CODIFICACION ESPECIAL DE SEDE:
# Pais = 1, 2, 3, 4 significa Sede_Ficticia_A, B, C, D. No son países reales.
# C_Colombia, C_Peru y C_Mexico son nombres HEREDADOS del esquema de regresión;
# aquí solo identifican los códigos ficticios 2, 3 y 4, respectivamente.
# =============================================================================

parametros_artificiales <- function() {
  # Una muestra equilibrada ilustrativa; N=800 no procede de la muestra real.
  # No cambie estos valores para imitar el ajuste o significación del artículo.
  list(
    seed = 20261010L,
    n_por_sede = 200L,
    sedes = paste0("Sede_Ficticia_", LETTERS[1:4]),
    prob_mujer = 0.52,
    prob_postgrado = 0.25,
    edad_media = 27,
    edad_sd = 8,
    edad_min = 18L,
    edad_max = 60L,
    prob_nse = c(0.10, 0.22, 0.36, 0.24, 0.08),
    csr_intercepto = -0.15,
    csr_efecto_postgrado = 0.35,
    csr_efecto_sede = c(-0.20, 0.10, 0.05, 0.00),
    # DGP continuo meramente didáctico: mezcla de ruido normal y dependencias
    # moderadas. Los nombres de componentes no implican validez psicométrica.
    st_efecto_csr = 0.20,
    se_efecto_csr = 0.05,
    sc_efecto_st = 0.30,
    sc_efecto_csr = 0.15,
    sbl_efecto_sc = 0.35,
    sbl_efecto_st = 0.20,
    psb_efecto_st = 0.25,
    psb_efecto_se = -0.10,
    psb_efecto_sc = 0.15,
    psb_efecto_sbl = 0.20,
    psb_efecto_csr = 0.10,
    carga_aps = 0.65,
    carga_pp = c(0.70, 0.65, 0.60),
    carga_valores = c(0.70, 0.65, 0.70, 0.65),
    cortes_aps = c(-1.20, -0.35, 0.40, 1.20),
    cortes_pp = c(-1.10, -0.25, 0.50, 1.35),
    cortes_valores = c(-1.30, -0.65, -0.05, 0.50, 1.20)
  )
}

esquema_artificial <- function() {
  # Los nombres y dominios se verificaron contra el código original y el
  # diccionario SPSS. No se importan esos archivos al ejecutar el generador.
  list(
    sc = paste0("APS", c(1:7, 9, 10, 14)),
    sbl = paste0("APS", c(11:13, 15, 16, 18:21)),
    pp = paste0("PP", 1:3),
    valores = c("Benevolencia", "Universalismo", "Poder", "Logro")
  )
}

ordinal_artificial <- function(componente, carga, cortes) {
  # Variable auxiliar continua = carga * componente + error nuevo independiente.
  # La discretización usa intervalos abiertos a la izquierda y cerrados a la
  # derecha. Cada llamada genera error nuevo. No extrae filas de ninguna base.
  ruido <- stats::rnorm(length(componente))
  z <- carga * componente + sqrt(1 - carga^2) * ruido
  as.integer(cut(z, breaks = c(-Inf, cortes, Inf), labels = FALSE,
                 right = TRUE, include.lowest = TRUE))
}

simular_base_artificial <- function(p = parametros_artificiales()) {
  # Fijar todos los componentes del RNG permite volver a obtener la misma base
  # con R moderno. El registro de sesión identifica la versión efectivamente usada.
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion",
          sample.kind = "Rejection")
  set.seed(p$seed)
  s <- esquema_artificial()
  n <- as.integer(p$n_por_sede * length(p$sedes))
  pais <- rep(seq_along(p$sedes), each = p$n_por_sede)

  # Demografía artificial y nombres explícitos que evitan confundir estos casos
  # con participantes. El identificador solo enumera casos del ejemplo.
  d <- data.frame(
    Tipo_dato = rep("SINTETICO_HIPOTETICO_NO_PARTICIPANTES", n),
    ID_sintetico = sprintf("SYN_%04d", seq_len(n)),
    Pais = as.integer(pais),
    Sede_ficticia = p$sedes[pais],
    Sexo = as.integer(stats::rbinom(n, 1, p$prob_mujer)),
    Estudios = as.integer(stats::rbinom(n, 1, p$prob_postgrado)),
    Edad = as.integer(pmax(p$edad_min, pmin(p$edad_max,
                        round(stats::rnorm(n, p$edad_media, p$edad_sd))))),
    NSE = as.integer(sample.int(5L, n, replace = TRUE, prob = p$prob_nse)),
    stringsAsFactors = FALSE
  )
  # El acotamiento de Edad es parte del ejemplo, no un tratamiento de valores
  # extremos aplicado al estudio. No se introducen ni imputan valores perdidos.

  # CSR1_rec y CSR1_bin: 0=No / 1=Sí.
  # CSR1 original: 1=Sí / 2=No. Restar 1 a CSR1 invertiría su significado.
  prob_csr <- stats::plogis(p$csr_intercepto +
                           p$csr_efecto_postgrado * d$Estudios +
                           p$csr_efecto_sede[d$Pais])
  d$CSR1_rec <- as.integer(stats::rbinom(n, 1, prob_csr))
  d$CSR1 <- ifelse(d$CSR1_rec == 1L, 1L, 2L)
  d$CSR1_bin <- as.integer(d$CSR1_rec)

  # Componentes latentes didácticos; no se exportan como estimaciones ni como
  # puntuaciones de factores. Toda dependencia procede de estas ecuaciones.
  st <- p$st_efecto_csr * d$CSR1_bin + stats::rnorm(n)
  se <- p$se_efecto_csr * d$CSR1_bin + stats::rnorm(n)
  sc <- p$sc_efecto_st * st + p$sc_efecto_csr * d$CSR1_bin + stats::rnorm(n)
  sbl <- p$sbl_efecto_sc * sc + p$sbl_efecto_st * st + stats::rnorm(n)
  psb <- p$psb_efecto_st * st + p$psb_efecto_se * se +
         p$psb_efecto_sc * sc + p$psb_efecto_sbl * sbl +
         p$psb_efecto_csr * d$CSR1_bin + stats::rnorm(n)

  # 19 ítems APS: mismas variables retenidas, respuesta ordinal de 1 a 5.
  for (v in s$sc) d[[v]] <- ordinal_artificial(sc, p$carga_aps, p$cortes_aps)
  for (v in s$sbl) d[[v]] <- ordinal_artificial(sbl, p$carga_aps, p$cortes_aps)
  # PP1–PP3: los TRES ítems que forman el resultado prosocial, de 1 a 5.
  for (j in seq_along(s$pp)) {
    d[[s$pp[j]]] <- ordinal_artificial(psb, p$carga_pp[j], p$cortes_pp)
  }
  # Cuatro valores observados de 1 a 6; no se generan otros ítems ni quijotismo.
  for (j in seq_along(s$valores)) {
    comp <- if (j <= 2L) st else se
    d[[s$valores[j]]] <- ordinal_artificial(comp, p$carga_valores[j],
                                            p$cortes_valores)
  }

  # Recodificaciones y controles reproducen las reglas verificadas del esquema.
  d$Sexo_num <- as.numeric(d$Sexo)
  d$Estudios_num <- as.numeric(d$Estudios)
  d$Edad_z <- as.numeric(scale(d$Edad))
  d$NSE_z <- as.numeric(scale(d$NSE))
  d$Edad_original <- as.numeric(d$Edad)
  d$NSE_original <- as.numeric(d$NSE)
  d$C_Colombia <- as.numeric(d$Pais == 2L) # aquí: Sede_Ficticia_B
  d$C_Peru <- as.numeric(d$Pais == 3L)     # aquí: Sede_Ficticia_C
  d$C_Mexico <- as.numeric(d$Pais == 4L)   # aquí: Sede_Ficticia_D
  # Medias explícitas: no son puntuaciones factoriales. na.rm=FALSE preserva la
  # regla completa aunque esta simulación no contiene faltantes.
  d$Prosocialidad <- rowMeans(d[, s$pp, drop = FALSE], na.rm = FALSE)
  d$ST_comp <- rowMeans(d[, c("Benevolencia", "Universalismo")], na.rm = FALSE)
  d$SE_comp <- rowMeans(d[, c("Poder", "Logro")], na.rm = FALSE)
  attr(d, "tipo_datos") <- "INTEGRAMENTE_ARTIFICIALES_NO_PARTICIPANTES"
  attr(d, "pais_codigos") <- stats::setNames(seq_along(p$sedes), p$sedes)
  attr(d, "semilla") <- p$seed
  d
}

diccionario_artificial <- function(d) {
  # El diccionario documenta cada columna y distingue creación aleatoria de
  # derivación determinista. No incluye valores o frecuencias de participantes.
  s <- esquema_artificial()
  x <- data.frame(variable = names(d), tipo_R = vapply(d, typeof, character(1)),
    origen = "simulado", dominio = "", descripcion_es = "",
    description_en = "", stringsAsFactors = FALSE)
  add <- function(vars, dominio, es, en, origen = "simulado") {
    ix <- match(vars, x$variable)
    stopifnot(!anyNA(ix))
    x$dominio[ix] <<- dominio; x$descripcion_es[ix] <<- es
    x$description_en[ix] <<- en; x$origen[ix] <<- origen
  }
  add("Tipo_dato", "SINTETICO_HIPOTETICO_NO_PARTICIPANTES", "Marca de ejemplo enteramente artificial", "Entirely artificial example marker", "metadato")
  add("ID_sintetico", sprintf("SYN_0001..SYN_%04d", nrow(d)), "Identificador secuencial artificial, sin correspondencia con personas", "Synthetic sequential identifier, unrelated to persons", "metadato")
  add("Pais", "1=A; 2=B; 3=C; 4=D (sedes ficticias)", "Código heredado Pais; no representa países reales", "Legacy Pais code for fictional sites, not real countries")
  add("Sede_ficticia", "Sede_Ficticia_A..Sede_Ficticia_D", "Etiqueta explícita de la sede artificial", "Explicit fictional site label", "derivado de Pais")
  add("Sexo", "0=Hombre; 1=Mujer", "Código binario del esquema original, generado artificialmente", "Artificial binary sex code matching original schema")
  add("Estudios", "0=Pregrado; 1=Postgrado", "Nivel de estudios generado artificialmente", "Artificial undergraduate/postgraduate code")
  add("Edad", "18..60, años enteros", "Edad artificial redondeada y acotada para el ejemplo", "Artificial rounded bounded age in years")
  add("NSE", "1=Bajo; 2=Medio-bajo; 3=Medio; 4=Medio-alto; 5=Alto", "Nivel socioeconómico ordinal artificial", "Artificial ordinal socioeconomic level")
  add("CSR1", "1=Sí; 2=No", "Curso de ética/RSE artificial, codificación original", "Artificial ethics/CSR course, original code", "derivado de CSR1_rec")
  add(c("CSR1_rec", "CSR1_bin"), "0=No; 1=Sí", "Curso de ética/RSE artificial, indicador recodificado", "Artificial ethics/CSR course, recoded indicator", "CSR1_rec simulado; CSR1_bin copia exacta")
  acuerdo <- "1=Totalmente en desacuerdo; 2=En desacuerdo; 3=Ni de acuerdo ni en desacuerdo; 4=De acuerdo; 5=Totalmente de acuerdo"
  semejanza <- "1=No se parece nada a mí; 2=No se parece a mí; 3=Se parece poco a mí; 4=Se parece algo a mí; 5=Se parece a mí; 6=Se parece mucho a mí"
  add(s$sc, acuerdo, "Respuesta artificial de ítem APS del bloque SC", "Artificial APS response in SC block")
  add(s$sbl, acuerdo, "Respuesta artificial de ítem APS del bloque SBL", "Artificial APS response in SBL block")
  add(s$pp, acuerdo, "Respuesta artificial de ítem de conducta prosocial autoinformada", "Artificial self-reported prosocial item response")
  add(s$valores, semejanza, "Respuesta artificial de valor observado, indicador único", "Artificial observed single-item value response")
  add("Sexo_num", "0..1", "Copia numérica exacta de Sexo", "Exact numeric copy of Sexo", "derivado")
  add("Estudios_num", "0..1", "Copia numérica exacta de Estudios", "Exact numeric copy of Estudios", "derivado")
  add(c("Edad_z", "NSE_z"), "real", "Estandarización de toda la base artificial: (x-media)/SD muestral; no usar antes de dividir entrenamiento", "Full synthetic sample z score; not training-only preprocessing", "derivado")
  add("Edad_original", "18..60", "Copia numérica de Edad, no es edad de participantes reales", "Numeric copy of artificial Edad", "derivado")
  add("NSE_original", "1..5", "Copia numérica del NSE artificial", "Numeric copy of artificial NSE", "derivado")
  add("C_Colombia", "0..1", "Nombre heredado: indicador Pais=2, Sede_Ficticia_B, sin inferencia sobre Colombia", "Legacy name: fictional site B indicator, unrelated to Colombia", "derivado")
  add("C_Peru", "0..1", "Nombre heredado: indicador Pais=3, Sede_Ficticia_C, sin inferencia sobre Perú", "Legacy name: fictional site C indicator, unrelated to Peru", "derivado")
  add("C_Mexico", "0..1", "Nombre heredado: indicador Pais=4, Sede_Ficticia_D, sin inferencia sobre México", "Legacy name: fictional site D indicator, unrelated to Mexico", "derivado")
  add("Prosocialidad", "1..5", "Media exacta de PP1, PP2 y PP3; no solo PP1 y PP2", "Exact mean of PP1, PP2 and PP3", "derivado")
  add("ST_comp", "1..6", "Media exacta de Benevolencia y Universalismo", "Exact mean of Benevolencia and Universalismo", "derivado")
  add("SE_comp", "1..6", "Media exacta de Poder y Logro", "Exact mean of Poder and Logro", "derivado")
  stopifnot(all(nzchar(x$descripcion_es)), all(nzchar(x$description_en)))
  x
}

validar_artificial <- function(d, replica, p) {
  # Pruebas del artefacto, sin estimar modelos ni comparar con resultados reales.
  s <- esquema_artificial()
  chequeos <- list()
  check <- function(nombre, condicion) chequeos[[nombre]] <<- isTRUE(condicion)
  check("tamano_y_cupos_segun_parametros", nrow(d) == p$n_por_sede * length(p$sedes) &&
        length(unique(d$Pais)) == length(p$sedes) && all(table(d$Pais) == p$n_por_sede))
  check("marcado_artificial", all(d$Tipo_dato == "SINTETICO_HIPOTETICO_NO_PARTICIPANTES"))
  check("identificadores_artificiales_unicos", !anyDuplicated(d$ID_sintetico) && all(grepl("^SYN_", d$ID_sintetico)))
  check("sedes_ficticias", identical(d$Sede_ficticia, p$sedes[d$Pais]))
  check("sin_NA", !anyNA(d))
  check("reproducible_misma_semilla", identical(d, replica))
  check("ordinales_1_a_5", all(as.matrix(d[, c(s$sc, s$sbl, s$pp)]) %in% 1:5))
  check("valores_1_a_6", all(as.matrix(d[, s$valores]) %in% 1:6))
  check("todas_categorias_en_ejemplo", all(vapply(d[, c(s$sc, s$sbl, s$pp)], function(x) identical(sort(unique(x)), 1:5), logical(1))))
  check("binarios", all(as.matrix(d[, c("Sexo", "Estudios", "CSR1_rec", "CSR1_bin")]) %in% 0:1))
  check("CSR1_sentido_correcto", identical(as.integer(d$CSR1 == 1L), d$CSR1_rec) && identical(d$CSR1_rec, d$CSR1_bin))
  check("edades_y_NSE_validos", all(d$Edad %in% 18:60) && all(d$NSE %in% 1:5))
  check("copias_controles", all(d$Sexo_num == d$Sexo) && all(d$Estudios_num == d$Estudios) && all(d$Edad_original == d$Edad) && all(d$NSE_original == d$NSE))
  check("dummies_sedes_ficticias", all(d$C_Colombia == (d$Pais == 2L)) && all(d$C_Peru == (d$Pais == 3L)) && all(d$C_Mexico == (d$Pais == 4L)))
  check("estandarizacion_exacta", identical(d$Edad_z, as.numeric(scale(d$Edad))) && identical(d$NSE_z, as.numeric(scale(d$NSE))))
  check("PP_tres_items", identical(d$Prosocialidad, rowMeans(d[, s$pp, drop = FALSE], na.rm = FALSE)))
  check("medias_valores", identical(d$ST_comp, rowMeans(d[, c("Benevolencia", "Universalismo")], na.rm = FALSE)) && identical(d$SE_comp, rowMeans(d[, c("Poder", "Logro")], na.rm = FALSE)))
  ans <- data.frame(prueba = names(chequeos), cumple = unlist(chequeos), row.names = NULL)
  if (!all(ans$cumple)) stop("Falló una comprobación de datos artificiales: ", paste(ans$prueba[!ans$cumple], collapse = ", "))
  ans
}

exportar_artificial <- function(salida) {
  dir.create(salida, recursive = TRUE, showWarnings = FALSE)
  p <- parametros_artificiales()
  d <- simular_base_artificial(p)
  replica <- simular_base_artificial(p)
  tests <- validar_artificial(d, replica, p)
  dic <- diccionario_artificial(d)
  csv <- function(x, nombre) utils::write.csv(x, file.path(salida, nombre), row.names = FALSE, fileEncoding = "UTF-8", na = "")
  csv(d, "Datos_SINTETICOS_V3.csv")
  saveRDS(d, file.path(salida, "Datos_SINTETICOS_V3.rds"), version = 3)
  csv(dic, "diccionario_SINTETICO_V3.csv")
  csv(tests, "comprobaciones.csv")
  # RDS conserva clases y los atributos de artificialidad. CSV es interoperable.
  # Las etiquetas SPSS se exportan solo si haven ya está instalado; no se instala
  # ningún paquete ni se utiliza una conexión a Internet.
  sav_ok <- FALSE
  if (requireNamespace("haven", quietly = TRUE)) {
    d_sav <- d
    labs <- list(
      Pais = stats::setNames(1:4, p$sedes),
      Sexo = c(Hombre = 0, Mujer = 1),
      Estudios = c(Pregrado = 0, Postgrado = 1),
      CSR1 = c(Si = 1, No = 2), CSR1_rec = c(No = 0, Si = 1),
      CSR1_bin = c(No = 0, Si = 1),
      NSE = c(Bajo = 1, Medio_bajo = 2, Medio = 3, Medio_alto = 4, Alto = 5))
    s <- esquema_artificial()
    likert <- c(Totalmente_en_desacuerdo = 1, En_desacuerdo = 2,
                Ni_de_acuerdo_ni_en_desacuerdo = 3, De_acuerdo = 4,
                Totalmente_de_acuerdo = 5)
    valor <- c(No_se_parece_nada_a_mi = 1, No_se_parece_a_mi = 2,
               Se_parece_poco_a_mi = 3, Se_parece_algo_a_mi = 4,
               Se_parece_a_mi = 5, Se_parece_mucho_a_mi = 6)
    for (v in c(s$sc, s$sbl, s$pp)) labs[[v]] <- likert
    for (v in s$valores) labs[[v]] <- valor
    for (v in names(d_sav)) attr(d_sav[[v]], "label") <- dic$descripcion_es[match(v, dic$variable)]
    for (v in names(labs)) d_sav[[v]] <- haven::labelled(d_sav[[v]], labels = labs[[v]], label = dic$descripcion_es[match(v, dic$variable)])
    haven::write_sav(d_sav, file.path(salida, "Datos_SINTETICOS_V3.sav"))
    sav_ok <- TRUE
  }
  # Tabla de TODOS los parámetros hipotéticos, legible sin paquetes adicionales.
  tabla_p <- data.frame(parametro = names(p),
    valor_hipotetico = vapply(p, function(x) paste(x, collapse = " | "), character(1)),
    procedencia = "Elección didáctica explícita; no estimada con datos reales",
    stringsAsFactors = FALSE)
  csv(tabla_p, "parametros_hipoteticos.csv")
  saveRDS(p, file.path(salida, "parametros_hipoteticos.rds"), version = 3)
  csv(data.frame(Pais = 1:4, Sede_ficticia = p$sedes, N_artificial = p$n_por_sede), "sedes_ficticias.csv")
  # La única frecuencia exportada pertenece al ejemplo artificial recién creado.
  s <- esquema_artificial()
  freq <- do.call(rbind, lapply(c(s$sc, s$sbl, s$pp, s$valores), function(v) {
    lev <- if (v %in% s$valores) 1:6 else 1:5
    data.frame(variable = v, categoria = lev,
               n_artificial = as.integer(table(factor(d[[v]], levels = lev))))
  }))
  csv(freq, "frecuencias_SOLO_SINTETICAS.csv")
  # Comprobar la ida y vuelta del RDS sin leer ningún archivo de participantes.
  stopifnot(identical(d, readRDS(file.path(salida, "Datos_SINTETICOS_V3.rds"))))
  vuelta <- utils::read.csv(file.path(salida, "Datos_SINTETICOS_V3.csv"), check.names = FALSE)
  stopifnot(identical(names(d), names(vuelta)), nrow(vuelta) == nrow(d))
  num <- vapply(d, is.numeric, logical(1))
  stopifnot(max(abs(as.matrix(d[, num]) - as.matrix(vuelta[, num]))) < 1e-12)
  writeLines(c("DATOS INTEGRAMENTE ARTIFICIALES: no contienen participantes.",
    paste("Semilla:", p$seed), paste("RNG:", paste(RNGkind(), collapse = "; ")),
    paste("Filas:", nrow(d)), paste("Columnas:", ncol(d)),
    paste("Pruebas internas correctas:", nrow(tests)),
    "Ida/vuelta CSV: valores numéricos coinciden con tolerancia 1e-12.",
    "Ida/vuelta RDS: objeto idéntico, incluidos atributos.",
    paste("Exportación SAV opcional:", sav_ok),
    "Se consultaron cero archivos de participantes y cero recursos de red.",
    "No se estimaron modelos ni se calcularon tablas empíricas del estudio.",
    "", capture.output(sessionInfo())), file.path(salida, "registro_generacion_y_sessionInfo.txt"), useBytes = TRUE)
  invisible(list(data = d, tests = tests, output = normalizePath(salida), sav = sav_ok))
}

if (sys.nframe() == 0L) {
  # Al abrir/source() en RStudio solo se definen funciones. Para generar desde
  # RStudio: exportar_artificial(file.path(getwd(), 'datos_sinteticos')).
  # Rscript resuelve por defecto una carpeta contigua a este mismo archivo.
  args <- commandArgs(trailingOnly = TRUE)
  archivo <- grep("^--file=", commandArgs(), value = TRUE)
  # Rscript en algunas versiones de macOS codifica espacios como ~+~ dentro
  # de --file. Se decodifica solo ese argumento de ruta antes de normalizarlo.
  ruta_script <- if (length(archivo)) gsub("~+~", " ", sub("^--file=", "", archivo[1]), fixed = TRUE) else NULL
  base <- if (length(ruta_script)) dirname(normalizePath(ruta_script)) else getwd()
  salida <- if (length(args)) args[1] else file.path(base, "datos_sinteticos")
  r <- exportar_artificial(salida)
  cat("Generados", nrow(r$data), "casos completamente artificiales en:\n", r$output, "\n")
}
