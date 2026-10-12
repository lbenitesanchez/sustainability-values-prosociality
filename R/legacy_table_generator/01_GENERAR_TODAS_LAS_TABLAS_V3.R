#!/usr/bin/env Rscript
# GENERACIÓN DOCUMENTADA DE LAS24 TABLAS DE LOS DOS WORDV3.
# Manuscript_V3_revised_contribution_Priale2024_verified_sources.docx
# Supplementary_V3_revised_contribution_Priale2024_verified_sources.docx
#
# ALCANCE: extrae estadísticas de modelos guardados y formatea exportaciones
# ya reproducidas desde la base en la auditoría. NO vuelve a estimar modelos.
# NO usa datos artificiales ni asegura que estos reproduzcan las cifras reales.
# Los RDS lavaan contienen datos privados y no forman parte del ZIP público.
#
# USO EN RSTUDIO (Source sólo define; no ejecuta sin autorización explícita):
# source("01_GENERAR_TODAS_LAS_TABLAS_V3.R", encoding="UTF-8")
# generar_todas_las_tablas_V3()
# Para rutas distintas, pase proyecto_fuente=, raiz_auditoria=,
# proyecto_referencia=, directorio_word=, directorio_modulos= y salida=.
#
# USO EN TERMINAL:
# Rscript --vanilla 01_GENERAR_TODAS_LAS_TABLAS_V3.R
# Rscript --vanilla 01_GENERAR_TODAS_LAS_TABLAS_V3.R --project=/ruta/proyecto --output=/ruta/nueva
# Opciones: --project --output --audit --reference --word --modules.
# Alias en castellano: --proyecto --salida --raiz-auditoria --referencia --word --modulos.
# Use comillas alrededor de cada argumento completo si su ruta tiene espacios.
#
# MAPA RÁPIDO (detalles y líneas exactas en mapa_tablas.csv):
# T1, T2,S1–S4→R/01; S5–S8→R/02; T3, S9→R/06;
# S10–S16, S21→R/07; S17–S20→R/08.

# Localiza este archivo tanto con Rscript como con source(), sin depender
# del directorio de trabajo. source() aporta el campo ofile de su marco.
.V3_SCRIPT_FILE <- local({
  ofiles <- lapply(sys.frames(), function(fr) fr$ofile)
  ofiles <- Filter(function(z) is.character(z) && length(z)==1L && nzchar(z), ofiles)
  if(length(ofiles)) normalizePath(tail(ofiles,1L)[[1L]],mustWork=TRUE) else {
    arg <- grep("^--file=",commandArgs(trailingOnly=FALSE),value=TRUE)
    if(length(arg)) normalizePath(gsub("~+~"," ",sub("^--file=","",arg[1L]),fixed=TRUE),mustWork=TRUE) else
      stop("No se puede localizar el script. Use source() o Rscript con su ruta.")
  }
})
.V3_PACKAGE_DIR <- dirname(.V3_SCRIPT_FILE)

generar_todas_las_tablas_V3 <- function(
    proyecto_fuente=file.path(raiz_auditoria,"proyecto"),
    salida=file.path(.V3_PACKAGE_DIR,paste0("tablas_generadas_",format(Sys.time(),"%Y%m%d_%H%M%S"))),
    raiz_auditoria=Sys.getenv("MA_AUDIT_ROOT", file.path(dirname(.V3_PACKAGE_DIR),"reproduccion-integral-V3-2026-10-10")),
    proyecto_referencia=file.path(raiz_auditoria,"referencia"),
    directorio_word=file.path(raiz_auditoria,"Word_referencia"),
    proyecto_codigo=file.path(raiz_auditoria,"proyecto"),
    directorio_modulos=file.path(.V3_PACKAGE_DIR,"R"),
    fallar_si_diferencias=TRUE) {
# [REQUISITOS] Esta operación exige el proyecto analítico y sus RDS/CSV
# privados ya reproducidos, la referencia congelada y los dos Word V3.
# El ZIP público de datos artificiales no contiene estas fuentes: es normal
# que no pueda reconstruir los resultados empíricos publicados sin acceso local.
  stopifnot(dir.exists(proyecto_fuente),dir.exists(proyecto_referencia),dir.exists(directorio_word))
  proyecto_fuente <- normalizePath(proyecto_fuente)
  proyecto_referencia <- normalizePath(proyecto_referencia)
  proyecto_codigo <- normalizePath(proyecto_codigo)
  salida <- normalizePath(salida,mustWork=FALSE)
  stopifnot(!dir.exists(salida),
    !startsWith(salida,paste0(proyecto_fuente,"/")),
    !startsWith(salida,paste0(proyecto_referencia,"/")),
    !startsWith(salida,paste0(proyecto_codigo,"/")))
# [PROTECCIÓN] La salida debe ser NUEVA y estar fuera de las fuentes.
# No se sobrescriben modelos, CSV ni documentos; tampoco se edita ningún Word.
  dir.create(salida,recursive=TRUE,showWarnings=FALSE)
  for(n in c("registro","resultados_base","tablas_finales"))dir.create(file.path(salida,n))
  original <- proyecto_fuente
  .libPaths(c(file.path(original,"revision-observaciones-2026-10-09/diagnosticos/R-library"),
    file.path(original,"revision-revista-2026-10-04/auditoria-codigos/R-library"),.libPaths()))
# [ENTORNO] No se instalan paquetes. Se exige lavaan/semTools para extraer
# parámetros y fiabilidad, jsonlite para registros y digest para SHA256.
  required <- c("lavaan","semTools","jsonlite","digest")
  stopifnot(all(vapply(required,requireNamespace,logical(1),quietly=TRUE)))
  inicio <- Sys.time()
# [AISLAMIENTO] Configura sólo rutas temporales para los módulos y restaura
# las variables de entorno al salir. Los módulos se cargan en entorno propio.
  ambiente_previo <- Sys.getenv(c("MA_PROJECT_DIR","MA_REVISION_DIR"),unset=NA)
  on.exit({for(n in names(ambiente_previo)) if(is.na(ambiente_previo[n]))Sys.unsetenv(n)else do.call(Sys.setenv,as.list(ambiente_previo[n]))},add=TRUE)
  Sys.setenv(MA_PROJECT_DIR=proyecto_fuente,MA_REVISION_DIR=salida)
  revision_fuente <- file.path(proyecto_fuente,"revision-observaciones-2026-10-09")
  revision_referencia <- file.path(proyecto_referencia,"revision-observaciones-2026-10-09")
  codigo <- normalizePath(directorio_modulos,mustWork=TRUE)
  codigos <- file.path(codigo,c("01_tablas_modelos.R","02_tablas_items_ordinal_prediccion.R",
    "06_TABLAS_PREDICCION_AGRUPADA.R","07_TABLAS_DIAGNOSTICOS_Y_SENSIBILIDAD.R",
    "08_TABLAS_REGULARIZACION_E_INCERTIDUMBRE.R"))
  stopifnot(all(file.exists(codigos)))
  sha <- function(p)digest::digest(p,file=TRUE,algo="sha256")
  escribir <- function(x,p)utils::write.csv(x,p,row.names=FALSE,na="",fileEncoding="UTF-8")
  exportar_matriz <- function(x,p)utils::write.table(x,p,sep=",",quote=TRUE,qmethod="double",row.names=FALSE,col.names=FALSE,fileEncoding="UTF-8",na="NA")
  manifest <- list();record <- function(path,role){
    manifest[[length(manifest)+1L]] <<- data.frame(path=path,role=role,sha256_before=sha(path),bytes=file.info(path)$size)
  }
  for(p in codigos)record(p,"formatting_code_read_only")
# [CONTRATO HISTÓRICO] Se copia sólo el contrato de cotejo en salida.
# No alimenta las estimaciones: cada módulo calcula sus tablas antes de leerlo.
  contrato <- file.path(revision_fuente,"registro/contrato_original.json")
  stopifnot(file.exists(contrato));record(contrato,"historical_contract_module02")
  file.copy(contrato,file.path(salida,"registro/tablas_word_esperadas.json"))
# [CONFIGURACIÓN CLAVE / T1, T2,S1–S8] prueba_sin_cv=TRUE desactiva
# definitivamente en esta ejecución la antigua CV de 175 ajustes. Se conservan
# los estadísticos de modelos existentes y se excluyen sus columnas históricas.
  cfg <- list(proyecto=proyecto_fuente,salida=file.path(salida,"resultados_base"),
    registro=file.path(salida,"registro"),prueba_sin_cv=TRUE)
  # Los módulos reciben el proyecto-fuente explícito. Esta etapa ejecuta
  # extractores de tablas; no vuelve a ejecutar el productor de modelos SEM.
  env <- new.env(parent=globalenv())
  sys.source(codigos[1],envir=env)
  sys.source(codigos[2],envir=env)
# [T1, T2,S1, S2, S3, S4] Ver R/01_tablas_modelos.R por fórmula y fuente.
# [S5, S6, S7, S8] Ver R/02_tablas_items_ordinal_prediccion.R.
# El módulo02 revisa además T3/S9 HISTÓRICAS; aquí se descartan ambas.
  a <- env$generar_tablas_modelos(cfg)
  b <- env$generar_tablas_items(cfg)
# [T1 ACTUAL] Columnas 1–4 y 7 = modelo, CFI, RMSEA, SRMR y npar.
# T2 y S1–S8 conservan el resultado completo de sus respectivos extractores.
  tablas_base <- c(list(T1=a$T1[,c(1:4,7),drop=FALSE],T2=a$T2),a[paste0("S",1:4)],b[paste0("S",5:8)])
  for(id in names(tablas_base))exportar_matriz(tablas_base[[id]],file.path(salida,"resultados_base",paste0(id,".csv")))
  # T3/S9 del módulo02 solo satisfacen su auditoría histórica. Se descartan
  # antes de formar las tablas finales, sustituidas por la predicción agrupada.
  # Copiamos agregados, no los RDS privados ni carpetas con ejecución activa.
  for(rel in c("evidencia/observaciones_2_3_6_7","diagnosticos/resultados","prediccion")) {
    origen <- file.path(revision_fuente,rel)
    archivos <- list.files(origen,pattern="\\.csv$",full.names=TRUE,recursive=TRUE)
    archivos <- archivos[!grepl("/privados/|/R-library/",archivos)]
    for(p in archivos) {
      sub <- substring(p,nchar(origen)+2L)
      target <- file.path(salida,rel,sub);dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
      record(p,"saved_aggregate_read_only")
      stopifnot(file.copy(p,target),identical(sha(p),sha(target)))
    }
  }
# [T3, S9] Módulo06: resultados predictivos AGRUPADOS definitivos.
# [S10, S11, S12, S13, S14, S15, S16, S21] Módulo07: diagnóstico/sensibilidad.
# [S17, S18, S19, S20] Módulo08: ridge, OOB, calibración y LISO.
# Todos formatean agregados ya reproducidos: no se ajustan modelos nuevos.
  for(i in 3:5)sys.source(codigos[i],envir=env)
  env$generar_prediccion_agrupada()
  env$formatear_evidencia();env$formatear_identificacion();env$formatear_probit()
  env$formatear_extension_predictiva()
# [MAPA DE COTEJO] Se mantiene el mapa congelado de la auditoría como
# contrato de rutas de salida. El mapa_tablas.csv de este paquete añade
# títulos literales, fuentes y líneas de los módulos documentados.
  mapa <- utils::read.csv(file.path(revision_referencia,"registro/MAPA_TABLAS_CODIGO.csv"),stringsAsFactors=FALSE)
  stopifnot(nrow(mapa)==24L,setequal(mapa$tabla,c(paste0("T",1:3),paste0("S",1:21))))
  read_matrix <- function(p)unname(as.matrix(utils::read.csv(p,header=FALSE,check.names=FALSE,
    colClasses="character",na.strings=NULL,stringsAsFactors=FALSE,fileEncoding="UTF-8")))
# [COTEJO] Normaliza sólo saltos de línea y UTF-8. No redondea de nuevo.
# Se admite únicamente la sustitución editorial documentada M0 Original/Baseline.
# Cualquier otra celda distinta hace fallar la comprobación.
  normal <- function(s)gsub("\r\n?","\n",enc2utf8(s),perl=TRUE)
  label_change <- function(a,b)(a=="M0 Original" & b=="M0 Baseline") | (a=="M0 Baseline" & b=="M0 Original")
  filas <- list();celdas <- list();nf<-0L;nc<-0L
  compare <- function(id,actual,esperada,tipo) {
    dims <- identical(dim(actual),dim(esperada));same<-labels<-mism<-0L
    if(dims)for(r in seq_len(nrow(actual)))for(c in seq_len(ncol(actual))) {
      x<-normal(actual[r,c]);y<-normal(esperada[r,c]);nc<<-nc+1L
      estado <- if(identical(x,y))"identical"else if(label_change(x,y))"intentional_M0_label"else"different"
      same<-same+(estado=="identical");labels<-labels+(estado=="intentional_M0_label");mism<-mism+(estado=="different")
      celdas[[nc]]<<-data.frame(table=id,comparison=tipo,row=r,column=c,computed=x,reference=y,status=estado)
    }
    nf<<-nf+1L;filas[[nf]]<<-data.frame(table=id,comparison=tipo,rows=nrow(actual),columns=ncol(actual),dimensions_match=dims,
      identical_cells=same,intentional_label_cells=labels,different_cells=if(dims)mism else NA_integer_,pass=dims && mism==0L)
  }
  # Word: solo tablas OOXML planas; no se requiere instalar paquetes.
  # La lectura conserva entidades XML y saltos. Las dimensiones y las celdas
  # se comparan además con los CSV, por lo que un fallo de lectura no se acepta.
  word_tables <- function(path) {
    record(path,"latest_Word_reference")
    tmp<-tempfile("word_xml_");dir.create(tmp);on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
    utils::unzip(path,files="word/document.xml",exdir=tmp)
    xml<-paste(readLines(file.path(tmp,"word/document.xml"),warn=FALSE,encoding="UTF-8"),collapse="")
    tokens<-function(pattern,x) regmatches(x,gregexpr(pattern,x,perl=TRUE))[[1L]]
    decode<-function(x) {
      for(pattern in c("&#x[0-9a-fA-F]+;","&#[0-9]+;")) {
        hits<-unique(tokens(pattern,x))
        for(h in hits) {
          hex<-startsWith(h,"&#x")
          n<-if(hex)strtoi(substr(h,4L,nchar(h)-1L),base=16L)else as.integer(substr(h,3L,nchar(h)-1L))
          x<-gsub(h,intToUtf8(n),x,fixed=TRUE)
        }
      }
      # Ampersand last prevents double decoding literal entity-looking text.
      for(h in c("&lt;","&gt;","&quot;","&apos;","&amp;"))
        x<-gsub(h,c("<",">","\"","'","&")[match(h,c("&lt;","&gt;","&quot;","&apos;","&amp;"))],x,fixed=TRUE)
      x
    }
    tabs<-tokens("(?s)<w:tbl(?:\\s[^>]*)?>.*?</w:tbl>",xml)
    lapply(tabs,function(t) {
      stopifnot(length(tokens("<w:tbl(?:\\s[^>]*)?>",t))==1L)
      rows<-tokens("(?s)<w:tr(?:\\s[^>]*)?>.*?</w:tr>",t)
      cells<-lapply(rows,function(row) {
        cc<-tokens("(?s)<w:tc(?:\\s[^>]*)?>.*?</w:tc>",row)
        vapply(cc,function(cell) {
          paragraphs<-tokens("(?s)<w:p(?:\\s[^>]*)?>.*?</w:p>",cell)
          paste(vapply(paragraphs,function(paragraph) {
            runs<-tokens("(?s)<w:t(?:\\s[^>]*)?>.*?</w:t>|<w:(?:br|cr|tab)(?:\\s[^>]*)?/>",paragraph)
            paste(vapply(runs,function(z) {
              if(grepl("^<w:(br|cr)",z))"\n"else if(startsWith(z,"<w:tab"))"\t"else
                decode(sub("</w:t>$","",sub("^<w:t(?:\\s[^>]*)?>","",z,perl=TRUE)))
            },character(1)),collapse="")
          },character(1)),collapse="\n")
        },character(1))
      })
      stopifnot(length(unique(lengths(cells)))==1L)
      unname(do.call(rbind,cells))
    })
  }
  wf<-list.files(directorio_word,pattern="\\.docx$",full.names=TRUE)
  mf<-wf[grepl("/Manuscript",wf)];sf<-wf[grepl("/Supplementary",wf)]
  stopifnot(length(mf)==1L,length(sf)==1L)
  wm<-word_tables(mf);ws<-word_tables(sf);stopifnot(length(wm)==3L,length(ws)==21L)
# [24 TABLAS / 48COTEJOS] Cada CSV generado se compara celda por celda
# con el CSV congelado y con su tabla enWord:3 del manuscrito + 21 del suplemento.
# La salida verificada conserva encabezados como primera fila del CSV.
  for(i in seq_len(nrow(mapa))) {
    id<-mapa$tabla[i];generated<-file.path(salida,mapa$csv[i]);reference<-file.path(revision_referencia,mapa$csv[i])
    stopifnot(file.exists(generated),file.exists(reference));record(reference,"reference_table_CSV")
    actual<-read_matrix(generated)
    compare(id,actual,read_matrix(reference),"reference_CSV")
    wt<-if(startsWith(id,"T"))wm[[as.integer(sub("T","",id))]]else ws[[as.integer(sub("S","",id))]]
    compare(id,actual,wt,"latest_Word")
    exportar_matriz(actual,file.path(salida,"tablas_finales",paste0(id,"_verificada.csv")))
  }
# [EVIDENCIA] Exporta resumen por tabla, detalle por celda, diferencias
# y hashes antes/después. sessionInfo registra R, plataforma y paquetes usados.
# Un cotejo completo acredita tablas, no sustituye la auditoría científica.
  summary<-do.call(rbind,filas);cells<-do.call(rbind,celdas)
  escribir(summary,file.path(salida,"registro/cotejo_24_tablas.csv"))
  escribir(cells,file.path(salida,"registro/cotejo_todas_celdas.csv"))
  escribir(cells[cells$status!="identical",,drop=FALSE],file.path(salida,"registro/diferencias_y_etiquetas.csv"))
  manifest<-unique(do.call(rbind,manifest));manifest$sha256_after<-vapply(manifest$path,sha,character(1));manifest$unchanged<-manifest$sha256_before==manifest$sha256_after
  escribir(manifest,file.path(salida,"registro/fuentes_leidas.csv"))
  # Los módulos guardan además hashes antes/después de los RDS que leen.
  comprobacion_modulo02 <- jsonlite::read_json(file.path(salida,"registro/02_comprobacion_reproduccion.json"),simplifyVector=TRUE)
  report<-list(status=if(all(summary$pass)&&all(manifest$unchanged))"EXTRACTION_AND_24_TABLE_COMPARISONS_PASSED"else"DIFFERENCES_REQUIRE_REVIEW",
    operation="Extraction from existing fitted RDS plus formatting saved aggregates; no model refitting in this step",
    source_project=proyecto_fuente,reference_project=proyecto_referencia,code_project=proyecto_codigo,documented_modules=codigo,output=salida,
    new_model_fits=0L,number_tables=24L,number_comparisons=nrow(summary),cells_compared=nrow(cells),
    identical_cells=sum(cells$status=="identical"),intentional_M0_label_cells=sum(cells$status=="intentional_M0_label"),
    differing_cells=sum(cells$status=="different"),all_compared_tables_pass=all(summary$pass),
    recorded_sources_unchanged=all(manifest$unchanged),module02_RDS_sources_unchanged=comprobacion_modulo02$original_sources_unchanged,
    old_T3_S9="Recalculated only to satisfy historical module audit; discarded, not used in final current24tables",
    limitations=c("No claim of refitting from raw SPSS by this step","Aggregates for T3/S9/S10-S21 are formatted, not reestimated here","Numerical narrative outside tables requires separate checking"),
    started=as.character(inicio),finished=as.character(Sys.time()),elapsed_seconds=as.numeric(difftime(Sys.time(),inicio,units="secs")),
    versions=setNames(lapply(required,function(p)as.character(utils::packageVersion(p))),required))
  jsonlite::write_json(report,file.path(salida,"registro/resultado_extraccion_y_cotejo.json"),pretty=TRUE,auto_unbox=TRUE)
  writeLines(capture.output(sessionInfo()),file.path(salida,"registro/sessionInfo.txt"))
  print(report[c("status","new_model_fits","number_tables","number_comparisons","cells_compared","intentional_M0_label_cells","differing_cells")])
  if(fallar_si_diferencias && (!all(summary$pass)||!all(manifest$unchanged)))stop("Hay diferencias: consultar registro/diferencias_y_etiquetas.csv")
  invisible(report)
}

# Rscript ejecuta; source() sólo deja disponible la función anterior.
if(sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  allowed <- c(project="proyecto_fuente",output="salida",audit="raiz_auditoria",
               reference="proyecto_referencia",word="directorio_word",modules="directorio_modulos",
               proyecto="proyecto_fuente",salida="salida",`raiz-auditoria`="raiz_auditoria",
               referencia="proyecto_referencia",modulos="directorio_modulos")
  options <- list()
  for(arg in args) {
    if(!grepl("^--[^=]+=.+$",arg)) stop("Formato esperado: --opcion=/ruta; recibido: ",arg)
    name <- sub("=.*$","",sub("^--","",arg))
    if(!name %in% names(allowed)) stop("Opción no reconocida: ",name)
    key <- unname(allowed[name])
    if(key %in% names(options)) stop("Opción duplicada: ",name)
    options[[key]] <- sub("^--[^=]+=","",arg)
  }
  do.call(generar_todas_las_tablas_V3,options)
}
