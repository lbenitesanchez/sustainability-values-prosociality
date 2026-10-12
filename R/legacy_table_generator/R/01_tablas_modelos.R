# MÓDULO 01 — TABLAS ACTUALES T1, T2 Y S1–S4.
# El principal llama generar_tablas_modelos(config) con prueba_sin_cv=TRUE.
# Así se extraen estadísticas de 28 ajustes guardados: NO se vuelve a ajustar
# ningún modelo. T1 actual conserva CFI, RMSEA, SRMR y npar; excluye las dos
# columnas de predicción históricas. T3 y S9 actuales se generan en módulo06.
# Las funciones de la antigua CV de 175 ajustes se conservan por trazabilidad,
# pero NO son parte de la ejecución V3. No llamarlas para generar estas24 tablas.
# Búsqueda en RStudio: [T1], [T2], [S1], [S2], [S3], [S4].
# Las cifras del Word sólo se leen después de calcular las tablas, para cotejo.

#### Preparación y funciones de apoyo ####

# [AUXILIAR / T1, T2,S1–S4] Huella SHA256: comprueba que leer los modelos
# no altera los archivos. La huella no se usa para calcular resultados.
.tab_modelos_sha <- function(path) digest::digest(path, algo="sha256", file=TRUE)
# [FORMATO / T1, T2,S1–S4] Redondea sólo la representación final (3 decimales
# por defecto). Los valores no finitos permanecen NA, nunca se imputan como cero.
.tab_modelos_num <- function(x, digits=3L) {
  # NA se muestra como NA: no se inventa un índice cuando el test no existe.
  ifelse(is.finite(x), sprintf(paste0("%.", digits, "f"), x), "NA")
}
# [FORMATO / T1, T2,S1–S4] Une encabezado y filas de texto, sin nombres de fila.
# Los CSV numéricos se guardan aparte antes de esta conversión a texto.
.tab_modelos_matriz <- function(header, rows) {
  ans <- rbind(as.character(header), do.call(rbind, rows))
  storage.mode(ans) <- "character"
  dimnames(ans) <- NULL
  ans
}
# [SALIDA / T1, T2,S1–S4] CSV numéricos con sus columnas, sin índice de filas.
.tab_modelos_exportar <- function(x, nombre, salida) {
  utils::write.csv(x, file.path(salida, nombre), row.names=FALSE, na="", fileEncoding="UTF-8")
}

# [LEGADO INACTIVO EN V3] Recupera definiciones para la CV histórica.
# El principal actual nunca entra en esta función ni en .tab_modelos_cv().
.tab_modelos_extraer_funciones_originales <- function(proyecto) {
  # Leer un script NO implica ejecutarlo entero. El script original contiene
  # ajustes y escrituras: se evalúan exclusivamente asignaciones autorizadas
  # que construyen los modelos y definen sus funciones. Se excluyen los bucles
  # originales de ajuste, source(), saveRDS(), writeLines() y demás efectos.
  script <- file.path(proyecto,"sem-alternativas-2026-10-08/codigo/comparar_modelos_SEM.R")
  environment <- new.env(parent=environment(.tab_modelos_extraer_funciones_originales))
  wanted <- c("sc","sbl","aps","pp","values","controls","join","measurement",
              "specs","structural","safe","fit_fun")
  expressions <- as.list(parse(script))
  for (name in wanted) {
    hit <- which(vapply(expressions, function(expr) {
      is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        is.symbol(expr[[2]]) && identical(as.character(expr[[2]]), name)
    }, logical(1)))
    stopifnot(length(hit)==1L)
    eval(expressions[[hit]], environment)
  }
  # La matriz target conserva el mismo ORDEN de 19 indicadores; sus ceros
  # orientan la rotación, mientras los unos dejan libres las cargas primarias.
  target <- matrix(0, nrow=length(environment$aps), ncol=2L,
                   dimnames=list(environment$aps,c("SC","SBL")))
  target[environment$sc,"SC"] <- 1
  target[environment$sbl,"SBL"] <- 1
  environment$target <- target
  for (id in names(environment$specs)) {
    spec <- environment$specs[[id]]
    spec$id <- id
    measurement <- spec$measurement
    if (grepl("_S1_",id))
      measurement <- sub(paste0("\n",spec$covariance),"",measurement,fixed=TRUE)
    spec$syntax <- paste(measurement,
                         environment$structural(spec$factors,spec$covariance),sep="\n")
    archived <- paste(readLines(file.path(proyecto,"sem-alternativas-2026-10-08/resultados",paste0(id,"_syntax.txt")),warn=FALSE),collapse="\n")
    stopifnot(identical(spec$syntax,archived))
    environment$specs[[id]] <- spec
  }
  # Este segundo script sólo define ayudantes y run_sem_cv. No se llama su
  # función de ejecución original, que escribiría en el directorio anterior.
  helper <- file.path(proyecto,"sem-alternativas-2026-10-08/codigo/cross_validate_models.R")
  sys.source(helper,envir=environment)
  environment
}

#### [T1-CV] Reestimación documentada de los 175 ajustes de entrenamiento ####

.tab_modelos_cv <- function(data, definitions, config, fingerprints) {
  # Se vuelven a estimar 7 modelos x 5 folds x 5 repeticiones. Los estimadores
  # de entrenamiento son ML sin SE/test, idénticos al protocolo original.
  # Esto REPRODUCE la evaluación previa; no es una nueva búsqueda de modelos.
  # Cada archivo privado guarda ajuste, momentos, filas y estados de RNG para
  # que una interrupción pueda reanudarse sin alterar las rotaciones ESEM.
  cache_dir <- file.path(config$privado,"T1_CV_175_ajustes")
  dir.create(cache_dir,recursive=TRUE,showWarnings=FALSE)
  protocol_file <- file.path(cache_dir,"protocolo_cache.rds")
  protocol <- list(version="T1-CV-reproduction-v1",N=653L,folds=5L,repeats=5L,
    seed=20261008L,models=names(definitions$specs),fingerprints=fingerprints,
    syntax=lapply(definitions$specs,`[[`,"syntax"),
    training="ML; se=none; test=none; fixed.x=TRUE; same target rotation rstarts20",
    target="Observed mean of PP1 PP2 PP3; remaining observed variables as predictors",
    R_version=R.version.string,lavaan_version=as.character(utils::packageVersion("lavaan")),
    RNG_kind=c("Mersenne-Twister","Inversion","Rejection"))
  if (file.exists(protocol_file)) {
    if (!identical(readRDS(protocol_file),protocol))
      stop("El cache privado T1 pertenece a otra configuración. No se reutiliza ni sobrescribe.")
  } else saveRDS(protocol,protocol_file)
  old_seed_exists <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  old_rng_kind <- RNGkind()
  if (old_seed_exists) old_seed <- get(".Random.seed",envir=.GlobalEnv)
  on.exit({
    do.call(RNGkind,as.list(old_rng_kind))
    if (old_seed_exists) assign(".Random.seed",old_seed,envir=.GlobalEnv)
    else if (exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv)
  },add=TRUE)
  # El Rscript original usó el generador predeterminado. Se fija aquí para
  # que un ajuste previo de RNGkind en RStudio no cambie las particiones.
  RNGkind("Mersenne-Twister","Inversion","Rejection")
  strata <- definitions$.sem_cv_strata(data)
  assignment <- vector("list",5L)
  # Mismo orden de países que los datos y mismo algoritmo exacto del original.
  for (r in 1:5) {
    set.seed(20261008L+r)
    a <- integer(nrow(data))
    for (s in unique(strata)) {
      rows <- which(strata==s)
      a[rows] <- sample(rep(1:5,length.out=length(rows)))
    }
    assignment[[r]] <- a
  }
  assignment_file <- file.path(cache_dir,"asignaciones_privadas.rds")
  if (file.exists(assignment_file)) stopifnot(identical(readRDS(assignment_file),assignment))
  else saveRDS(assignment,assignment_file)
  previous_assignments <- file.path(config$proyecto,"sem-valores-relativos-2026-10-08/resultados/privados/cv_assignments_private.rds")
  if (file.exists(previous_assignments))
    stopifnot(identical(assignment,readRDS(previous_assignments)))
  pp <- c("PP1","PP2","PP3")
  records <- list(); index <- 0L; fitted <- reused <- 0L
  shared_ov <- NULL
  for (r in 1:5) for (id in names(definitions$specs)) {
    for (k in 1:5) {
      train <- which(assignment[[r]]!=k); test <- which(assignment[[r]]==k)
      filename <- file.path(cache_dir,sprintf("%s_r%02d_f%02d.rds",id,r,k))
      before_rng <- get(".Random.seed",envir=.GlobalEnv)
      if (file.exists(filename)) {
        cached <- readRDS(filename)
        stopifnot(identical(cached$protocol_digest,digest::digest(protocol,algo="sha256")),
                  identical(cached$model,id),identical(cached$training_rows,train),
                  identical(cached$testing_rows,test),identical(cached$random_before,before_rng))
        fit <- cached$fit
        checked <- definitions$.sem_cv_fit_check(fit)
        stopifnot(checked$admissible,identical(cached$mean,checked$mean),
                  identical(cached$covariance,checked$covariance))
        assign(".Random.seed",cached$random_after,envir=.GlobalEnv)
        reused <- reused+1L
      } else {
        fit <- definitions$fit_fun(definitions$specs[[id]],data[train,,drop=FALSE],
                                   estimator="ML",fast=TRUE)
        checked <- definitions$.sem_cv_fit_check(fit)
        if (!checked$admissible)
          stop(sprintf("T1 CV no admisible: %s r%d f%d: %s",id,r,k,checked$reason))
        cached <- list(protocol_digest=digest::digest(protocol,algo="sha256"),model=id,
          training_rows=train,testing_rows=test,fit=fit,mean=checked$mean,
          covariance=checked$covariance,diagnostics=checked,
          random_before=before_rng,random_after=get(".Random.seed",envir=.GlobalEnv))
        saveRDS(cached,filename)
        fitted <- fitted+1L
      }
      if (is.null(shared_ov)) shared_ov <- checked$ov
      stopifnot(setequal(shared_ov,checked$ov))
      mu <- checked$mean[shared_ov]; covariance <- checked$covariance[shared_ov,shared_ov,drop=FALSE]
      observed <- definitions$.sem_cv_numeric_matrix(data[test,,drop=FALSE],shared_ov)
      # PP1–PP3 no entran en X. Sólo se usan después para evaluar el error.
      predictor_names <- setdiff(shared_ov,pp)
      predicted <- definitions$.sem_cv_conditional_pp(observed[,predictor_names,drop=FALSE],mu,covariance,pp)
      errors <- rowMeans(predicted)-rowMeans(observed[,pp,drop=FALSE])
      score <- definitions$.sem_cv_gaussian_logscore(observed,mu,covariance)
      stopifnot(all(is.finite(c(errors,score))))
      # No se exportan predicciones personales. El cache privado conserva los
      # momentos y filas suficientes para recalcularlas exactamente.
      index <- index+1L
      records[[index]] <- data.frame(model=id,repetition=r,fold=k,
        n_train=length(train),n_test=length(test),sum_squared_error=sum(errors^2),
        sum_absolute_error=sum(abs(errors)),sum_joint_logscore=sum(score),
        RMSE_PP_mean=sqrt(mean(errors^2)),MAE_PP_mean=mean(abs(errors)),
        mean_joint_logscore=mean(score),admissible=checked$admissible,
        warnings=paste(attr(fit,"captured_warnings"),collapse=" | "))
      writeLines(sprintf("T1 CV: %d/175; %s r%d f%d; nuevos=%d cache=%d",index,id,r,k,fitted,reused),
                 file.path(config$salida,"T1_CV_avance.txt"))
    }
    cat(sprintf("[T1-CV] Repetición %d/5, %s: 5/5 pliegues completos\n",r,id))
    flush.console()
  }
  folds <- do.call(rbind,records)
  repetitions <- do.call(rbind,lapply(1:5,function(r)do.call(rbind,lapply(names(definitions$specs),function(id) {
    z <- folds[folds$repetition==r & folds$model==id,]
    stopifnot(nrow(z)==5L,sum(z$n_test)==653L)
    data.frame(model=id,repetition=r,n_test=653L,
      RMSE_PP_mean=sqrt(sum(z$sum_squared_error)/653),
      MAE_PP_mean=sum(z$sum_absolute_error)/653,
      mean_joint_logscore=sum(z$sum_joint_logscore)/653)
  }))))
  overall <- do.call(rbind,lapply(names(definitions$specs),function(id) {
    z <- repetitions[repetitions$model==id,]
    data.frame(model=id,repetitions=5L,RMSE_PP_mean=mean(z$RMSE_PP_mean),
      MAE_PP_mean=mean(z$MAE_PP_mean),mean_joint_logscore=mean(z$mean_joint_logscore),
      RMSE_partition_SD=sd(z$RMSE_PP_mean))
  }))
  .tab_modelos_exportar(folds,"T1_CV_pliegues_recalculados.csv",config$salida)
  .tab_modelos_exportar(repetitions,"T1_CV_repeticiones_recalculadas.csv",config$salida)
  .tab_modelos_exportar(overall,"T1_CV_resumen_recalculado.csv",config$salida)
  list(overall=overall,new_fits=fitted,reused_fits=reused,complete=TRUE,cache_dir=cache_dir)
}

#### Función pública para el archivo principal de RStudio ####

# [ENTRADA / T1, T2,S1–S4] config identifica fuentes, salida y registro.
# Guardas: existen las fuentes, son objetos lavaan y no cambian sus hashes.
# El principal fija prueba_sin_cv=TRUE; cambiarlo ejecutaría el legado histórico.
generar_tablas_modelos <- function(config) {
  stopifnot(is.list(config),all(c("proyecto","salida","registro")%in%names(config)))
  # TRUE produce realmente sólo S1–S4: no calcula T1/T2 ni ejecuta CV.
  # Si se omite esta opción, continúa la reproducción completa original.
  solo_suplemento <- isTRUE(config$solo_suplemento)
  config$proyecto <- normalizePath(config$proyecto,mustWork=TRUE)
  if (is.null(config$privado)) config$privado <- file.path(config$salida,"privados")
  dir.create(config$salida,recursive=TRUE,showWarnings=FALSE)
  dir.create(config$privado,recursive=TRUE,showWarnings=FALSE)
  .libPaths(c(file.path(config$proyecto,"revision-revista-2026-10-04/auditoria-codigos/R-library"),
              path.expand("~/Library/R/arm64/4.6/library"),.libPaths()))
  suppressPackageStartupMessages({library(lavaan);library(semTools);library(digest);library(jsonlite)})
# [FUENTES / T1, T2,S1,S3,S4] Siete especificaciones M0–M6 × CFA/SEM
# × MLR/WLSMV. La estructura/estimador vienen del objeto, no se redefinen aquí.
  ids <- c("M0_original","M1_ASP_unico","M2_puentes2","M3_puentes4",
           "M4_S1_refSC","M5_S1_refSBL","M6_ESEM_target")
  input <- file.path(config$proyecto,"sem-alternativas-2026-10-08/resultados")
  # Los ajustes completos ya existen: se vuelven a extraer sus estadísticas,
  # sin estimarlos. El suplemento necesita 17 objetos y el modo completo 28.
  paths <- if (solo_suplemento) {
    c(file.path(input,paste0(ids,"_CFA_MLR.rds")),
      file.path(input,paste0(ids,"_SEM_WLSMV.rds")),
      file.path(input,paste0(c("M0_original","M4_S1_refSC","M5_S1_refSBL"),"_CFA_WLSMV.rds")))
  } else {
    unlist(lapply(ids,function(id)outer(c("CFA","SEM"),c("MLR","WLSMV"),
      function(kind,est)file.path(input,paste0(id,"_",kind,"_",est,".rds")))))
  }
  names(paths) <- NULL
  dictionary_file <- file.path(config$proyecto,"revision-revisores-2026-10-06/conceptual/instrument_literal_evidence.csv")
  data_file <- file.path(config$proyecto,"revision-integral-2026-10-04/analysis/outputs/canonical-full-20261004-02/01_audit/analytic_data.rds")
  p3_file <- file.path(config$proyecto,"revision-integral-2026-10-04/sensitivities/P3_country_partial_demographics.rds")
  source_script <- file.path(config$proyecto,"sem-alternativas-2026-10-08/codigo/comparar_modelos_SEM.R")
  helper_file <- file.path(config$proyecto,"sem-alternativas-2026-10-08/codigo/cross_validate_models.R")
  protected <- c(paths,dictionary_file,data_file,p3_file,source_script,helper_file)
  stopifnot(all(file.exists(protected)))
  before <- vapply(protected,.tab_modelos_sha,character(1))
  fits <- setNames(lapply(paths,readRDS),sub("\\.rds$","",basename(paths)))
  stopifnot(all(vapply(fits,inherits,logical(1),what="lavaan")))
# [AUXILIAR] Identifica un ajuste por modelo, tipo y estimador. Leer un RDS
# lavaan no ajusta otro modelo; los modelos guardados contienen datos privados.
  getfit <- function(id,kind,est)fits[[paste(id,kind,est,sep="_")]]
  # fitMeasures puede avisar del SEM ordinal M6 ya problemático. Se registra
  # el aviso sin rellenar su CFI/RMSEA ausente con otro índice diferente.
  fitmeasure <- function(fit,metrics) {
    captured <- character()
    ans <- withCallingHandlers(lavaan::fitMeasures(fit,metrics),warning=function(w) {
      captured <<- c(captured,conditionMessage(w));invokeRestart("muffleWarning")
    })
    attr(ans,"warnings") <- captured
    ans
  }

  # En modo suplemento ni siquiera se entra en los cálculos del artículo.
  # Estos ceros documentan que la CV no fue solicitada, no una CV fallida.
  cv <- list(new_fits=0L,reused_fits=0L,complete=FALSE)
  if (!solo_suplemento) {
  #### [T1] Artículo: ajuste estructural y predicción interna ####
  labels <- c("M0 Original","M1 One attitude","M2 Two bridges","M3 Four bridges",
              "M4 SC reference S1","M5 SBL reference S1","M6 ESEM")
  numeric_T1 <- do.call(rbind,lapply(ids,function(id) {
    fm <- fitmeasure(getfit(id,"SEM","MLR"),c("cfi.robust","rmsea.robust","srmr","npar"))
    data.frame(model=id,CFI_robust=unname(fm["cfi.robust"]),RMSEA_robust=unname(fm["rmsea.robust"]),
      SRMR=unname(fm["srmr"]),npar=unname(fm["npar"]))
  }))
  .tab_modelos_exportar(numeric_T1,"T1_ajuste_recalculado.csv",config$salida)
# [T1 ACTUAL] La rama siguiente evita los175 ajustes históricos. El principal
# elimina sus dos columnas vacías; no interpreta esos NA como predicción V3.
  if (isTRUE(config$prueba_sin_cv)) {
    # Sólo sirve para probar extracción de los ajustes completos sin ejecutar
    # los 175 entrenamientos. No se declara T1 verificada en este modo.
    cv <- list(overall=data.frame(model=ids,mean_joint_logscore=NA_real_,RMSE_PP_mean=NA_real_),
               new_fits=0L,reused_fits=0L,complete=FALSE)
  } else {
    definitions <- .tab_modelos_extraer_funciones_originales(config$proyecto)
    p3 <- readRDS(p3_file); original <- as.data.frame(readRDS(data_file))
    data <- as.data.frame(lavaan::lavInspect(p3,"data"))
    common <- intersect(names(data),names(original))
    stopifnot(nrow(data)==653L,max(abs(as.matrix(data[,common])-as.matrix(original[,common])))==0)
    data$Pais <- original$Pais
    cv <- .tab_modelos_cv(data,definitions,config,before)
  }
  T1 <- .tab_modelos_matriz(c("Model","CFI robust","RMSEA robust","SRMR","Joint log score","PP RMSE","npar"),
    lapply(seq_along(ids),function(i) {
      row <- numeric_T1[i,]; pred <- cv$overall[cv$overall$model==ids[i],]
      c(labels[i],.tab_modelos_num(row$CFI_robust),.tab_modelos_num(row$RMSEA_robust),
        .tab_modelos_num(row$SRMR),.tab_modelos_num(pred$mean_joint_logscore),
        .tab_modelos_num(pred$RMSE_PP_mean,4),as.character(as.integer(row$npar)))
    }))

  #### [T2] Artículo: coeficientes estandarizados e intervalos robustos ####
  # Holm se aplica a los 21 p NO estandarizados de cada modelo, antes de
  # seleccionar seis filas para la tabla. Los IC de beta siguen siendo
  # puntuales: que no crucen cero no implica superar la corrección Holm.
  focal_values <- c("Benevolencia","Universalismo","Poder","Logro")
  selected_ids <- c("M0_original","M2_puentes2","M6_ESEM_target")
  all_T2 <- list()
  for (id in selected_ids) {
    f <- getfit(id,"SEM","MLR")
    pe <- lavaan::parameterEstimates(f,ci=TRUE)
    ss <- lavaan::standardizedSolution(f,type="std.all",ci=TRUE)
    focal <- pe[pe$op=="~" & pe$lhs%in%c(focal_values,"SC","SBL","PSB") &
                  pe$rhs%in%c(focal_values,"SC","SBL","CSR1_bin"),]
    stopifnot(nrow(focal)==21L)
    focal$p_Holm_21 <- p.adjust(focal$pvalue,"holm")
    m <- match(paste(focal$lhs,focal$rhs),paste(ss$lhs,ss$rhs))
    focal$beta <- ss$est.std[m];focal$beta_CI_lower <- ss$ci.lower[m];focal$beta_CI_upper <- ss$ci.upper[m]
    focal$model <- id
    # ESEM añade una columna 'efa' a parameterEstimates; se conserva un
    # esquema común explícito para unir los tres modelos sin perder cifras.
    all_T2[[id]] <- as.data.frame(focal)[,c("lhs","op","rhs","label","est","se","z","pvalue",
      "ci.lower","ci.upper","p_Holm_21","beta","beta_CI_lower","beta_CI_upper","model")]
  }
  numeric_T2 <- do.call(rbind,all_T2)
  .tab_modelos_exportar(numeric_T2,"T2_familias_21_recalculadas.csv",config$salida)
  pairs <- data.frame(lhs=c("SC","SBL","PSB","PSB","SBL","PSB"),
                      rhs=c("Universalismo","Universalismo","SC","SBL","Poder","Benevolencia"))
  association <- c("Universalism to SC","Universalism to SBL","SC to PSB","SBL to PSB","Power to SBL","Benevolence to PSB")
  T2 <- .tab_modelos_matriz(c("Association","M0 Original","M2 Two bridges","M6 ESEM"),lapply(1:6,function(i) {
    cells <- vapply(selected_ids,function(id) {
      row <- numeric_T2[numeric_T2$model==id & numeric_T2$lhs==pairs$lhs[i] & numeric_T2$rhs==pairs$rhs[i],]
      stopifnot(nrow(row)==1L)
      paste0(.tab_modelos_num(row$beta),if(row$p_Holm_21<.05)"*"else"","\n[",
             .tab_modelos_num(row$beta_CI_lower),", ",.tab_modelos_num(row$beta_CI_upper),"]")
    },character(1))
    c(association[i],cells)
  }))
  }

  #### [S1] Suplemento: ajuste de los siete modelos de medición ####
  numeric_S1 <- do.call(rbind,lapply(ids,function(id) {
    f <- getfit(id,"CFA","MLR")
    fm <- fitmeasure(f,c("cfi.robust","tli.robust","rmsea.robust","srmr"))
    correlation <- lavaan::lavInspect(f,"cor.lv")
    r <- if(id=="M1_ASP_unico") NA_real_ else if(id%in%c("M4_S1_refSC","M5_S1_refSBL"))0 else correlation["SC","SBL"]
    data.frame(model=id,CFI_robust=unname(fm[1]),TLI_robust=unname(fm[2]),RMSEA_robust=unname(fm[3]),
      SRMR=unname(fm[4]),SC_SBL_r=r)
  }))
  .tab_modelos_exportar(numeric_S1,"S1_CFA_indices_recalculados.csv",config$salida)
  S1 <- .tab_modelos_matriz(c("Model","CFI robust","TLI robust","RMSEA robust","SRMR","SC–SBL r"),
    lapply(seq_along(ids),function(i) {
      z <- numeric_S1[i,]
      c(paste0("M",i-1L),.tab_modelos_num(z$CFI_robust),.tab_modelos_num(z$TLI_robust),
        .tab_modelos_num(z$RMSEA_robust),.tab_modelos_num(z$SRMR),
        if(ids[i]%in%c("M4_S1_refSC","M5_S1_refSBL"))"0 (fixed)"else .tab_modelos_num(z$SC_SBL_r))
    }))

  #### [S2] Suplemento: contenido literal y cargas de los dos ítems puente ####
  dictionary <- utils::read.csv(dictionary_file,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding="UTF-8")
  # El texto procede del diccionario literal SPSS cotejado anteriormente.
  # Las cargas y sus IC se extraen AHORA del objeto CFA M2 MLR, no del Word.
  loads <- lavaan::standardizedSolution(getfit("M2_puentes2","CFA","MLR"),type="std.all",ci=TRUE)
  loads <- loads[loads$op=="=~" & loads$rhs%in%c("APS2","APS13"),]
  .tab_modelos_exportar(loads,"S2_cargas_puentes_recalculadas.csv",config$salida)
  loading_cell <- function(item,factor) {
    z <- loads[loads$rhs==item & loads$lhs==factor,];stopifnot(nrow(z)==1L)
    paste0(factor,": ",.tab_modelos_num(z$est.std)," [",.tab_modelos_num(z$ci.lower),", ",.tab_modelos_num(z$ci.upper),"]")
  }
  S2 <- .tab_modelos_matriz(c("Item","Literal archived Spanish wording","Primary loading [95% CI]","Bridge loading [95% CI]"),
    lapply(c("APS2","APS13"),function(item) {
      label <- dictionary$variable_label[dictionary$variable==item];stopifnot(length(label)==1L)
      primary <- if(item=="APS2")"SC"else"SBL";bridge <- if(item=="APS2")"SBL"else"SC"
      c(item,label,loading_cell(item,primary),loading_cell(item,bridge))
    }))

  #### [S3] Suplemento: AVE simple y fiabilidad del componente específico ####
  # AVE = suma de varianzas explicadas / suma de varianzas de los indicadores.
  # semTools::AVE(obs.var=TRUE) usa la referencia observada; la versión WLSMV
  # usa respuestas ordinales latentes. Las dos escalas no son intercambiables.
  # En S–1, omega HS se obtiene del componente específico ortogonal: lambda
  # específica sumada al cuadrado por su varianza / varianza total del puntaje.
  # Esto conserva la definición matricial ARCHIVADA, no fiabilidad de todo SBL.
  numeric_S3 <- list(); index <- 0L
  for (est in c("MLR","WLSMV")) {
    f <- getfit("M0_original","CFA",est)
    ave <- semTools::AVE(f,obs.var=TRUE)
    for (factor in c("SC","SBL","PSB")) {
      index <- index+1L
      numeric_S3[[index]] <- data.frame(model="M0_original",domain=factor,quantity="AVE",estimator=est,value=unname(ave[factor]))
    }
    for (id in c("M4_S1_refSC","M5_S1_refSBL")) {
      f <- getfit(id,"CFA",est); m <- lavaan::lavInspect(f,"est")
      factor <- if(id=="M4_S1_refSC")"SBLspec"else"SCspec"
      domain_items <- if(factor=="SBLspec") c("APS11","APS12","APS13","APS15","APS16","APS18","APS19","APS20","APS21") else
        c("APS1","APS2","APS3","APS4","APS5","APS6","APS7","APS9","APS10","APS14")
      S <- m$lambda%*%m$psi%*%t(m$lambda)+m$theta
      omega <- sum(m$lambda[domain_items,factor])^2*m$psi[factor,factor]/sum(S[domain_items,domain_items])
      index <- index+1L
      numeric_S3[[index]] <- data.frame(model=id,domain=factor,quantity="omega_HS",estimator=est,value=omega)
    }
  }
  numeric_S3 <- do.call(rbind,numeric_S3)
  .tab_modelos_exportar(numeric_S3,"S3_AVE_omega_recalculadas.csv",config$salida)
  S3_specs <- data.frame(model=c(rep("M0_original",3),"M4_S1_refSC","M5_S1_refSBL"),domain=c("SC","SBL","PSB","SBLspec","SCspec"))
  S3 <- .tab_modelos_matriz(c("Model/domain","Quantity","MLR","WLSMV latent-response sensitivity"),lapply(1:5,function(i) {
    d <- S3_specs[i,];z <- numeric_S3[numeric_S3$model==d$model & numeric_S3$domain==d$domain,]
    c(paste0(if(i<=3)"M0"else if(i==4)"M4"else"M5"," / ",d$domain),
      if(i<=3)"AVE"else"ωHS",.tab_modelos_num(z$value[z$estimator=="MLR"]),.tab_modelos_num(z$value[z$estimator=="WLSMV"]))
  }))

  #### [S4] Suplemento: sensibilidad estructural ordinal ####
  numeric_S4 <- do.call(rbind,lapply(ids,function(id) {
    fm <- fitmeasure(getfit(id,"SEM","WLSMV"),c("cfi.scaled","rmsea.scaled","srmr"))
    data.frame(model=id,CFI_scaled=unname(fm[1]),RMSEA_scaled=unname(fm[2]),SRMR=unname(fm[3]),
               warnings=paste(attr(fm,"warnings"),collapse=" | "))
  }))
  .tab_modelos_exportar(numeric_S4,"S4_ordinal_indices_recalculados.csv",config$salida)
  S4 <- .tab_modelos_matriz(c("Model","CFI scaled","RMSEA scaled","SRMR"),lapply(seq_along(ids),function(i) {
    z <- numeric_S4[i,];c(paste0("M",i-1L),.tab_modelos_num(z$CFI_scaled),.tab_modelos_num(z$RMSEA_scaled),.tab_modelos_num(z$SRMR))
  }))
  tables <- if (solo_suplemento) list(S1=S1,S2=S2,S3=S3,S4=S4) else
    list(T1=T1,T2=T2,S1=S1,S2=S2,S3=S3,S4=S4)
  for (name in names(tables)) utils::write.table(tables[[name]],file.path(config$salida,paste0(name,"_tabla_formateada.csv")),
    sep=",",quote=TRUE,row.names=FALSE,col.names=FALSE,fileEncoding="UTF-8")

  #### Verificación final y protección de las fuentes ####
  after <- vapply(protected,.tab_modelos_sha,character(1));stopifnot(identical(before,after))
  # Sólo AQUÍ se lee el contrato extraído de los Word. No alimenta cálculos.
  contract_file <- file.path(config$registro,"tablas_word_esperadas.json")
  matches <- list()
  if (file.exists(contract_file)) {
    contract <- jsonlite::read_json(contract_file,simplifyVector=FALSE)
    for (name in names(tables)) {
      expected <- Filter(function(entry)identical(entry$key,name),contract)
      stopifnot(length(expected)==1L)
      wanted <- do.call(rbind,lapply(expected[[1]]$cells,unlist))
      dimnames(wanted) <- NULL;storage.mode(wanted) <- "character"
      # El modo de prueba omite sólo las dos columnas CV de T1.
      compare <- if(name=="T1" && !cv$complete)c(1:4,7)else seq_len(ncol(tables[[name]]))
      matches[[name]] <- identical(tables[[name]][,compare,drop=FALSE],wanted[,compare,drop=FALSE])
    }
    if (!all(unlist(matches))) stop("No coinciden todas las tablas recalculadas con el contrato Word: ",paste(names(matches)[!unlist(matches)],collapse=", "))
  }
  provenance <- list(sources_unchanged=TRUE,reused_full_sample_objects=length(fits),
    CV_requested=!solo_suplemento && !isTRUE(config$prueba_sin_cv),
    new_CV_training_fits=cv$new_fits,reused_CV_training_fits=cv$reused_fits,
    CV_complete=cv$complete,table_word_matches=matches,
    source_sha256=as.list(after),
    interpretation=if(solo_suplemento)
      "Recalculation of supplementary statistics from archived fits; no model estimation or CV requested" else
      "Recalculation of full-fit statistics plus reproduction of original CV; no new models, cases, item choices or selection",
    mode=if(solo_suplemento)"Supplement only: S1-S4; T1/T2 and 175 CV training fits not requested" else
      if(cv$complete)"full reproduction"else"Full-fit extraction smoke test; T1 CV deliberately not calculated")
  provenance_file <- if(solo_suplemento)"S1_S4_proveniencia.json"else"T1_T2_S1_S4_proveniencia.json"
  jsonlite::write_json(provenance,file.path(config$salida,provenance_file),
                       pretty=TRUE,auto_unbox=TRUE,digits=NA)
  attr(tables,"proveniencia") <- provenance
  tables
}
