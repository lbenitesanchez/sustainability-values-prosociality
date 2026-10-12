#!/usr/bin/env Rscript
# Public entry point for the completed M2 revision.
# Default: check exported files ONLY. This script never starts model fitting.
# See documentation/README_RERUN.md for the scope of each mode.
# source() only defines functions; call run_M2_package() explicitly in RStudio.

.M2_THIS_FILE <- local({
  f <- Filter(function(x) is.character(x) && length(x)==1L && nzchar(x),
              lapply(sys.frames(), function(fr) fr$ofile))
  if(length(f)) normalizePath(tail(f,1L)[[1L]], mustWork=TRUE) else {
    a <- grep('^--file=', commandArgs(FALSE), value=TRUE)
    if(!length(a)) stop('Use source() or Rscript to load R/00_master.R.')
    normalizePath(gsub('~+~',' ',sub('^--file=','',a[1L]),fixed=TRUE),mustWork=TRUE)
  }
})
.M2_PACKAGE <- dirname(dirname(.M2_THIS_FILE))

M2_expected_tables <- function() c(paste0('T',1:4),paste0('S',1:29))
M2_write_csv <- function(x,path) utils::write.csv(x,path,row.names=FALSE,na='',fileEncoding='UTF-8')
M2_record_environment <- function(out,mode) {
  writeLines(c(paste('Mode:',mode),paste('UTC:',format(Sys.time(),tz='UTC',usetz=TRUE)),
    paste('RNG:',paste(RNGkind(),collapse=' / ')),
    'The master did not set an empirical random seed or run model fitting.',
    capture.output(sessionInfo())),file.path(out,'sessionInfo.txt'),useBytes=TRUE)
}
M2_fresh_output <- function(output,root) {
  if(is.null(output)) output <- tempfile('M2_package_')
  output <- normalizePath(output,mustWork=FALSE)
  # Private/synthetic generated files must not silently enter the public package.
  if(identical(output,root) || startsWith(paste0(output,'/'),paste0(root,'/')))
    stop('Choose an output directory OUTSIDE the public package.')
  if(file.exists(output))stop('Output already exists; choose a new directory: ',output)
  if(!dir.create(output,recursive=TRUE))stop('Cannot create output directory: ',output)
  normalizePath(output,mustWork=TRUE)
}
M2_check_sources <- function(root,out) {
  index <- file.path(root,'documentation','source_index.csv')
  if(!file.exists(index))stop('Missing documentation/source_index.csv')
  z <- utils::read.csv(index,check.names=FALSE,stringsAsFactors=FALSE)
  stopifnot(all(c('public_path','md5','sha256')%in%names(z)))
  paths <- file.path(root,z$public_path)
  exists <- file.exists(paths)
  actual <- rep(NA_character_,length(paths))
  actual[exists] <- unname(tools::md5sum(paths[exists]))
  status <- data.frame(file=z$public_path,exists=exists,expected_md5=z$md5,
                       actual_md5=actual,identical=exists & !is.na(actual) & actual==z$md5)
  M2_write_csv(status,file.path(out,'source_hash_check.csv'))
  if(!all(status$identical))stop('One or more preserved source files are missing/changed; see source_hash_check.csv.')
  # SHA-256 provenance is preserved in the index. Base-R MD5 detects byte changes
  # here without making the default public-export check depend on a package.
  invisible(status)
}
M2_check_exports <- function(root,out) {
  mapfile <- file.path(root,'documentation','table_map.csv')
  if(!file.exists(mapfile))stop('Missing documentation/table_map.csv; exports have not been fully packaged.')
  z <- utils::read.csv(mapfile,check.names=FALSE,stringsAsFactors=FALSE)
  first_col <- function(candidates){a<-intersect(candidates,names(z));if(length(a))a[1]else NA_character_}
  idc <- first_col(c('table_id','table','tabla','id'))
  fc <- first_col(c('table_file','file','csv','output_file'))
  if(is.na(idc)||is.na(fc))stop('Map must have table_id and table_file columns (or documented equivalents).')
  id <- as.character(z[[idc]]); files <- as.character(z[[fc]])
  # The map may also include Figure 1. Validate figures separately below.
  keep <- id %in% M2_expected_tables()
  missing_ids <- setdiff(M2_expected_tables(),unique(id[keep]))
  unexpected <- setdiff(unique(id[grepl('^[TS][0-9]+$',id)]),M2_expected_tables())
  if(length(missing_ids)||length(unexpected))stop('Expected exactly T1-T4 and S1-S29. Missing: ',
      paste(missing_ids,collapse=', '),'; unexpected: ',paste(unexpected,collapse=', '))
  if(any(grepl('^/|^[A-Za-z]:',files[keep])))stop('Public table paths must be relative to package root.')
  paths <- file.path(root,files[keep]); tabids <- id[keep]
  rows <- lapply(seq_along(paths),function(i){
    p<-paths[i];ok<-file.exists(p);err<-'';nr<-nc<-NA_integer_
    if(ok){a<-tryCatch(utils::read.csv(p,header=FALSE,check.names=FALSE,
             colClasses='character',na.strings=NULL,stringsAsFactors=FALSE,fileEncoding='UTF-8'),error=function(e)e)
      if(inherits(a,'error'))err<-conditionMessage(a) else {nr<-nrow(a);nc<-ncol(a);ok<-nr>=2L&&nc>=2L}}
    data.frame(table_id=tabids[i],table_file=files[keep][i],exists=file.exists(p),
      readable_with_header_row=ok&&err=='',rows_including_header=nr,columns=nc,error=err)
  })
  status<-do.call(rbind,rows);M2_write_csv(status,file.path(out,'export_check.csv'))
  if(!all(status$readable_with_header_row))stop('Some mapped table/panel exports are missing or unreadable.')
  # The map must point to code actually shipped, not obsolete source locations.
  gc <- intersect(c('generator','script','code_file'),names(z))
  if(length(gc)){
    gp <- unique(as.character(z[[gc[1L]]]))
    gp <- gp[!is.na(gp) & nzchar(gp)]
    status_code <- data.frame(generator=gp,exists=file.exists(file.path(root,gp)))
    M2_write_csv(status_code,file.path(out,'mapped_code_check.csv'))
    if(!all(status_code$exists))stop('Map references missing generator files: ',paste(gp[!status_code$exists],collapse=', '))
  }
  fig <- files[grepl('^F(igure)?[ _]?1$',id,ignore.case=TRUE)]
  if(length(fig)&&!all(file.exists(file.path(root,fig))))stop('A mapped Figure 1 file is missing.')
  writeLines(c('PASS: all 33 table IDs are represented by readable exported table/panel CSVs.',
    paste('Mapped table-panel files:',nrow(status)),
    'This validates file completeness and preserved code; it does not re-estimate models,',
    'establish numerical equivalence to the Word documents, or reproduce from participant data.'),
    file.path(out,'STATUS.txt'))
  cat('PASS: exported files for 33 tables; no statistical analysis run.\n')
  invisible(status)
}
M2_private_extract <- function(root,out,private_project,private_reference,private_word,r_library=NULL) {
  required<-c(project=private_project,reference=private_reference,word=private_word)
  if(length(required)!=3L || any(!nzchar(required)) || any(!dir.exists(required)))
    stop('private_extract requires existing --private-project, --private-reference and --private-word directories.')
  required <- vapply(required,normalizePath,character(1),mustWork=TRUE)
  if(any(startsWith(paste0(out,'/'),paste0(required,'/'))))stop('Output must not be inside any private source directory.')
  if(!is.null(r_library)){
    if(!dir.exists(r_library))stop('Requested R library does not exist.')
    .libPaths(c(r_library,.libPaths()))
  }
  deps<-c('lavaan','semTools','jsonlite','digest','haven')
  absent<-deps[!vapply(deps,requireNamespace,logical(1),quietly=TRUE)]
  if(length(absent))stop('Required packages not installed in configured libraries: ',paste(absent,collapse=', '))
  # 1. Historical 24-table contract: extraction/formatting only; old CV is off.
  generator<-file.path(root,'R','legacy_table_generator','01_GENERAR_TODAS_LAS_TABLAS_V3.R')
  # source() supplies the ofile frame used by the preserved legacy loader.
  # sys.source() lacks that frame and failed in clean RStudio background jobs
  # (which also have no Rscript --file argument). Statistical code is unchanged.
  env<-new.env(parent=globalenv());source(generator,local=env,encoding='UTF-8')
  legacy<-file.path(out,'legacy_24_table_extraction')
  env$generar_todas_las_tablas_V3(proyecto_fuente=required[['project']],salida=legacy,
    raiz_auditoria=dirname(required[['project']]),proyecto_referencia=required[['reference']],
    directorio_word=required[['word']],proyecto_codigo=required[['project']],
    directorio_modulos=file.path(root,'R','legacy_table_generator','R'),fallar_si_diferencias=TRUE)
  # 2. Additional M2 extraction and sample summaries. The third argument names
  # the actual legacy item-path export, removing the original machine layout.
  extra_script<-file.path(root,'R','extensions','01_extract_complete_M2.R')
  previous_item_paths<-file.path(legacy,'resultados_base','02_S7_S8_caminos_numericos_IC_Holm18.csv')
  if(!file.exists(previous_item_paths))stop('Legacy item-path reference was not generated: ',previous_item_paths)
  extra<-file.path(out,'complete_M2_extraction')
  if(!is.null(r_library))Sys.setenv(M2_R_LIBRARY=normalizePath(r_library))
  cmd<-file.path(R.home('bin'),'Rscript')
  child_args<-c('--vanilla',extra_script,required[['project']],extra,previous_item_paths)
  rc<-system2(cmd,vapply(child_args,shQuote,character(1)),
    stdout=file.path(out,'complete_M2_execution.log'),stderr=file.path(out,'complete_M2_execution.log'))
  if(!identical(rc,0L))stop('Additional extraction failed; see complete_M2_execution.log.')
  # 3. Publication numbering and final labels. No numerical cell is changed.
  # Keep the historical extractor and its frozen comparisons byte-for-byte.
  # Only final T1/T2/S22 text is aligned with the current Word presentation:
  # M0 Baseline, abbreviated CSR1, and real LF breaks in institution headers.
  tableout<-file.path(out,'tables');dir.create(tableout)
  presentation_changes<-list()
  publication_labels<-function(id,path){
    if(!id%in%c('T1','T2','S22'))return(invisible(NULL))
    x<-unname(as.matrix(utils::read.csv(path,header=FALSE,check.names=FALSE,
        colClasses='character',na.strings=NULL,stringsAsFactors=FALSE,fileEncoding='UTF-8')))
    change<-function(row,column,old,new){
      before<-x[row,column]
      if(!before%in%c(old,new))stop('Unexpected publication label in ',id,
          ' [',row,',',column,']: ',before)
      x[row,column]<<-new
      presentation_changes[[length(presentation_changes)+1L]]<<-data.frame(
        table_id=id,row=row,column=column,before=before,after=new,
        changed=before!=new,operation='Text presentation only; no numerical changes')
    }
    if(id=='T1')change(2L,1L,'M0 Original','M0 Baseline')
    if(id=='T2')change(1L,2L,'M0 Original','M0 Baseline')
    if(id=='S22'){
      countries<-c('Argentina','Colombia','Peru','Mexico')
      for(i in seq_along(countries))change(1L,i+3L,
        paste('Institution in',countries[i]),paste0('Institution in\n',countries[i]))
      for(row in 17:19)change(row,1L,'CSR1: attended ethics/CSR course at university',
        'CSR1 course attendance')
    }
    utils::write.table(x,path,sep=',',quote=TRUE,qmethod='double',row.names=FALSE,
      col.names=FALSE,fileEncoding='UTF-8',na='NA')
    invisible(NULL)
  }
  map<-list();add<-function(id,panel,src,dest){
    if(!file.exists(src))stop('Missing generated export: ',src)
    if(!file.copy(src,file.path(tableout,dest)))stop('Cannot copy generated export.')
    publication_labels(id,file.path(tableout,dest))
    map[[length(map)+1L]]<<-data.frame(table_id=id,panel=panel,table_file=paste0('tables/',dest),
             source_file=src,operation='Existing-model extraction / aggregate formatting / final publication labels; no refitting')
  }
  for(id in c('T1','T2',paste0('S',1:21)))
    add(id,'',file.path(legacy,'tablas_finales',paste0(id,'_verificada.csv')),paste0(id,'.csv'))
  add('T4','',file.path(legacy,'tablas_finales','T3_verificada.csv'),'T4.csv')
  add('T3','',file.path(extra,'T3_items_MLR_formatted.csv'),'T3.csv')
  newfiles<-c(S22='S22_descriptives_formatted.csv',S23='S23_loadings_formatted.csv',
    S24='S24_loadings_formatted.csv',S25='S25_loadings_formatted.csv',S26='S26_loadings_formatted.csv',
    S27='S27_correlations_formatted.csv')
  for(id in names(newfiles))add(id,'',file.path(extra,newfiles[[id]]),paste0(id,'.csv'))
  for(est in c('MLR','WLSMV'))for(part in c('focal','controls')){
    panel<-paste0(if(est=='MLR')'A'else'B',if(part=='focal')'1'else'2')
    add('S28',panel,file.path(extra,paste0('S28_',est,'_',part,'_formatted.csv')),paste0('S28',panel,'.csv'))
  }
  for(est in c('MLR','WLSMV')){
    panel<-if(est=='MLR')'A'else'B'
    add('S29',panel,file.path(extra,paste0('S29_',est,'_indirects_formatted.csv')),paste0('S29',panel,'.csv'))
  }
  m<-do.call(rbind,map);stopifnot(setequal(unique(m$table_id),M2_expected_tables()),nrow(m)==37L)
  M2_write_csv(m,file.path(out,'generated_table_map.csv'))
  p<-do.call(rbind,presentation_changes);stopifnot(nrow(p)==9L)
  M2_write_csv(p,file.path(out,'publication_presentation_changes.csv'))
  writeLines(c('PASS: 33 canonical tables, 37 panel CSV files generated.',
    'Nine final label/header cells use the current Word presentation; numerical cells are unchanged.',
    'The original 24-table extractor used saved fits and aggregate exports; it did not refit models.',
    'Added tables used saved fits and new descriptive aggregation, not new statistical model fitting.',
    'This local output may contain private path information and is NOT automatically a public release.',
    'Participant-data reproduction and full historical training were not performed by this mode.'),file.path(out,'STATUS.txt'))
  cat('PASS: private extraction completed; no model training performed.\n')
  invisible(m)
}
run_M2_package <- function(mode='check_exports',output=NULL,package_root=.M2_PACKAGE,
                           private_project=NULL,private_reference=NULL,private_word=NULL,r_library=NULL) {
  mode<-match.arg(mode,c('check_exports','synthetic','private_extract'))
  root<-normalizePath(package_root,mustWork=TRUE)
  out<-M2_fresh_output(output,root)
  M2_record_environment(out,mode)
  M2_check_sources(root,out)
  if(mode=='check_exports')M2_check_exports(root,out)
  if(mode=='synthetic'){
    env<-new.env(parent=globalenv())
    sys.source(file.path(root,'R','synthetic','02_SIMULAR_DATOS_ARTIFICIALES.R'),envir=env)
    # Source defines functions only; call the preserved generator explicitly.
    z<-env$exportar_artificial(file.path(out,'synthetic_data'))
    stopifnot(nrow(z$data)==800L,ncol(z$data)==49L,all(z$tests$cumple))
    writeLines(c('PASS: 800 entirely artificial cases, 49 columns, generator checks passed.',
      'No participant records, fitted models, empirical frequencies or correlations were read.',
      'These data DO NOT reproduce any empirical table in the manuscript or supplement.'),file.path(out,'STATUS.txt'))
    cat('PASS: synthetic demonstration only; no empirical analysis run.\n')
  }
  if(mode=='private_extract')M2_private_extract(root,out,private_project,private_reference,private_word,r_library)
  cat('Output: ',out,'\n',sep='');invisible(out)
}
if(sys.nframe()==0L){
  args<-commandArgs(trailingOnly=TRUE)
  if(any(args%in%c('--help','-h'))){
    cat('Rscript --vanilla R/00_master.R [--mode=check_exports|synthetic|private_extract] [--output=/new/path]\n',
        'Private mode requires --private-project=... --private-reference=... --private-word=...\n',
        'Optional --r-library=/existing/library. No training/refitting mode is provided.\n',sep='')
  }else{
    valid<-c('mode','output','private-project','private-reference','private-word','r-library')
    getarg<-function(key,default=NULL){z<-args[startsWith(args,paste0('--',key,'='))];if(length(z)>1)stop('Duplicate argument: ',key);if(length(z))substring(z,nchar(key)+4L)else default}
    if(any(!grepl('^--[^=]+=',args)))stop('Use --name=value arguments; see --help.')
    keys<-sub('^--([^=]+)=.*$','\\1',args);if(any(!keys%in%valid))stop('Unknown argument: ',paste(keys[!keys%in%valid],collapse=', '))
    run_M2_package(mode=getarg('mode','check_exports'),output=getarg('output'),
      private_project=getarg('private-project'),private_reference=getarg('private-reference'),
      private_word=getarg('private-word'),r_library=getarg('r-library'))
  }
}
