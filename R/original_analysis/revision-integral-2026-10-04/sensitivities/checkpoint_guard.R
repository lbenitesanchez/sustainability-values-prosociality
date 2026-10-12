# Checkpoint schema and validation use only base R; no SEM is estimated here.
# Checkpoints, data hashes, resample indices and RDS bundles are private research files.
sens_hash_object <- function(x) {
  f <- tempfile("sens-hash-"); on.exit(unlink(f), add=TRUE)
  con <- file(f,"wb"); writeBin(serialize(x,NULL,version=2L),con); close(con)
  unname(tools::md5sum(f))
}
sens_checkpoint_meta <- function(input,dat,syntax,seed,B,resample_rows,point_names,
                                 lavaan_version=as.character(utils::packageVersion("lavaan"))) {
  if(!file.exists(input)) stop("Missing input SAV: ",input)
  stopifnot(length(resample_rows)==B,all(lengths(resample_rows)==nrow(dat)))
  list(schema_version=1L,statistical_spec_version="pooled-P3-observed-values-v1",
       source_md5=unname(tools::md5sum(input)),source_bytes=unname(file.info(input)$size),
       analytic_data_md5=sens_hash_object(dat),N=as.integer(nrow(dat)),
       syntax_md5=sens_hash_object(syntax),model_id="P3_country_partial_demographics",
       seed=as.integer(seed),bootstrap_seed=as.integer(seed+1L),B=as.integer(B),
       strata_counts=as.integer(table(dat$Pais)),strata_names=names(table(dat$Pais)),
       resample_rows_md5=sens_hash_object(resample_rows),
       resample_lengths=as.integer(lengths(resample_rows)),
       point_names=point_names,standardizations=c("std.all","std.nox"),
       RNGkind=RNGkind(),R_version=as.character(getRversion()),lavaan_version=lavaan_version)
}
sens_validate_results <- function(results,meta) {
  if(!is.list(results)||length(results)>meta$B) stop("Invalid bootstrap result list or B exceeded")
  if(!length(results)) return(invisible(TRUE))
  ids<-vapply(results,function(z) {
    if(!is.list(z)||length(z$i)!=1L||!is.numeric(z$i)||!is.finite(z$i)||z$i!=as.integer(z$i))
      stop("Invalid replicate entry")
    as.integer(z$i)
  },integer(1))
  if(!identical(ids,seq_len(length(results)))) stop("Replicate IDs must be unique consecutive 1..completed")
  for(z in results) {
    if(length(z$valid)!=1L||!is.logical(z$valid)||is.na(z$valid)) stop("Invalid validity flag")
    for(nm in c("values","values_nox")) {
      a<-z[[nm]]
      if(!is.numeric(a)||length(a)!=length(meta$point_names)) stop("Invalid effect-array dimension: ",nm)
      if(isTRUE(z$valid)&&(!identical(names(a),meta$point_names)||!all(is.finite(a))))
        stop("Valid replicate has nonfinite or mismatched effect-array names: ",nm)
    }
  }
  invisible(TRUE)
}
sens_validate_checkpoint <- function(checkpoint,expected_meta) {
  if(!is.list(checkpoint)||!identical(checkpoint$schema_version,1L)||
     !is.list(checkpoint$metadata)||!is.list(checkpoint$results))
    stop("Legacy or malformed checkpoint rejected; explicit verified migration is required")
  actual<-checkpoint$metadata
  if(!identical(names(actual),names(expected_meta))) stop("Checkpoint metadata fields differ")
  mismatch<-names(expected_meta)[!vapply(names(expected_meta),function(k)identical(actual[[k]],expected_meta[[k]]),logical(1))]
  if(length(mismatch)) stop("Stale checkpoint rejected; metadata mismatch: ",paste(mismatch,collapse=", "))
  sens_validate_results(checkpoint$results,expected_meta)
  invisible(TRUE)
}
sens_read_checkpoint <- function(path,expected_meta) {
  if(!file.exists(path)) stop("Requested recovery checkpoint does not exist: ",path)
  z<-readRDS(path);sens_validate_checkpoint(z,expected_meta);z$results
}
sens_write_checkpoint <- function(results,metadata,path,migration=NULL) {
  sens_validate_results(results,metadata)
  z<-list(schema_version=1L,metadata=metadata,results=results,
          written_UTC=format(Sys.time(),tz="UTC",usetz=TRUE),migration=migration)
  tmp<-paste0(path,".tmp-",Sys.getpid())
  on.exit(unlink(tmp),add=TRUE);saveRDS(z,tmp)
  if(!file.rename(tmp,path)) stop("Atomic checkpoint rename failed: ",path)
  invisible(path)
}
