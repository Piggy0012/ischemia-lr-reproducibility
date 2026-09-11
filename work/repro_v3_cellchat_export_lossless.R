# Preserve native double precision in decimal TSV; never recompute a model.
root <- normalizePath(getwd(),winslash="/")
.libPaths(c(file.path(root,"work/repro_decontx/rlibrary"),.libPaths()))
base <- file.path(root,"work/repro_v3_cellchat/controlled")
folders <- list.dirs(base,recursive=TRUE,full.names=TRUE)
folders <- folders[file.exists(file.path(folders,"run_telemetry.json"))]
sha <- function(x) digest::digest(file=x,algo="sha256",serialize=FALSE)
for (folder in folders) {
  telemetry <- jsonlite::read_json(file.path(folder,"run_telemetry.json"),simplifyVector=TRUE)
  if(telemetry$status!="complete")next
  rds <- file.path(folder,"native_network.rds")
  stopifnot(sha(rds)==telemetry$output_sha256[["native_network.rds"]])
  target <- file.path(folder,"full_network_lossless.tsv.gz")
  auditfile <- file.path(folder,"lossless_export_audit.json")
  if(file.exists(target)&&file.exists(auditfile)) {
    old <- jsonlite::read_json(auditfile,simplifyVector=TRUE)
    stopifnot(old$source_rds_sha256==sha(rds),old$output_tsv_sha256==sha(target))
    next
  }
  a <- readRDS(rds);dimn <- dimnames(a$net$prob)
  p <- as.vector(a$net$prob);pv <- as.vector(a$net$pval)
  grid <- expand.grid(source=dimn[[1]],target=dimn[[2]],interaction_name=dimn[[3]],stringsAsFactors=FALSE)
  grid$probability <- sprintf("%.17g",p);grid$pvalue <- sprintf("%.17g",pv)
  grid$dataset <- a$input_audit$dataset;grid$sample <- a$input_audit$sample;grid$condition <- a$input_audit$condition
  con <- gzfile(target,"wt");write.table(grid,con,sep="\t",quote=FALSE,row.names=FALSE);close(con)
  back <- read.delim(gzfile(target),check.names=FALSE,stringsAsFactors=FALSE)
  stopifnot(identical(as.double(back$probability),p),identical(as.double(back$pvalue),pv))
  original <- read.delim(gzfile(file.path(folder,"full_network.tsv.gz")),check.names=FALSE)
  report <- list(status="lossless_export_verified",source_rds_sha256=sha(rds),output_tsv_sha256=sha(target),
      original_tsv_sha256=sha(file.path(folder,"full_network.tsv.gz")),n_rows=length(p),
      n_original_tsv_probability_values_changed=sum(original$probability!=p),
      max_original_tsv_probability_absolute_error=max(abs(original$probability-p)),
      exact_probability_roundtrip=TRUE,exact_pvalue_roundtrip=TRUE,
      method="R native double -> sprintf('%.17g') character -> TSV -> R double, with identical assertions",
      no_model_rerun=TRUE,no_row_filter=TRUE,
      source_script_sha256=sha(file.path(root,"work/repro_v3_cellchat_export_lossless.R")))
  jsonlite::write_json(report,auditfile,auto_unbox=TRUE,pretty=TRUE,digits=17)
  cat("LOSSLESS_EXPORTED",a$input_audit$sample,length(p),"original rounded values",report$n_original_tsv_probability_values_changed,"\n")
  rm(a,grid,back,original,p,pv);gc()
}
