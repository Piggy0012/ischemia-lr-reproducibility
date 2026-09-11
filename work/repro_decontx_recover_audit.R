# Audit-only recovery of a successfully fitted model whose output metadata step was interrupted.
args <- commandArgs(trailingOnly=TRUE)
.libPaths(c(normalizePath(args[1],winslash="/"),.libPaths()))
dest <- normalizePath(args[2],winslash="/")
meta <- jsonlite::read_json(file.path(dest,"input_audit.json"),simplifyVector=TRUE)
saved <- readRDS(file.path(dest,"model_estimates.rds"))
stopifnot(saved$runParams$maxIter==500, saved$runParams$convergence==0.001,
          identical(saved$runParams$delta,c(10,10)), isTRUE(saved$runParams$estimateDelta))
cont <- read.delim(gzfile(file.path(dest,"contamination.tsv.gz")),check.names=FALSE)
cells <- read.delim(gzfile(file.path(dest,"cells.tsv.gz")),check.names=FALSE)
stopifnot(identical(cont$barcode,cells$barcode),identical(saved$runParams$z,cells$whole_brain_broad))
nz <- file.info(file.path(dest,"corrected_data.bin"))$size/8
stopifnot(file.info(file.path(dest,"corrected_indices.bin"))$size/4==nz,
          file.info(file.path(dest,"corrected_indptr.bin"))$size/4==meta$n_cells+1)
con <- file(file.path(dest,"corrected_data.bin"),"rb")
values <- readBin(con,numeric(),n=nz,size=8,endian="little");close(con)
stopifnot(all(is.finite(values)),all(values>=0),all(is.finite(cont$contamination)),
          all(cont$contamination>=0),all(cont$contamination<=1),
          abs(sum(values)-sum(cont$corrected_total))<1e-4,
          all(cont$corrected_total<=cont$raw_total+1e-6))
lines <- readLines(file.path(dest,"run.log"),warn=FALSE)
duration_line <- grep("Completed DecontX. Total time:",lines,value=TRUE)
stopifnot(length(duration_line)==1)
model_minutes <- as.numeric(sub(".*Total time: ([0-9.]+) mins.*","\\1",duration_line))
summary <- list(dataset=meta$dataset,sample=meta$sample,R=as.character(getRversion()),
  decontX=as.character(packageVersion("decontX")),Bioconductor="3.23",
  n_cells=meta$n_cells,n_genes=meta$n_genes,corrected_nnz=nz,
  started_model_log=grep("Analyzing all cells",lines,value=TRUE),
  finished=NULL,elapsed_seconds=NULL,reported_model_elapsed_seconds=model_minutes*60,
  requested_parameters=list(z_source="whole_brain_broad",background=NULL,batch=NULL,
    seed=20260911,delta=c(10,10),estimateDelta=TRUE,maxIter=500,
    convergence=0.001,iterLogLik=10,numerical_threads=1),
  runParams=saved$runParams,estimates_components=lapply(saved$estimates,names),
  batch_iterations=lapply(saved$estimates,function(e)e$iteration),
  estimated_delta=lapply(saved$estimates,function(e)e$delta),
  raw_total=sum(cont$raw_total),corrected_total=sum(values),
  contamination_quantiles=quantile(cont$contamination,c(0,.25,.5,.75,.9,.95,.99,1)),
  n_zero_corrected_total=sum(cont$corrected_total==0),
  fractional_values=any(abs(values-round(values))>1e-8),
  assumptions="Reference broad labels; filtered cells; background=NULL; computational sensitivity, no measured ambient profile",
  audit_recovery=list(reason="Rscript source was edited while a running R interpreter was reading the audit-only tail; the model completed, and all model outputs were saved before the parse error",
    original_process_exit_code=1,model_refitted=FALSE,audit_recovered_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
    raw_outputs_reused=c("corrected_data.bin","corrected_indices.bin","corrected_indptr.bin","corrected_counts.mtx","contamination.tsv.gz","model_estimates.rds","sessionInfo.txt","installed_packages.tsv"),
    timing_note="Original end-to-end R elapsed time was not saved; it is left null rather than fabricated. The model's own reported duration is retained."))
jsonlite::write_json(summary,file.path(dest,"r_audit.json"),auto_unbox=TRUE,pretty=TRUE,digits=16,null="null")
cat("AUDIT_RECOVERED_WITHOUT_MODEL_REFIT",meta$sample,"\n")
