Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
           VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1", RCPP_PARALLEL_NUM_THREADS="1")
args <- commandArgs(trailingOnly=TRUE)
.libPaths(c(normalizePath(args[1],winslash="/"), .libPaths()))
dest <- normalizePath(args[2],winslash="/")
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(decontX))
RcppParallel::setThreadOptions(numThreads=1)
meta <- jsonlite::read_json(file.path(dest,"input_audit.json"),simplifyVector=TRUE)
read_binary <- function(name, what, n, size) {
  con <- file(file.path(dest,name),"rb"); on.exit(close(con))
  readBin(con,what=what,n=n,size=size,endian="little")
}
nr <- meta$n_genes; nc <- meta$n_cells; nz <- meta$nnz
x <- new("dgCMatrix", i=read_binary("input_indices.bin",integer(),nz,4),
         p=read_binary("input_indptr.bin",integer(),nc+1L,4),
         x=as.numeric(read_binary("input_data.bin",integer(),nz,4)),
         Dim=as.integer(c(nr,nc)))
genes <- read.delim(gzfile(file.path(dest,"genes.tsv.gz")),check.names=FALSE)$symbol
cells <- read.delim(gzfile(file.path(dest,"cells.tsv.gz")),check.names=FALSE)
rownames(x) <- genes; colnames(x) <- cells$barcode
stopifnot(all(colSums(x)==cells$n_umis), length(cells$whole_brain_broad)==nc,
          !anyNA(cells$whole_brain_broad), anyDuplicated(cells$barcode)==0,
          anyDuplicated(genes)==0, all(x@x>=0), all(x@x==round(x@x)))
set.seed(20260911)
started <- Sys.time()
cat("START",meta$dataset,meta$sample,"genes",nr,"cells",nc,"nnz",nz,"\n")
fit <- decontX(x, z=cells$whole_brain_broad, background=NULL, batch=NULL,
              seed=20260911, delta=c(10,10), estimateDelta=TRUE,
              maxIter=500, convergence=0.001, iterLogLik=10,
              verbose=TRUE)
corrected <- as(fit$decontXcounts,"dgCMatrix")
stopifnot(identical(dim(corrected),dim(x)), identical(dimnames(corrected),dimnames(x)), all(is.finite(corrected@x)),
          all(corrected@x>=0), length(fit$contamination)==nc,
          all(is.finite(fit$contamination)), all(fit$contamination>=0),
          all(fit$contamination<=1))
corrected_total <- colSums(corrected)
stopifnot(all(corrected_total<=cells$n_umis+1e-6))
write_binary <- function(value,name,size) {
  con <- file(file.path(dest,name),"wb"); on.exit(close(con))
  writeBin(value,con,size=size,endian="little")
}
write_binary(corrected@x,"corrected_data.bin",8)
write_binary(corrected@i,"corrected_indices.bin",4)
write_binary(corrected@p,"corrected_indptr.bin",4)
writeMM(corrected,file.path(dest,"corrected_counts.mtx"))
contamination <- data.frame(barcode=cells$barcode, whole_brain_broad=cells$whole_brain_broad,
                           contamination=fit$contamination, raw_total=cells$n_umis,
                           corrected_total=corrected_total)
con <- gzfile(file.path(dest,"contamination.tsv.gz"),"wt")
write.table(contamination,con,sep="\t",row.names=FALSE,quote=FALSE); close(con)
saveRDS(list(runParams=fit$runParams, estimates=fit$estimates, z=fit$z),
        file.path(dest,"model_estimates.rds"),compress="gzip")
capture.output(sessionInfo(),file=file.path(dest,"sessionInfo.txt"))
pkg <- as.data.frame(installed.packages()[,c("Package","Version","LibPath")])
write.table(pkg,file.path(dest,"installed_packages.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
summary <- list(dataset=meta$dataset,sample=meta$sample,
                R=as.character(getRversion()),decontX=as.character(packageVersion("decontX")),
                Bioconductor=as.character(BiocManager::version()),
                n_cells=nc,n_genes=nr,corrected_nnz=length(corrected@x),
                started=as.character(started),finished=as.character(Sys.time()),
                elapsed_seconds=as.numeric(difftime(Sys.time(),started,units="secs")),
                requested_parameters=list(z_source="whole_brain_broad",background=NULL,batch=NULL,
                  seed=20260911,delta=c(10,10),estimateDelta=TRUE,maxIter=500,
                  convergence=0.001,iterLogLik=10,numerical_threads=1),
                runParams=fit$runParams,estimates_components=lapply(fit$estimates,names),
                batch_iterations=lapply(fit$estimates,function(e)e$iteration),
                estimated_delta=lapply(fit$estimates,function(e)e$delta),
                raw_total=sum(cells$n_umis),corrected_total=sum(corrected_total),
                contamination_quantiles=quantile(fit$contamination,c(0,.25,.5,.75,.9,.95,.99,1)),
                n_zero_corrected_total=sum(corrected_total==0),
                fractional_values=any(abs(corrected@x-round(corrected@x))>1e-8),
                assumptions="Reference broad labels; filtered cells; background=NULL; computational sensitivity, no measured ambient profile")
jsonlite::write_json(summary,file.path(dest,"r_audit.json"),auto_unbox=TRUE,pretty=TRUE,digits=16,null="null")
cat("COMPLETE",meta$sample,"median_contamination",median(fit$contamination),"\n")
