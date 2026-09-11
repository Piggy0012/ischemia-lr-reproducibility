args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
root <- normalizePath(getwd(),winslash="/")
out <- file.path(root,"work/repro_v3_cellchat")
.libPaths(c(file.path(out,"rlibrary"),file.path(root,"work/repro_decontx/rlibrary"),.libPaths()))
suppressPackageStartupMessages(library(CellChat))
future::plan("sequential")
input <- normalizePath(args[1],winslash="/",mustWork=TRUE)
meta <- read.delim(gzfile(file.path(input,"cells.tsv.gz")),check.names=FALSE,stringsAsFactors=FALSE)
genes <- read.delim(file.path(input,"genes.tsv"),check.names=FALSE,stringsAsFactors=FALSE)$symbol
x <- as(Matrix::readMM(gzfile(file.path(input,"logcp10k_signal.mtx.gz"))),"CsparseMatrix")
stopifnot(nrow(x)==length(genes),ncol(x)==nrow(meta),all(is.finite(x@x)),all(x@x>0))
rownames(x) <- genes;colnames(x) <- meta$barcode
rownames(meta) <- meta$barcode;meta$cell_type <- factor(meta$cell_type);meta$samples <- factor(meta$sample)
audit <- jsonlite::read_json(file.path(input,"input_audit.json"),simplifyVector=TRUE)
resource <- read.delim(file.path(out,"controlled_resource.tsv"),check.names=FALSE,stringsAsFactors=FALSE)
db <- CellChatDB.mouse
db$interaction <- db$interaction[match(resource$resource_row_id,rownames(db$interaction)),,drop=FALSE]
stopifnot(nrow(db$interaction)==1548,!anyNA(rownames(db$interaction)),setequal(CellChat::extractGene(db),genes))
obj <- CellChat::createCellChat(x,meta=meta,group.by="cell_type",datatype="RNA")
obj@DB <- db
obj <- CellChat::subsetData(obj)
stopifnot(nrow(obj@data.signaling)==length(genes),ncol(obj@data.signaling)==nrow(meta))
t0 <- proc.time()[3]
obj <- CellChat::computeCommunProb(obj,type="triMean",LR.use=db$interaction,raw.use=TRUE,population.size=FALSE,
                                 nboot=100,seed.use=20260911,Kh=0.5,n=1)
elapsed <- unname(proc.time()[3]-t0)
stopifnot(all(is.finite(obj@net$prob)),all(is.finite(obj@net$pval)),all(obj@net$prob>=0),all(obj@net$pval>=0 & obj@net$pval<=1))
saveRDS(list(net=obj@net,parameters=obj@options$parameter,LR=db$interaction,meta=meta,
             input_audit=audit,cellchat_version=as.character(packageVersion("CellChat"))),
        file.path(input,"native_network.rds"),compress="gzip")
d <- dimnames(obj@net$prob)
grid <- expand.grid(source=d[[1]],target=d[[2]],interaction_name=d[[3]],stringsAsFactors=FALSE)
grid$probability <- as.vector(obj@net$prob);grid$pvalue <- as.vector(obj@net$pval)
grid$dataset <- audit$dataset;grid$sample <- audit$sample;grid$condition <- audit$condition
stopifnot(nrow(grid)==prod(dim(obj@net$prob)),!anyDuplicated(grid[c("source","target","interaction_name")]))
con <- gzfile(file.path(input,"full_network.tsv.gz"),"wt")
oldopt <- options(digits=17)
write.table(grid,con,sep="\t",row.names=FALSE,quote=FALSE,na="NA")
close(con);options(oldopt)
report <- list(status="complete",arm="controlled_shared_resource",dataset=audit$dataset,sample=audit$sample,condition=audit$condition,
  parameters=obj@options$parameter,native_CellChat_version=as.character(packageVersion("CellChat")),
  native_RemoteSha=packageDescription("CellChat")$RemoteSha,actual_namespace_path=find.package("CellChat"),
  n_cells=nrow(meta),n_signal_genes=length(genes),n_context_types=nlevels(meta$cell_type),n_LR=nrow(db$interaction),
  network_rows=nrow(grid),zero_probability_rows=sum(grid$probability==0),
  target_rows=sum((grid$source=="Astrocyte" & grid$target=="Endothelial") | (grid$source=="Endothelial" & grid$target=="Astrocyte")),
  probability_min=min(grid$probability),probability_max=max(grid$probability),
  elapsed_compute_seconds=elapsed,all_probability_rows_preserved=TRUE,pvalue_filter_used=FALSE,
  inferential_unit_note="Native cell-label permutation P values are not animal-level disease P values.",
  rank_note="No native LR consensus rank is produced; any later strength percentile is analyst-derived.")
jsonlite::write_json(report,file.path(input,"model_audit.json"),auto_unbox=TRUE,pretty=TRUE,digits=17)
cat("CELLCHAT_COMPLETE",audit$sample,"seconds",elapsed,"rows",nrow(grid),"\n")
