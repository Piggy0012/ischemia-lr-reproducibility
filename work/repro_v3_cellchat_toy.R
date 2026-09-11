# Tiny smoke of verbatim official source; NOT a production CellChat package run.
root <- normalizePath(getwd(), winslash="/")
.libPaths(c(file.path(root,"work/repro_decontx/rlibrary"),.libPaths()))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))
probe <- file.path(root,"work/repro_v3_cellchat_probe")
args <- commandArgs(trailingOnly=TRUE)
stable <- length(args)>0 && args[1]=="stable"
src <- if(stable) file.path(probe,"stable_v2.1.2") else probe
source(file.path(src,"R/CellChat_class.R"))
source(file.path(src,"R/utilities.R"))
source(file.path(src,"R/modeling.R"))
load(file.path(src,"data/CellChatDB.mouse.rda"))
db <- CellChatDB.mouse
lr <- head(db$interaction[db$interaction$annotation=="Secreted Signaling",,drop=FALSE],6)
subunit_genes <- function(x,table) {
  if (is.na(x) || !nzchar(x)) return(character())
  if (x %in% rownames(table)) return(as.character(table[x,])[nzchar(as.character(table[x,]))])
  x
}
genes <- unique(unlist(c(lapply(c(lr$ligand,lr$receptor),subunit_genes,table=db$complex),
                        lapply(c(lr$agonist,lr$antagonist,lr$co_A_receptor,lr$co_I_receptor),subunit_genes,table=db$cofactor))))
stopifnot(length(genes)>2,length(genes)<100)
set.seed(20260911)
nc <- 120L
x <- matrix(rpois(length(genes)*nc,lambda=2),nrow=length(genes),dimnames=list(genes,paste0("toy",seq_len(nc))))
x[seq_len(ceiling(length(genes)/2)),1:40] <- x[seq_len(ceiling(length(genes)/2)),1:40]+5
x <- log1p(sweep(x,2,pmax(colSums(x),1),"/")*10000)
meta <- data.frame(cell_type=factor(rep(c("Astrocyte","Endothelial","Microglia"),each=40)),samples="toy",row.names=colnames(x))
obj <- createCellChat(x,meta=meta,group.by="cell_type",datatype="RNA")
obj@DB <- db
obj@data.signaling <- obj@data
future::plan("sequential")
t0 <- proc.time()[3]
a <- computeCommunProb(obj,type="triMean",LR.use=lr,raw.use=TRUE,population.size=FALSE,nboot=10,seed.use=20260911,Kh=0.5,n=1)
b <- computeCommunProb(obj,type="triMean",LR.use=lr,raw.use=TRUE,population.size=FALSE,nboot=10,seed.use=20260911,Kh=0.5,n=1)
stopifnot(identical(a@net$prob,b@net$prob),identical(a@net$pval,b@net$pval),all(is.finite(a@net$prob)),length(dim(a@net$prob))==3)
out <- list(status="source_engine_toy_smoke_pass",full_CellChat_package_installed=FALSE,production_qualified=FALSE,
  note="This only tests exact unmodified upstream R core functions on synthetic data. Formal biological runs require the complete native CellChat package.",
  source_commit=if(stable) "a0d3b2d231d46c8787177fffeac908270c253747" else "75253cd0c9e68410e6e721a6d3a0419a1d7e358f",genes=genes,input_shape=dim(x),n_cell_types=3,n_interactions=nrow(lr),
  complex_receptor_rows=sum(lr$receptor %in% rownames(db$complex)),cofactor_rows=sum(nzchar(lr$agonist)|nzchar(lr$antagonist)|nzchar(lr$co_A_receptor)|nzchar(lr$co_I_receptor)),
  score_shape=dim(a@net$prob),score_range=range(a@net$prob),positive_scores=sum(a@net$prob>0),exact_seed_repetition=TRUE,
  elapsed_seconds=unname(proc.time()[3]-t0),nboot=10,all_gene_normalization=FALSE,
  all_gene_normalization_note="Synthetic toy contains only signaling genes; not a biological input normalization strategy.")
suffix <- if(stable) "_v2_1_2" else "_v2_2_dev"
jsonlite::write_json(out,file.path(probe,paste0("source_toy_smoke",suffix,".json")),auto_unbox=TRUE,pretty=TRUE)
saveRDS(list(input=x,meta=meta,LR=lr,prob=a@net$prob,pval=a@net$pval),file.path(probe,paste0("source_toy_smoke",suffix,".rds")))
cat(jsonlite::toJSON(out,auto_unbox=TRUE,pretty=TRUE),"\n")
