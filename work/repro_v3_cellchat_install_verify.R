root <- normalizePath(getwd(),winslash="/")
out <- file.path(root,"work/repro_v3_cellchat")
lib <- file.path(out,"rlibrary")
dir.create(lib,recursive=TRUE,showWarnings=FALSE)
.libPaths(c(lib,file.path(root,"work/repro_decontx/rlibrary"),.libPaths()))
manifest <- jsonlite::read_json(file.path(out,"download_manifest.json"),simplifyVector=TRUE)
if (!file.exists(file.path(out,"installed_package_verified.json"))) {
  zips <- file.path(root,manifest$packages$local_path)
  install.packages(zips,repos=NULL,type="win.binary",lib=lib)
}
suppressPackageStartupMessages(library(CellChat))
stopifnot(as.character(packageVersion("CellChat"))=="2.2.0.9001")
desc <- packageDescription("CellChat")
stopifnot(desc$RemoteSha==manifest$remote_sha)
official <- new.env(parent=globalenv())
probe <- file.path(root,"work/repro_v3_cellchat_probe")
sys.source(file.path(probe,"R/modeling.R"),envir=official)
sys.source(file.path(probe,"R/utilities.R"),envir=official)
ns <- asNamespace("CellChat")
fn <- Filter(function(x) is.function(get(x,official)) && exists(x,ns,inherits=FALSE),ls(official,all.names=TRUE))
checks <- lapply(fn,function(x) {
  a <- get(x,official);b <- get(x,ns)
  list(function_name=x,body_matches=identical(deparse(body(a)),deparse(body(b))),formals_match=identical(formals(a),formals(b)))
})
stopifnot(all(vapply(checks,function(x) x$body_matches && x$formals_match,logical(1))))
expected <- new.env(parent=baseenv());actual <- new.env(parent=baseenv())
load(file.path(probe,"data/CellChatDB.mouse.rda"),envir=expected)
data(list="CellChatDB.mouse",package="CellChat",envir=actual)
stopifnot(identical(expected$CellChatDB.mouse,actual$CellChatDB.mouse))
db <- actual$CellChatDB.mouse
lr <- head(db$interaction[db$interaction$annotation=="Secreted Signaling",,drop=FALSE],6)
parts <- function(x,tab) {
  if(is.na(x)||!nzchar(x))return(character())
  if(x %in% rownames(tab))return(as.character(tab[x,])[nzchar(as.character(tab[x,]))])
  x
}
genes <- unique(unlist(c(lapply(c(lr$ligand,lr$receptor),parts,tab=db$complex),lapply(c(lr$agonist,lr$antagonist,lr$co_A_receptor,lr$co_I_receptor),parts,tab=db$cofactor))))
set.seed(20260911)
x <- matrix(rpois(length(genes)*120,2),nrow=length(genes),dimnames=list(genes,paste0("toy",seq_len(120))))
x[seq_len(ceiling(length(genes)/2)),1:40] <- x[seq_len(ceiling(length(genes)/2)),1:40]+5
x <- log1p(sweep(x,2,pmax(colSums(x),1),"/")*10000)
meta <- data.frame(cell_type=factor(rep(c("Astrocyte","Endothelial","Microglia"),each=40)),samples=factor(rep("toy",120)),row.names=colnames(x))
obj <- CellChat::createCellChat(x,meta=meta,group.by="cell_type",datatype="RNA")
obj@DB <- db; obj@data.signaling <- obj@data
future::plan("sequential")
t0 <- proc.time()[3]
a <- CellChat::computeCommunProb(obj,type="triMean",LR.use=lr,raw.use=TRUE,population.size=FALSE,nboot=10,seed.use=20260911,Kh=0.5,n=1)
b <- CellChat::computeCommunProb(obj,type="triMean",LR.use=lr,raw.use=TRUE,population.size=FALSE,nboot=10,seed.use=20260911,Kh=0.5,n=1)
stopifnot(identical(a@net$prob,b@net$prob),identical(a@net$pval,b@net$pval),all(is.finite(a@net$prob)),all(is.finite(a@net$pval)))
ip <- installed.packages()
report <- list(status="complete_native_package_verified",R_version=R.version.string,CellChat_version=as.character(packageVersion("CellChat")),
  RemoteSha=desc$RemoteSha,RemoteUrl=desc$RemoteUrl,library=lib,loaded_from=find.package("CellChat"),
  source_functions_verified=checks,source_function_count=length(checks),entire_native_mouse_DB_identical=TRUE,
  native_toy=list(status="pass",input_shape=dim(x),LR_rows=nrow(lr),output_shape=dim(a@net$prob),positive_scores=sum(a@net$prob>0),prob_range=range(a@net$prob),seed_repeat_exact=TRUE,elapsed_seconds=unname(proc.time()[3]-t0)),
  newly_installed_packages=unname(as.data.frame(ip[ip[,"LibPath"]==lib,c("Package","Version","LibPath")],stringsAsFactors=FALSE)),
  prior_library_modified=FALSE,scientific_run_uses_installed_namespace=TRUE)
jsonlite::write_json(report,file.path(out,"installed_package_verified.json"),auto_unbox=TRUE,pretty=TRUE,digits=17)
writeLines(capture.output(sessionInfo()),file.path(out,"installation_sessionInfo.txt"))
saveRDS(list(input=x,meta=meta,LR=lr,prob=a@net$prob,pval=a@net$pval),file.path(out,"native_package_toy.rds"))
cat("NATIVE_PACKAGE_VERIFIED",length(checks),"functions, complete DB, toy seconds",report$native_toy$elapsed_seconds,"\n")
