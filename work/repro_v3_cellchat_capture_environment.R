root <- normalizePath(getwd(),winslash="/")
out <- file.path(root,"work/repro_v3_cellchat")
.libPaths(c(file.path(out,"rlibrary"),file.path(root,"work/repro_decontx/rlibrary"),.libPaths()))
suppressPackageStartupMessages(library(CellChat))
ip <- installed.packages()
ip <- ip[!duplicated(ip[,"Package"]),,drop=FALSE]
cols <- intersect(c("Package","Version","LibPath","Built","Priority"),colnames(ip))
versions <- as.data.frame(ip[,cols,drop=FALSE],stringsAsFactors=FALSE)
versions$RemoteSha <- vapply(versions$Package,function(p) {
  x <- packageDescription(p)$RemoteSha;if(is.null(x))"" else x
},character(1))
write.table(versions,file.path(out,"R_packages_effective.tsv"),sep="\t",row.names=FALSE,quote=FALSE,na="")
writeLines(capture.output(sessionInfo()),file.path(out,"analysis_sessionInfo.txt"))
jsonlite::write_json(list(R_version=R.version.string,platform=R.version$platform,
  runtime_path=R.home(),libraries=.libPaths(),effective_package_count=nrow(versions),
  CellChat_version=as.character(packageVersion("CellChat")),CellChat_RemoteSha=packageDescription("CellChat")$RemoteSha,
  source="full installed native package; no LIANA CellChat wrapper; existing non-v3 libraries are read-only"),
  file.path(out,"environment_summary.json"),auto_unbox=TRUE,pretty=TRUE)
cat("ENVIRONMENT_CAPTURED",nrow(versions),"packages\n")
