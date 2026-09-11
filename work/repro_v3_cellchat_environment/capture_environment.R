# Capture only: no package installation, restore, update or model fitting.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
root <- normalizePath(args[1],winslash="/",mustWork=TRUE)
dest <- file.path(root,"work/repro_v3_cellchat_environment")
overlay <- file.path(root,"work/repro_v3_cellchat/rlibrary")
baseline <- file.path(root,"work/repro_decontx/rlibrary")
tooling <- file.path(root,"work/repro_v3_environment/renv_tooling")
stopifnot(dir.exists(overlay),dir.exists(baseline),dir.exists(tooling))
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
.libPaths(c(overlay,baseline,.Library))
suppressPackageStartupMessages(library(CellChat))
stopifnot(as.character(getRversion())=="4.6.1",
          as.character(packageVersion("CellChat"))=="2.2.0.9001",
          packageDescription("CellChat")$RemoteSha=="75253cd0c9e68410e6e721a6d3a0419a1d7e358f")
effective <- installed.packages(lib.loc=.libPaths(),fields=c("Repository","RemoteUrl","RemoteSha","biocViews"))
effective <- effective[!duplicated(effective[,"Package"]),,drop=FALSE]
effective <- effective[order(effective[,"Package"]),,drop=FALSE]
write.table(effective,file.path(dest,"R_packages_analysis_effective.tsv"),sep="\t",row.names=FALSE,quote=TRUE,qmethod="double",na="")
writeLines(capture.output(sessionInfo()),file.path(dest,"R_analysis_sessionInfo.txt"))
analysis_libpaths <- .libPaths()
loaded <- loadedNamespaces()
loaded_df <- data.frame(Package=loaded,Version=vapply(loaded,function(p)as.character(packageVersion(p)),character(1)),
                       Library=vapply(loaded,function(p)dirname(find.package(p)),character(1)))
write.table(loaded_df[order(loaded_df$Package),],file.path(dest,"R_loaded_namespaces.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
# The capture-only renv library is last and cannot shadow an analysis package.
.libPaths(c(analysis_libpaths,tooling))
suppressPackageStartupMessages(library(renv))
stopifnot(as.character(packageVersion("renv"))=="1.2.4")
ip <- installed.packages(lib.loc=.libPaths(),fields=c("Repository","RemoteUrl","RemoteSha","biocViews"))
ip <- ip[!duplicated(ip[,"Package"]),,drop=FALSE]
ip <- ip[order(ip[,"Package"]),,drop=FALSE]
stopifnot(setequal(setdiff(ip[,"Package"],effective[,"Package"]),"renv"))
write.table(ip,file.path(dest,"R_packages_lock_effective.tsv"),sep="\t",row.names=FALSE,quote=TRUE,qmethod="double",na="")
dir.create(file.path(dest,"R_DESCRIPTION"),showWarnings=FALSE)
descriptions <- list()
for (i in seq_len(nrow(ip))) {
  pkg <- ip[i,"Package"]
  desc <- file.path(ip[i,"LibPath"],pkg,"DESCRIPTION")
  stopifnot(file.copy(desc,file.path(dest,"R_DESCRIPTION",paste0(pkg,".DESCRIPTION")),overwrite=TRUE))
  descriptions[[pkg]] <- as.list(read.dcf(desc)[1,])
}
jsonlite::write_json(descriptions,file.path(dest,"R_package_descriptions.json"),pretty=TRUE,auto_unbox=TRUE)
# Same Bioconductor release as the inherited decontX environment. Repository and
# RemoteSha fields remain those inferred by renv from the installed packages.
options(repos=c(BioCsoft="https://bioconductor.org/packages/3.23/bioc",
                BioCann="https://bioconductor.org/packages/3.23/data/annotation",
                BioCexp="https://bioconductor.org/packages/3.23/data/experiment",
                BioCworkflows="https://bioconductor.org/packages/3.23/workflows",
                BioCbooks="https://bioconductor.org/packages/3.23/books",
                CRAN="https://cloud.r-project.org"),renv.bioconductor.version="3.23")
lock <- renv::lockfile_create(type="all",libpaths=.libPaths(),project=dest,prompt=FALSE)
renv::lockfile_write(lock,file=file.path(dest,"renv.lock"),project=dest)
readback <- renv::lockfile_read(file.path(dest,"renv.lock"),project=dest)
checks <- vapply(readback$Packages,function(x)identical(packageDescription(x$Package,fields="Version"),x$Version),logical(1))
stopifnot(all(checks))
baseline_ip <- installed.packages(lib.loc=c(baseline,.Library))
baseline_ip <- baseline_ip[!duplicated(baseline_ip[,"Package"]),,drop=FALSE]
new <- setdiff(effective[,"Package"],baseline_ip[,"Package"])
shared <- intersect(effective[,"Package"],baseline_ip[,"Package"])
changed <- shared[effective[shared,"Version"]!=baseline_ip[shared,"Version"]]
report <- list(status="actual_mixed_library_lock_captured",R_version=R.version.string,
  platform=R.version$platform,R_home=R.home(),analysis_libpaths=analysis_libpaths,
  capture_only_tooling_library=tooling,CellChat_version=as.character(packageVersion("CellChat")),
  CellChat_RemoteSha=packageDescription("CellChat")$RemoteSha,decontX_version=as.character(packageVersion("decontX")),
  analysis_effective_packages=nrow(effective),baseline_effective_packages=nrow(baseline_ip),
  overlay_new_packages=sort(new),overlay_existing_packages_with_changed_versions=sort(changed),
  installed_packages_in_lock_capture=nrow(ip),locked_nonbase_packages=length(readback$Packages),
  all_lock_versions_match_installed=all(checks),lock_generator="renv::lockfile_create(type='all')",
  capture_only_package="renv 1.2.4",full_fresh_restore_tested=FALSE,
  baseline_library_modified=FALSE,overlay_library_modified=FALSE,model_fit_performed=FALSE)
jsonlite::write_json(report,file.path(dest,"capture_audit.json"),pretty=TRUE,auto_unbox=TRUE)
cat("Captured",nrow(effective),"analysis packages;",length(readback$Packages),"nonbase lock entries;",
    length(new),"new and",length(changed),"changed baseline versions.\n")
