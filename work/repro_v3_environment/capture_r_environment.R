args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3L)
analysis_lib <- normalizePath(args[1], winslash="/", mustWork=TRUE)
tooling_lib <- normalizePath(args[2], winslash="/", mustWork=TRUE)
dest <- normalizePath(args[3], winslash="/", mustWork=TRUE)
.libPaths(c(analysis_lib, tooling_lib, .Library))
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(decontX))
suppressPackageStartupMessages(library(renv))
stopifnot(as.character(getRversion())=="4.6.1", as.character(packageVersion("decontX"))=="1.10.0")
options(repos=BiocManager::repositories(version="3.23"), renv.bioconductor.version="3.23")
ip <- installed.packages(lib.loc=.libPaths(), fields=c("Repository","RemoteType","RemoteUrl","RemoteSha","biocViews"))
ip <- ip[!duplicated(ip[,"Package"]),,drop=FALSE]
write.table(ip, file.path(dest,"R_installed_packages.tsv"), sep="\t", row.names=FALSE, quote=FALSE, na="")
descriptions <- list()
dir.create(file.path(dest,"R_DESCRIPTION"),showWarnings=FALSE)
for (i in seq_len(nrow(ip))) {
  name <- ip[i,"Package"]
  path <- file.path(ip[i,"LibPath"],name,"DESCRIPTION")
  file.copy(path,file.path(dest,"R_DESCRIPTION",paste0(name,".DESCRIPTION")),overwrite=TRUE)
  fields <- as.list(read.dcf(path)[1,])
  descriptions[[name]] <- fields
}
write_json(descriptions,file.path(dest,"R_package_descriptions.json"),pretty=TRUE,auto_unbox=TRUE)
lock <- renv::lockfile_create(type="all",libpaths=.libPaths(),project=dest,prompt=FALSE)
renv::lockfile_write(lock,file=file.path(dest,"renv.lock"),project=dest)
check <- renv::lockfile_read(file.path(dest,"renv.lock"),project=dest)
package_checks <- vapply(check$Packages,function(x) identical(packageDescription(x$Package,fields="Version"),x$Version),logical(1))
stopifnot(all(package_checks))
writeLines(capture.output(sessionInfo()),file.path(dest,"R_sessionInfo.txt"))
write_json(list(R_version=as.character(getRversion()),platform=R.version$platform,
  Bioconductor_version=as.character(BiocManager::version()),renv_version=as.character(packageVersion("renv")),
  decontX_version=as.character(packageVersion("decontX")),locked_packages=length(check$Packages),
  all_locked_versions_match_installed=all(package_checks),lockfile_generated_by="renv::lockfile_create(type='all')",
  original_analysis_library=analysis_lib,isolated_lock_tooling_library=tooling_lib,
  installed_package_count=nrow(ip),original_analysis_library_modified=FALSE,
  full_fresh_restore_tested=FALSE,model_fit_performed=FALSE),
  file.path(dest,"R_lock_capture_audit.json"),pretty=TRUE,auto_unbox=TRUE)
cat("Locked",length(check$Packages),"packages; all installed versions agree. No restore or model fit performed.\n")
