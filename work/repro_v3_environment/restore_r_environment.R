# Explicit opt-in full restore into this environment directory; never run by capture/verification.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
project <- normalizePath(args[1],winslash="/",mustWork=TRUE)
stopifnot(as.character(getRversion())=="4.6.1")
bootstrap <- file.path(project,"renv_bootstrap")
target <- file.path(project,"r-library")
dir.create(bootstrap,showWarnings=FALSE)
dir.create(target,showWarnings=FALSE)
.libPaths(c(bootstrap,.Library))
options(repos=c(CRAN="https://cloud.r-project.org"))
if (!requireNamespace("renv",quietly=TRUE)) {
  url <- "https://cran.r-project.org/src/contrib/Archive/renv/renv_1.2.4.tar.gz"
  # The pinned current version uses its current CRAN source endpoint.
  current <- available.packages()["renv","Version"]
  if (identical(current,"1.2.4")) url <- "https://cran.r-project.org/src/contrib/renv_1.2.4.tar.gz"
  install.packages(url,lib=bootstrap,repos=NULL,type="source")
}
stopifnot(as.character(packageVersion("renv"))=="1.2.4")
renv::restore(project=project,lockfile=file.path(project,"renv.lock"),library=target,prompt=FALSE)
.libPaths(c(target,bootstrap,.Library))
locked <- renv::lockfile_read(file.path(project,"renv.lock"))$Packages
stopifnot(all(vapply(locked,function(x) identical(packageDescription(x$Package,fields="Version"),x$Version),logical(1))))
cat("Full lock restore and version verification completed in",target,"\n")
