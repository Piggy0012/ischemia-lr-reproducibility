# Explicit user-run restore: installs the 42 archived packages into a NEW library.
# Requires a separately provisioned, version-matching baseline and R 4.6.1.
# First run: python environment/verify_bundle.py <extracted_bundle>
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3L,as.character(getRversion())=="4.6.1",.Platform$OS.type=="windows")
bundle <- normalizePath(args[1],winslash="/",mustWork=TRUE)
baseline <- normalizePath(args[2],winslash="/",mustWork=TRUE)
target <- normalizePath(args[3],winslash="/",mustWork=FALSE)
stopifnot(!identical(tolower(target),tolower(baseline)),!dir.exists(target))
env <- file.path(bundle,"environment")
lock <- read.delim(file.path(env,"package_lock.tsv"),check.names=FALSE,stringsAsFactors=FALSE)
bins <- read.delim(file.path(env,"offline_binary_manifest.tsv"),check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(nrow(bins)==42L)
.libPaths(c(baseline,.Library))
inherited <- lock[lock$scope=="baseline_or_R_runtime",]
before <- installed.packages(lib.loc=.libPaths())
before <- before[!duplicated(before[,"Package"]),,drop=FALSE]
stopifnot(all(inherited$Package %in% before[,"Package"]))
stopifnot(all(before[inherited$Package,"Version"]==inherited$Version))
paths <- file.path(bundle,bins$archive_path)
stopifnot(all(file.exists(paths)),all(unname(tools::md5sum(paths))==bins$md5))
dir.create(target,recursive=TRUE,showWarnings=FALSE)
install.packages(paths,lib=target,repos=NULL,type="win.binary")
.libPaths(c(target,baseline,.Library))
effective <- installed.packages(lib.loc=.libPaths())
effective <- effective[!duplicated(effective[,"Package"]),,drop=FALSE]
expected <- lock[lock$in_analysis_libpaths==1L,]
stopifnot(all(expected$Package %in% effective[,"Package"]),all(effective[expected$Package,"Version"]==expected$Version))
suppressPackageStartupMessages(library(CellChat))
stopifnot(packageDescription("CellChat")$RemoteSha=="75253cd0c9e68410e6e721a6d3a0419a1d7e358f")
stopifnot(normalizePath(find.package("CellChat"),winslash="/")==normalizePath(file.path(target,"CellChat"),winslash="/"))
writeLines(capture.output(sessionInfo()),paste0(target,"_sessionInfo.txt"))
jsonlite::write_json(list(status="offline_overlay_installed_and_namespace_loaded",R_version=R.version.string,
  installed_binary_packages=nrow(bins),effective_analysis_package_versions_checked=nrow(expected),
  CellChat_version=packageDescription("CellChat")$Version,CellChat_RemoteSha=packageDescription("CellChat")$RemoteSha,
  new_library=target,inherited_baseline_library=baseline,baseline_package_versions_checked=nrow(inherited),
  full_baseline_restored=FALSE,network_required=FALSE,model_fit_performed=FALSE),
  paste0(target,"_audit.json"),auto_unbox=TRUE,pretty=TRUE)
cat("42-package offline overlay restored; all effective versions agree; native namespace loaded.\n")
