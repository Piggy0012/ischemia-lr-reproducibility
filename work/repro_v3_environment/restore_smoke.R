# Bounded restoration smoke test: two infrastructure packages, no scientific fits.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) %in% c(2L,3L))
dest <- normalizePath(args[1],winslash="/",mustWork=TRUE)
tooling <- normalizePath(args[2],winslash="/",mustWork=TRUE)
.libPaths(c(tooling,.Library))
library(renv)
target <- file.path(dest,if (length(args)==3L) args[3] else "restore_smoke_library")
dir.create(target,showWarnings=FALSE)
stopifnot(length(list.files(target))==0L)
renv::restore(project=dest,lockfile=file.path(dest,"renv.lock"),library=target,
              packages=c("renv","jsonlite"),prompt=FALSE)
.libPaths(c(target,tooling,.Library))
stopifnot(packageDescription("renv",lib.loc=target,fields="Version")=="1.2.4",
          packageDescription("jsonlite",lib.loc=target,fields="Version")=="2.0.0")
jsonlite::write_json(list(status="passed",packages=list(renv="1.2.4",jsonlite="2.0.0"),
  restored_library_packages=as.data.frame(installed.packages(lib.loc=target)[,c("Package","Version"),drop=FALSE]),
  restored_to_empty_library=TRUE,full_232_package_restore_tested=FALSE,
  decontX_restore_tested=FALSE,original_analysis_library_modified=FALSE,
  model_fit_performed=FALSE),file.path(dest,"R_restore_smoke_audit.json"),pretty=TRUE,auto_unbox=TRUE)
cat("Two-package empty-library restoration passed; full scientific environment restore remains untested.\n")
