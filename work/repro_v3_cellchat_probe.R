args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if (length(args)) args[1] else getwd(), winslash="/", mustWork=TRUE)
probe <- file.path(root, "work/repro_v3_cellchat_probe")
.libPaths(c(file.path(root, "work/repro_decontx/rlibrary"), .libPaths()))
desc <- read.dcf(file.path(probe, "DESCRIPTION"))
parse_deps <- function(fields) {
  x <- paste(desc[1, intersect(fields, colnames(desc))], collapse=",")
  unique(trimws(sub("\\s*\\(.*\\)", "", strsplit(x, ",")[[1]])))
}
required <- setdiff(parse_deps(c("Depends", "Imports", "LinkingTo")), "R")
ip <- installed.packages()
versions <- setNames(lapply(required, function(p) if (p %in% rownames(ip)) unname(ip[p, "Version"]) else NULL), required)
load(file.path(probe, "data/CellChatDB.mouse.rda"))
db <- CellChatDB.mouse
for (nm in names(db)) {
  if (is.data.frame(db[[nm]])) {
    x <- db[[nm]]
    x <- data.frame(resource_row_id=rownames(x), x, check.names=FALSE)
    write.table(x, file.path(probe, paste0("CellChatDB.mouse_",nm,".tsv")), sep="\t", row.names=FALSE, quote=FALSE, na="")
  }
}
report <- list(
  status="dependency_and_resource_probe_complete",
  R_version=R.version.string,
  R_library_paths=.libPaths(),
  CellChat_installed="CellChat" %in% rownames(ip),
  CellChat_installed_version=if ("CellChat" %in% rownames(ip)) unname(ip["CellChat", "Version"]) else NULL,
  required_versions=versions,
  missing_required=required[!required %in% rownames(ip)],
  resource_components=lapply(db,dim),
  interaction_columns=colnames(db$interaction),
  annotation_counts=as.list(table(db$interaction$annotation)),
  interaction_version_counts=if("version" %in% colnames(db$interaction)) as.list(table(db$interaction$version)) else NULL,
  matrix_dense_conversion_warning="computeCommunProb calls as.matrix(object@data.signaling); estimate dense signaling-gene by context-cell allocation before full runs.",
  no_dependencies_installed=TRUE,
  no_biological_matrices_loaded=TRUE)
jsonlite::write_json(report,file.path(probe,"R_dependency_resource_probe.json"),pretty=TRUE,auto_unbox=TRUE,null="null")
cat(jsonlite::toJSON(report,pretty=TRUE,auto_unbox=TRUE,null="null"),"\n")
