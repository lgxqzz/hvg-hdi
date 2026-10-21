# Shared strict contracts. Source before computational modules.
validate_ids <- function(x, label) {
  if (!is.character(x) || is.null(x) || anyNA(x) || any(!nzchar(trimws(x))) || anyDuplicated(x) ||
      any(grepl("[\t\r\n]", x))) stop(label, " must be unique, nonempty IDs without control characters")
  invisible(x)
}
validate_matrix <- function(x) {
  if (!is.matrix(x) || !is.numeric(x) || any(!is.finite(x)) ||
      nrow(x) < 3L || ncol(x) < 2L) stop("Input must be a finite numeric genes-by-samples matrix (>=3 genes, >=2 samples)")
  validate_ids(rownames(x), "Gene IDs"); validate_ids(colnames(x), "Sample IDs")
  invisible(x)
}
align_metadata <- function(x, meta, id_col, cycle_col) {
  validate_matrix(x)
  if (!is.data.frame(meta) || !all(c(id_col, cycle_col) %in% names(meta))) stop("Missing metadata columns")
  if (anyDuplicated(names(meta))) stop("Duplicate metadata column names")
  validate_ids(as.character(meta[[id_col]]), "Metadata sample IDs")
  if (!setequal(colnames(x), meta[[id_col]])) stop("Matrix and metadata sample sets must match exactly")
  m <- meta[match(colnames(x), meta[[id_col]]), , drop=FALSE]
  if (anyNA(m[[cycle_col]]) || any(!nzchar(trimws(m[[cycle_col]])))) stop("Missing cycle labels")
  m
}
check_probability <- function(x, label) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x < 0 || x > 1) stop("Invalid probability: ", label)
}
write_tsv <- function(x, file) utils::write.table(x, file, sep="\t", row.names=FALSE, quote=TRUE, na="NA")
