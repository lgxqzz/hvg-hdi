# Raw counts are supported with the TMM profile.
normalize_input <- function(x, input_scale, profile) {
  validate_matrix(x)
  expected <- profile_settings(profile)$input_scale
  if (identical(input_scale, expected)) return(list(matrix=x, retained=rownames(x), normalization="provided_normalized"))
  if (input_scale != "raw_counts" || profile != "tmm_v1") stop("Input scale and profile incompatible; VST must be supplied with upstream design provenance")
  if (any(x < 0) || any(x != round(x)) || any(colSums(x) <= 0)) stop("Raw counts must be nonnegative integers with positive libraries")
  if (!requireNamespace("edgeR", quietly=TRUE)) stop("Install edgeR for raw-count normalization")
  d <- edgeR::DGEList(counts=x)
  keep <- rowMeans(edgeR::cpm(d) >= 1) >= .20
  if (sum(keep) < 3) stop("Too few genes after CPM filter")
  d <- d[keep,,keep.lib.sizes=FALSE]
  if (any(d$samples$lib.size <= 0)) stop("Empty library after filtering")
  d <- edgeR::calcNormFactors(d,method="TMM")
  list(matrix=edgeR::cpm(d,log=TRUE,prior.count=1), retained=rownames(d),
       normalization="CPM>=1 in >=20%; reset libraries; TMM; logCPM prior.count=1",
       factors=d$samples)
}
