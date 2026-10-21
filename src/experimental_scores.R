# Explicit experimental Hydroxymethylation Dispersion Index definitions.
# Ordered component references with explicit scoring definitions.
freeze_dispersion_reference <- function(raw, sample_ids) {
  raw <- as.matrix(raw)
  if (!is.numeric(raw) || nrow(raw)<2L || any(!is.finite(raw))) stop("Finite numeric reference with >=2 rows required")
  if (!identical(colnames(raw),c("IQR","MAD","ENT","Gini")) &&
      !identical(colnames(raw),c("IQR","MAD","ENT","Gini","BIM"))) stop("Ordered four or five components required")
  validate_ids(sample_ids,"Reference sample IDs")
  if (length(sample_ids)!=nrow(raw) || !identical(rownames(raw),sample_ids)) stop("Reference alignment mismatch")
  structure(list(center=colMeans(raw),sample_sd=apply(raw,2,sd),sample_ids=sample_ids,n=nrow(raw)),class="dispersion_reference")
}
score_experimental_dispersion <- function(raw, sample_ids, reference, definition, accepted=NULL) {
  if (length(definition)!=1L || is.na(definition) || !definition %in% c("HDI_4","HDI_5")) stop("Choose HDI_4 or HDI_5")
  components <- c("IQR","MAD","ENT","Gini",if(definition=="HDI_5") "BIM")
  raw <- as.matrix(raw)
  if (!inherits(reference,"dispersion_reference") || !identical(colnames(raw),components) ||
      !identical(names(reference$center),components) || !identical(names(reference$sample_sd),components)) stop("Definition/reference mismatch")
  validate_ids(sample_ids,"Scoring sample IDs")
  if (length(sample_ids)!=nrow(raw) || !identical(rownames(raw),sample_ids)) stop("Scoring alignment mismatch")
  if (any(!is.finite(raw)) || any(!is.finite(reference$sample_sd)|reference$sample_sd<=0)) stop("Undefined/constant component")
  if (definition=="HDI_5" && !identical(sample_ids,reference$sample_ids)) stop("HDI_5 requires the entire frozen scoring cohort")
  if (definition=="HDI_5" && (is.null(accepted) || !is.logical(accepted) || length(accepted)!=reference$n ||
      !identical(names(accepted),reference$sample_ids) || anyNA(accepted) || !all(accepted))) stop("HDI_5 requires accepted fits for entire frozen reference")
  z <- sweep(sweep(raw,2,reference$center,"-"),2,reference$sample_sd,"/")
  data.frame(sample_id=sample_ids,z,score=rowMeans(z),definition=definition,row.names=NULL)
}
