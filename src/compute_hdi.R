# Shared profile validation and unchanged four-component mathematical helpers.
profile_settings <- function(profile) {
  switch(profile,
    tmm_v1=list(input_scale="tmm_logcpm"),
    vst_sensitivity_v1=list(input_scale="deseq2_vst"),
    stop("Unknown profile: ", profile))
}
.shannon_entropy <- function(x, nbins=50L) {
  h <- graphics::hist(x, breaks=nbins, plot=FALSE)
  p <- h$density / sum(h$density)
  -sum(p[p > 0] * log(p[p > 0]))
}
.gini_coef <- function(x) {
  y <- as.numeric(scale(x))
  y <- sort(y - min(y) + 1e-6)
  n <- length(y)
  (2 * sum(seq_len(n) * y) - (n+1) * sum(y)) / (n * sum(y))
}
