# Baseline mean-window / MAD selection; inclusive boundaries retain ties.
compute_hvg_baseline <- function(counts, meta, id_col="sample_id", cycle_col="cycle",
                                 baseline="C0", mean_q_lo=.30, mean_q_hi=.70,
                                 mad_q=.60, out_file=NULL) {
  meta <- align_metadata(counts, meta, id_col, cycle_col)
  for (nm in c("mean_q_lo", "mean_q_hi", "mad_q")) check_probability(get(nm), nm)
  if (mean_q_lo >= mean_q_hi) stop("Mean quantile bounds must increase")
  if (length(baseline) != 1L || is.na(baseline) || !nzchar(baseline)) stop("Invalid baseline")
  ids <- meta[[id_col]][meta[[cycle_col]] == baseline]
  if (length(ids) < 3L) stop("At least three baseline samples required")
  x <- counts[, ids, drop=FALSE]
  means <- rowMeans(x)
  mads <- apply(x, 1L, stats::mad, constant=1.4826)
  limits <- stats::quantile(means, c(mean_q_lo, mean_q_hi), type=7)
  mid <- means >= limits[1] & means <= limits[2]
  if (!any(mid)) stop("Empty mean window")
  cutoff <- stats::quantile(mads[mid], mad_q, type=7)
  genes <- rownames(x)[mid & mads >= cutoff]
  if (length(genes) < 3L) stop("Fewer than three selected genes")
  tab <- data.frame(gene_id=genes, baseline_mean=means[genes], baseline_MAD=mads[genes], row.names=NULL)
  if (!is.null(out_file)) write_tsv(tab, out_file)
  attr(genes, "selection") <- list(mean_limits=unname(limits), mad_cutoff=unname(cutoff), baseline_ids=ids)
  genes
}
