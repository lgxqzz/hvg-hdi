#!/usr/bin/env Rscript
# Fixed synthetic success control, using the supported raw-count CLI unchanged.
script <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])
root <- dirname(dirname(normalizePath(script, mustWork=TRUE)))
source(file.path(root, "examples/generate_mixture_scenario.R"))
args <- commandArgs(TRUE)
if (length(args)>1L || any(startsWith(args, "--"))) {
  stop("Usage: Rscript examples/run_hdi5_demo.R [OUTPUT_PARENT]")
}
parent <- if (length(args)) args[1] else file.path(root, "examples/output/hdi5")
if (file.exists(parent) && !dir.exists(parent)) stop("Output parent must be a directory")
if (!dir.exists(parent) && !dir.create(parent, recursive=TRUE)) stop("Cannot create output parent")
out <- tempfile("separated3_", tmpdir=normalizePath(parent, mustWork=TRUE))
# Reserve the run directory before writing; a collision must never overwrite it.
if (!dir.create(out)) stop("Cannot create a fresh run directory")
message("HDI-5 demo output: ", out)
data_dir <- file.path(out, "inputs")
dir.create(data_dir)
scenario <- generate_mixture_scenario("separated3", 92301)
write.table(data.frame(gene_id=rownames(scenario$counts), scenario$counts, check.names=FALSE),
            file.path(data_dir, "counts.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(scenario$metadata, file.path(data_dir, "metadata.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
dput(scenario$parameters, file=file.path(data_dir, "scenario.R"))
status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(c(
  file.path(root, "src/run_pipeline.R"),
  "--counts", file.path(data_dir, "counts.tsv"),
  "--meta", file.path(data_dir, "metadata.tsv"),
  "--out_dir", file.path(out, "results"),
  "--input_scale", "raw_counts", "--profile", "tmm_v1", "--experimental_bim", "true")))
if (status != 0L) stop("HDI-5 demo pipeline failed; inspect output above")
h <- read.delim(file.path(out, "results/HDI_df.tsv"))
genes <- read.delim(file.path(out, "results/HVG_genes.tsv"))
message(sprintf("Samples: %d; HVGs: %d; finite HDI_4: %d/%d; finite HDI_5: %d/%d",
                nrow(h), nrow(genes), sum(is.finite(h$HDI_4)), nrow(h),
                sum(is.finite(h$HDI_5)), nrow(h)))
for (column in c("HDI_4_status", "HDI_5_status", "BIM_status")) {
  frequencies <- table(h[[column]], useNA="ifany")
  message(column, ": ", paste(names(frequencies), as.integer(frequencies), sep="=", collapse=", "))
}
if (!identical(h$sample_id, scenario$metadata$sample_id) ||
    !all(is.finite(h$HDI_4)) || !all(is.finite(h$HDI_5)) ||
    !isTRUE(all(h$HDI_4_status=="ok" & h$HDI_5_status=="ok" & h$BIM_status=="ok"))) {
  stop("HDI-5 success control did not produce valid scores for every sample; inspect retained outputs")
}
