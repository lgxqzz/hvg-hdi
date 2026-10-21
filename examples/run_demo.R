#!/usr/bin/env Rscript
script <- sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))][1])
root <- dirname(dirname(normalizePath(script)))
source(file.path(root,"examples/generate_synthetic_data.R"))
source(file.path(root,"src/plot_hdi_ridge.R"))
args <- commandArgs(TRUE)
figures <- "--figures" %in% args
args <- args[args!="--figures"]
if (length(args)>1L || any(startsWith(args,"--"))) stop("Usage: Rscript examples/run_demo.R [OUTPUT_PARENT] [--figures]")
parent <- if(length(args)) args[1] else file.path(root,"examples/output")
if (file.exists(parent) && !dir.exists(parent)) stop("Output parent must be a directory")
if (!dir.exists(parent) && !dir.create(parent,recursive=TRUE)) stop("Cannot create output parent")
out <- tempfile("synthetic_",tmpdir=parent)
data_dir <- file.path(out,"inputs")
generate_synthetic_data(data_dir)
status <- system2(file.path(R.home("bin"),"Rscript"),shQuote(c(file.path(root,"src/run_pipeline.R"),
  "--counts",file.path(data_dir,"counts.tsv"),"--meta",file.path(data_dir,"metadata.tsv"),
  "--out_dir",file.path(out,"results"),"--input_scale","raw_counts","--profile","tmm_v1")))
if (status != 0) stop("Demo pipeline failed")
h <- read.delim(file.path(out,"results/HDI_df.tsv")); m <- read.delim(file.path(data_dir,"metadata.tsv"))
if (figures) {
dir.create(file.path(root,"figures"),showWarnings=FALSE)
for (nm in c("longitudinal","components")) {
  plot <- if(nm=="longitudinal") plot_hdi_ridge(h,m) else plot_components(h)
  stem <- file.path(root,"figures",paste0("synthetic_",nm))
  # PNG for README viewing; vector PDF for print. No identifying PDF metadata.
  ggplot2::ggsave(paste0(stem,".png"),plot,device=grDevices::png,
                 width=7.2,height=4.2,dpi=300,bg="white")
  ggplot2::ggsave(paste0(stem,".pdf"),plot,device=grDevices::pdf,
                 width=7.2,height=4.2,bg="white",family="Helvetica",
                 useDingbats=FALSE,timestamp=FALSE,producer=FALSE,author="",
                 title=paste("Simulated",nm,"dispersion"))
}
}
message("Synthetic demo output: ",out)
