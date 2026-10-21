#!/usr/bin/env Rscript
# Locate modules relative to this script, independent of caller working directory.
main <- function() {
  script <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))][1])
  root <- dirname(dirname(normalizePath(script,mustWork=TRUE)))
  for (f in c("validation.R","compute_hvg.R","compute_hdi.R","normalize.R","bim_prototype.R","experimental_scores.R","candidate_scores.R")) source(file.path(root,"src",f),local=TRUE)
  a <- commandArgs(TRUE)
  if (identical(a,"--help")) { cat("Rscript src/run_pipeline.R --counts MATRIX.tsv --meta META.tsv --out_dir NEW_DIR --input_scale raw_counts|tmm_logcpm|deseq2_vst --profile tmm_v1|vst_sensitivity_v1 [--baseline C0 --id_col sample_id --cycle_col cycle --mean_q_lo .3 --mean_q_hi .7 --mad_q .6 --nbins 50 --experimental_bim false|true --bim_max_steps 5000 --bim_tolerance 1e-7]\n"); return(invisible(NULL)) }
  allowed <- c("counts","meta","out_dir","input_scale","profile","baseline","id_col","cycle_col","mean_q_lo","mean_q_hi","mad_q","nbins","experimental_bim","bim_max_steps","bim_tolerance")
  if (!length(a) || length(a) %% 2L != 0) stop("Arguments must be --key value pairs; use --help")
  keys <- sub("^--","",a[seq(1,length(a),2)])
  if (any(!startsWith(a[seq(1,length(a),2)],"--")) || any(!keys %in% allowed) || anyDuplicated(keys)) stop("Unknown, malformed or duplicate options")
  opts <- as.list(a[seq(2,length(a),2)]); names(opts) <- keys
  if (!all(c("counts","meta","out_dir","input_scale","profile") %in% keys)) stop("Missing required options; use --help")
  defaults <- list(baseline="C0",id_col="sample_id",cycle_col="cycle",mean_q_lo=.3,mean_q_hi=.7,mad_q=.6,nbins=50,experimental_bim="false",bim_max_steps=5000,bim_tolerance=1e-7)
  for (k in names(defaults)) if (is.null(opts[[k]])) opts[[k]] <- defaults[[k]]
  for (k in c("mean_q_lo","mean_q_hi","mad_q","nbins","bim_max_steps","bim_tolerance")) opts[[k]] <- as.numeric(opts[[k]])
  if (!opts$experimental_bim %in% c("true","false")) stop("experimental_bim must be true or false")
  if (dir.exists(opts$out_dir) || file.exists(opts$out_dir)) stop("Output path must be new to avoid stale or overwritten results")
  # Metadata are text: preserve IDs and extra fields without type inference.
  read_input <- function(f) read.delim(f,check.names=FALSE,stringsAsFactors=FALSE,
                                      colClasses="character",na.strings=c("NA",""))
  tab <- read.delim(opts$counts,check.names=FALSE,stringsAsFactors=FALSE,na.strings=c("NA",""),colClasses=c(gene_id="character"))
  if (ncol(tab) < 3 || names(tab)[1] != "gene_id") stop("Matrix TSV needs gene_id followed by numeric sample columns")
  validate_ids(tab[[1]],"Gene IDs")
  x <- as.matrix(tab[-1]); rownames(x) <- tab[[1]]
  meta <- align_metadata(x,read_input(opts$meta),opts$id_col,opts$cycle_col)
  norm <- normalize_input(x,opts$input_scale,opts$profile)
  genes <- compute_hvg_baseline(norm$matrix,meta,opts$id_col,opts$cycle_col,opts$baseline,opts$mean_q_lo,opts$mean_q_hi,opts$mad_q)
  result <- compute_hdi4(norm$matrix,genes,opts$profile,opts$nbins,
    experimental_bim=opts$experimental_bim=="true",max_steps=opts$bim_max_steps,tolerance=opts$bim_tolerance)
  reference <- attr(result,"reference")
  # Write only after all scientific calculations succeed.
  dir.create(opts$out_dir,recursive=TRUE)
  o <- function(f) file.path(opts$out_dir,f)
  write_tsv(data.frame(gene_id=as.character(genes)),o("HVG_genes.tsv"))
  write_tsv(result,o("HDI_df.tsv")); write_tsv(reference,o("standardization.tsv"))
  write_tsv(data.frame(sample_id=colnames(x),baseline=meta[[opts$cycle_col]]==opts$baseline),o("sample_roles.tsv"))
  write_tsv(meta,o("metadata.tsv"))
  write_tsv(attr(result,"bim_records"),o("BIM_fit_records.tsv"))
  write_tsv(data.frame(sample_id=attr(result,"reference_ids")),o("reference_samples.tsv"))
  write_tsv(data.frame(gene_id=norm$retained),o("retained_genes.tsv"))
  if (!is.null(norm$factors)) write_tsv(data.frame(sample_id=colnames(x),norm$factors,row.names=NULL),o("normalization_factors.tsv"))
  identities <- data.frame(role=c("counts","metadata"),md5=unname(tools::md5sum(c(opts$counts,opts$meta))))
  write_tsv(identities,o("input_identities.tsv"))
  # Exclude caller paths: configuration is portable and safe to share with synthetic data.
  config <- opts[setdiff(names(opts),c("counts","meta","out_dir"))]
  config$score_version <- "HDI4-v1"; config$optional_version <- "HDI5-experimental-v1"
  config$bim_protocol <- "returned-poststep-k23-v1"; config$bim_sd_floor <- 1e-6
  config$reference <- "all_input_samples_sample_SD"
  config$normalization <- norm$normalization; config$selection <- attr(genes,"selection")
  config$software <- lapply(c("edgeR","ggplot2"),function(p) if(requireNamespace(p,quietly=TRUE)) as.character(packageVersion(p)) else NA_character_)
  names(config$software) <- c("edgeR","ggplot2")
  dput(config,file=o("config.R"))
  write_tsv(data.frame(file=sort(list.files(file.path(root,"src"))),md5=unname(tools::md5sum(sort(list.files(file.path(root,"src"),full.names=TRUE))))),o("source_identities.tsv"))
  writeLines(capture.output(sessionInfo()),o("sessionInfo.txt"))
  writeLines(c("SUCCESS",paste("Samples",ncol(x)),paste("Baseline",sum(meta[[opts$cycle_col]]==opts$baseline)),
               paste("HVGs",length(genes)),paste("Profile",opts$profile),"HDI4-v1: four equal standardized components; complete scoring cohort reference",
               paste("HDI-5 status",unique(result$HDI_5_status))),o("analysis.log"))
  products <- sort(list.files(opts$out_dir,full.names=TRUE))
  write_tsv(data.frame(file=basename(products),md5=unname(tools::md5sum(products))),o("output_identities.tsv"))
  message("SUCCESS: ",ncol(x)," samples, ",length(genes)," HVGs; outputs in requested directory")
}
if (sys.nframe() == 0) main()
