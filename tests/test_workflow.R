library(testthat)
fixture <- function() {
  x <- matrix(rep(1:10,3),10,3) + outer(c(0,0,0:5,0,0),c(-1,0,1))
  x <- cbind(x,rep(100,10)); dimnames(x) <- list(paste0("g",1:10),paste0("s",1:4))
  list(x=x,m=data.frame(sample_id=colnames(x),cycle=c("C0","C0","C0","T1")))
}
test_that("HVG mean/MAD cutoffs have known values and exclude follow-up", {
  f <- fixture(); g <- compute_hvg_baseline(f$x,f$m,mean_q_lo=.2,mean_q_hi=.8)
  expect_identical(as.character(g),paste0("g",6:8))
  expect_equal(attr(g,"selection")$mean_limits,c(2.8,8.2))
  expect_equal(attr(g,"selection")$mad_cutoff,3*1.4826)
  f$x[,4] <- -1000
  expect_identical(g,compute_hvg_baseline(f$x,f$m,mean_q_lo=.2,mean_q_hi=.8))
  expect_identical(g,compute_hvg_baseline(f$x,f$m[4:1,],mean_q_lo=.2,mean_q_hi=.8))
})
test_that("ties and inclusive quantile boundaries retained", {
  f <- fixture(); f$x[,1:3] <- outer(1:10,c(1,1,1)) + rep(c(-1,0,1),each=10)
  expect_length(compute_hvg_baseline(f$x,f$m,mean_q_lo=0,mean_q_hi=1),10)
})
test_that("four raw components have independent numerical expectations", {
  v <- c(0,1,2,3)
  expect_equal(IQR(v),1.5)
  expect_equal(mad(v),1.4826)
  expect_equal(.shannon_entropy(v,2),-(.75*log(.75)+.25*log(.25)))
  y <- (v-mean(v))/sqrt(sum((v-mean(v))^2)/3)
  y <- y-min(y)+1e-6
  expect_equal(.gini_coef(v),sum(abs(outer(y,y,"-")))/(2*length(y)*sum(y)),tolerance=1e-12)
})
test_that("invalid inputs cannot be silently coerced, dropped or imputed", {
  f <- fixture()
  for (v in c(NA_real_,Inf,NaN)) {
    bad <- f$x; bad[1,1] <- v
    expect_error(compute_hvg_baseline(bad,f$m),"finite numeric")
  }
  bad <- f$x; rownames(bad)[2] <- rownames(bad)[1]; expect_error(validate_matrix(bad),"Gene IDs")
  bad <- f$x; colnames(bad)[2] <- colnames(bad)[1]; expect_error(validate_matrix(bad),"Sample IDs")
  bad <- f$x; rownames(bad)[1] <- ""; expect_error(validate_matrix(bad),"Gene IDs")
  expect_error(validate_matrix(matrix("1",3,3)),"finite numeric")
  expect_error(validate_ids(c("gene\t1","g2"),"test"),"control")
  m <- f$m; m$sample_id[1] <- "other"; expect_error(align_metadata(f$x,m,"sample_id","cycle"),"sets must match")
  m <- f$m; m$sample_id[1] <- m$sample_id[2]; expect_error(align_metadata(f$x,m,"sample_id","cycle"),"Metadata sample IDs")
  m <- f$m; m$cycle[1] <- NA; expect_error(align_metadata(f$x,m,"sample_id","cycle"),"Missing cycle")
  expect_error(compute_hvg_baseline(f$x,f$m,baseline="missing"),"three baseline")
  expect_error(compute_hvg_baseline(f$x,f$m,mean_q_lo=.9,mean_q_hi=.1),"increase")
  expect_error(compute_hvg_baseline(f$x,f$m,mad_q=1.1),"probability")
  expect_error(compute_hvg_baseline(f$x,f$m,mean_q_lo=.3,mean_q_hi=.31),"Empty mean window")
  expect_error(normalize_input(f$x,"raw_counts","vst_sensitivity_v1"),"incompatible")
  expect_error(normalize_input(f$x-.2,"raw_counts","tmm_v1"),"nonnegative integers")
})
test_that("raw-count normalization matches an independent edgeR recipe", {
  td <- tempfile(); generate_synthetic_data(td)
  t <- read.delim(file.path(td,"counts.tsv")); x <- as.matrix(t[-1]); rownames(x) <- t[[1]]
  d <- edgeR::DGEList(x); d <- d[rowMeans(edgeR::cpm(d)>=1)>=.2,,keep.lib.sizes=FALSE]
  d <- edgeR::calcNormFactors(d,method="TMM")
  expect_equal(normalize_input(x,"raw_counts","tmm_v1")$matrix,edgeR::cpm(d,log=TRUE,prior.count=1))
})
test_that("CLI is deterministic from root and outside cwd and writes genuine TSV", {
  td <- tempfile(); dir.create(td); data_dir <- file.path(td,"data"); generate_synthetic_data(data_dir)
  a <- c(file.path(root,"src/run_pipeline.R"),"--counts",file.path(data_dir,"counts.tsv"),"--meta",file.path(data_dir,"metadata.tsv"),"--input_scale","raw_counts","--profile","tmm_v1")
  run <- function(args) suppressWarnings(system2(file.path(R.home("bin"),"Rscript"),shQuote(args),stdout=TRUE,stderr=TRUE))
  old <- getwd(); on.exit(setwd(old)); setwd(root)
  out1 <- file.path(td,"out1"); log <- run(c(a,"--out_dir",out1))
  expect_null(attr(log,"status"))
  setwd(td); out2 <- file.path(td,"out2"); log <- run(c(a,"--out_dir",out2)); expect_null(attr(log,"status"))
  files <- c("HVG_genes.tsv","HDI_df.tsv","standardization.tsv","sample_roles.tsv","metadata.tsv","retained_genes.tsv","normalization_factors.tsv","input_identities.tsv","config.R","source_identities.tsv","sessionInfo.txt","analysis.log","output_identities.tsv")
  expect_true(all(file.exists(file.path(out1,files))))
  expect_identical(readLines(file.path(out1,"HDI_df.tsv")),readLines(file.path(out2,"HDI_df.tsv")))
  h <- read.delim(file.path(out1,"HDI_df.tsv")); expect_equal(nrow(h),48L)
  expect_true(all(is.finite(h$HDI_4))); expect_equal(ncol(h),20L)
  expect_true(grepl("\t",readLines(file.path(out1,"HDI_df.tsv"),n=1),fixed=TRUE))
  expect_equal(nrow(read.delim(file.path(out1,"HVG_genes.tsv"))),192L)
  expect_equal(nrow(read.delim(file.path(out1,"standardization.tsv"))),4L)
  expect_equal(attr(run(c(a,"--out_dir",out1)),"status"),1L)
  expect_equal(attr(run(c(a,"--out_dir",file.path(td,"bad"),"--unknown","1")),"status"),1L)
  expect_equal(attr(run(c(a,"--out_dir")),"status"),1L)
})

test_that("sensitivity profile is explicit and deterministic", {
  f <- fixture()
  expect_equal(normalize_input(f$x,"deseq2_vst","vst_sensitivity_v1")$matrix,f$x)
  expect_error(normalize_input(f$x,"deseq2_vst","tmm_v1"),"incompatible")
  expect_error(normalize_input(f$x*0,"raw_counts","tmm_v1"),"positive libraries")
})

test_that("CLI preserves text metadata IDs, ordering and numerical results", {
  td <- tempfile(); dir.create(td)
  data_dir <- file.path(td,"data"); generate_synthetic_data(data_dir)
  read_text <- function(path) read.delim(path,check.names=FALSE,colClasses="character")
  write_input <- function(x,path) write.table(x,path,sep="\t",row.names=FALSE,quote=TRUE,na="NA")
  run <- function(counts,meta,out,id_col="sample_id") {
    args <- c(file.path(root,"src/run_pipeline.R"),"--counts",counts,"--meta",meta,
              "--out_dir",out,"--input_scale","raw_counts","--profile","tmm_v1","--id_col",id_col)
    suppressWarnings(system2(file.path(R.home("bin"),"Rscript"),shQuote(args),stdout=TRUE,stderr=TRUE))
  }
  control <- file.path(td,"control")
  log <- run(file.path(data_dir,"counts.tsv"),file.path(data_dir,"metadata.tsv"),control)
  expect_null(attr(log,"status"))
  control_hdi <- read_text(file.path(control,"HDI_df.tsv"))
  counts <- read.delim(file.path(data_dir,"counts.tsv"),check.names=FALSE)
  ids <- sprintf("%03d",seq_len(ncol(counts)-1L))
  names(counts)[-1] <- ids
  count_path <- file.path(td,"padded-counts.tsv"); write_input(counts,count_path)
  for (id_col in c("sample_id","specimen_id")) {
    meta <- read_text(file.path(data_dir,"metadata.tsv"))
    meta$sample_id <- ids
    meta$participant_id <- sprintf("%03d",match(meta$participant_id,unique(meta$participant_id)))
    meta$external_id <- "00012345678901234567890"
    meta$numeric_annotation <- "01.50"
    names(meta)[names(meta)=="sample_id"] <- id_col
    meta_path <- file.path(td,paste0(id_col,".tsv"))
    # Reverse metadata rows to exercise ID matching rather than rowwise pairing.
    write_input(meta[nrow(meta):1,],meta_path)
    out <- file.path(td,id_col)
    log <- run(count_path,meta_path,out,id_col)
    expect_null(attr(log,"status"))
    saved <- read_text(file.path(out,"metadata.tsv"))
    hdi <- read_text(file.path(out,"HDI_df.tsv"))
    expect_identical(saved[[id_col]],ids)
    expect_identical(hdi$sample_id,ids)
    expect_identical(read_text(file.path(out,"sample_roles.tsv"))$sample_id,ids)
    expect_identical(saved,meta)
    expect_identical(saved$participant_id,meta$participant_id)
    expect_identical(saved$external_id,rep("00012345678901234567890",length(ids)))
    expect_identical(hdi[setdiff(names(hdi),"sample_id")],control_hdi[setdiff(names(control_hdi),"sample_id")])
    expect_identical(readLines(file.path(out,"HVG_genes.tsv")),readLines(file.path(control,"HVG_genes.tsv")))
    expect_identical(readLines(file.path(out,"standardization.tsv")),readLines(file.path(control,"standardization.tsv")))
  }
  # The text contract must not relax existing missing/blank/duplicate/mismatch checks.
  for (kind in c("duplicate","empty","missing","whitespace","unmatched")) {
    bad <- meta
    bad[[id_col]][1] <- switch(kind,duplicate=bad[[id_col]][2],empty="",missing=NA_character_,whitespace="   ",unmatched="999")
    bad_path <- file.path(td,paste0(kind,".tsv")); write_input(bad,bad_path)
    out <- file.path(td,paste0("invalid-",kind))
    log <- run(count_path,bad_path,out,id_col)
    expect_equal(attr(log,"status"),1L)
    expect_true(any(grepl(if(kind=="unmatched") "sets must match" else "Metadata sample IDs",log,fixed=TRUE)))
    expect_false(dir.exists(out))
  }
})

test_that("dependency setup is idempotent and fails after incomplete installs", {
  env <- new.env(parent=globalenv())
  source(file.path(root,"examples/install_dependencies.R"),local=env)
  required <- c("ggplot2","testthat","BiocManager","edgeR")
  unexpected_install <- function(packages) stop("Unexpected install")
  expect_message(value <- env$install_dependencies(function(package) TRUE,unexpected_install,unexpected_install),"Dependency setup complete")
  expect_identical(value,required)

  installed <- character(); calls <- list()
  available <- function(package) package %in% installed
  cran <- function(packages) { calls$cran <<- packages; installed <<- c(installed,packages) }
  bioc <- function(packages) { calls$bioc <<- packages; installed <<- c(installed,packages) }
  expect_message(env$install_dependencies(available,cran,bioc),"Dependency setup complete")
  expect_identical(calls$cran,c("ggplot2","testthat","BiocManager"))
  expect_identical(calls$bioc,"edgeR")
  expect_message(env$install_dependencies(available,unexpected_install,unexpected_install),"Dependency setup complete")

  expect_error(env$install_dependencies(function(package) package!="testthat",function(packages) NULL,unexpected_install),"missing or unloadable required packages: testthat")
  expect_error(env$install_dependencies(function(package) package!="edgeR",unexpected_install,function(packages) NULL),"missing or unloadable required packages: edgeR")
  expect_error(env$install_dependencies(function(package) FALSE,function(packages) NULL,unexpected_install),"ggplot2, testthat, BiocManager, edgeR")

  # Reproduce warning-only installation failure in a real Rscript process.
  # Injected operations leave installed package libraries unchanged.
  driver <- tempfile(fileext=".R")
  writeLines(c(
    sprintf("source(%s)",encodeString(file.path(root,"examples/install_dependencies.R"),quote='"')),
    'install_dependencies(available=function(package) package!="testthat", install_cran=function(packages) warning("simulated build failure"), install_bioc=function(packages) stop("Unexpected Bioconductor install"))'
  ),driver)
  log <- suppressWarnings(system2(file.path(R.home("bin"),"Rscript"),shQuote(driver),stdout=TRUE,stderr=TRUE))
  expect_equal(attr(log,"status"),1L)
  expect_true(any(grepl("simulated build failure",log,fixed=TRUE)))
  expect_true(any(grepl("missing or unloadable required packages: testthat",log,fixed=TRUE)))
  expect_false(any(grepl("Dependency setup complete",log,fixed=TRUE)))
})

test_that("routine Quick Start repeats with fresh outputs and preserves tracked figures", {
  td <- tempfile(); dir.create(td)
  assets <- list.files(file.path(root,"figures"),full.names=TRUE)
  before <- tools::md5sum(assets)
  run <- function() system2(file.path(R.home("bin"),"Rscript"),
    shQuote(c(file.path(root,"examples/run_demo.R"),td)),stdout=TRUE,stderr=TRUE)
  first <- run(); expect_null(attr(first,"status"))
  second <- run(); expect_null(attr(second,"status"))
  runs <- list.dirs(td,recursive=FALSE,full.names=TRUE)
  expect_length(runs,2L)
  expect_true(all(file.exists(file.path(runs,"inputs/counts.tsv"))))
  expect_true(all(file.exists(file.path(runs,"results/HDI_df.tsv"))))
  expect_identical(readLines(file.path(runs[1],"results/HDI_df.tsv")),
                   readLines(file.path(runs[2],"results/HDI_df.tsv")))
  expect_identical(tools::md5sum(assets),before)
})
