library(testthat)
candidate_fixture <- function() {
  set.seed(31); x <- matrix(rnorm(200*8),200,8)*rep(seq(.5,2,length.out=8),each=200)
  dimnames(x) <- list(paste0("g",1:200),paste0("s",1:8)); x
}
test_that("default is the exact four-component mean and never calls BIM", {
  x <- candidate_fixture()
  original <- fit_bim_prototype
  assign("fit_bim_prototype",function(...) stop("BIM called"),envir=.GlobalEnv)
  on.exit(assign("fit_bim_prototype",original,envir=.GlobalEnv))
  h <- compute_hdi4(x,rownames(x),"tmm_v1")
  expect_equal(h$HDI_4,unname(rowMeans(scale(raw4_components(x)))))
  expect_true(all(h$HDI_5_status=="not_requested"))
  expect_true(all(is.na(h$HDI_5)))
  expect_identical(h$sample_id,colnames(x))
  expect_false("HDI" %in% names(h))
  expect_true(all(h$HDI_4_version=="HDI4-v1"))
  expect_equal(nrow(attr(h,"bim_records")),0)
  expect_identical(attr(h,"reference_ids"),colnames(x))
  expect_error(compute_hdi4(x,c(rownames(x),"bad"),"tmm_v1"),"Missing selected")
  x[1,1] <- Inf; expect_error(compute_hdi4(x,rownames(x),"tmm_v1"),"finite numeric")
  x <- candidate_fixture(); x[,1] <- 1
  expect_error(compute_hdi4(x,rownames(x),"tmm_v1"),"Constant sample")
  x[,] <- 1:200
  expect_error(compute_hdi4(x,rownames(x),"tmm_v1"),"Undefined/constant")
})
test_that("optional scoring preserves all rows for each fit and reference failure", {
  x <- candidate_fixture(); base <- compute_hdi4(x,rownames(x),"tmm_v1")
  original <- fit_bim_prototype
  on.exit(assign("fit_bim_prototype",original,envir=.GlobalEnv))
  # Deterministic fit controls exercise integration, without a research fit sweep.
  for(kind in c("ok","nonconverged","boundary","constant")) {
    assign("fit_bim_prototype",local({type <- kind; function(v,k,...) {
      w <- if(type=="constant") .7 else .6+.2*plogis(mean(v))
      list(k=as.integer(k),finite=TRUE,criterion=k,weights=c(w,1-w),means=c(-1,1),sds=c(1,1),
        converged=type!="nonconverged",boundary=type=="boundary",steps=3L,
        reason=if(type %in% c("ok","constant")) "ok" else type)
    }}),envir=.GlobalEnv)
    h <- compute_hdi4(x,rownames(x),"tmm_v1",experimental_bim=TRUE)
    expect_identical(h$sample_id,base$sample_id)
    expect_identical(h$HDI_4,base$HDI_4)
    expect_equal(attr(h,"reference")[1:4,],attr(base,"reference"))
    expect_equal(nrow(attr(h,"bim_records")),16)
    if(kind=="ok") {
      expect_true(all(h$HDI_5_status=="ok"))
      expect_equal(h$HDI_5,rowMeans(cbind(as.matrix(h[paste0("z_",c("IQR","MAD","ENT","Gini"))]),h$z_BIM)))
    } else {
      expect_true(all(is.na(h$HDI_5)))
      expect_true(all(h$HDI_5_status==if(kind=="constant") "unavailable_constant_BIM" else "unavailable_rejected_fit"))
      if(kind!="constant") {
        expect_true(all(is.na(h$BIM)))
        expect_true(all(is.finite(attr(h,"bim_records")$diagnostic_BIM)))
        expect_true(is.na(tail(attr(h,"reference")$sample_sd,1)))
      }
    }
  }
})
test_that("actual optional one-step failure retains HDI4", {
  x <- candidate_fixture()
  h <- compute_hdi4(x,rownames(x),"tmm_v1",experimental_bim=TRUE,max_steps=1)
  expect_true(all(is.finite(h$HDI_4)))
  expect_true(all(is.na(h$HDI_5)))
  expect_true(all(h$HDI_5_status=="unavailable_rejected_fit"))
  expect_equal(nrow(h),8)
  expect_true(all(attr(h,"bim_records")$iterations==1))
})
test_that("CLI optional failure and provenance are complete and genuine", {
  td <- tempfile(); dir.create(td); generate_synthetic_data(file.path(td,"data"))
  out <- file.path(td,"out")
  args <- c(file.path(root,"src/run_pipeline.R"),"--counts",file.path(td,"data/counts.tsv"),
    "--meta",file.path(td,"data/metadata.tsv"),"--out_dir",out,"--input_scale","raw_counts","--profile","tmm_v1",
    "--experimental_bim","true","--bim_max_steps","1")
  log <- system2(file.path(R.home("bin"),"Rscript"),shQuote(args),stdout=TRUE,stderr=TRUE)
  expect_null(attr(log,"status"))
  h <- read.delim(file.path(out,"HDI_df.tsv"))
  expect_equal(nrow(h),48)
  expect_true(all(is.finite(h$HDI_4)))
  expect_true(all(is.na(h$HDI_5)))
  expect_true(all(h$HDI_5_status=="unavailable_rejected_fit"))
  records <- read.delim(file.path(out,"BIM_fit_records.tsv"))
  expect_equal(nrow(records),96)
  expect_identical(unique(records$sample_id),h$sample_id)
  expect_true(all(records$iterations==1))
  ref <- read.delim(file.path(out,"standardization.tsv"))
  expect_true(all(ref$n==48))
  expect_true(is.na(ref$sample_sd[5]))
  expect_identical(read.delim(file.path(out,"reference_samples.tsv"))$sample_id,h$sample_id)
  receipts <- read.delim(file.path(out,"output_identities.tsv"))
  expect_identical(unname(tools::md5sum(file.path(out,receipts$file))),receipts$md5)
  source_ids <- read.delim(file.path(out,"source_identities.tsv"))
  expect_identical(unname(tools::md5sum(file.path(root,"src",source_ids$file))),source_ids$md5)
  inputs <- read.delim(file.path(out,"input_identities.tsv"))
  expect_identical(unname(tools::md5sum(file.path(td,"data",c("counts.tsv","metadata.tsv")))),inputs$md5)
  conf <- dget(file.path(out,"config.R"))
  expect_identical(conf$score_version,"HDI4-v1")
  expect_identical(conf$bim_protocol,"returned-poststep-k23-v1")
  expect_identical(conf$bim_max_steps,1)
  expect_identical(conf$bim_sd_floor,1e-6)
})


test_that("actual accepted mixture controls produce separate HDI5 with unchanged HDI4", {
  # Same separated-distribution controls as the mixture API regression.
  set.seed(6102); two <- c(rnorm(144,-4,.5),rnorm(48,4,.8))
  set.seed(6103); three <- c(rnorm(96,-5,.4),rnorm(58,0,.6),rnorm(38,5,.5))
  x <- cbind(two=two,three=three); rownames(x) <- paste0("g",seq_len(192))
  base <- compute_hdi4(x,rownames(x),"tmm_v1")
  optional <- compute_hdi4(x,rownames(x),"tmm_v1",experimental_bim=TRUE)
  expect_identical(optional$sample_id,base$sample_id)
  expect_identical(optional$HDI_4,base$HDI_4)
  expect_equal(attr(optional,"reference")[1:4,],attr(base,"reference"))
  expect_true(all(optional$HDI_5_status=="ok"))
  expect_equal(optional$BIM,c(.25,.5),tolerance=1e-8)
  expect_equal(optional$HDI_5,(4*base$HDI_4+as.numeric(scale(optional$BIM)))/5)
  records <- attr(optional,"bim_records")
  expect_equal(nrow(records),4)
  expect_identical(records$k[records$accepted],c(2L,3L))
})


test_that("supported validation and reference preserve former meaningful checks", {
  x <- candidate_fixture()
  expect_error(compute_hdi4(x,rownames(x),"dip"),"Unknown profile")
  expect_error(compute_hdi4(x,rownames(x),"tmm_v1",nbins=1),"nbins")
  h <- compute_hdi4(x,rownames(x),"tmm_v1")
  raw <- as.matrix(h[c("IQR","MAD","ENT","Gini")])
  expect_equal(attr(h,"reference")$center,unname(colMeans(raw)))
  expect_equal(attr(h,"reference")$sample_sd,unname(apply(raw,2,sd)))
  expect_equal(unname(as.matrix(h[paste0("z_",colnames(raw))])),unname(sweep(sweep(raw,2,colMeans(raw),"-"),2,apply(raw,2,sd),"/")))
})
