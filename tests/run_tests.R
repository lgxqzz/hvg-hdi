#!/usr/bin/env Rscript
script <- sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))][1])
root <- dirname(dirname(normalizePath(script)))
if (!requireNamespace("testthat",quietly=TRUE)) stop("Install testthat")
for (f in c("validation.R","compute_hvg.R","compute_hdi.R","normalize.R","bim_prototype.R","experimental_scores.R","candidate_scores.R")) source(file.path(root,"src",f))
source(file.path(root,"examples/generate_synthetic_data.R"))
result <- testthat::test_file(file.path(root,"tests/test_workflow.R"),reporter="summary",stop_on_failure=TRUE)

prototype <- testthat::test_file(file.path(root,"tests/test_bim_prototype.R"),reporter="summary",stop_on_failure=TRUE)
experimental <- testthat::test_file(file.path(root,"tests/test_experimental_scores.R"),reporter="summary",stop_on_failure=TRUE)
candidate <- testthat::test_file(file.path(root,"tests/test_candidate_scores.R"),reporter="summary",stop_on_failure=TRUE)
d <- rbind(as.data.frame(candidate),as.data.frame(result),as.data.frame(prototype),as.data.frame(experimental))
cat(sprintf("RESULT: %d test blocks; %d expectations passed; %d failed; %d errors; %d warnings\n",nrow(d),sum(d$passed),sum(d$failed),sum(d$error),sum(d$warning)))
