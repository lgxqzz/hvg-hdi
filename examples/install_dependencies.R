# R >=4.3; install missing packages into your configured R library.
# Injectable operations allow failure tests without changing installed libraries.
install_dependencies <- function(
    available=function(package) requireNamespace(package,quietly=TRUE),
    install_cran=function(packages) utils::install.packages(packages,repos="https://cloud.r-project.org"),
    install_bioc=function(packages) BiocManager::install(packages,ask=FALSE,update=FALSE)) {
  cran <- c("ggplot2","testthat","BiocManager")
  required <- c(cran,"edgeR")
  missing_packages <- function(packages) packages[!vapply(packages,available,logical(1))]
  needed <- missing_packages(cran)
  if (length(needed)) install_cran(needed)
  if (!available("edgeR") && available("BiocManager")) install_bioc("edgeR")

  # install.packages can warn and return normally after a failed build.
  # Namespace loading checks runtime/demo/test imports, not just library folders.
  remaining <- missing_packages(required)
  if (length(remaining)) {
    stop("Dependency setup incomplete; missing or unloadable required packages: ",
         paste(remaining,collapse=", "),
         ". Resolve installation errors (including Linux system libraries) and rerun setup.",
         call.=FALSE)
  }
  message("Dependency setup complete: ",paste(required,collapse=", "))
  invisible(required)
}

if (sys.nframe()==0L) install_dependencies()
