# HDI4-v1 scoring with separate optional HDI5-experimental-v1.
raw4_components <- function(x, nbins=50L) {
  if (!is.matrix(x) || !is.numeric(x) || nrow(x)<3L || any(!is.finite(x))) stop("Finite numeric profiles with >=3 features required")
  if (length(nbins)!=1L || !is.finite(nbins) || nbins<2 || nbins!=as.integer(nbins)) stop("nbins must be an integer >=2")
  if (any(apply(x,2,sd)==0)) stop("Constant sample profile: Gini undefined")
  raw <- cbind(IQR=apply(x,2,IQR),MAD=apply(x,2,mad,constant=1.4826),
               ENT=apply(x,2,.shannon_entropy,nbins=nbins),Gini=apply(x,2,.gini_coef))
  if (any(!is.finite(raw))) stop("Undefined raw4 component")
  raw
}
.bim_record <- function(d, id) {
  do.call(rbind,lapply(d$candidates,function(f) {
    selected <- !is.null(d$fit) && identical(f$k,d$fit$k)
    data.frame(sample_id=id,k=f$k,selected=selected,finite=f$finite,
      criterion=f$criterion,converged=f$converged,iterations=f$steps,boundary=f$boundary,
      reason=f$reason,accepted=selected && d$accepted,
      diagnostic_BIM=if(isTRUE(f$finite)) 1-max(f$weights) else NA_real_,
      weights=if(isTRUE(f$finite)) paste(format(f$weights,digits=17),collapse=";") else NA_character_,
      means=if(isTRUE(f$finite)) paste(format(f$means,digits=17),collapse=";") else NA_character_,
      sds=if(isTRUE(f$finite)) paste(format(f$sds,digits=17),collapse=";") else NA_character_)
  }))
}
compute_hdi4 <- function(counts,fgenes,profile,nbins=50L,experimental_bim=FALSE,
                         max_steps=5000L,tolerance=1e-7,out_file=NULL) {
  validate_matrix(counts); profile_settings(profile); validate_ids(fgenes,"Selected genes")
  if(length(fgenes)<3L || !all(fgenes %in% rownames(counts))) stop("Missing selected features or fewer than three genes")
  if(!is.logical(experimental_bim) || length(experimental_bim)!=1L || is.na(experimental_bim)) stop("experimental_bim must be logical")
  x <- counts[fgenes,,drop=FALSE]; ids <- colnames(x)
  raw <- raw4_components(x,nbins)
  ref <- freeze_dispersion_reference(raw,ids)
  h4 <- score_experimental_dispersion(raw,ids,ref,"HDI_4")
  z <- as.matrix(h4[,c("IQR","MAD","ENT","Gini")]); colnames(z) <- paste0("z_",colnames(raw))
  ans <- data.frame(sample_id=ids,raw,z,HDI_4=h4$score,HDI_5=NA_real_,
    HDI_4_version="HDI4-v1",HDI_4_status="ok",HDI_5_version="HDI5-experimental-v1",
    HDI_5_status="not_requested",BIM=NA_real_,z_BIM=NA_real_,BIM_status="not_requested",
    BIM_method="not_requested",profile=profile,row.names=NULL)
  reference <- data.frame(component=names(ref$center),center=ref$center,sample_sd=ref$sample_sd,n=ref$n,row.names=NULL)
  # Empty, schema-stable table in four-only mode: never invoke a mixture fitter.
  records <- data.frame(sample_id=character(),k=integer(),selected=logical(),finite=logical(),
    criterion=numeric(),converged=logical(),iterations=integer(),boundary=logical(),reason=character(),
    accepted=logical(),diagnostic_BIM=numeric(),weights=character(),means=character(),sds=character())
  if(experimental_bim) {
    fits <- lapply(seq_len(ncol(x)),function(j) select_bim_prototype(lapply(2:3,function(k)
      fit_bim_prototype(x[,j],k,max_steps=max_steps,tolerance=tolerance))))
    records <- do.call(rbind,lapply(seq_along(ids),function(j) .bim_record(fits[[j]],ids[j])))
    accepted <- vapply(fits,function(d) isTRUE(d$accepted),logical(1))
    ans$BIM_status <- vapply(fits,`[[`,character(1),"reason")
    ans$BIM_method <- "returned-poststep-k23-v1"
    ans$BIM <- vapply(fits,function(d) if(d$accepted) d$value else NA_real_,numeric(1))
    # No partial reference: rejected-fit values appear only in diagnostic records.
    bcenter <- if(all(accepted)) mean(ans$BIM) else NA_real_
    bsd <- if(all(accepted)) sd(ans$BIM) else NA_real_
    reference <- rbind(reference,data.frame(component="BIM",center=bcenter,sample_sd=bsd,n=length(ids)))
    if(!all(accepted)) {
      ans$HDI_5_status <- "unavailable_rejected_fit"
    } else if(!is.finite(bsd) || bsd<=0) {
      ans$HDI_5_status <- "unavailable_constant_BIM"
    } else {
      ans$z_BIM <- (ans$BIM-bcenter)/bsd
      ans$HDI_5 <- (4*ans$HDI_4+ans$z_BIM)/5
      ans$HDI_5_status <- "ok"
    }
  }
  attr(ans,"reference") <- reference
  attr(ans,"reference_ids") <- ids
  attr(ans,"bim_records") <- records
  if(!is.null(out_file)) write_tsv(ans,out_file)
  ans
}
