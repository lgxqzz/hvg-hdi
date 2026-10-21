# R-only descriptive figures. Repeated observations are not independent inference.
plot_hdi_ridge <- function(hdi_df, metadata) {
  if (!requireNamespace("ggplot2",quietly=TRUE)) stop("Install ggplot2 for plots")
  if (anyDuplicated(metadata$sample_id) || !setequal(hdi_df$sample_id,metadata$sample_id)) stop("Plot metadata mismatch")
  d <- cbind(hdi_df,metadata[match(hdi_df$sample_id,metadata$sample_id),setdiff(names(metadata),"sample_id"),drop=FALSE])
  if (any(!is.finite(d$HDI_4)) || anyNA(d$cycle)) stop("Invalid plotting values")
  lev <- unique(as.character(d$cycle)); grid <- list()
  for (i in seq_along(lev)) {
    vals <- d$HDI_4[d$cycle==lev[i]]
    if (length(vals)<2 || sd(vals)==0) stop("Each ridge needs >=2 varying observations")
    den <- density(vals)
    grid[[i]] <- data.frame(x=den$x,y=i+den$y/max(den$y)*.72,base=i,cycle=lev[i])
  }
  g <- do.call(rbind,grid)
  pal <- setNames(rep(c("#527B99","#81A6AD","#B4C8B4"),length.out=length(lev)),lev)
  ggplot2::ggplot(g,ggplot2::aes(x=x,group=cycle,fill=cycle)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=base,ymax=y),colour="white",linewidth=.3) +
    ggplot2::geom_point(data=d,ggplot2::aes(x=HDI_4,y=match(cycle,lev)),inherit.aes=FALSE,shape=124,size=2,alpha=.7) +
    ggplot2::scale_fill_manual(values=pal) +
    ggplot2::scale_y_continuous(breaks=seq_along(lev),labels=vapply(lev,function(l) paste0(l,"  (n = ",sum(d$cycle==l),")"),character(1))) +
    ggplot2::labs(title="Simulated HDI-4 distributions",subtitle=paste(nrow(d),"simulated observations across",length(lev),"visits"),
                  x="HDI-4 (scoring-cohort standardization)",y=NULL,caption="Synthetic demonstration only. Ridges are descriptive densities;\nticks are observations.") +
    ggplot2::theme_classic(base_size=11) + ggplot2::theme(legend.position="none",plot.caption=ggplot2::element_text(hjust=0),plot.title=ggplot2::element_text(face="bold"))
}
plot_components <- function(hdi_df) {
  cols <- paste0("z_",c("IQR","MAD","ENT","Gini"))
  d <- data.frame(component=rep(sub("z_","",cols),each=nrow(hdi_df)),z=unlist(hdi_df[cols],use.names=FALSE))
  d$component <- factor(d$component,levels=sub("z_","",cols))
  ggplot2::ggplot(d,ggplot2::aes(component,z)) +
    ggplot2::geom_hline(yintercept=0,colour="grey80",linewidth=.4) +
    ggplot2::geom_violin(fill="#B3C8D6",colour="#527B99",linewidth=.4,trim=FALSE) +
    ggplot2::geom_boxplot(width=.12,outlier.shape=NA,fill="white",linewidth=.35) +
    ggplot2::labs(title="Four standardized components",subtitle=paste(nrow(hdi_df),"simulated samples; equal four-component mean"),
                  x=NULL,y="Component z score",caption="Synthetic demonstration only. Boxes: median and quartiles;\nwhiskers: 1.5 IQR. Each z score contributes z / 4 to HDI-4.") +
    ggplot2::theme_classic(base_size=11) + ggplot2::theme(plot.caption=ggplot2::element_text(hjust=0),plot.title=ggplot2::element_text(face="bold"))
}
