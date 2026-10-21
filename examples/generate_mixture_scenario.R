# Count-level synthetic controls with explicit latent memberships.
generate_mixture_scenario <- function(regime,seed) {
  regimes <- c("unimodal","balanced2","separated3","unbalanced")
  if(length(regime)!=1L || !regime %in% regimes) stop("Unknown scenario")
  if(length(seed)!=1L || !is.finite(seed)) stop("Finite seed required")
  set.seed(seed)
  r <- match(regime,regimes)
  probabilities <- list(1,c(.5,.5),rep(1/3,3),c(.9,.1))[[r]]
  offsets <- list(0,c(-2,2),c(-3,0,3),c(-2,2))[[r]]
  ng <- 1200L; np <- 16L
  meta <- expand.grid(participant_id=sprintf("sim_%02d",seq_len(np)),cycle=c("C0","T1","T2"),stringsAsFactors=FALSE)
  meta$sample_id <- sprintf("sim_sample_%03d",seq_len(nrow(meta)))
  meta$time <- match(meta$cycle,c("C0","T1","T2"))-1L
  membership <- sample(seq_along(probabilities),ng,replace=TRUE,prob=probabilities)
  abundance <- rnorm(ng,5.2,.25); gene_sd <- runif(ng,.15,.35)
  person <- matrix(rnorm(ng*np),ng)*gene_sd
  latent <- counts <- matrix(NA_real_,ng,nrow(meta),dimnames=list(sprintf("sim_gene_%04d",seq_len(ng)),meta$sample_id))
  for(j in seq_len(nrow(meta))) {
    p <- match(meta$participant_id[j],sprintf("sim_%02d",seq_len(np)))
    innovation <- rnorm(ng,0,gene_sd*(.8+.12*meta$time[j]))
    library <- runif(1,.8,1.3)
    latent[,j] <- abundance+offsets[membership]+person[,p]+innovation+log(library)
    counts[,j] <- rnbinom(ng,mu=exp(latent[,j]),size=12)
  }
  list(counts=counts,latent=latent,metadata=meta,
       genes=data.frame(gene_id=rownames(counts),group=membership),
       parameters=list(regime=regime,seed=seed,probabilities=probabilities,offsets=offsets,
                       log_abundance_sd=.25,gene_sd=c(.15,.35),NB_size=12))
}
