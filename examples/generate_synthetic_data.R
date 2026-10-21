# Simulated gene-body counts only: no measurements or parameters from patients.
generate_synthetic_data <- function(out_dir, seed=20261020L) {
  set.seed(seed)
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  n_genes <- 1200L; n_people <- 16L
  meta <- expand.grid(participant_id=sprintf("sim_%02d",seq_len(n_people)), cycle=c("C0","T1","T2"),stringsAsFactors=FALSE)
  meta$sample_id <- sprintf("sim_sample_%03d",seq_len(nrow(meta)))
  meta$time <- match(meta$cycle,c("C0","T1","T2"))-1L
  abundance <- exp(rnorm(n_genes,5.2,1))
  gene_sd <- runif(n_genes,.15,.9)
  person <- matrix(rnorm(n_genes*n_people),nrow=n_genes)*gene_sd
  x <- vapply(seq_len(nrow(meta)),function(j) {
    person_j <- match(meta$participant_id[j],sprintf("sim_%02d",seq_len(n_people)))
    spread <- runif(1,.6,1.5) + .12*meta$time[j]
    mu <- abundance * exp(person[,person_j] + rnorm(n_genes,0,gene_sd*spread)) * runif(1,.8,1.3)
    rnbinom(n_genes,mu=mu,size=12)
  },numeric(n_genes))
  colnames(x) <- meta$sample_id
  write.table(data.frame(gene_id=sprintf("sim_gene_%04d",seq_len(n_genes)),x,check.names=FALSE),file.path(out_dir,"counts.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  write.table(meta,file.path(out_dir,"metadata.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  invisible(meta)
}
