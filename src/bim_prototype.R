# Optional returned-parameter mixture fitting and convergence-checked selection.
.mixture_state_v1 <- function(x, means, sds, weights) {
  k <- length(weights)
  ld <- vapply(seq_len(k),function(j) log(weights[j]) + dnorm(x,means[j],sds[j],log=TRUE),numeric(length(x)))
  top <- if(k==2L) pmax(ld[,1],ld[,2]) else pmax(ld[,1],ld[,2],ld[,3])
  total <- top + log(rowSums(exp(ld-top)))
  list(loglik=sum(total),responsibility=exp(ld-total))
}
fit_bim_prototype <- function(x,k,max_steps=5000L,tolerance=1e-7,span=c(.15,.85)) {
  if(length(k)!=1L || !k %in% 2:3) stop("k must be 2 or 3")
  if(length(max_steps)!=1L || !is.finite(max_steps) || max_steps<1 || max_steps>5000 || max_steps!=as.integer(max_steps)) stop("max_steps must be an integer from 1 to 5000")
  if(length(tolerance)!=1L || !is.finite(tolerance) || tolerance<=0) stop("tolerance must be positive and finite")
  if(length(span)!=2L || any(!is.finite(span)) || any(span<0|span>1) || span[1]>=span[2]) stop("Invalid quantile span")
  invalid <- function(reason,steps=0L) list(k=as.integer(k),finite=FALSE,reason=reason,
    converged=FALSE,boundary=FALSE,steps=steps,criterion=Inf,loglik=NA_real_)
  if(!is.numeric(x) || length(x)<3L || any(!is.finite(x))) return(invalid("invalid_input"))
  if(length(unique(x))<k || sd(x)==0) return(invalid("degenerate_input"))
  means <- unname(quantile(x,seq(span[1],span[2],length.out=k),type=7))
  sds <- rep(sd(x),k); weights <- rep(1/k,k)
  state <- .mixture_state_v1(x,means,sds,weights)
  if(!is.finite(state$loglik)) return(invalid("invalid_initial_objective"))
  initial <- state$loglik; converged <- FALSE; min_change <- Inf
  for(step in seq_len(max_steps)) {
    r <- state$responsibility; mass <- colSums(r)
    if(any(!is.finite(mass)|mass<=0)) return(invalid("invalid_component_mass",step-1L))
    weights <- mass/length(x); means <- colSums(r*x)/mass
    sds <- pmax(1e-6,sqrt(colSums(r*outer(x,means,"-")^2)/mass))
    next_state <- .mixture_state_v1(x,means,sds,weights)
    if(!is.finite(next_state$loglik) || any(!is.finite(next_state$responsibility))) return(invalid("invalid_updated_objective",step))
    change <- next_state$loglik-state$loglik; min_change <- min(min_change,change)
    state <- next_state
    if(abs(change)<tolerance) { converged <- TRUE; break }
  }
  boundary <- any(sds<=1e-6)
  list(k=as.integer(k),finite=TRUE,weights=weights,means=means,sds=sds,
    initial_loglik=initial,loglik=state$loglik,objective_change=change,min_objective_change=min_change,
    criterion=(3*k-1)*log(length(x))-2*state$loglik,steps=as.integer(step),
    converged=converged,boundary=boundary,
    reason=if(!converged && boundary) "nonconverged_boundary" else if(!converged) "nonconverged" else if(boundary) "boundary" else "ok")
}
select_bim_prototype <- function(fits) {
  eligible <- which(vapply(fits,function(f) isTRUE(f$finite)&&is.finite(f$criterion),logical(1)))
  if(!length(eligible)) return(list(accepted=FALSE,reason="all_invalid",value=NA_real_,fit=NULL,candidates=fits))
  # Select first, then reject; a rejected winner never triggers alternate-k rescue.
  best <- fits[[eligible[which.min(vapply(fits[eligible],function(f) f$criterion,numeric(1)))]]]
  accepted <- isTRUE(best$converged) && !isTRUE(best$boundary)
  list(accepted=accepted,reason=best$reason,value=1-max(best$weights),fit=best,candidates=fits)
}
