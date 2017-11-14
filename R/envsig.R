envsig <- function(envlist, adjustment = c("holm", "hochberg", "hommel",
                                            "bonferroni", "BH", "BY", "fdr", "none")){
  if (!is.list(envlist)) stop(
    "The method 'envsig' must be applied to an object from the output of envelope()")
  if (is.null(envlist$variogram0)) stop(
    "The method 'envsig' must be applied to an object from the output of envelope()")
  switch(class(envlist$variogram0)[1],
         "gstatVariogram"={
           dims <- dim(envlist$variogram)
           pvals <- sapply(1:dims[1], function(i){
             # number of more extreme semivariances
             min(table(envlist$variogram0$gamma[i] < envlist$variogram[i,]))/dims[2]
           })
         },
         "variogram"={
           dims <- dim(envlist$variogram)
           pvals <- sapply(1:dims[1], function(i){
             min(table(envlist$variogram0$v[i] < envlist$variogram[i,]))/dims[2]
           })
         })
  # minimum p-value is 1/nsim
  pvals[pvals==1 | pvals==0] <- 1/dims[2]
  pvals.adj <- stats::p.adjust(pvals, method = adjustment, n = length(pvals))
  # Fisher's method
  pvals.all <- pchisq( -2*sum(log(pvals.adj)), df = 2*length(pvals.adj), lower.tail=FALSE)
  return(list(p.pointwise = pvals.adj, p.overall = pvals.all))
}
