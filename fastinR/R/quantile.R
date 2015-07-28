#' @S3method quantiles ind_props
# Returns data frame with one row for each predator/prey combination, 
# for each chain.  
# The first n.prey rows in each chain refer to the population; 
# subsequent rows refer to individual predators
object <- ind.siandfa.pelagicfishandcephalopods
quantile.ind_props <- function(object,...){
  
  prey.names <- object$prey.names
  n.prey <- length(prey.names)
  
  for (l in 1:object$nChains) { # for each chain (l)
    MCMCout <- object[[l]]
    n.preds <- (ncol(MCMCout)-n.prey)/n.prey #number of predators
    
    # calculate pop proportions
    popix <- grep('pop',colnames(MCMCout))
    pop.res <- apply(MCMCout[,popix],2,function(x){quantile(x,c(0.025,0.5,0.975))})
    
    # add prey type and chain num
    pop.res <- cbind.data.frame(pred=rep(NA, n.prey), prey=prey.names, 
                                chain=rep(l, n.prey), t(pop.res))
    rownames(pop.res) <- paste0(rownames(pop.res), l) #add chain to row name to avoid duplicate rownames
    
    if(l==1) res <- pop.res else res <- rbind.data.frame(res, pop.res)
    rm(pop.res)
    
    # calculate ind proportions
    indix <- 1:ncol(MCMCout) 
    indix=indix[-which(indix%in%popix)]
    
    k=0;
    for (a in 1:(length(indix)/n.prey)){ # for each predator (a)
      k=k+1
      ind.res <- apply(MCMCout[,indix[seq(k,(n.p-1)*np+k,np)]],2,function(x){quantile(x,c(0.025,0.5,0.975))})
      ind.res <- cbind.data.frame(pred=rep(a, n.prey), prey=prey.names, 
                                  chain=rep(l, n.prey), t(ind.res))
      rownames(ind.res) <- paste0(rownames(ind.res), l) #add chain to row name to avoid duplicate rownames
      res <- rbind.data.frame(res, ind.res)
      rm(ind.res)
    } # for a
  } #for l
  return(res)
} #function

object <- readRDS("~/Dropbox/Rutgers-Dropbox/Magdalena Bay not shared/Talia code (graphs and stats)/MagBay analysis/data-to import/data for fastin/pop_siandfa.RDS")

#' @S3method quantile pop_props
# Returns data frame with one row for each pop/prey combination, for each chain.  
quantile.pop_props <- function(object,...){
  
  prey.names <- object$prey.names
  n.prey <- length(prey.names)
  
  for (l in 1:object$nChains){
    
    MCMCout <- object[[l]]
    chain.res <- apply(MCMCout,2,function(x){quantile(x,c(0.025,0.5,0.975))})
    
    chain.res <- cbind.data.frame(prey=prey.names, 
                                  chain=rep(l, n.prey), t(chain.res))
    rownames(chain.res) <- paste0(rownames(chain.res), l) #add chain to row name to avoid duplicate rownames
    
    if(l==1) res <- chain.res else res <- rbind.data.frame(res, chain.res)
    rm(chain.res)
  } # for l  
  return(res)
} # function