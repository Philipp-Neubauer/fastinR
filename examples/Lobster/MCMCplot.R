MCMCplot <- function(res){
  
  class(res) <- 'mcmc.list'
  plot(res)
  
}