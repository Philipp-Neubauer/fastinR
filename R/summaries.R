summary.pop_props <- function(MCMCout){
  res <- apply(MCMCout$MCMC,2,function(x){quantile(x,c(0.025,0.5,0.975))})
  #if (ls(.GlobalEnv)[grep('preys.names',ls(.GlobalEnv))[1]] == "preys.names"){
  prey.names <- unique(guiGetSafe('datas')$prey.ix)
  colnames(res) <- prey.names
  cat('\n','\n','\n',"population diet proportions",'\n','\n','\n')
  print(res)
}
summary.ind_props <- function(MCMCout){
  
  # print pop proportions
  popix <- grep('pop',colnames(MCMCout$MCMC))
  res <- apply(MCMCout$MCMC[,popix],2,function(x){quantile(x,c(0.025,0.5,0.975))})
    
  prey.names <- unique(guiGetSafe('datas')$prey.ix)
  colnames(res) <- prey.names
  
  cat('\n','\n','\n',"population diet proportions",'\n','\n','\n')
  print(res)
  
  nms <- colnames(res)
  n.p=length(nms)
  # print ind proportions
 indix <- 1:ncol(MCMCout$MCMC) 
  indix=indix[-which(indix%in%popix)]
  
  k=1;cc=0;
  for (a in 1:(length(indix)/n.p)){
  cc=cc+1
  res <- apply(MCMCout$MCMC[,indix[k:(k+n.p-1)]],2,function(x){quantile(x,c(0.025,0.5,0.975))})
  
  colnames(res) <-  nms 
  cat('\n','\n',"diet proportions for predator ",paste(cc),'\n','\n')
  print(res)
  
  k=k+3
  }
}
summary.cov_props <- function(MCMCout){}