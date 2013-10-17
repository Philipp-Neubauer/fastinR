summary.pop_props <- function(object, ...){
 
  for (l in 1:length(object)){
      
      cat('\n','\n',"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",'\n')
      cat("Printing results for MCMC chain",l)
      cat('\n',"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",'\n')
      
      MCMCout <- object[[l]]
      res <- apply(MCMCout,2,function(x){quantile(x,c(0.025,0.5,0.975))})
                                       
      prey.names <- unique(guiGetSafe('datas')$prey.ix)
      colnames(res) <- prey.names
      cat('\n','\n','\n',"population diet proportions",'\n','\n','\n')
      print(res)
  }
}
summary.ind_props <- function(object, ...){
  
   for (l in 1:length(object)){
       cat('\n','\n',"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",'\n')
       cat("Printing results for MCMC chain",l)
       cat('\n',"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",'\n')
       MCMCout <- object[[l]]
                                        # print pop proportions
       popix <- grep('pop',colnames(MCMCout))
       res <- apply(MCMCout[,popix],2,function(x){quantile(x,c(0.025,0.5,0.975))})
    
       prey.names <- unique(guiGetSafe('datas')$prey.ix)
       colnames(res) <- prey.names
  
       cat('\n','\n','\n',"Posterior percentiles for",'\n',"population diet proportions",'\n','\n','\n')
       print(res)
  
       nms <- colnames(res)
       n.p=length(nms)
       np <- (ncol(MCMCout)-n.p)/n.p
  # print ind proportions
       indix <- 1:ncol(MCMCout) 
       indix=indix[-which(indix%in%popix)]
  
       k=0;
       for (a in 1:(length(indix)/n.p)){
           k=k+1
           res <- apply(MCMCout[,indix[seq(k,(n.p-1)*np+k,np)]],2,function(x){quantile(x,c(0.025,0.5,0.975))})
  
           colnames(res) <- nms 
           cat('\n','\n',"Posterior percentiles for",'\n',"diet proportions of predator ",paste(k),'\n','\n')
           print(res)
       }   
   }
}
summary.cov_props <- function(object, ...){

    for (l in 1:length(object)){
        cat('\n','\n',"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",'\n')
        cat("Printing results for MCMC chain",l)
        cat('\n',"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",'\n')
        MCMCout <- object[[l]]
                                        # print pop proportions
      datas <- guiGetSafe('datas')

      prey.names <- unique(datas$prey.ix)
  
      Covs <- guiGetSafe('Covs') 
      cidx <- apply(Covs,2,function(x){any(x!=0 & x!=1)})

                                        #number of groups and covariates
      nGr <- sum(cidx==F)
      nCs <- sum(cidx)    
      Csidx <- which(cidx)
      Gridx <- which(cidx==F)
      
      covnames <- colnames(Covs[Csidx])
      grnames <- colnames(Covs[Gridx])
      
      betaix <- grep('beta',colnames(MCMCout))
      
    # if there's covariates -> show effects (beta)
      if(nCs>0){
          eff <- apply(MCMCout[,betaix[(1:(nCs*(datas$n.preys-1))+nGr*(datas$n.preys-1))]],2,function(x){quantile(x,c(0.025,0.5,0.975))})
          k=0
          for (i in 1:nCs){

              cat('\n','\n',"Posterior percentiles for",'\n',"effect of ",covnames[i]," on odds ratios", '\n','\n')            
              this.eff <- eff[,(k+1):(k+datas$n.preys-1)]
              for (n in 1:ncol(this.eff)) colnames(this.eff)[n] <- paste(prey.names[n],'/',prey.names[datas$n.preys])
              print(this.eff)

              k=k+(datas$n.preys-1)
              
          }
      }

                                        # idem for Groups

      if(nGr>0){
          eff <- apply(MCMCout[,betaix[(1:(nGr*(datas$n.preys-1)))]],2,function(x){quantile(x,c(0.025,0.5,0.975))})
          k=(datas$n.preys-1)
          for (i in 2:nGr){
              
              cat('\n','\n',"Posterior percentiles for odds ratios",'\n',"of group (anova) contrast",grnames[i],' / ',grnames[1],'\n','\n')            
              this.eff <- eff[,(k+1):(k+datas$n.preys-1)]
              for (n in 1:ncol(this.eff)) colnames(this.eff)[n] <- paste(prey.names[n],'/',prey.names[datas$n.preys])
              print(this.eff)
              
              k=k+(datas$n.preys-1)
          }
          
                                        # do again to give group proportions
          k=0
          popix <- grep('pop',colnames(MCMCout))
          res <- apply(MCMCout[,popix],2,function(x){quantile(x,c(0.025,0.5,0.975))})
          for (i in 1:nGr){          
              
              this.eff <- res[,(k+1):(k+datas$n.preys)]
              for (n in 1:ncol(this.eff)) colnames(this.eff)[n] <- paste(prey.names[n])
              
              cat('\n','\n','\n',"Posterior percentiles for",'\n',"diet proportions of group ",i,'\n',
                  '(at mean values of continuous co-variates)','\n','\n')
              print(this.eff)
              k=k+(datas$n.preys)
          }
      }
      
 
  # print ind proportions
      indix <- (max(popix)+1):ncol(MCMCout)

      for (k in 1:datas$n.preds){
          
          res <- apply(MCMCout[,indix[seq(k,(datas$n.preys-1)*datas$n.preds+k,datas$n.preds)]],2,function(x){quantile(x,c(0.025,0.5,0.975))})
          
          colnames(res) <- prey.names
          cat('\n','\n',"Posterior percentiles for",'\n',"diet proportions of predator ",paste(k),'\n','\n')
          print(res)
          
      }
  }
}
