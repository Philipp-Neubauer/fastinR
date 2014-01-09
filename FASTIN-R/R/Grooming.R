plotvarselect <- function(prey_mat,prey.ix){
  
  n.fats=ncol(prey_mat)
  distan <- adist(t(prey_mat))
  PR.RDA <- vegan::capscale(as.dist(distan)~as.factor(prey.ix),comm=(prey_mat))
  X11()
  par(mfcol=c(2,1))
  sv = sort(clo(rowSums(sqrt((t(t(PR.RDA$CCA$v)*PR.RDA$CCA$eig)^2)))),decreasing =T,index.return=T)
  plot(cumsum(sort(clo(rowSums(sqrt(t(t(PR.RDA$CCA$v)*PR.RDA$CCA$eig)^2))),decreasing =T)),axes=F,xlab='',ylab='Cumulative source separation',ylim=c(0,1))
  axis(2)
  axis(1,at=1:n.fats,labels=F)
  text(1:n.fats, par("usr")[3] - 0.2, srt = 45, adj = 1,
       labels = colnames(prey_mat)[sv$ix], xpd = TRUE)
  
  ks <- kappas(prey_mat,sv$ix)
  plot(ks,axes=F,ylab='Prey matrix condition number',xlab='',ylim=c(0,max(ks)))
  axis(2)
  axis(1,at=1:n.fats,labels=F)
  text(1:n.fats, par("usr")[3] - max(ks)/5 , srt = 45, adj = 1,
       labels = colnames(prey_mat)[sv$ix], xpd = TRUE)
  
  
  cumsums <- as.data.frame(matrix(,n.fats,1))
  cumsums[,1] <- cumsum(sort(clo(rowSums(t(t(cbind(PR.RDA$CCA$v,PR.RDA$CA$v))*c(PR.RDA$CCA$eig,PR.RDA$CA$eig))^2)),decreasing =T))
  names(cumsums) <- 'Cumulative Proportion'
  rownames(cumsums) <-  colnames(prey_mat)[sv$ix]
  print(cumsums)
  
  par(mfcol=c(1,1))
  return(sv)
}

selectvars <- function(datas){
  
  # check if GUI is being used
  if(exists('GUI',envir=.GlobalEnv)){
    GUI <- get('GUI',envir=.GlobalEnv)
    if (GUI) datas <- guiGetSafe('datas')
  } else {
    GUI=F
  }
  
  prey.mat <- datas$datas.FA$preys
  
  n.fats=nrow(prey.mat)
  
  sv <- plotvarselect(prey.mat,datas$prey.ix)
  
  six=1:n.fats
  nv <- select.list(title='please choose the fatty acids to use (at least 3)',choices = colnames(prey.mat)[sv$ix],graphics=T,multiple=T)
  
  #nv <- readline(prompt = "please enter number of variables for analysis \n")
  six <- match(nv,colnames(prey.mat))
  
  n.fats =length(six)
  list(selecta=six,n.fats=n.fats)
  
  datas$datas.FA$preys <- clo(datas$datas.FA$preys[,six])
  datas$datas.FA$preds <- alr(datas$datas.FA$preds.FA[,six])
  datas$datas.FA$n.fats <- n.fats
  datas$datas.FA$preym <- clo(datas$datas.FA$preym[,six])
  datas$datas.FA$preds <- clo(datas$datas.FA$preds[,six])
  
  ifelse(GUI,guiSet('datas',datas),return(datas))
  
}