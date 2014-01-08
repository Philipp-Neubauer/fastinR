plotvarselect <- function(prey_mat){
  
  n.fats=nrow(prey_mat)
  distan <- adist(prey_mat)
  RLR.RDA <- vegan::capscale(as.dist(distan)~as.factor(group),comm=t(prey_mat))
  par(mfcol=c(2,1))
  sv = sort(clo(rowSums(sqrt((t(t(RLR.RDA$CCA$v)*RLR.RDA$CCA$eig)^2)))),decreasing =T,index.return=T)
  plot(cumsum(sort(clo(rowSums(sqrt(t(t(RLR.RDA$CCA$v)*RLR.RDA$CCA$eig)^2))),decreasing =T)),axes=F,xlab='',ylab='Cumulative source separation',ylim=c(0,1))
  axis(2)
  axis(1,at=1:n.fats,labels=F)
  text(1:n.fats, par("usr")[3] - 0.2, srt = 45, adj = 1,
       labels = rownames(prey_mat)[sv$ix], xpd = TRUE)
  
  ks <- kappas(prey_mat,sv$ix)
  plot(ks,axes=F,ylab='Prey matrix condition number',xlab='',ylim=c(0,max(ks)))
  axis(2)
  axis(1,at=1:n.fats,labels=F)
  text(1:n.fats, par("usr")[3] - max(ks)/5 , srt = 45, adj = 1,
       labels = rownames(prey_mat)[sv$ix], xpd = TRUE)
  
  print(data.frame(cum_prop=cumsum(sort(clo(rowSums(sqrt(t(t(RLR.RDA$CCA$v)*RLR.RDA$CCA$eig)^2))),decreasing =T))))
  par(mfcol=c(1,1))
  return(sv)
}

selectvars <- function(prey.mat){
  n.fats=nrow(prey.mat)
  
  sv <- plotvarselect(prey.mat)
  
  six=1:n.fats
  nv <- select.list(title='please choose the fatty acids to use (at least 3)',choices = rownames(prey.mat)[sv$ix],graphics=T,multiple=T)
  
  #nv <- readline(prompt = "please enter number of variables for analysis \n")
  six <- match(nv,rownames(prey.mat))
  
  n.fats =length(six)
  list(selecta=six,n.fats=n.fats)
}