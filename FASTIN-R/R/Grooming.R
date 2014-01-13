#' Plot cosine distance and matrix condition number for variable selection
#' 
#' Cosine distance and matrix condition number are graphed as a function of Fatty Acids, where fatty acids are ordered by their relative contribution to Canonical axes in a CAP (Canonical Analysis of Principal Coordinates).
#' 
#' @param prey_mat (row-wise) dataframe of prey compositions
#' @param prey.ix Indexes rows of \code{prey_mat}
#' 
#' @return Two plots with 1) relative contributions of Fatty Acids to source separation and 2) matrix condition number.
#' @author Philipp Neubauer
#' @references Neubauer.P. and Jensen, O.P. (in prep)
#' @export
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

#' Select Fatty Acids from imported data (should really be a method for subset...)
#' 
#' @param datas A data structure produced by \code{\link{addSI}} and \code{\link{addFA}}
#' @param ix Index of variables to retain, if not supplied, \code{\link{plotvarselect}} is called to itneractively select variables based on cosine distance and prey matrix condition.
#' 
#' @return a data structure of the same form as datas, with Fatty Acids selected by ix.
#' @author Philipp Neubauer
#' @references Neubauer.P. and Jensen, O.P. (in prep)
#' @export
selectvars <- function(datas,ix=NULL){
  
  # check if GUI is being used
  if(exists('GUI',envir=.GlobalEnv)){
    GUI <- get('GUI',envir=.GlobalEnv)
    if (GUI) datas <- guiGetSafe('datas')
  } else {
    GUI=F
  }
  
  prey.mat <- datas$datas.FA$preys
  
  n.fats=nrow(prey.mat)
  prey.ix <- datas$prey.ix
  
  if(is.null(ix)){
    sv <- plotvarselect(prey.mat,prey.ix)
    nv <- select.list(title='please choose the fatty acids to use (at least 3)',choices = colnames(prey.mat)[sv$ix],graphics=T,multiple=T)
  } else if(is.numeric(ix)){
    nv <- colnames(prey.mat)[ix]
  } else {nv <- ix}
  
  #nv <- readline(prompt = "please enter number of variables for analysis \n")
  six <- match(nv,colnames(prey.mat))
  
  n.fats  <- length(six)
  m.fats <- n.fats-1
  
  datas$datas.FA$preys <- clo(datas$datas.FA$preys[,six])
  mprey <- aggregate(datas$datas.FA$preys,list(prey.ix),gmean)[,2:(n.fats+1)]
  datas$datas.FA$preym <- alr(mprey)
  
  datas$datas.FA$preds.FA <- datas$datas.FA$preds.FA[,six]
  datas$datas.FA$preds <- alr(datas$datas.FA$preds.FA)
  datas$datas.FA$n.fats <- n.fats
  datas$datas.FA$m.fats <- m.fats
  
  datas$datas.FA$mean_c <- datas$datas.FA$mean_c[,six]
  datas$datas.FA$tau_c <- datas$datas.FA$tau_c[,six]
  
  n.preys <- datas$n.preys
  R <- array(,c(m.fats,m.fats,n.preys))
  ni<-rep(NA,n.preys)
  for (i in 1:n.preys){
    ni[i] <- max(n.fats+1,sum(prey.ix==unique(prey.ix)[i])-1)
    R[,,i]=cov(alr(prey.mat[prey.ix==unique(prey.ix)[i],six]))*ni[i]
  }
  
  datas$datas.FA$R <- R
  datas$datas.FA$Rnot <- datas$datas.FA$Rnot[six[1:m.fats],six[1:m.fats]]
  datas$datas.FA$ni <- ni
  
  if(GUI & dev.cur()!=1) dev.off()
  
  ifelse(GUI,guiSet('datas',datas),return(datas))
  
}