#' Plot cosine distance and matrix condition number for variable selection
#' 
#' Cosine distance based separation on CAP axes and matrix condition number are graphed as a function of Fatty Acids, where fatty acids are ordered by their relative contribution to Canonical axes in a CAP (Canonical Analysis of Principal Coordinates).
#' 
#' @param prey.mat (row-wise) dataframe of prey compositions
#' @param prey.ix Indexes rows of \code{prey.mat}
#' 
#' @return Two plots with 1) relative contributions of Fatty Acids to source separation and 2) matrix condition number.
#' @author Philipp Neubauer
#' @references Neubauer.P. and Jensen, O.P. (in prep)
#' #' @seealso \code{\link{add_FA}},\code{\link{select_vars}},\code{\link{run_MCMC}},\code{\link{simulation}}
#' @examples  \dontrun{
#' # load simulated example
#' data('Sim')
#' preys <- datas$datas.FA$preys
#' prey.ix <- datas$prey.ix
#' var_select_plot(preys,prey.ix)
#' }
#' @export
var_select_plot <- function(prey.mat,prey.ix){
  
  n.fats=ncol(prey.mat)
  distan <- cdist(prey.mat)
  PR.RDA <- vegan::capscale(as.dist(distan)~as.factor(prey.ix),comm=(prey.mat))

  par(mfcol=c(2,1),mar=c(6,6,4,4))
  sv = sort(clo(rowSums(sqrt(t(t(PR.RDA$CCA$v)*PR.RDA$CCA$eig)^2))),decreasing =T,index.return=T)
  plot(cumsum(sort(clo(rowSums(sqrt(t(t(PR.RDA$CCA$v)*PR.RDA$CCA$eig)^2))),decreasing =T)),axes=F,xlab='',ylab=expression(Delta),ylim=c(0,1),pch=16)
  axis(2,at=seq(0,1,0.2))
  axis(1,at=1:n.fats,labels=F)
  text(1:n.fats, par("usr")[3] - 0.2, srt = 45, adj = 1,
       labels = colnames(prey.mat)[sv$ix], xpd = TRUE)
  text(n.fats/2, par("usr")[3] - 0.6,labels='Fatty Acids', xpd = TRUE)
  
  
  ks <- kappas(prey.mat,sv$ix)
  plot(log(ks),axes=F,ylab=expression(log(kappa)),xlab='',ylim=c(0,max(log(ks))),pch=16)
  axis(2)
  axis(1,at=1:n.fats,labels=F)
  text(1:n.fats, par("usr")[3] - max(log(ks))/5 , srt = 45, adj = 1,
       labels = colnames(prey.mat)[sv$ix], xpd = TRUE)
  text(n.fats/2, par("usr")[3] - 4,labels='Fatty Acids', xpd = TRUE)
    
  cumsums <- as.data.frame(matrix(,n.fats,1))
  cumsums[,1] <- cumsum(sort(clo(rowSums(sqrt((t(t(PR.RDA$CCA$v)*PR.RDA$CCA$eig)^2)))),decreasing =T))
  cumsums[,2] <- kappas(prey.mat,sv$ix)
  names(cumsums) <- c('Cumulative Separation','Prey Matrix Condition')
  rownames(cumsums) <-  colnames(prey.mat)[sv$ix]
  print(cumsums)
  
  par(mfcol=c(1,1))
  return(sv)
}

#' Combine MCMC chains for prey items
#' 
#' Prey items can be combined into groups post hoc by summing their individual MCMC draws
#' 
#' @param mcmc output from \code{\link{run_MCMC}}
#' @param index A vector of integer indices or prey names to combine.
#' 
#' @return An object of the same class as the input with combined MCMC chains.
#' 
#' @author Philipp Neubauer
#' @references Neubauer.P. and Jensen, O.P. (in prep)
#' @seealso \code{\link{run_MCMC}}
#' @export
pcombine <- function(mcmc,index){
  if(is.character(index[1])){
    index <- match(index,mcmc$prey.names)
  }
  
  nindex <- (1:length(mcmc$prey.names))[-index]
  
  for (i in 1:mcmc$nChains){
    mcmc[[i]] <- cbind(rowSums(mcmc[[i]][,index]),mcmc[[i]][,nindex])
  }
  
  mcmc$prey.names <- c('Combined prey',mcmc$prey.names[nindex])
  
  return(mcmc)
}

#' Prints cosine distance and matrix condition number for variable selection
#' 
#' Cosine distance based separation on CAP axes and matrix condition number are printed as a function of Fatty Acids, where fatty acids are ordered by their relative contribution to Canonical axes in a CAP (Canonical Analysis of Principal Coordinates). The original indices of these fatty acids in the data are given to facilitate selection
#' 
#' @param prey.mat (row-wise) dataframe of prey compositions
#' @param prey.ix Indexes rows of \code{prey.mat}
#' 
#' @return Prints a matrix in the console with columns being 1) relative contributions of Fatty Acids to source separation, 2) matrix condition number and 3) original index used to select the fatty acids.
#' 
#' @author Philipp Neubauer
#' @references Neubauer.P. and Jensen, O.P. (in prep)
#' #' @seealso \code{\link{add_FA}},\code{\link{select_vars}},\code{\link{run_MCMC}},\code{\link{simulation}}
#' @examples  \dontrun{
#' # load simulated example
#' data('Sim')
#' preys <- datas$datas.FA$preys
#' prey.ix <- datas$prey.ix
#' print_var_list(preys,prey.ix)
#' }
#' @export
print_var_list <- function(prey.mat,prey.ix){
  
  n.fats=ncol(prey.mat) 

  distan <- cdist(prey.mat)
  PR.RDA <- vegan::capscale(distan~as.factor(prey.ix),comm=(prey.mat))
  
  sv = sort(clo(rowSums(sqrt((t(t(PR.RDA$CCA$v)*PR.RDA$CCA$eig)^2)))),decreasing =T,index.return=T)
        
  cumsums <- as.data.frame(matrix(,n.fats,2))
  cumsums[,1] <- cumsum(sort(clo(rowSums(sqrt((t(t(PR.RDA$CCA$v)*PR.RDA$CCA$eig)^2)))),decreasing =T))
  cumsums[,2] <- kappas(prey.mat,sv$ix)
  names(cumsums) <- c('Cumulative Separation','Prey Matrix Condition')
  rownames(cumsums) <-  colnames(prey.mat)[sv$ix]
  print(cumsums)

  return(sv)
}


#' Select Fatty Acids from imported data (should really be a method for subset...)
#' 
#' @param datas A data structure produced by \code{\link{add_SI}} and \code{\link{add_FA}}
#' @param ix Index of variables to retain, if not supplied, \code{\link{var_select_plot}} is called to itneractively select variables based on cosine distance and prey matrix condition.
#' @param plot Logical (T/F): should selection be itneractive? Needs tcl-tk package(and X-quartz on Mac)
#' @return a data structure of the same form as datas, with Fatty Acids selected by ix.
#' @author Philipp Neubauer
#' @references Neubauer.P. and Jensen, O.P. (in prep)
#' @seealso \code{\link{add_FA}},\code{\link{var_select_plot}},\code{\link{run_MCMC}},\code{\link{simulation}}
#' @examples \dontrun{
#' # load simulated example
#' data('Sim')
#' select_vars(datas,4:1)
#' # or select visualy
#' select_vars(datas)
#' }
#' @export
select_vars <- function(datas,ix=NULL,plot=T){
  
  # check if GUI is being used
  if(exists('GUI',envir=.GlobalEnv)){
    GUI <- get('GUI',envir=.GlobalEnv)
    if (GUI) datas <- guiGetSafe('datas')
  } else {
    GUI=F
  }
  
  colnames(datas$datas.FA$preys) <- sub('X(*)','\\1',colnames(datas$datas.FA$preys))
  
  colnames(datas$datas.FA$tau_c) <- colnames(datas$datas.FA$preys)
  
  colnames(datas$datas.FA$mean_c)  <- colnames(datas$datas.FA$tau_c)
  
  prey.mat <- clo(datas$datas.FA$preys*((datas$datas.FA$mean_c/datas$datas.FA$tau_c)[datas$prey.ix,]))
  
  n.fats=nrow(prey.mat)
  prey.ix <- datas$prey.ix
  
  if(is.null(ix)){
    if(plot==T){
      sv <- var_select_plot(prey.mat,prey.ix)
    } else {
      sv <- print_var_list(prey.mat,prey.ix)
    }
    nv <- select.list(title='please choose the fatty acids to use (at least 3)',choices = colnames(prey.mat)[sv$ix],graphics=F,multiple=T)
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
  datas$datas.FA$ni <- ni
  
  if(GUI & dev.cur()!=1) dev.off()
  
  ifelse(GUI,guiSet('datas',datas),return(datas))
  
}