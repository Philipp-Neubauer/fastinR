# some convenience functions
kappas <- function(mat,vars){
  ks <- vector(,length(vars))
  for (i in 1:length(vars)){
    ks[i] <- kappa(mat[vars[1:i],],exact=T)
  }
  return(ks)
}

clo <- function(x){
  
  if (is.null(dim(x))){xc <- x/sum(x)} else {xc <- t(apply(x,1,function(y){y/sum(y)}))}
  return(xc)
}

alr <- function(x){
  x<-clo(x)
  if (is.null(dim(x))){xc <- log(x[1:(length(x)-1)]/x[length(x)])} else {xc <- t(apply(x,1,function(y){log(y[1:(length(y)-1)]/y[length(y)])}))}
  return(xc)
}

clr <- function(x){
  x<-clo(x)
  if (is.null(dim(x))){xc <- log(x)-mean(log(x))} else { t(apply(x,1,function(y){log(y)-mean(log(y))}))}
}

adist <- function(mat){
  
  dims <- dim(mat)
  dists <- matrix(,dims[2],dims[2])
  for (i in 1:(dims[2]-1)){
    for (j in (i+1):dims[2]){
      dists[j,i] <- 1-(t(mat[,i]) %*% mat[,j])/(sqrt(t(mat[,i]) %*% mat[,i])*sqrt(t(mat[,j]) %*% mat[,j]))
    }}
  dista <- as.dist(dists)
  return(dista)
}

gmean <- function(x) if (is.null(dim(x))) {exp(mean(log(x)))} else { t(apply(x,1,function(y){exp(mean(log(x)))}))}