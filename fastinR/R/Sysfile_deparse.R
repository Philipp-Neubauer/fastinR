##### Pop props -----

#' @export
.Poppropanalysis <- function(datas) UseMethod(".Poppropanalysis", datas)

#' @S3method .Poppropanalysis FA
.Poppropanalysis.FA <- function(datas)
{
  
  sysfile <- paste(system.file("exec",package = "fastinR"),"/Pop.prop.analysis.FA.bugs",sep='')
  return(sysfile)
  
}

#' @S3method .Poppropanalysis SI
.Poppropanalysis.SI <- function(datas)
{
  
  sysfile <- paste(system.file("exec",package = "fastinR"),"/Pop.prop.analysis.SI.bugs",sep='')
  
  return(sysfile)
}

#' @S3method .Poppropanalysis combined
.Poppropanalysis.combined <- function(datas)
{
    
  sysfile <- paste(system.file("exec",package = "fastinR"),"/Pop.prop.analysis.combined.bugs",sep='')
  
  return(sysfile)
  
}

#### Pop and Ind props ------
#' @export
.PopandIndprops <- function(datas) UseMethod(".PopandIndprops", datas)

#' @S3method .PopandIndprops FA
.PopandIndprops.FA <- function(datas)
{
  
  sysfile <- paste(system.file("exec",package = "fastinR"),"/Pop.and.Ind.props.FA.bugs",sep='')
  return(sysfile)
  
}

#' @S3method .PopandIndprops SI
.PopandIndprops.SI <- function(datas)
{
  
  sysfile <- paste(system.file("exec",package = "fastinR"),"/Pop.and.Ind.props.SI.bugs",sep='')
  return(sysfile)
}

#' @S3method .PopandIndprops combined
.PopandIndprops.combined <- function(datas)
{
  
  sysfile <- paste(system.file("exec",package = "fastinR"),"/Pop.and.Ind.props.combined.bugs",sep='')
  return(sysfile)
  
}

####### Analysis with Cov -------------
#' @export
.AnalysiswithCov <- function(datas) UseMethod(".AnalysiswithCov",datas)

#' @S3method .AnalysiswithCov FA 
.AnalysiswithCov.FA <- function(datas)
{
  sysfile <- paste(system.file("exec",package = "fastinR"),"/Analysis.with.Cov.FA.bugs",sep='')
  return(sysfile)
  
}

#' @S3method .AnalysiswithCov SI
.AnalysiswithCov.SI <- function(datas)
{
  sysfile <- paste(system.file("exec",package = "fastinR"),"/Analysis.with.Cov.SI.bugs",sep='')
  return(sysfile)
}

#' @S3method .AnalysiswithCov combined
.AnalysiswithCov.combined <- function(datas)
{
  
  sysfile <- paste(system.file("exec",package = "fastinR"),"/Analysis.with.Cov.combined.bugs",sep='')
  return(sysfile)
  
}