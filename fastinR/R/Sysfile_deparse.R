##### Pop props -----

#' @export
.Poppropanalysis <- function(datas) UseMethod(".Poppropanalysis", datas)

#' @export
.Poppropanalysis.FA <- function(datas)
{
  
  sysfile <- 'Pop_prop_analysis_FA'
  return(sysfile)
  
}

#' @export
.Poppropanalysis.SI <- function(datas)
{
  
  sysfile <- 'Pop_prop_analysis_SI'
  
  return(sysfile)
}

#' @export
.Poppropanalysis.combined <- function(datas)
{
    
  sysfile <- paste(system.file("exec",package = "fastinR"),"/Pop.prop.analysis.combined.bugs",sep='')
  
  return(sysfile)
  
}

#### Pop and Ind props ------
#' @export
.PopandIndprops <- function(datas) UseMethod(".PopandIndprops", datas)

#' @export
.PopandIndprops.FA <- function(datas)
{
  
  sysfile <- paste(system.file("exec",package = "fastinR"),"/Pop.and.Ind.props.FA.bugs",sep='')
  return(sysfile)
  
}

#' @export
.PopandIndprops.SI <- function(datas)
{
  
  sysfile <- paste(system.file("exec",package = "fastinR"),"/Pop.and.Ind.props.SI.bugs",sep='')
  return(sysfile)
}

#' @export
.PopandIndprops.combined <- function(datas)
{
  
  sysfile <- paste(system.file("exec",package = "fastinR"),"/Pop.and.Ind.props.combined.bugs",sep='')
  return(sysfile)
  
}

####### Analysis with Cov -------------
#' @export
.AnalysiswithCov <- function(datas) UseMethod(".AnalysiswithCov",datas)

#' @export
.AnalysiswithCov.FA <- function(datas)
{
  sysfile <- paste(system.file("exec",package = "fastinR"),"/Analysis.with.Cov.FA.bugs",sep='')
  return(sysfile)
  
}

#' @export
.AnalysiswithCov.SI <- function(datas)
{
  sysfile <- paste(system.file("exec",package = "fastinR"),"/Analysis.with.Cov.SI.bugs",sep='')
  return(sysfile)
}

#' @export
.AnalysiswithCov.combined <- function(datas)
{
  
  sysfile <- paste(system.file("exec",package = "fastinR"),"/Analysis.with.Cov.combined.bugs",sep='')
  return(sysfile)
  
}