##### Pop props -----

.Poppropanalysis <- function(datas) UseMethod(".Poppropanalysis", datas)

.Poppropanalysis.FA <- function(datas)
{
  
  sysfile <- paste(system.file("exec",package = "FASTIN"),"/Pop.prop.analysis.FA.bugs",sep='')
  return(sysfile)
  
}

.Poppropanalysis.SI <- function(datas)
{
  
  sysfile <- paste(system.file("exec",package = "FASTIN"),"/Pop.prop.analysis.SI.bugs",sep='')
  
  return(sysfile)
}

.Poppropanalysis.combined <- function(datas)
{
    
  sysfile <- paste(system.file("exec",package = "FASTIN"),"/Pop.prop.analysis.combined.bugs",sep='')
  
  return(sysfile)
  
}

#### Pop and Ind props ------

.PopandIndprops <- function(datas) UseMethod(".PopandIndprops", datas)

.PopandIndprops.FA <- function(datas)
{
  
  sysfile <- paste(system.file("exec",package = "FASTIN"),"/Pop.and.Ind.props.FA.bugs",sep='')
  return(sysfile)
  
}

.PopandIndprops.SI <- function(datas)
{
  
  sysfile <- paste(system.file("exec",package = "FASTIN"),"/Pop.and.Ind.props.SI.bugs",sep='')
  return(sysfile)
}

.PopandIndprops.combined <- function(datas)
{
  
  sysfile <- paste(system.file("exec",package = "FASTIN"),"/Pop.and.Ind.props.combined.bugs",sep='')
  return(sysfile)
  
}

####### Analysis with Cov -------------

.AnalysiswithCov <- function(datas) UseMethod(".AnalysiswithCov",datas)

.AnalysiswithCov.FA <- function(datas)
{
  sysfile <- paste(system.file("exec",package = "FASTIN"),"/Analysis.with.Cov.FA.bugs",sep='')
  return(sysfile)
  
}

.AnalysiswithCov.SI <- function(datas)
{
  sysfile <- paste(system.file("exec",package = "FASTIN"),"/Analysis.with.Cov.SI.bugs",sep='')
  return(sysfile)
}

.AnalysiswithCov.combined <- function(datas)
{
  
  sysfile <- paste(system.file("exec",package = "FASTIN"),"/Analysis.with.Cov.combined.bugs",sep='')
  return(sysfile)
  
}