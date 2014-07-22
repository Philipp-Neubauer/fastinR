#' @S3method print Fatty_Acid_Profiles
print.Fatty_Acid_Profiles <- function(x, ...){
   cat('Fatty Acid profile data with:','\n')
   cat(datas$n.preds,' Predators','\n')
   cat(datas$n.preys,' Prey species','\n')
   cat('Profiles containing ',datas$datas.FA$n.fats,' Fatty Acids','\n')
}

#' @S3method print Stable_Isotopes
print.Stable_Isotopes <- function(x, ...){
  cat('Stable Isotope data with:','\n')
  cat(datas$n.preds,' Predators','\n')
  cat(datas$n.preys,' Prey species','\n')
  cat(datas$datas.SI$isos,' Stable Isotopes','\n')
}


#' @S3method print Combined_Markers
print.Combined_Markers <- function(x, ...){
  cat('Combined marker data with:','\n')
  cat(datas$n.preds,' Predators','\n')
  cat(datas$n.preys,' Prey species','\n')
  cat('Profiles containing ',datas$datas.FA$n.fats,' Fatty Acids','\n')
  cat('and ',datas$datas.SI$isos,' Stable Isotopes','\n')
  
}