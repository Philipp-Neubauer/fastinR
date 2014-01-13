#' Fatty Acids and Stable Isotopes in R
#'
#' This packages provides standalone functions and gui facilities to perfom Bayesian estiamtion of diet proportions based on fatty acid profiles, stable isotopes and their combination. Bayesian estiamtion is currently performed in JAGS using the R package rjags. The package can handle covariates and grouped (e.g., regional) data as predictors of diet proportions, and will estimate both individual and popualtion level diet proportions. Functionality for simulating, plotting and extracting/writing summary statistics to file is provided as well.
#' @author Philipp Neubauer
#' @references  Neubauer.P. and Jensen, O.P. (in prep)
#' @importFrom MCMCpack rdirichlet
#' @importFrom reshape melt.data.frame
#' @import vegan fgui tcltk rjags lattice MASS ggplot2 grid
#' @docType package
#' @name FASTIN-package
NULL