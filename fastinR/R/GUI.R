#' @name fastinR_GUI
#' @title GUI for the fastinR set of functions
#' @description Provides a gui to input data to be analysed with the \code{\link{run_MCMC}} function. The gui rolls the functionality of most of the functions contained within fastinR into a single interface, but some options (plottin options mainly) are not available through the gui.
#' @details Please consult the help individual functions to obtain help on a particular topic
#' @references  Neubauer.P. and Jensen, O.P. (in prep)
#' @author Philipp Neubauer
#' @section Warning: The Tcl/Tk gui interface is very unpredictable, leading to odd errors like internal functions not being found when calling FASTIN functions from the console after having used the GUI. Problems seem to be related to environments. It's usually best to completely restart R when this happens.
#' @seealso \code{\link{simulation}},\code{\link{add_FA}},\code{\link{add_SI}},\code{\link{run_MCMC}}
#' @examples
#' \dontrun{fastin_GUI()}
#' @export
fastinR_GUI <- function(){
  #require(tcltk)   # this is needed - but leads to crashes...
  #require(fgui)
  GUI <<- TRUE
  options(warn = -2)
  # gui helper functions
  
  .fastin <- function(FA.data=NULL,SI.data=NULL,add.covs,varselect,plotdatas,Load.Data=NULL,Save.Data=NULL,MCMC,disp.diags=NULL,Plot.Outputs=NULL,Display.Summaries=NULL,Save.Outputs=NULL){}
  
  dispsummaries <- function(Display.Summary=NULL){output <- guiGetSafe('MCMCout') ; 
                                                  if(class(output)=='pop_props'|class(output)=='ind_props'|class(output)=='cov_props')
                                                  {summary(output)}
  }
  
  dispdiags <- function(accuracy=0.01,proba=0.95,quant=0.025){
    output <- guiGetSafe('MCMCout') 
    if(class(output)=='pop_props'|class(output)=='ind_props'|class(output)=='cov_props')
    {
      diags(MCMCout=output,accuracy,proba,quant)
    }
  }
  
  plotoutputs<- function(){output <- guiGetSafe('MCMCout') 
                           if(class(output)=='pop_props'|class(output)=='ind_props'|class(output)=='cov_props')
                           {plot(output)}
  }
  
  saveoutputs <- function(Path="MCMCout.RData"){output <- guiGetSafe('MCMCout');save(output,file=Path)}
  
  SaveData <- function(Path="datas.RData"){datas <- guiGetSafe('datas');stopifnot(length(datas)>1);save(datas,file=Path)}
  LoadData <- function(Path=NULL){ load(Path); guiSet('datas',datas)}
  
  guiSet( "LIST_WIDTH", 50)
  guiSet( "ENTRY_WIDTH", 10)
  
  #Gui - tried some meaningful indentation here, but still not quite right...
  output <- gui(.fastin, title = 'FASTIN main menu',
                argCommand=list(add.covs=guiNestedF(addCovs,"add.covs",  
                                                    argFilter=list(Groups="{{} {.csv}}",Covariates="{{} {.csv}}"),
                                                    argText=c(Covariates = "Add Covariates (optional)",Groups = "Add Groups (optional)"),
                                                    cancelButton=F,
                                                    exec='Add'),
                                Save.Outputs=guiNestedF(saveoutputs,"Save.Outputs",argText = list(Path='Choose filename'),cancelButton=F,exec='save'),
                                Save.Data=guiNestedF(SaveData,"Save.Data",argText = list(Path='Choose filename'),cancelButton=F,exec='save'),
                                Load.Data=guiNestedF(LoadData,"Load.Data",argFilter=list(Path= "{{} {.RData}}"),cancelButton=F,exec='load'),
                                Display.Summaries=dispsummaries,
                                plotdatas = dataplot,
                                varselect = select_vars,
                                Plot.Outputs=plotoutputs,
                                SI.data = guiNestedF(add_SI,"SI.data",
                                                     argFilter=list(SI.predators= "{{} {.csv}}",SI.preys= "{{} {.csv}}",Frac.Coeffs.mean="{{} {.csv}}",Frac.Coeffs.var="{{} {.csv}}"), 
                                                     title = 'Stable Isotope data entry form',
                                                     exec = "Add Stable Isotope data", closeOnExec = TRUE,
                                                     output = NULL,
                                                     argType=list(datas='i'),
                                                     argText=c(SI.predators="Load predator(s) stable isotope data (csv)",
                                                               SI.preys="Load prey stable isotope data (csv)",
                                                               Frac.Coeffs.mean="Load prey and stable isotope specific fractionation means (csv)",
                                                               Frac.Coeffs.var="Load prey and stable isotope specific fractionation variance (csv)",
                                                               FC.mean="Mean fractionation coefficients (use R's c() notation for multiple isoptopes)",
                                                               FC.var="Variance of fractionation coefficients (use R's c() notation for multiple isoptopes)",
                                                               R.diag.SI = "Diagonal of the prior for predator (co)-variance matrix (>0, smaller value is less informative)"),
                                                     cancelButton=F), 
                                FA.data = guiNestedF(add_FA,"FA.data",
                                                     argFilter=list(FA.predators= "{{} {.csv}}",FA.preys= "{{} {.csv}}",fat.conts = "{{} {.csv}}",Conv.Coeffs.mean="{{} {.csv}}",Conv.Coeffs.var="{{} {.csv}}"), 
                                                     argType=list(datas='i'),
                                                     #argEdit = list(CC.mean=NULL,CC.var=NULL,R.diag=NULL),
                                                     title = 'Fatty Acid Profile data entry form',
                                                     argText=c(FA.predators="Load predator(s) fatty acid data (csv)",
                                                               FA.preys="Load prey fatty acid data (csv)",
                                                               Conv.Coeffs.mean="Load prey and fatty acid specific conversion means (csv)",
                                                               Conv.Coeffs.var="Load prey and fatty acid specific conversion variance (csv)",
                                                               fat.conts = "Prey fat content (csv)",
                                                               FC.mean="Mean fat content (use R's c() notation)",
                                                               FC.var="Variance of fat content (use R's c() notation)", 
                                                               CC.mean="Mean conversion coefficients (use R's c() notation for multiple FAs)",
                                                               CC.var="SD of conversion coefficients (use R's c() notation for multiple FAs)",
                                                               R.diag = "Diagonal of the prior for predator (co)-variance matrix (>0, smaller value is less informative)"),
                                                     exec = "Add Fatty Acid data", cancelButton=F,closeOnExec = TRUE,output = NULL),  
                                MCMC = guiNestedF(run_MCMC,"MCMC",
                                                  title = 'Run FASTIN analysis',
                                                  exec = "Run MCMC", closeOnExec = TRUE,
                                                  output = NULL,
                                                  
                                                  argOption = list(Data.Type = c('Fatty.Acid.Profiles','Stable.Isotopes','Combined.Analysis'), defaultChoice=1,spawn=c(FALSE,TRUE),defaultChoice=1,Analysis.Type = c('Population.proportions','Individual.proportions','Analysis.with.Covariates'), defaultChoice=1),    
                                                  argType=list(datas='i',plott='i',Covs='i'),
                                                  argText=list(nIter = 'Number of MCMC iterations',
                                                               nBurnin = 'Number of MCMC iterations to discard (burn-in)',
                                                               nChains = 'Number of Markov Chains',
                                                               nThin = 'Thinning interval of MCMC chains',
                                                               even= 'Prior eveness of proportions',
                                                               Rnot= 'Prior diagonal predator FA covariance matrix',
                                                               Rnot.SI= 'Prior diagonal predator SI covariance matrix'
                                                  ),cancelButton=F),
                                disp.diags = guiNestedF(dispdiags,"disp.diags",
                                                        title = 'Show Convergence Diagnositcs', 
                                                        exec = "Display diagnostics in the R console",
                                                        argText=c(proba = 'probability',
                                                                  quant='quantile to estiamte'
                                                        ),cancelButton=F)
                ),                                 
                argText=c(SI.data='Add Stable Isotope data', 
                          FA.data='Add Fatty Acid profiles',
                          plotdatas = 'Plot data',
                          varselect= 'Select a subset of Fatty Acids',
                          Save.Data= "Save data",
                          Load.Data = 'Load previously saved dataset',
                          add.covs='Add covariates and/or groups',
                          MCMC = 'Run Bayesian analysis (MCMC)',
                          disp.diags = 'Show convergence diagnostics'
                ),
                exec=NULL,output=NULL,argGridOrder=c(1,1,1,2,2,3,3,4,4,4,5,5), argGridSticky=rep("w",length(formals(.fastin)))
  )
  GUI <<- FALSE
  return(output)
  
}