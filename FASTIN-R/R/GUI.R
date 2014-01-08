.fastin <- function(SI.data=NULL,FA.data=NULL,groupings,Load.Data=NULL,Save.Data=NULL,add.covs,MCMC,Save.Outputs=NULL,Display.Summaries=NULL,Plot.Outputs=NULL,disp.diags=NULL){}

resetsc <- function(datas=NULL){if (is.null(datas)){datas <- guiGetSafe('datas')};if(length(datas)>1){datas$SC=F;guiSet('datas',datas);return(datas)}}

FASTIN <- function(){
  #require(tcltk)   # this is needed - but leads to crashes...
  #require(fgui)
  GUI <- TRUE
  # gui helper functions
    
    pnorm_even <- function(even=0.1){p=2*(1-pnorm(log(95)/2,0,sqrt(1/even)));return(p)}
  
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
  
    saveoutputs <- function(Path="MCMCout.Rdata"){output <- guiGetSafe('MCMCout');save(output,file=Path)}
  
    SaveData <- function(Path="datas.Rdata"){datas <- guiGetSafe('datas');stopifnot(length(datas)>1);save(datas,file=Path)}
    LoadData <- function(Path=NULL){ load(Path); guiSet('datas',datas)}
  
    guiSet( "LIST_WIDTH", 50)
    guiSet( "ENTRY_WIDTH", 10)
  
  #Gui - tried some meaningful indentation here, but still not quite right...
    output <- gui(.fastin, title = 'FASTIN main menu',
                 argCommand=list(add.covs=guiNestedF(addcovs,"add.covs",  argFilter=list(Groups="{{} {.csv}}",Covariates="{{} {.csv}}"),
                                                     argText=c(Covariates = "Add Covariates (optional)",Groups = "Add Groups (optional)"),cancelButton=F,exec='Add'),
                     Save.Outputs=guiNestedF(saveoutputs,"Save.Outputs",argText = list(Path='Choose filename'),cancelButton=F,exec='save'),
                     Save.Data=guiNestedF(SaveData,"Save.Data",argText = list(Path='Choose filename'),cancelButton=F,exec='save'),
                     Load.Data=guiNestedF(LoadData,"Load.Data",argFilter=list(Path= "{{} {.Rdata}}"),cancelButton=F,exec='load'),
                     Display.Summaries=dispsummaries,
                     groupings = resetsc,
                     Plot.Outputs=plotoutputs,
                     SI.data = guiNestedF(addSI,"SI.data",
                                                      argFilter=list(predators.SI= "{{} {.csv}}",preys.SI= "{{} {.csv}}",Frac.Coeffs.mean="{{} {.csv}}",Frac.Coeffs.var="{{} {.csv}}"), 
                                                      title = 'Stable Isotope data entry form',
                                                      exec = "Add Stable Isotope data", closeOnExec = TRUE,
                                                      output = NULL,
                                                      argText=c(predators.SI="Load predator(s) stable isotope data (csv)",
                                                                preys.SI="Load prey stable isotope data (csv)",
                                                                Frac.Coeffs.mean="Load prey and stable isotope specific fractionation means (csv)",
                                                                Frac.Coeffs.var="Load prey and stable isotope specific fractionation sd (csv)",
                                                                FC.mean="Mean fractionation coefficients (use R's c() notation)",
                                                                FC.var="SD of fractionation coefficients (use R's c() notation)",
                                                                R.diag.SI = "Diagonal of the prior for predator (co)-variance matrix (>0, smaller value is less informative)"),
                                                      cancelButton=F), 
                     FA.data = guiNestedF(addFA,"FA.data",argFilter=list(predators.FA= "{{} {.csv}}",preys.FA= "{{} {.csv}}",fat.conts = "{{} {.csv}}",Conv.Coeffs.mean="{{} {.csv}}",Conv.Coeffs.var="{{} {.csv}}"), 
                                                      #argEdit = list(CC.mean=NULL,CC.var=NULL,R.diag=NULL),
                                                      title = 'Fatty Acid Profile data entry form',
                                                      argText=c(predators.FA="Load predator(s) fatty acid data (csv)",
                                                                preys.FA="Load prey fatty acid data (csv)",
                                                                Conv.Coeffs.mean="Load prey and fatty acid specific conversion means (csv)",
                                                                Conv.Coeffs.var="Load prey and fatty acid specific conversion sd (csv)",
                                                                fat.conts = "Prey fat content (csv)",
                                                                FC.mean="Mean fat content (use R's c() notation)",
                                                                FC.var="SD of fat content (use R's c() notation)", 
                                                                CC.mean="Mean conversion coefficients (use R's c() notation)",
                                                                CC.var="SD of conversion coefficients (use R's c() notation)",
                                                                R.diag = "Diagonal of the prior for predator (co)-variance matrix (>0, smaller value is less informative)"),
                                                      exec = "Add Fatty Acid data", cancelButton=F,closeOnExec = TRUE,output = NULL),  
                     MCMC = guiNestedF(run_MCMC,"MCMC",
                                                   title = 'Run FASTIN analysis',
                                                   exec = "Run MCMC", closeOnExec = TRUE,
                                                   output = NULL,
                                                   argCommand=list(even = guiNestedF(pnorm_even,'even',exec='Set new prior',argSlider=list(even=c(0.001,2,0.001),output=c(0,1,0.001)),callback=guiExec,argText=c(output = "Probability that prop(prey x) = 0.95*prop(preys other than x)"))),
                                                   argOption = list(Data.Type = c('Fatty.Acid.Profiles','Stable.Isotopes','Combined.Analysis'), defaultChoice=1,Analysis.Type = c('Population.proportions','Individual.proportions','Analysis.with.Covariates'), defaultChoice=1),                                                   
                                                   argText=c(nIter = 'Number of MCMC iterations',
                                                             nBurnin = 'Number of MCMC iterations to discard (burn-in)',
                                                             nChains = 'Number of Markov Chains',
                                                             nThin = 'Thinning interval of MCMC chains',
                                                             even= 'Prior eveness of proportions'
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
                           groupings = 'Reset previous prey grouping',
                           Save.Data= "Save data",
                           Load.Data = 'Load previously saved dataset',
                           add.covs='Add covariates and/or groups',
                           MCMC = 'Run Bayesian analysis (MCMC)',
                     disp.diags = 'Show convergence diagnostics'
                 ),
                 exec=NULL,output=NULL,argGridOrder=c(1,1,1,2,2,2,3,3,4,4,4), argGridSticky=rep("w",length(formals(.fastin)))
  )
  
  return(output)
  
}
