<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{fastinR: Fatty Acids and Stable Isotopes in R}
-->

Estimating diet proportions form Fatty Acids and Stable Isotopes
========================================================

DISCLAIMER: This is an evolving package and vignette, please read instructions [here](https://github.com/Philipp-Neubauer/fastinR/blob/master/README.md) for installation and to see what dependencies are required. If you find a bug or want so suggest improvements, please [submit an issue on github](https://github.com/Philipp-Neubauer/fastinR/issues?state=open), collaborations and contributions are very welcome!

We start with a very simple simulated example. The easiest way to simulate relevant data is through the built in simulation GUI, which is called from the command line once fastinR has been loaded:




```r
library(fastinR, warn.conflicts = F)
```

```
## Loading required package: rjags
## Loading required package: coda
## Loading required package: lattice
## Linked to JAGS 3.4.0
## Loaded modules: basemod,bugs
```

```
simulation()
```
The simulation window displays, and the user is invited to select parameters for simulated data, which can be plotted to visualize source overlap etc. For more details call `?simulation`.

The package comes with a simulated dataset which can be called with , which loads an object called _datas_ into the current environment.

Loading data from files
-----------------------

Raw data should be stored in .csv files in a prescribed format, see `?add_FA`, `?add_SI` and `?add_Covs` for details on formatting. You can also inspect saved output from `simulation()` to get a better idea of the correct file format.

Files are imported with the help of three constructor functions: `add_FA` to import fatty acid data, `add_SI` to import Stable Isotope data and `add_Covs` to import data on covariates or groupings that may be influencing predator diets. 

The package comes with simulated data in .csv files, which can be found using the `system.file` function. The `add_SI` constructor takes separate files for predator and prey Stable Isotope data, as well as files for additive fractionation coefficient means and variances - these are optional and can be specified manually in the function call, see the function help for details.


```r
SI.predators <- system.file("extdata", "Simdata_SI_preds.csv", package = "fastinR")
SI.preys <- system.file("extdata", "Simdata_SI_preys.csv", package = "fastinR")
Frac.Coeffs.mean <- system.file("extdata", "Simdata_SI_fc_means.csv", package = "fastinR")
Frac.Coeffs.var <- system.file("extdata", "Simdata_SI_fc_var.csv", package = "fastinR")

dats <- add_SI(SI.predators, SI.preys, Frac.Coeffs.mean, Frac.Coeffs.var)
```


We can now visualize the data on an Non-Metric Dimensional Scaling plot like so:


```r
dataplot(dats)
```

```
## Run 0 stress 0 
## Run 1 stress 0.4121 
## Run 2 stress 8.945e-05 
## ... procrustes: rmse 8.395e-05  max resid 0.0004449 
## *** Solution reached
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 


The dats object now has a set of data in the right format for further analysis. We'll add fatty acids before we proceed. As with `add_SI`, the `add_FA` constructor takes separate files for predator and prey fatty acid data, conversion coefficient means and variances as well as fat content - these are again optional and can be specified manually in the function call, see the function help for details.


```r
FA.predators <- system.file("extdata", "Simdata_FA_preds.csv", package = "fastinR")
FA.preys <- system.file("extdata", "Simdata_FA_preys.csv", package = "fastinR")
Conv.Coeffs.mean <- system.file("extdata", "Simdata_FA_cc_means.csv", package = "fastinR")
Conv.Coeffs.var <- system.file("extdata", "Simdata_FA_cc_var.csv", package = "fastinR")
fat.conts <- system.file("extdata", "Simdata_fat_cont.csv", package = "fastinR")

dats <- add_FA(FA.predators, FA.preys, fat.conts, Conv.Coeffs.mean, Conv.Coeffs.var, 
    datas = dats)
```



*Note* the last argument in `add_FA` now contains a reference to the data object that was constructed from the Stable Isotopes beforehand. `add_FA` just adds fatty acid data to the mix - the same would apply if data were added the other way around (only `add_Covs` works separately, it needs it's own object).

Plotting the joint dataset:


```r
dataplot(dats)
```

```
## Run 0 stress 0.06002 
## Run 1 stress 0.06002 
## ... procrustes: rmse 3.045e-05  max resid 0.0002547 
## *** Solution reached
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 


Data Grooming
-------------

The grooming step in this case corresponds to selecting relevant variables for the analysis. In practice, one could measure lots of Fatty Acids, and researchers often choose arbitrary cutoff points in proportions to reduce the dataset. However, even Fatty Acids occurring in low proportions may be useful for source discrimination and diet estimation if they introduce systematic differences between sources (prey items). 
The `select_vars` function provides a graphical way to choose a subset of Fatty Acids according to their contribution to prey separation and reduction in collinearity in the prey matrix. Using `select_vars` in the simulated example shows that only really 2 Fatty Acids contribute to source separation - the plot of cumulative separation converges to one after adding the two most important fatty acid proportions (see also cumulative separation scores displayed in the console): all remaining FAs only add to collinearity in the source matrix. We choose F Atty Acid no 3,5 and 7 for analysis (the second call does this directly without the graphical interface).

```
dats.subset <- select_vars(dats)
```
or

```r
dats.subset <- select_vars(dats, c(3, 2, 6))

# inspecting
dats.subset$datas.FA$n.fats
```

```
## [1] 3
```


The data object `dats` is now ready for analysis, we can plot the dataset with it's new subset of fatty acids as before:


```r
dataplot(dats.subset)
```

```
## Run 0 stress 0.02032 
## Run 1 stress 0.02032 
## ... procrustes: rmse 0.001026  max resid 0.009586 
## *** Solution reached
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8.png) 


Apart from a rotation, the overall configuration should be similar.

Bayesian Analysis
-----------------

The actual analysis can be performed locally (in the active R session) or in a distributed way, using as many R sessions as Markov Chains. A good strategy is usually to run a few chains locally for short runs, and then run 3 or more chains in parallel in a distributed way once a satisfactory set of priors, thinning interval etc has been found. The `run_MCMC` functions sets up the MCMC runs and takes `spawn = T` or `spawn = F` as parameter to run chains in R slave processes or locally, respectively. 

For a combined analysis of stable isotopes and fatty acids, it is often useful to run the two datatypes separately to assure that good priors can be found for both datasets independently, and then combining them for a final analysis. The analysis type is chosen in the appropriate option in the function call.

### Estiamting population proportions

#### From Stable Isotopes alone

Lets start with an analysis of the stable isotopes, estimating only global (population) level diets. We will use the default prior on the predator covariance matrix, and will adjust this prior subsequently. *WARNING* This might take a while depending on your resources, the size of the dataset and the parameters used for the MCMC.


```r
Pop.SI <- run_MCMC(datas = dats.subset, nIter = 10000, nBurnin = 1000, nChains = 3, 
    nThin = 10, Data.Type = "Stable.Isotopes", Analysis.Type = "Population.proportions", 
    Rnot_SI = 1, plott = F, spawn = F)
```

```
## Compiling model graph
##    Resolving undeclared variables
##    Allocating nodes
##    Graph Size: 96
## 
## Initializing model
## 
## 
##  proceeding to burn-in phase 
## 
##  sampling from parameter distributions
```


Plotting the MCMC is the easiest way to ensure that the sampler is mixing - meaning that the chain explores the posterior distribution of each parameter efficiently.


```r
MCMCplot(Pop.SI)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10.png) 


The mixing isn't great, `diags` function gives more information:


```r
diags(Pop.SI)
```

```
## 
##  
## ################################# 
## Raftery-Lewis diagnostics 
## ################################# 
##  
## [[1]]
## 
## Quantile (q) = 0.025
## Accuracy (r) = +/- 0.01
## Probability (s) = 0.95 
##                                                
##          Burn-in  Total Lower bound  Dependence
##          (M)      (N)   (Nmin)       factor (I)
##  prop[1] 30       10530 937          11.2      
##  prop[2] 30       10530 937          11.2      
##  prop[3] 20       9690  937          10.3      
## 
## 
## [[2]]
## 
## Quantile (q) = 0.025
## Accuracy (r) = +/- 0.01
## Probability (s) = 0.95 
##                                                
##          Burn-in  Total Lower bound  Dependence
##          (M)      (N)   (Nmin)       factor (I)
##  prop[1] 30       10530 937          11.2      
##  prop[2] 30       10530 937          11.2      
##  prop[3] 30       11430 937          12.2      
## 
## 
## [[3]]
## 
## Quantile (q) = 0.025
## Accuracy (r) = +/- 0.01
## Probability (s) = 0.95 
##                                                
##          Burn-in  Total Lower bound  Dependence
##          (M)      (N)   (Nmin)       factor (I)
##  prop[1] 30       11430 937          12.2      
##  prop[2] 30       11430 937          12.2      
##  prop[3] 40       11610 937          12.4      
## 
## 
## 
##  Based on these diagnostics you should repeat the MCMC with  11610  iterations 
##  and a thinning interval of  12  ,if these values are higher than the values 
##  used to produce these diagnostics 
##  
## 
##  
## ################################# 
## Gelman-Rubin diagnostics 
## ################################# 
##  
## Potential scale reduction factors:
## 
##         Point est. Upper C.I.
## prop[1]          1       1.00
## prop[2]          1       1.01
## prop[3]          1       1.02
## 
## Multivariate psrf
## 
## 1.01
## 
##  Both univariate upper C.I. and multivariate psrf 
##  should be close to 1 if the chains converged 
##  
## 
```


Based on the diagnostics, it seems that we're doing OK, but that the Stable Isotopes don't give much certainty about diet proportions:


```r
summary(Pop.SI)
```

```
## 
##  
##  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
## Printing results for MCMC chain 1
##  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
## 
##  
##  
##  population diet proportions 
##  
##  
##        Prey_1  Prey_2   Prey_3
## 2.5%  0.02529 0.01992 0.009298
## 50%   0.44921 0.40645 0.118759
## 97.5% 0.93022 0.86187 0.385118
## 
##  
##  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
## Printing results for MCMC chain 2
##  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
## 
##  
##  
##  population diet proportions 
##  
##  
##        Prey_1  Prey_2  Prey_3
## 2.5%  0.02319 0.02022 0.01166
## 50%   0.43297 0.41608 0.11863
## 97.5% 0.91325 0.84226 0.40128
## 
##  
##  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
## Printing results for MCMC chain 3
##  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
## 
##  
##  
##  population diet proportions 
##  
##  
##       Prey_1  Prey_2  Prey_3
## 2.5%  0.0241 0.01951 0.01127
## 50%   0.4443 0.40226 0.13686
## 97.5% 0.9269 0.82588 0.41910
```


```r
plot(Pop.SI, save = F)
```

```
## Using  as id variables
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-131.png) ![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-132.png) 


The credible intervals are very large for Prey items 1 and 2, and the correlation plot suggests that there is a reasonably strong posterior correlation between estimates of prey items 1 & 2.

Increasing the default prior on the predator covariance matrix sometimes leads to better mixing. We'll also run more iterations and set `spawn=T` for this.


```r
Pop.SI2 <- run_MCMC(datas = dats.subset, nIter = 30000, nBurnin = 1000, nChains = 3, 
    nThin = 30, Data.Type = "Stable.Isotopes", Analysis.Type = "Population.proportions", 
    Rnot_SI = 5, plott = F, spawn = T)
```

```
## +++++
```



```r
MCMCplot(Pop.SI2)
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15.png) 



```r
diags(Pop.SI2)
```

```
## 
##  
## ################################# 
## Raftery-Lewis diagnostics 
## ################################# 
##  
## [[1]]
## 
## Quantile (q) = 0.025
## Accuracy (r) = +/- 0.01
## Probability (s) = 0.95 
##                                                
##          Burn-in  Total Lower bound  Dependence
##          (M)      (N)   (Nmin)       factor (I)
##  prop[1] 60       26790 937          28.6      
##  prop[2] 60       26790 937          28.6      
##  prop[3] 60       29070 937          31.0      
## 
## 
## [[2]]
## 
## Quantile (q) = 0.025
## Accuracy (r) = +/- 0.01
## Probability (s) = 0.95 
##                                                
##          Burn-in  Total Lower bound  Dependence
##          (M)      (N)   (Nmin)       factor (I)
##  prop[1] 60       29070 937          31.0      
##  prop[2] 60       29070 937          31.0      
##  prop[3] 60       26790 937          28.6      
## 
## 
## [[3]]
## 
## Quantile (q) = 0.025
## Accuracy (r) = +/- 0.01
## Probability (s) = 0.95 
##                                                
##          Burn-in  Total Lower bound  Dependence
##          (M)      (N)   (Nmin)       factor (I)
##  prop[1] 60       29070 937          31.0      
##  prop[2] 60       26790 937          28.6      
##  prop[3] 60       29190 937          31.2      
## 
## 
## 
##  Based on these diagnostics you should repeat the MCMC with  29190  iterations 
##  and a thinning interval of  31  ,if these values are higher than the values 
##  used to produce these diagnostics 
##  
## 
##  
## ################################# 
## Gelman-Rubin diagnostics 
## ################################# 
##  
## Potential scale reduction factors:
## 
##         Point est. Upper C.I.
## prop[1]          1       1.00
## prop[2]          1       1.01
## prop[3]          1       1.01
## 
## Multivariate psrf
## 
## 1
## 
##  Both univariate upper C.I. and multivariate psrf 
##  should be close to 1 if the chains converged 
##  
## 
```


Mixing doesn't seem to get much better, indicating that this is probably as good as it gets with stable isotopes alone.

#### Fatty Acids alone

Lets repeat this with Fatty Acids, again starting with the default prior. Note that `spawn=T` this time to save some time.


```r
Pop.FA <- run_MCMC(datas = dats.subset, nIter = 10000, nBurnin = 1000, nChains = 3, 
    nThin = 10, Data.Type = "Fatty.Acid.Profiles", Analysis.Type = "Population.proportions", 
    Rnot = 0.2, plott = F, spawn = T)
```

```
## ++++++++
```



```r
MCMCplot(Pop.FA)
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-18.png) 


Once again the mixing isn't great:


```r
diags(Pop.FA)
```

```
## 
##  
## ################################# 
## Raftery-Lewis diagnostics 
## ################################# 
##  
## [[1]]
## 
## Quantile (q) = 0.025
## Accuracy (r) = +/- 0.01
## Probability (s) = 0.95 
##                                                
##          Burn-in  Total Lower bound  Dependence
##          (M)      (N)   (Nmin)       factor (I)
##  prop[1] 100      25880 937          27.6      
##  prop[2] 50       13530 937          14.4      
##  prop[3] 80       21170 937          22.6      
## 
## 
## [[2]]
## 
## Quantile (q) = 0.025
## Accuracy (r) = +/- 0.01
## Probability (s) = 0.95 
##                                                
##          Burn-in  Total Lower bound  Dependence
##          (M)      (N)   (Nmin)       factor (I)
##  prop[1] 140      36460 937          38.9      
##  prop[2] 100      20860 937          22.3      
##  prop[3] 90       23360 937          24.9      
## 
## 
## [[3]]
## 
## Quantile (q) = 0.025
## Accuracy (r) = +/- 0.01
## Probability (s) = 0.95 
##                                                
##          Burn-in  Total Lower bound  Dependence
##          (M)      (N)   (Nmin)       factor (I)
##  prop[1] 140      50160 937          53.5      
##  prop[2] 60       16080 937          17.2      
##  prop[3] 100      25880 937          27.6      
## 
## 
## 
##  Based on these diagnostics you should repeat the MCMC with  50160  iterations 
##  and a thinning interval of  54  ,if these values are higher than the values 
##  used to produce these diagnostics 
##  
## 
##  
## ################################# 
## Gelman-Rubin diagnostics 
## ################################# 
##  
## Potential scale reduction factors:
## 
##         Point est. Upper C.I.
## prop[1]       1.01       1.04
## prop[2]       1.01       1.05
## prop[3]       1.02       1.02
## 
## Multivariate psrf
## 
## 1.02
## 
##  Both univariate upper C.I. and multivariate psrf 
##  should be close to 1 if the chains converged 
##  
## 
```


Again not too bad according to the diagnostics, although the chains should be run for longer (see Raftery-Lewis diagnostics ) - but lets try again with `Rnot = 1`


```r
Pop.FA2 <- run_MCMC(datas = dats.subset, nIter = 10000, nBurnin = 1000, nChains = 3, 
    nThin = 10, Data.Type = "Fatty.Acid.Profiles", Analysis.Type = "Population.proportions", 
    Rnot = 1, plott = F, spawn = T)
```

```
## ++++++++
```



```r
MCMCplot(Pop.FA2)
```

![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-21.png) 



```r
diags(Pop.FA2)
```

```
## 
##  
## ################################# 
## Raftery-Lewis diagnostics 
## ################################# 
##  
## [[1]]
## 
## Quantile (q) = 0.025
## Accuracy (r) = +/- 0.01
## Probability (s) = 0.95 
##                                                
##          Burn-in  Total Lower bound  Dependence
##          (M)      (N)   (Nmin)       factor (I)
##  prop[1] 110      28820 937          30.8      
##  prop[2] 30       10530 937          11.2      
##  prop[3] 40       12430 937          13.3      
## 
## 
## [[2]]
## 
## Quantile (q) = 0.025
## Accuracy (r) = +/- 0.01
## Probability (s) = 0.95 
##                                                
##          Burn-in  Total Lower bound  Dependence
##          (M)      (N)   (Nmin)       factor (I)
##  prop[1] 40       12430 937          13.3      
##  prop[2] 30       10530 937          11.2      
##  prop[3] 70       19260 937          20.6      
## 
## 
## [[3]]
## 
## Quantile (q) = 0.025
## Accuracy (r) = +/- 0.01
## Probability (s) = 0.95 
##                                                
##          Burn-in  Total Lower bound  Dependence
##          (M)      (N)   (Nmin)       factor (I)
##  prop[1] 80       21170 937          22.6      
##  prop[2] 30       10530 937          11.2      
##  prop[3] 50       13530 937          14.4      
## 
## 
## 
##  Based on these diagnostics you should repeat the MCMC with  28820  iterations 
##  and a thinning interval of  31  ,if these values are higher than the values 
##  used to produce these diagnostics 
##  
## 
##  
## ################################# 
## Gelman-Rubin diagnostics 
## ################################# 
##  
## Potential scale reduction factors:
## 
##         Point est. Upper C.I.
## prop[1]       1.01       1.03
## prop[2]       1.00       1.01
## prop[3]       1.02       1.05
## 
## Multivariate psrf
## 
## 1.02
## 
##  Both univariate upper C.I. and multivariate psrf 
##  should be close to 1 if the chains converged 
##  
## 
```


Looks very similar, so lets do a final run with 30k iterations as suggested by the `diags` output and a thinning interval of 30 should be better. (Note that values suggested by the diagnostics may vary from one run to the next. Also note that we already ran 30k as we ran 3 parallel chains here...).


```r
Pop.FA3 <- run_MCMC(datas = dats.subset, nIter = 30000, nBurnin = 1000, nChains = 3, 
    nThin = 30, Data.Type = "Fatty.Acid.Profiles", Analysis.Type = "Population.proportions", 
    Rnot = 1, plott = F, spawn = T)
```

```
## ++++++++++++++++++++
```



```r
MCMCplot(Pop.FA3)
```

![plot of chunk unnamed-chunk-24](figure/unnamed-chunk-24.png) 



```r
diags(Pop.FA3)
```

```
## 
##  
## ################################# 
## Raftery-Lewis diagnostics 
## ################################# 
##  
## [[1]]
## 
## Quantile (q) = 0.025
## Accuracy (r) = +/- 0.01
## Probability (s) = 0.95 
##                                                
##          Burn-in  Total Lower bound  Dependence
##          (M)      (N)   (Nmin)       factor (I)
##  prop[1] 90       34290 937          36.6      
##  prop[2] 60       29070 937          31.0      
##  prop[3] 60       26790 937          28.6      
## 
## 
## [[2]]
## 
## Quantile (q) = 0.025
## Accuracy (r) = +/- 0.01
## Probability (s) = 0.95 
##                                                
##          Burn-in  Total Lower bound  Dependence
##          (M)      (N)   (Nmin)       factor (I)
##  prop[1] 90       31590 937          33.7      
##  prop[2] 60       29070 937          31.0      
##  prop[3] 90       31590 937          33.7      
## 
## 
## [[3]]
## 
## Quantile (q) = 0.025
## Accuracy (r) = +/- 0.01
## Probability (s) = 0.95 
##                                                
##          Burn-in  Total Lower bound  Dependence
##          (M)      (N)   (Nmin)       factor (I)
##  prop[1] 120      37290 937          39.8      
##  prop[2] 90       31590 937          33.7      
##  prop[3] 60       29070 937          31.0      
## 
## 
## 
##  Based on these diagnostics you should repeat the MCMC with  37290  iterations 
##  and a thinning interval of  40  ,if these values are higher than the values 
##  used to produce these diagnostics 
##  
## 
##  
## ################################# 
## Gelman-Rubin diagnostics 
## ################################# 
##  
## Potential scale reduction factors:
## 
##         Point est. Upper C.I.
## prop[1]          1       1.00
## prop[2]          1       1.01
## prop[3]          1       1.00
## 
## Multivariate psrf
## 
## 1
## 
##  Both univariate upper C.I. and multivariate psrf 
##  should be close to 1 if the chains converged 
##  
## 
```

```r
summary(Pop.FA3)
```

```
## 
##  
##  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
## Printing results for MCMC chain 1
##  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
## 
##  
##  
##  population diet proportions 
##  
##  
##       Prey_1 Prey_2  Prey_3
## 2.5%  0.1503 0.0173 0.01181
## 50%   0.6302 0.2145 0.08448
## 97.5% 0.9299 0.7408 0.49537
## 
##  
##  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
## Printing results for MCMC chain 2
##  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
## 
##  
##  
##  population diet proportions 
##  
##  
##       Prey_1  Prey_2  Prey_3
## 2.5%  0.1543 0.02095 0.01150
## 50%   0.6294 0.21638 0.09061
## 97.5% 0.9285 0.76039 0.51392
## 
##  
##  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
## Printing results for MCMC chain 3
##  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
## 
##  
##  
##  population diet proportions 
##  
##  
##       Prey_1  Prey_2  Prey_3
## 2.5%  0.1639 0.01747 0.01271
## 50%   0.6251 0.20583 0.10026
## 97.5% 0.9331 0.69805 0.52857
```


```r
plot(Pop.FA3, save = F)
```

```
## Using  as id variables
```

![plot of chunk unnamed-chunk-26](figure/unnamed-chunk-261.png) ![plot of chunk unnamed-chunk-26](figure/unnamed-chunk-262.png) 

There is still substantial uncertainty about the diet proportions: it looks as though Prey 1 is a dominant source, but the credible intervals are large. Let's try combining stable isotopes and fatty acids:

#### Combining fatty acids and stable isotopes


```r
Pop.Combined <- run_MCMC(datas = dats.subset, nIter = 30000, nBurnin = 1000, 
    nChains = 3, nThin = 30, Data.Type = "Combined.Analysis", Analysis.Type = "Population.proportions", 
    Rnot = 1, Rnot_SI = 1, plott = F, spawn = T)
```

```
## ++++++++++++++++++++++++
```



```r
MCMCplot(Pop.Combined)
```

![plot of chunk unnamed-chunk-28](figure/unnamed-chunk-28.png) 


```r
diags(Pop.Combined)
```

```
## 
##  
## ################################# 
## Raftery-Lewis diagnostics 
## ################################# 
##  
## [[1]]
## 
## Quantile (q) = 0.025
## Accuracy (r) = +/- 0.01
## Probability (s) = 0.95 
##                                                
##          Burn-in  Total Lower bound  Dependence
##          (M)      (N)   (Nmin)       factor (I)
##  prop[1] 90       31590 937          33.7      
##  prop[2] 90       34290 937          36.6      
##  prop[3] 60       29070 937          31.0      
## 
## 
## [[2]]
## 
## Quantile (q) = 0.025
## Accuracy (r) = +/- 0.01
## Probability (s) = 0.95 
##                                                
##          Burn-in  Total Lower bound  Dependence
##          (M)      (N)   (Nmin)       factor (I)
##  prop[1] 60       29070 937          31.0      
##  prop[2] 60       29070 937          31.0      
##  prop[3] 150      40590 937          43.3      
## 
## 
## [[3]]
## 
## Quantile (q) = 0.025
## Accuracy (r) = +/- 0.01
## Probability (s) = 0.95 
##                                                
##          Burn-in  Total Lower bound  Dependence
##          (M)      (N)   (Nmin)       factor (I)
##  prop[1] 60       29070 937          31.0      
##  prop[2] 90       34290 937          36.6      
##  prop[3] 150      44190 937          47.2      
## 
## 
## 
##  Based on these diagnostics you should repeat the MCMC with  44190  iterations 
##  and a thinning interval of  47  ,if these values are higher than the values 
##  used to produce these diagnostics 
##  
## 
##  
## ################################# 
## Gelman-Rubin diagnostics 
## ################################# 
##  
## Potential scale reduction factors:
## 
##         Point est. Upper C.I.
## prop[1]       1.00       1.01
## prop[2]       1.01       1.02
## prop[3]       1.00       1.00
## 
## Multivariate psrf
## 
## 1.01
## 
##  Both univariate upper C.I. and multivariate psrf 
##  should be close to 1 if the chains converged 
##  
## 
```

```r
summary(Pop.Combined)
```

```
## 
##  
##  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
## Printing results for MCMC chain 1
##  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
## 
##  
##  
##  population diet proportions 
##  
##  
##       Prey_1 Prey_2  Prey_3
## 2.5%  0.2330 0.0247 0.01435
## 50%   0.5848 0.2979 0.08575
## 97.5% 0.9173 0.6539 0.29892
## 
##  
##  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
## Printing results for MCMC chain 2
##  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
## 
##  
##  
##  population diet proportions 
##  
##  
##       Prey_1  Prey_2  Prey_3
## 2.5%  0.2153 0.03118 0.01344
## 50%   0.5851 0.29470 0.08594
## 97.5% 0.9027 0.67988 0.28856
## 
##  
##  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
## Printing results for MCMC chain 3
##  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
## 
##  
##  
##  population diet proportions 
##  
##  
##       Prey_1  Prey_2  Prey_3
## 2.5%  0.2258 0.03145 0.01568
## 50%   0.6097 0.27467 0.09186
## 97.5% 0.9134 0.67428 0.28766
```


```r
plot(Pop.Combined, save = F)
```

```
## Using  as id variables
```

![plot of chunk unnamed-chunk-30](figure/unnamed-chunk-301.png) ![plot of chunk unnamed-chunk-30](figure/unnamed-chunk-302.png) 


The combined analysis reduces uncertainty slightly, especially for source 3, but sources remain uncertain with somewhat large credible intervals. 

To compare the three approaches explicitly, we can use `multiplot`. For `multiplot` the three results need to be combined in a list:


```r
Pop.list <- list(Stable.Isotopes = Pop.SI2, Fatty.Acids = Pop.FA3, Combined = Pop.Combined)

multiplot(Pop.list, save = F)
```

```
## Using V1 as id variables
```

![plot of chunk unnamed-chunk-31](figure/unnamed-chunk-31.png) 


The original data was simulated to include groups of predators with rather different proportions (visible in the dataplots above). Lets have a look at the proportions that were used to simulate the dataset:


```r
proportions.path <- system.file("extdata", "Simdata_props.csv", package = "fastinR")
proportions <- read.csv(proportions.path, header = F, row.names = 1)

colnames(proportions) <- unique(dats$prey.ix)

proportions
```

```
##    Prey_1  Prey_2  Prey_3
## 1  0.1315 0.78796 0.08052
## 2  0.8168 0.03774 0.14550
## 3  0.1349 0.76823 0.09684
## 4  0.2429 0.64427 0.11285
## 5  0.1494 0.72239 0.12826
## 6  0.7965 0.07278 0.13076
## 7  0.7307 0.08529 0.18398
## 8  0.8189 0.07136 0.10973
## 9  0.8169 0.07488 0.10824
## 10 0.1727 0.74464 0.08268
```

```r

# overall mean proportions
colMeans(proportions)
```

```
## Prey_1 Prey_2 Prey_3 
## 0.4811 0.4010 0.1179
```


Looking at the proportions, it is clear that there are two groups of predators (e.g., juveniles and adults), one group preys preferentially on Prey 1 while the other preys mostly on Prey 2. Let's see if we can pick this up by estimating individual predator proportions:

### Estiamting individual proportions

Lets try this with fatty acids first - given that there was limited information about diets in the stable isotopes alone, they can't be expected SI alone to provide the extra information we're after. 

We now need to deal with one extra prior to be set: the prior for the variance of diet proportions. After some exploratory runs, `even=2` seems like a reasonable compromise between our ability to detect differences (facilitated for smaller values of `even`) and the ability of the Markov Chains to mix properly (easier for larger `even` - see the documentation of `run_MCMC`).

#### From Fatty Acids


```r
Ind.FA <- run_MCMC(datas = dats.subset, nIter = 10000, nBurnin = 1000, nChains = 3, 
    nThin = 10, Data.Type = "Fatty.Acid.Profiles", Analysis.Type = "Individual.proportions", 
    Rnot = 1, even = 2, plott = F, spawn = T)
```

```
## +++++++++
```


We won't show the output of the next commands anymore since it is far too long. They are the same as for a population proportion analysis:

```
MCMCplot(Ind.FA)

diags(Ind.FA)

summary(Ind.FA)
```

```r
plot(Ind.FA, save = F)
```

```
## Using  as id variables
```

![plot of chunk unnamed-chunk-34](figure/unnamed-chunk-341.png) 

```
## Using rep(i, nrow(outs)) as id variables
```

![plot of chunk unnamed-chunk-34](figure/unnamed-chunk-342.png) ![plot of chunk unnamed-chunk-34](figure/unnamed-chunk-343.png) 


Despite this MCMC run being an order of magnitude too short, we can see that the grouped pattern emerges quite clearly in the last plot. Nevertheless, posterior correlations remain for prey items 1&2 remain strong (we'll soon see why!). First, let's see if the grouping becomes more pronounced with the inclusion of stable isotopes:

#### Combined analysis


```r
Ind.Combined <- run_MCMC(datas = dats.subset, nIter = 10000, nBurnin = 1000, 
    nChains = 3, nThin = 10, Data.Type = "Combined.Analysis", Analysis.Type = "Individual.proportions", 
    Rnot = 1, Rnot_SI = 1, even = 2, plott = F, spawn = T)
```

```
## ++++++++++
```

Again skipping diagnostics output:

```
MCMCplot(Ind.Combined)

diags(Ind.Combined)

summary(Ind.Combined)
```

```r
plot(Ind.Combined, save = F)
```

```
## Using  as id variables
```

![plot of chunk unnamed-chunk-36](figure/unnamed-chunk-361.png) 

```
## Using rep(i, nrow(outs)) as id variables
```

![plot of chunk unnamed-chunk-36](figure/unnamed-chunk-362.png) ![plot of chunk unnamed-chunk-36](figure/unnamed-chunk-363.png) 


The uncertainty for predators preying predominantly on Prey 2 is significantly reduced (but once again the analysis should be run for an order of magnitude longer with a larger thinning interval to ensure that the estimates are reliable).

While the patterns here offer great insights into individual diet proportions, they do not provide a means to estimate the population distribution of diet proportions for these two groups. This can be achieved in an anova type linear model:

### Estiamting group proportions

We now need to add covariates in the form of group membership. The reasoning and procedure is the same for continuous covariates (e.g., size). To add the covariates, we use the `add_Covs` constructor:


```r
group.path <- system.file("extdata", "Simdata_groups.csv", package = "fastinR")
Covs <- add_Covs(Groups = group.path)


Cov.Combined <- run_MCMC(datas = dats.subset, Covs = Covs, nIter = 10000, nBurnin = 1000, 
    nChains = 3, nThin = 10, Data.Type = "Combined.Analysis", Analysis.Type = "Analysis.with.Covariates", 
    Rnot = 1, Rnot_SI = 1, even = 2, plott = F, spawn = T)
```

```
## ++++++++++
```


```
MCMCplot(Cov.Combined)

diags(Cov.Combined)

summary(Cov.Combined)
```

```r
plot(Cov.Combined, save = F)
```

```
## Using  as id variables
## Using  as id variables
```

![plot of chunk unnamed-chunk-38](figure/unnamed-chunk-381.png) 

```
## Using rep(i, nrow(outs)) as id variables
```

![plot of chunk unnamed-chunk-38](figure/unnamed-chunk-382.png) ![plot of chunk unnamed-chunk-38](figure/unnamed-chunk-383.png) 


We can now clearly see the difference between the two group's diet proportions!
