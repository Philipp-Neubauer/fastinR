Analysis of Stowasser et al. 2006 data
========================================================

Here we demonstrate the validity of our methods by re-analysing the experimental dataset of Stowasser et al. The original paper investigated the evolution of stable isotope (C/N) and fatty acid signatures for squid fed either fish, crustacean or mixed diets, as well as switched diet regimes. Our specific aim here is to estiamte diet proportions of mixed diet individuals, and to show the advantage of using both markers in concert over using a single marker for diet analysis. 

We begin the analysis using just SI. Hussey et al 2014 (Ecology Letters Volume 17, Issue 2, pages 239–250, February 2014), following earlier analyses by Caut et al and others, showed that discimintation is prey N dependent, and used a meta-analytic model to estiamte discrimination coefficients. 

We initially estiamte prey specific discimination in stable isotopes following Hussey et al 2014, and develop analogous methods for fatty acid signatures. Discrimination is estiamted for both types of marker from animals that were fed on diets containing a single prey type (Fish or Crustaceans).

Stable Isotope analysis
------------------------

We first read in the data tables for SI:


```
## Loading required package: dplyr
## 
## Attaching package: 'dplyr'
## 
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## 
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```



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

```r
require(dplyr)

prey.table <- read.csv("Prey_SI.csv", header = T, stringsAsFactors = F)
pred.table <- read.csv("Predator_SI.csv", header = T, stringsAsFactors = F)

# use only mixed feed individuals for analysis - using dplyr
pred.table <- pred.table %.% filter(Feed == "Mixed") %.% select(X_13C, X_15N)
```

We estimate discrimintation coefficients (also termed fractionation or conversion for SI and FA, respectively) using Bayesian models. The analysis is done in a separate file (Discrim.Analysis.R) and produces files discr.means.csv/discr.var.csv for SI, and corresponding cc_FA.csv and cc_FA_var.csv for FA.The estiamtion isn't trivial, so we just pass over it here. We instead proceed with the analysis by adding the necessary items for the final analysis to the workspace, using the add_SI function in fastin-R. This function just adds the data, and puts it into a specific format.


```r
squid.data <- add_SI(SI.predators = pred.table, SI.preys = prey.table, Frac.Coeffs.mean = "discr.means.csv", 
    Frac.Coeffs.var = "discr.var.csv")
```


This is all we need for a first analysis using SI data. We can now visualize the data on an Non-Metric Dimensional Scaling plot like so:


```r
dataplot(squid.data)
```

```
## Run 0 stress 0 
## Run 1 stress 9.944e-05 
## ... procrustes: rmse 0.002492  max resid 0.007928 
## *** Solution reached
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 


SIs seem to provide some separation among diet items, even within fish, but there is still strong overlap. One way to find out if the analysis of SI alone yields anything interpretable.

The run_MCMC function is the main workhorse for any type of analysis, it uses methods dispatch (=different models) based on function arguments and data types. Notice that the `Data.Type` argument is set to `'Stable.Isotopes'` and the `Analysis.Type` is `'Population.proportions'`. Lets not worry about the rest for now.


```r
Squid.SI.run1 <- run_MCMC(datas = squid.data, nIter = 10000, nBurnin = 1000, 
    nChains = 3, nThin = 10, Data.Type = "Stable.Isotopes", Analysis.Type = "Population.proportions", 
    Rnot_SI = 0.01, plott = F, spawn = F)
```

```
## Compiling model graph
##    Resolving undeclared variables
##    Allocating nodes
##    Graph Size: 114
## 
## Initializing model
## 
## 
##  proceeding to burn-in phase 
## 
##  sampling from parameter distributions
```


Plotting the MCMC is the easiest way to ensure that the sampler is mixing - meaning that the Markov chain explores the posterior distribution of each parameter efficiently.


```r
MCMCplot(Squid.SI.run1)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 


This looks rather uninformative, let's see if this is due to chains not running long enough. The `diags` function gives more information, it dispalys two types of diagnostics (from coda package):


```r
diags(Squid.SI.run1)
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
##  prop[1] 20       9690  937          10.3      
##  prop[2] 30       10110 937          10.8      
##  prop[3] 20       9690  937          10.3      
##  prop[4] 30       10530 937          11.2      
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
##  prop[1] 40       12430 937          13.30     
##  prop[2] 20       8930  937           9.53     
##  prop[3] 30       11430 937          12.20     
##  prop[4] 20       8930  937           9.53     
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
##  prop[1] 30       10530 937          11.2      
##  prop[2] 30       10620 937          11.3      
##  prop[3] 30       11430 937          12.2      
##  prop[4] 30       11430 937          12.2      
## 
## 
## 
##  Based on these diagnostics you should repeat the MCMC with  12430  iterations 
##  and a thinning interval of  13  ,if these values are higher than the values 
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
## prop[1]          1       1.01
## prop[2]          1       1.01
## prop[3]          1       1.01
## prop[4]          1       1.00
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


Lets run it again for a bit longer, this time using spawn=T, which 'spawns' invisible R processes that do the calculations


```r
Squid.SI.run2 <- run_MCMC(datas = squid.data, nIter = 1e+05, nBurnin = 10000, 
    nChains = 3, nThin = 100, Data.Type = "Stable.Isotopes", Analysis.Type = "Population.proportions", 
    Rnot_SI = 0.01, plott = F, spawn = T)

diags(Squid.SI.run2)
```



```r
MCMCplot(Squid.SI.run2)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 

Looks like that's about as good as it'll get, so let's look at a summary and produce a nicer plot for the outcome.


```r
summary(Squid.SI.run2)
```



```r
plot(Squid.SI.run2, save = F)
```

```
## Using  as id variables
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-111.png) ![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-112.png) 


Based on posterior medians, we'd say that Mixed diet treatment squid were mostly feeding on fish (at the level of the predator population), but there's considerable uncertainty. Does this reflect variability in the individual diets or uncertainty due to lack of signal in the data? We can set the `Analysis.Type` to `'Individual.proportions'` to see how individual diets vary and get a better idea.


```r
Squid.SI.run3 <- run_MCMC(datas = squid.data, nIter = 1e+06, nBurnin = 10000, 
    nChains = 3, nThin = 1000, Data.Type = "Stable.Isotopes", Analysis.Type = "Individual.proportions", 
    Rnot_SI = 10, plott = F, spawn = T)

diags(Squid.SI.run3)
```



```r
MCMCplot(Squid.SI.run3)
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-131.png) ![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-132.png) ![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-133.png) ![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-134.png) ![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-135.png) ![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-136.png) 



```r
summary(Squid.SI.run3)
```


Scolling through this summary, you should see both individual psoterior sumamries and population level estimates for all chains. 


```r
plot(Squid.SI.run3, save = F)
```

```
## Using  as id variables
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-151.png) 

```
## Using rep(i, nrow(outs)) as id variables
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-152.png) ![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-153.png) 


This still looks very uncertain, and we could play with Rnot to see if we can get a better prior that improves mixing and perhaps the estiamtes. I tried that allready, and it doesn't work, which means that the uncertainty in the SI analysis is irreducible.

Let's switch to fatty acids and see what that data can reveal independently of SIs.

FAtty Acid Analysis
-------------------

First, we need to calculate the moments for fat content of prey types. The FA prey data has 3 profiles for striped mullet - since these are similar to other fish and were fed in very low proportions, we'll exclude them here. We do not have the original fat content data (could ask for it?), so we just use the empirical mean and variance to calculate priors for a log-normal model of fat content:


```r

# load fatty acid data tables

# prey
prey.ix <- t(read.csv("Prey_FA.csv", header = F, stringsAsFactors = F, row.names = 1))[, 
    1]
prey.table.FA <- t(read.csv("Prey_FA.csv", header = T, stringsAsFactors = F, 
    row.names = 1))

mullets <- which(prey.ix == "Striped Mullet")
prey.ix <- prey.ix[-mullets]
prey.table.FA <- prey.table.FA[-mullets, ]

# need to replace 0 proportions in dataset by min for that FA (ad-hoc..), no
# need to worry if it doesn't sum to 100 anymore
for (i in 1:ncol(prey.table.FA)) prey.table.FA[prey.table.FA[, i] == 0, i] = min(prey.table.FA[prey.table.FA[, 
    i] > 0, i])

# cumbersome, but need to add column of prey.ix
prey.table.FA <- data.frame(prey.ix, prey.table.FA)

# predators - subset to Mixed feed only taht have SI data (to be quicker)
pred.table.FA <- t(read.csv("Predator_FA.csv", header = T, stringsAsFactors = F, 
    row.names = 1))
pred.table.FA <- pred.table.FA[grep("Mixed", rownames(pred.table.FA))[c(1, 2, 
    3, 13, 14)], ]
# same here, replace
for (i in 1:ncol(pred.table.FA)) pred.table.FA[pred.table.FA[, i] == 0, i] = min(pred.table.FA[pred.table.FA[, 
    i] > 0, i])

Silverside <- 2.27
Silverside.sd <- 1.02

Sailfin.Molly <- 3.97
Sailfin.Molly.sd <- 1.25

Sheepshead.Minnow <- 3.59
Sheepshead.Minnow.sd <- 1.45

fish.means <- c(Silverside, Sailfin.Molly, Sheepshead.Minnow)

fish.vars <- c(Silverside.sd, Sailfin.Molly.sd, Sheepshead.Minnow.sd)^2
```


Assuming a log normal model for fat content - this is done internally in fastin-R, so the next bit is for illustrative purposes only...


```r
fish.ln.vars = log(fish.vars + fish.means^2) - 2 * log(fish.means)
fish.ln.means = log(fish.means) - fish.ln.vars/2

hist(rlnorm(10000, meanlog = fish.ln.means[1], sdlog = sqrt(fish.ln.vars[1])), 
    30)
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-171.png) 

```r

shrimp.mean <- 2.43
shrimp.var <- 0.66^2

shrimp.ln.var = log(shrimp.var + shrimp.mean^2) - 2 * log(shrimp.mean)
shrimp.ln.mean = log(shrimp.mean) - shrimp.ln.var/2

hist(rlnorm(10000, meanlog = shrimp.ln.mean, sdlog = sqrt(shrimp.ln.var)), 30)
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-172.png) 

```r

fat.cont <- rbind(cbind(fish.ln.means, fish.ln.vars), c(shrimp.ln.mean, shrimp.ln.var))
colnames(fat.cont) <- NULL
rownames(fat.cont) <- unique(prey.ix)

write.table(fat.cont, file = "fat.cont.csv", col.names = F, sep = ",")
```


We estimate discrimintation coefficients (also termed fractionation or conversion for SI and FA, respectively) using Bayesian models. The analysis is done in a separate file (Discrim.Analysis.R) and produces files discr.means.csv/discr.var.csv for SI, and corresponding cc_FA.csv and cc_FA_var.csv for FA.The estiamtion isn't trivial for FA, so we just pass over it here. We instead proceed with the analysis by adding the necessary items for the FA analysis to the workspace.


```r
# note the datas argument here taking our initial SI object
squid.data.FA <- add_FA(FA.predators = pred.table.FA, FA.preys = prey.table.FA, 
    fat.conts = "fat.cont.csv", Conv.Coeffs.mean = "cc_FA.csv", Conv.Coeffs.var = "cc_FA_var.csv", 
    LN.par = T)
```


Plotting the new dataset in MDS-scaled space:


```r
dataplot(squid.data.FA)
```

```
## Run 0 stress 0.1333 
## Run 1 stress 0.1467 
## Run 2 stress 0.1475 
## Run 3 stress 0.1475 
## Run 4 stress 0.1326 
## ... New best solution
## ... procrustes: rmse 0.0217  max resid 0.08105 
## Run 5 stress 0.1262 
## ... New best solution
## ... procrustes: rmse 0.07705  max resid 0.5188 
## Run 6 stress 0.1262 
## ... procrustes: rmse 0.0001232  max resid 0.0008824 
## *** Solution reached
```

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19.png) 


MDS sugegsts that FA alone may not allow us to accurately estiamte diet proportions either. To confirm this, let's select a subset of FA for analysis first to cut down on computation time. This is done using the ```select_vars``` function.


```r
squid.subset <- select_vars(squid.data.FA)
```

```
##           Cumulative Proportion
## X22.6n.3                 0.6450
## X18.1n.9                 0.7593
## X20.5n.3                 0.8406
## X16.1n.7                 0.8899
## X18.2n.6                 0.9282
## X22.5n.3                 0.9547
## X18.1n.7                 0.9743
## X18.0                    0.9877
## X18.3n.3                 0.9921
## X16.0                    0.9955
## X24.1n.9                 0.9967
## X20.4n.6                 0.9975
## X20.4n.3                 0.9981
## X20.1n.9                 0.9987
## X18.4n.3                 0.9991
## X16.3n.6                 0.9994
## X14.0                    0.9995
## X18.3n.6                 0.9997
## X15.0                    0.9998
## X16.4n.3                 0.9999
## X21.5n.3                 0.9999
## X20.1n.11                1.0000
## X22.1n.9                 1.0000
## X12.0                    1.0000
## X16.2n.6                 1.0000
```

```
## Error: select.list() cannot be used non-interactively
```

![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-20.png) 


The two plots (and console output) suggest that using FAs up to 16:0 will preserve 99.5% of the total cumulative separation (based on all FAs) between preys, and keeps the prey matrix-condition number low, meaning less collinearity among items. Great!


```r
dataplot(squid.subset)
```

```
## Error: object 'squid.subset' not found
```


This doesn't look much better, so let's see what a first stab at the analysis spits out...


```r
Squid.subset.analysis <- run_MCMC(datas = squid.subset, nIter = 30000, nBurnin = 10000, 
    nChains = 3, nThin = 30, Data.Type = "Fatty.Acid.Profiles", Analysis.Type = "Population.proportions", 
    Rnot = 0.2, plott = F, spawn = T)
```

```
## Error: object 'squid.subset' not found
```

```r

diags(Squid.subset.analysis)
```

```
## Error: object 'Squid.subset.analysis' not found
```



```r
MCMCplot(Squid.subset.analysis)
```

```
## Error: object 'Squid.subset.analysis' not found
```



```r
summary(Squid.subset.analysis)
```

```
## Error: object 'Squid.subset.analysis' not found
```



```r
plot(Squid.subset.analysis, save = F)
```

```
## Error: object 'Squid.subset.analysis' not found
```


FAs suggest a greater importance of silversides in the diets of mixed diet squid, but there's still lots of uncvertainty. Let's compare more explicitly using the ```multiplot``` function:


```r
Pop.list <- list(Stable.Isotopes = Squid.SI.run2, Fatty.Acids = Squid.subset.analysis)
```

```
## Error: object 'Squid.subset.analysis' not found
```

```r

multiplot(Pop.list, save = F)
```

```
## Error: object 'Pop.list' not found
```


Calculating individual proportions with FA data alone works in the same way as before. The R.not prior is often key to get well bahving chains.


```r
Squid.subset.ind <- run_MCMC(datas = squid.subset, nIter = 30000, nBurnin = 10000, 
    nChains = 3, nThin = 30, Data.Type = "Fatty.Acid.Profiles", Analysis.Type = "Individual.proportions", 
    Rnot = 1, plott = F, spawn = T)
```

```
## Error: object 'squid.subset' not found
```

```r

diags(Squid.subset.ind)
```

```
## Error: object 'Squid.subset.ind' not found
```



```r
MCMCplot(Squid.subset.ind)
```

```
## Error: object 'Squid.subset.ind' not found
```



```r
summary(Squid.subset.ind)
```

```
## Error: object 'Squid.subset.ind' not found
```



```r
plot(Squid.subset.ind, save = F)
```

```
## Error: object 'Squid.subset.ind' not found
```

Looking at the individual (still preliminary!) psoterior of proportions, it seems that some of the uncertainty at the population level may stem from the dominance of Sheepshead minnow for squid 1,2 and 5, whereas Silverside dominantes for squid no 3 and 4. NEvertheless, a lot of uncertainty remains.

Lets see is a combined analysis can improve on this. We need to combine the two data objects first.

Combined analysis for population proportions
--------------------------------------------


```r
# note the datas argument here taking our initial SI object
squid.data.comb <- add_FA(FA.predators = pred.table.FA, FA.preys = prey.table.FA, 
    fat.conts = "fat.cont.csv", Conv.Coeffs.mean = "cc_FA.csv", Conv.Coeffs.var = "cc_FA_var.csv", 
    datas = squid.data, LN.par = T)
```


And take the subset again to run a combined analysis:


```r
squid.subset.comb <- select_vars(squid.data.comb)
```

```
## Error: select.list() cannot be used non-interactively
```

![plot of chunk unnamed-chunk-32](figure/unnamed-chunk-32.png) 

```r

Squid.combined.analysis <- run_MCMC(datas = squid.subset.comb, nIter = 30000, 
    nBurnin = 10000, nChains = 3, nThin = 30, Data.Type = "Combined.Analysis", 
    Analysis.Type = "Population.proportions", plott = F, spawn = T)
```

```
## Error: object 'squid.subset.comb' not found
```

```r

diags(Squid.combined.analysis)
```

```
## Error: object 'Squid.combined.analysis' not found
```



```r
MCMCplot(Squid.combined.analysis)
```

```
## Error: object 'Squid.combined.analysis' not found
```



```r
summary(Squid.combined.analysis)
```

```
## Error: object 'Squid.combined.analysis' not found
```



```r
plot(Squid.combined.analysis, save = F)
```

```
## Error: object 'Squid.combined.analysis' not found
```

The results seem to make sense: the feed for the original study included 60 silversides, 45 Sheepshead minnows and 50 Sailfin mollys, with 76 shrimp over all treatments. Assuming equal allocation suggests that the inferred proportions are quite likely to be close to the actual diet proportions. Lets see how different they are between markers and the combined analysis:


```r
Pop.list <- list(Stable.Isotopes = Squid.SI.run2, Fatty.Acids = Squid.subset.analysis, 
    Combined.Analysis = Squid.combined.analysis)
```

```
## Error: object 'Squid.subset.analysis' not found
```

```r

multiplot(Pop.list, save = F)
```

```
## Error: object 'Pop.list' not found
```


Combined analysis for individual proportions
--------------------------------------------


```r
Squid.combined.ind <- run_MCMC(datas = squid.subset.comb, nIter = 30000, nBurnin = 10000, 
    nChains = 3, nThin = 30, Data.Type = "Combined.Analysis", Analysis.Type = "Individual.proportions", 
    Rnot = 0.2, Rnot_SI = 10, plott = F, spawn = T)
```

```
## Error: object 'squid.subset.comb' not found
```

```r

diags(Squid.combined.ind)
```

```
## Error: object 'Squid.combined.ind' not found
```



```r
MCMCplot(Squid.combined.ind)
```

```
## Error: object 'Squid.combined.ind' not found
```



```r
summary(Squid.combined.ind)
```

```
## Error: object 'Squid.combined.ind' not found
```



```r
plot(Squid.combined.ind, save = F)
```

```
## Error: object 'Squid.combined.ind' not found
```

