fastinR
========

Fatty Acids and Stable Isotopes in R

Disclaimer
----------
This repository is still in development and will be updated frequently to eliminate bugs and add improvements. Please [file an issue](https://github.com/Philipp-Neubauer/fastinR/issues?milestone=1&state=open) if you find a bug or have a suggestion, that way it is visible to other users/contributors and progress on the bug/issue can be traced.

Install
-------

The package uses jags for Bayesian computations, JAGS needs to be installed manually from [here](http://sourceforge.net/projects/mcmc-jags/files/JAGS/3.x/).

Requirements: Some R package requirements, these should install automatically. If not, use R's ```install.package``` to install dependencies manually.

On MAC, Xquarz is needed to support the gui and plotting from the gui. Not needed iif plotting is manually turned off for all functions (this usually involves setting save=F manually in the function arguments).

To install directly from github you'll also need git (get it [here](http://git-scm.com/)) and the devtools package for R, and run ```require('devtools')``` once the devtools package is installed.

Please install development versions of the package directly from github using 
```R
install_github("philipp-neubauer/fastinR/fastinR")
```

You can also [clone the package](https://help.github.com/articles/fork-a-repo) and [install from source](http://stackoverflow.com/questions/1474081/how-do-i-install-an-r-package-from-source). This is the preferred way to contribute to code development.

Compiled versions will be released at some point too...

Running
-------

Calling ```fastinR_GUI()``` opens a gui that lets the user input relevant data and choose analysis options. Calling ```simulation()``` will open a gui to simulate data that can be input through the ```fastinR_GUI()```. Standalone functions are available for all options and are documented.

A tutorial using simulated data is available [here](http://figshare.com/articles/Estimating_diet_proportions_from_fatty_acids_and_stable_isotopes_the_fastinR_package/900392) 

Current Limitations
-------------------

 - plotting SI in dataplot function should just give SI as axes if n<=2 isotopes. 
 - multiplot only works for population proportion settings, not when estiamting individual proportions.
 - the violin multiplot only works for lists of length 3 (i.e, for comparing three methods). Should be useful in otehr circumstances too, so needs to be generalised.
 - the default prior setup works ok in some cases but needs lots of manual tweaking in others. Better default options should be possible, but will require further work.
 - The jags backend is very inefficient for data with lots of dimensions (```predators*preys*fatty acids```), a custom backend would probably be better but would be significantly more work
