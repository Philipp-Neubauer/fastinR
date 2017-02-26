fastinR
========

[![DOI](https://zenodo.org/badge/16589/Philipp-Neubauer/fastinR.svg)](https://zenodo.org/badge/latestdoi/16589/Philipp-Neubauer/fastinR)


Fatty Acids and Stable Isotopes in R - the published paper is available @ [peerj](https://peerj.com/articles/920/)

Disclaimer
----------
This repository is still in development and will be updated **in**frequently to eliminate bugs and add improvements. Please [file an issue](https://github.com/Philipp-Neubauer/fastinR/issues?milestone=1&state=open) if you find a bug or have a suggestion, that way it is visible to other users/contributors and progress on the bug/issue can be traced.

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

You can also [fork and clone the package](https://help.github.com/articles/fork-a-repo) and [install from source](http://stackoverflow.com/questions/1474081/how-do-i-install-an-r-package-from-source). This is the preferred way to contribute to code development.

Compiled versions will be released at some point too...

Running
-------

Calling ```fastinR_GUI()``` opens a gui that lets the user input relevant data and choose analysis options. Calling ```simulation()``` will open a gui to simulate data that can be input through the ```fastinR_GUI()```. NOTE that the GUI is depricated and will no longer be maintained - the TCL/TK GUI interface is far too much of a pain to work with to make maintenance worth anyone's while. That said, it should work in general, but might produce some obscure messages and such.

Standalone functions are available for all options and are documented.

A tutorial using simulated data accompanies the original paper [as supplement 1](https://peerj.com/articles/920/#supp-1), and code for all analyses perfomed for the paper are [available from peerJ as well](https://peerj.com/articles/920/#supplemental-information).

Current Limitations
-------------------

Minor
**********

 - plotting SI in dataplot function should just give SI as axes if n<=2 isotopes. 
 - multiplot only works for population proportion settings, not when estiamting individual proportions.
 - the violin multiplot only works for lists of length 3 (i.e, for comparing all three methods - SI, FA and combined). Should be useful in other circumstances too, so needs to be generalised.
 
More Substantial
*************

 - the fully Bayesian approach mean strong constraints on the number of FA x prey x predator that can be run. The Empirical Bayes approach used in most SI models should mitigate this (at a price, of course). First priority...
 - the default prior setup works ok in some cases but needs lots of manual tweaking in others. Better default options should be possible, but will require further work.
 - The jags backend is very inefficient for data with lots of dimensions (```predators*preys*fatty acids```), a custom backend would probably be better but would be significantly more work
