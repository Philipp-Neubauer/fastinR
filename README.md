fastinR
========

Fatty Acids and Stable Isotopes in R

Install
=======

Requirements: Some R package requirements, these should install automatically. On MAC, Xquarz is needed (I think...) to support the gui and plotting from the gui. Not needed iif plotting is manually turned off for all functions.

To install directly from github you'll also need git (get it here: http://git-scm.com/) and the devtools package for R.

Please install development versions of the package directly from github using 
```R
install_github("philipp-neubauer/fastinR/fastinR")
```

Running
===========

Calling ```fastinR_GUI()``` opens a gui that lets the user input relevant data and choose analysis options. Calling ```simulation()``` will open a gui to simulate data that can be input through the ```fastinR_GUI()```. Standalone functions are available for all options and are documented, vignette to come.
