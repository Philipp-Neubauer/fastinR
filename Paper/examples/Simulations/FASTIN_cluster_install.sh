git clone https://github.com/Philipp-Neubauer/FASTIN-R.git

cd FASTIN-R

echo 'install.packages(c("vegan", "fgui", "rjags", "robCompositions", "compositions", "denstrip"),repos="http://cran.stat.auckland.ac.nz/",dependencies=T);install.packages("FASTIN-R",repos=NULL)' > inst.pkg.R

R CMD BATCH inst.pkg.R
