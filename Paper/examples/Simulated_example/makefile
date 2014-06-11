all: Simulated_Example.pdf

Simulated_example.pdf: Simulated example.tex
	pdflatex Simulated example.tex

Simulated example.tex:Simulated example.Rnw
	Rscript --vanilla -e 'library(knitr);knit("Simulated_example.Rnw")'	
