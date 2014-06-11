all: Simulated_example.pdf

Simulated_example.pdf: Simulated_example.tex
	grep -vwE "(\\#\\#)" Simulated_example_pre.tex > Simulated_example.tex;pdflatex Simulated example.tex

Simulated_example.tex:Simulated_example.Rnw
	Rscript --vanilla -e 'library(knitr);knit("Simulated_example.Rnw")';mv Simulated_example.tex Simulated_example_pre.tex	
