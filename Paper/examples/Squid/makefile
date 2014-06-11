all: Squid_re-analysis.pdf

Squid_re-analysis.pdf: Squid_re-analysis_pre.tex
	grep -vwE "(\\#\\#)" Squid_re-analysis_pre.tex > Squid_re-analysis.tex; pdflatex Squid_re-analysis.tex

Squid_re-analysis_pre.tex: Squid_re-analysis.Rnw
	Rscript -e 'library(knitr);knit("Squid_re-analysis.Rnw")';mv Squid_re-analysis.tex Squid_re-analysis_pre.tex	
