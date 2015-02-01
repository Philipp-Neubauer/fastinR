all: Manuscript_draft.pdf S1-Simulation.pdf S2-Sensitivity.pdf  S3-Variable_selection.pdf S4-Discrimination_analysis.pdf S5-Squid_re-analysis.pdf

Manuscript_draft.pdf: Main_paper.tex 
	xelatex Main_paper.tex; biber Main_paper; xelatex Main_paper.tex; xelatex Main_paper.tex; mv Main_paper.pdf > Manuscript_draft.pdf

# get simualtion tex from relevant dir, need to run locally to repeat analysis
S1-Simulation.pdf: ../../examples/Simulated_example.tex
	 pdflatex ../../examples/Simulated_example.tex ; mv  ../../examples/Simulated_example.pdf > S1-Simulation.pdf

S2-Sensitivity.pdf: Appendix_1.tex
	xelatex Appendix_1.tex; biber Appendix_1; xelatex Appendix_1.tex; xelatex Appendix_1.tex; mv Appendix_1.pdf S2-Sensitivity.pdf

S3-Variable_selection.pdf: Appendix_2.tex
	xelatex Appendix_2.tex; biber Appendix_2; xelatex Appendix_2.tex; xelatex Appendix_2.tex; mv Appendix_2.pdf S3-Variable_selection

S4-Discrimination_analysis.pdf: ../../examples/Discrimination_analysis.tex
	pdflatex ../../examples/Discrimination_analysis.tex ; mv  ../../examples/Discrimination_analysis.pdf > S4 -Discrimination_analysis.pdf

# get Squid analysis tex from relevant dir, need to run locally to repeat analysis
S5-Squid_re-analysis.pdf: ../../examples/Squid_re-analysis.tex	
	pdflatex ../../examples/Squid_re-analysis.tex ; mv  ../../examples/Squid_re-analysis.pdf > S1-Simulation.pdf

clean:
	rm -fv *.aux *.log *.toc *.blg *.bbl *.synctex.gz
	rm -fv *.out *.bcf *blx.bib *.run.xml
	rm -fv *.fdb_latexmk *.fls
