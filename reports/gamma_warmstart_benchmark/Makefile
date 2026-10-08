# Makefile for gamma_warmstart_benchmark report

all: plots report.pdf

plots:
	gnuplot plots/plot_runtime.gp
	gnuplot plots/plot_speedup.gp
	gnuplot plots/plot_gamma.gp
	gnuplot plots/plot_correlations.gp

report.pdf: report.tex plots
	pdflatex -interaction=nonstopmode report.tex
	pdflatex -interaction=nonstopmode report.tex

clean:
	rm -f *.aux *.log *.out *.toc report.pdf

.PHONY: all plots clean
