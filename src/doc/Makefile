### This directory contains the sources for the manual "Multi-state
### modelling with R: the msm package" which is included in the source
### package as inst/doc/msm-manual.pdf, and in the installed package as
### doc/msm-manual.pdf.
### 
### To build msm-manual.pdf, first the msm package must be
### installed. Then type
### 
### make 
###
### Or equivalently, open R in this directory and run
### 
### Sweave("msm-manual.Rnw")
### 
### Exit from R.  The resulting LaTeX file, msm-manual.tex, can be
### compiled to PDF with the shell command 
### 
### texi2dvi msm-manual.tex 

LATEX=msm-manual.tex
PDF=msm-manual.pdf

all: $(LATEX) $(PDF)

$(LATEX):
	R --vanilla --slave -e Sweave\(\'msm-manual.Rnw\'\)

$(PDF):
	texi2dvi msm-manual.tex 

clean:
	rm -f msm-manual*.pdf msm-manual.tex msm-manual*.eps *.dvi *.toc *.aux *.bbl *.blg *.log *~ Rplots.pdf
