# Build the manuscript and the supplementary material.
#
#   make        both documents
#   make paper  main manuscript
#   make supp   supplementary material
#   make clean  remove LaTeX by-products

LATEX  := pdflatex -interaction=nonstopmode -halt-on-error
BIBTEX := bibtex

.PHONY: all paper supp clean distclean

all: paper supp

paper: pdf/paper.pdf

supp: pdf/supplementary.pdf

pdf/paper.pdf: paper.tex referencias.bib $(wildcard figures/*.pdf)
	@mkdir -p pdf
	$(LATEX) paper.tex
	$(BIBTEX) paper
	$(LATEX) paper.tex
	$(LATEX) paper.tex
	@mv paper.pdf pdf/paper.pdf
	@echo "==> pdf/paper.pdf"

pdf/supplementary.pdf: supplementary.tex referencias.bib $(wildcard figures/*.pdf)
	@mkdir -p pdf
	$(LATEX) supplementary.tex
	$(BIBTEX) supplementary
	$(LATEX) supplementary.tex
	$(LATEX) supplementary.tex
	@mv supplementary.pdf pdf/supplementary.pdf
	@echo "==> pdf/supplementary.pdf"

clean:
	rm -f *.aux *.log *.bbl *.blg *.out *.spl *.toc *.synctex.gz

distclean: clean
	rm -rf pdf cache cache_dF
