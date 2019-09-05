TEXINPUTS:=.

LEC= \
	pdf/background.pdf \
	pdf/2019-08-30.pdf \
	pdf/2019-09-04.pdf \
	pdf/2019-09-06.pdf 

.PHONY: all copy
all: $(LEC)

pdf/%.pdf: notes/%.tex
	xelatex -output-dir pdf $<
	xelatex -output-dir pdf $<

pdf/2019-09-04.pdf: notes/2019-09-04.tex
	xelatex -output-dir pdf $<
	bibtex pdf/2019-09-04
	xelatex -output-dir pdf $<
	xelatex -output-dir pdf $<

pdf/2019-08-30.pdf: notes/2019-08-30.tex
	xelatex -output-dir pdf $<
	bibtex pdf/2019-08-30
	xelatex -output-dir pdf $<
	xelatex -output-dir pdf $<
	
clean:
	rm -f pdf/*.aux pdf/*.log pdf/*.out
	rm -f pdf/*.nav pdf/*.snm pdf/*.toc 
	rm -f pdf/*.vrb

distclean: clean
	rm -f pdf/*.pdf

copy:
	mkdir -p ../web/lec
	cp $(LEC) ../web/lec
