TEXINPUTS:=.

LEC= \
	pdf/background.pdf \
	pdf/2019-08-30.pdf \
	pdf/2019-09-04.pdf \
	pdf/2019-09-06.pdf \
	pdf/2019-09-09.pdf \
	pdf/2019-09-11.pdf \
	pdf/2019-09-13.pdf \
	pdf/2019-09-16.pdf \
	pdf/2019-09-18.pdf \
	pdf/2019-09-20.pdf \
	pdf/2019-09-23.pdf \
	pdf/2019-09-25.pdf \
	pdf/2019-09-27.pdf \
	pdf/2019-09-30.pdf \
	pdf/2019-10-02.pdf \
	pdf/2019-10-04.pdf \
	pdf/2019-10-07.pdf \
	pdf/2019-10-09.pdf \
	pdf/2019-10-11.pdf \
	pdf/2019-10-16.pdf \
	pdf/2019-10-18.pdf \
	pdf/2019-10-21.pdf \
	pdf/2019-10-23.pdf \
	pdf/2019-10-25.pdf \
	pdf/2019-10-28.pdf \
	pdf/2019-10-30.pdf \
	pdf/2019-11-01.pdf \
	pdf/2019-11-04.pdf \
	pdf/2019-11-06.pdf \
	pdf/2019-11-08.pdf \
	pdf/2019-11-13.pdf \
	pdf/2019-11-15.pdf 

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
