TEXINPUTS:=.

HW=pdf/hw1.pdf \
   pdf/hw2.pdf \
   pdf/hw3.pdf \
   pdf/hw4.pdf \
   pdf/hw5.pdf \
   pdf/mt.pdf \
   pdf/hw6.pdf \
   pdf/hw7.pdf \
   pdf/hw8.pdf \
   pdf/hw9.pdf \
   pdf/hw10.pdf \
   pdf/final.pdf

.PHONY: all copy
all: $(HW)

pdf/%.pdf: tex/%.tex
	xelatex -output-dir pdf $<
	xelatex -output-dir pdf $<

clean:
	rm -f pdf/*.aux pdf/*.log pdf/*.out
	rm -f pdf/*.nav pdf/*.snm pdf/*.toc
	rm -f pdf/*.vrb

distclean: clean
	rm -f pdf/*.pdf

copy: $(HW)
	mkdir -p ../web/hw
	cp $(HW) ../web/hw
