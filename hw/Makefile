TEXINPUTS:=.

HW=pdf/hw1.pdf 

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
	#mkdir -p ../web/hw
	#cp $(HW) ../web/hw
