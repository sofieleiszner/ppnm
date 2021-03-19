CFLAGS = -Wall -O1 -std=gnu11
LDLIBS = -lm

ms.pdf: ms.tex fig-pyxplot.pdf
	pdflatex $<
	pdflatex $<
    
out.txt: main 
	./$< > $@
    
main: main.o myexp.o


fig-pyxplot.pdf: fig-pyxplot.gpi out.txt 
	pyxplot $<
    
    
clean: 
	$(RM) *.o fig-pyxplot.pdf out.txt ms.aux ms.pdf ms.txt ms.log