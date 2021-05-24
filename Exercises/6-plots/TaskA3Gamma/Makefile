CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)


gamma.pyxplot.png: data.txt Makefile gammatabulated.txt out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Gamma function";\
	plot \
	 "$<" using 1:2 with line title "Gamma function from math.h"\
	,"$<" using 1:3 with line title "GSL gamma function"\
	,"$<" using 1:4 with line title "Stirling Approximation in gamma.c"\
	,"gammatabulated.txt" using 1:2 with points title "Tabulated gamma values"\
	' | pyxplot


out.txt: main
	./$< > $@

data.txt: main

main: main.o gamma.o

.PHONEY:clean
clean:
	$(RM) *.o main *.png data.txt *.Identifier *.save 
