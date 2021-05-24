CFLAGS = $(shell gsl-config --cflags) 
LDLIBS = $(shell gsl-config --libs)

erf.pyxplot.png: data.txt Makefile tabulated.data.txt out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function";\
	plot \
	 "$<" using 1:2 with line title "Error function from math.h"\
	,"$<" using 1:3 with line title "Error function from gsl"\
	,"$<" using 1:4 with line title "Error Function approximation in erf.c"\
	,"tabulated.data.txt" using 1:2 with points title "Tabulated data"\
	' | pyxplot


out.txt: main 
	./$< > $@
    
data.txt: main

main: main.o erf.o 

.PHONEY: clean 
clean: 
	$(RM) *.o main *.png data.txt *.Identifier *.save
