CFLAGS = $(shell gsl-config --cflags) 
LDLIBS = $(shell gsl-config --libs)

erf.pyxplot.png: data.txt Makefile tabulated.data.txt 
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function";\
	plot \
	 "$<" using 1:2 with line title "Error Function from math.h"\
	,"$<" using 1:3 with line title "Error Funktion from gsl: gsl-sf-erf"\
	,"$<" using 1:4 with line title "Homemade Error Function"\
	,"tabulated.data.txt" using 1:2 with points title "Tabulated data for test"\
	' | pyxplot


data.txt: main 
	./main > data.txt

main: main.o erf.o 

.PHONEY: clean 
clean: 
	$(RM) *.o main *.png data.txt 
