CFLAGS = -Wall -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

least.pyxplot.png: datafile.txt data.txt Makefile main out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "Time (days)";\
	set ylabel "ln(Activity)";\
	set title "Least squares fit";\
	plot \
	"datafile.txt" using 1:2 with lines title "Our fit"\
	,"datafile.txt" using 1:3 with lines title "Our fit (+δck)"\
	,"datafile.txt" using 1:4 with lines title "Our fit (-δck)"\
    ,"data.txt" using 1:4:5 with yerrorbars title "Data with errorbars"\
	' | tee tmp.gpi | pyxplot

out.txt: main 
	./$< > $@

datafile.txt: main

main: main.o GS_decomp.o  GS_solve.o inverse.o leastsq.o

main.o GS_decomp.o MeasureDecompTime.o: GS_decomp.c GS_solve.c inverse.c leastsq.c
        
    
clean: 
	$(RM) main main.o GS_decomp.o out.txt least.pyxplot.png datafile.txt