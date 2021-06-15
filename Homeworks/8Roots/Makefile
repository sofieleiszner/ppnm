CFLAGS = -Wall -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

default: Hydrogen.pyxplot.png C1.pyxplot.png C2.pyxplot.png C1AndC2.pyxplot.png

Hydrogen.pyxplot.png: hydrogenatom.txt main Makefile out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Hydrogenatom";\
	plot \
	"hydrogenatom.txt" using 1:2 with points title "Calculated"\
    , x*exp(-x) with lines title "Exact Result"\
	' | tee tmp.gpi | pyxplot
    
    
C1.pyxplot.png: hydrogenatom.txt main Makefile out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "rmax";\
	set ylabel "Energy";\
	set title "C1";\
	plot \
	"C1.txt" using 1:2 with points title "Calculated"\
    , -0.5 with lines title "Exact Result"\
	' | tee tmp.gpi | pyxplot

    
C2.pyxplot.png: hydrogenatom.txt main Makefile out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "rmax";\
	set ylabel "Energy";\
	set title "C2";\
	plot \
	"C2.txt" using 1:2 with points title "Calculated"\
    , -0.5 with lines title "Exact Result"\
	' | tee tmp.gpi | pyxplot

C1AndC2.pyxplot.png: hydrogenatom.txt main Makefile out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "rmax";\
	set ylabel "Energy";\
	set title "Comparisson of C1 and C2";\
	plot \
	"C1.txt" using 1:2 with points title "C1"\
	,"C2.txt" using 1:2 with points title "C2"\
    , -0.5 with lines title "Exact Result"\
	' | tee tmp.gpi | pyxplot



out.txt: main 
	./$< > $@
    
hydrogenatom.txt: main    

    
main: ODEfunctions.o functions.o main.o  functionsLinEq.o
   
    
clean: 
	$(RM) main *.o *.txt *.png *.gpi
    
    