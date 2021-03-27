CFLAGS = -Wall -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

default: TheWay.pyxplot.png TheWaySIR.pyxplot.png TheWayNewton.pyxplot.png

TheWay.pyxplot.png: TheWay.txt main Makefile out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Sine function example";\
	plot \
	"TheWay.txt" using 1:2 with points title "u"\
	,"TheWay.txt" using 1:3 with points title "u prime"\
	' | tee tmp.gpi | pyxplot
    

TheWaySIR.pyxplot.png: TheWaySIR.txt main Makefile out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "Time (days)";\
	set ylabel "People";\
	set title "SIR";\
	plot \
	"TheWaySIR.txt" using 1:2 with points title "Susceptible"\
	,"TheWaySIR.txt" using 1:3 with points title "Infectious"\
	,"TheWaySIR.txt" using 1:4 with points title "Recovered or dead"\
	' | tee tmp.gpi | pyxplot

TheWayNewton.pyxplot.png: TheWayNewton.txt main Makefile out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "Time";\
	set ylabel "Distance";\
	set title "SIR";\
	plot \
	"TheWayNewton.txt" using 1:2 with points title "Susceptible"\
	,"TheWayNewton.txt" using 1:3 with points title "Infectious"\
	,"TheWayNewton.txt" using 1:4 with points title "Recovered or dead"\
	' | tee tmp.gpi | pyxplot


out.txt: main 
	./$< > $@
    
main: functions.o main.c 
    
TheWay.txt: main    

TheWaySIR.txt: main    
    
clean: 
	$(RM) main main.o functions.o *.txt *.png *.gpi
    
    
    