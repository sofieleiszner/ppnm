CFLAGS = -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

default: TaskA.png TaskB.png
TaskA.png: data.txt main Makefile out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Neural network interpolates cos(5*x-1)*exp(-abs(x))";\
	plot \
	"data.txt" using 1:2 with points title "Exact y"\
	,"data.txt" using 1:3 with points title "Estimated y"\
	' | tee tmp.gpi | pyxplot

TaskB.png: data.txt main Makefile out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Neural network interpolates cos(5*x-1)*exp(-abs(x))";\
	plot \
	"dataB.txt" using 1:2 with lines title "Derivative"\
	,"dataB.txt" using 1:3 with lines title "Estimated y"\
	,"dataB.txt" using 1:4 with lines title "Anti-derivative"\
	,"dataB.txt" using 1:5 with lines title "Exact y"\
	' | tee tmp.gpi | pyxplot

out.txt: main 
	./$< > $@ 2>log
    
main: main.o functions.o lineqfunctions.o ann.o
    
clean: 
	$(RM) main *.o out.txt *.png *.gpi log