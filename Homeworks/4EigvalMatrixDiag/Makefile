CFLAGS = -Wall -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

default: eigenfunctions.png timediagonalization.png out.txt

eigenfunctions.png: main eigenfunctions.txt Makefile out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "N";\
	set ylabel "Time (s)";\
	set title "Our function";\
	plot \
	"eigenfunctions.txt" using 1:2 with line title "n=1 eigenfunction"\
	,"eigenfunctions.txt" using 1:3 with line title "n=2 eigenfunction"\
	,"eigenfunctions.txt" using 1:4 with line title "n=3 eigenfunction"\
	,"exacteigenfunctions.txt" using 1:2 with line title "n=1 exact eigenfunction"\
	,"exacteigenfunctions.txt" using 1:3 with line title "n=2 exact eigenfunction"\
	,"exacteigenfunctions.txt" using 1:4 with line title "n=3 exact eigenfunction"\
	' | pyxplot

timediagonalization.png: out_time_jacobi.txt out_time_gsl.txt out_time_jacobiop.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "N (Matrix of size NxN)";\
	set ylabel "Time (s)";\
	set title "Performance of our function compared to GSL: Jacobi diagonalization";\
	a=0.1; b=5; f(x)=a*x**b ;\
	c=0.1; d=3; g(x)=c*x**d ;\
	e=0.1; f=5; h(x)=e*x**f ;\
    fit f(x) "out_time_jacobi.txt" using 1:2 via a,b ;\
    fit g(x) "out_time_gsl.txt" using 1:2 via c,d ;\
    fit h(x) "out_time_jacobiop.txt" using 1:2 via e,f ;\
    plot \
	"out_time_jacobi.txt" using 1:2 with points title "Our implementation"\
    ,"out_time_gsl.txt" using 1:2 with points title "GSL"\
	,"out_time_jacobiop.txt" using 1:2 with points title "Our optimized implementation"\
    ,f(x) with lines title sprintf("%3.1f*10^-^7*n^{%3.1f}",a*10000000,b)\
    ,g(x) with lines title sprintf("%3.1f*10^-^1^1*n^{%3.1f}",c*100000000000,d)\
    ,h(x) with lines title sprintf("%3.1f*10^-^7*n^{%3.1f}",e*10000000,f)\
    '| gnuplot

eigenfunctions.txt: main

out.txt: main 
	./$< > $@

main: main.o functions.o

clean:
	$(RM) main *.o *.png *.txt *.gpi