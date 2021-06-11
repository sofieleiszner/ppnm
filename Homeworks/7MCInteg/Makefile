CFLAGS = -Wall -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

default: Error_Integral1.png Error_Integral2.png

Error_Integral1.png: main Makefile A.txt B.txt out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "N";\
	set ylabel "Error";\
	set title "Error comparison for integral from 0 to 1 of sqrt(x)";\
	a=1; b=-0.5; f(x)=a*x**b ;\
	c=1; d=-0.5; g(x)=c*x**d ;\
	fit f(x) "A.txt" using 1:3 via a,b ;\
	fit g(x) "B.txt" using 1:3 via c,d ;\
	plot \
	"A.txt" using 1:3 with points title "Plain MC error "\
	,"B.txt" using 1:3 with points title "Quasi-random MC error"\
	, f(x) with lines title sprintf("%3.1f*n^{%3.1f}",a,b)\
	, g(x) with lines title sprintf("%3.1f*n^{%3.1f}",c,d)\
	' |gnuplot
	
Error_Integral2.png: main Makefile AA.txt BB.txt out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "N";\
	set ylabel "Error";\
	set title "Error Comparison for integral given in task A description ";\
	e=1; f=-0.5; i(x)=e*x**f ;\
	g=1; h=-0.2; j(x)=g*x**h ;\
	fit i(x) "AA.txt" using 1:3 via e,f ;\
	fit j(x) "BB.txt" using 1:3 via g,h ;\
	plot \
	"AA.txt" using 1:3 with points title "Plain MC error "\
	,"BB.txt" using 1:3 with points title "Quasi-random MC error"\
	, i(x) with lines title sprintf("%3.1f*n^{%3.1f}",e,f)\
	, j(x) with lines title sprintf("%3.1f*n^{%3.1f}",g,h)\
	' |gnuplot

A.txt B.txt AA.txt BB.txt:main

out.txt: main 
	./$< > $@
    
main: functions.c main.c 
    
clean: 
	$(RM) main main.o functions.o *.txt *.png *.gpi *.log