CFLAGS=-Wall -Werror -O1
CFLAGS+=$(shell gsl-config --cflags)
LDLIBS+=$(shell gsl-config --libs)
CC=gcc

.PHONEY:default
default: fig-linterp-pyxplot.png fig-quadterp-pyxplot.png fig-cubicterp-pyxplot.png

fig-linterp-pyxplot.png: linterp.dat input.dat main out.txt Makefile
	echo '\
	set term png;\
	set width 8;\
	set xlabel "$x$";\
	set ylabel "$y$";\
	set fontsize .5;\
	set key right;\
	set output "$@";\
	plot \
	"linterp.dat" using 1:2 with lines title "Linear interpolation" \
	,"input.dat" using 1:2 with points title "original points" \
	,"linterp.dat" using 1:3 with lines title "linear interpolation integration" \
	,"linterp.dat" using 1:7 with point title "GSL linear interpolation" \
	,"linterp.dat" using 1:8 with point title "GSL linear integration" \
	' | tee tmp.gpi | pyxplot

fig-quadterp-pyxplot.png: linterp.dat input.dat main out.txt Makefile
	echo '\
	set terminal png;\
	set width 8;\
	set xlabel "$x$";\
	set ylabel "$y$";\
	set fontsize .5;\
	set key right;\
	set output "$@";\
	plot \
	"input.dat" using 1:2 with points title "original points" \
	,"linterp.dat" using 1:4 with lines title "quadratic interpolation" \
	,"linterp.dat" using 1:5 with lines title "quadratic integration" \
	,"linterp.dat" using 1:6 with lines title "quadratic differential" \
	' | tee tmp.gpi | pyxplot
	
fig-cubicterp-pyxplot.png: linterp.dat input.dat main out.txt Makefile
	echo '\
	set terminal png;\
	set width 8;\
	set xlabel "$x$";\
	set ylabel "$y$";\
	set fontsize .5;\
	set key right;\
	set output "$@";\
	plot \
	"input.dat" using 1:2 with points title "original points" \
	,"linterp.dat" using 1:9 with lines title "cubic interpolation" \
	,"linterp.dat" using 1:10 with lines title "cubic differential" \
	,"linterp.dat" using 1:11 with lines title "cubic integration" \
	,"linterp.dat" using 1:12 with point title "cubic GSL interpolation" \
	,"linterp.dat" using 1:13 with point title "cubic GSL integration" \
	,"linterp.dat" using 1:14 with point title "cubic GSL derivative" \
	' | tee tmp.gpi | pyxplot
	
out.txt: main 
	./$< > $@

main: main.o linterp.o linterp_integ.o binsearch.o cubic.o quad.o quad_integ.o quad_diff.o

linterp.dat: main

clean:
	$(RM) main *.o linterp.dat *.png *.gpi


