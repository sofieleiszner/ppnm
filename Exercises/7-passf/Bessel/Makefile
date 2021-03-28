CFLAGS = -Wall -Werror -O1
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
CC = gcc

.PHONEY: default
default: erf.pyxplot.png 

erf.pyxplot.png: out.txt Makefile 
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "J1(x)";\
	set title "My Bessel Function";\
	plot \
	 "$<" using 1:2 with line title "My Bessel Function"\
     ' | pyxplot
     
out.txt: main
	./$< > $@

.PHONEY:clean
clean:
	$(RM) *.o main *.png out.txt
