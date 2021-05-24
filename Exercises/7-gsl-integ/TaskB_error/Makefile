CFLAGS = -Wall -Werror -O1
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
CC = gcc

.PHONEY:default
default:  erf.pyxplot.png 

erf.pyxplot.png: data.txt Makefile out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "err(x)";\
	set title "Error function";\
	plot \
	 "$<" using 1:2 with line title "Error Function using integral representation"\
     ' | pyxplot
     
     
out.txt: main
	./$< > $@

main: main.o function.o

data.txt: main

.PHONEY:clean
clean:
	$(RM) *.o main *.png out.txt
