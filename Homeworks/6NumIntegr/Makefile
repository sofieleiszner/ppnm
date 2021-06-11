CFLAGS = -Wall -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

#default: TheWay.pyxplot.png TheWaySIR.pyxplot.png TheWayNewton.pyxplot.png Figure8_Newton.pyxplot.png

out.txt: main 
	./$< > $@
    
main: functions.c main.c 
    
clean: 
	$(RM) main main.o functions.o *.txt *.png *.gpi