CFLAGS = -Wall -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

out.txt: main 
	./$< > $@
    
main: functions.o main.o lineqfunctions.o
    
clean: 
	$(RM) main *.o out.txt *.png *.gpi