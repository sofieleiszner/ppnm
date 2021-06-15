CFLAGS = -Wall -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

all: out.txt
	@cat  out.txt


out.txt: main
	./$< > $@

    
main: main.o functions.o

    
    
clean: 
	$(RM) main *.o out.txt
    
    
test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)
