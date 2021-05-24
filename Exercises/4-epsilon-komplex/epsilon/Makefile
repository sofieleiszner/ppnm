CC = gcc 
CFLAGS = -Wall -fwrapv 
LDLIBS = -lm

out.txt: main 
	./$< > $@

main: equal.o main.o
    
clean: 
	$(RM) main main.o equal.o  