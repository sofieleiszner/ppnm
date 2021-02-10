CC = gcc
CFLAGS = -std=gnu99 -O -Wall
LDLIBS = -lm 

out.txt: math
	./math > out.txt

math.o: math.c 
	$(CC) $(CFLAGS) -c math.c -o math.o 

math: math.o
	$(CC) -o math math.o $(LDLIBS)

clean: 
	$(RM) *.o out* math 
