CC = gcc 
CFLAGS = -std=gnu11 -Wall
LDLIBS = -lm

out.txt: epsilon 
	./epsilon > out.txt

epsilon:  epsilon.o equal.o
	$(CC) -o epsilon epsilon.o $(LDLIBS)

epsilon.o: epsilon.c 
	$(CC) $(CFLAGS) -c epsilon.c 

equal.o: equal.c 
	$(CC) $(CFLAGS) -c $< -o $@


clean: 
	$(RM) epsilon.o epsilon out.txt
