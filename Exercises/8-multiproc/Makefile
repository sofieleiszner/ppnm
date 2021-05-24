CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += -lpthread
LDLIBS = -lm 
CC = gcc 


out.txt: main_thread
	./$< > $@

main_thread: main_thread.o 
	$(CC) -o main_thread main_thread.o $(LDLIBS) -lpthread

main_thread.o: main_thread.c
	$(CC) $(CFLAGS) -c main_thread.c -lpthread
    
#Det virker at skrive: 
#gcc main_thread.c -o masin_thread -lpthread -lm
#Og derefter køre programmet
