CFLAGS = -Wall -Werror -O1
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
CC = gcc

out.txt: main
	./$< > $@

main: main.o functions.o