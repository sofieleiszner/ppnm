CFLAGS = -Wall -O1 -std=gnu11 -Werror
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS = $(shell /usr/bin/gsl-config --libs)

all: output ; cat $<
output: main ; ./$< > $@
clean: ; $(RM) main output.txt 

