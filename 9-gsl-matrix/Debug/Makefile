CFLAGS = -Wall -O1 -std=gnu11 -Werror
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS = $(shell /usr/bin/gsl-config --libs)

all: output.txt
	cat output.txt

output.txt: main
	./main > output.txt

clean: 
	$(RM) main output.txt 

