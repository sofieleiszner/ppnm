CC = gcc 
CFLAGS = -O -Wall -std=gnu11  
LDLIBS = -lm

    
all: output.txt
	cat $<

output.txt: read-cmdline read-stdin read-file input.txt
	./read-cmdline $$(cat input.txt)  > output.txt
	@echo                            >> output.txt
	./read-stdin   <      input.txt  >> output.txt
	@echo                            >> output.txt
	./read-file           input.txt     output.txt


input.txt:
	echo $$(seq 1 1 5) > input.txt



#Compiles the three files
read-cmdline: read-cmdline.o
read-cmdline.o: read-cmdline.c

read-stdin: read-stdin.o
read-stdin.o: read-stdin.c

read-file: read-file.o
read-file.o: read-file.c

    
clean: 
	$(RM) output.txt read-stdin read-cmdline read-file main 

