CFLAGS=-g -O2 -Wall -std=c99


PROGRAM_NAME= factor
OBJS = main.o 
CC = mpicc

$(PROGRAM_NAME): $(OBJS)
	$(CC) -o $@ $? -lgmp

clean:
	rm  *.o $(PROGRAM_NAME) *~
