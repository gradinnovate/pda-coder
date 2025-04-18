export TMPDIR := $(CURDIR)
CC = gcc
CFLAGS = -O3 -std=c11 -Wall -g

all: hyperml

hyperml: hyperml.c hyperml.h
	$(CC) $(CFLAGS) hyperml.c -o hyperml

clean:
	rm -f hyperml