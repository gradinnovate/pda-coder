CC = gcc
CFLAGS = -O2 -std=c99 -Wall -Wextra
TARGET = hyperfm
OBJS = hyperfm.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS)

hyperfm.o: hyperfm.c
	$(CC) $(CFLAGS) -c hyperfm.c

clean:
	rm -f $(TARGET) $(OBJS)
