CC = gcc
CFLAGS = -std=c99 -O2 -Wall
TARGET = floorplan
OBJS = floorplan.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) -lm

floorplan.o: floorplan.c floorplan.h
	$(CC) $(CFLAGS) -c floorplan.c

clean:
	rm -f $(OBJS) $(TARGET)

gen-large-sample:
	bash gen_large_sample.sh > sample_large.txt

.PHONY: all clean gen-large-sample