CC ?= cc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -pedantic
LDLIBS = -lm

.PHONY: all clean example

all: percolate

percolate: src/percolate.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

example: all
	./percolate --er 1000 4 1 nbc inf 2 results.csv sequence.txt

clean:
	rm -f percolate results.csv sequence.txt
