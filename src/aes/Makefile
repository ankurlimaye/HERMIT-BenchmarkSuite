
CC = gcc

CFILES = $(wildcard *.c)
INCFILES = $(wildcard include/*.c) $(wildcard include/*.h)
XFILES = aes

CFLAGS = -Wall -static
LDFLAGS = -lm

$(XFILES):
	$(CC) $(CFILES) $(CFLAGS) -o $@

all: $(XFILES)

.PHONY: clean
clean:
	rm -f $(XFILES)