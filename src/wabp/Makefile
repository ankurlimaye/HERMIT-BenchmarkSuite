
CC = gcc

CFILES = $(wildcard *.c)
INCFILES = $(wildcard include/*.c) $(wildcard include/*.h)
XFILES = $(CFILES:.c=)

CFLAGS = -static
LDFLAGS = -lm

$(XFILES):
	$(CC) $(CFILES) $(CFLAGS) -o $@ $(LDFLAGS)

all: $(XFILES)

.PHONY: clean
clean:
	rm -f $(XFILES) *.wabp