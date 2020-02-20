
CC = gcc

CFILES = $(wildcard *.c)
INCFILES = $(wildcard include/*.c) $(wildcard include/*.h)
XFILES = $(CFILES:.c=)

CFLAGS = -static

$(XFILES):
	$(CC) $(CFILES) $(CFLAGS) -o $@

all: $(XFILES)

.PHONY: clean
clean:
	rm -f $(XFILES) *.qrs