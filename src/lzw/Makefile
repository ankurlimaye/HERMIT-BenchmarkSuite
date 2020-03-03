
CC = gcc

CFILES = lzw.c
INCFILES = $(wildcard include/*.c) $(wildcard include/*.h)
XFILES = lzw

CFLAGS = -static

$(XFILES):
	$(CC) $(CFILES) $(CFLAGS) -o $@

all: $(XFILES)

.PHONY: clean
clean:
	rm -f $(XFILES) *.lzw

