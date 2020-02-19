
CC = gcc

CFILES = $(wildcard *.c)
INCFILES = $(wildcard include/*.c) $(wildcard include/*.h)
XFILES = $(CFILES:.c=)

$(XFILES):
	$(CC) $(CFILES) -o $@

all: $(XFILES)

.PHONY: clean
clean:
	rm -f $(XFILES) *.qrs