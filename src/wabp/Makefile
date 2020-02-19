
CC = gcc

CFILES = $(wildcard *.c)
INCFILES = $(wildcard include/*.c) $(wildcard include/*.h)
XFILES = $(CFILES:.c=)

LDFLAGS = -lm

$(XFILES):
	$(CC) $(CFILES) -o $@ $(LDFLAGS)

all: $(XFILES)

.PHONY: clean
clean:
	rm -f $(XFILES) *.wabp