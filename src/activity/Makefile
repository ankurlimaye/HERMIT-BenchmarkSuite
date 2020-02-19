
CC = gcc

CFILES = $(wildcard *.c)
XFILES = $(CFILES:.c=)

LDFLAGS = -lm

$(XFILES):
	$(CC) $(CFILES) -o $@ $(LDFLAGS)

all: $(XFILES)

.PHONY: clean
clean:
	rm -f $(XFILES) *.txt