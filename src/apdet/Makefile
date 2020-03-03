# Makefile		Joe Mietus		Dec 19 2001

CC = gcc
BINDIR = /usr/local/bin
ALL = av detruns filt ht htavsd htmedfilt ldetrend linsamp mm rrlist smooth

all :		$(ALL)

help:
	@echo "Type 'make all' to compile everything in this directory."
	@echo "Type 'make install' to install executables into" $(BINDIR)
	@echo "Type 'make clean' to remove binaries and backups."

install :	all
	-mkdir $(BINDIR)/apdet
	install -s $(ALL) $(BINDIR)/apdet
	install get_apdet $(BINDIR)

clean :
	rm -f $(ALL) *.o *~ core

av :		av.c
	$(CC) -O -o $@ av.c -lm

detruns :	detruns.c
	$(CC) -O -o $@ detruns.c

filt :		filt.c
	$(CC) -O -o $@ filt.c

ht :		ht.c
	$(CC) -O -o $@ ht.c -lm

htavsd :	htavsd.c
	$(CC) -O -o $@ htavsd.c -lm

htmedfilt :	htmedfilt.c
	$(CC) -O -o $@ htmedfilt.c

ldetrend :	ldetrend.c
	$(CC) -O -o $@ ldetrend.c

linsamp :	linsamp.c
	$(CC) -O -o $@ linsamp.c

mm :		mm.c
	$(CC) -O -o $@ mm.c

rrlist :	rrlist.c
	$(CC) -O -o $@ rrlist.c

smooth :	smooth.c
	$(CC) -O -o $@ smooth.c

