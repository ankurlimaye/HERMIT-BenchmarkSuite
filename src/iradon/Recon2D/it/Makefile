CC = gcc

SOURCEDIR = ../sourcedir

EXE   = it

SOURCES  = $(SOURCEDIR)/it.c \
           $(SOURCEDIR)/em.c \
           $(SOURCEDIR)/art.c \
           $(SOURCEDIR)/cg.c \
           $(SOURCEDIR)/amatrix.c \
           $(SOURCEDIR)/itreadargs.c \
           $(SOURCEDIR)/sparse.c \
           $(SOURCEDIR)/misc.c \
           $(SOURCEDIR)/eval.c \
           $(SOURCEDIR)/imgtools.c \
           $(SOURCEDIR)/calc.c \
           $(SOURCEDIR)/egif_lib.c \
           $(SOURCEDIR)/dgif_lib.c \
           $(SOURCEDIR)/gif_err.c \
           $(SOURCEDIR)/gif_hash.c

IDIR      = -I../include

OBJS        = $(SOURCES:.c=.o)

CFLAGS     = -O3 -finline-functions -Winline -Wall -m486 -malign-loops=2 -malign-jumps=2 -malign-functions=2 -Wstrict-prototypes

#CFLAGS      = -g -Wstrict-prototypes -Winline -Wall

LFLAGS      = -lm


$(EXE) : $(OBJS)
	$(CC) $(CFLAGS) -o $(EXE) $(OBJS) $(LFLAGS)

%.o : %.c
	$(CC) -c $(IDIR) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJS) $(EXE)
