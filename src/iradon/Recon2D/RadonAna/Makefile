CC = gcc

SOURCEDIR= ../sourcedir

EXE  = RadonAna

SOURCES  = $(SOURCEDIR)/radonana.c \
           $(SOURCEDIR)/ran1.c \
           $(SOURCEDIR)/gasdev.c \
           $(SOURCEDIR)/misc.c \
           $(SOURCEDIR)/eval.c \
           $(SOURCEDIR)/imgtools.c \
           $(SOURCEDIR)/calc.c \
           $(SOURCEDIR)/egif_lib.c \
           $(SOURCEDIR)/dgif_lib.c \
           $(SOURCEDIR)/gif_err.c \
           $(SOURCEDIR)/gif_hash.c

IDIR     = -I../include
      
OBJS     = $(SOURCES:.c=.o)

CFLAGS   = -O3 -finline-functions -Winline -Wall -Wstrict-prototypes
#CFLAGS      = -g

LFLAGS   = -lm 

all     :       $(EXE)

$(EXE): $(OBJS)
	$(CC) $(CFLAGS) -o $(EXE) $(OBJS) $(LFLAGS)

%.o : %.c
	$(CC) -c $(IDIR) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJS) $(EXE)
