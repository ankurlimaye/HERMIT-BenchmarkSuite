CC = gcc

sd= ../sourcedir


EXE   = Poisson

SOURCES  = poisson.c \
           $(sd)/ran1.c \
           $(sd)/gammln.c \
           $(sd)/poidev.c \
           $(sd)/misc.c \
           $(sd)/egif_lib.c \
           $(sd)/dgif_lib.c \
           $(sd)/gif_err.c \
           $(sd)/gif_hash.c \
           $(sd)/eval.c \
           $(sd)/imgtools.c

IDIR      = -I../include

OBJS        = $(SOURCES:.c=.o)
CFLAGS     = -O2 -finline-functions -Winline -Wall
#CFLAGS      = -g
LFLAGS      = -lm


$(EXE) : $(OBJS)
	$(CC) $(CFLAGS) -o $(EXE) $(OBJS) $(LFLAGS)

%.o : %.c
	$(CC) -c $(IDIR) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJS) $(EXE)
