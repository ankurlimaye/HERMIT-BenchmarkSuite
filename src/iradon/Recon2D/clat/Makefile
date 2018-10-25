CC = gcc

OUTPUT   = clat

sources  = clat.c

OBJS       = $(sources:.c=.o)
# CFLAGS     = -g
LFLAGS     = -lm

$(OUTPUT): $(OBJS)
	$(CC) $(CFLAGS) -o $(OUTPUT) $(OBJS) $(LFLAGS)
