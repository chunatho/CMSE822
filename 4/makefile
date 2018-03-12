
CC := gcc

CFLAGS := -std=gnu99 -Wall -Wextra -Wpedantic -O3 -march=ivybridge -fopenmp

PI_SRCDIR := src/pi
PI_CFILES := $(filter-out ./$(PI_SRCDIR)/main.c,$(wildcard ./$(PI_SRCDIR)/*.c))
PI_OFILES := $(PI_CFILES:.c=.o)

CS_SRCDIR := src/countsort
CS_CFILES := $(filter-out ./$(CS_SRCDIR)/main.c,$(wildcard ./$(CS_SRCDIR)/*.c))
CS_OFILES := $(CS_CFILES:.c=.o)

LIB := -lm


all: pi countsort
	@echo "Compilation complete."


pi: $(PI_OFILES)
	$(CC) -o pi.x ./$(PI_SRCDIR)/main.c $(PI_OFILES) $(CFLAGS) $(LIB)


$(PI_SRCDIR)/%.o: $(PI_SRCDIR)/%.c
	$(CC) -c -o $@ $< $(CFLAGS) $(LIB)


countsort: $(CS_OFILES)
	$(CC) -o cs.x ./$(CS_SRCDIR)/main.c $(CS_OFILES) $(CFLAGS) $(LIB)

$(CS_SRCDIR)/%.o: $(CS_SRCDIR)/%.c
	$(CC) -c -o $@ $< $(CFLAGS) $(LIB)


clean:
	-rm *.x $(PI_OFILES) $(CS_OFILES)
