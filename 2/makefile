
CC := gcc
TAUCC := tau_cc.sh

CFLAGS := -std=gnu99 -Wall -Wextra -Wpedantic -O3

SRCDIR := src
CFILES := $(filter-out ./$(SRCDIR)/main.c,$(wildcard ./$(SRCDIR)/*.c))
OFILES := $(CFILES:.c=.o)
TAUOFILES := $(CFILES:.c=.tau.o)

LIB := -lm


transpose: $(OFILES)
	$(CC) -o transpose.x ./$(SRCDIR)/main.c $(OFILES) $(CFLAGS) $(LIB)

$(SRCDIR)/%.o: $(SRCDIR)/%.c
	$(CC) -c -o $@ $< $(CFLAGS) $(LIB)


transpose-tau: $(TAUOFILES)
	$(TAUCC) -o transpose-tau.x ./$(SRCDIR)/main.c $(TAUOFILES) $(CFLAGS) $(LIB)

$(SRCDIR)/%.tau.o: $(SRCDIR)/%.c
	$(TAUCC) -c -o $@ $< $(CFLAGS) $(LIB)


clean:
	-rm *.x $(OFILES) $(TAUOFILES) main.o
