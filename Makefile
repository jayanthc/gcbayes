#
# Makefile for gcbayes
#
# Created by Jayanth Chennamangalam
#

# C compilers and flags
CC = gcc
ifeq ($(OPT_DEBUG), yes)
CFLAGS = -g -pedantic -Wall -std=gnu99
else
CFLAGS = -O3 -pedantic -Wall -std=gnu99
endif

# linker flags
LFLAGS = -lm -lpgplot -lcpgplot

# directories
SRCDIR = ./src
IDIR = ./src
BINDIR = ./bin

# command definitions
DELCMD = rm

all: colourmap.o \
	 gcbayes_func.o \
	 gcbayes \
	 tags

colourmap.o: $(SRCDIR)/colourmap.c \
	         $(SRCDIR)/colourmap.h
	$(CC) -c $(CFLAGS) $< -o $(IDIR)/$@

gcbayes_func.o: $(SRCDIR)/gcbayes_func.c \
	            $(SRCDIR)/gcbayes.h
	$(CC) -c $(CFLAGS) $< -o $(IDIR)/$@

gcbayes: $(SRCDIR)/gcbayes.c \
         $(SRCDIR)/gcbayes.h \
         $(IDIR)/gcbayes_func.o \
         $(IDIR)/colourmap.o
	$(CC) $(CFLAGS) $(SRCDIR)/gcbayes.c $(IDIR)/gcbayes_func.o \
        $(IDIR)/colourmap.o $(LFLAGS) -o $(BINDIR)/$@

tags: $(SRCDIR)/*.c $(SRCDIR)/*.h
	ctags $(SRCDIR)/*.c $(SRCDIR)/*.h

install:
	@echo Nothing to install.

clean:
	$(DELCMD) $(IDIR)/gcbayes_func.o
	$(DELCMD) $(IDIR)/colourmap.o

