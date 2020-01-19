#
# Makefile for gcbayes
#
# Created by Jayanth Chennamangalam
#

# C compilers and flags
CC = gcc
ifeq ($(shell uname), Darwin)
CFLAGS_INC_PGPLOT =-I/opt/local/include# define if needed (as -I[...])
else
CFLAGS_INC_PGPLOT =# define if needed (as -I[...])
endif
ifeq ($(OPT_DEBUG), yes)
CFLAGS = -g -pedantic -Wall -std=gnu99 $(CFLAGS_INC_PGPLOT)
else
CFLAGS = -O3 -pedantic -Wall -std=gnu99 $(CFLAGS_INC_PGPLOT)
endif

# linker flags
ifeq ($(shell uname), Darwin)
LFLAGS_PGPLOT_DIR =-L/opt/local/lib# define if not in $PATH (as -L[...])
else
LFLAGS_PGPLOT_DIR =# define if not in $PATH (as -L[...])
endif
LFLAGS = -lm -lpgplot -lcpgplot $(LFLAGS_PGPLOT_DIR)

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

