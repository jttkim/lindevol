# Makefile for the LindEvol 1.23 program

.PHONY : install all

PREFIX	= $(HOME)
BINDIR	= $(PREFIX)/bin
INCLUDEPATH = $(PREFIX)/include
LIBPATH = $(PREFIX)/lib
LND1LIB = liblnd1.a

ALL 	= lnd1v23 lnd1c23 lnd1f23 stplot1
CC 	= gcc
LD	= gcc
INCLUDE	= -I$(INCLUDEPATH)
STDLIBS	= -lm
LIBS 	= -llnd1 -lgenome -lgnlib -ljklib $(STDLIBS)
# OPT 	= -Wall -O2 -fomit-frame-pointer
OPT 	= -Wall -g
CFLAGS 	= $(OPT) $(DEFINE) $(INCLUDE)
LDFLAGS	= -L$(LIBPATH) -L.

include objects 

lnd1v23: lnd1v23.o $(LND1LIB)
	$(LD) $(LDFLAGS) lnd1v23.o -o lnd1v23 $(LIBS)

all : $(ALL)
	@echo 'all lnd1v23 programs updated'

install : $(ALL)
	cp $(ALL) $(BINDIR)

lnd1c23: lnd1c23.o $(LND1LIB)
	$(LD) $(LDFLAGS) lnd1c23.o -o lnd1c23 $(LIBS)

lnd1f23: lnd1f23.o $(LND1LIB)
	$(LD) $(LDFLAGS) lnd1f23.o -o lnd1f23 $(LIBS)

lnd1c23.o : lnd1v23.c
	$(CC) $(CFLAGS) -c -DRANDOM_FITNESS lnd1v23.c -o lnd1c23.o

lnd1f23.o : lnd1v23.c
	$(CC) $(CFLAGS) -c -DFUJI_FITNESS lnd1v23.c -o lnd1f23.o

stplot1 : stplot1.o $(LND1LIB)
	$(LD) $(LDFLAGS) stplot1.o -o stplot1 $(LIBS)

$(LND1LIB) : $(LIBOBJS)
	rm -rf $(LND1LIB)
	ar -rsc $(LND1LIB) $(LIBOBJS)

clean :
	rm -f lnd1v23 lnd1c23 lnd1f23 core $(OBJS) $(LIBOBJS) lnd1v23.o lnd1c23.o lnd1f23.o mdepend

mdepend :
	$(CC) $(INCLUDE) -MM *.c > mdepend

include mdepend

