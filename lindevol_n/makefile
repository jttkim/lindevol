# Makefile for the LindEvol program

LIBPREFIX  = lib
LIBPOSTFIX = .a

BINDIR	= $(HOME)/bin

PRG 	= lnd6v01
PRGS	= $(PRG) stplot6 xlmovie
LNDLIB	= lndlib
LNDLIBNAME = $(LIBPREFIX)$(LNDLIB)$(LIBPOSTFIX)
CC 	= gcc
CCXX 	= g++
AR	= ar -rsc
STDLIBS	= -lXext -lX11 -lm
# OPT 	= -Wall -O6
OPT 	= -Wall -g
# OPT	= -fullwarn -O2 -pedantic

QTDIR	= /usr/local/qt
QTLIBS	= -lqt
# QTINCLUDE = -I/usr/include/qt
# QTLIBDIR  = /usr/lib
QTCINCLUDE = -I$(QTDIR)/include
QTLIBDIR = $(QTDIR)/lib

include objects

JKLIBS 	= -lgfaport -lgnlib -lptlib -lgenome -ljklib
JKINCL  = -I$(HOME)/include

INCLUDE = $(JKINCL) $(QTINCLUDE)
CFLAGS 	= $(OPT) $(DBGDEF) $(DEFINE) $(INCLUDE)
LDFLAGS	= $(OPT) -L. -L$(HOME)/lib -L$(QTLIBDIR) -L/usr/lib/X11
LIBS	= -llndlib $(JKLIBS) $(QTLIBS) $(STDLIBS)

$(PRG) : $(OBJS) $(LNDLIBNAME)
	$(CC) $(LDFLAGS) $(OBJS) -o $(PRG) $(LIBS)

.PHONY : all install

all : $(PRGS)
	@echo
	@echo '***********************************************************'
	@echo '* All programs of the lindevol project are now up to date *'
	@echo '***********************************************************'

install : $(PRG)
	cp $(PRG) $(BINDIR)

$(LNDLIBNAME) : $(LNDLIBOBJS)
	rm -f $(LNDLIBNAME)
	$(AR) $(LNDLIBNAME) $(LNDLIBOBJS)

stplot6: $(LNDLIBNAME) $(STPLOTOBJS)
	$(CC) $(LDFLAGS) $(STPLOTOBJS) -o stplot6 $(LIBS)

xlmovie : xlmovie.o
	$(CC) $(LDFLAGS) xlmovie.o -o xlmovie -lwhndlr -lXext -lX11 $(LIBS)

readtree : readtree.o
	$(CC) $(LDFLAGS) readtree.o -o readtree $(LIBS)

cmptree : cmptree.o
	$(CC) $(LDFLAGS) cmptree.o -o cmptree $(LIBS)

jfdiff: jfdiff.o
	$(CC) $(LDFLAGS) jfdiff.o -o jfdiff $(LIBS)

plottree: plottree.o
	$(CC) $(LDFLAGS) plottree.o -o plottree $(LIBS)

plotcons: plotcons.o
	$(CC) $(LDFLAGS) plotcons.o -o plotcons $(LIBS)

rndtrdst : rndtrdst.o
	$(CC) $(LDFLAGS) rndtrdst.o -o rndtrdst $(LIBS)

xmsf : xmsf.o $(LNDLIBNAME)
	$(CC) $(LDFLAGS) $(CLFAGS) xmsf.o -o xmsf $(LIBS)

counttrs : counttrs.o $(LNDLIBNAME)
	$(CC) $(LDFLAGS) $(CLFAGS) counttrs.o -o counttrs $(LIBS)

treechk : treechk.o $(LNDLIBNAME)
	$(CC) $(LDFFAGS) treechk.o -o treechk $(LIBS)

rtrchk : rtrchk.o $(LNDLIBNAME)
	$(CC) $(LDLFAGS) rtrchk.o -o rtrchk $(LIBS)

tree2dmt : tree2dmt.o $(LNDLIBNAME)
	$(CC) $(LDLFAGS) tree2dmt.o -o tree2dmt $(LIBS)

distcorr : distcorr.o $(LNDLIBNAME)
	$(CC) $(LDLFAGS) distcorr.o -o distcorr $(LIBS)

treecorr : treecorr.o $(LNDLIBNAME)
	$(CC) $(LDFAGS) treecorr.o -o treecorr $(LIBS)

clean :
	rm -f $(PRG) $(LNDLIBOBJS) $(OBJS) $(LNDLIBNAME) mdepend

mdepend :
	gcc $(INCLUDE) -MM *.c > mdepend
	g++ $(INCLUDE) -MM *.cxx > mdepend

include mdepend

%.o : %.cxx
	$(CCXX) -c $(CFLAGS) $*.cxx -o $*.o

mlndwdgt.cxx : lndwidgt.h
	moc -o mlndwdgt.cxx lndwidgt.h

mctrlwdg.cxx : ctrlwdgt.h
	moc -o mctrlwdg.cxx ctrlwdgt.h

