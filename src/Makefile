# compiler choice
FF = gfortran
FLAG = -O3

# equation solver
LINEQ = 4

# build directory
BDIR = .

# optional libraries
LIBS =

# Architecture for the static build
MACHINE := $(shell uname -m)

# Self-contained ("make static") build
STATICLIBS = -static-libgfortran -static-libgcc -static-libquadmath

# Oldest macOS the static binary will load on
MINOS = -mmacosx-version-min=11.0

# Module files written by the compiler.  They are listed explicitly because the code also
# writes its LSD model spectrum to <prefix>.mod, which must not be removed by "make clean".
MODS = lsddata.mod lsdutils.mod sparse.mod lsdcalc.mod

# object files of the static build, kept apart from those of the ordinary one
SOBJ = ilsd_static.o lineq$(LINEQ)_static.o

.PHONY: all static clean

all: $(BDIR)/ilsd

$(BDIR)/ilsd: ilsd.o lineq$(LINEQ).o Makefile
	$(FF) $(FLAG) -o $(BDIR)/ilsd ilsd.o lineq$(LINEQ).o $(LIBS)

ilsd.o $(MODS): ilsd.f90 Makefile
	$(FF) -c $(FLAG) ilsd.f90

lineq$(LINEQ).o: lineq$(LINEQ).f Makefile
	$(FF) -c $(FLAG) lineq$(LINEQ).f

# Self-contained binary that runs without a Fortran installation
static: ilsd.f90 lineq$(LINEQ).f Makefile
	$(FF) -c $(FLAG) $(MINOS) ilsd.f90 -o ilsd_static.o
	$(FF) -c $(FLAG) $(MINOS) lineq$(LINEQ).f -o lineq$(LINEQ)_static.o
	$(FF) $(FLAG) $(MINOS) -o $(BDIR)/ilsd_$(MACHINE) $(SOBJ) $(STATICLIBS) $(LIBS)
	rm -f $(SOBJ)

clean:
	rm -f ilsd.o lineq$(LINEQ).o $(SOBJ) $(MODS) $(BDIR)/ilsd $(BDIR)/ilsd_$(MACHINE)
