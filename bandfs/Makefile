# Makefile for the templated Fortran procedures bandf and bands
# Designed to be launched from the root project folder, e.g.
#
#   make -f bandfs/Makefile
#

FYPP=fypp

FFILES = src/bandf.f src/bands.f

all: $(FFILES)

src/bandf.f: bandfs/bandf.fpp
	$(FYPP) $< $@
src/bands.f: bandfs/bands.fpp
	$(FYPP) $< $@

FC=gfortran
FFLAGS=-Wall

banded_example: src/bandf.f src/bands.f bandfs/banded_example.f90
	$(FC) $(FFLAGS) -I$(CURDIR)/bandfs/ -o $@ $^ -llapack

.PHONY: clean
clean:
	$(RM) $(FFILES)