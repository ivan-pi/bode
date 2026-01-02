# Makefile for the templated Fortran procedures bandf and bands
# Designed to be launched from the root project folder, e.g.
#
#   make -f bandfs/Makefile
#

FYPP=fypp

# Directory Definitions
SRCDIR = src
FPPDIR = bandfs

TARGETS = src/bandf.f src/bands.f

all: $(TARGETS)

# The Generalized Pattern Rule
$(SRCDIR)/%.f: $(FPPDIR)/%.fpp
	$(FYPP) $< $@

FC=gfortran
FFLAGS=-Wall

banded_example: $(TARGETS) bandfs/banded_example.f90
	$(FC) $(FFLAGS) -I$(CURDIR)/bandfs/ -o $@ $^ -llapack

bandfs_demo: $(TARGETS) bandfs/bandfs_demo.f90
	$(FC) $(FFLAGS) -I$(CURDIR)/bandfs/ -o $@ $^

.PHONY: clean
clean:
	$(RM) $(TARGETS)