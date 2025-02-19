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

.PHONY: clean
clean:
	$(RM) $(FFILES)