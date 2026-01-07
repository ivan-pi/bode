FYPP = fypp
FC = gfortran
FFLAGS = -Wall -Og -g

# Directory Definitions
TEMPLATE_DIR = templates
GEN_DIR = .
DIST_NAME = bandfs
DATE = $(shell date +%Y-%m-%d)

# Generated Files
TARGET_NAMES = bandf.f bands.f bndmv.f90 bndsl.f90
TARGETS = $(addprefix $(GEN_DIR)/, $(TARGET_NAMES))

# Example Programs
EXAMPLES = banded_example bandfs_demo bandfs_reaction_diffusion \
		   bandfs_fem_1d bandfs_fem_1d_quadratic

# --- Target-Specific Variables ---
# This appends -llapack only when building banded_example
banded_example: LIBS = -llapack

all: examples

.PHONY: generate
generate: $(TARGETS)

# --- Template Rules ---

$(GEN_DIR)/%.f: $(TEMPLATE_DIR)/%.fpp
	$(FYPP) $< $@

$(GEN_DIR)/%.f90: $(TEMPLATE_DIR)/%.fpp
	$(FYPP) $< $@

# --- Merged Build Rule for Examples ---

examples: $(EXAMPLES)

$(EXAMPLES): %: $(TARGETS) %.f90
	$(FC) $(FFLAGS) -I./include -o $@ $^ $(LIBS)

# --- Distribution Rule ---

dist: clean
	mkdir -p $(DIST_NAME)
	cp -r templates complex banded_example.f90 bandfs_demo.f90 bandfs.fi Makefile README.md $(DIST_NAME)/
	tar -czvf $(DIST_NAME)-$(DATE).tar.gz $(DIST_NAME)
	$(RM) -r $(DIST_NAME)

.PHONY: clean dist
clean:
	$(RM) $(EXAMPLES)
	$(RM) *.tar.gz