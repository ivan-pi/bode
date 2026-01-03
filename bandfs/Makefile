FYPP = fypp
FC = gfortran
FFLAGS = -Wall -Og -g

# Directory Definitions
TEMPLATE_DIR = templates
GEN_DIR = .

# List of templates to be generated
TARGET_NAMES = bandf.f bands.f bndmv.f90

# Variable Expansion
TARGETS = $(addprefix $(GEN_DIR)/, $(TARGET_NAMES))

all: $(TARGETS)

# --- Pattern Rules ---

$(GEN_DIR)/%.f: $(TEMPLATE_DIR)/%.fpp
	$(FYPP) $< $@

$(GEN_DIR)/%.f90: $(TEMPLATE_DIR)/%.fpp
	$(FYPP) $< $@

# --- Build Rules ---

examples: banded_example bandfs_demo

banded_example: $(TARGETS) banded_example.f90
	$(FC) $(FFLAGS) -I. -o $@ $^ -llapack

bandfs_demo: $(TARGETS) bandfs_demo.f90
	$(FC) $(FFLAGS) -I. -o $@ $^

.PHONY: clean
clean:
	$(RM) $(TARGETS) banded_example bandfs_demo