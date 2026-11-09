# ALATDYN -- build everything into DESTDIR (make.inc: ~/ALADYN/bin)
#
#   make  (= make all)    FOCEX + utilities, SC-SNAPS, THERMACOND, ANFOMOD, links to the Python tools
#   make focex            FOCEX (focex-<ver>.x) and its utilities (read_outcar.x, read_qe.x, ...)
#   make sc_snaps         SC-SNAPS (sc_snaps.x, poscar2xyz.py)
#   make thermacond       THERMACOND with MPI (tcond-<ver>.x); make thermacond-serial: no MPI
#   make anfomod          ANFOMOD
#   make python-tools     links focex-agent, focex-diagnose (FOCEX/AGENT), bandplot and allplots (FOCEX/PLOTTING), aladyn-gui (GUI)
#   make scop8            SCOP8 (not in "all": it does not compile at present)
#   make clean            remove the object files of every part
# Another place: make DESTDIR=/some/where.  THERMACOND needs mpif90.

include make.inc

# macOS with Homebrew: mpif90 must call the gfortran that compiled Open MPI's mpi.mod (Homebrew's),
# not another one earlier on the PATH ("mpi.mod ... created by a different version of GNU Fortran").
# Exported to the sub-makes; values already set in the environment are kept.
# OMPI_FC and SDKROOT come from make.inc (configure); without them, this guess:
ifeq ($(shell uname -s),Darwin)
  ifneq ($(wildcard /opt/homebrew/bin/gfortran),)
    OMPI_FC ?= /opt/homebrew/bin/gfortran
    SDKROOT ?= $(shell xcrun --show-sdk-path 2>/dev/null)
  endif
endif
ifneq ($(OMPI_FC),)
  export OMPI_FC
endif
ifneq ($(SDKROOT),)
  export SDKROOT
endif
# the compilers of make.inc (configure), passed to every part
FF ?= gfortran
MPIFC ?= mpif90
$(shell mkdir -p $(DESTDIR))
# SC_SNAPS's Makefile quotes its BINDIR, where the shell does not expand ~
BINDIR_ABS := $(patsubst ~%,$(HOME)%,$(DESTDIR))

.PHONY: all focex sc_snaps thermacond thermacond-serial anfomod scop8 python-tools clean

all: python-tools focex sc_snaps thermacond thermacond-serial anfomod

focex:
	$(MAKE) -C FOCEX DESTDIR=$(DESTDIR) FF=$(FF)
	$(MAKE) -C FOCEX/utility DESTDIR=$(DESTDIR)

# SC_SNAPS is a git subtree of github.com/KeivanS/SC_SNAPS; its default target starts the GUI
sc_snaps:
	$(MAKE) -C SC_SNAPS compile BINDIR=$(BINDIR_ABS) FF=$(FF)

thermacond:
	$(MAKE) -C THERMACOND DESTDIR=$(DESTDIR) MPIFC=$(MPIFC)

thermacond-serial:
	$(MAKE) -C THERMACOND serial DESTDIR=$(DESTDIR) MPIFC=$(MPIFC)

anfomod:
	$(MAKE) -C ANFOMOD DESTDIR=$(DESTDIR)

scop8:
	$(MAKE) -C SCOP8

# the Python tools stay in the repository (they find their templates next to themselves). DESTDIR gets a
# launcher for each that starts it with PYTHON (make.local, written by ./configure: the Python with numpy,
# spglib, pymatgen), so they work whatever python3 comes first on PATH (a batch job's 'module load python'
# puts the module's own python3 first). Without PYTHON: symlinks, run by the python3 on PATH.
PYTHON ?=
PYTOOLS = FOCEX/AGENT/focex-agent.py:focex-agent FOCEX/AGENT/diagnose.py:focex-diagnose \
          FOCEX/PLOTTING/allplots:allplots FOCEX/PLOTTING/kappaplots:kappaplots \
          FOCEX/utility/shell-optimizer:shell-optimizer GUI/aladyn-gui.py:aladyn-gui
python-tools:
	mkdir -p $(DESTDIR)
	ln -sf $(CURDIR)/FOCEX/PLOTTING/bandplot $(DESTDIR)/bandplot
	@for t in $(PYTOOLS); do \
	    src=$(CURDIR)/$${t%%:*}; dst=$(DESTDIR)/$${t##*:}; rm -f "$$dst"; \
	    if [ -n "$(PYTHON)" ]; then \
	        printf '#!/bin/sh\n# launcher written by make python-tools: %s with PYTHON of make.local\nexec "%s" "%s" "$$@"\n' \
	            "$$src" "$(PYTHON)" "$$src" > "$$dst" && chmod 755 "$$dst" && echo "$$dst: launcher, $(PYTHON)"; \
	    else ln -s "$$src" "$$dst" && echo "$$dst -> $$src"; fi; \
	done

# the parts' clean rules fail when there is nothing to remove; carry on regardless
clean:
	-$(MAKE) -C ANFOMOD clean
	-$(MAKE) -C FOCEX clean
	-$(MAKE) -C FOCEX/utility clean
	-$(MAKE) -C SC_SNAPS clean
	-$(MAKE) -C SCOP8 clean
	-$(MAKE) -C THERMACOND clean
