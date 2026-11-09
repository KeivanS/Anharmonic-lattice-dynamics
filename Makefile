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
# OMPI_FC and SDKROOT come from make.inc (configure.sh); without them, this guess:
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
# the compilers of make.inc (configure.sh), passed to every part
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

# the Python tools stay in the repository (they find their templates next to themselves);
# DESTDIR gets symlinks to them. Run them with the Python that has numpy, spglib, pymatgen.
python-tools:
	ln -sf $(CURDIR)/FOCEX/AGENT/focex-agent.py $(DESTDIR)/focex-agent
	ln -sf $(CURDIR)/FOCEX/AGENT/diagnose.py $(DESTDIR)/focex-diagnose
	ln -sf $(CURDIR)/FOCEX/PLOTTING/bandplot $(DESTDIR)/bandplot
	ln -sf $(CURDIR)/FOCEX/PLOTTING/allplots $(DESTDIR)/allplots
	ln -sf $(CURDIR)/FOCEX/utility/shell-optimizer $(DESTDIR)/shell-optimizer
	ln -sf $(CURDIR)/GUI/aladyn-gui.py $(DESTDIR)/aladyn-gui

# the parts' clean rules fail when there is nothing to remove; carry on regardless
clean:
	-$(MAKE) -C ANFOMOD clean
	-$(MAKE) -C FOCEX clean
	-$(MAKE) -C FOCEX/utility clean
	-$(MAKE) -C SC_SNAPS clean
	-$(MAKE) -C SCOP8 clean
	-$(MAKE) -C THERMACOND clean
