# ALATDYN -- build everything into DESTDIR (make.inc: ~/ALADYN/bin)
#
#   make  (= make all)    FOCEX + utilities, SC-SNAPS, THERMACOND, ANFOMOD, links to the Python tools
#   make focex            FOCEX (focex-<ver>.x) and its utilities (read_outcar.x, read_qe.x, ...)
#   make sc_snaps         SC-SNAPS (sc_snaps.x, poscar2xyz.py)
#   make thermacond       THERMACOND with MPI (tcond-<ver>.x); make thermacond-serial: no MPI
#   make anfomod          ANFOMOD
#   make python-tools     links focex-agent, focex-diagnose to FOCEX/AGENT
#   make scop8            SCOP8 (not in "all": it does not compile at present)
#   make clean            remove the object files of every part
# Another place: make DESTDIR=/some/where.  THERMACOND needs mpif90 (on macOS with Homebrew
# Open MPI: export OMPI_FC=/opt/homebrew/bin/gfortran SDKROOT=$(xcrun --show-sdk-path)).

include make.inc
$(shell mkdir -p $(DESTDIR))
# SC_SNAPS's Makefile quotes its BINDIR, where the shell does not expand ~
BINDIR_ABS := $(patsubst ~%,$(HOME)%,$(DESTDIR))

.PHONY: all focex sc_snaps thermacond thermacond-serial anfomod scop8 python-tools clean

all: python-tools focex sc_snaps thermacond anfomod

focex:
	$(MAKE) -C FOCEX DESTDIR=$(DESTDIR)
	$(MAKE) -C FOCEX/utility DESTDIR=$(DESTDIR)

# SC_SNAPS is a git subtree of github.com/KeivanS/SC_SNAPS; its default target starts the GUI
sc_snaps:
	$(MAKE) -C SC_SNAPS compile BINDIR=$(BINDIR_ABS)

thermacond:
	$(MAKE) -C THERMACOND DESTDIR=$(DESTDIR)

thermacond-serial:
	$(MAKE) -C THERMACOND serial DESTDIR=$(DESTDIR)

anfomod:
	$(MAKE) -C ANFOMOD DESTDIR=$(DESTDIR)

scop8:
	$(MAKE) -C SCOP8

# the Python tools stay in the repository (they find their templates next to themselves);
# DESTDIR gets symlinks to them. Run them with the Python that has numpy, spglib, pymatgen.
python-tools:
	ln -sf $(CURDIR)/FOCEX/AGENT/focex-agent.py $(DESTDIR)/focex-agent
	ln -sf $(CURDIR)/FOCEX/AGENT/diagnose.py $(DESTDIR)/focex-diagnose

# the parts' clean rules fail when there is nothing to remove; carry on regardless
clean:
	-$(MAKE) -C ANFOMOD clean
	-$(MAKE) -C FOCEX clean
	-$(MAKE) -C FOCEX/utility clean
	-$(MAKE) -C SC_SNAPS clean
	-$(MAKE) -C SCOP8 clean
	-$(MAKE) -C THERMACOND clean
