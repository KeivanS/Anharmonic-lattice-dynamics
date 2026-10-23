include make.inc
$(shell mkdir -p $(DESTDIR))

SUBDIRS := ANFOMOD FOCEX FOCEX/utility THERMACOND SCOP8   # SCOP8 last: it does not compile now

.PHONY: all $(SUBDIRS) clean python-tools

all: python-tools $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) -C $@ DESTDIR=$(DESTDIR)

# the Python tools stay in the repository (they find their templates next to themselves);
# DESTDIR gets symlinks to them. Run them with the Python that has numpy, spglib, pymatgen.
python-tools:
	ln -sf $(CURDIR)/FOCEX/AGENT/focex-agent.py $(DESTDIR)/focex-agent
	ln -sf $(CURDIR)/FOCEX/AGENT/diagnose.py $(DESTDIR)/focex-diagnose

clean:
	$(MAKE) -C ANFOMOD clean
	$(MAKE) -C FOCEX clean
	$(MAKE) -C FOCEX/utility clean
	$(MAKE) -C SCOP8 clean
	$(MAKE) -C THERMACOND clean
