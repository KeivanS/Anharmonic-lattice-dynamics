include make.inc
$(shell mkdir -p $(DESTDIR))

SUBDIRS := ANFOMOD FOCEX SCOP8 THERMACOND

.PHONY: all $(SUBDIRS) clean

all: $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) -C $@

clean:
	$(MAKE) -C ANFOMOD clean
	$(MAKE) -C FOCEX clean
	$(MAKE) -C SCOP8 clean
	$(MAKE) -C THERMACOND clean
