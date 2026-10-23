# SC-Snaps GUI — Makefile
#
BINDIR ?= $(HOME)/BIN

# Use the local virtual environment if one exists (see "make venv"), otherwise
# fall back to the system interpreter. Override with e.g. make run PYTHON=python3.12
ifneq ($(wildcard .venv/bin/python),)
PYTHON ?= .venv/bin/python
else
PYTHON ?= python3
endif

# bare "make" starts the GUI and opens it in the browser
.DEFAULT_GOAL := run

.PHONY: run compile venv clean help

run:
	$(PYTHON) sc-snaps-gui.py

compile:
	@mkdir -p "$(BINDIR)"
	gfortran -O2 -o sc_snaps.x sc_snaps.f90
	mv -f sc_snaps.x "$(BINDIR)/sc_snaps.x"
	cp -f poscar2xyz.py "$(BINDIR)/poscar2xyz.py"
	@echo "Compiled and moved sc_snaps.x to $(BINDIR)/sc_snaps.x"

# One-time setup of the Python dependencies. --system-site-packages keeps any
# numpy/matplotlib already provided by the system or by Homebrew, so only what
# is missing gets installed. Needed on Python installs marked EXTERNALLY-MANAGED
# (PEP 668, e.g. Homebrew), where a plain "pip install" is refused.
venv:
	python3 -m venv --system-site-packages .venv
	./.venv/bin/python -m pip install --upgrade pip
	./.venv/bin/python -m pip install flask numpy matplotlib scipy
	@echo "Virtual environment ready — 'make run' will now use it."

clean:
	@find . -name "*.pyc" -delete
	@find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
	@echo "Clean."

help:
	@echo "make         — same as 'make run' (default target)"
	@echo "make venv    — one-time: create .venv and install the Python packages"
	@echo "make compile — compile sc_snaps.f90 and move sc_snaps.x to $(BINDIR)"
	@echo "make run     — start SC-Snaps GUI on http://localhost:5050"
	@echo "make clean   — remove .pyc / __pycache__"
	@echo ""
	@echo "Current settings:"
	@echo "  PYTHON = $(PYTHON)"
	@echo "  BINDIR = $(BINDIR)"
