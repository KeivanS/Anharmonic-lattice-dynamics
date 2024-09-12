include ../../make.inc

# List of source files
SOURCES = $(wildcard *.f90)

# combine many of them into one name "ALL" ; used for backup
ALL=$(SOURCES) cell.inp snaps.inp supercell.inp Makefile 

# List of executable files (replace .f90 with .x)
EXECUTABLES = $(SOURCES:.f90=.x)

# Default target: build all executables
all: $(EXECUTABLES)
	mv $(EXECUTABLES) $(DESTDIR)

# Rule to compile .f90 files into .x executables
%.x: %.f90
	$(FF) $(FFLAGS) $< -o $@

# tests , cleanup and archiving
write:
	echo $(ALL)

clean:
	rm -rf *.o *.mod
