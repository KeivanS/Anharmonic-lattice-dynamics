# THERMACOND
#   make            parallel build with mpif90 -> $(DESTDIR)/tcond-$(ver).x
#                   (run: tcond-2.0.x, or mpirun -np N tcond-2.0.x)
#   make serial     serial build, no MPI library, same compiler as mpif90 -> $(DESTDIR)/tcond-$(ver)-serial.x
#   make DEBUG=1    either of the above with -g -fbounds-check -fbacktrace instead of -O2
#   make clean
# Objects and .mod files go to obj_mpi/ or obj_serial/, not next to the sources.
# macOS with Homebrew Open MPI: if mpif90 picks an older gfortran than the one Open MPI was
# built with ("created by a different version of GNU Fortran"), set
#   export OMPI_FC=/opt/homebrew/bin/gfortran   and   export SDKROOT=$(xcrun --show-sdk-path)

ver = 2.0
DESTDIR ?= ~/ALADYN/bin
EXE = tcond-$(ver).x
EXE_SERIAL = tcond-$(ver)-serial.x

MPIFC ?= mpif90
# serial: the compiler mpif90 uses (Open MPI's --showme:command), so that serial and MPI builds give
# the same numbers; gfortran if that does not answer.  make serial FC=... overrides it.
FC    := $(or $(shell $(MPIFC) --showme:command 2>/dev/null | head -1),gfortran)
ifeq ($(DEBUG),1)
  OPT = -g -fbounds-check -fbacktrace
else
  OPT = -O2
endif
FLAGS = $(OPT) -ffree-line-length-0

# modules.f90 first (it holds the modules); kappa.f90 (the main program) last
SRC  = modules.f90 zhegv.f force_constants.f90 sub.f90 rest.f90 extratools.f90 svd.f90 \
       tetrahedron.f90 cg.f90 emissivity.f90
OBJS = $(addsuffix .o,$(basename $(SRC)))

$(EXE): $(addprefix obj_mpi/,$(OBJS) kappa.o)
	$(MPIFC) $(FLAGS) $^ -o $@
	mkdir -p $(DESTDIR)
	mv $@ $(DESTDIR)/

serial: $(EXE_SERIAL)
$(EXE_SERIAL): $(addprefix obj_serial/,mpi_stub.o $(OBJS) kappa.o)
	$(FC) $(FLAGS) $^ -o $@
	mkdir -p $(DESTDIR)
	mv $@ $(DESTDIR)/

# every source uses the modules of modules.f90; the main program uses all of them
$(addprefix obj_mpi/,$(filter-out modules.o,$(OBJS))):    obj_mpi/modules.o
$(addprefix obj_serial/,$(filter-out modules.o,$(OBJS))): obj_serial/modules.o
obj_mpi/kappa.o:    $(addprefix obj_mpi/,$(OBJS))
obj_serial/kappa.o: $(addprefix obj_serial/,mpi_stub.o $(OBJS))

obj_mpi/%.o: %.f90
	@mkdir -p obj_mpi
	$(MPIFC) $(FLAGS) -Jobj_mpi -c $< -o $@
obj_mpi/%.o: %.f
	@mkdir -p obj_mpi
	$(MPIFC) $(FLAGS) -Jobj_mpi -c $< -o $@
obj_serial/%.o: %.f90
	@mkdir -p obj_serial
	$(FC) $(FLAGS) -Jobj_serial -c $< -o $@
obj_serial/%.o: %.f
	@mkdir -p obj_serial
	$(FC) $(FLAGS) -Jobj_serial -c $< -o $@

clean:
	rm -rf obj_mpi obj_serial $(EXE) $(EXE_SERIAL)

# source archive with today's date
archive:
	tar -czf thermacond.`date "+%F"`.tar.gz Makefile *.f90 *.f

.PHONY: serial clean archive
