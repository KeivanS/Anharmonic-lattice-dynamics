# THERMACOND
#   make            parallel build with mpif90 -> kap8            (run: kap8, or mpirun -np N kap8)
#   make serial     serial build with gfortran and no MPI library -> kap8_serial
#   make DEBUG=1    either of the above with -g -fbounds-check -fbacktrace instead of -O2
#   make clean
# Objects and .mod files go to obj_mpi/ or obj_serial/, not next to the sources.
# macOS with Homebrew Open MPI: if mpif90 picks an older gfortran than the one Open MPI was
# built with ("created by a different version of GNU Fortran"), set
#   export OMPI_FC=/opt/homebrew/bin/gfortran   and   export SDKROOT=$(xcrun --show-sdk-path)

MPIFC ?= mpif90
FC     = gfortran
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

kap8: $(addprefix obj_mpi/,$(OBJS) kappa.o)
	$(MPIFC) $(FLAGS) $^ -o $@

serial: kap8_serial
kap8_serial: $(addprefix obj_serial/,mpi_stub.o $(OBJS) kappa.o)
	$(FC) $(FLAGS) $^ -o $@

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
	rm -rf obj_mpi obj_serial kap8 kap8_serial

# source archive with today's date
archive:
	tar -czf thermacond.`date "+%F"`.tar.gz Makefile *.f90 *.f

.PHONY: serial clean archive
