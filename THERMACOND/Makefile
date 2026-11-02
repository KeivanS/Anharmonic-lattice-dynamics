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

# modules shared with FOCEX (../common, one file per module) first, then modules.f90 (the
# THERMACOND modules); kappa.f90 (the main program) last
vpath %.f90 ../common
COMMON = constants.f90 params.f90 geometry.f90 ios_mod.f90 born.f90
# shared sources that use THERMACOND's own modules (lattice, atoms_force_constants): after them
COMMON2 = tools.f90 ewald.f90 nonanal.f90 dielectric.f90 phonons.f90 gruneisen.f90
SRC  = $(COMMON) modules.f90 $(COMMON2) zhegv.f force_constants.f90 sub.f90 rest.f90 extratools.f90 svd.f90 \
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

# the common modules in their order of use, then modules.f90; every other source uses them all
COMMONO = $(COMMON:.f90=.o)
COMMON2O = $(COMMON2:.f90=.o)
NONCOMMON = $(filter-out $(COMMONO) $(COMMON2O),$(OBJS))
define chain
$(1)/params.o: $(1)/constants.o
$(1)/geometry.o: $(1)/params.o
$(1)/ios_mod.o: $(1)/geometry.o
$(1)/atoms_force_constants.o: $(1)/ios_mod.o
$(1)/lattice.o: $(1)/atoms_force_constants.o
$(1)/svd_stuff.o: $(1)/lattice.o
$(1)/linalgb.o: $(1)/svd_stuff.o
$(1)/born.o: $(1)/ios_mod.o
$(1)/modules.o: $(addprefix $(1)/,$(COMMONO))
$(addprefix $(1)/,$(COMMON2O)): $(1)/modules.o
$(1)/nonanal.o: $(1)/ewald.o
$(1)/gruneisen.o: $(1)/phonons.o
endef
$(eval $(call chain,obj_mpi))
$(eval $(call chain,obj_serial))
$(addprefix obj_mpi/,$(filter-out modules.o,$(NONCOMMON))):    obj_mpi/modules.o $(addprefix obj_mpi/,$(COMMON2O))
$(addprefix obj_serial/,$(filter-out modules.o,$(NONCOMMON))): obj_serial/modules.o $(addprefix obj_serial/,$(COMMON2O))
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
