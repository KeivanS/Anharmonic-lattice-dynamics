PROG = scop8

SRCS =	check.f90 ConjugateGradient.f90 constants.f90 extratools.f90 \
	force_constants.f90 force_update.f90 kp_1d.f90 MatrixDiagonalize.f90 \
	mods9.f90 modules_tetra.f90 others3_nshells.f90 Broyden.f95 \
	DFT_force_constants.f95 Fourier_force_constants.f95 geometry.f95 \
	Iteration_parameters.f95 main.f95 Structure_info.f95 VA_math.f95 \
	zhegv.f mpi_params.f95

OBJS =	check.o ConjugateGradient.o constants.o extratools.o \
	force_constants.o force_update.o kp_1d.o MatrixDiagonalize.o mods9.o \
	modules_tetra.o others3_nshells.o Broyden.o DFT_force_constants.o \
	Fourier_force_constants.o geometry.o Iteration_parameters.o main.o \
	Structure_info.o VA_math.o zhegv.o mpi_params.o

LIBS =	

OPT = -O
FC = mpif90
FFLAGS = $(OPT) -ffpe-summary=invalid,zero,overflow 
F90 = mpif90
F90FLAGS = $(OPT) -ffpe-summary=invalid,zero,overflow 
LDFLAGS = 
all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)
	mv $(PROG) ~/ALADYN/bin/

.PHONY: clean
clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f .f90 .F90 .f95
.SUFFIXES: $(SUFFIXES) .c .cpp .cxx

.f90.o .f95.o .F90.o:
	$(F90) $(F90FLAGS) -c $<

.f.o:
	$(FC) $(FFLAGS) -c $<

check.o: Broyden.o force_update.o
ConjugateGradient.o: VA_math.o check.o
extratools.o: MatrixDiagonalize.o constants.o geometry.o mods9.o
force_constants.o: constants.o geometry.o
force_update.o: VA_math.o
kp_1d.o: Structure_info.o
MatrixDiagonalize.o: mods9.o
mods9.o: constants.o geometry.o 
modules_tetra.o: mods9.o
others3_nshells.o: DFT_force_constants.o constants.o geometry.o
Broyden.o: Iteration_parameters.o
DFT_force_constants.o: Structure_info.o mods9.o
Fourier_force_constants.o: DFT_force_constants.o kp_1d.o mods9.o
geometry.o: constants.o
Iteration_parameters.o: DFT_force_constants.o geometry.o mods9.o
main.o: Broyden.o ConjugateGradient.o VA_math.o check.o force_update.o mpi_params.o
Structure_info.o: constants.o geometry.o mods9.o
VA_math.o: Broyden.o Fourier_force_constants.o Iteration_parameters.o \
	MatrixDiagonalize.o constants.o geometry.o kp_1d.o mods9.o \
	modules_tetra.o mpi_params.o



# this is how to attach today's date to the backup .tar file
A=`date "+%F-%k_%M"`
C=scp.`date "+%F|%T"`.tar

# the use of the source file itself in the dependencies (1st line) is also necessary
oldscp: $(SUBS)  $(MAIN)
	$(FF) $(FLAGS) $(SUBS) $(MAIN)  -o $@
	mv oldscp ~/BIN

scp: $(SUBS)
	ar rcs libAll.a *.o
	$(FF) libAll.a $(MAIN) -o $@
	mv scp ~/ALADYN/bin/
	make clean

write:
	echo $(SUBS)
	echo $(ALL)

archive:
	rm -rf LASTBACKUP scp.tar.gz
	mkdir LASTBACKUP
	cp $(ALL) LASTBACKUP
	tar -cvf scp.tar LASTBACKUP
	mv scp.tar $(C)
	gzip $(C)

%.o : %.mod

# the above line is added because of a 'm2c' error
# $@ is the name of the file to be made
# $? is the names of the changed dependents
# $< the name of the related file that caused the action
# $* the prefix shared by target and dependent files

