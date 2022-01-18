FF=gfortran
FLAGS= #-O3 #-C # to check everything O3 #g -p #for profiling with gprof
DFLAGS= #-fbounds-check -fbacktrace -ffpe-trap=ivalid,zero,overflow,underflow #-debug extended #g -inline-debug-info #-p #for profiling with gprof

#MOD=modules_tetra.f90 geometry.f95   
#SUB1=mods9.f90 constants.f90 geometry.f95 force_constants.f90 zhegv.f MatrixDiagonalize.f90 extratools.f90
#SUB2= Structure_info.f95 DFT_force_constants.f95 kp_1d.f90 others3_nshells.f90
#SUB3=Fourier_force_constants.f95 Iteration_parameters.f95 VA_math.f95 Broyden.f95
#SUB4= check.f90 ConjugateGradient.f90 force_update.f90

SUB1=constants.f90 geometry.f95 force_constants.f90 zhegv.f 
SUB2=mods9.f90 MatrixDiagonalize.f90 modules_tetra.f90 extratools.f90
SUB3=Structure_info.f95 DFT_force_constants.f95 kp_1d.f90 
SUB4=Fourier_force_constants.f95 Iteration_parameters.f95 VA_math.f95 Broyden.f95 force_update.f90 others3_nshells.f90 check.f90 ConjugateGradient.f90 
MAIN=main.f95

SUBALL= $(SUB1)  $(SUB2)  $(SUB3)  $(SUB4)
# combine many of them into one name "ALL"
#ALL=$(MOD) $(SUB1)  $(SUB2)  $(SUB3)  $(SUB4) $(MAIN)  params.inp Makefile
ALL=$(SUBALL) $(MAIN)  params.inp params.phon params.born kpbs.in q_mesh.in iteration_parameters.in Makefile

# declare the list of object files associated with the .f90 files
#MODS=$(addsuffix .o, $(basename $(MOD)))
SUBS=$(addsuffix .o, $(basename $(SUBALL)))
OBJS =$(MODS) $(SUBS) 

#gfortran -o test main.f95 libAll.a 
#rm *.o *.mod


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
	mv scp ~/BIN
	make clean

# rules to make the "simple" object files from sources
.SUFFIXES: .f90 .f95
.f95.o:
	$(FF) $(FLAGS) -c $<
.f95 :
	$(FF) $(FLAGS) $< -o $@
.f90.o:
	$(FF) $(FLAGS) -c $<
.f90 :
	$(FF) $(FLAGS) $< -o $@
.f.o:
	$(FF) $(FLAGS) -c $<
.f :
	$(FF) $(FLAGS) $< -o $@

# explicit rules to make objects if there are some dependencies 
#force_constants8.o: force_constants8.f modules.o
#	$(FF) $(FLAGS) -c force_constants8.f


# tests , cleanup and archiving
write:
	echo $(SUBS)
	echo $(ALL)
clean:
	rm  *.o *.mod 
archive:
	rm -rf LASTBACKUP scp.tar.gz
	mkdir LASTBACKUP
	cp $(ALL) LASTBACKUP
	tar -cvf scp.tar LASTBACKUP
	mv scp.tar $(C)
	gzip $(C)

# $@ is the name of the file to be made
# $? is the names of the changed dependents 
# $< the name of the related file that caused the action
# $* the prefix shared by target and dependent files
