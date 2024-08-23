FF=mpif90  -fbounds-check -g #gfortran #ifort
FLAGS= -ffree-line-length-0 -g #-fdec-math #-fbounds-check -fbacktrace #g -inline-debug-info #-p #for profiling with gprof

MOD=modules.f90  
SUBKAP=zhegv.f force_constants.f90 sub.f90 rest.f90  extratools.f90 svd.f90 tetrahedron.f90 cg.f90
KAPA=kappa.f90

# combine many of them into one name "ALL"
ALL=$(SUBKAP) $(KAPA) structure.params dielectric.params latdyn.params kpbs.params 

# declare the list of object files associated with the .f90 files
MODS=$(addsuffix .o, $(basename $(MOD)))
SUBKAPS=$(addsuffix .o, $(basename $(SUBKAP)))


# this is how to attach today's date to the backup .tar file
version=8
A=`date "+%F-%k_%M"`
C=kap$(version).`date "+%F|%T"`
#MYBIN=/Users/ke4c/BIN/


kap$(version): $(MODS) $(SUBKAPS) $(KAPA) 
	$(FF) $(FLAGS) $(MODS) $(SUBKAPS) $(KAPA) -o $@ 
#	mv kap$(version) $(MYBIN)

mfpband: mfpthermal.f90
	$(FF) mfpthermal.f90 -o $@
#	mv mfpband $(MYBIN)

# rules to make the "simple" object files from sources
.SUFFIXES: .f90
.f90.o:
	$(FF) $(FLAGS) -c $<
.f90 :
	$(FF) $(FLAGS) $< -o $@
.f.o:
	$(FF) $(FLAGS) -c $<
.f :
	$(FF) $(FLAGS) $< -o $@

# tests , cleanup and archiving
ec:
#	echo $(MYBIN)
write:
	echo $(ALL)
clean:
	rm  *.o *.mod 
archive:
	mkdir BACKUP-$(C)
	cp $(ALL) Makefile BACKUP-$(C)
	tar -cvf all.tar BACKUP-$(C)
	mv all.tar $(C).tar
	gzip $(C).tar

tt:
	rm -rf TEST test.tar.gz test.tar
	mkdir TEST
	cp $(MOD) $(MOD2) $(SUB3) Makefile params.inp TEST
	tar -cvf test.tar TEST
	gzip test.tar

# $@ is the name of the file to be made
# $? is the names of the changed dependents 
# $< the name of the related file that caused the action
# $* the prefix shared by target and dependent files
