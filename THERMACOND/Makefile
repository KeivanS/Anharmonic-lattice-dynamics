FF=gfortran
ZFLAGS=   #-C # to check everything O3 #g -p #for profiling with gprof
FLAGS=-ffree-line-length-250 #-O3 -xHost  #-parallel # -pedantic -fbounds-check -fbacktrace #g -inline-debug-info #-p #for profiling with gprof

MOD=modules.f90  extratools.f90 force_constants.f90 
SUBKAP=mods2.f90 tetrahedron.f90  rest.f90   sub.f90  zhegv.f svd.f90 
KAPA=kappa.f90

# combine many of them into one name "ALL"
ALL=$(SUBKAP) $(KAPA) $(MOD) 

# declare the list of object files associated with the .f90 files
MODS=$(addsuffix .o, $(basename $(MOD)))
SUBKAPS=$(addsuffix .o, $(basename $(SUBKAP)))

OBJ2=$(MODS) $(SUBKAPS) #force_constants8.o 

# this is how to attach today's date to the backup .tar file
A=`date "+%F-%k_%M"`
C=kap.`date "+%F|%T"`.tar

# the use of the source file itself in the dependencies (1st line) is also necessary

kaptet: $(OBJ2) $(KAPA) 
	$(FF) $(FLAGS) $(OBJ2) $(KAPA) -o $@ 
#	mv kaptet ../


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

# explicit rules to make objects if there are some dependencies 
force_constants8.o: force_constants8.f modules.o
	$(FF) $(FLAGS) -c force_constants8.f

extratools.o: extratools.f90 modules.o mods2.o
	$(FF) $(FLAGS) -c extratools.f90

zhegv.o: zhegv.f
	$(FF) $(ZFLAGS) -c zhegv.f 

# tests , cleanup and archiving
write:
	echo $(ALL)
clean:
	rm  *.o *.mod 
archive:
	rm -rf LASTBACKUP kap.tar.gz
	mkdir LASTBACKUP
	cp $(ALL) Makefile LASTBACKUP  ##params.inp params.phon log.dat 
	tar -cvf kap.tar LASTBACKUP
	mv kap.tar $(C)
	gzip $(C)


# $@ is the name of the file to be made
# $? is the names of the changed dependents 
# $< the name of the related file that caused the action
# $* the prefix shared by target and dependent files
