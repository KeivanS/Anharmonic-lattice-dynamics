FF=gfortran #ifort
DFLAGS= #-O3  #-C # to check everything O3 #g -p #for profiling with gprof
FLAGS= #-fbounds-check -fbacktrace #g -inline-debug-info #-p #for profiling with gprof

MOD=modules.f90  svd.f90 extratools.f90
MOD9=modules9.f90 mods9.f90 modules_tetra.f90
SUBKAP=zhegv.f  force_constants.f90 subs.f90 cg.f90  extratools.f90 svd.f90
KAPA9=thermalsplit.f90
MAIN3=fc234.f90
FIT2=$(MOD) force_constants9.f others3.f90 

# combine many of them into one name "ALL"
ALL=$(SUBKAP) $(KAPA9)  mods9.f90 modules9.f90 modules_tetra.f90

# declare the list of object files associated with the .f90 files
MODS=$(addsuffix .o, $(basename $(MOD)))
MOD9S=$(addsuffix .o, $(basename $(MOD9)))
SUBKAPS=$(addsuffix .o, $(basename $(SUBKAP)))
FIT2S=$(addsuffix .o, $(basename $(FIT2)))

#OBJ2=modules.o force_constants8.o svd.o extratools.o mods2.o zhegv.o rest_bandsort.o
#OBJ9=modules9.o mods9.o modules_tetra.o force_constants9.o extratools.o zhegv.o subs.o

# this is how to attach today's date to the backup .tar file
version=10
A=`date "+%F-%k_%M"`
C=therm-$(version).`date "+%F|%T"`
MYBIN=/Users/keivan/BIN/


# the use of the source file itself in the dependencies (1st line) is also necessary
#fc234: $(FIT2S)  $(MAIN3)
#	$(FF) $(FLAGS) $(FIT2S) $(MAIN3) -o $@ 
#	mv fc234 $(MYBIN)

thermsplit-$(version): $(MOD9S) $(SUBKAPS) $(KAPA9) 
	$(FF) $(FLAGS) $(MOD9S) $(SUBKAPS) $(KAPA9) -o $@ 
#	mv thermsplit-$(version) $(MYBIN)

#thermsplit2: $(OBJ9) $(KAPA9) 
#	$(FF) $(FLAGS) $(OBJ9) $(KAPA9) -o $@ 
#	mv thermsplit2 $(MYBIN)

test: $(OBJ9) test_tetrahedron.f90
	$(FF) $(FLAGS) $(OBJ9) test_tetrahedron.f90  -o $@ 

coll: collect.f90
	$(FF) collect.f90 -o $@
	mv coll $(MYBIN)

mfpband: mfpthermal.f90
	$(FF) mfpthermal.f90 -o $@
	mv mfpband $(MYBIN)

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
	echo $(MYBIN)
write:
	echo $(ALL)
clean:
	rm  *.o *.mod 
archive:
	mkdir BACKUP-$(C)
	cp $(ALL) Makefile params.inp params.phon params.born kpbs.in BACKUP-$(C)
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
