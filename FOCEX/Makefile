FF=gfortran
DFLAGS=  #-C # to check everything O3 #g -p #for profiling with gprof
FLAGS= -fcheck=all -fdec-math #-Wall# -fbounds-check -fbacktrace #-ffpe-trap=invalid,zero,overflow,underflow,inexact #g -inline-debug-info #-p #for profiling with gprof -pedantic -std=f95 
FLAGS= -fcheck=all -fdec-math#-Wall# -fbounds-check -fbacktrace #-ffpe-trap=invalid,zero,overflow,underflow,inexact #g -inline-debug-info #-p #for profiling with gprof -pedantic -std=f95 

MOD=modules9.f90  svd.f90 extratools.f90 ios.f90
MOD9=modules9.f90 mods9.f90 modules_tetra.f90
SUBKAP=zhegv.f  force_constants.f90 extratools.f90 svd.f90
MAIN3=fc234.f90
FIT2=$(MOD) force_constants.f90 others3_nshells.f90 

# combine many of them into one name "ALL"
ALL=$(SUBKAP) $(MAIN3) $(FIT2)  modules9.f90 

# declare the list of object files associated with the .f90 files
MODS=$(addsuffix .o, $(basename $(MOD)))
MOD9S=$(addsuffix .o, $(basename $(MOD9)))
SUBKAPS=$(addsuffix .o, $(basename $(SUBKAP)))
FIT2S=$(addsuffix .o, $(basename $(FIT2)))

#OBJ2=modules.o force_constants8.o svd.o extratools.o mods2.o zhegv.o rest_bandsort.o
#OBJ9=modules9.o mods9.o modules_tetra.o force_constants9.o extratools.o zhegv.o subs.o

# this is how to attach today's date to the backup .tar file
version=13
A=`date "+%F-%k_%M"`
C=fit-$(version).`date "+%F|%T"`
MYBIN=/Users/keivan/BIN/


# the use of the source file itself in the dependencies (1st line) is also necessary
fc234-$(version): $(FIT2S)  $(MAIN3)
	$(FF) $(FLAGS) $(FIT2S) $(MAIN3) -o $@ 
#	mv fc234-$(version) $(MYBIN)

fcborn: $(FIT2S)  $(BORN)
	$(FF) $(FLAGS) $(FIT2S) $(BORN) -o $@ 
	mv fcborn $(MYBIN)

thermsplit: $(MOD9S) $(SUBKAPS) $(KAPA9) 
	$(FF) $(FLAGS) $(MOD9S) $(SUBKAPS) $(KAPA9) -o $@ 
	mv thermsplit $(MYBIN)

thermsplit2: $(OBJ9) $(KAPA9) 
	$(FF) $(FLAGS) $(OBJ9) $(KAPA9) -o $@ 
	mv thermsplit2 $(MYBIN)

test: $(OBJ9) test_tetrahedron.f90
	$(FF) $(FLAGS) $(OBJ9) test_tetrahedron.f90  -o $@ 

coll: collect.f90
	$(FF) collect.f90 -o $@
	mv coll $(MYBIN)

mfp: mfpthermal.f90
	$(FF) mfpthermal.f90 -o $@
	mv mfp $(MYBIN)

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
#	rm -rf LASTBACKUP all.tar.gz

arcfc234:
	rm -rf ARCH_FC234 arch-fc234.tar.gz arch-fc234.tar
	mkdir ARCH_FC234
	cp $(FIT2) $(MAIN3) Makefile params.inp ARCH_FC234
	tar -cvf arch-fc234.tar ARCH_FC234
	gzip arch-fc234.tar

# $@ is the name of the file to be made
# $? is the names of the changed dependents 
# $< the name of the related file that caused the action
# $* the prefix shared by target and dependent files
