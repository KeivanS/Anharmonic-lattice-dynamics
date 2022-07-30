FF=gfortran
DFLAGS=  #-C # to check everything O3 #g -p #for profiling with gprof
FLAGS= -fcheck=all -fdec-math #-Wall# -fbounds-check -fbacktrace #-ffpe-trap=invalid,zero,overflow,underflow,inexact #g -inline-debug-info #-p #for profiling with gprof -pedantic -std=f95 

MOD=modules_fcs.f90  
MAIN1=fc234.f90
FIT2=$(MOD) force_constants.f90 others3_nshells.f90 svd.f90 extratools.f90 ios.f90 

# combine many of them into one name "ALL"
ALL=$(MAIN1) $(FIT2)   

# declare the list of object files associated with the .f90 files
MODS=$(addsuffix .o, $(basename $(MOD)))
FIT2S=$(addsuffix .o, $(basename $(FIT2)))


# this is how to attach today's date to the backup .tar file
version=14
A=`date "+%F-%k_%M"`
C=fit-$(version).`date "+%F|%T"`
MYBIN=/Users/keivan/BIN/


# the use of the source file itself in the dependencies (1st line) is also necessary
fc234-$(version): $(FIT2S)  $(MAIN1)
	$(FF) $(FLAGS) $(FIT2S) $(MAIN1) -o $@ 
#	mv fc234-$(version) $(MYBIN)


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

arcfc234:
	rm -rf ARCH_FC234 arch-fc234.tar.gz arch-fc234.tar
	mkdir ARCH_FC234
	cp $(FIT2) $(MAIN1) Makefile params.inp ARCH_FC234
	tar -cvf arch-fc234.tar ARCH_FC234
	gzip arch-fc234.tar

# $@ is the name of the file to be made
# $? is the names of the changed dependents 
# $< the name of the related file that caused the action
# $* the prefix shared by target and dependent files
