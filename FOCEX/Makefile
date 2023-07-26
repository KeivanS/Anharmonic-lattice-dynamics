FF=/usr/local/gfortran/bin/gfortran #~/LIB/forcomp #gfortran #nagfor #/opt/homebrew/bin/gfortran #ifort
DFLAGS= # -O3  #-C # to check everything O3 #g -p #for profiling with gprof
FLAGS=  #-fbounds-check # -fbacktrace #g -inline-debug-info #-p #for profiling with gprof

#SUB=modules9.f90 kpoints.f90 dyn_modules.f90 fourier_stuff.f90 ewald_born.f90 svd.f90 extratools.f90 ios.f90 others3_nshells.f90 new_amatrix.f90 force_constants_new.f90  dynamical.f90
SUB=fc_modules.f90 kpoints.f90 ios.f90 fourier_stuff.f90 ewald_born.f90 svd.f90 extratools.f90 others.f90 setup_amatrix.f90 force_constant_symmetries.f90
MAIN=fc234main.f90

# combine many of them into one name "ALL"
ALL=$(SUB) $(MAIN) structure.params dielectric.params Makefile POSCAR1 FORCEDISP1

# declare the list of object files associated with the .f90 files
SUBS=$(addsuffix .o, $(basename $(SUB)))


# this is how to attach today's date to the backup .tar file
A=`date "+%F-%k_%M"`
C=fcx`date "+%F_%T"`.tar
#MYBIN=/Users/keivan/BIN/
backup=FCX_$(A)
ver=5
exe=fcx_$(ver).x 


# the use of the source file itself in the dependencies (1st line) is also necessary
$(exe): $(SUBS)  $(MAIN)
	$(FF) $(FLAGS) $(SUBS) $(MAIN) -o $@ 
	echo 'executable was made on' $A
	mv $(exe)  ~/BIN

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
	echo $(A) $(C) $(backup) 
clean:
	rm  *.o *.mod 
archive:
#	make clean
#	rm -rf LASTBACKUP
	mkdir $(backup)
	cp $(ALL) $(backup)
	tar -cvf all.tar $(backup)
	gzip all.tar
	mv all.tar.gz $(backup)


# $@ is the name of the file to be made
# $? is the names of the changed dependents 
# $< the name of the related file that caused the action
# $* the prefix shared by target and dependent files
