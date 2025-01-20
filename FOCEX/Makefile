include ../make.inc

SUBF90=fc_modules.f90 dyn_modules.f90 kpoints.f90 modules_tetra.f90 fourier_stuff.f90 ios.f90 ewald_born.f90 svd_ridge.f90 extratools.f90 others.f90 setup_amatrix.f90 force_constant_symmetries.f90 dynamical.f90 zhegv.f 
MAIN=focex.f90 
PH=ph.f90

# declare the list of object files associated with the .f90 files
SUBO=$(addsuffix .o, $(basename $(SUBF90)))
MAINO=$(addsuffix .o, $(basename $(MAIN)))
PHO=$(addsuffix .o, $(basename $(PH)))

# combine many of them into one name "ALL" ; used for backup
ALL=$(SUBF90) $(MAIN) $(PH)  structure.params dielectric.params default.params kpbs.params latdyn.params Makefile 


# this is how to attach today's date to the backup .tar file
# contains time up to minutes
M=`date "+%F-%k_%M"`       
# below contains seconds as well as minutes and everything else
A=`date "+%F_%T"`
B=`date "+%F"`
ver=8.9
backup=FCX_$(ver)_$(M)
arch=focex_$(ver)_$(M)
exe=fcx_$(ver).x

$(exe): $(SUBO) $(MAINO) 
	$(FF) $(FLAGS)  -o $@  $^  
	echo 'executable was made on' $A
	mv $(exe)  $(DESTDIR)

#ph.x: $(SUBO) $(PHO)
#	$(FF) $(FLAGS)  -o $@  $^  
#	echo 'executable was made on' $A
#	mv ph.x  $(DESTDIR)

# rules to make the "simple" object files from sources
%.o: %.f90
	$(FF) $(FLAGS) -c $<
%.o: %.f
	$(FF) $(FLAGS) -c $<


# tests , cleanup and archiving
write:
	echo $(A) $(B) $(M) $(arch) $(backup) 
clean:
	rm  *.o *.mod 
backup:
	mkdir $(backup)
	cp $(ALL) $(backup)
#	tar -cvf focex.tar $(backup)
#	gzip focex.tar

archive:
	tar -cvf $(arch) $(ALL)
	gzip $(arch)

# $@ is the name of the file to be made
# $? is the names of the changed dependents 
# $< the name of the related file that caused the action
# $* the prefix shared by target and dependent files
