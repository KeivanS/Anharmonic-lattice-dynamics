include ../make.inc

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


kap$(version): $(MODS) $(SUBKAPS) $(KAPA) 
	$(F90) $(FLAGS) $(MODS) $(SUBKAPS) $(KAPA) -o $@ 
	mv kap$(version) $(DESTDIR)

mfpband: mfpthermal.f90
	$(F90) mfpthermal.f90 -o $@
	mv mfpband $(DESTDIR)

# rules to make the "simple" object files from sources
.SUFFIXES: .f90
%.o : %.mod
.f90.o:
	$(F90) $(FLAGS) -c $<
.f90 :
	$(F90) $(FLAGS) $< -o $@
.f.o:
	$(F90) $(FLAGS) -c $<
.f :
	$(F90) $(FLAGS) $< -o $@

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
