CC=gfortran
FOUTCAR=readoutcar.f90
FPWSCF=readpfpwscf.f90

all: readoutcar readpfpwscf

readoutcar: $(FOUTCAR)
	$(CC) $(FOUTCAR) -o readoutcar.x

readpfpwscf: $(FPWSCF)
	$(CC) $(FPWSCF) -o readpfpwscf.x
