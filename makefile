all: aout

aout: mod_spectra_heom.o spectra_heom.o
	ifort -o aout mod_spectra_heom.o spectra_heom.o -O2 -llapack -lblas

%.o: %.f90
	ifort -c -O2 $<

clean:
	rm *.o aout

