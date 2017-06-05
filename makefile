FC=gfortran
FFLAG=-lfftw3 -I /home/junxie/opt/fftw/include -L /home/junxie/opt/fftw/lib -fbounds-check
FFLAG=-lfftw3 -I /usr/include -L /usr/lib64 -fbounds-check
objects=pre_processing.o filter.o do_norm1.o do_norm2.o do_norm3.o smooth.o smoothf.o taperf.o do_whiten1.o do_whiten2.o do_whiten3.o sacio.o globe_data.o
all:sacio.mod globe_data.mod ../bin/pre_processing
.f.o:
	$(FC) $(FFLAG) $< -c
%.o:%.f90
	$(FC) $(FFLAG) $< -c 
sacio.mod:sacio.f90
	$(FC) $< -c
globe_data.mod:globe_data.f90
	$(FC) $< -c
../bin/pre_processing:$(objects)
	$(FC) $(FFLAG) $^ -o $@ 
clean:
	-rm *.o *.mod 
