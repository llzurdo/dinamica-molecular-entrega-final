MAKEFILE = Makefile
exe = simple
fcomp = gfortran
flags = -O3

OBJS = ziggurat.o globals.o md_module.o md_integrator.o simple.o

.SUFFIXES:
.SUFFIXES: .f90 .o

.f90.o:
	$(fcomp) -c $(flags) $<

$(exe): $(OBJS) Makefile
	$(fcomp) $(flags) -o $(exe) $(OBJS)

clean:
	rm -f *.o *.mod

# Dependencias
simple.o: simple.f90 ziggurat.o globals.o md_module.o md_integrator.o
md_integrator.o: md_integrator.f90 globals.o ziggurat.o
md_module.o: md_module.f90 globals.o
ziggurat.o: ziggurat.f90
globals.o: globals.f90

