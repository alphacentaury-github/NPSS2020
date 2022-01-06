# October 2015
# AE Lovell for LJ Titus' local transfer code

FFLAGS = -std=legacy -fmax-errors=0 	# specific for gfortran 
FC=gfortran
F90C=$(FC)

COMMON = front_end.o main.o polpak.o bound_state.o scattering_state.o tmatrix.o local_binding_potental.o \
	 nmout.o usewitt.o nmin.o calc_binding_source.o coulfg.o local_scattering_potential.o nm.o \
	 nonlocal_adiabatic.o calc_source.o local_adiabatic.o cgamma.o gaussquad.o cg.o \
	 nonlocal_binding_potential.o mcsphjy.o nonlocal_scattering_potential.o CoulombPhase.o diffCS.o racah2.o

# name the executable and put it wherever you want	 
target: NLAT
install: NLAT
#	cp -fp NLAT $(HOME)/bin 	#or wherever you want the executable
	mv NLAT ../

NLAT: $(COMMON)
	$(FC) $(FFLAGS) -o NLAT $(COMMON)

# allows you to compile all of the types of fortran files
%.o: %.f90
	$(FC) $(FORM) -c $(FFLAGS) $<

%.o: %.f
	$(FC) $(FORM) -c $(FFLAGS) $<

%.o: %.for
	$(FC) $(FORM) -c $(FFLAGS) $<

# remove the .o files	
clean:
	\rm -f *.o

# to compile > make install clean
