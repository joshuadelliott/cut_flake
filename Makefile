FC = gfortran
FFLAGS = -fcheck=bounds 

all: flake protonate 

trash: 
	rm -r *.mod *.o *.x

cell_data.o: cell_data.f90
	$(FC) $(FFLAGS) -c cell_data.f90

xyz_types.o: xyz_types.f90
	$(FC) $(FFLAGS) -c xyz_types.f90

xyz_io.o: xyz_io.f90
	$(FC) $(FFLAGS) -c xyz_io.f90

xyz_ops.o: xyz_ops.f90
	$(FC) $(FFLAGS) -c xyz_ops.f90


flake: xyz_types.o cell_data.o xyz_io.o xyz_ops.o
	$(FC) $(FFLAGS) -o flake.x flake.f90 \
	cell_data.o \
	xyz_types.o \
	xyz_io.o \
	xyz_ops.o 

protonate: xyz_types.o cell_data.o xyz_io.o xyz_ops.o
	$(FC) $(FFLAGS) -o protonate.x protonate.f90 \
	cell_data.o \
	xyz_types.o \
	xyz_io.o \
	xyz_ops.o 
