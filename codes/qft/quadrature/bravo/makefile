# Find all source files, create a list of corresponding object files
PROGRAM = quadrature

EXT = f08
# Find all source files, create a list of corresponding object files
SRCS=$(wildcard *.$(EXT))
OBJS=$(patsubst %.$(EXT),%.o,$(SRCS))

# Ditto for mods (They will be in both lists)
MODS=$(wildcard mod*.$(EXT))
MOD_OBJS=$(patsubst %.$(EXT),%.o,$(MODS))

# Compiler/Linker settings
# gfortran - gnu (execution: ./a.out) (flags: gflags in ~/.bashrc)
# mpif90 - gnu, MPI Mac Desktop (execution: mpirun -np 4 ./a.out)
FC = ${fortran_compiler_mpi_gnu}
FLFLAGS = -g
FCFLAGS = -c ${gflags}
PRG_OBJ = $(PROGRAM).o

# Clean the suffixes
.SUFFIXES:

# Set the suffixes we are interested in
.SUFFIXES: .$(EXT) .o

# make without parameters will make first target found.
default : $(PROGRAM)

# Compiler steps for all objects
$(OBJS) : %.o : %.$(EXT)
	$(FC) $(FCFLAGS) -o $@ $<

# Linker
$(PROGRAM) : $(OBJS)
	$(FC) $(FLFLAGS) -o $@ $^

debug:
	@echo "${bold}PROGRAM${normal}  = $(PROGRAM)"
	@echo "${bold}PRG_OBJ${normal}  = $(PRG_OBJ)"
	@echo "${bold}SRCS${normal}     = $(SRCS)"
	@echo "${bold}OBJS${normal}     = $(OBJS)"
	@echo "${bold}MODS${normal}     = $(MODS)"
	@echo "${bold}MOD_OBJS${normal} = $(MOD_OBJS)"

clean:
	rm -rf $(OBJS) $(PROGRAM) $(patsubst %.o,%.mod,$(MOD_OBJS))
	rm -rf *.mod *.smod *.o *.dSYM

.PHONY: debug default clean

# Dependencies

# Main program depends on all modules
$(PRG_OBJ) : $(MOD_OBJS)

# Blocks and allocations depends on shared
mod_constants.o               : mod_set_precision.o
mod_cuboid.o                  : mod_set_precision.o mod_parameters.o mod_point_in_space.o
mod_domain.o                  : mod_set_precision.o
mod_mpfi_sgi.o                :
mod_parameters.o              : mod_set_precision.o
mod_partition.o               : mod_set_precision.o mod_constants.o mod_parameters.o
mod_point_in_space.o          : mod_parameters.o
mod_set_precision.o           :
mod_time_step.o               :

# rditldmt@ITLDMT-MD-O2034:bravo $ date
# Tue Jan 10 17:52:40 CST 2017
# rditldmt@ITLDMT-MD-O2034:bravo $ pwd
# /Users/rditldmt/hpc/fortran/projects/qft/quadrature/bravo
# rditldmt@ITLDMT-MD-O2034:bravo $ make debug
# PROGRAM  = quadrature
# PRG_OBJ  = quadrature.o
# SRCS     = mod_constants.f08 mod_mpif_sgi.f08 mod_parameters.f08 mod_partition.f08 mod_point_in_space.f08 mod_set_precision.f08 quadrature.f08
# OBJS     = mod_constants.o mod_mpif_sgi.o mod_parameters.o mod_partition.o mod_point_in_space.o mod_set_precision.o quadrature.o
# MODS     = mod_constants.f08 mod_mpif_sgi.f08 mod_parameters.f08 mod_partition.f08 mod_point_in_space.f08 mod_set_precision.f08
# MOD_OBJS = mod_constants.o mod_mpif_sgi.o mod_parameters.o mod_partition.o mod_point_in_space.o mod_set_precision.o

# rditldmt@ITLDMT-MD-O2034:bravo $ echo $gflags
# -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only

# rditldmt@ITLDMT-MD-O2034:bravo $ mpichversion
# MPICH Version:    	3.2
# MPICH Release date:	Wed Nov 11 22:06:48 CST 2015
# MPICH Device:    	ch3:nemesis
# MPICH configure: 	--prefix=/opt/local --disable-dependency-tracking --disable-silent-rules --enable-base-cache --enable-cache --enable-cxx --enable-fast=O2 --enable-shared --enable-smpcoll --with-device=ch3:nemesis --with-pm=hydra --with-thread-package=posix --enable-versioning F90FLAGS= F90= --enable-timer-type=mach_absolute_time --libdir=/opt/local/lib/mpich-gcc5 --sysconfdir=/opt/local/etc/mpich-gcc5 --program-suffix=-mpich-gcc5 --includedir=/opt/local/include/mpich-gcc5 --enable-f77 --enable-fc
# MPICH CC: 	/opt/local/bin/gcc-mp-5 -pipe -m64   -O2
# MPICH CXX: 	/opt/local/bin/g++-mp-5 -pipe -m64  -O2
# MPICH F77: 	/opt/local/bin/gfortran-mp-5 -pipe  -O2
# MPICH FC: 	/opt/local/bin/gfortran-mp-5 -pipe -m64  -O2
