# Find all source files, create a list of corresponding object files
PROGRAM = AaM0

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
#FC = ${fortran_compiler_mpi_gnu}
FC      = gfortran
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
# mod_allocator.o               : mod_set_precision.o
mod_constants.o               : mod_set_precision.o mod_constants.o
mod_extents.o                 : mod_set_precision.o
mod_fields.o                  : mod_set_precision.o mod_constants.o mod_extents.o mod_masses.o mod_file_handling.o mod_inputs.o
mod_file_handling.o           :
mod_inputs.o                  : mod_set_precision.o
mod_masses.o                  : mod_set_precision.o mod_constants.o
mod_randoms.o                 : mod_set_precision.o
mod_set_precision.o           :
mod_time_step.o               :

# dantopa@Muntz-Szasz.attlocal.net:alpha $ date
# Sun Jan 22 15:03:03 CST 2017
# dantopa@Muntz-Szasz.attlocal.net:alpha $ pwd
# /Users/dantopa/Documents/GitHub_desktop/UNM-QFT/codes/qft/greens/alpha
# dantopa@Muntz-Szasz.attlocal.net:alpha $ echo $gflags
# -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only
# dantopa@Muntz-Szasz.attlocal.net:alpha $ gcc --version
# gcc (MacPorts gcc7 7-20170115_0) 7.0.0 20170115 (experimental)
# Copyright (C) 2017 Free Software Foundation, Inc.
# This is free software; see the source for copying conditions.  There is NO
# warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# dantopa@Muntz-Szasz.attlocal.net:alpha $ make debug
# PROGRAM  = AaM0
# PRG_OBJ  = AaM0.o
# SRCS     = AaM0.f08 mod_constants.f08 mod_extents.f08 mod_fields.f08 mod_file_handling.f08 mod_inputs.f08 mod_masses.f08 mod_parameter_sets.f08 mod_randoms.f08 mod_set_precision.f08 mod_time_stamp.f08
# OBJS     = AaM0.o mod_constants.o mod_extents.o mod_fields.o mod_file_handling.o mod_inputs.o mod_masses.o mod_parameter_sets.o mod_randoms.o mod_set_precision.o mod_time_stamp.o
# MODS     = mod_constants.f08 mod_extents.f08 mod_fields.f08 mod_file_handling.f08 mod_inputs.f08 mod_masses.f08 mod_parameter_sets.f08 mod_randoms.f08 mod_set_precision.f08 mod_time_stamp.f08
# MOD_OBJS = mod_constants.o mod_extents.o mod_fields.o mod_file_handling.o mod_inputs.o mod_masses.o mod_parameter_sets.o mod_randoms.o mod_set_precision.o mod_time_stamp.o
# dantopa@Muntz-Szasz.attlocal.net:alpha $ make
# gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_set_precision.o mod_set_precision.f08
# make: Circular mod_constants.o <- mod_constants.o dependency dropped.
# gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_constants.o mod_constants.f08
# gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_extents.o mod_extents.f08
# gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_masses.o mod_masses.f08
# gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_file_handling.o mod_file_handling.f08
# gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_inputs.o mod_inputs.f08
# gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_fields.o mod_fields.f08
# mod_fields.f08:66:0:
#
#      subroutine write_f_sub ( me, myInputs )
#
# Warning: 'write_f_sub' defined but not used [-Wunused-function]
# mod_fields.f08:170:0:
#
#      subroutine extrema_sub ( me )
#
# Warning: 'extrema_sub' defined but not used [-Wunused-function]
# gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_parameter_sets.o mod_parameter_sets.f08
# gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_randoms.o mod_randoms.f08
# gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_time_stamp.o mod_time_stamp.f08
# gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o AaM0.o AaM0.f08
# gfortran -g -o AaM0 AaM0.o mod_constants.o mod_extents.o mod_fields.o mod_file_handling.o mod_inputs.o mod_masses.o mod_parameter_sets.o mod_randoms.o mod_set_precision.o mod_time_stamp.o
