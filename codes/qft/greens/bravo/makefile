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

# Blocks and allocations depends on shared# mod_allocator.o               : mod_set_precision.o
mod_constants.o               : mod_set_precision.o
mod_extents.o                 : mod_set_precision.o
mod_fields.o                  : mod_set_precision.o mod_constants.o mod_extents.o mod_masses.o mod_file_handling.o mod_inputs.o
mod_file_handling.o           :
mod_inputs.o                  : mod_set_precision.o
mod_masses.o                  : mod_set_precision.o mod_constants.o
mod_randoms.o                 : mod_set_precision.o mod_constants.o
mod_set_precision.o           :
mod_time_step.o               :

# rditldmt@ITLDMT-MD-O2034:alpha $ date
# Mon Jan 23 18:19:39 CST 2017
# rditldmt@ITLDMT-MD-O2034:alpha $ pwd
# /Users/rditldmt/Documents/GitHub_Desktop/UNM-QFT/codes/qft/greens/alpha

# rditldmt@ITLDMT-MD-O2034:alpha $ make debug
# PROGRAM  = AaM0
# PRG_OBJ  = AaM0.o
# SRCS     = AaM0.f08 mod_constants.f08 mod_extents.f08 mod_fields.f08 mod_fields_allocate.f08 mod_fields_output.f08 mod_fields_utilities.f08 mod_file_handling.f08 mod_inputs.f08 mod_masses.f08 mod_parameter_sets.f08 mod_randoms.f08 mod_set_precision.f08 mod_time_stamp.f08
# OBJS     = AaM0.o mod_constants.o mod_extents.o mod_fields.o mod_fields_allocate.o mod_fields_output.o mod_fields_utilities.o mod_file_handling.o mod_inputs.o mod_masses.o mod_parameter_sets.o mod_randoms.o mod_set_precision.o mod_time_stamp.o
# MODS     = mod_constants.f08 mod_extents.f08 mod_fields.f08 mod_fields_allocate.f08 mod_fields_output.f08 mod_fields_utilities.f08 mod_file_handling.f08 mod_inputs.f08 mod_masses.f08 mod_parameter_sets.f08 mod_randoms.f08 mod_set_precision.f08 mod_time_stamp.f08
# MOD_OBJS = mod_constants.o mod_extents.o mod_fields.o mod_fields_allocate.o mod_fields_output.o mod_fields_utilities.o mod_file_handling.o mod_inputs.o mod_masses.o mod_parameter_sets.o mod_randoms.o mod_set_precision.o mod_time_stamp.o
# rditldmt@ITLDMT-MD-O2034:alpha $

# rditldmt@ITLDMT-MD-O2034:alpha $ echo $gflags
# -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only
# rditldmt@ITLDMT-MD-O2034:alpha $

# rditldmt@ITLDMT-MD-O2034:alpha $ make
# gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_set_precision.o mod_set_precision.f08
# gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_constants.o mod_constants.f08
# gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_extents.o mod_extents.f08
# gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_masses.o mod_masses.f08
# gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_file_handling.o mod_file_handling.f08
# gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_inputs.o mod_inputs.f08
# gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_fields.o mod_fields.f08
# gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_fields_allocate.o mod_fields_allocate.f08
# mod_fields_allocate.f08:2:19:
#
#  submodule ( mFields ) smFieldsAllocate
#                    1
# Warning: USE statement at (1) has no ONLY qualifier [-Wuse-without-only]
# gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_fields_output.o mod_fields_output.f08
# mod_fields_output.f08:2:19:
#
#  submodule ( mFields ) smFieldsOutput
#                    1
# Warning: USE statement at (1) has no ONLY qualifier [-Wuse-without-only]
# gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_fields_utilities.o mod_fields_utilities.f08
# mod_fields_utilities.f08:2:19:
#
#  submodule ( mFields ) smFieldsUtilities
#                    1
# Warning: USE statement at (1) has no ONLY qualifier [-Wuse-without-only]
# gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_parameter_sets.o mod_parameter_sets.f08
# gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_randoms.o mod_randoms.f08
# gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_time_stamp.o mod_time_stamp.f08
# gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o AaM0.o AaM0.f08
# gfortran -g -o AaM0 AaM0.o mod_constants.o mod_extents.o mod_fields.o mod_fields_allocate.o mod_fields_output.o mod_fields_utilities.o mod_file_handling.o mod_inputs.o mod_masses.o mod_parameter_sets.o mod_randoms.o mod_set_precision.o mod_time_stamp.o
