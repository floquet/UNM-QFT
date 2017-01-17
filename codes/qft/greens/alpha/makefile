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
mod_fields.o                  : mod_set_precision.o mod_constants.o mod_extents.o
mod_file_handling.o           :
mod_randoms.o                 : mod_set_precision.o
mod_set_precision.o           :
mod_time_step.o               :

# rditldmt@ITLDMT-MD-O2034:alpha $ date
# Mon Jan 16 18:20:17 CST 2017
# rditldmt@ITLDMT-MD-O2034:alpha $ pwd
# /Users/rditldmt/hpc/fortran/projects/qft/greens/alpha
# rditldmt@ITLDMT-MD-O2034:alpha $
# rditldmt@ITLDMT-MD-O2034:alpha $ make debug
# PROGRAM  = AaM0
# PRG_OBJ  = AaM0.o
# SRCS     = AaM0.f08 mod_constants.f08 mod_extents.f08 mod_fields.f08 mod_file_handling.f08 mod_randoms.f08 mod_set_precision.f08 mod_time_stamp.f08
# OBJS     = AaM0.o mod_constants.o mod_extents.o mod_fields.o mod_file_handling.o mod_randoms.o mod_set_precision.o mod_time_stamp.o
# MODS     = mod_constants.f08 mod_extents.f08 mod_fields.f08 mod_file_handling.f08 mod_randoms.f08 mod_set_precision.f08 mod_time_stamp.f08
# MOD_OBJS = mod_constants.o mod_extents.o mod_fields.o mod_file_handling.o mod_randoms.o mod_set_precision.o mod_time_stamp.o
# rditldmt@ITLDMT-MD-O2034:alpha $
# rditldmt@ITLDMT-MD-O2034:alpha $ echo $gflags
# -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only
