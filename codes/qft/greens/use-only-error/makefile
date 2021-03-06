# Find all source files, create a list of corresponding object files
PROGRAM = demonstration

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
#FLFLAGS = -g
#FCFLAGS = -c ${gflags}
FCFLAGS = -c -Wuse-without-only
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
# mod_constants.o               : mod_set_precision.o
# mod_extents.o                 : mod_set_precision.o
# mod_fields.o                  : mod_set_precision.o mod_constants.o mod_extents.o mod_masses.o mod_file_handling.o mod_inputs.o
# mod_file_handling.o           :
# mod_inputs.o                  : mod_set_precision.o
# mod_masses.o                  : mod_set_precision.o mod_constants.o
# mod_randoms.o                 : mod_set_precision.o mod_constants.o
# mod_set_precision.o           :
# mod_time_step.o               :
