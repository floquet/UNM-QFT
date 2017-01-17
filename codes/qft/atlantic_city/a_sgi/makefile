# Find all source files, create a list of corresponding object files
PROGRAM = makesAC

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
	rm -rf $(OBJS) $(patsubst %.o,%.mod,$(MOD_OBJS))
	rm -f *.mod *.smod *.o *.dSYM

.PHONY: debug default clean

# Dependencies

# Main program depends on all modules
$(PRG_OBJ) : $(MOD_OBJS)

# Blocks and allocations depends on shared
mod_file_handling.o           :
mod_library_functions.o       : mod_set_precision.o
mod_mpif_sgi.o                :
mod_set_precision.o           :

# dantopa@thunder08:a_sgi $ date
# Mon Jan  9 13:52:07 EST 2017
# dantopa@thunder08:a_sgi $ pwd
# /home/dantopa/hpc/fortran/projects/qft/atlantic_city/a_sgi
# dantopa@thunder08:a_sgi $ make debug
# PROGRAM  = makesAC
# PRG_OBJ  = makesAC.o
# SRCS     = makesAC.f08 mod_file_handling.f08 mod_library_functions.f08 mod_mpif_sgi.f08 mod_set_precision.f08
# OBJS     = makesAC.o mod_file_handling.o mod_library_functions.o mod_mpif_sgi.o mod_set_precision.o
# MODS     = mod_file_handling.f08 mod_library_functions.f08 mod_mpif_sgi.f08 mod_set_precision.f08
# MOD_OBJS = mod_file_handling.o mod_library_functions.o mod_mpif_sgi.o mod_set_precision.o
# dantopa@thunder08:a_sgi $ make
# mpif90 -c -g -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -O -pedantic -fcheck=bounds -fmax-errors=5 -o mod_file_handling.o mod_file_handling.f08
# mpif90 -c -g -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -O -pedantic -fcheck=bounds -fmax-errors=5 -o mod_set_precision.o mod_set_precision.f08
# mpif90 -c -g -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -O -pedantic -fcheck=bounds -fmax-errors=5 -o mod_library_functions.o mod_library_functions.f08
# mpif90 -c -g -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -O -pedantic -fcheck=bounds -fmax-errors=5 -o mod_mpif_sgi.o mod_mpif_sgi.f08
# mpif90 -c -g -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -O -pedantic -fcheck=bounds -fmax-errors=5 -o makesAC.o makesAC.f08
# makesAC.f08:46:15:
# 
#          proc = myproc
#                1
# Warning: Conversion from INTEGER(4) to INTEGER(8) at (1) [-Wconversion-extra]
# makesAC.f08:47:18:
# 
#          numproc = mynumproc
#                   1
# Warning: Conversion from INTEGER(4) to INTEGER(8) at (1) [-Wconversion-extra]
# mpif90 -g -o makesAC makesAC.o mod_file_handling.o mod_library_functions.o mod_mpif_sgi.o mod_set_precision.o
# dantopa@thunder08:a_sgi $ make clean
# rm -rf makesAC.o mod_file_handling.o mod_library_functions.o mod_mpif_sgi.o mod_set_precision.o mod_file_handling.mod mod_library_functions.mod mod_mpif_sgi.mod mod_set_precision.mod
# rm -f *.mod *.smod *.o *.dSYM

