RM = rm -f
CC = gcc
FC = gfortran

CCFLAGS = -Wall
FCFLAGS = -Wall

# solvers/ directory that can be zipped and downloaded from netlib
ASL_DIR = ../solvers
# directory holding the compiled MA28 library
MA28_DIR = ../ma28-1.0.0/src

LIBASL = $(ASL_DIR)/amplsolver.a
LIBMA28 = $(MA28_DIR)/libma28.a

DEPS = -ldl -lm -lgfortran

BIN_DIR = ..

newton: newton.c iter_log.c convergence.c linalg.c
	$(CC) $(CCFLAGS) -I$(ASL_DIR) $^ $(LIBASL) $(LIBMA28) $(DEPS) -o $(BIN_DIR)/$@

clean:
	$(RM) $(BIN_DIR)/newton
