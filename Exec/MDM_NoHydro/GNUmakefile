# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..

AMREX_HOME ?= ../../subprojects/amrex
FFTW_DIR ?= /opt/homebrew/Cellar/fftw/3.3.10_3

USE_HDF5=FALSE

# compilation options
COMP    = gcc

USE_MPI = TRUE
USE_OMP = TRUE
USE_CUDA = FALSE
NO_MPI_CHECKING=FALSE

# Cosmo

PRECISION = DOUBLE
DEBUG     = FALSE

# physics
DIM      = 3
USE_GRAV = TRUE
NO_HYDRO = TRUE
USE_HYDRO = FALSE
USE_FDM = TRUE
USE_PARTICLES = TRUE
AMREX_PARTICLES = TRUE
AXIONYX = TRUE
USE_AXCOMPLEX = FALSE
USE_AXREAL = FALSE
CONST_SPECIES = TRUE

USE_CGS = FALSE

Bpack := ./Make.package
Blocs := .

USE_FFT := TRUE
# AMReX native FFT uses FFTW3 locally on each rank (no fftw3_mpi needed)
INCLUDE_LOCATIONS += $(FFTW_DIR)/include
LIBRARIES         += -L$(FFTW_DIR)/lib

include $(TOP)/Exec/Make.Nyx  

# Add near the top of GNUmakefile
GIT_HASH  := $(shell git rev-parse --short HEAD 2>/dev/null || echo unknown)
GIT_BRANCH := $(shell git rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown)
GIT_DATE  := $(shell git log -1 --format=%cd --date=short 2>/dev/null || echo unknown)

$(shell echo '#define AXIONYX_GIT_HASH   "$(GIT_HASH)"'  > axionyx_version.H)
$(shell echo '#define AXIONYX_GIT_BRANCH "$(GIT_BRANCH)"' >> axionyx_version.H)
$(shell echo '#define AXIONYX_GIT_DATE   "$(GIT_DATE)"'   >> axionyx_version.H)                            
