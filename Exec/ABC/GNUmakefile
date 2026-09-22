# Cluster build  (GCC + system FFTW + MPI)
# Usage:  make -j<N>

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..

AMREX_HOME ?= ../../subprojects/amrex
FFTW_DIR ?= /sw/rev/25.04/cascadelake_opx_rocky8/linux-rocky8-cascadelake/gcc-11.5.0/fftw-3.3.10-5572zu2ox4k6xsy3yfkrblzykhmenkgi/
USE_HDF5 = FALSE

# compilation options
COMP = gcc

USE_MPI  = TRUE
USE_OMP  = TRUE
USE_CUDA = FALSE
NO_MPI_CHECKING = FALSE

PRECISION = DOUBLE
DEBUG     = FALSE

# physics
DIM      = 3
USE_GRAV = TRUE
NO_HYDRO = FALSE
USE_HYDRO = TRUE
USE_HEATCOOL = TRUE
USE_FDM = TRUE
USE_PARTICLES = TRUE
AMREX_PARTICLES = TRUE
AXIONYX = TRUE
USE_AXCOMPLEX = FALSE
USE_AXREAL = FALSE
USE_CONST_SPECIES = TRUE

USE_CGS = FALSE

# Heating/cooling network (USE_HEATCOOL) needs the real CVode integrator,
# which needs SUNDIALS. Local build, not a cluster module -- lib64, not
# lib, is where this particular CMake install put the .so/.a files.
USE_SUNDIALS = TRUE
SUNDIALS_ROOT ?= $(HOME)/sundials/install
SUNDIALS_LIB_DIR = $(SUNDIALS_ROOT)/lib64
LDFLAGS += -Wl,-rpath,$(SUNDIALS_LIB_DIR)

Bpack := ./Make.package
Blocs := .

USE_FFT = TRUE

# AMReX native FFT uses FFTW3 locally on each rank (no fftw3_mpi needed)
INCLUDE_LOCATIONS += $(FFTW_DIR)/include
LIBRARIES         += -L$(FFTW_DIR)/lib

include $(TOP)/Exec/Make.Nyx

# This ld (binutils 2.43) does not copy transitive DT_NEEDED entries, so
# libsundials_cvode.so needing libsundials_core.so is not enough on its
# own -- core has to be listed explicitly, after cvode/arkode above.
LIBRARIES += -lsundials_core

# Same story for libmpi.so needing opal_common_ofi_* from libopen-pal.so --
# it's already on mpicxx's own -L search path, just never requested with
# -l. No new path to hardcode; the directory is whatever mpicxx supplies.
LIBRARIES += -lopen-pal

GIT_HASH   := $(shell git rev-parse --short HEAD 2>/dev/null || echo unknown)
GIT_BRANCH := $(shell git rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown)
GIT_DATE   := $(shell git log -1 --format=%cd --date=short 2>/dev/null || echo unknown)

$(shell echo '#define AXIONYX_GIT_HASH   "$(GIT_HASH)"'   > axionyx_version.H)
$(shell echo '#define AXIONYX_GIT_BRANCH "$(GIT_BRANCH)"' >> axionyx_version.H)
$(shell echo '#define AXIONYX_GIT_DATE   "$(GIT_DATE)"'   >> axionyx_version.H)
