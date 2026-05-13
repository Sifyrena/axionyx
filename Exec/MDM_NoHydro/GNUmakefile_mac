# macOS / Apple Silicon build  (Clang + Homebrew FFTW + LLVM OpenMP)
# Usage:  make -f GNUmakefile_mac -j$(sysctl -n hw.logicalcpu)

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..

AMREX_HOME ?= ../../subprojects/amrex
FFTW_DIR   ?= /opt/homebrew/opt/fftw

USE_HDF5 = FALSE

# Homebrew LLVM clang — AMReX maps COMP=llvm to llvm.mak which handles -fopenmp/libomp
COMP = llvm

USE_MPI  = TRUE
USE_OMP  = FALSE
USE_CUDA = FALSE
NO_MPI_CHECKING = FALSE

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

USE_FFT = TRUE
INCLUDE_LOCATIONS += $(FFTW_DIR)/include
LIBRARIES         += -L$(FFTW_DIR)/lib

include $(TOP)/Exec/Make.Nyx

# ── macOS OpenMP / FFTW note ──────────────────────────────────────────────────
# libfftw3_omp must be patched once after brew install to point at LLVM's libomp:
#   install_name_tool -change \
#       /opt/homebrew/opt/libomp/lib/libomp.dylib \
#       /opt/homebrew/opt/llvm/lib/libomp.dylib \
#       /opt/homebrew/opt/fftw/lib/libfftw3_omp.3.dylib
#   install_name_tool -change \
#       /opt/homebrew/opt/libomp/lib/libomp.dylib \
#       /opt/homebrew/opt/llvm/lib/libomp.dylib \
#       /opt/homebrew/opt/fftw/lib/libfftw3f_omp.3.dylib
# Without this, two libomp copies load → "OMP: Error #15" crash.
# ─────────────────────────────────────────────────────────────────────────────

GIT_HASH   := $(shell git rev-parse --short HEAD 2>/dev/null || echo unknown)
GIT_BRANCH := $(shell git rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown)
GIT_DATE   := $(shell git log -1 --format=%cd --date=short 2>/dev/null || echo unknown)

$(shell echo '#define AXIONYX_GIT_HASH   "$(GIT_HASH)"'   > axionyx_version.H)
$(shell echo '#define AXIONYX_GIT_BRANCH "$(GIT_BRANCH)"' >> axionyx_version.H)
$(shell echo '#define AXIONYX_GIT_DATE   "$(GIT_DATE)"'   >> axionyx_version.H)

# Apple Silicon: Homebrew LLVM ld doesn't add an ad-hoc signature automatically.
# macOS refuses to run unsigned arm64 binaries.  Sign after every link.
all_targets: $(executable)
	codesign --sign - --force $(executable)
