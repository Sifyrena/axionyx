// fftw_omp_stubs.cpp
//
// No-op stubs for FFTW's thread-init entry points.
//
// Homebrew's libfftw3_omp.dylib links /opt/homebrew/opt/libomp/lib/libomp.dylib,
// while AMReX's Clang (llvm.mak) links /opt/homebrew/opt/llvm/lib/libomp.dylib.
// Two different dylib paths → two OMP runtime initializations → "OMP: Error #15".
//
// Solution: skip linking libfftw3_omp entirely and satisfy the three symbols that
// amrex::FFT::detail::Initialize() calls with harmless no-ops.  AMReX manages
// threading above the FFTW layer via its own OpenMP/MPI model; FFTW's internal
// thread pool is never actually needed.

extern "C" {
    int  fftw_init_threads (void)    { return 1; }
    void fftw_plan_with_nthreads(int) {}
    int  fftwf_init_threads (void)   { return 1; }
    void fftwf_plan_with_nthreads(int) {}
}
