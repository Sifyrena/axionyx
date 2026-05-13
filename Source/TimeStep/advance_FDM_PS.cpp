// This version makes use of AMReX native FFT (C2C via FFT::R2C<Real,both,true>)

#ifdef FDM
#ifdef BL_USE_MPI
#include <AMReX_BLProfiler.H>
#include <Nyx.H>
#include <AMReX_MultiFab.H>
#ifdef GRAVITY
#include <Gravity.H>
#endif
#include <AMReX_FFT.H>
#include <AMReX_GpuComplex.H>

#include <string>

using namespace amrex;

using cMF      = FabArray<BaseFab<GpuComplex<Real>>>;
using AmrexFFT = FFT::R2C<Real, FFT::Direction::both, true>;

// Cached FFT object — FFTW plan creation is expensive; construct once per domain.
static std::unique_ptr<AmrexFFT> s_fft;
static Box                       s_fft_domain;

static AmrexFFT& get_fft(Box const& domain)
{
    if (!s_fft || s_fft_domain != domain) {
        s_fft        = std::make_unique<AmrexFFT>(domain);
        s_fft_domain = domain;
        // Register a one-time cleanup so s_fft is destroyed *inside*
        // amrex::Finalize(), before AMReX tears down MPI and vtables.
        // Without this, the static unique_ptr destructs after Finalize()
        // and triggers "pure virtual function called".
        static bool registered = false;
        if (!registered) {
            amrex::ExecOnFinalize([]{ s_fft.reset(); });
            registered = true;
        }
    }
    return *s_fft;
}

void drift(AmrexFFT& fft, cMF& psi, MultiFab& Ax_new,
           Box const& domain, Real dt, Real h, Real a_half, Real hbaroverm);

void fdm_timestep(AmrexFFT& fft, cMF& psi,
                  MultiFab& Ax_new, MultiFab& phi, Gravity* gravity, Geometry& geom,
                  int level, Real h, Real dt_c, Real a_c,
                  Real dt_d, Real a_d, Real a_new, Real hbaroverm);

void fdm_timestep_NG(AmrexFFT& fft, cMF& psi,
                     MultiFab& Ax_new, MultiFab& phi, Geometry& geom,
                     int level, Real h, Real dt_c, Real a_c,
                     Real dt_d, Real a_d, Real a_new, Real hbaroverm);

void Nyx::advance_FDM_PS(amrex::Real time,
                         amrex::Real dt,
                         amrex::Real a_old,
                         amrex::Real a_new)
{
  BL_PROFILE("Nyx::advance_FDM_PS_NEW()");

  // *****************************************
  // define constants
  // *****************************************
  const Real h = geom.CellSize(0);

  Real w0, w1, w2, w3;
  Real c1, c2, c3, c4, c5, c6, c7, c8;
  Real d1, d2, d3, d4, d5, d6, d7, d8;
  Real a_c1, a_c2, a_c3, a_c4, a_c5, a_c6, a_c7, a_c8;
  Real a_d1, a_d2, a_d3, a_d4, a_d5, a_d6, a_d7, a_d8;
  Real a_time;

  //******************************************
  // define weights for time steps (see Levkov et. al. 2018)
  //******************************************

  if (order == 6)
  {
    w1 = -1.17767998417887;
    w2 = 0.235573213359359;
    w3 = 0.784513610477560;
    w0 = 1.0 - 2.0 * (w1 + w2 + w3);

    c1 = w3 / 2.0 * dt;
    c2 = (w2 + w3) / 2.0 * dt;
    c3 = (w1 + w2) / 2.0 * dt;
    c4 = (w0 + w1) / 2.0 * dt;
    c5 = (w0 + w1) / 2.0 * dt;
    c6 = (w1 + w2) / 2.0 * dt;
    c7 = (w2 + w3) / 2.0 * dt;
    c8 = w3 / 2.0 * dt;

    d1 = w3 * dt;
    d2 = w2 * dt;
    d3 = w1 * dt;
    d4 = w0 * dt;
    d5 = w1 * dt;
    d6 = w2 * dt;
    d7 = w3 * dt;
    d8 = 0.0;

    a_time = state[Axion_Type].prevTime() + 0.5 * c1;
    a_c1 = get_comoving_a(a_time);
    a_time += 0.5 * (c1 + c2);
    a_c2 = get_comoving_a(a_time);
    a_time += 0.5 * (c2 + c3);
    a_c3 = get_comoving_a(a_time);
    a_time += 0.5 * (c3 + c4);
    a_c4 = get_comoving_a(a_time);
    a_time += 0.5 * (c4 + c5);
    a_c5 = get_comoving_a(a_time);
    a_time += 0.5 * (c5 + c6);
    a_c6 = get_comoving_a(a_time);
    a_time += 0.5 * (c6 + c7);
    a_c7 = get_comoving_a(a_time);
    a_time += 0.5 * (c7 + c8);
    a_c8 = get_comoving_a(a_time);

    a_time = state[Axion_Type].prevTime() + 0.5 * d1;
    a_d1 = get_comoving_a(a_time);
    a_time += 0.5 * (d1 + d2);
    a_d2 = get_comoving_a(a_time);
    a_time += 0.5 * (d2 + d3);
    a_d3 = get_comoving_a(a_time);
    a_time += 0.5 * (d3 + d4);
    a_d4 = get_comoving_a(a_time);
    a_time += 0.5 * (d4 + d5);
    a_d5 = get_comoving_a(a_time);
    a_time += 0.5 * (d5 + d6);
    a_d6 = get_comoving_a(a_time);
    a_time += 0.5 * (d6 + d7);
    a_d7 = get_comoving_a(a_time);
    a_time += 0.5 * (d7 + d8);
    a_d8 = get_comoving_a(a_time);
  }
  else if (order == 2)
  {
    c1 = 0.5 * dt;
    c2 = 0.5 * dt;

    d1 = dt;
    d2 = 0.0;

    a_time = state[Axion_Type].prevTime() + 0.5 * c1;
    a_c1 = get_comoving_a(a_time);
    a_time += 0.5 * (c1 + c2);
    a_c2 = get_comoving_a(a_time);

    a_time = state[Axion_Type].prevTime() + 0.5 * d1;
    a_d1 = get_comoving_a(a_time);
    a_time += 0.5 * (d1 + d2);
    a_d2 = get_comoving_a(a_time);
  }
  else
    amrex::Error("Order of algorithm is not implemented! Use 2 or 6 Please!");

  // *****************************************
  // Get Axion State
  // *****************************************
  MultiFab& Ax_old = get_old_data(Axion_Type);
  MultiFab& Ax_new = get_new_data(Axion_Type);

  if (Ax_old.contains_nan(0, Ax_old.nComp(), 0))
  {
    for (int i = 0; i < Ax_old.nComp(); i++)
    {
      if (Ax_old.contains_nan(i, 1, 0))
      {
        std::cout << "Testing component i for NaNs: " << i << std::endl;
        amrex::Abort("Ax_old has NaNs in this component::advance_FDM_FFT()");
      }
    }
  }

  // *****************************************
  // Define Potential
  // *****************************************
  MultiFab phi(Ax_new.boxArray(), Ax_new.DistributionMap(), 1, 1);
  phi.setVal(0.0);
#ifdef GRAVITY
  MultiFab& grav_phi = get_old_data(PhiGrav_Type);
  phi.ParallelCopy(grav_phi, 0, 0, 1, 1, 1, parent->Geom(level).periodicity());
#endif

  // *****************************************
  // Construct FFT and load wavefunction into complex MultiFab
  // *****************************************
  AmrexFFT& fft = get_fft(geom.Domain());
  cMF psi(Ax_old.boxArray(), Ax_old.DistributionMap(), 1, 0);

  for (MFIter mfi(psi); mfi.isValid(); ++mfi)
  {
    auto const ax = Ax_old.const_array(mfi);
    auto p = psi[mfi].array();
    amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
      p(i,j,k) = GpuComplex<Real>(ax(i,j,k,Nyx::AxRe), ax(i,j,k,Nyx::AxIm));
    });
  }

  //  *******************************************
  //  higher order time steps
  //  *******************************************
  if (order == 6)
  {
    fdm_timestep(fft, psi, Ax_new, phi, gravity, geom, level, h, c1, a_c1, d1, a_d1, a_new, hbaroverm);
    fdm_timestep(fft, psi, Ax_new, phi, gravity, geom, level, h, c2, a_c2, d2, a_d2, a_new, hbaroverm);
    fdm_timestep(fft, psi, Ax_new, phi, gravity, geom, level, h, c3, a_c3, d3, a_d3, a_new, hbaroverm);
    fdm_timestep(fft, psi, Ax_new, phi, gravity, geom, level, h, c4, a_c4, d4, a_d4, a_new, hbaroverm);
    fdm_timestep(fft, psi, Ax_new, phi, gravity, geom, level, h, c5, a_c5, d5, a_d5, a_new, hbaroverm);
    fdm_timestep(fft, psi, Ax_new, phi, gravity, geom, level, h, c6, a_c6, d6, a_d6, a_new, hbaroverm);
    fdm_timestep(fft, psi, Ax_new, phi, gravity, geom, level, h, c7, a_c7, d7, a_d7, a_new, hbaroverm);
    fdm_timestep(fft, psi, Ax_new, phi, gravity, geom, level, h, c8, a_c8, d8, a_d8, a_new, hbaroverm);
  }
  else if (order == 2)
  {
    fdm_timestep(fft, psi, Ax_new, phi, gravity, geom, level, h, c1, a_c1, d1, a_d1, a_new, hbaroverm);
    fdm_timestep(fft, psi, Ax_new, phi, gravity, geom, level, h, c2, a_c2, d2, a_d2, a_new, hbaroverm);
  }

  //  *******************************************
  //  copy everything back to Axion_State
  //  *******************************************
  for (MFIter mfi(Ax_new); mfi.isValid(); ++mfi)
  {
    auto const p    = psi[mfi].const_array();
    auto arr        = Ax_new[mfi].array();
    auto const axold = Ax_old[mfi].const_array();
    amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
      Real re = p(i,j,k).real();
      Real im = p(i,j,k).imag();
      arr(i,j,k,Nyx::AxRe)   = re;
      arr(i,j,k,Nyx::AxIm)   = im;
      arr(i,j,k,Nyx::AxDens) = re*re + im*im;
      arr(i,j,k,Nyx::AxPhas) = std::atan2(im, re);
      int N = int((axold(i,j,k,Nyx::AxPhas) - arr(i,j,k,Nyx::AxPhas)) / (2.0*M_PI));
      if (std::abs(axold(i,j,k,Nyx::AxPhas) - arr(i,j,k,Nyx::AxPhas) - 2.0*M_PI*N) < M_PI)
      {
        arr(i,j,k,Nyx::AxPhas) += 2*M_PI*N;
      }
      else
      {
        arr(i,j,k,Nyx::AxPhas) += 2*M_PI*(N + int((0 < axold(i,j,k,Nyx::AxPhas)) - (axold(i,j,k,Nyx::AxPhas) < 0)));
      }
    });
  }
  Ax_new.FillBoundary(geom.periodicity());

  if (Ax_new.contains_nan(0, Ax_new.nComp(), 0))
  {
    for (int i = 0; i < Ax_new.nComp(); i++)
    {
      if (Ax_new.contains_nan(i, 1, 0))
      {
        std::cout << "Testing component i for NaNs: " << i << std::endl;
        amrex::Abort("Ax_new has NaNs in this component::advance_FDM_FFT()");
      }
    }
  }
}

void Nyx::advance_FDM_PS_NG(amrex::Real time,
                             amrex::Real dt,
                             amrex::Real a_old,
                             amrex::Real a_new)
{
  BL_PROFILE("Nyx::advance_FDM_PS_NEW()");

  // *****************************************
  // define constants
  // *****************************************
  const Real h = geom.CellSize(0);

  Real w0, w1, w2, w3;
  Real c1, c2, c3, c4, c5, c6, c7, c8;
  Real d1, d2, d3, d4, d5, d6, d7, d8;
  Real a_c1, a_c2, a_c3, a_c4, a_c5, a_c6, a_c7, a_c8;
  Real a_d1, a_d2, a_d3, a_d4, a_d5, a_d6, a_d7, a_d8;
  Real a_time;

  //******************************************
  // define weights for time steps (see Levkov et. al. 2018)
  //******************************************

  if (order == 6)
  {
    w1 = -1.17767998417887;
    w2 = 0.235573213359359;
    w3 = 0.784513610477560;
    w0 = 1.0 - 2.0 * (w1 + w2 + w3);

    c1 = w3 / 2.0 * dt;
    c2 = (w2 + w3) / 2.0 * dt;
    c3 = (w1 + w2) / 2.0 * dt;
    c4 = (w0 + w1) / 2.0 * dt;
    c5 = (w0 + w1) / 2.0 * dt;
    c6 = (w1 + w2) / 2.0 * dt;
    c7 = (w2 + w3) / 2.0 * dt;
    c8 = w3 / 2.0 * dt;

    d1 = w3 * dt;
    d2 = w2 * dt;
    d3 = w1 * dt;
    d4 = w0 * dt;
    d5 = w1 * dt;
    d6 = w2 * dt;
    d7 = w3 * dt;
    d8 = 0.0;

    a_time = state[Axion_Type].prevTime() + 0.5 * c1;
    a_c1 = get_comoving_a(a_time);
    a_time += 0.5 * (c1 + c2);
    a_c2 = get_comoving_a(a_time);
    a_time += 0.5 * (c2 + c3);
    a_c3 = get_comoving_a(a_time);
    a_time += 0.5 * (c3 + c4);
    a_c4 = get_comoving_a(a_time);
    a_time += 0.5 * (c4 + c5);
    a_c5 = get_comoving_a(a_time);
    a_time += 0.5 * (c5 + c6);
    a_c6 = get_comoving_a(a_time);
    a_time += 0.5 * (c6 + c7);
    a_c7 = get_comoving_a(a_time);
    a_time += 0.5 * (c7 + c8);
    a_c8 = get_comoving_a(a_time);

    a_time = state[Axion_Type].prevTime() + 0.5 * d1;
    a_d1 = get_comoving_a(a_time);
    a_time += 0.5 * (d1 + d2);
    a_d2 = get_comoving_a(a_time);
    a_time += 0.5 * (d2 + d3);
    a_d3 = get_comoving_a(a_time);
    a_time += 0.5 * (d3 + d4);
    a_d4 = get_comoving_a(a_time);
    a_time += 0.5 * (d4 + d5);
    a_d5 = get_comoving_a(a_time);
    a_time += 0.5 * (d5 + d6);
    a_d6 = get_comoving_a(a_time);
    a_time += 0.5 * (d6 + d7);
    a_d7 = get_comoving_a(a_time);
    a_time += 0.5 * (d7 + d8);
    a_d8 = get_comoving_a(a_time);
  }
  else if (order == 2)
  {
    c1 = 0.5 * dt;
    c2 = 0.5 * dt;

    d1 = dt;
    d2 = 0.0;

    a_time = state[Axion_Type].prevTime() + 0.5 * c1;
    a_c1 = get_comoving_a(a_time);
    a_time += 0.5 * (c1 + c2);
    a_c2 = get_comoving_a(a_time);

    a_time = state[Axion_Type].prevTime() + 0.5 * d1;
    a_d1 = get_comoving_a(a_time);
    a_time += 0.5 * (d1 + d2);
    a_d2 = get_comoving_a(a_time);
  }
  else
    amrex::Error("Order of algorithm is not implemented! Use 2 or 6 Please!");

  // *****************************************
  // Get Axion State
  // *****************************************
  MultiFab& Ax_old = get_old_data(Axion_Type);
  MultiFab& Ax_new = get_new_data(Axion_Type);

  if (Ax_old.contains_nan(0, Ax_old.nComp(), 0))
  {
    for (int i = 0; i < Ax_old.nComp(); i++)
    {
      if (Ax_old.contains_nan(i, 1, 0))
      {
        std::cout << "Testing component i for NaNs: " << i << std::endl;
        amrex::Abort("Ax_old has NaNs in this component::advance_FDM_FFT()");
      }
    }
  }

  // *****************************************
  // Define Potential
  // *****************************************
  MultiFab phi(Ax_new.boxArray(), Ax_new.DistributionMap(), 1, 1);
  phi.setVal(0.0);
#ifdef GRAVITY
  MultiFab& grav_phi = get_old_data(PhiGrav_Type);
  phi.ParallelCopy(grav_phi, 0, 0, 1, 1, 1, parent->Geom(level).periodicity());
#endif

  // *****************************************
  // Construct FFT and load wavefunction into complex MultiFab
  // *****************************************
  AmrexFFT& fft = get_fft(geom.Domain());
  cMF psi(Ax_old.boxArray(), Ax_old.DistributionMap(), 1, 0);

  for (MFIter mfi(psi); mfi.isValid(); ++mfi)
  {
    auto const ax = Ax_old.const_array(mfi);
    auto p = psi[mfi].array();
    amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
      p(i,j,k) = GpuComplex<Real>(ax(i,j,k,Nyx::AxRe), ax(i,j,k,Nyx::AxIm));
    });
  }

  //  *******************************************
  //  higher order time steps
  //  *******************************************
  if (order == 6)
  {
    fdm_timestep_NG(fft, psi, Ax_new, phi, geom, level, h, c1, a_c1, d1, a_d1, a_new, hbaroverm);
    fdm_timestep_NG(fft, psi, Ax_new, phi, geom, level, h, c2, a_c2, d2, a_d2, a_new, hbaroverm);
    fdm_timestep_NG(fft, psi, Ax_new, phi, geom, level, h, c3, a_c3, d3, a_d3, a_new, hbaroverm);
    fdm_timestep_NG(fft, psi, Ax_new, phi, geom, level, h, c4, a_c4, d4, a_d4, a_new, hbaroverm);
    fdm_timestep_NG(fft, psi, Ax_new, phi, geom, level, h, c5, a_c5, d5, a_d5, a_new, hbaroverm);
    fdm_timestep_NG(fft, psi, Ax_new, phi, geom, level, h, c6, a_c6, d6, a_d6, a_new, hbaroverm);
    fdm_timestep_NG(fft, psi, Ax_new, phi, geom, level, h, c7, a_c7, d7, a_d7, a_new, hbaroverm);
    fdm_timestep_NG(fft, psi, Ax_new, phi, geom, level, h, c8, a_c8, d8, a_d8, a_new, hbaroverm);
  }
  else if (order == 2)
  {
    fdm_timestep_NG(fft, psi, Ax_new, phi, geom, level, h, c1, a_c1, d1, a_d1, a_new, hbaroverm);
    fdm_timestep_NG(fft, psi, Ax_new, phi, geom, level, h, c2, a_c2, d2, a_d2, a_new, hbaroverm);
  }

  //  *******************************************
  //  copy everything back to Axion_State
  //  *******************************************
  for (MFIter mfi(Ax_new); mfi.isValid(); ++mfi)
  {
    auto const p     = psi[mfi].const_array();
    auto arr         = Ax_new[mfi].array();
    auto const axold = Ax_old[mfi].const_array();
    amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
      Real re = p(i,j,k).real();
      Real im = p(i,j,k).imag();
      arr(i,j,k,Nyx::AxRe)   = re;
      arr(i,j,k,Nyx::AxIm)   = im;
      arr(i,j,k,Nyx::AxDens) = re*re + im*im;
      arr(i,j,k,Nyx::AxPhas) = std::atan2(im, re);
      int N = int((axold(i,j,k,Nyx::AxPhas) - arr(i,j,k,Nyx::AxPhas)) / (2.0*M_PI));
      if (std::abs(axold(i,j,k,Nyx::AxPhas) - arr(i,j,k,Nyx::AxPhas) - 2.0*M_PI*N) < M_PI)
      {
        arr(i,j,k,Nyx::AxPhas) += 2*M_PI*N;
      }
      else
      {
        arr(i,j,k,Nyx::AxPhas) += 2*M_PI*(N + int((0 < axold(i,j,k,Nyx::AxPhas)) - (axold(i,j,k,Nyx::AxPhas) < 0)));
      }
    });
  }
  Ax_new.FillBoundary(geom.periodicity());

  if (Ax_new.contains_nan(0, Ax_new.nComp(), 0))
  {
    for (int i = 0; i < Ax_new.nComp(); i++)
    {
      if (Ax_new.contains_nan(i, 1, 0))
      {
        std::cout << "Testing component i for NaNs: " << i << std::endl;
        amrex::Abort("Ax_new has NaNs in this component::advance_FDM_FFT()");
      }
    }
  }
}


inline void fdm_timestep(AmrexFFT& fft, cMF& psi,
                         MultiFab& Ax_new, MultiFab& phi, Gravity* gravity, Geometry& geom,
                         int level, Real h, Real dt_c, Real a_c,
                         Real dt_d, Real a_d, Real a_new, Real hbaroverm)
{
  //  *******************************************
  //  drift by dt_c
  //  *******************************************
  drift(fft, psi, Ax_new, geom.Domain(), dt_c, h, a_c, hbaroverm);
  Ax_new.FillBoundary(geom.periodicity());

  if (!dt_d)
    return;

  //  *******************************************
  //  re-calculate potential
  //  *******************************************
  int fill_interior = 0;
  int grav_n_grow = 1;
  gravity->solve_for_new_phi(level, phi,
                             gravity->get_grad_phi_curr(level),
                             fill_interior, grav_n_grow);

  //  *******************************************
  //  kick by dt_d
  //  *******************************************
  for (MFIter mfi(psi); mfi.isValid(); ++mfi)
  {
    auto p             = psi[mfi].array();
    auto const phiarr  = phi.const_array(mfi);
    amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
      Real arg = phiarr(i,j,k) * a_new / a_d / hbaroverm * dt_d;
      p(i,j,k) *= amrex::exp(GpuComplex<Real>(Real(0.0), arg));
    });
  }
}

inline void fdm_timestep_NG(AmrexFFT& fft, cMF& psi,
                            MultiFab& Ax_new, MultiFab& phi, Geometry& geom,
                            int level, Real h, Real dt_c, Real a_c,
                            Real dt_d, Real a_d, Real a_new, Real hbaroverm)
{
  //  *******************************************
  //  drift by dt_c
  //  *******************************************
  drift(fft, psi, Ax_new, geom.Domain(), dt_c, h, a_c, hbaroverm);
  Ax_new.FillBoundary(geom.periodicity());

  if (!dt_d)
    return;

  //  *******************************************
  //  kick by dt_d (phi fixed from initial setup — no gravity re-solve)
  //  *******************************************
  for (MFIter mfi(psi); mfi.isValid(); ++mfi)
  {
    auto p             = psi[mfi].array();
    auto const phiarr  = phi.const_array(mfi);
    amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
      Real arg = phiarr(i,j,k) * a_new / a_d / hbaroverm * dt_d;
      p(i,j,k) *= amrex::exp(GpuComplex<Real>(Real(0.0), arg));
    });
  }
}

inline void drift(AmrexFFT& fft, cMF& psi, MultiFab& Ax_new,
                  Box const& domain, Real dt, Real h, Real a_half, Real hbaroverm)
{
  if (!dt)
    return;

  const Real tpi = 2.0 * Real(4 * std::atan(1.0));
  const int  nx  = domain.length(0);
  const int  ny  = domain.length(1);
  const int  nz  = domain.length(2);
  // FFTW convention: forward then backward scales by N — divide here
  const Real norm = Real(1.0) / Real(nx * ny * nz);

  fft.forwardThenBackward(psi, psi,
    [=] AMREX_GPU_DEVICE (int gx, int gy, int gz, GpuComplex<Real>& val)
    {
      // wrap 0..N-1 indices to signed wavenumbers -N/2..N/2
      int ix = (gx > nx/2) ? gx - nx : gx;
      int iy = (gy > ny/2) ? gy - ny : gy;
      int iz = (gz > nz/2) ? gz - nz : gz;
      Real kx = tpi * Real(ix) / Real(nx);
      Real ky = tpi * Real(iy) / Real(ny);
      Real kz = tpi * Real(iz) / Real(nz);
      Real k2 = (kx*kx + ky*ky + kz*kz) / (h*h);
      val *= amrex::exp(GpuComplex<Real>(Real(0.0),
                        -hbaroverm * k2 / (a_half*a_half) / Real(2.0) * dt)) * norm;
    });

  // Update AxDens in Ax_new so FillBoundary and gravity solve see current density
  for (MFIter mfi(Ax_new); mfi.isValid(); ++mfi)
  {
    auto const p = psi[mfi].const_array();
    auto arr     = Ax_new[mfi].array();
    amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
      Real re = p(i,j,k).real();
      Real im = p(i,j,k).imag();
      arr(i,j,k,Nyx::AxDens) = re*re + im*im;
    });
  }
}
#endif
#endif // FDM
