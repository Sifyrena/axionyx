#ifdef AXIONYX
#include <AMReX_BLProfiler.H>
#include <Nyx.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

#ifdef AXCOMPLEX
void Nyx::axionyxEstimateTimeStep(Real& est_dt)
{
  const Real* dx = geom.CellSize();
  Real cur_time = state[Axion_Type].curTime();
  Real dt = string_dt*dx[0]/std::sqrt(msa*msa+string_massterm*pow(cur_time/string_time1,string_n+2.0)/string_time1/string_time1+12);
  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "...estdt from strings :  "<< dt <<'\n';
  est_dt = std::min(est_dt, dt);
}
#endif

#ifdef AXREAL
void Nyx::axionyxEstimateTimeStep(Real& est_dt)
{
  const Real* dx = geom.CellSize();
  Real cur_time = state[Axion_Type].curTime();
  Real dt = ax_dt*dx[0]/std::sqrt(pow(cur_time/ax_time1,ax_n)*pow(cur_time,3)+12);
  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "...estdt from axion only :  "<< dt <<'\n';
  est_dt = std::min(est_dt, dt);
}
#endif

#ifdef FDM
void Nyx::axionyxEstimateTimeStep(Real& est_dt)
{
  if (vonNeumann_dt >0 && (levelmethod[level]==FDlevel || levelmethod[level]==PSlevel) ){
    Real cur_time = state[Axion_Type].curTime();
    Real a = get_comoving_a(cur_time);
    const Real* dx = geom.CellSize();
    Real dt = dx[0]*dx[0]*a*a/6.0/hbaroverm;
#ifdef GRAVITY
    const MultiFab& phi = get_new_data(PhiGrav_Type);
    Real phi_max = std::abs(phi.max(0)-phi.min(0));
    dt = std::min(dt,hbaroverm/phi_max);
#endif
    if (levelmethod[level]==PSlevel)
      dt *= vonNeumann_dt;
    else if (vonNeumann_dt<1.0)
      dt *= vonNeumann_dt;
    if (verbose && ParallelDescriptor::IOProcessor())
      std::cout << "...estdt from von Neumann stability :  "<< dt <<'\n';
    est_dt = std::min(est_dt, dt);
  }
}
#endif

#endif
