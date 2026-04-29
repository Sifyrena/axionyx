#ifdef AXREAL
 
#include "Nyx.H"
#include <AMReX_MultiFab.H>
#include "compute_laplace_operator.H"

using namespace amrex;

void
Nyx::advance_AxReal_FD (Real time, Real dt)
{
    if (verbose && ParallelDescriptor::IOProcessor() ){
	std::cout << "Advancing the axions at level " << level <<  "...\n";
    }
    const Real* dx      = geom.CellSize();
    const Real invdeltasq  = 1.0 / dx[0] / dx[0];
    const Real dt_half = 0.5*dt;

    MultiFab&  Ax_old = get_old_data(Axion_Type);
    MultiFab&  Ax_new = get_new_data(Axion_Type);

#ifdef DEBUG
    if (Ax_old.contains_nan(0, Ax_old.nComp(), 0))
      {
        for (int i = 0; i < Ax_old.nComp(); i++)
          {
            if (Ax_old.contains_nan(i,1,0))
              {
		std::cout << "Testing component i for NaNs: " << i << std::endl;
		amrex::Abort("Ax_old has NaNs in this component::advance_AxReal_FD()");
              }
          }
      }
#endif

    kick_AxReal_FD(time, dt_half, Ax_old, Ax_new, invdeltasq, neighbors);
    drift_AxReal_FD(dt, Ax_old, Ax_new);
    kick_AxReal_FD(time+dt, dt_half, Ax_new, Ax_new, invdeltasq, neighbors);

#ifdef DEBUG
    if (Ax_new.contains_nan(0, Ax_new.nComp(), 0))
      {
        for (int i = 0; i < Ax_new.nComp(); i++)
          {
            if (Ax_new.contains_nan(i,1,0))
              {
		std::cout << "Testing component i for NaNs: " << i << std::endl;
		amrex::Abort("Ax_new has NaNs in this component::advance_AxReal_FD()");
              }
          }
      }
#endif
}

void 
Nyx::kick_AxReal_FD (Real time, Real dt, MultiFab&  mf_old, MultiFab&  mf_new, const Real invdeltasq, int neighbors)
{
  for (FillPatchIterator 
	 fpi(*this, mf_old, neighbors, time, Axion_Type, Nyx::Psi, 1);
       fpi.isValid(); ++fpi)
    {
      const Box& bx  = fpi.validbox();
      auto const arr_in   = fpi().array();
      auto const arr_old  = mf_old[fpi].array();
      auto const arr_new  = mf_new[fpi].array();
      
      ParallelFor(bx,
		  [=] AMREX_GPU_DEVICE (int i, int j, int k)
		  {
		    arr_new(i,j,k,Nyx::PsiPrime) = arr_old(i,j,k,Nyx::PsiPrime)+dt
		      *(compute_laplace_operator(arr_in,i,j,k,Nyx::Psi,invdeltasq,neighbors)
			-pow(time/Nyx::ax_time1,Nyx::ax_n)*pow(time,3)*std::sin(arr_old(i,j,k,Nyx::Psi)/time));
		  });
    }
}

void 
Nyx::drift_AxReal_FD (Real dt, MultiFab&  mf_old, MultiFab&  mf_new)
{
  for (MFIter mfi(mf_old,TilingIfNotGPU()); mfi.isValid(); ++mfi){
    auto const arr_old  = mf_old.array(mfi);
    auto const arr_new = mf_new.array(mfi);
    const Box& bx = mfi.tilebox();
    ParallelFor(bx,
		[=] AMREX_GPU_DEVICE (int i, int j, int k)
		{
		  arr_new(i,j,k,Nyx::Psi) = arr_old(i,j,k,Nyx::Psi) + arr_new(i,j,k,Nyx::PsiPrime)*dt;		  
		});
  }
}
#endif
