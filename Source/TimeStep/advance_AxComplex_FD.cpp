#ifdef AXCOMPLEX
 
#include "Nyx.H"
#include <AMReX_MultiFab.H>
#include "compute_laplace_operator.H"

using namespace amrex;

void
Nyx::advance_AxComplex_FD (Real time, Real dt)
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
		amrex::Abort("Ax_old has NaNs in this component::advance_AxComplex_FD()");
              }
          }
      }
#endif

    if(jaxions_ic){
      advance_initial_conditions(time+dt, Ax_old, Ax_new, neighbors);
    }else{
      rescale_field(Ax_old);
      kick_AxComplex_FD(time, dt_half, Ax_old, Ax_new, invdeltasq, neighbors);
      project_velocity(Ax_old, Ax_new);
      drift_AxComplex_FD(dt, Ax_old, Ax_new);
      project_velocity(Ax_new, Ax_new);
      kick_AxComplex_FD(time+dt, dt_half, Ax_new, Ax_new, invdeltasq, neighbors);
      rescale_field(Ax_new);
    }

#ifdef DEBUG
    if (Ax_new.contains_nan(0, Ax_new.nComp(), 0))
      {
        for (int i = 0; i < Ax_new.nComp(); i++)
          {
            if (Ax_new.contains_nan(i,1,0))
              {
		std::cout << "Testing component i for NaNs: " << i << std::endl;
		amrex::Abort("Ax_new has NaNs in this component::advance_AxComplex_FD()");
              }
          }
      }
#endif
}


void Nyx::kick_AxComplex_FD(Real time, Real dt, MultiFab&  mf_old, MultiFab&  mf_new, const Real invdeltasq, int neighbors)
{
  for (FillPatchIterator 
	 fpi(*this, mf_old, neighbors, time, Axion_Type, 0, 2);
       fpi.isValid(); ++fpi)
    {
      const Box& bx  = fpi.validbox();
      auto const arr_in   = fpi().array();
      auto const arr_old  = mf_old[fpi].array();
      auto const arr_new  = mf_new[fpi].array();
      Real lprs = 1.0;
      if(Nyx::prsstring)
	lprs = pow(Nyx::msa/time,2)/2.0;
      else{
	if(Nyx::string_stop_time<=0.0) amrex::Abort("prsstring needs string_stop_time>0");
	lprs = pow(Nyx::msa/Nyx::string_stop_time,2)/2.0;
      }
      
      if(Nyx::fixmodulus[level])
	ParallelFor(bx,
		  [=] AMREX_GPU_DEVICE (int i, int j, int k)
		  {
		    arr_new(i,j,k,Nyx::AxvRe) = arr_old(i,j,k,Nyx::AxvRe)*(time-dt)/(time+dt)-time*dt/(time+dt)
		      *(arr_old(i,j,k,Nyx::AxRe)*compute_laplace_operator(arr_in,i,j,k,Nyx::AxIm,invdeltasq,neighbors)
			-arr_old(i,j,k,Nyx::AxIm)*compute_laplace_operator(arr_in,i,j,k,Nyx::AxRe,invdeltasq,neighbors)		    
			-string_massterm*pow(time/Nyx::string_time1,Nyx::string_n+2.0)/Nyx::string_time1/Nyx::string_time1
                        *arr_in(i,j,k,Nyx::AxIm))*arr_old(i,j,k,Nyx::AxIm);

		    arr_new(i,j,k,Nyx::AxvIm) = arr_old(i,j,k,Nyx::AxvIm)*(time-dt)/(time+dt)+time*dt/(time+dt)
		      *(arr_old(i,j,k,Nyx::AxRe)*compute_laplace_operator(arr_in,i,j,k,Nyx::AxIm,invdeltasq,neighbors)
			-arr_old(i,j,k,Nyx::AxIm)*compute_laplace_operator(arr_in,i,j,k,Nyx::AxRe,invdeltasq,neighbors)
			-string_massterm*pow(time/Nyx::string_time1,Nyx::string_n+2.0)/Nyx::string_time1/Nyx::string_time1
                        *arr_in(i,j,k,Nyx::AxIm))*arr_old(i,j,k,Nyx::AxRe);
		  });
      else
	ParallelFor(bx,
		  [=] AMREX_GPU_DEVICE (int i, int j, int k)
		  {
		    arr_new(i,j,k,Nyx::AxvRe) = arr_old(i,j,k,Nyx::AxvRe)*(time-dt)/(time+dt)+time*dt/(time+dt)
		      *(compute_laplace_operator(arr_in,i,j,k,Nyx::AxRe,invdeltasq,neighbors)
			-lprs*time*time*(arr_in(i,j,k,Nyx::AxRe)*arr_in(i,j,k,Nyx::AxRe)+arr_in(i,j,k,Nyx::AxIm)*arr_in(i,j,k,Nyx::AxIm)-1)*arr_in(i,j,k,Nyx::AxRe)
			+2.0*string_massterm*pow(time/Nyx::string_time1,Nyx::string_n+2.0)/Nyx::string_time1/Nyx::string_time1*arr_in(i,j,k,Nyx::AxIm)*arr_in(i,j,k,Nyx::AxIm)
			*pow(arr_in(i,j,k,Nyx::AxRe)*arr_in(i,j,k,Nyx::AxRe)+arr_in(i,j,k,Nyx::AxIm)*arr_in(i,j,k,Nyx::AxIm),-1.5));
		    
		    arr_new(i,j,k,Nyx::AxvIm) = arr_old(i,j,k,Nyx::AxvIm)*(time-dt)/(time+dt)+time*dt/(time+dt)
		      *(compute_laplace_operator(arr_in,i,j,k,Nyx::AxIm,invdeltasq,neighbors)
			-lprs*time*time*(arr_in(i,j,k,Nyx::AxRe)*arr_in(i,j,k,Nyx::AxRe)+arr_in(i,j,k,Nyx::AxIm)*arr_in(i,j,k,Nyx::AxIm)-1)*arr_in(i,j,k,Nyx::AxIm)
			-2.0*string_massterm*pow(time/Nyx::string_time1,Nyx::string_n+2.0)/Nyx::string_time1/Nyx::string_time1*arr_in(i,j,k,Nyx::AxRe)*arr_in(i,j,k,Nyx::AxIm)
			*pow(arr_in(i,j,k,Nyx::AxRe)*arr_in(i,j,k,Nyx::AxRe)+arr_in(i,j,k,Nyx::AxIm)*arr_in(i,j,k,Nyx::AxIm),-1.5));
		  });
      }
}


void Nyx::drift_AxComplex_FD(Real dt, MultiFab&  mf_old, MultiFab&  mf_new)
{
  for (MFIter mfi(mf_new,TilingIfNotGPU()); mfi.isValid(); ++mfi){
    auto const arr_old = mf_old.array(mfi);
    auto const arr_new = mf_new.array(mfi);
    const Box& bx = mfi.tilebox();
    if(Nyx::fixmodulus[level])
      ParallelFor(bx,
		  [=] AMREX_GPU_DEVICE (int i, int j, int k)
		  {
		    Real phase = (arr_old(i,j,k,Nyx::AxRe)*arr_new(i,j,k,Nyx::AxvIm)-arr_old(i,j,k,Nyx::AxIm)*arr_new(i,j,k,Nyx::AxvRe))*dt;
		    arr_new(i,j,k,Nyx::AxRe) = arr_old(i,j,k,Nyx::AxRe)*std::cos(phase)-arr_old(i,j,k,Nyx::AxIm)*std::sin(phase);
		    arr_new(i,j,k,Nyx::AxIm) = arr_old(i,j,k,Nyx::AxRe)*std::sin(phase)+arr_old(i,j,k,Nyx::AxIm)*std::cos(phase);
		  });
    else
      ParallelFor(bx,
		  [=] AMREX_GPU_DEVICE (int i, int j, int k)
		  {
		    arr_new(i,j,k,Nyx::AxRe) = arr_old(i,j,k,Nyx::AxRe) + arr_new(i,j,k,Nyx::AxvRe)*dt;
		    arr_new(i,j,k,Nyx::AxIm) = arr_old(i,j,k,Nyx::AxIm) + arr_new(i,j,k,Nyx::AxvIm)*dt;
		  });
  }
}


void Nyx::rescale_field(MultiFab& mf)
{

  if(!Nyx::fixmodulus[level])
    return;

  for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi){
    auto const arr = mf.array(mfi);
    const Box& bx = mfi.tilebox();
    ParallelFor(bx,
		[=] AMREX_GPU_DEVICE (int i, int j, int k)
		{
		  Real amp = std::sqrt(arr(i,j,k,Nyx::AxRe)*arr(i,j,k,Nyx::AxRe)+arr(i,j,k,Nyx::AxIm)*arr(i,j,k,Nyx::AxIm));
		  arr(i,j,k,Nyx::AxRe) /= amp;
		  arr(i,j,k,Nyx::AxIm) /= amp;
		});
  }
}


void Nyx::project_velocity(MultiFab& mf_old, MultiFab& mf_new)
{

  if(!Nyx::fixmodulus[level])
    return;

  for (MFIter mfi(mf_old,TilingIfNotGPU()); mfi.isValid(); ++mfi){
    auto const arr_old = mf_old.array(mfi);
    auto const arr_new = mf_new.array(mfi);
    const Box& bx = mfi.tilebox();
    ParallelFor(bx,
		[=] AMREX_GPU_DEVICE (int i, int j, int k)
		{
		  arr_new(i,j,k,Nyx::AxvRe) = -(arr_old(i,j,k,Nyx::AxRe)*arr_new(i,j,k,Nyx::AxvIm)-arr_old(i,j,k,Nyx::AxIm)*arr_new(i,j,k,Nyx::AxvRe))*arr_old(i,j,k,Nyx::AxIm);
		  arr_new(i,j,k,Nyx::AxvIm) =  (arr_old(i,j,k,Nyx::AxRe)*arr_new(i,j,k,Nyx::AxvIm)-arr_old(i,j,k,Nyx::AxIm)*arr_new(i,j,k,Nyx::AxvRe))*arr_old(i,j,k,Nyx::AxRe);
		});
  }
}



void Nyx::advance_initial_conditions(Real time, MultiFab&  mf_old, MultiFab&  mf_new, int neighbors)
{
  for (FillPatchIterator 
	 fpi(*this, mf_old, neighbors, time, Axion_Type,   0, Nyx::NUM_AX);
       fpi.isValid(); ++fpi)
    {
      const Box& bx  = fpi.validbox();
      auto const arr_in   = fpi().array();
      auto const arr_out  = mf_new[fpi].array();
      Real diff = 0.01;
      
      ParallelFor(bx,
		  [=] AMREX_GPU_DEVICE (int i, int j, int k)
		  {
		    arr_out(i,j,k,Nyx::AxRe) = arr_in(i,j,k,Nyx::AxRe)+compute_laplace_operator(arr_in,i,j,k,Nyx::AxRe,diff,neighbors);
		    arr_out(i,j,k,Nyx::AxIm) = arr_in(i,j,k,Nyx::AxIm)+compute_laplace_operator(arr_in,i,j,k,Nyx::AxIm,diff,neighbors);
		    
		    Real phase = atan2(arr_out(i,j,k,Nyx::AxIm),arr_out(i,j,k,Nyx::AxRe));
		    arr_out(i,j,k,Nyx::AxvRe ) = cos(phase);
		    arr_out(i,j,k,Nyx::AxvIm ) = sin(phase);
		    arr_out(i,j,k,Nyx::AxRe  ) = time*arr_out(i,j,k,Nyx::AxvRe);
		    arr_out(i,j,k,Nyx::AxIm  ) = time*arr_out(i,j,k,Nyx::AxvIm);
		  }); 
    }
}
#endif
