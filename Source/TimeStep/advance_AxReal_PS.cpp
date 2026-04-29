#ifdef AXREAL
#ifdef BL_USE_MPI
#include <AMReX_BLProfiler.H>
#include <Nyx.H>
#include <AMReX_MultiFab.H>
#include <Distribution.H>
#include <AlignedAllocator.h>
#include <Dfft.H>

#include <LLUTs.h>

#include <string>

#define ALIGN 16

using namespace amrex;

void fdm_timestep(hacc::Dfft &dfft, MultiFab &Ax_new, MultiFab &Ax_old, Geometry &geom,
		  int const level, int const gridsize, Real const h, Real const dt_c,  
		  Real const dt_d, 
		  Real const time, std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> >* a,std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> >* b);

void Nyx::advance_AxReal_PS(amrex::Real time,
			 amrex::Real dt)
{

  BL_PROFILE("Nyx::advance_AxReal_PS()");

    // *****************************************
    //define constants
    // *****************************************
    const Real h = geom.CellSize(0);

    Real w0,w1,w2,w3;
    Real c1,c2,c3,c4,c5,c6,c7,c8;
    Real d1,d2,d3,d4,d5,d6,d7,d8;

    //******************************************
    //define weights for time steps (see Levkov et. al. 2018)
    //******************************************

    if(order==6){
      w1 = -1.17767998417887;
      w2 = 0.235573213359359;
      w3 = 0.784513610477560;
      w0 = 1.0-2.0*(w1+w2+w3);

      c1 = w3/2.0*dt;
      c2 = (w2+w3)/2.0*dt;
      c3 = (w1+w2)/2.0*dt;
      c4 = (w0+w1)/2.0*dt;
      c5 = (w0+w1)/2.0*dt;
      c6 = (w1+w2)/2.0*dt;
      c7 = (w2+w3)/2.0*dt;
      c8 = w3/2.0*dt;

      d1 = w3*dt;
      d2 = w2*dt;
      d3 = w1*dt;
      d4 = w0*dt;
      d5 = w1*dt;
      d6 = w2*dt;
      d7 = w3*dt;
      d8 = 0.0;

    }else if(order==2){

    c1 = 0.5*dt;
    c2 = 0.5*dt;

    d1 = dt;
    d2 = 0.0;

    }else
      amrex::Error("Order of algorithm not implemented!");

    // *****************************************
    // Get Axion State
    // *****************************************
    MultiFab& Ax_old = get_old_data(Axion_Type);
    MultiFab& Ax_new = get_new_data(Axion_Type);
    const BoxArray& ba = Ax_old.boxArray();
    const DistributionMapping& dm = Ax_old.DistributionMap();
    
    if (Ax_old.contains_nan(0, Ax_old.nComp(), 0))                                                                                                                                                      
      {                                                                                                                                                                                                   
	for (int i = 0; i < Ax_old.nComp(); i++)                                                                                                                                                        
	  {                                                                                                                                                                                               
	    if (Ax_old.contains_nan(i,1,0))                                                                                                                                                             
	      {                                                                                                                                                                                           
		std::cout << "Testing component i for NaNs: " << i << std::endl;                                                                                                                        
		amrex::Abort("Ax_old has NaNs in this component::advance_AxReal_FFT()");
	      }                                                                                                                                                                                           
	  }                                                                                                                                                                                               
      } 
    
    // *****************************************
    // We assume that all grids have the same size hence
    // we have the same nx,ny,nz on all ranks
    // *****************************************
    int nx = ba[0].size()[0];
    int ny = ba[0].size()[1];
    int nz = ba[0].size()[2];
    int gridsize = nx*ny*nz;

    Box domain(geom.Domain());

    int nbx = domain.length(0) / nx;
    int nby = domain.length(1) / ny;
    int nbz = domain.length(2) / nz;
    int nboxes = nbx * nby * nbz;

    if (nboxes != ba.size())
        amrex::Error("NBOXES NOT COMPUTED CORRECTLY");

    // *****************************************
    // This unfortunately seems neccessary (cf. amrex/Src/Extern/SWFFT/README)
    // *****************************************
    if(ParallelDescriptor::NProcs() < nboxes)
        amrex::Error("Number of MPI ranks has to be larger or equal to the number of root level grids!");

    Vector<int> rank_mapping;
    rank_mapping.resize(nboxes);

    for (int ib = 0; ib < nboxes; ++ib)
    {
        int i = ba[ib].smallEnd(0) / nx;
        int j = ba[ib].smallEnd(1) / ny;
        int k = ba[ib].smallEnd(2) / nz;
        int local_index = i*nbx*nby + j*nbx + k;

        rank_mapping[local_index] = dm[ib];

        if (verbose)
            amrex::Print() << "LOADING RANK NUMBER " << dm[ib] << " FOR GRID NUMBER " << ib
                 << " WHICH IS LOCAL NUMBER " << local_index << std::endl;
    }
    // *****************************************
    // Assume for now that nx = ny = nz
    // *****************************************
    int Ndims[3] = { nbz, nby, nbx };
    int     n[3] = {domain.length(2), domain.length(1), domain.length(0)};
    hacc::Distribution d(MPI_COMM_WORLD,n,Ndims,&rank_mapping[0]);
    hacc::Dfft dfft(d);

    //  *******************************************
    //  prepare fft
    //  *******************************************

    std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > a;
    std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > b;
    a.resize(gridsize);
    b.resize(gridsize);
    dfft.makePlans(&a[0],&b[0],&a[0],&b[0]);

    for (MFIter mfi(Ax_old,false); mfi.isValid(); ++mfi){
      auto const arr = Ax_old.array(mfi);
      const Box& bx = mfi.validbox();
      const Dim3 lo = amrex::lbound(bx);
      const Dim3 hi = amrex::ubound(bx);
      const Dim3 w ={hi.x-lo.x+1,hi.y-lo.y+1,hi.z-lo.z+1};
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(size_t k=0; k<(size_t)w.z; k++) {
        for(size_t j=0; j<(size_t)w.y; j++) {
	  AMREX_PRAGMA_SIMD
	  for(size_t i=0; i<(size_t)w.x; i++) {
            size_t local_indx_threaded = (size_t)w.y*(size_t)w.z*i+(size_t)w.z*j+k;
            complex_t temp(arr(i+lo.x,j+lo.y,k+lo.z,Nyx::Psi),0.0);
            a[local_indx_threaded] = temp;
          }}}
    }


    //  *******************************************
    //  higher order time steps
    //  *******************************************
    if(order==6){
      fdm_timestep(dfft, Ax_new, Ax_old, geom, level, gridsize, h, c1, d1, time, &a, &b);
      fdm_timestep(dfft, Ax_new, Ax_old, geom, level, gridsize, h, c2, d2, time, &a, &b);
      fdm_timestep(dfft, Ax_new, Ax_old, geom, level, gridsize, h, c3, d3, time, &a, &b);
      fdm_timestep(dfft, Ax_new, Ax_old, geom, level, gridsize, h, c4, d4, time, &a, &b);
      fdm_timestep(dfft, Ax_new, Ax_old, geom, level, gridsize, h, c5, d5, time, &a, &b);
      fdm_timestep(dfft, Ax_new, Ax_old, geom, level, gridsize, h, c6, d6, time, &a, &b);
      fdm_timestep(dfft, Ax_new, Ax_old, geom, level, gridsize, h, c7, d7, time, &a, &b);
      fdm_timestep(dfft, Ax_new, Ax_old, geom, level, gridsize, h, c8, d8, time, &a, &b);
    }else if(order==2){
      fdm_timestep(dfft, Ax_new, Ax_old, geom, level, gridsize, h, c1, d1, time, &a, &b);
      fdm_timestep(dfft, Ax_new, Ax_old, geom, level, gridsize, h, c2, d2, time+dt, &a, &b);
    }

    Ax_new.FillBoundary(geom.periodicity());

    if (Ax_new.contains_nan(0, Ax_new.nComp(), 0))                                                                                                                                                      
      {                                                                                                                                                                                                   
	for (int i = 0; i < Ax_new.nComp(); i++)                                                                                                                                                        
	  {                                                                                                                                                                                               
	    if (Ax_new.contains_nan(i,1,0))                                                                                                                                                             
	      {                                                                                                                                                                                           
		std::cout << "Testing component i for NaNs: " << i << std::endl;                                                                                                                        
		amrex::Abort("Ax_new has NaNs in this component::advance_AxReal_FFT()");
	      }                                                                                                                                                                                           
	  }                                                                                                                                                                                               
      } 

}

inline void fdm_timestep(hacc::Dfft &dfft, MultiFab &Ax_new, MultiFab &Ax_old, Geometry &geom,
			 int const level, int const gridsize, Real const h, Real const dt_c,  
			 Real const dt_d, Real const time, 
			 std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> >* a,std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> >* b)
{
  const Real pi = 4 * std::atan(1.0);
  const Real tpi = 2 * pi;
  const Real hsq = h*h;
  const std::complex<double> imagi(0.0,1.0);
  const int *self = dfft.self_kspace();
  const int *local_ng = dfft.local_ng_kspace();
  const int *global_ng = dfft.global_ng();
  size_t local_indx = 0;

  dfft.forward(&((*a)[0]));    
#ifdef _OPENMP
#pragma omp parallel for       
#endif
    for(size_t i=0; i<(size_t)local_ng[0]; i++) {
      int global_i = local_ng[0]*self[0] + i;
      if (global_i > global_ng[0]/2.){
	global_i = global_i - global_ng[0];
      }
	
	
      for(size_t j=0; j<(size_t)local_ng[1]; j++) {
	int global_j = local_ng[1]*self[1] + j;
	if (global_j > global_ng[1]/2.){
	  global_j = global_j - global_ng[1];
	}
	  
        AMREX_PRAGMA_SIMD	  
	for(size_t k=0; k<(size_t)local_ng[2]; k++) {
	  int global_k = local_ng[2]*self[2] + k;
          size_t local_indx_threaded = (size_t)local_ng[1]*(size_t)local_ng[2]*i+(size_t)local_ng[2]*j+k;
	  if (global_k > global_ng[2]/2.){
	    global_k = global_k - global_ng[2];
	  }
	    
	  double kx = tpi * double(global_i)/double(global_ng[0]);
	  double ky = tpi * double(global_j)/double(global_ng[1]);
	  double kz = tpi * double(global_k)/double(global_ng[2]);
	  double k2 = (kx*kx + ky*ky + kz*kz)/h/h;
	
	  // (*a)[local_indx_threaded] *= std::exp(-k2*dt_c);
	  (*a)[local_indx_threaded] *= -k2*dt_c;
          (*a)[local_indx_threaded] /= dfft.global_size();	
	}
      }
    }
  dfft.backward(&((*a)[0]));
  
  for (MFIter mfi(Ax_old,false); mfi.isValid(); ++mfi){
    auto const axnew = Ax_new[mfi].array();
    auto const axold = Ax_old[mfi].array();
    const Box& bx = mfi.validbox();
    const Dim3 lo = amrex::lbound(bx);
    const Dim3 hi = amrex::ubound(bx);
    const Dim3 w ={hi.x-lo.x+1,hi.y-lo.y+1,hi.z-lo.z+1};

    amrex::ParallelFor(bx,
		       [=] AMREX_GPU_DEVICE (int i, int j, int k)
		       {
			 size_t local_indx_threaded = (size_t)w.y*(size_t)w.z*(i-lo.x)+(size_t)w.z*(j-lo.y)+(k-lo.z);
			 complex_t temp = (*a)[local_indx_threaded];
			 axnew(i,j,k,Nyx::PsiPrime) += std::real(temp)-pow(time/Nyx::ax_time1,Nyx::ax_n)*pow(time,3)*std::sin(axnew(i,j,k,Nyx::Psi)/time)*dt_c;
			 axnew(i,j,k,Nyx::Psi)  += axnew(i,j,k,Nyx::PsiPrime)*dt_d;
			 (*a)[local_indx_threaded] = complex_t(axnew(i,j,k,Nyx::Psi),0.0);
		       });
  }
}

#endif
#endif
