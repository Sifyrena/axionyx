#ifdef INFGRAV
#ifdef BL_USE_MPI
#include <AMReX_BLProfiler.H>
#include <Nyx.H>
#include <AMReX_MultiFab.H>
#ifdef GRAVITY
#include <Gravity.H>
#endif
#include <constants_cosmo.H>
#include <Distribution.H>
#include <AlignedAllocator.h>
#include <Dfft.H>

#include <string>

#define ALIGN 16

using namespace amrex;

inline Real compute_Oijlm(int i, int j, int l, int m, const Real arr[3]){
  Real P_ij, P_il, P_jm, P_lm;
  if (i==j)
    P_ij = 1.0 - arr[i-1]*arr[j-1];
  else
    P_ij = 0.0 - arr[i-1]*arr[j-1];

  if (i==l)
    P_il = 1.0 - arr[i-1]*arr[l-1];
  else
    P_il = 0.0 - arr[i-1]*arr[l-1];

  if (j==m)
    P_jm = 1.0 - arr[j-1]*arr[m-1];
  else
    P_jm = 0.0 - arr[j-1]*arr[m-1];

  if (l==m)
    P_lm = 1.0 - arr[l-1]*arr[m-1];
  else
    P_lm = 0.0 - arr[l-1]*arr[m-1];

  return P_il * P_jm - 0.5*P_ij * P_lm;
  
}

inline Real compute_tau(Real t){
  Real t_init = 1.26813e-05;
  Real a_init = 0.00247875;
  return 3.0*pow(t_init,2.0/3.0)/a_init * (pow(t+t_init,1.0/3.0) - pow(t_init,1.0/3.0)); // pow(t+t_init,1.0/3.0): t_init required here because t starts at 0 while t_start=t_init in runlog. 
  
}

inline Real compute_a(Real t){
  Real t_init = 1.26813e-05;
  Real a_init = 0.00247875;
  return a_init*std::pow((t+t_init)/t_init,2./3.);
}

hacc::Distribution Nyx::initialize_FFT(){

  BL_PROFILE("Nyx::initialize_FFT()");
  
  MultiFab& Ax_new = get_level(level).get_new_data(Axion_Type); //get_level(0)...

  const BoxArray& ba = Ax_new.boxArray();
  const DistributionMapping& dm = Ax_new.DistributionMap();

  int nx = ba[0].size()[0];
  int ny = ba[0].size()[1];
  int nz = ba[0].size()[2];
  int gridsize = nx*ny*nz;

  Box domain(get_level(0).geom.Domain());

  int nbx = domain.length(0) / nx;
  int nby = domain.length(1) / ny;
  int nbz = domain.length(2) / nz;
  int nboxes = nbx * nby * nbz;

  if (nboxes != ba.size())
    amrex::Error("NBOXES NOT COMPUTED CORRECTLY");

  // *****************************************
  // This unfortunately seems neccessary (cf. amrex/Src/Extern/SWFFT/README)
  // *****************************************
  
  if(ParallelDescriptor::NProcs() < nboxes){
    //Print() << "Nboxes" << nboxes << std::endl;
    //Print()<< "nprocs" << ParallelDescriptor::NProcs() << std::endl;
    amrex::Error("Number of MPI ranks has to be larger or equal to the number of root level grids!");
  }

  Vector<int> rank_mapping;
  rank_mapping.resize(nboxes);

  for (int ib = 0; ib < nboxes; ++ib)
    {
      int i = ba[ib].smallEnd(0) / nx;
      int j = ba[ib].smallEnd(1) / ny;
      int k = ba[ib].smallEnd(2) / nz;
      int local_index = i*nbx*nby + j*nbx + k;

      rank_mapping[local_index] = dm[ib];

      //if (verbose)
	// amrex::Print() << "LOADING RANK NUMBER " << dm[ib] << " FOR GRID NUMBER " << ib
	// 	       << " WHICH IS LOCAL NUMBER " << local_index << std::endl;
    }
  // *****************************************
  // Assume for now that nx = ny = nz
  // *****************************************
  int Ndims[3] = { nbz, nby, nbx };
  int     n[3] = {domain.length(2), domain.length(1), domain.length(0)};
  hacc::Distribution d(MPI_COMM_WORLD,n,Ndims,&rank_mapping[0]);
  //hacc::Dfft dfft(d);

  return d;

}

void Nyx::fill_tensor_TT(){

  BL_PROFILE("Nyx::fill_tensor_TT()");

  MultiFab& TRe = get_level(level).get_new_data(TensorRe_Type); // get_level(0)...                                                                                                                               
  MultiFab& TIm = get_level(level).get_new_data(TensorIm_Type); // get_level(0)...                                                                                                                               
 
  std::unique_ptr<MultiFab> vel_x = get_level(level).derive("fdm_particle_x_velocity", state[Axion_Type].curTime(), 0);
  std::unique_ptr<MultiFab> vel_y = get_level(level).derive("fdm_particle_y_velocity", state[Axion_Type].curTime(), 0);
  std::unique_ptr<MultiFab> vel_z = get_level(level).derive("fdm_particle_z_velocity", state[Axion_Type].curTime(), 0); // in units of l_u/t_u ["km/s"]                                                           

  MultiFab& Ax_new = get_level(level).get_new_data(Axion_Type); //get_level(0)...                                                                                                                                 
  MultiFab& IcosRe = get_level(level).get_new_data(CosRe_Type);
  MultiFab& IcosIm = get_level(level).get_new_data(CosIm_Type);
  MultiFab& IsinRe = get_level(level).get_new_data(SinRe_Type);
  MultiFab& IsinIm = get_level(level).get_new_data(SinIm_Type);

  const BoxArray& ba = Ax_new.boxArray();
  const DistributionMapping& dm = Ax_new.DistributionMap();

  int nx = ba[0].size()[0];
  int ny = ba[0].size()[1];
  int nz = ba[0].size()[2];
  int gridsize = nx*ny*nz;

  const Real vec[3] = {1.0/std::sqrt(2.0), 1.0/std::sqrt(2.0), 0}; // six directions of k-vector: (1,1,0), (1,0,1), (0,1,1), (-1,1,0), (-1,0,1), (0,-1,1)
  const Real h = get_level(level).geom.CellSize(0);  
  hacc::Distribution d = initialize_FFT();
  hacc::Dfft dfft1(d); 
  hacc::Dfft dfft2(d); 
  hacc::Dfft dfft3(d); 
  hacc::Dfft dfft4(d); 
  hacc::Dfft dfft5(d); 
  hacc::Dfft dfft6(d); 

  const int *self = dfft1.self_kspace();
  const int *local_ng = dfft1.local_ng_kspace(); 
  const int *global_ng = dfft1.global_ng(); 

  std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > a_11, a_12, a_13, a_22, a_23, a_33;
  std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > b_11, b_12, b_13, b_22, b_23, b_33;
  a_11.resize(gridsize);
  a_12.resize(gridsize);
  a_13.resize(gridsize);
  a_22.resize(gridsize);
  a_23.resize(gridsize);
  a_33.resize(gridsize);
  b_11.resize(gridsize);
  b_12.resize(gridsize);
  b_13.resize(gridsize);
  b_22.resize(gridsize);
  b_23.resize(gridsize);
  b_33.resize(gridsize);

  dfft1.makePlans(&a_11[0],&b_11[0],&a_11[0],&b_11[0]);
  dfft2.makePlans(&a_12[0],&b_12[0],&a_12[0],&b_12[0]);
  dfft3.makePlans(&a_13[0],&b_13[0],&a_13[0],&b_13[0]);
  dfft4.makePlans(&a_22[0],&b_22[0],&a_22[0],&b_22[0]);
  dfft5.makePlans(&a_23[0],&b_23[0],&a_23[0],&b_23[0]);
  dfft6.makePlans(&a_33[0],&b_33[0],&a_33[0],&b_33[0]);

  for (MFIter mfi(Ax_new,false); mfi.isValid(); ++mfi){

    const Box& bx = mfi.validbox();

    Array4<Real> const& fab_ax = Ax_new.array(mfi);
    Array4<Real> const& fab_vx = (*vel_x).array(mfi);
    Array4<Real> const& fab_vy = (*vel_y).array(mfi);
    Array4<Real> const& fab_vz = (*vel_z).array(mfi);

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
            //amrex::Print() << "INDEX THREADED: " << local_indx_threaded << std::endl; // running from 0-1023. Makes sense because of tiling: 16 x 16 x 4 [32 root grid, 8 grids with size 16 x 16 x 16]      

            complex_t temp_11(fab_ax(i+lo.x,j+lo.y,k+lo.z,Nyx::AxDens)*fab_vx(i+lo.x,j+lo.y,k+lo.z)*fab_vx(i+lo.x,j+lo.y,k+lo.z),0); // include ,0 in fab_vx etc. ? 
            complex_t temp_12(fab_ax(i+lo.x,j+lo.y,k+lo.z,Nyx::AxDens)*fab_vx(i+lo.x,j+lo.y,k+lo.z)*fab_vy(i+lo.x,j+lo.y,k+lo.z),0);
            complex_t temp_13(fab_ax(i+lo.x,j+lo.y,k+lo.z,Nyx::AxDens)*fab_vx(i+lo.x,j+lo.y,k+lo.z)*fab_vz(i+lo.x,j+lo.y,k+lo.z),0);
            complex_t temp_22(fab_ax(i+lo.x,j+lo.y,k+lo.z,Nyx::AxDens)*fab_vy(i+lo.x,j+lo.y,k+lo.z)*fab_vy(i+lo.x,j+lo.y,k+lo.z),0);
            complex_t temp_23(fab_ax(i+lo.x,j+lo.y,k+lo.z,Nyx::AxDens)*fab_vy(i+lo.x,j+lo.y,k+lo.z)*fab_vz(i+lo.x,j+lo.y,k+lo.z),0);
            complex_t temp_33(fab_ax(i+lo.x,j+lo.y,k+lo.z,Nyx::AxDens)*fab_vz(i+lo.x,j+lo.y,k+lo.z)*fab_vz(i+lo.x,j+lo.y,k+lo.z),0);

            a_11[local_indx_threaded] = temp_11;
            a_12[local_indx_threaded] = temp_12;
            a_13[local_indx_threaded] = temp_13;
            a_22[local_indx_threaded] = temp_22;
            a_23[local_indx_threaded] = temp_23;
            a_33[local_indx_threaded] = temp_33;
          }
      }
    }

    
  } // end MFIter

  dfft1.forward(&(a_11[0]));
  dfft2.forward(&(a_12[0]));
  dfft3.forward(&(a_13[0]));
  dfft4.forward(&(a_22[0]));
  dfft5.forward(&(a_23[0]));
  dfft6.forward(&(a_33[0]));

  d.redistribute_2_to_3(&a_11[0],&b_11[0],2);
  d.redistribute_2_to_3(&a_12[0],&b_12[0],2);
  d.redistribute_2_to_3(&a_13[0],&b_13[0],2);
  d.redistribute_2_to_3(&a_22[0],&b_22[0],2);
  d.redistribute_2_to_3(&a_23[0],&b_23[0],2);
  d.redistribute_2_to_3(&a_33[0],&b_33[0],2);

  for (MFIter mfi(Ax_new,false); mfi.isValid(); ++mfi){
    const Box& bx = mfi.validbox();//mfi.tilebox();                                                                                                                                                                

    Array4<Real> const& fab_tRe = TRe.array(mfi);
    Array4<Real> const& fab_tIm = TIm.array(mfi);
    //Array4<Real> const& fab_cos = Icos.array(mfi);                                                                                                                                                               
    //Array4<Real> const& fab_sin = Isin.array(mfi);                                                                                                                                                               
    const Dim3 lo = amrex::lbound(bx);
    const Dim3 hi = amrex::ubound(bx);
    const Dim3 w ={hi.x-lo.x+1,hi.y-lo.y+1,hi.z-lo.z+1};

    amrex::ParallelFor(bx,
		       [=] AMREX_GPU_DEVICE (int i, int j, int k)
		       {
			 size_t local_indx_threaded = (size_t)w.y*(size_t)w.z*(i-lo.x)+(size_t)w.z*(j-lo.y)+(k-lo.z);			 		       

			 // returns correct (initial) expressions for real and imag part of TT tensor components!! 
			 fab_tRe(i,j,k,Nyx::T_11) = compute_Oijlm(1,1,1,1,vec)*b_11[local_indx_threaded].real()*std::pow(h,3) + compute_Oijlm(1,1,1,2,vec)* b_12[local_indx_threaded].real()*std::pow(h,3)
			   + compute_Oijlm(1,1,1,3,vec)* b_13[local_indx_threaded].real()*std::pow(h,3) + compute_Oijlm(1,1,2,1,vec)* b_12[local_indx_threaded].real()*std::pow(h,3)
			   + compute_Oijlm(1,1,2,2,vec)* b_22[local_indx_threaded].real()*std::pow(h,3) + compute_Oijlm(1,1,2,3,vec)* b_23[local_indx_threaded].real()*std::pow(h,3)
			   + compute_Oijlm(1,1,3,1,vec)* b_13[local_indx_threaded].real()*std::pow(h,3) + compute_Oijlm(1,1,3,2,vec)* b_23[local_indx_threaded].real()*std::pow(h,3)
			   + compute_Oijlm(1,1,3,3,vec)* b_33[local_indx_threaded].real()*std::pow(h,3);

			 fab_tRe(i,j,k,Nyx::T_12) = compute_Oijlm(1,2,1,1,vec)*b_11[local_indx_threaded].real()*std::pow(h,3) + compute_Oijlm(1,2,1,2,vec)* b_12[local_indx_threaded].real()*std::pow(h,3)
			   + compute_Oijlm(1,2,1,3,vec)* b_13[local_indx_threaded].real()*std::pow(h,3) + compute_Oijlm(1,2,2,1,vec)* b_12[local_indx_threaded].real()*std::pow(h,3)
			   + compute_Oijlm(1,2,2,2,vec)* b_22[local_indx_threaded].real()*std::pow(h,3) + compute_Oijlm(1,2,2,3,vec)* b_23[local_indx_threaded].real()*std::pow(h,3)
			   + compute_Oijlm(1,2,3,1,vec)* b_13[local_indx_threaded].real()*std::pow(h,3) + compute_Oijlm(1,2,3,2,vec)* b_23[local_indx_threaded].real()*std::pow(h,3)
			   + compute_Oijlm(1,2,3,3,vec)* b_33[local_indx_threaded].real()*std::pow(h,3);

			 fab_tRe(i,j,k,Nyx::T_13) = compute_Oijlm(1,3,1,1,vec)*b_11[local_indx_threaded].real()*std::pow(h,3) + compute_Oijlm(1,3,1,2,vec)* b_12[local_indx_threaded].real()*std::pow(h,3)
			   + compute_Oijlm(1,3,1,3,vec)* b_13[local_indx_threaded].real()*std::pow(h,3) + compute_Oijlm(1,3,2,1,vec)* b_12[local_indx_threaded].real()*std::pow(h,3)
			   + compute_Oijlm(1,3,2,2,vec)* b_22[local_indx_threaded].real()*std::pow(h,3) + compute_Oijlm(1,3,2,3,vec)* b_23[local_indx_threaded].real()*std::pow(h,3)
			   + compute_Oijlm(1,3,3,1,vec)* b_13[local_indx_threaded].real()*std::pow(h,3) + compute_Oijlm(1,3,3,2,vec)* b_23[local_indx_threaded].real()*std::pow(h,3)
			   + compute_Oijlm(1,3,3,3,vec)* b_33[local_indx_threaded].real()*std::pow(h,3);

			 fab_tRe(i,j,k,Nyx::T_22) = compute_Oijlm(2,2,1,1,vec)*b_11[local_indx_threaded].real()*std::pow(h,3) + compute_Oijlm(2,2,1,2,vec)* b_12[local_indx_threaded].real()*std::pow(h,3)
			   + compute_Oijlm(2,2,1,3,vec)* b_13[local_indx_threaded].real()*std::pow(h,3) + compute_Oijlm(2,2,2,1,vec)* b_12[local_indx_threaded].real()*std::pow(h,3)
			   + compute_Oijlm(2,2,2,2,vec)* b_22[local_indx_threaded].real()*std::pow(h,3) + compute_Oijlm(2,2,2,3,vec)* b_23[local_indx_threaded].real()*std::pow(h,3)
			   + compute_Oijlm(2,2,3,1,vec)* b_13[local_indx_threaded].real()*std::pow(h,3) + compute_Oijlm(2,2,3,2,vec)* b_23[local_indx_threaded].real()*std::pow(h,3)
			   + compute_Oijlm(2,2,3,3,vec)* b_33[local_indx_threaded].real()*std::pow(h,3);

			 fab_tRe(i,j,k,Nyx::T_23) = compute_Oijlm(2,3,1,1,vec)*b_11[local_indx_threaded].real()*std::pow(h,3) + compute_Oijlm(2,3,1,2,vec)* b_12[local_indx_threaded].real()*std::pow(h,3)
			   + compute_Oijlm(2,3,1,3,vec)* b_13[local_indx_threaded].real()*std::pow(h,3) + compute_Oijlm(2,3,2,1,vec)* b_12[local_indx_threaded].real()*std::pow(h,3)
			   + compute_Oijlm(2,3,2,2,vec)* b_22[local_indx_threaded].real()*std::pow(h,3) + compute_Oijlm(2,3,2,3,vec)* b_23[local_indx_threaded].real()*std::pow(h,3)
			   + compute_Oijlm(2,3,3,1,vec)* b_13[local_indx_threaded].real()*std::pow(h,3) + compute_Oijlm(2,3,3,2,vec)* b_23[local_indx_threaded].real()*std::pow(h,3)
			   + compute_Oijlm(2,3,3,3,vec)* b_33[local_indx_threaded].real()*std::pow(h,3);

			 fab_tRe(i,j,k,Nyx::T_33) = compute_Oijlm(3,3,1,1,vec)*b_11[local_indx_threaded].real()*std::pow(h,3) + compute_Oijlm(3,3,1,2,vec)* b_12[local_indx_threaded].real()*std::pow(h,3)
			   + compute_Oijlm(3,3,1,3,vec)* b_13[local_indx_threaded].real()*std::pow(h,3) + compute_Oijlm(3,3,2,1,vec)* b_12[local_indx_threaded].real()*std::pow(h,3)
			   + compute_Oijlm(3,3,2,2,vec)* b_22[local_indx_threaded].real()*std::pow(h,3) + compute_Oijlm(3,3,2,3,vec)* b_23[local_indx_threaded].real()*std::pow(h,3)
			   + compute_Oijlm(3,3,3,1,vec)* b_13[local_indx_threaded].real()*std::pow(h,3) + compute_Oijlm(3,3,3,2,vec)* b_23[local_indx_threaded].real()*std::pow(h,3)
			   + compute_Oijlm(3,3,3,3,vec)* b_33[local_indx_threaded].real()*std::pow(h,3);


			 fab_tIm(i,j,k,Nyx::T_11) = compute_Oijlm(1,1,1,1,vec)*b_11[local_indx_threaded].imag()*std::pow(h,3) + compute_Oijlm(1,1,1,2,vec)* b_12[local_indx_threaded].imag()*std::pow(h,3)
			   + compute_Oijlm(1,1,1,3,vec)* b_13[local_indx_threaded].imag()*std::pow(h,3) + compute_Oijlm(1,1,2,1,vec)* b_12[local_indx_threaded].imag()*std::pow(h,3)
			   + compute_Oijlm(1,1,2,2,vec)* b_22[local_indx_threaded].imag()*std::pow(h,3) + compute_Oijlm(1,1,2,3,vec)* b_23[local_indx_threaded].imag()*std::pow(h,3)
			   + compute_Oijlm(1,1,3,1,vec)* b_13[local_indx_threaded].imag()*std::pow(h,3) + compute_Oijlm(1,1,3,2,vec)* b_23[local_indx_threaded].imag()*std::pow(h,3)
			   + compute_Oijlm(1,1,3,3,vec)* b_33[local_indx_threaded].imag()*std::pow(h,3);

			 fab_tIm(i,j,k,Nyx::T_12) = compute_Oijlm(1,2,1,1,vec)*b_11[local_indx_threaded].imag()*std::pow(h,3) + compute_Oijlm(1,2,1,2,vec)* b_12[local_indx_threaded].imag()*std::pow(h,3)
			   + compute_Oijlm(1,2,1,3,vec)* b_13[local_indx_threaded].imag()*std::pow(h,3) + compute_Oijlm(1,2,2,1,vec)* b_12[local_indx_threaded].imag()*std::pow(h,3)
			   + compute_Oijlm(1,2,2,2,vec)* b_22[local_indx_threaded].imag()*std::pow(h,3) + compute_Oijlm(1,2,2,3,vec)* b_23[local_indx_threaded].imag()*std::pow(h,3)
			   + compute_Oijlm(1,2,3,1,vec)* b_13[local_indx_threaded].imag()*std::pow(h,3) + compute_Oijlm(1,2,3,2,vec)* b_23[local_indx_threaded].imag()*std::pow(h,3)
			   + compute_Oijlm(1,2,3,3,vec)* b_33[local_indx_threaded].imag()*std::pow(h,3);

			 fab_tIm(i,j,k,Nyx::T_13) = compute_Oijlm(1,3,1,1,vec)*b_11[local_indx_threaded].imag()*std::pow(h,3) + compute_Oijlm(1,3,1,2,vec)* b_12[local_indx_threaded].imag()*std::pow(h,3)
			   + compute_Oijlm(1,3,1,3,vec)* b_13[local_indx_threaded].imag()*std::pow(h,3) + compute_Oijlm(1,3,2,1,vec)* b_12[local_indx_threaded].imag()*std::pow(h,3)
			   + compute_Oijlm(1,3,2,2,vec)* b_22[local_indx_threaded].imag()*std::pow(h,3) + compute_Oijlm(1,3,2,3,vec)* b_23[local_indx_threaded].imag()*std::pow(h,3)
			   + compute_Oijlm(1,3,3,1,vec)* b_13[local_indx_threaded].imag()*std::pow(h,3) + compute_Oijlm(1,3,3,2,vec)* b_23[local_indx_threaded].imag()*std::pow(h,3)
			   + compute_Oijlm(1,3,3,3,vec)* b_33[local_indx_threaded].imag()*std::pow(h,3);

			 fab_tIm(i,j,k,Nyx::T_22) = compute_Oijlm(2,2,1,1,vec)*b_11[local_indx_threaded].imag()*std::pow(h,3) + compute_Oijlm(2,2,1,2,vec)* b_12[local_indx_threaded].imag()*std::pow(h,3)
			   + compute_Oijlm(2,2,1,3,vec)* b_13[local_indx_threaded].imag()*std::pow(h,3) + compute_Oijlm(2,2,2,1,vec)* b_12[local_indx_threaded].imag()*std::pow(h,3)
			   + compute_Oijlm(2,2,2,2,vec)* b_22[local_indx_threaded].imag()*std::pow(h,3) + compute_Oijlm(2,2,2,3,vec)* b_23[local_indx_threaded].imag()*std::pow(h,3)
			   + compute_Oijlm(2,2,3,1,vec)* b_13[local_indx_threaded].imag()*std::pow(h,3) + compute_Oijlm(2,2,3,2,vec)* b_23[local_indx_threaded].imag()*std::pow(h,3)
			   + compute_Oijlm(2,2,3,3,vec)* b_33[local_indx_threaded].imag()*std::pow(h,3);

			 fab_tIm(i,j,k,Nyx::T_23) = compute_Oijlm(2,3,1,1,vec)*b_11[local_indx_threaded].imag()*std::pow(h,3) + compute_Oijlm(2,3,1,2,vec)* b_12[local_indx_threaded].imag()*std::pow(h,3)
			   + compute_Oijlm(2,3,1,3,vec)* b_13[local_indx_threaded].imag()*std::pow(h,3) + compute_Oijlm(2,3,2,1,vec)* b_12[local_indx_threaded].imag()*std::pow(h,3)
			   + compute_Oijlm(2,3,2,2,vec)* b_22[local_indx_threaded].imag()*std::pow(h,3) + compute_Oijlm(2,3,2,3,vec)* b_23[local_indx_threaded].imag()*std::pow(h,3)
			   + compute_Oijlm(2,3,3,1,vec)* b_13[local_indx_threaded].imag()*std::pow(h,3) + compute_Oijlm(2,3,3,2,vec)* b_23[local_indx_threaded].imag()*std::pow(h,3)
			   + compute_Oijlm(2,3,3,3,vec)* b_33[local_indx_threaded].imag()*std::pow(h,3);

			 fab_tIm(i,j,k,Nyx::T_33) = compute_Oijlm(3,3,1,1,vec)*b_11[local_indx_threaded].imag()*std::pow(h,3) + compute_Oijlm(3,3,1,2,vec)* b_12[local_indx_threaded].imag()*std::pow(h,3)
			   + compute_Oijlm(3,3,1,3,vec)* b_13[local_indx_threaded].imag()*std::pow(h,3) + compute_Oijlm(3,3,2,1,vec)* b_12[local_indx_threaded].imag()*std::pow(h,3)
			   + compute_Oijlm(3,3,2,2,vec)* b_22[local_indx_threaded].imag()*std::pow(h,3) + compute_Oijlm(3,3,2,3,vec)* b_23[local_indx_threaded].imag()*std::pow(h,3)
			   + compute_Oijlm(3,3,3,1,vec)* b_13[local_indx_threaded].imag()*std::pow(h,3) + compute_Oijlm(3,3,3,2,vec)* b_23[local_indx_threaded].imag()*std::pow(h,3)
			   + compute_Oijlm(3,3,3,3,vec)* b_33[local_indx_threaded].imag()*std::pow(h,3);

		       });

  } // end MFIter

   if (TRe.contains_nan(0, TRe.nComp(), 0)){ // T.nGrow()
    for (int i = 0; i < TRe.nComp(); i++)                                                                                                                                                                 
      {                                                                                                
	if (ParallelDescriptor::IOProcessor())                                                                                                                                                           
	  std::cout << "initialize tensor: testing component " << i << " for NaNs" << std::endl;                                                                                                          
	if (TRe.contains_nan(i,1,0))                                                                                                                                                                 
	  {                                                                                                                                                                                                   	    amrex::Abort("TRe has NaNs in this component");                                                                                                                                                
	  }                                                                                                                                                                                               
      }                                                                                                                                                                                                   
    amrex::Abort("TRe has NaNs before the first call.");
  }

   amrex::Print() << "Done filling tensor!" << std::endl;

}

void Nyx::compute_integrals_1(Real t_prev){

  // Another Tensor State for cos and sin integral required? Is there another way?

  // old and new states (not required if FDM_SLIM_STATE is set to TRUE in Makefile...)
  //MultiFab& T_old = get_level(level).get_old_data(Tensor_Type);
  MultiFab& TRe = get_level(level).get_old_data(TensorRe_Type);
  MultiFab& TIm = get_level(level).get_old_data(TensorIm_Type);
  MultiFab& int_cosRe = get_level(level).get_old_data(CosRe_Type);
  MultiFab& int_cosIm = get_level(level).get_old_data(CosIm_Type);
  MultiFab& int_sinRe = get_level(level).get_old_data(SinRe_Type);
  MultiFab& int_sinIm = get_level(level).get_old_data(SinIm_Type);

  Box domain(get_level(0).geom.Domain());

  Real cur_time  = state[TensorRe_Type].curTime();

  amrex::Print() << "TIMES at compute_integrals_1: " << cur_time << " " << t_prev << std::endl;
  amrex::Print() << "TIME COMPARISON: " << cur_time << " " << state[TensorRe_Type].prevTime() << std::endl;

  Real a_old = compute_a(t_prev);  //get_comoving_a(prev_time);
  Real a_new = compute_a(cur_time); //get_comoving_a(cur_time);

  Real tau_new = compute_tau(cur_time);
  Real tau_old = compute_tau(t_prev);

  amrex::Print() << "tau old/new: " << tau_old << " " << tau_new << ", diff: " << tau_new-tau_old << std::endl;
  amrex::Print() << "a_old/new in integrals1: " << a_old << " " << a_new << std::endl;

  const Real pi = 4.0 * std::atan(1.0);
  const Real tpi = 2.0 * pi;
  const Real c_u = 142948.34557923174; // speed of light in units of l_u/t_u
  const Real L = get_level(level).geom.ProbHi(0);

  for (MFIter mfi(TRe,false); mfi.isValid(); ++mfi){

    const Box& bx = mfi.validbox(); //mfi.tilebox();
    Array4<Real> const& fab_tRe = TRe[mfi].array();
    Array4<Real> const& fab_tIm = TIm[mfi].array();
    Array4<Real> const& fab_icosRe = int_cosRe[mfi].array();
    Array4<Real> const& fab_icosIm = int_cosIm[mfi].array();
    Array4<Real> const& fab_isinRe = int_sinRe[mfi].array();
    Array4<Real> const& fab_isinIm = int_sinIm[mfi].array();

    const Dim3 lo = amrex::lbound(bx);
    const Dim3 hi = amrex::ubound(bx);
    const Dim3 w ={hi.x-lo.x+1,hi.y-lo.y+1,hi.z-lo.z+1};

    double kx_, ky_, kz_, ksq_;


    AMREX_PARALLEL_FOR_3D ( bx, i, j, k, 
			    {
			      if (i < domain.length(0)/2.)
				kx_ = tpi*double(i)/L;
			      else
				kx_ = tpi*double(i - domain.length(0))/L;

			      if (j < domain.length(1)/2.)
				ky_ = tpi*double(j)/L;
			      else
				ky_ = tpi*double(j - domain.length(1))/L;

			      if (k < domain.length(2)/2.)
				kz_ = tpi*double(k)/L;
			      else
				kz_ = tpi*double(k - domain.length(2))/L;

			      ksq_ = std::sqrt(kx_*kx_ + ky_*ky_ + kz_*kz_);
			      size_t local_indx_threaded = (size_t)w.y*(size_t)w.z*(i-lo.x)+(size_t)w.z*(j-lo.y)+(k-lo.z);

			      fab_icosRe(i,j,k,Nyx::Cos_11) += 0.5*std::cos(ksq_*tau_old*c_u/tpi)*a_old*fab_tRe(i,j,k,Nyx::T_11)*(tau_new-tau_old);
			      fab_icosRe(i,j,k,Nyx::Cos_12) += 0.5*std::cos(ksq_*tau_old*c_u/tpi)*a_old*fab_tRe(i,j,k,Nyx::T_12)*(tau_new-tau_old);
			      fab_icosRe(i,j,k,Nyx::Cos_13) += 0.5*std::cos(ksq_*tau_old*c_u/tpi)*a_old*fab_tRe(i,j,k,Nyx::T_13)*(tau_new-tau_old);
			      fab_icosRe(i,j,k,Nyx::Cos_22) += 0.5*std::cos(ksq_*tau_old*c_u/tpi)*a_old*fab_tRe(i,j,k,Nyx::T_22)*(tau_new-tau_old);
			      fab_icosRe(i,j,k,Nyx::Cos_23) += 0.5*std::cos(ksq_*tau_old*c_u/tpi)*a_old*fab_tRe(i,j,k,Nyx::T_23)*(tau_new-tau_old);
			      fab_icosRe(i,j,k,Nyx::Cos_33) += 0.5*std::cos(ksq_*tau_old*c_u/tpi)*a_old*fab_tRe(i,j,k,Nyx::T_33)*(tau_new-tau_old);

			      fab_isinRe(i,j,k,Nyx::Sin_11) += 0.5*std::sin(ksq_*tau_old*c_u/tpi)*a_old*fab_tRe(i,j,k,Nyx::T_11)*(tau_new-tau_old);
			      fab_isinRe(i,j,k,Nyx::Sin_12) += 0.5*std::sin(ksq_*tau_old*c_u/tpi)*a_old*fab_tRe(i,j,k,Nyx::T_12)*(tau_new-tau_old);
			      fab_isinRe(i,j,k,Nyx::Sin_13) += 0.5*std::sin(ksq_*tau_old*c_u/tpi)*a_old*fab_tRe(i,j,k,Nyx::T_13)*(tau_new-tau_old);
			      fab_isinRe(i,j,k,Nyx::Sin_22) += 0.5*std::sin(ksq_*tau_old*c_u/tpi)*a_old*fab_tRe(i,j,k,Nyx::T_22)*(tau_new-tau_old);
			      fab_isinRe(i,j,k,Nyx::Sin_23) += 0.5*std::sin(ksq_*tau_old*c_u/tpi)*a_old*fab_tRe(i,j,k,Nyx::T_23)*(tau_new-tau_old);
			      fab_isinRe(i,j,k,Nyx::Sin_33) += 0.5*std::sin(ksq_*tau_old*c_u/tpi)*a_old*fab_tRe(i,j,k,Nyx::T_33)*(tau_new-tau_old);

			      fab_icosIm(i,j,k,Nyx::Cos_11) += 0.5*std::cos(ksq_*tau_old*c_u/tpi)*a_old*fab_tIm(i,j,k,Nyx::T_11)*(tau_new-tau_old);
			      fab_icosIm(i,j,k,Nyx::Cos_12) += 0.5*std::cos(ksq_*tau_old*c_u/tpi)*a_old*fab_tIm(i,j,k,Nyx::T_12)*(tau_new-tau_old);
			      fab_icosIm(i,j,k,Nyx::Cos_13) += 0.5*std::cos(ksq_*tau_old*c_u/tpi)*a_old*fab_tIm(i,j,k,Nyx::T_13)*(tau_new-tau_old);
			      fab_icosIm(i,j,k,Nyx::Cos_22) += 0.5*std::cos(ksq_*tau_old*c_u/tpi)*a_old*fab_tIm(i,j,k,Nyx::T_22)*(tau_new-tau_old);
			      fab_icosIm(i,j,k,Nyx::Cos_23) += 0.5*std::cos(ksq_*tau_old*c_u/tpi)*a_old*fab_tIm(i,j,k,Nyx::T_23)*(tau_new-tau_old);
			      fab_icosIm(i,j,k,Nyx::Cos_33) += 0.5*std::cos(ksq_*tau_old*c_u/tpi)*a_old*fab_tIm(i,j,k,Nyx::T_33)*(tau_new-tau_old);

			      fab_isinIm(i,j,k,Nyx::Sin_11) += 0.5*std::sin(ksq_*tau_old*c_u/tpi)*a_old*fab_tIm(i,j,k,Nyx::T_11)*(tau_new-tau_old);
			      fab_isinIm(i,j,k,Nyx::Sin_12) += 0.5*std::sin(ksq_*tau_old*c_u/tpi)*a_old*fab_tIm(i,j,k,Nyx::T_12)*(tau_new-tau_old);
			      fab_isinIm(i,j,k,Nyx::Sin_13) += 0.5*std::sin(ksq_*tau_old*c_u/tpi)*a_old*fab_tIm(i,j,k,Nyx::T_13)*(tau_new-tau_old);
			      fab_isinIm(i,j,k,Nyx::Sin_22) += 0.5*std::sin(ksq_*tau_old*c_u/tpi)*a_old*fab_tIm(i,j,k,Nyx::T_22)*(tau_new-tau_old);
			      fab_isinIm(i,j,k,Nyx::Sin_23) += 0.5*std::sin(ksq_*tau_old*c_u/tpi)*a_old*fab_tIm(i,j,k,Nyx::T_23)*(tau_new-tau_old);
			      fab_isinIm(i,j,k,Nyx::Sin_33) += 0.5*std::sin(ksq_*tau_old*c_u/tpi)*a_old*fab_tIm(i,j,k,Nyx::T_33)*(tau_new-tau_old);

			      // use this as approximation without imaginary parts in cos/sin
			      // fab_icosRe(i,j,k,Nyx::Cos_11) += 0.5*std::cos(ksq_*tau_old*c_u/tpi)*a_old*std::sqrt(std::pow(fab_tRe(i,j,k,Nyx::T_11),2)+std::pow(fab_tIm(i,j,k,Nyx::T_11),2))*(tau_new-tau_old);
			      // fab_icosRe(i,j,k,Nyx::Cos_12) += 0.5*std::cos(ksq_*tau_old*c_u/tpi)*a_old*std::sqrt(std::pow(fab_tRe(i,j,k,Nyx::T_12),2)+std::pow(fab_tIm(i,j,k,Nyx::T_12),2))*(tau_new-tau_old);
			      // fab_icosRe(i,j,k,Nyx::Cos_13) += 0.5*std::cos(ksq_*tau_old*c_u/tpi)*a_old*std::sqrt(std::pow(fab_tRe(i,j,k,Nyx::T_13),2)+std::pow(fab_tIm(i,j,k,Nyx::T_13),2))*(tau_new-tau_old);
			      // fab_icosRe(i,j,k,Nyx::Cos_22) += 0.5*std::cos(ksq_*tau_old*c_u/tpi)*a_old*std::sqrt(std::pow(fab_tRe(i,j,k,Nyx::T_22),2)+std::pow(fab_tIm(i,j,k,Nyx::T_22),2))*(tau_new-tau_old);
			      // fab_icosRe(i,j,k,Nyx::Cos_23) += 0.5*std::cos(ksq_*tau_old*c_u/tpi)*a_old*std::sqrt(std::pow(fab_tRe(i,j,k,Nyx::T_23),2)+std::pow(fab_tIm(i,j,k,Nyx::T_23),2))*(tau_new-tau_old);
			      // fab_icosRe(i,j,k,Nyx::Cos_33) += 0.5*std::cos(ksq_*tau_old*c_u/tpi)*a_old*std::sqrt(std::pow(fab_tRe(i,j,k,Nyx::T_33),2)+std::pow(fab_tIm(i,j,k,Nyx::T_33),2))*(tau_new-tau_old);
			      
			      // fab_isinRe(i,j,k,Nyx::Sin_11) += 0.5*std::sin(ksq_*tau_old*c_u/tpi)*a_old*std::sqrt(std::pow(fab_tRe(i,j,k,Nyx::T_11),2)+std::pow(fab_tIm(i,j,k,Nyx::T_11),2))*(tau_new-tau_old);
			      // fab_isinRe(i,j,k,Nyx::Sin_12) += 0.5*std::sin(ksq_*tau_old*c_u/tpi)*a_old*std::sqrt(std::pow(fab_tRe(i,j,k,Nyx::T_12),2)+std::pow(fab_tIm(i,j,k,Nyx::T_12),2))*(tau_new-tau_old);
			      // fab_isinRe(i,j,k,Nyx::Sin_13) += 0.5*std::sin(ksq_*tau_old*c_u/tpi)*a_old*std::sqrt(std::pow(fab_tRe(i,j,k,Nyx::T_13),2)+std::pow(fab_tIm(i,j,k,Nyx::T_13),2))*(tau_new-tau_old);
			      // fab_isinRe(i,j,k,Nyx::Sin_22) += 0.5*std::sin(ksq_*tau_old*c_u/tpi)*a_old*std::sqrt(std::pow(fab_tRe(i,j,k,Nyx::T_22),2)+std::pow(fab_tIm(i,j,k,Nyx::T_22),2))*(tau_new-tau_old);
			      // fab_isinRe(i,j,k,Nyx::Sin_23) += 0.5*std::sin(ksq_*tau_old*c_u/tpi)*a_old*std::sqrt(std::pow(fab_tRe(i,j,k,Nyx::T_23),2)+std::pow(fab_tIm(i,j,k,Nyx::T_23),2))*(tau_new-tau_old);
			      // fab_isinRe(i,j,k,Nyx::Sin_33) += 0.5*std::sin(ksq_*tau_old*c_u/tpi)*a_old*std::sqrt(std::pow(fab_tRe(i,j,k,Nyx::T_33),2)+std::pow(fab_tIm(i,j,k,Nyx::T_33),2))*(tau_new-tau_old);
			    });

  }

}

void Nyx::compute_integrals_2(Real t_prev){

  MultiFab& TRe = get_level(level).get_new_data(TensorRe_Type);
  MultiFab& TIm = get_level(level).get_new_data(TensorIm_Type);
  MultiFab& int_cosRe = get_level(level).get_new_data(CosRe_Type);
  MultiFab& int_cosIm = get_level(level).get_new_data(CosIm_Type);
  MultiFab& int_sinRe = get_level(level).get_new_data(SinRe_Type);
  MultiFab& int_sinIm = get_level(level).get_new_data(SinIm_Type);

  Box domain(get_level(0).geom.Domain());

  Real cur_time  = state[TensorRe_Type].curTime();

  amrex::Print() << "TIMES at compute_integrals_2: " << cur_time << " " << t_prev << std::endl;

  Real a_old = compute_a(t_prev); //get_comoving_a(prev_time);
  Real a_new = compute_a(cur_time); //get_comoving_a(cur_time);

  Real tau_new = compute_tau(cur_time);
  Real tau_old = compute_tau(t_prev);

  const Real pi = 4.0 * std::atan(1.0);
  const Real tpi = 2.0 * pi;
  const Real c_u = 142948.34557923174; // speed of light in units of l_u/t_u
  const Real L = get_level(level).geom.ProbHi(0);

  for (MFIter mfi(TRe,false); mfi.isValid(); ++mfi){
    const Box& bx = mfi.validbox(); //mfi.tilebox();
    Array4<Real> const& fab_tRe = TRe[mfi].array();
    Array4<Real> const& fab_tIm = TIm[mfi].array();
    Array4<Real> const& fab_icosRe = int_cosRe[mfi].array();
    Array4<Real> const& fab_icosIm = int_cosIm[mfi].array();
    Array4<Real> const& fab_isinRe = int_sinRe[mfi].array();
    Array4<Real> const& fab_isinIm = int_sinIm[mfi].array();

    const Dim3 lo = amrex::lbound(bx);
    const Dim3 hi = amrex::ubound(bx);
    const Dim3 w ={hi.x-lo.x+1,hi.y-lo.y+1,hi.z-lo.z+1};

    double kx_, ky_, kz_, ksq_;


    AMREX_PARALLEL_FOR_3D ( bx, i, j, k, 
			    {
			      if (i < domain.length(0)/2.)
				kx_ = tpi*double(i)/L;
			      else
				kx_ = tpi*double(i - domain.length(0))/L;

			      if (j < domain.length(1)/2.)
				ky_ = tpi*double(j)/L;
			      else
				ky_ = tpi*double(j - domain.length(1))/L;

			      if (k < domain.length(2)/2.)
				kz_ = tpi*double(k)/L;
			      else
				kz_ = tpi*double(k - domain.length(2))/L;

			      ksq_ = std::sqrt(kx_*kx_ + ky_*ky_ + kz_*kz_);

			      size_t local_indx_threaded = (size_t)w.y*(size_t)w.z*(i-lo.x)+(size_t)w.z*(j-lo.y)+(k-lo.z);
			      
			      fab_icosRe(i,j,k,Nyx::Cos_11) += 0.5*std::cos(ksq_*tau_new*c_u/tpi)*a_new*fab_tRe(i,j,k,Nyx::T_11)*(tau_new-tau_old);
			      fab_icosRe(i,j,k,Nyx::Cos_12) += 0.5*std::cos(ksq_*tau_new*c_u/tpi)*a_new*fab_tRe(i,j,k,Nyx::T_12)*(tau_new-tau_old);
			      fab_icosRe(i,j,k,Nyx::Cos_13) += 0.5*std::cos(ksq_*tau_new*c_u/tpi)*a_new*fab_tRe(i,j,k,Nyx::T_13)*(tau_new-tau_old);
			      fab_icosRe(i,j,k,Nyx::Cos_22) += 0.5*std::cos(ksq_*tau_new*c_u/tpi)*a_new*fab_tRe(i,j,k,Nyx::T_22)*(tau_new-tau_old);
			      fab_icosRe(i,j,k,Nyx::Cos_23) += 0.5*std::cos(ksq_*tau_new*c_u/tpi)*a_new*fab_tRe(i,j,k,Nyx::T_23)*(tau_new-tau_old);
			      fab_icosRe(i,j,k,Nyx::Cos_33) += 0.5*std::cos(ksq_*tau_new*c_u/tpi)*a_new*fab_tRe(i,j,k,Nyx::T_33)*(tau_new-tau_old);

			      fab_isinRe(i,j,k,Nyx::Sin_11) += 0.5*std::sin(ksq_*tau_new*c_u/tpi)*a_new*fab_tRe(i,j,k,Nyx::T_11)*(tau_new-tau_old);
			      fab_isinRe(i,j,k,Nyx::Sin_12) += 0.5*std::sin(ksq_*tau_new*c_u/tpi)*a_new*fab_tRe(i,j,k,Nyx::T_12)*(tau_new-tau_old);
			      fab_isinRe(i,j,k,Nyx::Sin_13) += 0.5*std::sin(ksq_*tau_new*c_u/tpi)*a_new*fab_tRe(i,j,k,Nyx::T_13)*(tau_new-tau_old);
			      fab_isinRe(i,j,k,Nyx::Sin_22) += 0.5*std::sin(ksq_*tau_new*c_u/tpi)*a_new*fab_tRe(i,j,k,Nyx::T_22)*(tau_new-tau_old);
			      fab_isinRe(i,j,k,Nyx::Sin_23) += 0.5*std::sin(ksq_*tau_new*c_u/tpi)*a_new*fab_tRe(i,j,k,Nyx::T_23)*(tau_new-tau_old);
			      fab_isinRe(i,j,k,Nyx::Sin_33) += 0.5*std::sin(ksq_*tau_new*c_u/tpi)*a_new*fab_tRe(i,j,k,Nyx::T_33)*(tau_new-tau_old);

			      fab_icosIm(i,j,k,Nyx::Cos_11) += 0.5*std::cos(ksq_*tau_new*c_u/tpi)*a_new*fab_tIm(i,j,k,Nyx::T_11)*(tau_new-tau_old);
			      fab_icosIm(i,j,k,Nyx::Cos_12) += 0.5*std::cos(ksq_*tau_new*c_u/tpi)*a_new*fab_tIm(i,j,k,Nyx::T_12)*(tau_new-tau_old);
			      fab_icosIm(i,j,k,Nyx::Cos_13) += 0.5*std::cos(ksq_*tau_new*c_u/tpi)*a_new*fab_tIm(i,j,k,Nyx::T_13)*(tau_new-tau_old);
			      fab_icosIm(i,j,k,Nyx::Cos_22) += 0.5*std::cos(ksq_*tau_new*c_u/tpi)*a_new*fab_tIm(i,j,k,Nyx::T_22)*(tau_new-tau_old);
			      fab_icosIm(i,j,k,Nyx::Cos_23) += 0.5*std::cos(ksq_*tau_new*c_u/tpi)*a_new*fab_tIm(i,j,k,Nyx::T_23)*(tau_new-tau_old);
			      fab_icosIm(i,j,k,Nyx::Cos_33) += 0.5*std::cos(ksq_*tau_new*c_u/tpi)*a_new*fab_tIm(i,j,k,Nyx::T_33)*(tau_new-tau_old);

			      fab_isinIm(i,j,k,Nyx::Sin_11) += 0.5*std::sin(ksq_*tau_new*c_u/tpi)*a_new*fab_tIm(i,j,k,Nyx::T_11)*(tau_new-tau_old);
			      fab_isinIm(i,j,k,Nyx::Sin_12) += 0.5*std::sin(ksq_*tau_new*c_u/tpi)*a_new*fab_tIm(i,j,k,Nyx::T_12)*(tau_new-tau_old);
			      fab_isinIm(i,j,k,Nyx::Sin_13) += 0.5*std::sin(ksq_*tau_new*c_u/tpi)*a_new*fab_tIm(i,j,k,Nyx::T_13)*(tau_new-tau_old);
			      fab_isinIm(i,j,k,Nyx::Sin_22) += 0.5*std::sin(ksq_*tau_new*c_u/tpi)*a_new*fab_tIm(i,j,k,Nyx::T_22)*(tau_new-tau_old);
			      fab_isinIm(i,j,k,Nyx::Sin_23) += 0.5*std::sin(ksq_*tau_new*c_u/tpi)*a_new*fab_tIm(i,j,k,Nyx::T_23)*(tau_new-tau_old);
			      fab_isinIm(i,j,k,Nyx::Sin_33) += 0.5*std::sin(ksq_*tau_new*c_u/tpi)*a_new*fab_tIm(i,j,k,Nyx::T_33)*(tau_new-tau_old);

			      // use this as approximation without imaginary parts in cos/sin
			      // fab_icosRe(i,j,k,Nyx::Cos_11) += 0.5*std::cos(ksq_*tau_new*c_u/tpi)*a_new*std::sqrt(std::pow(fab_tRe(i,j,k,Nyx::T_11),2)+std::pow(fab_tIm(i,j,k,Nyx::T_11),2))*(tau_new-tau_old);
			      // fab_icosRe(i,j,k,Nyx::Cos_12) += 0.5*std::cos(ksq_*tau_new*c_u/tpi)*a_new*std::sqrt(std::pow(fab_tRe(i,j,k,Nyx::T_12),2)+std::pow(fab_tIm(i,j,k,Nyx::T_12),2))*(tau_new-tau_old);
			      // fab_icosRe(i,j,k,Nyx::Cos_13) += 0.5*std::cos(ksq_*tau_new*c_u/tpi)*a_new*std::sqrt(std::pow(fab_tRe(i,j,k,Nyx::T_13),2)+std::pow(fab_tIm(i,j,k,Nyx::T_13),2))*(tau_new-tau_old);
			      // fab_icosRe(i,j,k,Nyx::Cos_22) += 0.5*std::cos(ksq_*tau_new*c_u/tpi)*a_new*std::sqrt(std::pow(fab_tRe(i,j,k,Nyx::T_22),2)+std::pow(fab_tIm(i,j,k,Nyx::T_22),2))*(tau_new-tau_old);
			      // fab_icosRe(i,j,k,Nyx::Cos_23) += 0.5*std::cos(ksq_*tau_new*c_u/tpi)*a_new*std::sqrt(std::pow(fab_tRe(i,j,k,Nyx::T_23),2)+std::pow(fab_tIm(i,j,k,Nyx::T_23),2))*(tau_new-tau_old);
			      // fab_icosRe(i,j,k,Nyx::Cos_33) += 0.5*std::cos(ksq_*tau_new*c_u/tpi)*a_new*std::sqrt(std::pow(fab_tRe(i,j,k,Nyx::T_33),2)+std::pow(fab_tIm(i,j,k,Nyx::T_33),2))*(tau_new-tau_old);
			      
			      // fab_isinRe(i,j,k,Nyx::Sin_11) += 0.5*std::sin(ksq_*tau_new*c_u/tpi)*a_new*std::sqrt(std::pow(fab_tRe(i,j,k,Nyx::T_11),2)+std::pow(fab_tIm(i,j,k,Nyx::T_11),2))*(tau_new-tau_old);
			      // fab_isinRe(i,j,k,Nyx::Sin_12) += 0.5*std::sin(ksq_*tau_new*c_u/tpi)*a_new*std::sqrt(std::pow(fab_tRe(i,j,k,Nyx::T_12),2)+std::pow(fab_tIm(i,j,k,Nyx::T_12),2))*(tau_new-tau_old);
			      // fab_isinRe(i,j,k,Nyx::Sin_13) += 0.5*std::sin(ksq_*tau_new*c_u/tpi)*a_new*std::sqrt(std::pow(fab_tRe(i,j,k,Nyx::T_13),2)+std::pow(fab_tIm(i,j,k,Nyx::T_13),2))*(tau_new-tau_old);
			      // fab_isinRe(i,j,k,Nyx::Sin_22) += 0.5*std::sin(ksq_*tau_new*c_u/tpi)*a_new*std::sqrt(std::pow(fab_tRe(i,j,k,Nyx::T_22),2)+std::pow(fab_tIm(i,j,k,Nyx::T_22),2))*(tau_new-tau_old);
			      // fab_isinRe(i,j,k,Nyx::Sin_23) += 0.5*std::sin(ksq_*tau_new*c_u/tpi)*a_new*std::sqrt(std::pow(fab_tRe(i,j,k,Nyx::T_23),2)+std::pow(fab_tIm(i,j,k,Nyx::T_23),2))*(tau_new-tau_old);
			      // fab_isinRe(i,j,k,Nyx::Sin_33) += 0.5*std::sin(ksq_*tau_new*c_u/tpi)*a_new*std::sqrt(std::pow(fab_tRe(i,j,k,Nyx::T_33),2)+std::pow(fab_tIm(i,j,k,Nyx::T_33),2))*(tau_new-tau_old);
			    });

  }

}


void Nyx::compute_Sk(){
  
  amrex::Print() << "Start compute_Sk" << std::endl;
  const Real c_u = 142948.34557923174; // speed of light in units of l_u/t_u
  const Real pi = 4.0 * std::atan(1.0);
  const Real tpi = 2.0 * pi;

  //MultiFab& Ax_new = get_level(level).get_new_data(Axion_Type); 
  MultiFab& int_cosRe = get_level(level).get_new_data(CosRe_Type);
  MultiFab& int_cosIm = get_level(level).get_new_data(CosIm_Type);
  MultiFab& int_sinRe = get_level(level).get_new_data(SinRe_Type);
  MultiFab& int_sinIm = get_level(level).get_new_data(SinIm_Type);
  MultiFab S_k(int_cosRe.boxArray(), int_cosRe.DistributionMap(), 1, 0);
  S_k.setVal(0.);

  amrex::Print() << "Defined Multifabs in compute_Sk" << std::endl;

  const BoxArray& ba = int_cosRe.boxArray();
  Box domain(get_level(0).geom.Domain());
  const Real L = get_level(level).geom.ProbHi(0); 
  const int Nx = ba[0].size()[0];
  //double kx, ky, kz, k2;

  const int nbins = std::ceil(std::sqrt(3)*std::floor(domain.length(0)/2));
  //  amrex::Print() << "NUMBER OF BINS FOR SPECTRUM: " << nbins << " " << std::ceil(std::floor(domain.length(0)/2)) <<std::endl;
  Vector<Real> spectrum(nbins,0.0);
  Vector<Real> mode(nbins,0.0);
  Vector<int>  nmode(nbins,0);

  for (MFIter mfi(int_cosRe,false); mfi.isValid(); ++mfi){
    const Box& bx = mfi.validbox(); //mfi.tilebox();
    Array4<Real> const& fab_icosRe = int_cosRe[mfi].array();
    Array4<Real> const& fab_icosIm = int_cosIm[mfi].array();
    Array4<Real> const& fab_isinRe = int_sinRe[mfi].array();
    Array4<Real> const& fab_isinIm = int_sinIm[mfi].array();
    Array4<Real> const& fab_Sk = S_k[mfi].array();

    const Dim3 lo = amrex::lbound(bx);
    const Dim3 hi = amrex::ubound(bx);
    const Dim3 w ={hi.x-lo.x+1,hi.y-lo.y+1,hi.z-lo.z+1};
    
    double kx_, ky_, kz_, ksq_;

    AMREX_PARALLEL_FOR_3D ( bx, i, j, k, 
			    {
			      if (i < domain.length(0)/2.)
				kx_ = tpi*double(i)/L;
			      else
				kx_ = tpi*double(i - domain.length(0))/L;

			      if (j < domain.length(1)/2.)
				ky_ = tpi*double(j)/L;
			      else
				ky_ = tpi*double(j - domain.length(1))/L;

			      if (k < domain.length(2)/2.)
				kz_ = tpi*double(k)/L;
			      else
				kz_ = tpi*double(k - domain.length(2))/L;

			      ksq_ = std::sqrt(kx_*kx_ + ky_*ky_ + kz_*kz_);
			      
			      size_t local_indx_threaded = (size_t)w.y*(size_t)w.z*(i-lo.x)+(size_t)w.z*(j-lo.y)+(k-lo.z);

			      fab_Sk(i,j,k) = 4.0*pi*Gconst/(L*L*L*std::pow(c_u,4))*std::pow(ksq_,3)*  //std::pow(std::sqrt(k2[local_indx_threaded]),3)*
			      	( std::pow(fab_icosRe(i,j,k,Cos_11),2) + std::pow(fab_icosIm(i,j,k,Cos_11),2)
			      	  + std::pow(fab_icosRe(i,j,k,Cos_22),2) + std::pow(fab_icosIm(i,j,k,Cos_22),2)
			      	  + std::pow(fab_icosRe(i,j,k,Cos_33),2) + std::pow(fab_icosIm(i,j,k,Cos_33),2)
			      	  + 2.0*std::pow(fab_icosRe(i,j,k,Cos_12),2) + 2.0*std::pow(fab_icosIm(i,j,k,Cos_12),2)
			      	  + 2.0*std::pow(fab_icosRe(i,j,k,Cos_13),2) + 2.0*std::pow(fab_icosIm(i,j,k,Cos_13),2)
			      	  + 2.0*std::pow(fab_icosRe(i,j,k,Cos_23),2) + 2.0*std::pow(fab_icosIm(i,j,k,Cos_23),2)
			      	  + std::pow(fab_isinRe(i,j,k,Sin_11),2) + std::pow(fab_isinIm(i,j,k,Sin_11),2)
			      	  + std::pow(fab_isinRe(i,j,k,Sin_22),2) + std::pow(fab_isinIm(i,j,k,Sin_22),2)
			      	  + std::pow(fab_isinRe(i,j,k,Sin_33),2) + std::pow(fab_isinIm(i,j,k,Sin_33),2)
			      	  + 2.0*std::pow(fab_isinRe(i,j,k,Sin_12),2) + 2.0*std::pow(fab_isinIm(i,j,k,Sin_12),2)
			      	  + 2.0*std::pow(fab_isinRe(i,j,k,Sin_13),2) + 2.0*std::pow(fab_isinIm(i,j,k,Sin_13),2)
			      	  + 2.0*std::pow(fab_isinRe(i,j,k,Sin_23),2) + 2.0*std::pow(fab_isinIm(i,j,k,Sin_23),2)
			      	  );

			      // use this as approximation without imaginary cos/sin integral parts
			      // fab_Sk(i,j,k) = 4.0*pi*Gconst/(L*L*L*std::pow(c_u,4))*std::pow(ksq_,3)*			       
			      // 	( pow(std::abs(fab_icosRe(i,j,k,Nyx::Cos_11)),2) + pow(std::abs(fab_isinRe(i,j,k,Nyx::Sin_11)),2)
			      // 	  + pow(std::abs(fab_icosRe(i,j,k,Nyx::Cos_22)),2) + pow(std::abs(fab_isinRe(i,j,k,Nyx::Sin_22)),2)
			      // 	  + pow(std::abs(fab_icosRe(i,j,k,Nyx::Cos_33)),2) + pow(std::abs(fab_isinRe(i,j,k,Nyx::Sin_33)),2)
			      // 	  + 2.0*pow(std::abs(fab_icosRe(i,j,k,Nyx::Cos_12)),2) + 2.0*pow(std::abs(fab_isinRe(i,j,k,Nyx::Sin_12)),2)
			      // 	  + 2.0*pow(std::abs(fab_icosRe(i,j,k,Nyx::Cos_13)),2) + 2.0*pow(std::abs(fab_isinRe(i,j,k,Nyx::Sin_13)),2)
			      // 	  + 2.0*pow(std::abs(fab_icosRe(i,j,k,Nyx::Cos_23)),2) + 2.0*pow(std::abs(fab_isinRe(i,j,k,Nyx::Sin_23)),2)
			      // 	  );

			      const int bin = std::floor(std::sqrt((i-domain.length(0)/2)*(i-domain.length(0)/2) + (j-domain.length(1)/2)*(j-domain.length(1)/2) 
			      					   + (k-domain.length(2)/2)*(k-domain.length(2)/2)));
			      
			      nmode[bin]++;
			      spectrum[bin] += fab_Sk(i,j,k);
			      mode[bin] += ksq_; //std::sqrt(k2[local_indx_threaded]);
			    });

  }

    // output 1D spectrum

    //if (parent->levelSteps(level)%Nyx::gw_nstep_spectrum == 0){ // introduce gw_nstep_spectrum in Nyx.cpp
  
      int nstep = parent->levelSteps(level);

      ParallelDescriptor::ReduceIntSum(nmode.dataPtr(), nmode.size(),
    				     ParallelDescriptor::IOProcessorNumber());
      ParallelDescriptor::ReduceRealSum(spectrum.dataPtr(), spectrum.size(),
    				      ParallelDescriptor::IOProcessorNumber());
      ParallelDescriptor::ReduceRealSum(mode.dataPtr(), mode.size(),
    				      ParallelDescriptor::IOProcessorNumber());

      if (ParallelDescriptor::IOProcessor()) {
	if ( ! amrex::UtilCreateDirectory("./spectra", 0755)) {
	  amrex::CreateDirectoryFailed("./spectra");
	}
	std::string fileName(gw_log+"/Sk_"+std::to_string(1+nstep));
	std::stringstream stream;
	//stream << std::fixed << std::setprecision(5) << time;
	//fileName += stream.str();
	std::ofstream spectrum_out(fileName.c_str());
	
	for(int bin=0; bin < spectrum.size(); ++bin) {
	  spectrum_out << mode[bin] << " " << spectrum[bin] << " " << nmode[bin] << std::endl;
	}
	spectrum_out.close();
      }

      //}
 
  amrex::Print() << "Finished computing S_k" << std::endl;
}

void Nyx::compute_grav(){
  
  compute_integrals_1(gw_prevtime);
  amrex::Print() << "Done with computing integrals_1" << std::endl;
  fill_tensor_TT();
  compute_integrals_2(gw_prevtime);
  amrex::Print() << "Done with computing integrals_2" << std::endl;
  compute_Sk();
  amrex::Print() << "Done with computing S_k" << std::endl;
  gw_prevtime = state[State_for_Time].curTime();
}

#endif
#endif // INFGRAV
