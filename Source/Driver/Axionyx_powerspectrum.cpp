#ifdef AXCOMPLEX
#ifdef BL_USE_MPI
#include <AMReX_BLProfiler.H>
#include <Nyx.H>
#include <AMReX_MultiFab.H>
#ifdef GRAVITY
#include <Gravity.H>
#endif
#include <Distribution.H>
#include <AlignedAllocator.h>
#include <Dfft.H>

#include <string>

#define ALIGN 16

using namespace amrex;

void Nyx::prepare_and_compute_powerspectrum(amrex::Real time)
{
  
  BL_PROFILE("Nyx::compute_FDM_powerspectrum()");

  MultiFab mf = construct_regular_mf(time);
  hacc::Distribution d = initialize_fft_distribution(mf); 
  hacc::Dfft dfft(d);
  std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > a = initialize_fft_input_field(mf);
  std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > b;
  b.resize(a.size());
  dfft.makePlans(&a[0],&b[0],&a[0],&b[0]);
  dfft.forward(&(a[0]));    
  actually_compute_powerspectrum(dfft,a,time);
}    

MultiFab Nyx::construct_regular_mf(amrex::Real time){
  BoxArray ba(Domain());
  int nprocs = ParallelDescriptor::NProcs();
  int useful_nprocs = 1;
  int grid_length = Domain().length(0);
  while(nprocs>=8){
    useful_nprocs *=8;
    nprocs /=8;
    grid_length /= 2;
  }
  if(useful_nprocs<nprocs)
    amrex::Print()<<"WARNING: Since fft needs regularly spaced grids, only "<<useful_nprocs<<" are used in fft!"<<std::endl;
  ba.maxSize(grid_length);
  Vector<int> pmap(useful_nprocs);
  for(int i =0; i<pmap.size(); i++) pmap[i]=i;
  DistributionMapping dm(pmap);
  MultiFab fft(ba, dm, 1, 0);
#ifdef AXCOMPLEX
  const auto& mf = get_level(0).derive("axion_velocity", time, 0);
  fft.ParallelCopy(*mf, 0, 0, fft.nComp(), (*mf).nGrow(), fft.nGrow(),parent->Geom(level).periodicity(), FabArrayBase::COPY);
#else
  MultiFab& mf = get_level(0).get_new_data(Axion_Type);
  fft.ParallelCopy(mf, 0, 0, fft.nComp(), mf.nGrow(), fft.nGrow(),parent->Geom(level).periodicity(), FabArrayBase::COPY);
#endif
  return fft;
}


const hacc::Distribution Nyx::initialize_fft_distribution(MultiFab& mf){
  
  const BoxArray& ba = mf.boxArray();
  const DistributionMapping& dm = mf.DistributionMap();
  
  Vector<int> rank_mapping;
  rank_mapping.resize(ba.size());
  
  Vector<int> smallxEnds;
  for (int ib = 0; ib < ba.size(); ++ib){
    smallxEnds.push_back(ba[ib].smallEnd(0));
  }
  std::sort(smallxEnds.begin(), smallxEnds.end());
  smallxEnds.erase( unique( smallxEnds.begin(), smallxEnds.end() ), smallxEnds.end() );
  
  Vector<int> smallyEnds;
  for (int ib = 0; ib < ba.size(); ++ib){
    smallyEnds.push_back(ba[ib].smallEnd(1));
  }
  std::sort(smallyEnds.begin(), smallyEnds.end());
  smallyEnds.erase( unique( smallyEnds.begin(), smallyEnds.end() ), smallyEnds.end() );
  
  Vector<int> smallzEnds;
  for (int ib = 0; ib < ba.size(); ++ib){
    smallzEnds.push_back(ba[ib].smallEnd(2));
  }
  std::sort(smallzEnds.begin(), smallzEnds.end());
  smallzEnds.erase( unique( smallzEnds.begin(), smallzEnds.end() ), smallzEnds.end() );
  
  if ((smallxEnds.size()*smallyEnds.size()*smallzEnds.size()) != ba.size())
    amrex::Error("GRIDS NEED TO BE DISTRIBUTED MORE REGULARLY");
  if(ParallelDescriptor::NProcs() != ba.size())
    amrex::Error("Number of MPI ranks has to be equal to the number of root level grids!");
  
  for (int ib = 0; ib < ba.size(); ++ib)
    {
      int i = (std::find(smallxEnds.begin(), smallxEnds.end(), ba[ib].smallEnd(0))-smallxEnds.begin());
      int j = (std::find(smallyEnds.begin(), smallyEnds.end(), ba[ib].smallEnd(1))-smallyEnds.begin());
      int k = (std::find(smallzEnds.begin(), smallzEnds.end(), ba[ib].smallEnd(2))-smallzEnds.begin());
      int local_index = i*smallzEnds.size()*smallyEnds.size() + j*smallzEnds.size() + k;  
      rank_mapping[local_index] = dm[ib];
    }
  int Ndims[3] = { smallzEnds.size(), smallyEnds.size(), smallxEnds.size() };
  int     n[3] = {Domain().length(2), Domain().length(1), Domain().length(0)};
  
  return hacc::Distribution(MPI_COMM_WORLD,n,Ndims,&rank_mapping[0]);
}


 std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > Nyx::initialize_fft_input_field(amrex::MultiFab &mf){
  
  std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > a;

  for (MFIter mfi(mf,false); mfi.isValid(); ++mfi){
    Array4<Real> const& arr = mf.array(mfi);
    const Box& bx = mfi.validbox();
    const Dim3 lo = amrex::lbound(bx);
    const Dim3 hi = amrex::ubound(bx);
    const Dim3 w ={hi.x-lo.x+1,hi.y-lo.y+1,hi.z-lo.z+1};
    const int gridsize = (hi.x-lo.x+1)*(hi.y-lo.y+1)*(hi.z-lo.z+1);
    a.resize(gridsize);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(size_t k=0; k<(size_t)w.z; k++) {
      for(size_t j=0; j<(size_t)w.y; j++) {
	AMREX_PRAGMA_SIMD
	  for(size_t i=0; i<(size_t)w.x; i++) {

	    size_t local_indx_threaded = (size_t)w.y*(size_t)w.z*i+(size_t)w.z*j+k;	    
	    complex_t temp(arr(i+lo.x,j+lo.y,k+lo.z,0),0);
	    a[local_indx_threaded] = temp;		
	  }}}
  }
  return a;
}


void Nyx::actually_compute_powerspectrum(hacc::Dfft& dfft, std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> >& a, amrex::Real time){
  const Real pi = 4 * std::atan(1.0);
  const Real tpi = 2 * pi;
  const Real h = get_level(0).geom.CellSize(0);
  const int *self = dfft.self_kspace();
  const int *local_ng = dfft.local_ng_kspace();
  const int *global_ng = dfft.global_ng();
  const int nbins = std::ceil(std::sqrt(3)*std::floor(Domain().length(0)/2));
  Vector<Real> spectrum(nbins,0.0);
  Vector<Real> mode(nbins,0.0);
  Vector<int>  nmode(nbins,0);
  
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
	      
	      const int bin = std::floor(std::sqrt(global_i*global_i+global_j*global_j+global_k*global_k));

	      nmode[bin]++;
	      spectrum[bin] += std::pow(std::abs(a[local_indx_threaded]),2);
	      mode[bin] += std::sqrt(k2);
	    }
	}
      }
      
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
	std::string fileName("spectra/spectrum_");
	std::stringstream stream;
	stream << std::fixed << std::setprecision(5) << time;
	fileName += stream.str();
	std::ofstream spectrum_out(fileName.c_str());
	for(int bin=0; bin < spectrum.size(); ++bin) {
	  spectrum[bin] *= time*time/(geom.ProbHi(0)-geom.ProbLo(0))/(geom.ProbHi(1)-geom.ProbLo(1))/(geom.ProbHi(2)-geom.ProbLo(2))*pow(h,6)/2.0;
	  spectrum_out << mode[bin] << " " << spectrum[bin] << " " << nmode[bin] << std::endl;
	}
	spectrum_out.close();
      }
}
#endif
#endif
