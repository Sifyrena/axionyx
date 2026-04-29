#ifdef AXREAL
#include <AMReX_BLProfiler.H>
#include <Nyx.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

void Nyx::initialize_axion_only_simulation (std::string readin_ics_fname)
{

  MultiFab& Ax_new = get_new_data(Axion_Type);

  if (level == 0){

    Nyx::ax_init_time = getInitialTimeFromHeader(readin_ics_fname);
    setInitialTime();
    MultiFab mf;
    icReadAndPrepareFab(readin_ics_fname, 0, mf);
    convertComplexAxionToRealAxion(mf, Ax_new);
    mf.clear();
    mend_axions();
    if (ParallelDescriptor::IOProcessor())
      std::cout << "Loading Jaxions initial conditions on root grid ...done\n";

  }else{

    FillCoarsePatch(Ax_new, 0, Nyx::ax_init_time, Axion_Type, 0, Ax_new.nComp());
    if (ParallelDescriptor::IOProcessor())
      std::cout << "Interpolating Jaxions initial conditions onto level "<<level<<" ...done\n";

  }
}


void Nyx::convertComplexAxionToRealAxion (amrex::MultiFab& mf_src, amrex::MultiFab& mf_state)
{

  BL_PROFILE("Nyx::initialize_axion_only_simulation()");

  int AxRe = 0;
  int AxIm = 1;
  int AxvRe = 2;
  int AxvIm = 3;
  amrex::Real time = Nyx::ax_init_time;

  for (MFIter mfi(mf_state,TilingIfNotGPU()); mfi.isValid(); ++mfi){
    Array4<Real> const& arr_src = mf_src.array(mfi);
    Array4<Real> const& arr_state = mf_state.array(mfi);
    const Box& bx = mfi.tilebox();
    ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                  arr_state(i,j,k,Nyx::Psi) = time*std::atan2(arr_src(i,j,k,AxIm),arr_src(i,j,k,AxRe));
                  arr_state(i,j,k,Nyx::PsiPrime) = time*(arr_src(i,j,k,AxRe)*arr_src(i,j,k,AxvIm)-arr_src(i,j,k,AxIm)*arr_src(i,j,k,AxvRe))
		    / (arr_src(i,j,k,AxRe)*arr_src(i,j,k,AxRe)+arr_src(i,j,k,AxIm)*arr_src(i,j,k,AxIm))
		    + std::atan2(arr_src(i,j,k,AxIm),arr_src(i,j,k,AxRe));
                });
  }
}


Real Nyx::getInitialTimeFromHeader (std::string readin_ics_fname)
{
  std::string mfDirName(readin_ics_fname);
  int LineInWhichTimeIsStored = 8;
  std::ifstream file(mfDirName+"/Header");
  std::string str_time;
  for (int i = 1; i <= LineInWhichTimeIsStored; i++)
    std::getline(file, str_time);
  return float(std::floor(10000.0*stof(str_time)))/10000.0;
}


void Nyx::setInitialTime ()
{
  parent->setCumTime(Nyx::ax_init_time);
  old_a_time = Nyx::ax_init_time;
  new_a_time = Nyx::ax_init_time;
  // We set dt to be large for this new level to avoid screwing up computeNewDt.
  Real dummy_dt = 1.e100;
  setTimeLevel(parent->cumTime(), dummy_dt, dummy_dt);
}


void Nyx::mend_axions ()
{

  MultiFab& Ax_new = get_level(0).get_new_data(Axion_Type);
  amrex::Real time = Nyx::ax_init_time;
  const Real pi = 4 * std::atan(1.0);
  const Real tpi = 2 * pi;

  for (FillPatchIterator fpi(*this, Ax_new, 1, time, Axion_Type, 0, 1); fpi.isValid(); ++fpi)
    {
      const Box& bx  = fpi.validbox();
      Array4<Real> const& arr_in   = fpi().array();
      Array4<Real> const& arr_out  = Ax_new[fpi].array();
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);
      for (int k = lo.z; k <= hi.z; ++k) {
	for (int j = lo.y; j <= hi.y; ++j) {
	  for (int i = lo.x; i <= hi.x; ++i) {
	    if(arr_in(i-1,j,k,0)-arr_in(i,j,k,0)>5.5) {
	      arr_out(i,j,k,0) += time*tpi;
	      arr_out(i,j,k,1) += tpi;
	    }
	  }
	}
      }
      
    }
}
#endif
