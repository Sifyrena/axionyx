#ifdef AXCOMPLEX
#include <AMReX_BLProfiler.H>
#include <Nyx.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

void Nyx::initialize_axion_string_simulation (std::string readin_ics_fname)
{

  MultiFab& Ax_new = get_new_data(Axion_Type);

  if (level == 0){

    Nyx::string_init_time = getInitialTimeFromHeader(readin_ics_fname);
    setInitialTime();
    icReadAndPrepareFab(readin_ics_fname, 0, Ax_new);
    if (ParallelDescriptor::IOProcessor())
      std::cout << "Loading Jaxions initial conditions on root grid ...done\n";

  }else{

    FillCoarsePatch(Ax_new, 0, Nyx::string_init_time, Axion_Type, 0, Ax_new.nComp());
    if (ParallelDescriptor::IOProcessor())
      std::cout << "Interpolating Jaxions initial conditions onto level "<<level<<" ...done\n";    

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
  parent->setCumTime(Nyx::string_init_time);
  old_a_time = Nyx::string_init_time;
  new_a_time = Nyx::string_init_time;
  // We set dt to be large for this new level to avoid screwing up computeNewDt.
  Real dummy_dt = 1.e100;
  setTimeLevel(parent->cumTime(), dummy_dt, dummy_dt);
}
#endif
