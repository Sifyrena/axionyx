
#include <Nyx.H>
#ifdef GRAVITY
#include <Gravity.H>
#endif
#include <constants_cosmo.H>
#include <Forcing.H>

using namespace amrex;

using std::string;

Real
Nyx::advance (Real time,
              Real dt,
              int  iteration,
              int  ncycle)

  // Arguments:
  //    time      : the current simulation time
  //    dt        : the timestep to advance (e.g., go from time to
  //                time + dt)
  //    iteration : where we are in the current AMR subcycle.  Each
  //                level will take a number of steps to reach the
  //                final time of the coarser level below it.  This
  //                counter starts at 1
  //    ncycle    : the number of subcycles at this level

{

  const std::string region_name = nStep() > 1 ? "R::Nyx::advance" : "R::Nyx::advance::STEP1";
  BL_PROFILE_REGION(region_name);
  MultiFab::RegionTag amrlevel_tag("AmrLevel_Level_" + std::to_string(level));

#ifndef NO_HYDRO // If Hydro!
  MultiFab& S_new = get_new_data(State_Type);
  if (S_new.contains_nan(0, 1, 0))                                                                
    {amrex::Abort("Watchdog: New Hydro State has NaNs");
    } 
#endif
    MultiFab& Ax_new = get_new_data(Axion_Type);

  {
    const char* axnames[] = {"AxDens","AxRe","AxIm","AxPhas"};
    bool found = false;
    for (int c = 0; c < Ax_new.nComp(); ++c) {
        if (Ax_new.contains_nan(c, 1, 0)) {
            amrex::Print() << "Watchdog: NaN in Axion_Type component " << c
                           << " (" << axnames[c] << ") at level " << level
                           << " step " << nStep() << "\n";
            found = true;
        }
    }
    if (found) amrex::Abort("Watchdog: New FDM State has NaNs.");
  }

#ifndef NO_HYDRO
return advance_hydro_axionyx(time, dt, iteration, ncycle);
#else
if (!do_hydro)
{
  return advance_particles_only(time, dt, iteration, ncycle);
}
else
{
return advance_axionyx(time, dt, iteration, ncycle);
}
#endif
}