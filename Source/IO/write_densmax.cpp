#ifdef FDM
#include <iomanip>
#include <Nyx.H>

using namespace amrex;

void
Nyx::write_densmax ()
{
    int ndatalogs = parent->NumDataLogs();
    Real time_unit = 1.0;// 3.0856776e19 / 31557600.0; // conversion to Julian years

    int rlp = Nyx::runlog_precision;
    //    std::cout << "ndatalogs" << ndatalogs << std::endl;
    if (ndatalogs > 1)
    {

#ifdef FDM
      Real max_dens = 0.0; // , max_dens_nb = 0.0;
      int lev_finest = parent->finestLevel();
      if (lev_finest){
          max_dens = get_level(lev_finest).get_new_data(Axion_Type).max(Nyx::AxDens);
      }

#endif



#ifdef NO_HYDRO
#ifdef FDM
#ifdef GRAVITY
	Real time = state[PhiGrav_Type].curTime();
#else
	Real time = state[Axion_Type].curTime();
#endif
#else
	Real time = state[PhiGrav_Type].curTime();
#endif
#else
        Real time  = state[State_Type].curTime();
#endif
        Real dt    = parent->dtLevel(0);
        int  nstep = parent->levelSteps(0);

        if (ParallelDescriptor::IOProcessor())
        {
            std::ostream& data_loga = parent->DataLog(1);

            if (time == 0.0)
            {
                data_loga << std::setw( 8) <<  "#  nstep";
                data_loga << std::setw(14) <<  "       time    ";
                data_loga << std::setw(14) <<  "       dt      ";
                data_loga << std::setw(14) <<  "         z     ";
#ifdef FDM
		data_loga << std::setw(22) <<  " max_dens";
#endif

                data_loga << '\n';

                Real old_z = (1. / old_a) - 1.;
                data_loga << std::setw( 8) <<  nstep;
                data_loga << std::setw(14) <<  std::setprecision(rlp) <<  (initial_time+time) * time_unit;
                data_loga << std::setw(14) <<  std::setprecision(rlp) <<    dt * time_unit;
                data_loga << std::setw(14) <<  std::setprecision(rlp) << old_z;
#ifdef FDM
		data_loga << std::setw(22) <<  std::setprecision(6) << max_dens;
#endif
                data_loga << '\n';
            }
            else
            {
                const Real new_z = (1. / new_a) - 1.;
                data_loga << std::setw( 8) <<  nstep;
                data_loga << std::setw(14) <<  std::setprecision(rlp) <<  (initial_time+time) * time_unit;
                data_loga << std::setw(14) <<  std::setprecision(rlp) <<    dt * time_unit;
                data_loga << std::setw(14) <<  std::setprecision(rlp) << new_z;
#ifdef FDM
		data_loga << std::setw(22) <<  std::setprecision(10) << max_dens;

#endif

                data_loga << std::endl;
            }
        }
    }
}
#endif
