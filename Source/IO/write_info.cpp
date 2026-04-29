#include <iomanip>
#include <Nyx.H>

using namespace amrex;

void
Nyx::write_info ()
{
    int ndatalogs = parent->NumDataLogs();
    Real time_unit = 1.0; //3.0856776e19 / 31557600.0; // conversion to Julian years

    int rlp  = Nyx::runlog_precision;
    int rlpt = amrex::min(rlp,Nyx::runlog_precision_terse);

    if (ndatalogs > 0)
    {
#ifndef NO_HYDRO
        MultiFab& S_new = get_new_data(State_Type);
        MultiFab& D_new = get_new_data(DiagEOS_Type);
        Real      max_t = 0;

        Real rho_T_avg=0.0, T_avg=0.0, Tinv_avg=0.0, T_meanrho=0.0;
        Real whim_mass_frac, whim_vol_frac, hh_mass_frac, hh_vol_frac, igm_mass_frac, igm_vol_frac;

        if (do_hydro)
        {
            // Removed reset internal energy before call to compute_temp, still compute new temp
            compute_new_temp(S_new,D_new);
            max_t = D_new.norm0(Temp_comp);
            compute_rho_temp(rho_T_avg, T_avg, Tinv_avg, T_meanrho);
            compute_gas_fractions(1.0e5, 120.0, whim_mass_frac, whim_vol_frac,
                                  hh_mass_frac, hh_vol_frac, igm_mass_frac, igm_vol_frac);
        }
#endif

#ifdef AXCOMPLEX
        Real stringdens;
        Vector<Vector<Real> > test_val;
        compute_axionyx_quantities(stringdens, test_val);
        MultiFab&  Ax_new = get_level(0).get_new_data(Axion_Type);
#endif
#ifdef FDM
	MultiFab&  Ax_new = get_level(0).get_new_data(Axion_Type);
        Real mass, epot, ekinrho, ekinv, ekin, etot, angmom_x, angmom_y, angmom_z, max_dens;
        compute_axionyx_quantities(mass, epot, ekinrho, ekinv, ekin, etot, angmom_x, angmom_y, angmom_z, max_dens);
#endif
#ifndef NO_HYDRO
        Real Bmass, Bepot, Bekin, Beint, Betot;
        compute_baryon_quantities(Bmass, Bepot, Bekin, Beint, Betot);
#endif
        
        Real time  = state[State_for_Time].curTime();
        Real dt    = parent->dtLevel(0);
        int  nstep = parent->levelSteps(0);

        if (ParallelDescriptor::IOProcessor())
        {
            std::ostream& data_loga = parent->DataLog(0);
#ifdef AXCOMPLEX
            if (time == Nyx::string_init_time)
#elif AXREAL
	    if (time == Nyx::ax_init_time)
#else
            if (time == 0.0)
#endif
            {
                data_loga << std::setw(8) <<  "Step#"; 
                data_loga << std::setw(25) <<  "TIME";
                data_loga << std::setw(25) <<  "dt";
#if !defined(AXCOMPLEX) && !defined(AXREAL)
                data_loga << std::setw(8) <<  "z";
#endif
#ifdef AXCOMPLEX
                data_loga << std::setw(40) <<  " string_dens";
                for(int i=0; i<Nyx::test_posx.size(); i++){
                  data_loga << std::setw(25) <<  "    test_mRe   ";
                  data_loga << std::setw(25) <<  "    test_mIm   ";
                  data_loga << std::setw(25) <<  "    test_vRe   ";
                  data_loga << std::setw(25) <<  "    test_vIm   ";
                }
#endif
#ifdef FDM
                data_loga << std::setw(25) <<  "A_Mass";
                data_loga << std::setw(25) <<  "A_Epot";
                data_loga << std::setw(25) <<  "A_Eqp";
                data_loga << std::setw(25) <<  "A_Ekinv";
                data_loga << std::setw(25) <<  "A_Ekin";
                data_loga << std::setw(25) <<  "A_Etot";
                data_loga << std::setw(40) <<  "A_rhomax";
                data_loga << std::setw(40) <<  "A_GBrmax";
                if(fdm_halo)
                  {
                    data_loga << std::setw(40) <<  "fdm_halo_x";
                    data_loga << std::setw(40) <<  "fdm_halo_y";
                    data_loga << std::setw(40) <<  "fdm_halo_z";
                  }
#endif
#if !defined(AXCOMPLEX) && !defined(AXREAL)
                data_loga << std::setw(25) <<  "      a        ";
#endif
#ifndef NO_HYDRO // YES HYDRO
                if (do_hydro == 1)
                {
                  data_loga << std::setw(25) <<  "T_max";
                  data_loga << std::setw(25) <<  "<T>_rho";
                  data_loga << std::setw(25) <<  "<T>_V";
                  data_loga << std::setw(25) <<  "T @ <rho>";
                  data_loga << std::setw(25) <<  "T(21cm)";
                  data_loga << std::setw(25) <<  "B_Eadiab.";
                  data_loga << std::setw(25) <<  "WHIM_m";
                  data_loga << std::setw(25) <<  "WHIM_v";
                  data_loga << std::setw(25) <<  "HH_m";
                  data_loga << std::setw(25) <<  "HH_v";
                  data_loga << std::setw(25) <<  "IGM_m";
                  data_loga << std::setw(25) <<  "IGM_v";
                  data_loga << std::setw(25) <<  "B_Mass";
                  data_loga << std::setw(25) <<  "B_Epot";
                  data_loga << std::setw(25) <<  "B_Ekin";
                  data_loga << std::setw(25) <<  "B_Eint";
                  data_loga << std::setw(25) <<  "B_Etot";
                }
#endif
                data_loga << '\n';

                Real old_z = (1. / old_a) - 1.;
                data_loga << std::setw( 8) <<  nstep;
                data_loga << std::setw(25) <<  std::setprecision(rlp) <<  (initial_time+time) * time_unit;
                data_loga << std::setw(25) <<  std::setprecision(rlp) <<    dt * time_unit;
#if !defined(AXCOMPLEX) && !defined(AXREAL)
                data_loga << std::setw(8) <<  std::setprecision(rlp) << old_z;
#endif
#ifdef AXCOMPLEX
                data_loga << std::setw(40) <<  std::setprecision(6) << stringdens;
                for(int i=0; i<Nyx::test_posx.size(); i++){
                  for(int n=0 ; n<Ax_new.nComp(); n++)
                    data_loga << std::setw(25) <<  std::setprecision(6) << test_val[i][n];
                }
#endif
#ifdef FDM
                data_loga << std::setw(25) <<  std::setprecision(rlp) << mass;
                data_loga << std::setw(25) <<  std::setprecision(rlp) << epot;
                data_loga << std::setw(25) <<  std::setprecision(rlp) << ekinrho;
                data_loga << std::setw(25) <<  std::setprecision(rlp) << ekinv;
                data_loga << std::setw(25) <<  std::setprecision(rlp) << ekin;
                data_loga << std::setw(25) <<  std::setprecision(rlp) << etot;
                data_loga << std::setw(40) <<  std::setprecision(rlp) << max_dens;
                if(fdm_halo)
                  {
                    data_loga << std::setw(40) <<  std::setprecision(6) << fdm_halo_pos[0];
                    data_loga << std::setw(40) <<  std::setprecision(6) << fdm_halo_pos[1];
                    data_loga << std::setw(40) <<  std::setprecision(6) << fdm_halo_pos[2];
                  }
#endif
#if !defined(AXCOMPLEX) && !defined(AXREAL)
                data_loga << std::setw(25) <<  std::setprecision(rlp) << old_a;
#endif
#ifndef NO_HYDRO
                if (do_hydro == 1)
                {
                   data_loga << std::setw(25) <<  std::setprecision(rlpt) << max_t;
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << rho_T_avg;
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << T_avg;
                   data_loga << std::setw(25) <<  std::setprecision(rlpt) << T_meanrho;
                   data_loga << std::setw(25) <<  std::setprecision(rlpt) << 1.0/Tinv_avg;
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << 0.021*(1.0+old_z)*(1.0+old_z);
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << whim_mass_frac;
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << whim_vol_frac;
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << hh_mass_frac;
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << hh_vol_frac;
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << igm_mass_frac;
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << igm_vol_frac;
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << Bmass;
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << Bepot;
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << Bekin;
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << Beint;
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << Betot;
                }
#endif
                data_loga << '\n';
            }
            else
            {
                const Real new_z = (1. / new_a) - 1.;
                data_loga << std::setw( 8) <<  nstep;
                data_loga << std::setw(25) <<  std::setprecision(rlp) <<  (initial_time+time) * time_unit;
                data_loga << std::setw(25) <<  std::setprecision(rlp) <<    dt * time_unit;
#if !defined(AXCOMPLEX) && !defined(AXREAL)
                data_loga << std::setw(8) <<  std::setprecision(rlp) << new_z;
#endif
#ifdef AXCOMPLEX
                data_loga << std::setw(40) <<  std::setprecision(6) << stringdens;
                for(int i=0; i<Nyx::test_posx.size(); i++){
                  for(int n=0 ; n<Ax_new.nComp(); n++)
                    data_loga << std::setw(25) <<  std::setprecision(6) << test_val[i][n];
                }
#endif
#ifdef FDM
                data_loga << std::setw(25) <<  std::setprecision(rlp) << mass;
                data_loga << std::setw(25) <<  std::setprecision(rlp) << epot;
                data_loga << std::setw(25) <<  std::setprecision(rlp) << ekinrho;
                data_loga << std::setw(25) <<  std::setprecision(rlp) << ekinv;
                data_loga << std::setw(25) <<  std::setprecision(rlp) << ekin;
                data_loga << std::setw(25) <<  std::setprecision(rlp) << etot;
                data_loga << std::setw(40) <<  std::setprecision(rlp) << max_dens;
                if(fdm_halo)
                  {
                    data_loga << std::setw(40) <<  std::setprecision(6) << fdm_halo_pos[0];
                    data_loga << std::setw(40) <<  std::setprecision(6) << fdm_halo_pos[1];
                    data_loga << std::setw(40) <<  std::setprecision(6) << fdm_halo_pos[2];
                  }
#endif
#if !defined(AXCOMPLEX) && !defined(AXREAL)
                data_loga << std::setw(25) <<  std::setprecision(rlp) << new_a;
#endif
#ifndef NO_HYDRO
                if (do_hydro == 1)
                {
                   data_loga << std::setw(25) <<  std::setprecision(rlpt) << max_t;
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << rho_T_avg;
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << T_avg;
                   data_loga << std::setw(25) <<  std::setprecision(rlpt) << T_meanrho;
                   data_loga << std::setw(25) <<  std::setprecision(rlpt) << 1.0/Tinv_avg;
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << 0.021*(1.0+new_z)*(1.0+new_z);
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << whim_mass_frac;
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << whim_vol_frac;
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << hh_mass_frac;
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << hh_vol_frac;
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << igm_mass_frac;
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << igm_vol_frac;
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << Bmass;
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << Bepot;
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << Bekin;
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << Beint;
                   data_loga << std::setw(25) <<  std::setprecision(rlp) << Betot;
                }
#endif
                data_loga << std::endl;
            }
        }
    }
}
