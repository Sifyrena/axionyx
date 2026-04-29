#include <iomanip>
#include <Nyx.H>
#include <Prob.H>
#ifdef GRAVITY
#include <Gravity.H>
#endif
using namespace amrex;

#include <fstream>
#include <vector>

amrex::Vector<double> read_binary_to_buffer(const std::string& filename, const int N)
{
    amrex::Vector<double> buffer(N, 0.0); // Allocate buffer on all ranks

   if (amrex::ParallelDescriptor::MyProc() == 0) {
        //Only Rank 0 reads the file
        std::ifstream file(filename, std::ios::binary);
        if (!file.is_open()) {
        amrex::Print() << "Error: Unable to open file " << filename << std::endl;
        } else {
        file.read(reinterpret_cast<char*>(buffer.data()), N * sizeof(double));
        if (!file) {
            amrex::Print() << "Error: Failed to read the entire binary file." << std::endl;
            }
            file.close();
            amrex::Print() << "Praise Bodo." << std::endl;
        }
    } // Do not need anymore!

    // Broadcast the buffer to all ranks
    amrex::ParallelDescriptor::Bcast(buffer.data(), N, 0);
    amrex::Print()<<"Size = "<<N<<std::endl;
    amrex::ParallelDescriptor::Bcast(buffer.data(), N, 0);
    return buffer;
}


namespace
{
    std::string ascii_particle_file;
    std::string binary_particle_file;
    std::string    sph_particle_file;

#ifdef AGN
    std::string agn_particle_file;
#endif

#ifdef NEUTRINO_PARTICLES
    std::string neutrino_particle_file;
#endif
}

static int  do_santa_barbara = 0;
static int  init_sb_vels     = 1;
static int  do_readin_ics    = 0;
static int  fix_random_seed  = 0;
std::string readin_ics_fname;

#ifndef NO_HYDRO
void Nyx::set_small_values_given_average (amrex::Real average_dens, amrex::Real average_temp, amrex::Real a,
                                          amrex::Real & small_dens_inout, amrex::Real & small_temp_inout, 
                                          amrex::Real &small_pres_inout, Real gamma_minus_1_in, 
                                          Real h_species_in)
{

    Real lsmall_dens, lsmall_temp, frac;
    frac = 1e-6;
    if (small_dens_inout <= 0.e0)
        lsmall_dens = frac * average_dens;
    else
        lsmall_dens = small_dens_inout;

    if (small_temp_inout <= 0.e0)
        lsmall_temp = frac * average_temp;
    else
        lsmall_temp = small_temp_inout;

    ParmParse pp_nyx("nyx");

    // We approximate Ne = 1.0 for consistency with iterate_ne in EOS
    const amrex::Real typical_Ne = 1.e0;
    Real lsmall_pres=0.0;
    Real dummy_small_pres=0.0;
    Real eint=0.0;
    auto atomic_rates = atomic_rates_glob;

    if (!(pp_nyx.query("small_pres", dummy_small_pres)))
        nyx_eos_given_RT(atomic_rates, gamma_minus_1_in, h_species_in, &eint, 
                         &lsmall_pres, lsmall_dens, lsmall_temp, typical_Ne, a);
    else
        lsmall_pres = small_pres_inout;

    small_dens_inout = lsmall_dens;
    small_temp_inout = lsmall_temp;
    small_pres_inout = lsmall_pres;
}
#endif // No Hydro

void
Nyx::read_init_params ()
{
    BL_PROFILE("Nyx::read_init_params()");
    ParmParse pp("nyx");

    pp.query("do_santa_barbara", do_santa_barbara);
    pp.query("init_sb_vels", init_sb_vels);

    pp.query("fix_random_seed", fix_random_seed);
    // Note that the value of 1024UL is not significant -- the point here is just to set the
    //      same seed for all MPI processes for the purpose of regression testing
    if (fix_random_seed)
        amrex::InitRandom(1024UL);

    if (do_hydro == 0 && do_santa_barbara == 1)
           amrex::Error("Nyx::cant have do_hydro == 0 and do_santa_barbara == 1");
    if (do_santa_barbara == 0 && init_with_sph_particles == 1)
           amrex::Error("Nyx::cant have do_santa_barbara == 0 and init_with_sph_particles == 1");
    if (do_santa_barbara == 0 && init_sb_vels == 1)
    {
       init_sb_vels = 0;
       if (ParallelDescriptor::IOProcessor())
           std::cout << "Nyx::setting init_sb_vels to 0 since do_santa_barbara = 0\n";
    }

    pp.query("do_readin_ics",       do_readin_ics);
    pp.query("readin_ics_fname", readin_ics_fname);
#ifdef AMREX_PARTICLES  
    pp.query("ascii_particle_file", ascii_particle_file);

    // Input error check
    if (do_dm_particles && !ascii_particle_file.empty() && particle_init_type != "AsciiFile")
    {
        if (ParallelDescriptor::IOProcessor())
            std::cerr << "ERROR::particle_init_type is not AsciiFile but you specified ascii_particle_file" << std::endl;
        amrex::Error();
    }

    pp.query("sph_particle_file", sph_particle_file);

    // Input error check
    if (init_with_sph_particles != 1 && !sph_particle_file.empty())
    {
        if (ParallelDescriptor::IOProcessor())
            std::cerr << "ERROR::init_with_sph_particles is not 1 but you specified sph_particle_file" << std::endl;
        amrex::Error();
    }

    // Input error check
    if (init_with_sph_particles == 1 && sph_particle_file.empty())
    {
        if (ParallelDescriptor::IOProcessor())
            std::cerr << "ERROR::init_with_sph_particles is 1 but you did not specify sph_particle_file" << std::endl;
        amrex::Error();
    }

    pp.query("binary_particle_file", binary_particle_file);

    // Input error check
    if (!binary_particle_file.empty() && (particle_init_type != "BinaryFile" &&
                                          particle_init_type != "BinaryMetaFile" && 
                                          particle_init_type != "BinaryMortonFile"))
    {
        if (ParallelDescriptor::IOProcessor())
            std::cerr << "ERROR::particle_init_type is not BinaryFile, BinaryMetaFile, or BinaryMortonFile but you specified binary_particle_file" << std::endl;
        amrex::Error();
    }

#ifdef AGN
    pp.query("agn_particle_file", agn_particle_file);
    if (!agn_particle_file.empty() && particle_init_type != "AsciiFile")
    {
        if (ParallelDescriptor::IOProcessor())
            std::cerr << "ERROR::particle_init_type is not AsciiFile but you specified agn_particle_file" << std::endl;
        amrex::Error();
    }
#endif

#ifdef NEUTRINO_PARTICLES
    pp.query("neutrino_particle_file", neutrino_particle_file);
    if (!neutrino_particle_file.empty() && particle_init_type != "AsciiFile" &&
                                           particle_init_type != "BinaryMetaFile" && 
                                           particle_init_type != "BinaryFile")
    {
        if (ParallelDescriptor::IOProcessor())
            std::cerr << "ERROR::particle_init_type is not AsciiFile or BinaryFile but you specified neutrino_particle_file" << std::endl;
        amrex::Error();
    }
#endif
#endif
}

void
Nyx::init_zhi ()
{
    BL_PROFILE("Nyx::init_zhi()");

    if (ParallelDescriptor::IOProcessor()) std::cout << "Reading z_HI from file...";

    const int file_res = inhomo_grid;
    const int prob_res = geom.Domain().longside();
    const int ratio = prob_res / file_res;

    BL_ASSERT(ratio >= 1);
    
#ifndef NO_HYDRO
    MultiFab& D_new = get_new_data(DiagEOS_Type);

    const BoxArray& my_ba = D_new.boxArray();
    const DistributionMapping& my_dmap = D_new.DistributionMap();

    BL_ASSERT(my_ba.coarsenable(ratio));
    BoxArray coarse_ba = my_ba;
    coarse_ba.coarsen(ratio);
    MultiFab zhi(coarse_ba, my_dmap, 1, 0);

    MultiFab zhi_from_file;
    VisMF::Read(zhi_from_file, inhomo_zhi_file);
    zhi.ParallelCopy(zhi_from_file, geom.periodicity());
    if(D_new.nComp()>2)
    {
        int l_Zhi_comp = Zhi_comp;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(D_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const auto fab_zhi=zhi.array(mfi);
                const auto fab_D_new=D_new.array(mfi);

                amrex::ParallelFor(
                               bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                                   fab_D_new(i, j ,k, l_Zhi_comp) = fab_zhi(i/ratio,j/ratio,k/ratio);
                               });

            }
    }
#endif
    if (ParallelDescriptor::IOProcessor()) std::cout << "done.\n";
}

void 
Nyx::initData ()
{
    BL_PROFILE("Nyx::initData()");
    
    // Here we initialize the grid data and the particles from a plotfile.
    if (!parent->theRestartPlotFile().empty())
    {
        amrex::Abort("AmrData requires fortran");
        return;
    }

#ifndef NO_HYDRO    
    MultiFab&   S_new  = get_new_data(State_Type);
    MultiFab&   D_new  = get_new_data(DiagEOS_Type);
    // We need this because otherwise we might operate on uninitialized data.
    S_new.setVal(0.0);
    D_new.setVal(0.0);
#endif

#ifdef AXIONYX
    MultiFab&   Ax_new   = get_new_data(Axion_Type);
#endif

    // If you run a pure N-body simulation and Nyx segfaults here, then
    // please check Prob_3d.f90... You might set other variables then density...
    if (verbose && ParallelDescriptor::IOProcessor())
        amrex::Print() << "Initializing the data at level " << level << '\n';

    const auto dx = geom.CellSizeArray();
    const auto geomdata = geom.data();

    // Make sure dx = dy = dz -- that's all we guarantee to support
    const Real SMALL = 1.e-13;
    if ( (fabs(dx[0] - dx[1]) > SMALL) || (fabs(dx[0] - dx[2]) > SMALL) )
        amrex::Abort("We don't support dx != dy != dz");

#ifndef NO_HYDRO    
#ifdef AMREX_PARTICLES
    if ( (do_santa_barbara == 0) && (do_readin_ics == 0) && (particle_init_type != "Cosmological") )
#else
    if ( (do_santa_barbara == 0) && (do_readin_ics == 0) )
#endif
    {
        if (do_hydro == 1) 
        {
            MultiFab&   D_new    = get_new_data(DiagEOS_Type);
            D_new.setVal(0., Temp_comp);
            D_new.setVal(0.,   Ne_comp);

    const Box& domain = geomdata.Domain();

    // Get the number of cells in each dimension for the full domain
    const int nx = domain.length(0); 
    const int ny = domain.length(1); 
    const int nz = domain.length(2); 
    
    GpuArray<amrex::Real,max_prob_param> prob_param;

    amrex::Vector<double> axionRebuffer = read_binary_to_buffer("AxionRe.bin", nx * ny * nz);
    amrex::Vector<double> axionImbuffer = read_binary_to_buffer("AxionIm.bin", nx * ny * nz);
    
    #ifdef _OPENMP
    #pragma omp parallel if (Gpu::notInLaunchRegion())
    #endif
        for (MFIter mfi(Ax_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<amrex::Real> const& axion = Ax_new.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
            int idx = (j * (nx) * (nx) + i * (nx) + k);
            axion(i,j,k, Nyx::AxRe) = axionRebuffer[idx];
            axion(i,j,k, Nyx::AxIm) = axionImbuffer[idx];
            
            axion(i,j,k, Nyx::AxDens) = axion(i,j,k, Nyx::AxIm) * axion(i,j,k, Nyx::AxIm) + axion(i,j,k, Nyx::AxRe) * axion(i,j,k, Nyx::AxRe);
            
            axion(i,j,k, Nyx::AxPhas) = atan2(axion(i,j,k, Nyx::AxIm), axion(i,j,k, Nyx::AxRe));
            });
        };
    amrex::Print()<<"Axion Data Loaded from Files."<<std::endl;
// The above block is now handled in Prob.cpp
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();
    prob_initdata_Basic_on_box(
        bx, S_new.array(mfi), D_new.array(mfi), Ax_new.array(mfi),
        geomdata, prob_param
    );
}


            if (inhomo_reion) init_zhi();

            // First reset internal energy before call to compute_temp
            MultiFab reset_e_src(S_new.boxArray(), S_new.DistributionMap(), 1, NUM_GROW);
            reset_e_src.setVal(0.0);

            reset_internal_energy(S_new,D_new,reset_e_src);
            compute_new_temp     (S_new,D_new);

            // Define (rho E) given (rho e) and the momenta
            enforce_consistent_e(S_new);


            if (S_new.contains_nan(0, 1, 0))                                                                                                                                                      
    {                                                                                                                                                                                                                                                                                                                                            
		amrex::Print() << "FW Init: New Hydro State has NaNs." << std::endl;                                                                                                                        
		amrex::Abort("S_new has NaNs in this component::advance_FDM_FFT()");                                                                                                                   
    } else {
        amrex::Print() << "FW Init: New Hydro State looks OK." << std::endl;   
    }     


        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const auto fab_S_new=S_new.array(mfi); 
                const auto fab_D_new=D_new.array(mfi);

#ifdef AXIONYX
                const auto fab_Ax_new=Ax_new.array(mfi);
#endif
                GpuArray<amrex::Real,max_prob_param> prob_param;
                prob_param_fill(prob_param);
                prob_param_special_fill(prob_param);
		comoving_type=int(std::round(prob_param[comoving_type_comp]));
#ifdef AXIONYX
		prob_initdata_state_on_box(bx, fab_S_new, fab_D_new, fab_Ax_new, geomdata, prob_param);
#else
         prob_initdata_state_on_box(bx, fab_S_new, fab_D_new, geomdata, prob_param); // Use Nyx for non=Axion Runs. the overheads might not be clean enough.
#endif
            }
        }
    }
#else // No Hydro Part
#ifdef AXIONYX
#ifdef AMREX_PARTICLES
    if ( (do_santa_barbara == 0) && (do_readin_ics == 0) && (particle_init_type != "Cosmological") )
#else
    if ( (do_santa_barbara == 0) && (do_readin_ics == 0) )
#endif
    {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
#ifdef EIGENIC
        static bool eigen_func_IC = false; 
        ParmParse ppn("nyx");
        ppn.query("eigen_func_IC", eigen_func_IC); // check if initialisation done with eigenfunctions 

        if(eigen_func_IC){
            amrex::Real pi = std::atan(1.0)*4.0;

            static int Nmax=0;
            static int Nr=0;

            // TODO: read these parameters from the header of the anl file!  
            amrex::Print() << "Initializing the data at level " << level << '\n';
            ppn.query("Nmax", Nmax); // get the Number of eigenfucntionsform the inputs
            amrex::Print() << "Eigenfucntion ICs: Number of eigenfucntions: " << Nmax <<std::endl;
            ppn.query("Nr", Nr); // get the Number of points in radius form the inputs
            amrex::Print() << "Number of points in radius: " << Nr << std::endl;

            amrex::Vector<amrex::Real> r;
            r.resize(Nr);
            amrex::Vector<amrex::Vector<int>> nl; 
            nl.resize(Nmax, amrex::Vector<int>(2)); 
            amrex::Vector<amrex::Real> anl;
            anl.resize(Nmax);
            amrex::Vector<amrex::Vector<amrex::Real>> Rnlr; 
            Rnlr.resize(Nmax, amrex::Vector<amrex::Real>(Nr)); 
            amrex::Vector<amrex::Vector<amrex::Real>> randnlm; 
            randnlm.resize(Nmax, amrex::Vector<amrex::Real>(2*Nmax +1)); 

            read_FDMhalo_eigen(Nmax, Nr, r, nl, anl, Rnlr);

            amrex::Print() << "Read eigenmodes Rnlrs and factors anl from a file.\n";

            int l;
            srand48(12345);
            for(int p = 0; p<Nmax; p++){
                l =nl[p][1];
                for(int m = -l; m<l+1; m++){
                    randnlm[p][m+l] = drand48()*2*pi;
                }
            }
            // amrex::Print() << "After random number \n";
        

            for (MFIter mfi(Ax_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const auto fab_Ax_new=Ax_new.array(mfi);
                GpuArray<amrex::Real,max_prob_param> prob_param;
                prob_param_fill(prob_param);
                prob_param_special_fill(prob_param);
                comoving_type=int(std::round(prob_param[comoving_type_comp]));
                prob_initdata_state_on_box_eigen(bx, fab_Ax_new, geomdata, prob_param, Nmax, Nr, r, nl, anl, Rnlr, randnlm);
            }
        }
        else
#endif //EIGENIC        
        { GpuArray<amrex::Real,max_prob_param> prob_param;
	       	const Box& domain = geomdata.Domain();

    // Get the number of cells in each dimension for the full domain
    const int nx = domain.length(0); 
    const int ny = domain.length(1); 
    const int nz = domain.length(2); 
    
    // Only loading Axion Data? We have that already .... 
amrex::Print() << "Initialising Data from Prob.cpp";	
{
            for (MFIter mfi(Ax_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const auto fab_Ax_new=Ax_new.array(mfi);
                GpuArray<amrex::Real,max_prob_param> prob_param;
                prob_param_fill(prob_param);
                prob_param_special_fill(prob_param);
                comoving_type=int(std::round(prob_param[comoving_type_comp]));
                prob_initdata_state_on_box(bx, fab_Ax_new, geomdata, prob_param);
            }
        }

	}
}
#endif
#endif // end NO_HYDRO

#ifdef GRAVITY
    if (!do_grav)
    {
        //
        // Set these to zero so they're defined for the plotfile.
        //
        MultiFab& G_new = get_new_data(Gravity_Type);
        G_new.setVal(0);
    }
    else 
    {
        //
        // Initialize this to zero before first solve.
        //
        MultiFab& Phi_new = get_new_data(PhiGrav_Type);
        Phi_new.setVal(0.);
    }
#endif

#ifndef NO_HYDRO

#ifdef SDC
    //
    // Initialize this to zero before we use it in advance
    //
    if (do_hydro)
    {
        MultiFab& IR_new = get_new_data(SDC_IR_Type);
        IR_new.setVal(0.0);
    }
#endif

    //
    // Read in initial conditions from a file.
    // By now only for fixed grid ics.
    // Layout and units have to be as in \vec U.
    //
    if (do_readin_ics)
    {
        std::string mfDirName(readin_ics_fname);

        MultiFab mf;
        mfDirName.append("/Level_0/Cell");

        VisMF::Read(mf, mfDirName.c_str());

        MultiFab& S_new_crse = get_level(0).get_new_data(State_Type);
        
        S_new_crse.MultiFab::ParallelCopy(mf, 0, 0, 6, 0, 0);
#ifndef CONST_SPECIES
        S_new_crse.MultiFab::ParallelCopy(mf, 0, FirstSpec_comp, 1, 0, 0);
#endif

        if (do_hydro == 1) 
        {
            MultiFab&  D_new = get_new_data(DiagEOS_Type);
            D_new.setVal(0.);
        }

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Readin stuff...done\n";
    }
#endif

#ifdef AXCOMPLEX
    if (do_readin_ics)
	initialize_axion_string_simulation (readin_ics_fname);
#endif
#ifdef AXREAL
    if (do_readin_ics)
	initialize_axion_only_simulation (readin_ics_fname);
#endif

    amrex::Gpu::Device::synchronize();

#ifdef AMREX_PARTICLES
    if (level == 0)
        init_particles();
#endif 

    amrex::Gpu::Device::synchronize();

#if defined(AMREX_PARTICLES) || defined(FDM)
    if ( particle_init_type == "Cosmological")
      initcosmo();
#ifdef NO_HYDRO
    if ( particle_init_type == "FWMDM")
      initFDMcosmo();
#endif
#endif

    //
    // Must redistribute particles before calling `init_santa_barbara` so that
    // the particles already live on the higher level when we go to put some of
    // the mass onto the grid.
    //
#ifdef AMREX_PARTICLES
    if (level > 0)
        particle_redistribute();
#endif

    //
    // With this call we define the initial data on the current level but we
    // also may need to modify the data on the level below since particles
    // previous at the coarser level may now live at the finer level and
    // distribute their mass differently.
    //
#ifndef NO_HYDRO
    if (do_santa_barbara == 1)
        init_santa_barbara(init_sb_vels);

    if (do_hydro)
        check_initial_species();

    //
    // Need to compute this in case we want to use overdensity for regridding.
    //
    if (level == 0) 
        compute_average_density();
#endif

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Done initializing the level " << level << " data\n";
}

#ifdef AXIONYX
#ifdef EIGENIC 
void Nyx::read_FDMhalo_eigen(int Nmax, 
                        int Nr, 
                        amrex::Vector<amrex::Real> &r, 
                        amrex::Vector<amrex::Vector<int>> &nl, 
                        amrex::Vector<amrex::Real> &anl, 
                        amrex::Vector<amrex::Vector<amrex::Real>> &Rnlr)
{

    amrex::Real tempvar;
    amrex::Real n;
    amrex::Real l;

    std::ifstream infile;
    std::string eigen_func_input_name = "";
    ParmParse pp("nyx");
    pp.query("eigen_func_input_name", eigen_func_input_name); // get the name of the file form the inputs
    amrex::Print() << eigen_func_input_name << std::endl;
    infile.open(eigen_func_input_name);//open the text file

    if (!infile) {
        amrex::Error("Can't open the eigenfunctions data file in Nyx::read_FDMhalo_eigen.");
    }

    for (int i = 0; i < Nr; i++){
        infile >> tempvar;
        //read r and convert from kpc to Mpc. 
        r[i]= tempvar*0.001;
        // amrex::Print() << r[i] << '\n';
    }
    // for (int n = 0; n < Nmax; n++){
    //     for (int l=0; l < Nmax; l++){

    for(int p = 0; p<Nmax; p++){
        infile >> n;
        nl[p][0] = int(n);
        infile >> l;
        nl[p][1] = int(l);
        infile >> tempvar;
        anl[p] = tempvar; 
        // std::cout << "n, l anl" << n << " " << l << " " << anl[p]<< std::endl;
        for (int ir = 0; ir < Nr; ir++){
            infile >> tempvar;
            Rnlr[p][ir] = tempvar;
        }
    }
    infile.close();
}
#endif //EIGENIC
#endif //AXIONYX

#if 0
void
Nyx::init_from_plotfile ()
{
    BL_PROFILE("Nyx::init_from_plotfile()");
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << " " << std::endl; 
        std::cout << "Initializing the data from " << parent->theRestartPlotFile() << std::endl;
    }

    if (parent->maxLevel() > 0)
        amrex::Abort("We can only restart from single-level plotfiles");

    // Make sure to read in "a" before we call ReadPlotFile since we will use a 
    //      when we construct e from T.

    bool is_checkpoint = false;

    // Now read in the time as well as the grid and particle data
    bool first = true;
    bool rhoe_infile;
    ReadPlotFile(first,parent->theRestartPlotFile(),rhoe_infile);

    // This is just a dummy value so that we can set the current time of the StateData
    Real dummy_dt = 1.e100;
    setTimeLevel(parent->cumTime(), dummy_dt, dummy_dt);

    comoving_a_post_restart(parent->theRestartPlotFile());

#ifndef NO_HYDRO
    for (int lev = 0; lev <= parent->finestLevel(); ++lev)
    {
        Nyx& nyx_lev = get_level(lev);
        MultiFab& S_new = nyx_lev.get_new_data(State_Type);
        MultiFab& D_new = nyx_lev.get_new_data(DiagEOS_Type);
        int ns = S_new.nComp();
        int nd = D_new.nComp();

        // Construct internal energy given density, temperature and species
        if (! rhoe_infile)
           init_e_from_T(old_a);

        // Define (rho E) given (rho e) and the momenta
        nyx_lev.enforce_consistent_e(S_new);
    }

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Done initializing the grid data from " << parent->theRestartPlotFile() << std::endl;
#endif

    // Now read the particles from the plotfile
    particle_post_restart(parent->theRestartPlotFile(),is_checkpoint);

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << "Done initializing the particles from the plotfile " << std::endl;
        std::cout << " " << std::endl; 
    }
}
#endif

#ifndef NO_HYDRO
void
Nyx::check_initial_species ()
{
    // Verify that the sum of (rho X)_i = rho at every cell.
#ifdef CONST_SPECIES
#ifndef AMREX_USE_FLOAT
    if (amrex::Math::abs(1.0 - h_species - he_species) > 1.e-8)
        amrex::Abort("Error:: Failed check of initial species summing to 1, A");
#else
    if (amrex::Math::abs(1.0 - h_species - he_species) > 1.e-6)
        amrex::Abort("Error:: Failed check of initial species summing to 1, B");
#endif
#else
    ReduceOps<ReduceOpMax> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    int iden  = Density_comp;
    if (FirstSpec_comp > 0)
    {
        int iufs = FirstSpec_comp;
        int nspec = NumSpec;
        MultiFab&   S_new    = get_new_data(State_Type);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          const auto state_fab = S_new.array(mfi);

          const Box& tbx = mfi.tilebox();

          reduce_op.eval(tbx, reduce_data, [state_fab,nspec,iden,iufs]
          AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
          {
               Real sum = state_fab(i,j,k,iufs);
               for (int n = 1; n < nspec; n++)
                  sum += state_fab(i,j,k,iufs+n);

               sum /= state_fab(i,j,k,iden);

               Real x = amrex::Math::abs(amrex::Math::abs(sum) - 1.);
               return x;
          });
        }

        ReduceTuple hv = reduce_data.value();
        ParallelDescriptor::ReduceRealMax(amrex::get<0>(hv));
#ifndef AMREX_USE_FLOAT
        std::cout << "Pre-FW INITCONSTSP HERE "<< amrex::Math::abs(get<0>(hv) - 1) * 1e8 << std::endl;
        // if (amrex::Math::abs(get<0>(hv) - 1) > 1.e-8)
        //    std::cout << "FW INITCONSTSP HERE "<< get<0>(hv) << std::endl;
        //    std::cout << std::endl;    
        //   amrex::Abort("Error:: Failed check of initial species summing to 1, C");
#else
        if (amrex::Math::abs(get<0>(hv) - 1) > 1.e-6)
            amrex::Abort("Error:: Failed check of initial species summing to 1, D");
#endif
    }
#endif
}

void
Nyx::init_e_from_T (const Real& a)
{

    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& D_new = get_new_data(DiagEOS_Type);
    Real h_species_in=h_species;
    Real gamma_minus_1_in=gamma - 1.0;
    auto atomic_rates = atomic_rates_glob;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        const auto s_arr    = S_new.array(mfi);
        const auto d_arr = D_new.array(mfi);

        amrex::Print() << "FW 2026's internal energy debug!" << std::endl;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
             Real rho = s_arr(i,j,k,Density_comp);
             Real T   = d_arr(i,j,k,Temp_comp);
             Real ne  = d_arr(i,j,k,  Ne_comp);

             Real eint;

             // Call EOS to get the internal energy
             Real dummy_pres=0.0;
             // Set temp to small_temp and compute corresponding internal energy
             nyx_eos_given_RT(atomic_rates, gamma_minus_1_in, h_species_in, &eint, &dummy_pres, rho, T, ne, a);

             s_arr(i,j,k,Eint_comp) = s_arr(i,j,k,Density_comp) * eint;

             s_arr(i,j,k,Eden_comp) = s_arr(i,j,k,Eint_comp) + 0.5 * ( 
                                 s_arr(i,j,k,Xmom_comp)*s_arr(i,j,k,Xmom_comp) +
                                 s_arr(i,j,k,Ymom_comp)*s_arr(i,j,k,Ymom_comp) +
                                 s_arr(i,j,k,Zmom_comp)*s_arr(i,j,k,Zmom_comp) ) / s_arr(i,j,k,Density_comp);
        });
        amrex::Gpu::Device::streamSynchronize();

        // ---- HOST DEBUG PRINT ----
        if (bx.contains(IntVect(0,0,0))) {
            const auto& s = S_new.const_array(mfi);
            const auto& d = D_new.const_array(mfi);
    
            amrex::Print()
                << "DEBUG (0,0,0): "
                << "rho="  << s(0,0,0,Density_comp)
                << " T="   << d(0,0,0,Temp_comp)
                << " ne="  << d(0,0,0,Ne_comp)
                << " Eint="<< s(0,0,0,Eint_comp)
                << " Eden="<< s(0,0,0,Eden_comp)
                << "\n";
        }
    }
}
#endif
