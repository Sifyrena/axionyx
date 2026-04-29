#include <iomanip>

#ifdef GRAVITY
#include <Gravity.H>
#endif

#include <Nyx.H>
#ifdef AMREX_PARTICLES
#include <NyxParticleContainer.H>
#endif
#include <constants_cosmo.H>

using namespace amrex;

// #ifdef GRAVITY
//  This function reads a multifab on disk.
void Nyx::icReadAndPrepareFab(std::string mfDirName, int nghost, MultiFab &mf)
{
  if (level > 0 && nghost > 0)
  {
    std::cout << "Are sure you want to do what you are doing?" << std::endl;
    amrex::Abort();
  }

  //
  // Read from the directory into mf
  //
  MultiFab mf_read;

  if (!mfDirName.empty() && mfDirName[mfDirName.length() - 1] != '/')
    mfDirName += '/';
  std::string Level = amrex::Concatenate("Level_", level, 1);
  mfDirName.append(Level);
  mfDirName.append("/Cell");

  VisMF::Read(mf_read, mfDirName.c_str());

  if (ParallelDescriptor::IOProcessor())
    std::cout << "mf read" << '\n';

  if (mf_read.contains_nan())
  {
    for (int i = 0; i < mf_read.nComp(); i++)
    {
      if (mf_read.contains_nan(i, 1))
      {
        std::cout << "Found NaNs in read_mf in component " << i << ". Going ahead. " << std::endl;
       // amrex::Abort("Nyx::init_particles: Your initial conditions contain NaNs!");
      }
    }
  }

  const auto &ba = parent->boxArray(level);
  const auto &dm = parent->DistributionMap(level);
  const auto &ba_read = mf_read.boxArray();
  int nc = mf_read.nComp();

  // if we don't use a cic scheme for the initial conditions,
  // we can safely set the number of ghost cells to 0
  // for multilevel ICs we can't use ghostcells
  mf.define(ba, dm, nc, nghost);

  mf.MultiFab::ParallelCopy(mf_read, 0, 0, nc, 0, 0);

  if (!((ba.contains(ba_read) && ba_read.contains(ba))))
  {
    if (ParallelDescriptor::IOProcessor())
    {
      std::cout << "ba      :" << ba << std::endl;
      std::cout << "ba_read :" << ba_read << std::endl;
      std::cout << "Read mf and hydro mf DO NOT cover the same domain!"
                << std::endl;
    }
    ParallelDescriptor::Barrier();
    if (ParallelDescriptor::IOProcessor())
    {
      amrex::Abort();
    }
  }

  mf_read.clear();

  mf.FillBoundary();
  mf.EnforcePeriodicity(geom.periodicity());

  // FIXME
  // mf.setVal(0);

  if (mf.contains_nan())
  {
    for (int i = 0; i < mf.nComp(); i++)
    {
      if (mf.contains_nan(i, 1, nghost))
      {
        std::cout << "Found NaNs in component " << i << ". Initialized! " << std::endl;
        //amrex::Abort("Nyx::init_particles: Your initial conditions contain NaNs!");
      }
    }
  }
}
// #endif

// This is for if the code was compiled with Hydro
#ifndef NO_HYDRO
void Nyx::initcosmo()
{

  //     if(parent->useFixedUpToLevel()<level)
  //       return;
  if (ParallelDescriptor::IOProcessor())
    std::cout << "Calling InitCosmo for level " << level << std::endl;

#ifdef GRAVITY
  Real comoving_OmL;
  const Real len[AMREX_SPACEDIM] = {geom.ProbLength(0), geom.ProbLength(1), geom.ProbLength(2)};

  Real particleMass;
  std::string mfDirName;

  Real redshift = -1;
  Vector<int> n_part(AMREX_SPACEDIM);

  if (level > parent->useFixedUpToLevel())
  {
    std::cout << "You have more refinement than grids, there might be a problem with your refinement criterion..." << std::endl;

    MultiFab &S_new = get_level(level).get_new_data(State_Type);
    MultiFab &D_new = get_level(level).get_new_data(DiagEOS_Type);

    FillCoarsePatch(S_new, 0, 0, State_Type, 0, S_new.nComp());
    FillCoarsePatch(D_new, 0, 0, DiagEOS_Type, 0, D_new.nComp());
    return;
  }

  //
  // Read the init directory name and particle mass from the inputs file
  //
  ParmParse pp("cosmo");
  pp.get("initDirName", mfDirName);
#ifdef FDM
  pp.query("mixed_cosmology", mixed_cosmology);
  std::string FDMmfDirName;
  if (mixed_cosmology)
  {
    if (ParallelDescriptor::IOProcessor())
      std::cout << "Using mixed cosmology" << std::endl;
    pp.get("FDMinitDirName", FDMmfDirName);
  }
#endif
#ifdef NUFLUID
  std::string nuMfDirName;
  pp.get("nuInitDirName", nuMfDirName);
#endif
  ParmParse pp2("nyx");
  pp2.get("initial_z", redshift);
  pp2.getarr("n_particles", n_part, 0, AMREX_SPACEDIM);

#ifdef NUFLUID
  Real comoving_OmNu;
  fort_get_omnu(&comoving_OmNu);
#endif

  // We now define this here instead of reading it.
  comoving_OmL = 1. - comoving_OmM;

  // //For whatever reason OmB is divided by OmM, which I need to undo here...
  // Real realOmB = comoving_OmB*comoving_OmM;

  // compute \Omega_K - just for checking flatness...
  Real comoving_OmK = 1 - comoving_OmM - comoving_OmL;
  if (std::abs(comoving_OmK) > 1e-3)
  {
    if (ParallelDescriptor::IOProcessor())
    {
      std::cout << "You assume a non-flat universe - \\Omega_K = "
                << comoving_OmK << std::endl;
    }
    amrex::Abort();
  }

  // compute \rho_{baryon}=3H_0^2\Omega_{baryon}/(8\pi G)
  Real rhoB = 3 * comoving_h * 100 * comoving_h * 100 * comoving_OmB / (8 * M_PI * Gconst);
  if (ParallelDescriptor::IOProcessor())
  {
    std::cout << "Mean baryonic matter density is " << rhoB
              << " M_sun/Mpc^3." << '\n';
  }

#ifdef NUFLUID
  // compute \rho_{\nu}=3H_0^2\Omega_{\nu}/(8\pi G)
  Real rhoNu = 3 * comoving_h * 100 * comoving_h * 100 * comoving_OmNu / (8 * M_PI * Gconst);
  if (ParallelDescriptor::IOProcessor())
  {
    std::cout << "Mean neutrino matter density is " << rhoNu
              << " M_sun/Mpc^3." << std::endl;
  }
#endif

  // compute \rho_{dark matter}=3H_0^2\Omega_{dark matter}/(8\pi G)
#ifdef NUFLUID
  Real comoving_OmD = comoving_OmM - comoving_OmB - comoving_OmNu;
#else
  Real comoving_OmD = comoving_OmM - comoving_OmB;
#endif
  Real rhoD = 3 * comoving_h * 100 * comoving_h * 100 * comoving_OmD / (8 * M_PI * Gconst);
  if (ParallelDescriptor::IOProcessor())
  {
    std::cout << "Mean dark matter density is " << rhoD
              << " M_sun/Mpc^3." << '\n';
  }

  if (!do_hydro && comoving_OmB > 0)
  {
    if (ParallelDescriptor::IOProcessor())
    {
      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << "You chose a baryonic matter content greater than 0, "
                << "but chose not to do hydro" << std::endl;
      std::cout << "I will be using \\rho_M for the dark matter density..." << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;
    }
    rhoD += rhoB;
  }

  // Reads the mf and checks the data...
  MultiFab mf;
  icReadAndPrepareFab(mfDirName, 0, mf);
#ifdef FDM
  MultiFab FDMmf;
  if (mixed_cosmology)
    icReadAndPrepareFab(FDMmfDirName, 0, FDMmf);
#endif
#ifdef NUFLUID
  MultiFab nuMf;
  icReadAndPrepareFab(nuMfDirName, 0, nuMf);
#endif

  // we have to calculate the initial a on our own
  // as there is no code path to get_comoving_a
  Real comoving_a = 1 / (1 + redshift);

  std::string icSource;
  pp.get("ic-source", icSource);
  int baryon_den, baryon_vx;
  int part_dx, part_vx;
  Real vel_fac[AMREX_SPACEDIM], dis_fac[AMREX_SPACEDIM];
  const Real *dx = geom.CellSize();

#ifdef NUFLUID
  int nu_den, nu_vx, nu_vy, nu_vz;
#endif

  if (icSource == "MUSIC")
  {
    if (do_hydro)
    {
      baryon_den = 0;
      baryon_vx = 1;
      part_dx = 4;
      part_vx = 7;
#ifdef NUFLUID
      nu_den = -1;
      nu_vx = -1;
#endif
    }
    else
    {
      baryon_den = -1;
      baryon_vx = -1;
      part_dx = 0;
      part_vx = 3;
#ifdef NUFLUID
      nu_den = -1;
      nu_vx = -1;
#endif

#ifdef FDM
      baryon_den = 0;
      baryon_vx = 1;
      part_dx = 4;
      part_vx = 7;
#endif
    }
    for (int n = 0; n < AMREX_SPACEDIM; n++)
    {
      // vel_fac[n] = comoving_a * len[n]/comoving_h;
      // vel_fac[n] = len[n]/comoving_h;
      vel_fac[n] = len[n] * comoving_h; // we need the boxlength in Mpc/h
      amrex::Print()<<"Velocity Prefactor is "<<vel_fac[n]<<std::endl;
      dis_fac[n] = len[n];
    }
    particleMass = rhoD * dx[0] * dx[1] * dx[2];
#ifdef FDM
    if (mixed_cosmology)
      particleMass *= (1.0 - ratio_fdm);
#endif
    if (ParallelDescriptor::IOProcessor())
    {
      std::cout << "Particle mass at level " << level
                << " is " << particleMass
                << " M_sun." << '\n';
    }
  }
  else if (icSource == "nyx")
  {
    baryon_den = 3;
    baryon_vx = 0;
    part_dx = 0;
    part_vx = 0;
#ifdef NUFLUID
    nu_den = 3;
    nu_vx = 0;
#endif
    // Velocities are proportional to displacements by
    // LBox [MPc] * a * H [km/s/MPc]
    for (int n = 0; n < AMREX_SPACEDIM; n++)
    {
      vel_fac[n] = len[n] * comoving_a * std::sqrt(comoving_OmM / pow(comoving_a, 3) + comoving_OmL) * comoving_h * 100;
      dis_fac[n] = len[n];
    }
    // Compute particle mass
    Real simulationVolume = len[0] * len[1] * len[2];
    Real numberOfParticles = n_part[0] * n_part[1] * n_part[2];
    particleMass = rhoD * simulationVolume / numberOfParticles;
#ifdef FDM
    if (mixed_cosmology)
      particleMass *= (1.0 - ratio_fdm);
#endif
    if (ParallelDescriptor::IOProcessor())
    {
      std::cout << "Particle mass is " << particleMass
                << " M_sun." << '\n';
    }
  }
  else
  {
    std::cout << "No clue from which code the initial conditions originate..." << std::endl
              << "Aborting...!" << std::endl;
    amrex::Abort();
  }

  BoxArray baWhereNot;
  if (level < parent->initialBaLevels())
  {
    //       baWhereNot = parent->boxArray(level+1);
    //       baWhereNot.coarsen(parent->refRatio()[level]);
    baWhereNot = parent->initialBa(level + 1);
    // std::cout << "Don't generate particles in region: " << baWhereNot << std::endl;
  }
  // initialize particles
  // std::cout << "FIXME!!!!!!!!!!!!!!!!!!" << std::endl;
  BoxArray myBaWhereNot(baWhereNot.size());
  for (int i = 0; i < baWhereNot.size(); i++)
    myBaWhereNot.set(i, baWhereNot[i]);
  if (level < parent->initialBaLevels())
    myBaWhereNot.coarsen(parent->refRatio(level));

  // we also have to restrict the creation of particles to non refined parts of the domain.
  // if (level == 0)
  //   Nyx::theDMPC()->InitCosmo1ppcMultiLevel(mf, dis_fac, vel_fac,
  // 					      particleMass,
  // 					      part_dx, part_vx,
  // 					      myBaWhereNot,
  // 					      level, parent->initialBaLevels()+1);
#ifdef AMREX_PARTICLES
  if (do_dm_particles)
    Nyx::theDMPC()->InitCosmo1ppcMultiLevel(mf, dis_fac, vel_fac,
                                            particleMass,
                                            part_dx, part_vx,
                                            myBaWhereNot,
                                            level, parent->initialBaLevels() + 1);
#endif

#ifdef FDM

  /*Copy density & velocity from baryon to axion state*/
  MultiFab &Ax_new = get_level(level).get_new_data(Axion_Type);
  Ax_new.setVal(0.);
  if (mixed_cosmology)
    Ax_new.ParallelCopy(FDMmf, baryon_den, 0, 4, 0, Ax_new.nGrow(), parent->Geom(level).periodicity(), FabArrayBase::COPY);
  else
    Ax_new.ParallelCopy(mf, baryon_den, 0, 4, 0, Ax_new.nGrow(), parent->Geom(level).periodicity(), FabArrayBase::COPY);
  Ax_new.plus(1, 0, 1, Ax_new.nGrow());
  Real meanfdmdens;
  if (mixed_cosmology)
    meanfdmdens = ratio_fdm * rhoD;
  else
    meanfdmdens = rhoD;
  Ax_new.mult(meanfdmdens, 0, 1, Ax_new.nGrow());
  Ax_new.mult(vel_fac[0], 1, 1, Ax_new.nGrow());
  Ax_new.mult(vel_fac[1], 2, 1, Ax_new.nGrow());
  Ax_new.mult(vel_fac[2], 3, 1, Ax_new.nGrow());
  /*Once velocity is defined on all levels, construct phase from velocity divergence*/
  if (level == (parent->initialBaLevels()))
  {

    Vector<std::unique_ptr<MultiFab>> div(parent->initialBaLevels() + 1);
    Vector<MultiFab *> phase(parent->initialBaLevels() + 1);
    Vector<Vector<std::unique_ptr<MultiFab>>> grad_phase(parent->initialBaLevels() + 1);
    Vector<std::array<MultiFab *, AMREX_SPACEDIM>> grad_phase_aa;

    for (int lev = 0; lev <= parent->initialBaLevels(); lev++)
    {

      const BoxArray &ba = parent->boxArray()[lev];
      const DistributionMapping &dmap = parent->DistributionMap()[lev];
      div[lev].reset(new MultiFab(ba, dmap, 1, 0));
      div[lev]->setVal(0.0);
      phase[lev] = new MultiFab(ba, dmap, 1, 1);
      phase[lev]->setVal(0.0);

      MultiFab &Axvel = get_level(lev).get_new_data(Axion_Type);
      // for (MFIter mfi(Axvel,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      for (FillPatchIterator
               fpi(*this, Axvel, 1, state[Axion_Type].curTime(), Axion_Type, 0, Nyx::NUM_AX);
           fpi.isValid(); ++fpi)
      {
        const Box &bx = fpi.tilebox();
        Array4<Real> const &velarr = fpi().array(); // Axvel.array(fpi);
        Array4<Real> const &divarr = (*(div[lev]))[fpi].array();
        AMREX_PARALLEL_FOR_3D(bx, i, j, k,
                              {
                                divarr(i, j, k, 0) = (velarr(i + 1, j, k, 1) - velarr(i - 1, j, k, 1)) / dx[0] * comoving_a / 2.0 + (velarr(i, j + 1, k, 2) - velarr(i, j - 1, k, 2)) / dx[1] * comoving_a / 2.0 + (velarr(i, j, k + 1, 3) - velarr(i, j, k - 1, 3)) / dx[2] * comoving_a / 2.0;
                              });
      }
      grad_phase[lev].resize(BL_SPACEDIM);
      for (int n = 0; n < BL_SPACEDIM; ++n)
        grad_phase[lev][n].reset(new MultiFab(BoxArray(ba).surroundingNodes(n), dmap, 1, 1));
      grad_phase_aa.push_back({AMREX_D_DECL(GetVecOfVecOfPtrs(grad_phase)[lev][0],
                                            GetVecOfVecOfPtrs(grad_phase)[lev][1],
                                            GetVecOfVecOfPtrs(grad_phase)[lev][2])});
    }

    const MultiFab *crse_bcdata = nullptr;
    Real rel_eps = 1.e-10;
    Real abs_eps = 0.;
    gravity->solve_with_MLMG(0, parent->initialBaLevels(), phase, amrex::GetVecOfConstPtrs(div),
                             grad_phase_aa, crse_bcdata, rel_eps, abs_eps);
    /*Fill FDM particle container and axion state with initial information*/
    for (int lev = 0; lev <= parent->initialBaLevels(); lev++)
    {

      BoxArray baWhereNotfdm;
      if (lev < parent->initialBaLevels())
        baWhereNotfdm = parent->initialBa(lev + 1);
      BoxArray myBaWhereNotfdm(baWhereNotfdm.size());
      for (int i = 0; i < baWhereNotfdm.size(); i++)
        myBaWhereNotfdm.set(i, baWhereNotfdm[i]);
      if (lev < parent->initialBaLevels())
        myBaWhereNotfdm.coarsen(parent->refRatio(lev));

      MultiFab &Ax_new = get_level(lev).get_new_data(Axion_Type);

      Ax_new.ParallelCopy((*phase[lev]), 0, Nyx::AxPhas, 1, 0, Ax_new.nGrow(), parent->Geom(lev).periodicity(), FabArrayBase::COPY);
      const Real moverhbar = 1.0 / hbaroverm;
      Ax_new.mult(moverhbar, Nyx::AxPhas, 1, Ax_new.nGrow());

      for (MFIter mfi(Ax_new, TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
        const Box &bx = mfi.tilebox();
        Array4<Real> const &fab = Ax_new.array(mfi);
        AMREX_PARALLEL_FOR_3D(bx, i, j, k,
                              {
                                fab(i, j, k, Nyx::AxRe) = std::sqrt(fab(i, j, k, Nyx::AxDens)) * std::cos(fab(i, j, k, Nyx::AxPhas));
                                fab(i, j, k, Nyx::AxIm) = std::sqrt(fab(i, j, k, Nyx::AxDens)) * std::sin(fab(i, j, k, Nyx::AxPhas));
                              });
      }
      Ax_new.FillBoundary(geom.periodicity());
    }
  }
#endif

  MultiFab &S_new = get_level(level).get_new_data(State_Type);
  MultiFab &D_new = get_level(level).get_new_data(DiagEOS_Type);

#ifdef NUFLUID
  // copy density
  S_new.MultiFab::ParallelCopy(nuMf, nu_den, NuDens, 1, 0, 0);
  S_new.plus(1, NuDens, 1, S_new.nGrow());
  S_new.mult(rhoNu, NuDens, 1, S_new.nGrow());

  // copy velocities...
  S_new.MultiFab::ParallelCopy(nuMf, nu_vx, NuXMom, 3, 0, 0);

  //...and "transform" to momentum
  S_new.mult(vel_fac[0], NuXMom, 1, S_new.nGrow());
  S_new.mult(vel_fac[1], NuYMom, 1, S_new.nGrow());
  S_new.mult(vel_fac[2], NuZMom, 1, S_new.nGrow());
  MultiFab::Multiply(S_new, S_new, NuDens, NuXMom, 1, S_new.nGrow());
  MultiFab::Multiply(S_new, S_new, NuDens, NuYMom, 1, S_new.nGrow());
  MultiFab::Multiply(S_new, S_new, NuDens, NuZMom, 1, S_new.nGrow());
#endif

  // initialize hydro
  // units for velocity should be okay
  if (do_hydro)
  {
    if (ParallelDescriptor::IOProcessor())
    {
      std::cout << "Performing hydro init..." << '\n';
    }

    // Fill everything with old data...should only affect ghostzones, but
    // seems to have no effect...
    if (level > 0)
    {
      FillCoarsePatch(S_new, 0, 0, State_Type, 0, S_new.nComp());
      FillCoarsePatch(D_new, 0, 0, DiagEOS_Type, 0, D_new.nComp());
    }

    // copy density
    S_new.setVal(0.);
    S_new.MultiFab::ParallelCopy(mf, baryon_den, Density_comp, 1, 0, 0);
    S_new.plus(1, Density_comp, 1, S_new.nGrow());
    S_new.mult(rhoB, Density_comp, 1, S_new.nGrow());

    //      //This block assigns "the same" density for the baryons as for the dm.
    //      Vector<std::unique_ptr<MultiFab> > particle_mf;
    //      Nyx::theDMPC()->AssignDensity(particle_mf);
    //      particle_mf[0]->mult(comoving_OmB / comoving_OmD);
    //      S_new.copy(*particle_mf[0], 0, Density_comp, 1);

    // copy velocities...
    S_new.MultiFab::ParallelCopy(mf, baryon_vx, Xmom_comp, 3, 0, 0);

    //...and "transform" to momentum
    S_new.mult(vel_fac[0], Xmom_comp, 1, S_new.nGrow());
    MultiFab::Multiply(S_new, S_new, Density_comp, Xmom_comp, 1, S_new.nGrow());
    S_new.mult(vel_fac[1], Ymom_comp, 1, S_new.nGrow());
    MultiFab::Multiply(S_new, S_new, Density_comp, Ymom_comp, 1, S_new.nGrow());
    S_new.mult(vel_fac[2], Zmom_comp, 1, S_new.nGrow());
    MultiFab::Multiply(S_new, S_new, Density_comp, Zmom_comp, 1, S_new.nGrow());

    Real tempInit = 0.021 * (1.0 + redshift) * (1.0 + redshift);

    D_new.setVal(tempInit, Temp_comp);
    D_new.setVal(0.0, Ne_comp);
    if (inhomo_reion > 0)
      D_new.setVal(0.0, Zhi_comp);

    init_e_from_T(old_a);

#ifndef CONST_SPECIES
    // Convert X_i to (rho X)_i
    {
      S_new.setVal(0.75, FirstSpec_comp);
      S_new.setVal(0.25, FirstSpec_comp + 1);

      for (int i = 0; i < NumSpec; i++)
      {
        MultiFab::Multiply(S_new, S_new, Density_comp, FirstSpec_comp + i, 1, 0);
      }
    }
#endif
  }
  else
  {
    S_new.setVal(0.0, Density_comp);
    D_new.setVal(-23, Temp_comp);
    D_new.setVal(-42, Ne_comp);
  }

  mf.clear();

#endif // ifdef GRAVITY
}
#else
#ifdef FDM // FDM, NoHydro
void Nyx::initcosmo()
{
  if (ParallelDescriptor::IOProcessor())
    std::cout << "No Hydro! Calling Frank Wang's New InitCosmo for level " << level << std::endl;

#ifdef GRAVITY
  Real comoving_OmL;
  const Real len[BL_SPACEDIM] = {geom.ProbLength(0), geom.ProbLength(1), geom.ProbLength(2)};

  Real particleMass;
  std::string mfDirName;

  Real redshift = -1;
  bool clear_phase = false;
  Vector<int> n_part(BL_SPACEDIM);

  if (level > parent->useFixedUpToLevel())
  {
    std::cout << "You have more refinement than grids, there might be a problem with your refinement criterion..." << std::endl;
    return;
  }

  ParmParse pp("cosmo");
  pp.get("initDirName", mfDirName);
#ifdef FDM
  pp.query("mixed_cosmology", mixed_cosmology);
  std::string FDMmfDirName;
  if (mixed_cosmology)
  {
    if (ParallelDescriptor::IOProcessor())
      std::cout << "Using mixed dark matter (MDM) cosmology" << std::endl;
    pp.get("FDMinitDirName", FDMmfDirName);
  }
#endif
  ParmParse pp2("nyx");
  pp2.get("initial_z", redshift);
  pp2.getarr("n_particles", n_part, 0, BL_SPACEDIM);
  pp2.query("clear_phase", clear_phase);
  amrex::Print() << "Initial Axion Phase Clear?" << clear_phase << std::endl;

  comoving_OmL = 1. - comoving_OmM;

  Real rhoD = 3 * comoving_h * 100 * comoving_h * 100 * comoving_OmM / (8 * M_PI * Gconst);
  if (ParallelDescriptor::IOProcessor())
  {
    std::cout << "Mean dark matter density is " << rhoD
              << " M_sun/Mpc^3." << '\n';
  }

  MultiFab mf;
  icReadAndPrepareFab(mfDirName, 0, mf);
#ifdef FDM
  MultiFab FDMmf;
  if (mixed_cosmology)
    icReadAndPrepareFab(FDMmfDirName, 0, FDMmf);
#endif
  Real comoving_a = 1 / (1 + redshift);

  std::string icSource;
  pp.get("ic-source", icSource);
  int baryon_den, baryon_vx;
  int part_dx, part_vx;
  Real vel_fac[BL_SPACEDIM], dis_fac[BL_SPACEDIM];
  const Real *dx = geom.CellSize();

  if (icSource == "MUSIC")
  {
    baryon_den = 0;
    baryon_vx = 1;
    part_dx = 4;
    part_vx = 7;

    for (int n = 0; n < BL_SPACEDIM; n++)
    {
      amrex::Print()<<"The Length Thing is "<< len[n] <<std::endl;
      vel_fac[n] = len[n] * comoving_h;
      amrex::Print()<<"The Velocity Transformation Factor is " << vel_fac[n] << std::endl;
      dis_fac[n] = len[n];
    }
    particleMass = rhoD * dx[0] * dx[1] * dx[2];
#ifdef FDM
    if (mixed_cosmology)
      particleMass *= (1.0 - ratio_fdm);
#endif
    if (ParallelDescriptor::IOProcessor())
    {
      std::cout << "The particle mass at level " << level
                << " is " << particleMass
                << " M_sun." << '\n';
    }
  }
  else if (icSource == "nyx")
  {
    baryon_den = 3;
    baryon_vx = 0;
    part_dx = 0;
    part_vx = 0;

    for (int n = 0; n < BL_SPACEDIM; n++)
    {
      vel_fac[n] = len[n] * comoving_a * std::sqrt(comoving_OmM / pow(comoving_a, 3) + comoving_OmL) * comoving_h * 100;
      dis_fac[n] = len[n];
    }
    Real simulationVolume = len[0] * len[1] * len[2];
    Real numberOfParticles = n_part[0] * n_part[1] * n_part[2];

    amrex::Print() << "There are " << numberOfParticles << " particles here!\n";
    particleMass = rhoD * simulationVolume / numberOfParticles;
#ifdef FDM
    if (mixed_cosmology)
      particleMass *= (1.0 - ratio_fdm);
      amrex::Print()<<"Using Mixed Cosmo, Scaling Down N Body Particle Mass\n";
#endif
    if (ParallelDescriptor::IOProcessor())
    {
      std::cout << "Particle mass is " << particleMass
                << " M_sun." << '\n';
    }
  }
  else
  {
    std::cout << "We have no clue from which code the initial coniditions originate..." << std::endl
              << "Aborting...!" << std::endl;
    amrex::Abort();
  }

  BoxArray baWhereNot;
  if (level < parent->initialBaLevels())
    baWhereNot = parent->initialBa(level + 1);
  BoxArray myBaWhereNot(baWhereNot.size());
  for (int i = 0; i < baWhereNot.size(); i++)
    myBaWhereNot.set(i, baWhereNot[i]);
  if (level < parent->initialBaLevels())
    myBaWhereNot.coarsen(parent->refRatio(level));

#ifdef AMREX_PARTICLES
  if (do_dm_particles)
  {
    Nyx::theDMPC()->InitCosmo1ppcMultiLevel(mf, dis_fac, vel_fac,
                                            particleMass,
                                            part_dx, part_vx,
                                            myBaWhereNot,
                                            level, parent->initialBaLevels() + 1);
  }
#endif
#ifdef FDM

  /*Copy density & velocity from baryon to axion state*/
  MultiFab &Ax_new = get_level(level).get_new_data(Axion_Type);
  Ax_new.setVal(0.);
  if (mixed_cosmology)
    Ax_new.ParallelCopy(FDMmf, baryon_den, 0, 4, 0, Ax_new.nGrow(), parent->Geom(level).periodicity(), FabArrayBase::COPY);
  else
    Ax_new.ParallelCopy(mf, baryon_den, 0, 4, 0, Ax_new.nGrow(), parent->Geom(level).periodicity(), FabArrayBase::COPY);
  Ax_new.plus(1, 0, 1, Ax_new.nGrow());
  Real meanfdmdens;
  if (mixed_cosmology)
    meanfdmdens = ratio_fdm * rhoD;
  else
    meanfdmdens = rhoD;
  Ax_new.mult(meanfdmdens, 0, 1, Ax_new.nGrow());
  Ax_new.mult(vel_fac[0], 1, 1, Ax_new.nGrow());
  Ax_new.mult(vel_fac[1], 2, 1, Ax_new.nGrow());
  Ax_new.mult(vel_fac[2], 3, 1, Ax_new.nGrow());

  /*Once velocity is defined on all levels, construct phase from velocity divergence*/
  if (level == (parent->initialBaLevels()))
  {

    Vector<std::unique_ptr<MultiFab>> div(parent->initialBaLevels() + 1);
    Vector<MultiFab *> phase(parent->initialBaLevels() + 1);
    Vector<Vector<std::unique_ptr<MultiFab>>> grad_phase(parent->initialBaLevels() + 1);
    Vector<std::array<MultiFab *, AMREX_SPACEDIM>> grad_phase_aa;

    for (int lev = 0; lev <= parent->initialBaLevels(); lev++)
    {

      const BoxArray &ba = parent->boxArray()[lev];
      const DistributionMapping &dmap = parent->DistributionMap()[lev];
      div[lev].reset(new MultiFab(ba, dmap, 1, 0));
      div[lev]->setVal(0.0);
      phase[lev] = new MultiFab(ba, dmap, 1, 1);
      phase[lev]->setVal(0.0);

      MultiFab &Axvel = get_level(lev).get_new_data(Axion_Type);

      /* This is still experimental
      MultiFab cdm_density(ba, dmap, 1, 1); // 1 component, 1 ghost
      cdm_density.setVal(0.0);

      NyxParticleContainer::ParticleLevel &plev = DarkMatterContainer->GetParticles(lev);
      DarkMatterContainer->AssignDensity(cdm_density, lev, 0);
      */

      // for (MFIter mfi(Axvel,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      for (FillPatchIterator
               fpi(*this, Axvel, 1, state[Axion_Type].curTime(), Axion_Type, 0, Nyx::NUM_AX);
           fpi.isValid(); ++fpi)
      {
        const Box &bx = fpi.tilebox();
        Array4<Real> const &velarr = fpi().array(); // Axvel.array(fpi);
        Array4<Real> const &divarr = (*(div[lev]))[fpi].array();
        AMREX_PARALLEL_FOR_3D(bx, i, j, k,
                              {
                                divarr(i, j, k, 0) = (velarr(i + 1, j, k, 1) - velarr(i - 1, j, k, 1)) / dx[0] * comoving_a / 2.0 + (velarr(i, j + 1, k, 2) - velarr(i, j - 1, k, 2)) / dx[1] * comoving_a / 2.0 + (velarr(i, j, k + 1, 3) - velarr(i, j, k - 1, 3)) / dx[2] * comoving_a / 2.0;
                              });
      }

      grad_phase[lev].resize(BL_SPACEDIM);
      for (int n = 0; n < BL_SPACEDIM; ++n)
        grad_phase[lev][n].reset(new MultiFab(BoxArray(ba).surroundingNodes(n), dmap, 1, 1));
      grad_phase_aa.push_back({AMREX_D_DECL(GetVecOfVecOfPtrs(grad_phase)[lev][0],
                                            GetVecOfVecOfPtrs(grad_phase)[lev][1],
                                            GetVecOfVecOfPtrs(grad_phase)[lev][2])});
    }
    const MultiFab *crse_bcdata = nullptr;
    Real rel_eps = 1.e-10;
    Real abs_eps = 0.;
    gravity->solve_with_MLMG(0, parent->initialBaLevels(), phase, amrex::GetVecOfConstPtrs(div),
                             grad_phase_aa, crse_bcdata, rel_eps, abs_eps);

    /*Fill FDM particle container and axion state with initial information*/
    for (int lev = 0; lev <= parent->initialBaLevels(); lev++)
    {

      BoxArray baWhereNotfdm;
      if (lev < parent->initialBaLevels())
        baWhereNotfdm = parent->initialBa(lev + 1);
      BoxArray myBaWhereNotfdm(baWhereNotfdm.size());
      for (int i = 0; i < baWhereNotfdm.size(); i++)
        myBaWhereNotfdm.set(i, baWhereNotfdm[i]);
      if (lev < parent->initialBaLevels())
        myBaWhereNotfdm.coarsen(parent->refRatio(lev));

      MultiFab &Ax_new = get_level(lev).get_new_data(Axion_Type);

      Ax_new.ParallelCopy((*phase[lev]), 0, Nyx::AxPhas, 1, 0, Ax_new.nGrow(), parent->Geom(lev).periodicity(), FabArrayBase::COPY);
      if (clear_phase)
      {
        amrex::Print() << "FW Debug: Setting Axion Phase to Zero\n";
        for (MFIter mfi(Ax_new, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          const Box &bx = mfi.tilebox();
          Array4<Real> const &fab = Ax_new.array(mfi);
          AMREX_PARALLEL_FOR_3D(bx, i, j, k,
                                {
                                  fab(i, j, k, Nyx::AxRe) = std::sqrt(fab(i, j, k, Nyx::AxDens));
                                  fab(i, j, k, Nyx::AxIm) = 0;
                                  fab(i, j, k, Nyx::AxPhas) = 0;
                                });
        }
      }
      else
      {
        amrex::Print() << "Using MUSIC GENERATED AXION PHASE\n";
        for (MFIter mfi(Ax_new, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          const Box &bx = mfi.tilebox();
          Array4<Real> const &fab = Ax_new.array(mfi);
          AMREX_PARALLEL_FOR_3D(bx, i, j, k,
                                {
                                  fab(i, j, k, Nyx::AxRe) = std::sqrt(fab(i, j, k, Nyx::AxDens)) * std::cos(fab(i, j, k, Nyx::AxPhas) / hbaroverm);
                                  fab(i, j, k, Nyx::AxIm) = std::sqrt(fab(i, j, k, Nyx::AxDens)) * std::sin(fab(i, j, k, Nyx::AxPhas) / hbaroverm);
                                });
        }
      }
      Ax_new.FillBoundary(geom.periodicity());
    }
  }
#endif

#ifdef INFGRAV
  MultiFab &TRe = get_level(level).get_new_data(TensorRe_Type);
  MultiFab &TIm = get_level(level).get_new_data(TensorIm_Type);
  MultiFab &IcosRe = get_level(level).get_new_data(CosRe_Type);
  MultiFab &IcosIm = get_level(level).get_new_data(CosIm_Type);
  MultiFab &IsinRe = get_level(level).get_new_data(SinRe_Type);
  MultiFab &IsinIm = get_level(level).get_new_data(SinIm_Type);

  TRe.setVal(0.0); // at the beginning/initialization set tensor components to zero
  TIm.setVal(0.0);
  IcosRe.setVal(0.0);
  IcosIm.setVal(0.0);
  IsinRe.setVal(0.0);
  IsinIm.setVal(0.0);

  fill_tensor_TT();
  gw_prevtime = 0.0;

#endif
  mf.clear();
#endif // ifdef GRAVITY
}

void Nyx::initFDMcosmo()
{
  
  amrex::Print() << "No Hydro! Calling the New InitCosmo for Axions Only at level " << level << std::endl;

#ifdef GRAVITY
  Real comoving_OmL;
  const Real len[BL_SPACEDIM] = {geom.ProbLength(0), geom.ProbLength(1), geom.ProbLength(2)};

  Real particleMass;
  std::string mfDirName;

  Real redshift = -1;
  bool clear_phase = false;
  Vector<int> n_part(BL_SPACEDIM);

  if (level > parent->useFixedUpToLevel())
  {
    std::cout << "You have more refinement than grids, there might be a problem with your refinement criterion..." << std::endl;
    return;
  }

  ParmParse pp("cosmo");
  pp.get("initDirName", mfDirName);
#ifdef FDM
  pp.query("mixed_cosmology", mixed_cosmology);
  std::string FDMmfDirName;
  if (mixed_cosmology)
  {
    if (ParallelDescriptor::IOProcessor())
      std::cout << "Using mixed cosmology" << std::endl;
    pp.get("FDMinitDirName", FDMmfDirName);
  }
#endif
  ParmParse pp2("nyx");
  pp2.get("initial_z", redshift);
  pp2.getarr("n_particles", n_part, 0, BL_SPACEDIM);
  pp2.query("clear_phase", clear_phase);
  amrex::Print() << "Initial Axion Phase Clear?" << clear_phase << std::endl;

  comoving_OmL = 1. - comoving_OmM;

  Real rhoD = 3 * comoving_h * 100 * comoving_h * 100 * comoving_OmM / (8 * M_PI * Gconst);
  amrex::Print() << "Mean dark matter density is " << rhoD
              << " M_sun/Mpc^3." << '\n';

  MultiFab mf;
  icReadAndPrepareFab(mfDirName, 0, mf);
#ifdef FDM
  MultiFab FDMmf;
  if (mixed_cosmology)
    icReadAndPrepareFab(FDMmfDirName, 0, FDMmf);
#endif
  Real comoving_a = 1 / (1 + redshift);

  std::string icSource;
  pp.get("ic-source", icSource);
  int baryon_den, baryon_vx;
  int part_dx, part_vx;
  Real vel_fac[BL_SPACEDIM], dis_fac[BL_SPACEDIM];
  const Real *dx = geom.CellSize();

  if (icSource == "MUSIC")
  {
    baryon_den = 0;
    baryon_vx = 1;
    part_dx = 4;
    part_vx = 7;

    for (int n = 0; n < BL_SPACEDIM; n++)
    {
      amrex::Print()<<"The Length Thing is "<< len[n] <<std::endl;
      vel_fac[n] = len[n] * comoving_h;
      amrex::Print()<<"The Velocity Transformation Factor is " << vel_fac[n] << std::endl;
      dis_fac[n] = len[n];
    }
    particleMass = rhoD * dx[0] * dx[1] * dx[2];
#ifdef FDM
    if (mixed_cosmology)
      particleMass *= (1.0 - ratio_fdm);
#endif
    if (ParallelDescriptor::IOProcessor())
    {
      std::cout << "The particle mass at level " << level
                << " is " << particleMass
                << " M_sun." << '\n';
    }
  }
  else 
  {
    std::cout << "This function takes only MUSIC input." << std::endl
              << "Aborting...!" << std::endl;
    amrex::Abort();
  }

  BoxArray baWhereNot;
  if (level < parent->initialBaLevels())
    baWhereNot = parent->initialBa(level + 1);
  BoxArray myBaWhereNot(baWhereNot.size());
  for (int i = 0; i < baWhereNot.size(); i++)
    myBaWhereNot.set(i, baWhereNot[i]);
  if (level < parent->initialBaLevels())
    myBaWhereNot.coarsen(parent->refRatio(level));


#ifdef FDM

  /*Copy density & velocity from baryon to axion state*/
  MultiFab &Ax_new = get_level(level).get_new_data(Axion_Type);
  Ax_new.setVal(0.);
  if (mixed_cosmology)
    Ax_new.ParallelCopy(FDMmf, baryon_den, 0, 4, 0, Ax_new.nGrow(), parent->Geom(level).periodicity(), FabArrayBase::COPY);
  else
    Ax_new.ParallelCopy(mf, baryon_den, 0, 4, 0, Ax_new.nGrow(), parent->Geom(level).periodicity(), FabArrayBase::COPY);
  Ax_new.plus(1, 0, 1, Ax_new.nGrow());
  Real meanfdmdens;
  if (mixed_cosmology)
    meanfdmdens = ratio_fdm * rhoD;
  else
    meanfdmdens = rhoD;
  Ax_new.mult(meanfdmdens, 0, 1, Ax_new.nGrow());
  Ax_new.mult(vel_fac[0], 1, 1, Ax_new.nGrow());
  Ax_new.mult(vel_fac[1], 2, 1, Ax_new.nGrow());
  Ax_new.mult(vel_fac[2], 3, 1, Ax_new.nGrow());

  /*Once velocity is defined on all levels, construct phase from velocity divergence*/
  if (level == (parent->initialBaLevels()))
  {
    Vector<std::unique_ptr<MultiFab>> div(parent->initialBaLevels() + 1);
    Vector<MultiFab *> phase(parent->initialBaLevels() + 1);
    Vector<Vector<std::unique_ptr<MultiFab>>> grad_phase(parent->initialBaLevels() + 1);
    Vector<std::array<MultiFab *, AMREX_SPACEDIM>> grad_phase_aa;

    for (int lev = 0; lev <= parent->initialBaLevels(); lev++)
    {

      const BoxArray &ba = parent->boxArray()[lev];
      const DistributionMapping &dmap = parent->DistributionMap()[lev];
      div[lev].reset(new MultiFab(ba, dmap, 1, 0));
      div[lev]->setVal(0.0);
      phase[lev] = new MultiFab(ba, dmap, 1, 1);
      phase[lev]->setVal(0.0);

      MultiFab &Axvel = get_level(lev).get_new_data(Axion_Type);

      for (FillPatchIterator
               fpi(*this, Axvel, 1, state[Axion_Type].curTime(), Axion_Type, 0, Nyx::NUM_AX);
           fpi.isValid(); ++fpi)
      {
        const Box &bx = fpi.tilebox();
        Array4<Real> const &velarr = fpi().array(); // Axvel.array(fpi);
        Array4<Real> const &divarr = (*(div[lev]))[fpi].array();
        AMREX_PARALLEL_FOR_3D(bx, i, j, k,
                              {
                                divarr(i, j, k, 0) = (velarr(i + 1, j, k, 1) - velarr(i - 1, j, k, 1)) / dx[0] * comoving_a / 2.0 + (velarr(i, j + 1, k, 2) - velarr(i, j - 1, k, 2)) / dx[1] * comoving_a / 2.0 + (velarr(i, j, k + 1, 3) - velarr(i, j, k - 1, 3)) / dx[2] * comoving_a / 2.0;
                              });
      }

      grad_phase[lev].resize(BL_SPACEDIM);
      for (int n = 0; n < BL_SPACEDIM; ++n)
        grad_phase[lev][n].reset(new MultiFab(BoxArray(ba).surroundingNodes(n), dmap, 1, 1));
      grad_phase_aa.push_back({AMREX_D_DECL(GetVecOfVecOfPtrs(grad_phase)[lev][0],
                                            GetVecOfVecOfPtrs(grad_phase)[lev][1],
                                            GetVecOfVecOfPtrs(grad_phase)[lev][2])});
    }
    const MultiFab *crse_bcdata = nullptr;
    Real rel_eps = 1.e-10;
    Real abs_eps = 0.;
    gravity->solve_with_MLMG(0, parent->initialBaLevels(), phase, amrex::GetVecOfConstPtrs(div),
                             grad_phase_aa, crse_bcdata, rel_eps, abs_eps);

    /*Fill FDM particle container and axion state with initial information*/
    for (int lev = 0; lev <= parent->initialBaLevels(); lev++)
    {

      BoxArray baWhereNotfdm;
      if (lev < parent->initialBaLevels())
        baWhereNotfdm = parent->initialBa(lev + 1);
      BoxArray myBaWhereNotfdm(baWhereNotfdm.size());
      for (int i = 0; i < baWhereNotfdm.size(); i++)
        myBaWhereNotfdm.set(i, baWhereNotfdm[i]);
      if (lev < parent->initialBaLevels())
        myBaWhereNotfdm.coarsen(parent->refRatio(lev));

      MultiFab &Ax_new = get_level(lev).get_new_data(Axion_Type);

      Ax_new.ParallelCopy((*phase[lev]), 0, Nyx::AxPhas, 1, 0, Ax_new.nGrow(), parent->Geom(lev).periodicity(), FabArrayBase::COPY);
      if (clear_phase)
      {
        amrex::Print() << "FW Debug: Setting Axion Phase to Zero\n";
        for (MFIter mfi(Ax_new, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          const Box &bx = mfi.tilebox();
          Array4<Real> const &fab = Ax_new.array(mfi);
          AMREX_PARALLEL_FOR_3D(bx, i, j, k,
                                {
                                  fab(i, j, k, Nyx::AxRe) = std::sqrt(fab(i, j, k, Nyx::AxDens));
                                  fab(i, j, k, Nyx::AxIm) = 0;
                                  fab(i, j, k, Nyx::AxPhas) = 0;
                                });
        }
      }
      else
      {
        amrex::Print() << "Using MUSIC GENERATED AXION PHASE\n";
        for (MFIter mfi(Ax_new, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          const Box &bx = mfi.tilebox();
          Array4<Real> const &fab = Ax_new.array(mfi);
          AMREX_PARALLEL_FOR_3D(bx, i, j, k,
                                {
                                  fab(i, j, k, Nyx::AxRe) = std::sqrt(fab(i, j, k, Nyx::AxDens)) * std::cos(fab(i, j, k, Nyx::AxPhas) / hbaroverm);
                                  fab(i, j, k, Nyx::AxIm) = std::sqrt(fab(i, j, k, Nyx::AxDens)) * std::sin(fab(i, j, k, Nyx::AxPhas) / hbaroverm);
                                });
        }
      }
      Ax_new.FillBoundary(geom.periodicity());
    }
  }
#endif

#ifdef INFGRAV
  MultiFab &TRe = get_level(level).get_new_data(TensorRe_Type);
  MultiFab &TIm = get_level(level).get_new_data(TensorIm_Type);
  MultiFab &IcosRe = get_level(level).get_new_data(CosRe_Type);
  MultiFab &IcosIm = get_level(level).get_new_data(CosIm_Type);
  MultiFab &IsinRe = get_level(level).get_new_data(SinRe_Type);
  MultiFab &IsinIm = get_level(level).get_new_data(SinIm_Type);

  TRe.setVal(0.0); // at the beginning/initialization set tensor components to zero
  TIm.setVal(0.0);
  IcosRe.setVal(0.0);
  IcosIm.setVal(0.0);
  IsinRe.setVal(0.0);
  IsinIm.setVal(0.0);

  fill_tensor_TT();
  gw_prevtime = 0.0;

#endif
  mf.clear();
#endif // ifdef GRAVITY
}


#else
void Nyx::initcosmo()
{
  amrex::Abort("Nyx::initcosmo is not defined for NO_HYDRO=TRUE");
}
#endif

#endif
