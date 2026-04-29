#ifdef AXIONYX
#include <Nyx.H>
#ifdef GRAVITY
#include <Gravity.H>
#endif

using namespace amrex;

Real Nyx::advance_axionyx(Real time,
                          Real dt,
                          int iteration,
                          int ncycle)
{
  BL_PROFILE("Nyx::advance_fdm()");

  if (level == 0 || iteration > 1)
  {
    amrex::Gpu::LaunchSafeGuard lsg(true);
    const int finest_level = parent->finestLevel();
    for (int lev = level; lev <= finest_level; lev++)
    {
      Real dt_lev = parent->dtLevel(lev);
      for (int k = 0; k < NUM_STATE_TYPE; k++)
      {
        get_level(lev).state[k].allocOldData();
        get_level(lev).state[k].swapTimeLevels(dt_lev);
      }
    }
  }

#ifdef GRAVITY
  if (do_grav)
  {

    const int finest_level = parent->finestLevel();
    BL_PROFILE_VAR("solve_for_old_phi", solve_for_old_phi);
    if (level == 0 || iteration > 1)
    {
      MultiFab::RegionTag amrGrav_tag("Gravity_" + std::to_string(level));
      // fix fluxes on finer grids

      if (do_reflux)
      {
        for (int lev = level; lev < finest_level; lev++)
        {
          gravity->zero_phi_flux_reg(lev + 1);
        }
      }
      // swap grav data

      for (int lev = level; lev <= finest_level; lev++)
        get_level(lev).gravity->swap_time_levels(lev);

      // Solve for phi using the previous phi as a guess.

      int use_previous_phi_as_guess = 1;
      int ngrow_for_solve = 1; // iteration + stencil_deposition_width;
      gravity->multilevel_solve_for_old_phi(level, finest_level,
                                            ngrow_for_solve,
                                            use_previous_phi_as_guess);
    }
    BL_PROFILE_VAR_STOP(solve_for_old_phi);
  }

#endif

  BL_PROFILE_VAR("just_the_fdm", just_the_fdm);

#ifdef FDM
  const Real prev_time = state[Axion_Type].prevTime();
  const Real cur_time = state[Axion_Type].curTime();
  const Real a_old = get_comoving_a(prev_time);
  const Real a_new = get_comoving_a(cur_time);

  if (levelmethod[level] == FDlevel)
    get_level(level).advance_FDM_FD(time, dt, a_old, a_new);
#ifdef BL_USE_MPI
  else if (levelmethod[level] == PSlevel)
  
    if (do_grav){
    get_level(level).advance_FDM_PS(time, dt, a_old, a_new);}else{
    get_level(level).advance_FDM_PS_NG(time, dt, a_old, a_new);
    }
#endif
  else
    amrex::Abort("Nyx::advance_fdm -- can only do FD or PS");
#endif

#ifdef AXCOMPLEX
  if (levelmethod[level] == FDlevel)
    get_level(level).advance_AxComplex_FD(time, dt);
#ifdef BL_USE_MPI
  else if (levelmethod[level] == PSlevel)
    get_level(level).advance_AxComplex_PS(time, dt);
#endif
  else
    amrex::Abort("Nyx::advance_AxComplex -- can only do FD or PS");
#endif

#ifdef AXREAL
  if (levelmethod[level] == FDlevel)
    get_level(level).advance_AxReal_FD(time, dt);
#ifdef BL_USE_MPI
  else if (levelmethod[level] == PSlevel)
    get_level(level).advance_AxReal_PS(time, dt);
#endif
  else
    amrex::Abort("Nyx::advance_AxComplex -- can only do FD or PS");
#endif

#ifdef GRAVITY
  if (do_grav)
  {
    MultiFab::Copy(parent->getLevel(level).get_new_data(PhiGrav_Type),
                   parent->getLevel(level).get_old_data(PhiGrav_Type),
                   0, 0, 1, 0);

    // Solve for new Gravity

    BL_PROFILE_VAR("solve_for_new_phi", solve_for_new_phi);
    int use_previous_phi_as_guess = 1;

    MultiFab::RegionTag amrGrav_tag("Gravity_" + std::to_string(level));
    int fill_interior = 0;
    int grav_n_grow = 1;
    gravity->solve_for_new_phi(level, get_new_data(PhiGrav_Type),
                               gravity->get_grad_phi_curr(level),
                               fill_interior, grav_n_grow);

    BL_PROFILE_VAR_STOP(solve_for_new_phi);
  }
#endif

  BL_PROFILE_VAR_STOP(just_the_fdm);

  return dt;
}
#endif

/* ##########
 */
#ifndef NO_HYDRO // FW PhD
Real Nyx::advance_hydro_axionyx(
    Real time,
    Real dt,
    int iteration,
    int ncycle)
{

  //const Real prev_time = state[State_Type].prevTime();
  //const Real cur_time = state[State_Type].curTime();

  //const Real a_old = get_comoving_a(prev_time);
  //const Real a_new = get_comoving_a(cur_time);

  // Begin updating hydro and particles

  int stencil_deposition_width = 1;

  int stencil_interpolation_width = 1;

  int ghost_width = ncycle + stencil_deposition_width;

  int where_width = ghost_width + (1 - iteration) - 1;

  int grav_n_grow = ghost_width + (1 - iteration) + (iteration - 1) +
                    stencil_interpolation_width;


  // Sanity checks
  if (!do_hydro)
    amrex::Print()<<"In `advance_hydro_plus_particles` but `do_hydro` not true\n";

  //if (Nyx::theActiveParticles().size() <= 0)
  //  amrex::Print()<<"In `advance_hydro_plus_particles` but no active particles"<<std::endl;
    
  if (!do_grav)
        amrex::Print()<<"`do_grav` not set to true\n";

  const int finest_level = parent->finestLevel();
  int finest_level_to_advance;
  bool nosub = !parent->subCycle();

  if (nosub)
  {
    if (level > 0)
      return dt;

    finest_level_to_advance = finest_level;
  }
  else
  {
    // This level was advanced by a previous multilevel advance.
    if (level > 0 && ncycle == 1)
      return dt;

    // Find the finest level to advance
    int lev = level;
    while (lev < finest_level && parent->nCycle(lev + 1) == 1)
      lev++;
    finest_level_to_advance = lev;

    // We must setup virtual and Ghost Particles
    //
    // Setup the virtual particles that represent finer level particles
    //
    setup_virtual_particles();
    //
    // Setup ghost particles for use in finer levels. Note that Ghost particles
    // that will be used by this level have already been created, the
    // particles being set here are only used by finer levels.

    //
    for (int lev = level; lev <= finest_level_to_advance && lev < finest_level; lev++)
    {
      get_level(lev).setup_ghost_particles(ghost_width);
    }
  }

  Real dt_lev;

  if (level == 0 || iteration > 1)
  {

    for (int lev = level; lev <= finest_level; lev++)
    {
      dt_lev = parent->dtLevel(lev);
      for (int k = 0; k < NUM_STATE_TYPE; k++)
      {
        get_level(lev).state[k].allocOldData();
        get_level(lev).state[k].swapTimeLevels(dt_lev);
      }
    }
  }

  const Real prev_time = state[State_Type].prevTime();
  const Real cur_time = state[State_Type].curTime();
  
  // amrex::Print()<<"TIMES " << prev_time << ", " << cur_time << std::endl;
  
  const Real a_old = get_comoving_a(prev_time);
  const Real a_new = get_comoving_a(cur_time);

  //
  // Move current data to previous, clear current.
  // Don't do this if a coarser level has done this already.
  //

  if (do_grav)
  {
    BL_PROFILE_VAR("solve_for_old_phi", solve_for_old_phi);
    if (level == 0 || iteration > 1)
    {

      MultiFab::RegionTag amrGrav_tag("Gravity_" + std::to_string(level));

      if (do_reflux)
      {
        for (int lev = level; lev < finest_level; lev++)
        {
          gravity->zero_phi_flux_reg(lev + 1);
        }
      }

      // swap grav data
      for (int lev = level; lev <= finest_level; lev++)
        get_level(lev).gravity->swap_time_levels(lev);

      //
      // Solve for phi using the previous phi as a guess.
      //
      int use_previous_phi_as_guess = 1;
      int ngrow_for_solve = iteration + stencil_deposition_width;
      gravity->multilevel_solve_for_old_phi(level, finest_level,
                                            ngrow_for_solve,
                                            use_previous_phi_as_guess);
    }

  } // if (do_grav)
  //
  // Call the hydro advance at each level to be advanced
  //

  MultiFab::RegionTag amrhydro_tag("Hydro_" + std::to_string(level));

  for (int lev = level; lev <= finest_level_to_advance; lev++)
  {
    get_level(lev).strang_hydro(time, dt, a_old, a_new);
    // Original Axionyx

    if (levelmethod[level] == FDlevel)
      get_level(lev).advance_FDM_FD(time, dt, a_old, a_new);

    else if (levelmethod[level] == PSlevel)
      if (do_grav){
      get_level(level).advance_FDM_PS(time, dt, a_old, a_new);}else{
      get_level(level).advance_FDM_PS_NG(time, dt, a_old, a_new);
      }

    else
      amrex::Abort("Nyx::advance_fdm -- can only do FD or PS");
  }

  if (do_reflux)
  {
    for (int lev = level; lev < finest_level_to_advance; lev++)
    {
      get_level(lev).reflux();
    }
  }
  // Always average down the new state from finer to coarser.
  for (int lev = finest_level_to_advance - 1; lev >= level; lev--)
  {
    get_level(lev).average_down(State_Type);
    get_level(lev).average_down(DiagEOS_Type);
    get_level(lev).average_down(Axion_Type);
  }

  // Advancing the Particles Stitch.


  if (Nyx::theActiveParticles().size() > 0)
  {
    // Advance the particle velocities to the half-time and the positions to the new time
    // We use the cell-centered gravity to correctly interpolate onto particle locations
    const Real a_half = 0.5 * (a_old + a_new);

    if (particle_verbose && ParallelDescriptor::IOProcessor())
      amrex::Print() << "NBody: Particles x and xdot updated with MoveKickDrift!\n";

    MultiFab::RegionTag amrMoveKickDrift_tag("MoveKickDrift_" + std::to_string(level));

    for (int lev = level; lev <= finest_level_to_advance; lev++)
    {

      // We need grav_n_grow grow cells to track boundary particles
      const auto &ba = get_level(lev).get_new_data(PhiGrav_Type).boxArray();
      const auto &dm = get_level(lev).get_new_data(PhiGrav_Type).DistributionMap();

      MultiFab grav_vec_old(ba, dm, AMREX_SPACEDIM, grav_n_grow);

      get_level(lev).gravity->get_old_grav_vector(lev, grav_vec_old, time);

      amrex::Real acc_min = grav_vec_old.min(0, lev);
      amrex::Real acc_max = grav_vec_old.max(0, lev);

      // Print to the standard output
      if (particle_verbose && ParallelDescriptor::IOProcessor())
      {
        amrex::Print() << "FW is here to check. Acceleration MultiFab at level " << lev << ":"
                       << " Min = " << acc_min
                       << ", Max = " << acc_max << std::endl;
      }

      for (int i = 0; i < Nyx::theActiveParticles().size(); i++)
        Nyx::theActiveParticles()[i]->moveKickDrift(grav_vec_old, lev, dt, a_old, a_half, where_width);

      // Only need the coarsest virtual particles here.
      if (lev == level && level < finest_level)
        for (int i = 0; i < Nyx::theVirtualParticles().size(); i++)
          Nyx::theVirtualParticles()[i]->moveKickDrift(grav_vec_old, lev, dt, a_old, a_half, where_width);

      // Miiiight need all Ghosts
      // amrex::Print()<<"A Ghost! "<<Nyx::theGhostParticles().size() <<"\n";
      for (int i = 0; i < Nyx::theGhostParticles().size(); i++)
        Nyx::theGhostParticles()[i]->moveKickDrift(grav_vec_old, lev, dt, a_old, a_half, where_width);

    }
  }

#ifdef FDM
  for (int lev = level; lev <= finest_level_to_advance; lev++)
    if (levelmethod[lev] == FDlevel)
    {
      get_level(lev).advance_FDM_FD(time, dt, a_old, a_new);
      if (lev == parent->finestLevel())
        write_densmax();

      if (lev < parent->finestLevel() && false)
      {
        MultiFab &Ax_new = get_level(lev).get_new_data(Axion_Type);
        Nyx *nyx_level = dynamic_cast<Nyx *>(&(parent->getLevel(lev)));
        Real cdm_mass = nyx_level->vol_weight_sum("fdm_particle_mass_density", cur_time, true);
        Real fdm_mass = nyx_level->vol_weight_sum("AxDens", cur_time, true);
        Real offset = cdm_mass / fdm_mass;
        amrex::Print() << "FDM mass offset in advance_particles FDlevel " << lev << " =  "
                       << offset << " = " << cdm_mass << " / " << fdm_mass << '\n';

        for (MFIter mfi(Ax_new, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          const Box &bx = mfi.growntilebox();
          Array4<Real> const &fab = Ax_new.array(mfi);
          AMREX_PARALLEL_FOR_3D(bx, i, j, k,
                                {
                                  fab(i, j, k, Nyx::AxRe) *= std::sqrt(offset);
                                  fab(i, j, k, Nyx::AxIm) *= std::sqrt(offset);
                                  fab(i, j, k, Nyx::AxDens) *= offset;
                                });
        }
      }
    }
    else if (levelmethod[lev] == PSlevel)
      if (do_grav){
      get_level(level).advance_FDM_PS(time, dt, a_old, a_new);}else{
      get_level(level).advance_FDM_PS_NG(time, dt, a_old, a_new);
      }

  for (int lev = finest_level_to_advance - 1; lev >= level; lev--)
    get_level(lev).average_down(Axion_Type);
#endif

  //

  // Here we use the "old" phi from the current time step as a guess for this
  // solve
  //
  if (do_grav)
  {

    // std::cout << "Evaluating New Gravity Here" << std::endl;
    for (int lev = level; lev <= finest_level_to_advance; lev++)
    {
      MultiFab::Copy(parent->getLevel(lev).get_new_data(PhiGrav_Type),
                     parent->getLevel(lev).get_old_data(PhiGrav_Type),
                     0, 0, 1, 0);
    }

    // Solve for new Gravity
    int use_previous_phi_as_guess = 1;
    if (finest_level_to_advance > level)
    {
      MultiFab::RegionTag amrGrav_tag("Gravity_" + std::to_string(level));
      // The particle may be as many as "iteration" ghost cells out
      int ngrow_for_solve = iteration + stencil_deposition_width;
      gravity->multilevel_solve_for_new_phi(level, finest_level_to_advance,
                                            ngrow_for_solve,
                                            use_previous_phi_as_guess);
    }
    else
    {
      MultiFab::RegionTag amrGrav_tag("Gravity_" + std::to_string(level));
      int fill_interior = 0;
      gravity->solve_for_new_phi(level, get_new_data(PhiGrav_Type),
                                 gravity->get_grad_phi_curr(level),
                                 fill_interior, grav_n_grow);
    }
    // Reflux
    if (do_reflux)
    {
      MultiFab::RegionTag amrGrav_tag("Gravity_" + std::to_string(level));
      for (int lev = level; lev <= finest_level_to_advance; lev++)
      {
        gravity->add_to_fluxes(lev, iteration, ncycle);
      }
    }

    //
    // Now do corrector part of source term update
    //
    for (int lev = level; lev <= finest_level_to_advance; lev++)
    {
      // std::cout << "Evaluating Source Corrections" << std::endl;
      MultiFab::RegionTag amrGrav_tag("Gravity_" + std::to_string(lev));

      // Now do corrector part of source term update
      correct_gsrc(lev, time, prev_time, cur_time, dt);

      MultiFab &S_new = get_level(lev).get_new_data(State_Type);
      MultiFab &D_new = get_level(lev).get_new_data(DiagEOS_Type);
      MultiFab &Ax_new = get_level(lev).get_new_data(Axion_Type);

      // First reset internal energy before call to compute_temp
      MultiFab reset_e_src(S_new.boxArray(), S_new.DistributionMap(), 1, NUM_GROW);
      reset_e_src.setVal(0.0);
      get_level(lev).reset_internal_energy(S_new, D_new, reset_e_src);
      get_level(lev).compute_new_temp(S_new, D_new);
    }

    // Must average down again after doing the gravity correction;
    //      always average down from finer to coarser.
    // Here we average down both the new state and phi and gravity.
    for (int lev = finest_level_to_advance - 1; lev >= level; lev--)
      get_level(lev).average_down();

    if (Nyx::theActiveParticles().size() > 0)
    {
      // Advance the particle velocities by dt/2 to the new time. We use the
      // cell-centered gravity to correctly interpolate onto particle
      // locations.
      MultiFab::RegionTag amrMoveKickDrift_tag("MoveKick_" + std::to_string(level));
      const Real a_half = 0.5 * (a_old + a_new);

      if (particle_verbose && ParallelDescriptor::IOProcessor())
        amrex::Print() << "NBody: Particles xdot updated with MoveKick\n";

      for (int lev = level; lev <= finest_level_to_advance; lev++)
      {
        const auto &ba = get_level(lev).get_new_data(State_Type).boxArray();
        const auto &dm = get_level(lev).get_new_data(State_Type).DistributionMap();
        MultiFab grav_vec_new(ba, dm, AMREX_SPACEDIM, grav_n_grow);
        get_level(lev).gravity->get_new_grav_vector(lev, grav_vec_new, cur_time);

        for (int i = 0; i < Nyx::theActiveParticles().size(); i++)
          Nyx::theActiveParticles()[i]->moveKick(grav_vec_new, lev, dt, a_new, a_half);

        // Virtual particles will be recreated, so we need not kick them.

        // Ghost particles need to be kicked except during the final iteration.
        if (iteration != ncycle)
          for (int i = 0; i < Nyx::theGhostParticles().size(); i++)
            Nyx::theGhostParticles()[i]->moveKick(grav_vec_new, lev, dt, a_new, a_half);
      }
    }

  } // do_grav

  //
  // Synchronize Energies
  //
  for (int lev = level; lev <= finest_level_to_advance; lev++)
  {
    MultiFab::RegionTag amrReset_tag("Reset_" + std::to_string(lev));
    MultiFab &S_new = get_level(lev).get_new_data(State_Type);
    MultiFab &D_new = get_level(lev).get_new_data(DiagEOS_Type);
    MultiFab &Ax_new = get_level(lev).get_new_data(Axion_Type);

    get_level(lev).reset_internal_energy_nostore(S_new, D_new);
  }

  return dt;
}
#endif