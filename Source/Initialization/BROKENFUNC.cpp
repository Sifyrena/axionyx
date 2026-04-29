
Nyx::check_initial_species ()
{

std::
    // Verify that the sum of (rho X)_i = rho at every cell.
#ifdef CONST_SPECIES
#ifndef AMREX_USE_FLOAT
    if (amrex::Math::abs(1.0 - h_species - he_species) > 1.e-8)
        amrex::Abort("Error:: Failed check of initial species summing to 1");
#else
    if (amrex::Math::abs(1.0 - h_species - he_species) > 1.e-6)
        amrex::Abort("Error:: Failed check of initial species summing to 1");
#endif
#else
    amrex::ReduceOps<ReduceOpMax> reduce_op;
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
        if (get<0>(hv) > 1.e-8)
            amrex::Abort("Error:: Failed check of initial species summing to 1");
#else
        if (get<0>(hv) > 1.e-6)
            amrex::Abort("Error:: Failed check of initial species summing to 1");
#endif
    }
#endif
}