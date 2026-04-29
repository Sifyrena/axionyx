#ifdef FDM
 
#include "Nyx.H"
#include <AMReX_MultiFab.H>
#ifdef GRAVITY
#include "Gravity.H"
#endif

using namespace amrex;

void
Nyx::advance_FDM_FD (Real time,
                      Real dt,
                      Real a_old,
                      Real a_new)
{

	Print() << "Finite Difference Advance was called!" << std::endl;  
    if (verbose && ParallelDescriptor::IOProcessor() ){
	std::cout << "Advancing the axions at level " << level <<  "...\n";
    }
    const Real a_half = get_comoving_a(time+0.5*dt);
    const Real* dx      = geom.CellSize();
    const Real invdeltasq_old  = 1.0 / ( a_old  * dx[0] ) / ( a_old  * dx[0] );                                                                                                       
    const Real invdeltasq_half = 1.0 / ( a_half * dx[0] ) / ( a_half * dx[0] );                                                                                                       
    const Real invdeltasq_new  = 1.0 / ( a_new  * dx[0] ) / ( a_new  * dx[0] );

    MultiFab&  Ax_old = get_old_data(Axion_Type);
    MultiFab&  Ax_new = get_new_data(Axion_Type);

#ifdef DEBUG
    if (Ax_old.contains_nan(0, Ax_old.nComp(), 0))
      {
        for (int i = 0; i < Ax_old.nComp(); i++)
          {
            if (Ax_old.contains_nan(i,1,0))
              {
		std::cout << "Testing component i for NaNs: " << i << std::endl;
		amrex::Abort("Ax_old has NaNs in this component::advance_FDM_FD()");
              }
          }
      }
#endif

#ifdef GRAVITY
    MultiFab& Phi_old = get_old_data(PhiGrav_Type);
    for (FillPatchIterator 
	   fpi(*this, Ax_old, 4, time, Axion_Type,   0, Nyx::NUM_AX),
	   pfpi(*this, Phi_old, 3, time, PhiGrav_Type, 0, 1);
	 fpi.isValid() && pfpi.isValid(); ++fpi,++pfpi)
#else
    for (FillPatchIterator 
	   fpi(*this, Ax_old, 4, time, Axion_Type,   0, Nyx::NUM_AX);
	 fpi.isValid(); ++fpi)
#endif
      {
	const Box& bxfour  = fpi.validbox();
	const Box& bxthree = grow(bxfour,1);
	const Box& bxtwo   = grow(bxthree,1);
	const Box& bxone   = grow(bxtwo,1);
	auto const ax_in   = fpi().array();
	auto const ax_out  = Ax_new[fpi].array();
#ifdef GRAVITY
	auto const phi     = pfpi().array();
#endif
	FArrayBox kr_one_fab(bxone);
	Array4<Real> const kr_one = kr_one_fab.array();
	FArrayBox ki_one_fab(bxone);
	Array4<Real> const ki_one = ki_one_fab.array();

	ParallelFor(bxone,
                           [=] AMREX_GPU_DEVICE (int i, int j, int k)
			   {
			     kr_one(i,j,k) =-hbaroverm*(6.0*ax_in(i+1,j,k,Nyx::AxIm)
							+6.0*ax_in(i-1,j,k,Nyx::AxIm)
							+6.0*ax_in(i,j+1,k,Nyx::AxIm)
							+6.0*ax_in(i,j-1,k,Nyx::AxIm)
							+6.0*ax_in(i,j,k+1,Nyx::AxIm)
							+6.0*ax_in(i,j,k-1,Nyx::AxIm)
							+3.0*ax_in(i+1,j+1,k,Nyx::AxIm)
							+3.0*ax_in(i+1,j-1,k,Nyx::AxIm)
							+3.0*ax_in(i-1,j+1,k,Nyx::AxIm)
							+3.0*ax_in(i-1,j-1,k,Nyx::AxIm)
							+3.0*ax_in(i+1,j,k+1,Nyx::AxIm)
							+3.0*ax_in(i+1,j,k-1,Nyx::AxIm)
							+3.0*ax_in(i-1,j,k+1,Nyx::AxIm)
							+3.0*ax_in(i-1,j,k-1,Nyx::AxIm)
							+3.0*ax_in(i,j+1,k+1,Nyx::AxIm)
							+3.0*ax_in(i,j+1,k-1,Nyx::AxIm)
							+3.0*ax_in(i,j-1,k+1,Nyx::AxIm)
							+3.0*ax_in(i,j-1,k-1,Nyx::AxIm)
							+2.0*ax_in(i+1,j+1,k+1,Nyx::AxIm)
							+2.0*ax_in(i+1,j+1,k-1,Nyx::AxIm)
							+2.0*ax_in(i+1,j-1,k+1,Nyx::AxIm)
							+2.0*ax_in(i+1,j-1,k-1,Nyx::AxIm)
							+2.0*ax_in(i-1,j+1,k+1,Nyx::AxIm)
							+2.0*ax_in(i-1,j+1,k-1,Nyx::AxIm)
							+2.0*ax_in(i-1,j-1,k+1,Nyx::AxIm)
							+2.0*ax_in(i-1,j-1,k-1,Nyx::AxIm)
							-88.0*ax_in(i,j,k,Nyx::AxIm))
			       *invdeltasq_old/52.0
#ifdef GRAVITY
			       - phi(i,j,k)/hbaroverm*(ax_in(i,j,k,Nyx::AxIm))
#endif
			       ;

			     ki_one(i,j,k) = hbaroverm*(6.0*ax_in(i+1,j,k,Nyx::AxRe)
							+6.0*ax_in(i-1,j,k,Nyx::AxRe)
							+6.0*ax_in(i,j+1,k,Nyx::AxRe)
							+6.0*ax_in(i,j-1,k,Nyx::AxRe)
							+6.0*ax_in(i,j,k+1,Nyx::AxRe)
							+6.0*ax_in(i,j,k-1,Nyx::AxRe)
							+3.0*ax_in(i+1,j+1,k,Nyx::AxRe)
							+3.0*ax_in(i+1,j-1,k,Nyx::AxRe)
							+3.0*ax_in(i-1,j+1,k,Nyx::AxRe)
							+3.0*ax_in(i-1,j-1,k,Nyx::AxRe)
							+3.0*ax_in(i+1,j,k+1,Nyx::AxRe)
							+3.0*ax_in(i+1,j,k-1,Nyx::AxRe)
							+3.0*ax_in(i-1,j,k+1,Nyx::AxRe)
							+3.0*ax_in(i-1,j,k-1,Nyx::AxRe)
							+3.0*ax_in(i,j+1,k+1,Nyx::AxRe)
							+3.0*ax_in(i,j+1,k-1,Nyx::AxRe)
							+3.0*ax_in(i,j-1,k+1,Nyx::AxRe)
							+3.0*ax_in(i,j-1,k-1,Nyx::AxRe)
							+2.0*ax_in(i+1,j+1,k+1,Nyx::AxRe)
							+2.0*ax_in(i+1,j+1,k-1,Nyx::AxRe)
							+2.0*ax_in(i+1,j-1,k+1,Nyx::AxRe)
							+2.0*ax_in(i+1,j-1,k-1,Nyx::AxRe)
							+2.0*ax_in(i-1,j+1,k+1,Nyx::AxRe)
							+2.0*ax_in(i-1,j+1,k-1,Nyx::AxRe)
							+2.0*ax_in(i-1,j-1,k+1,Nyx::AxRe)
							+2.0*ax_in(i-1,j-1,k-1,Nyx::AxRe)
							-88.0*ax_in(i,j,k,Nyx::AxRe))
			       *invdeltasq_old/52.0
#ifdef GRAVITY
			       + phi(i,j,k)/hbaroverm*(ax_in(i,j,k,Nyx::AxRe))
#endif
			       ;
			       
			      });

	FArrayBox kr_two_fab(bxtwo);
	Array4<Real> const kr_two = kr_two_fab.array();
	FArrayBox ki_two_fab(bxtwo);
	Array4<Real> const ki_two = ki_two_fab.array();

	ParallelFor(bxtwo,
 			   [=] AMREX_GPU_DEVICE (int i, int j, int k)
			   {
			     kr_two(i,j,k) =-hbaroverm*(6.0*(ax_in(i+1,j,k,Nyx::AxIm)+ki_one(i+1,j,k)*dt/2.0)
							+6.0*(ax_in(i-1,j,k,Nyx::AxIm)+ki_one(i-1,j,k)*dt/2.0)
							+6.0*(ax_in(i,j+1,k,Nyx::AxIm)+ki_one(i,j+1,k)*dt/2.0)
							+6.0*(ax_in(i,j-1,k,Nyx::AxIm)+ki_one(i,j-1,k)*dt/2.0)
							+6.0*(ax_in(i,j,k+1,Nyx::AxIm)+ki_one(i,j,k+1)*dt/2.0)
							+6.0*(ax_in(i,j,k-1,Nyx::AxIm)+ki_one(i,j,k-1)*dt/2.0)
							+3.0*(ax_in(i+1,j+1,k,Nyx::AxIm)+ki_one(i+1,j+1,k)*dt/2.0)
							+3.0*(ax_in(i+1,j-1,k,Nyx::AxIm)+ki_one(i+1,j-1,k)*dt/2.0)
							+3.0*(ax_in(i-1,j+1,k,Nyx::AxIm)+ki_one(i-1,j+1,k)*dt/2.0)
							+3.0*(ax_in(i-1,j-1,k,Nyx::AxIm)+ki_one(i-1,j-1,k)*dt/2.0)
							+3.0*(ax_in(i+1,j,k+1,Nyx::AxIm)+ki_one(i+1,j,k+1)*dt/2.0)
							+3.0*(ax_in(i+1,j,k-1,Nyx::AxIm)+ki_one(i+1,j,k-1)*dt/2.0)
							+3.0*(ax_in(i-1,j,k+1,Nyx::AxIm)+ki_one(i-1,j,k+1)*dt/2.0)
							+3.0*(ax_in(i-1,j,k-1,Nyx::AxIm)+ki_one(i-1,j,k-1)*dt/2.0)
							+3.0*(ax_in(i,j+1,k+1,Nyx::AxIm)+ki_one(i,j+1,k+1)*dt/2.0)
							+3.0*(ax_in(i,j+1,k-1,Nyx::AxIm)+ki_one(i,j+1,k-1)*dt/2.0)
							+3.0*(ax_in(i,j-1,k+1,Nyx::AxIm)+ki_one(i,j-1,k+1)*dt/2.0)
							+3.0*(ax_in(i,j-1,k-1,Nyx::AxIm)+ki_one(i,j-1,k-1)*dt/2.0)
							+2.0*(ax_in(i+1,j+1,k+1,Nyx::AxIm)+ki_one(i+1,j+1,k+1)*dt/2.0)
							+2.0*(ax_in(i+1,j+1,k-1,Nyx::AxIm)+ki_one(i+1,j+1,k-1)*dt/2.0)
							+2.0*(ax_in(i+1,j-1,k+1,Nyx::AxIm)+ki_one(i+1,j-1,k+1)*dt/2.0)
							+2.0*(ax_in(i+1,j-1,k-1,Nyx::AxIm)+ki_one(i+1,j-1,k-1)*dt/2.0)
							+2.0*(ax_in(i-1,j+1,k+1,Nyx::AxIm)+ki_one(i-1,j+1,k+1)*dt/2.0)
							+2.0*(ax_in(i-1,j+1,k-1,Nyx::AxIm)+ki_one(i-1,j+1,k-1)*dt/2.0)
							+2.0*(ax_in(i-1,j-1,k+1,Nyx::AxIm)+ki_one(i-1,j-1,k+1)*dt/2.0)
							+2.0*(ax_in(i-1,j-1,k-1,Nyx::AxIm)+ki_one(i-1,j-1,k-1)*dt/2.0)
							-88.0*(ax_in(i,j,k,Nyx::AxIm)+ki_one(i,j,k)*dt/2.0))
			       *invdeltasq_half/52.0
#ifdef GRAVITY
			       - phi(i,j,k)/hbaroverm*(ax_in(i,j,k,Nyx::AxIm)+ki_one(i,j,k)*dt/2.0)
#endif
			       ;

			     ki_two(i,j,k) = hbaroverm*(6.0*(ax_in(i+1,j,k,Nyx::AxRe)+kr_one(i+1,j,k)*dt/2.0)
							+6.0*(ax_in(i-1,j,k,Nyx::AxRe)+kr_one(i-1,j,k)*dt/2.0)
							+6.0*(ax_in(i,j+1,k,Nyx::AxRe)+kr_one(i,j+1,k)*dt/2.0)
							+6.0*(ax_in(i,j-1,k,Nyx::AxRe)+kr_one(i,j-1,k)*dt/2.0)
							+6.0*(ax_in(i,j,k+1,Nyx::AxRe)+kr_one(i,j,k+1)*dt/2.0)
							+6.0*(ax_in(i,j,k-1,Nyx::AxRe)+kr_one(i,j,k-1)*dt/2.0)
							+3.0*(ax_in(i+1,j+1,k,Nyx::AxRe)+kr_one(i+1,j+1,k)*dt/2.0)
							+3.0*(ax_in(i+1,j-1,k,Nyx::AxRe)+kr_one(i+1,j-1,k)*dt/2.0)
							+3.0*(ax_in(i-1,j+1,k,Nyx::AxRe)+kr_one(i-1,j+1,k)*dt/2.0)
							+3.0*(ax_in(i-1,j-1,k,Nyx::AxRe)+kr_one(i-1,j-1,k)*dt/2.0)
							+3.0*(ax_in(i+1,j,k+1,Nyx::AxRe)+kr_one(i+1,j,k+1)*dt/2.0)
							+3.0*(ax_in(i+1,j,k-1,Nyx::AxRe)+kr_one(i+1,j,k-1)*dt/2.0)
							+3.0*(ax_in(i-1,j,k+1,Nyx::AxRe)+kr_one(i-1,j,k+1)*dt/2.0)
							+3.0*(ax_in(i-1,j,k-1,Nyx::AxRe)+kr_one(i-1,j,k-1)*dt/2.0)
							+3.0*(ax_in(i,j+1,k+1,Nyx::AxRe)+kr_one(i,j+1,k+1)*dt/2.0)
							+3.0*(ax_in(i,j+1,k-1,Nyx::AxRe)+kr_one(i,j+1,k-1)*dt/2.0)
							+3.0*(ax_in(i,j-1,k+1,Nyx::AxRe)+kr_one(i,j-1,k+1)*dt/2.0)
							+3.0*(ax_in(i,j-1,k-1,Nyx::AxRe)+kr_one(i,j-1,k-1)*dt/2.0)
							+2.0*(ax_in(i+1,j+1,k+1,Nyx::AxRe)+kr_one(i+1,j+1,k+1)*dt/2.0)
							+2.0*(ax_in(i+1,j+1,k-1,Nyx::AxRe)+kr_one(i+1,j+1,k-1)*dt/2.0)
							+2.0*(ax_in(i+1,j-1,k+1,Nyx::AxRe)+kr_one(i+1,j-1,k+1)*dt/2.0)
							+2.0*(ax_in(i+1,j-1,k-1,Nyx::AxRe)+kr_one(i+1,j-1,k-1)*dt/2.0)
							+2.0*(ax_in(i-1,j+1,k+1,Nyx::AxRe)+kr_one(i-1,j+1,k+1)*dt/2.0)
							+2.0*(ax_in(i-1,j+1,k-1,Nyx::AxRe)+kr_one(i-1,j+1,k-1)*dt/2.0)
							+2.0*(ax_in(i-1,j-1,k+1,Nyx::AxRe)+kr_one(i-1,j-1,k+1)*dt/2.0)
							+2.0*(ax_in(i-1,j-1,k-1,Nyx::AxRe)+kr_one(i-1,j-1,k-1)*dt/2.0)
							-88.0*(ax_in(i,j,k,Nyx::AxRe)+kr_one(i,j,k)*dt/2.0))
			       *invdeltasq_half/52.0
#ifdef GRAVITY
			       + phi(i,j,k)/hbaroverm*(ax_in(i,j,k,Nyx::AxRe)+kr_one(i,j,k)*dt/2.0)
#endif
			       ;

			       });

	FArrayBox kr_three_fab(bxthree);
	Array4<Real> const kr_three = kr_three_fab.array();
	FArrayBox ki_three_fab(bxthree);
	Array4<Real> const ki_three = ki_three_fab.array();

	ParallelFor(bxthree,
 			   [=] AMREX_GPU_DEVICE (int i, int j, int k)
			   {
			     kr_three(i,j,k) =-hbaroverm*(6.0*(ax_in(i+1,j,k,Nyx::AxIm)+ki_two(i+1,j,k)*dt/2.0)
							+6.0*(ax_in(i-1,j,k,Nyx::AxIm)+ki_two(i-1,j,k)*dt/2.0)
							+6.0*(ax_in(i,j+1,k,Nyx::AxIm)+ki_two(i,j+1,k)*dt/2.0)
							+6.0*(ax_in(i,j-1,k,Nyx::AxIm)+ki_two(i,j-1,k)*dt/2.0)
							+6.0*(ax_in(i,j,k+1,Nyx::AxIm)+ki_two(i,j,k+1)*dt/2.0)
							+6.0*(ax_in(i,j,k-1,Nyx::AxIm)+ki_two(i,j,k-1)*dt/2.0)
							+3.0*(ax_in(i+1,j+1,k,Nyx::AxIm)+ki_two(i+1,j+1,k)*dt/2.0)
							+3.0*(ax_in(i+1,j-1,k,Nyx::AxIm)+ki_two(i+1,j-1,k)*dt/2.0)
							+3.0*(ax_in(i-1,j+1,k,Nyx::AxIm)+ki_two(i-1,j+1,k)*dt/2.0)
							+3.0*(ax_in(i-1,j-1,k,Nyx::AxIm)+ki_two(i-1,j-1,k)*dt/2.0)
							+3.0*(ax_in(i+1,j,k+1,Nyx::AxIm)+ki_two(i+1,j,k+1)*dt/2.0)
							+3.0*(ax_in(i+1,j,k-1,Nyx::AxIm)+ki_two(i+1,j,k-1)*dt/2.0)
							+3.0*(ax_in(i-1,j,k+1,Nyx::AxIm)+ki_two(i-1,j,k+1)*dt/2.0)
							+3.0*(ax_in(i-1,j,k-1,Nyx::AxIm)+ki_two(i-1,j,k-1)*dt/2.0)
							+3.0*(ax_in(i,j+1,k+1,Nyx::AxIm)+ki_two(i,j+1,k+1)*dt/2.0)
							+3.0*(ax_in(i,j+1,k-1,Nyx::AxIm)+ki_two(i,j+1,k-1)*dt/2.0)
							+3.0*(ax_in(i,j-1,k+1,Nyx::AxIm)+ki_two(i,j-1,k+1)*dt/2.0)
							+3.0*(ax_in(i,j-1,k-1,Nyx::AxIm)+ki_two(i,j-1,k-1)*dt/2.0)
							+2.0*(ax_in(i+1,j+1,k+1,Nyx::AxIm)+ki_two(i+1,j+1,k+1)*dt/2.0)
							+2.0*(ax_in(i+1,j+1,k-1,Nyx::AxIm)+ki_two(i+1,j+1,k-1)*dt/2.0)
							+2.0*(ax_in(i+1,j-1,k+1,Nyx::AxIm)+ki_two(i+1,j-1,k+1)*dt/2.0)
							+2.0*(ax_in(i+1,j-1,k-1,Nyx::AxIm)+ki_two(i+1,j-1,k-1)*dt/2.0)
							+2.0*(ax_in(i-1,j+1,k+1,Nyx::AxIm)+ki_two(i-1,j+1,k+1)*dt/2.0)
							+2.0*(ax_in(i-1,j+1,k-1,Nyx::AxIm)+ki_two(i-1,j+1,k-1)*dt/2.0)
							+2.0*(ax_in(i-1,j-1,k+1,Nyx::AxIm)+ki_two(i-1,j-1,k+1)*dt/2.0)
							+2.0*(ax_in(i-1,j-1,k-1,Nyx::AxIm)+ki_two(i-1,j-1,k-1)*dt/2.0)
							-88.0*(ax_in(i,j,k,Nyx::AxIm)+ki_two(i,j,k)*dt/2.0))
			       *invdeltasq_half/52.0
#ifdef GRAVITY
			       - phi(i,j,k)/hbaroverm*(ax_in(i,j,k,Nyx::AxIm)+ki_two(i,j,k)*dt/2.0)
#endif
			       ;

			     ki_three(i,j,k) = hbaroverm*(6.0*(ax_in(i+1,j,k,Nyx::AxRe)+kr_two(i+1,j,k)*dt/2.0)
							+6.0*(ax_in(i-1,j,k,Nyx::AxRe)+kr_two(i-1,j,k)*dt/2.0)
							+6.0*(ax_in(i,j+1,k,Nyx::AxRe)+kr_two(i,j+1,k)*dt/2.0)
							+6.0*(ax_in(i,j-1,k,Nyx::AxRe)+kr_two(i,j-1,k)*dt/2.0)
							+6.0*(ax_in(i,j,k+1,Nyx::AxRe)+kr_two(i,j,k+1)*dt/2.0)
							+6.0*(ax_in(i,j,k-1,Nyx::AxRe)+kr_two(i,j,k-1)*dt/2.0)
							+3.0*(ax_in(i+1,j+1,k,Nyx::AxRe)+kr_two(i+1,j+1,k)*dt/2.0)
							+3.0*(ax_in(i+1,j-1,k,Nyx::AxRe)+kr_two(i+1,j-1,k)*dt/2.0)
							+3.0*(ax_in(i-1,j+1,k,Nyx::AxRe)+kr_two(i-1,j+1,k)*dt/2.0)
							+3.0*(ax_in(i-1,j-1,k,Nyx::AxRe)+kr_two(i-1,j-1,k)*dt/2.0)
							+3.0*(ax_in(i+1,j,k+1,Nyx::AxRe)+kr_two(i+1,j,k+1)*dt/2.0)
							+3.0*(ax_in(i+1,j,k-1,Nyx::AxRe)+kr_two(i+1,j,k-1)*dt/2.0)
							+3.0*(ax_in(i-1,j,k+1,Nyx::AxRe)+kr_two(i-1,j,k+1)*dt/2.0)
							+3.0*(ax_in(i-1,j,k-1,Nyx::AxRe)+kr_two(i-1,j,k-1)*dt/2.0)
							+3.0*(ax_in(i,j+1,k+1,Nyx::AxRe)+kr_two(i,j+1,k+1)*dt/2.0)
							+3.0*(ax_in(i,j+1,k-1,Nyx::AxRe)+kr_two(i,j+1,k-1)*dt/2.0)
							+3.0*(ax_in(i,j-1,k+1,Nyx::AxRe)+kr_two(i,j-1,k+1)*dt/2.0)
							+3.0*(ax_in(i,j-1,k-1,Nyx::AxRe)+kr_two(i,j-1,k-1)*dt/2.0)
							+2.0*(ax_in(i+1,j+1,k+1,Nyx::AxRe)+kr_two(i+1,j+1,k+1)*dt/2.0)
							+2.0*(ax_in(i+1,j+1,k-1,Nyx::AxRe)+kr_two(i+1,j+1,k-1)*dt/2.0)
							+2.0*(ax_in(i+1,j-1,k+1,Nyx::AxRe)+kr_two(i+1,j-1,k+1)*dt/2.0)
							+2.0*(ax_in(i+1,j-1,k-1,Nyx::AxRe)+kr_two(i+1,j-1,k-1)*dt/2.0)
							+2.0*(ax_in(i-1,j+1,k+1,Nyx::AxRe)+kr_two(i-1,j+1,k+1)*dt/2.0)
							+2.0*(ax_in(i-1,j+1,k-1,Nyx::AxRe)+kr_two(i-1,j+1,k-1)*dt/2.0)
							+2.0*(ax_in(i-1,j-1,k+1,Nyx::AxRe)+kr_two(i-1,j-1,k+1)*dt/2.0)
							+2.0*(ax_in(i-1,j-1,k-1,Nyx::AxRe)+kr_two(i-1,j-1,k-1)*dt/2.0)
							-88.0*(ax_in(i,j,k,Nyx::AxRe)+kr_two(i,j,k)*dt/2.0))
			       *invdeltasq_half/52.0
#ifdef GRAVITY
			       + phi(i,j,k)/hbaroverm*(ax_in(i,j,k,Nyx::AxRe)+kr_two(i,j,k)*dt/2.0)
#endif
			       ;

			       });

	FArrayBox kr_four_fab(bxfour);
	Array4<Real> const kr_four = kr_four_fab.array();
	FArrayBox ki_four_fab(bxfour);
	Array4<Real> const ki_four = ki_four_fab.array();

	ParallelFor(bxfour,
 			   [=] AMREX_GPU_DEVICE (int i, int j, int k)
			   {
			     kr_four(i,j,k) =-hbaroverm*(6.0*(ax_in(i+1,j,k,Nyx::AxIm)+ki_three(i+1,j,k)*dt)
							+6.0*(ax_in(i-1,j,k,Nyx::AxIm)+ki_three(i-1,j,k)*dt)
							+6.0*(ax_in(i,j+1,k,Nyx::AxIm)+ki_three(i,j+1,k)*dt)
							+6.0*(ax_in(i,j-1,k,Nyx::AxIm)+ki_three(i,j-1,k)*dt)
							+6.0*(ax_in(i,j,k+1,Nyx::AxIm)+ki_three(i,j,k+1)*dt)
							+6.0*(ax_in(i,j,k-1,Nyx::AxIm)+ki_three(i,j,k-1)*dt)
							+3.0*(ax_in(i+1,j+1,k,Nyx::AxIm)+ki_three(i+1,j+1,k)*dt)
							+3.0*(ax_in(i+1,j-1,k,Nyx::AxIm)+ki_three(i+1,j-1,k)*dt)
							+3.0*(ax_in(i-1,j+1,k,Nyx::AxIm)+ki_three(i-1,j+1,k)*dt)
							+3.0*(ax_in(i-1,j-1,k,Nyx::AxIm)+ki_three(i-1,j-1,k)*dt)
							+3.0*(ax_in(i+1,j,k+1,Nyx::AxIm)+ki_three(i+1,j,k+1)*dt)
							+3.0*(ax_in(i+1,j,k-1,Nyx::AxIm)+ki_three(i+1,j,k-1)*dt)
							+3.0*(ax_in(i-1,j,k+1,Nyx::AxIm)+ki_three(i-1,j,k+1)*dt)
							+3.0*(ax_in(i-1,j,k-1,Nyx::AxIm)+ki_three(i-1,j,k-1)*dt)
							+3.0*(ax_in(i,j+1,k+1,Nyx::AxIm)+ki_three(i,j+1,k+1)*dt)
							+3.0*(ax_in(i,j+1,k-1,Nyx::AxIm)+ki_three(i,j+1,k-1)*dt)
							+3.0*(ax_in(i,j-1,k+1,Nyx::AxIm)+ki_three(i,j-1,k+1)*dt)
							+3.0*(ax_in(i,j-1,k-1,Nyx::AxIm)+ki_three(i,j-1,k-1)*dt)
							+2.0*(ax_in(i+1,j+1,k+1,Nyx::AxIm)+ki_three(i+1,j+1,k+1)*dt)
							+2.0*(ax_in(i+1,j+1,k-1,Nyx::AxIm)+ki_three(i+1,j+1,k-1)*dt)
							+2.0*(ax_in(i+1,j-1,k+1,Nyx::AxIm)+ki_three(i+1,j-1,k+1)*dt)
							+2.0*(ax_in(i+1,j-1,k-1,Nyx::AxIm)+ki_three(i+1,j-1,k-1)*dt)
							+2.0*(ax_in(i-1,j+1,k+1,Nyx::AxIm)+ki_three(i-1,j+1,k+1)*dt)
							+2.0*(ax_in(i-1,j+1,k-1,Nyx::AxIm)+ki_three(i-1,j+1,k-1)*dt)
							+2.0*(ax_in(i-1,j-1,k+1,Nyx::AxIm)+ki_three(i-1,j-1,k+1)*dt)
							+2.0*(ax_in(i-1,j-1,k-1,Nyx::AxIm)+ki_three(i-1,j-1,k-1)*dt)
							-88.0*(ax_in(i,j,k,Nyx::AxIm)+ki_three(i,j,k)*dt))
			       *invdeltasq_new/52.0
#ifdef GRAVITY
			       - phi(i,j,k)/hbaroverm*(ax_in(i,j,k,Nyx::AxIm)+ki_three(i,j,k)*dt)
#endif
			       ;

			     ki_four(i,j,k) = hbaroverm*(6.0*(ax_in(i+1,j,k,Nyx::AxRe)+kr_three(i+1,j,k)*dt)
							+6.0*(ax_in(i-1,j,k,Nyx::AxRe)+kr_three(i-1,j,k)*dt)
							+6.0*(ax_in(i,j+1,k,Nyx::AxRe)+kr_three(i,j+1,k)*dt)
							+6.0*(ax_in(i,j-1,k,Nyx::AxRe)+kr_three(i,j-1,k)*dt)
							+6.0*(ax_in(i,j,k+1,Nyx::AxRe)+kr_three(i,j,k+1)*dt)
							+6.0*(ax_in(i,j,k-1,Nyx::AxRe)+kr_three(i,j,k-1)*dt)
							+3.0*(ax_in(i+1,j+1,k,Nyx::AxRe)+kr_three(i+1,j+1,k)*dt)
							+3.0*(ax_in(i+1,j-1,k,Nyx::AxRe)+kr_three(i+1,j-1,k)*dt)
							+3.0*(ax_in(i-1,j+1,k,Nyx::AxRe)+kr_three(i-1,j+1,k)*dt)
							+3.0*(ax_in(i-1,j-1,k,Nyx::AxRe)+kr_three(i-1,j-1,k)*dt)
							+3.0*(ax_in(i+1,j,k+1,Nyx::AxRe)+kr_three(i+1,j,k+1)*dt)
							+3.0*(ax_in(i+1,j,k-1,Nyx::AxRe)+kr_three(i+1,j,k-1)*dt)
							+3.0*(ax_in(i-1,j,k+1,Nyx::AxRe)+kr_three(i-1,j,k+1)*dt)
							+3.0*(ax_in(i-1,j,k-1,Nyx::AxRe)+kr_three(i-1,j,k-1)*dt)
							+3.0*(ax_in(i,j+1,k+1,Nyx::AxRe)+kr_three(i,j+1,k+1)*dt)
							+3.0*(ax_in(i,j+1,k-1,Nyx::AxRe)+kr_three(i,j+1,k-1)*dt)
							+3.0*(ax_in(i,j-1,k+1,Nyx::AxRe)+kr_three(i,j-1,k+1)*dt)
							+3.0*(ax_in(i,j-1,k-1,Nyx::AxRe)+kr_three(i,j-1,k-1)*dt)
							+2.0*(ax_in(i+1,j+1,k+1,Nyx::AxRe)+kr_three(i+1,j+1,k+1)*dt)
							+2.0*(ax_in(i+1,j+1,k-1,Nyx::AxRe)+kr_three(i+1,j+1,k-1)*dt)
							+2.0*(ax_in(i+1,j-1,k+1,Nyx::AxRe)+kr_three(i+1,j-1,k+1)*dt)
							+2.0*(ax_in(i+1,j-1,k-1,Nyx::AxRe)+kr_three(i+1,j-1,k-1)*dt)
							+2.0*(ax_in(i-1,j+1,k+1,Nyx::AxRe)+kr_three(i-1,j+1,k+1)*dt)
							+2.0*(ax_in(i-1,j+1,k-1,Nyx::AxRe)+kr_three(i-1,j+1,k-1)*dt)
							+2.0*(ax_in(i-1,j-1,k+1,Nyx::AxRe)+kr_three(i-1,j-1,k+1)*dt)
							+2.0*(ax_in(i-1,j-1,k-1,Nyx::AxRe)+kr_three(i-1,j-1,k-1)*dt)
							-88.0*(ax_in(i,j,k,Nyx::AxRe)+kr_three(i,j,k)*dt))
			       *invdeltasq_new/52.0
#ifdef GRAVITY
			       + phi(i,j,k)/hbaroverm*(ax_in(i,j,k,Nyx::AxRe)+kr_three(i,j,k)*dt)
#endif
			       ;

			       });

			     ParallelFor(bxfour,
 			   [=] AMREX_GPU_DEVICE (int i, int j, int k)
			   {
			     ax_out(i,j,k,Nyx::AxRe) = ax_in(i,j,k,Nyx::AxRe) + dt*(kr_one(i,j,k)+2.0*kr_two(i,j,k)+2.0*kr_three(i,j,k)+kr_four(i,j,k))/6.0;			     
			     ax_out(i,j,k,Nyx::AxIm) = ax_in(i,j,k,Nyx::AxIm) + dt*(ki_one(i,j,k)+2.0*ki_two(i,j,k)+2.0*ki_three(i,j,k)+ki_four(i,j,k))/6.0;			     	
			     ax_out(i,j,k,Nyx::AxDens) = ax_out(i,j,k,Nyx::AxRe)*ax_out(i,j,k,Nyx::AxRe)+ax_out(i,j,k,Nyx::AxIm)*ax_out(i,j,k,Nyx::AxIm);
			     ax_out(i,j,k,Nyx::AxPhas) = std::atan2(ax_out(i,j,k,Nyx::AxIm),ax_out(i,j,k,Nyx::AxRe));
			   });
      }

    Ax_new.FillBoundary(geom.periodicity());

#ifdef DEBUG
    if (Ax_new.contains_nan(0, Ax_new.nComp(), 0))
      {
        for (int i = 0; i < Ax_new.nComp(); i++)
          {
            if (Ax_new.contains_nan(i,1,0))
              {
		std::cout << "Testing component i for NaNs: " << i << std::endl;
		amrex::Abort("Ax_new has NaNs in this component::advance_FDM_FD()");
              }
          }
      }
#endif
}
#endif
