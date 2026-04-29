#include <AMReX_REAL.H>

#include <Nyx.H>
#include "compute_laplace_operator.H"

#ifdef AXCOMPLEX
#include <LLUTs.h>
#endif

using namespace amrex;

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef AXCOMPLEX
  inline void stringcorre( amrex::Real *da, amrex::Real &re);

  void derstring(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
		 const FArrayBox& datfab, const Geometry& geomdata,
		 Real time, const int* /*bcrec*/, int level)
  {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      // Here dat contains (AxRe,AxIm)

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
	der(i,j,k,0) = 0.0;
	int hand = 0;

	if ( dat(i,j,k,1)*dat(i+1,j,k,1) < 0)
	  if( (dat(i,j,k,1)*dat(i+1,j,k,0)-dat(i,j,k,0)*dat(i+1,j,k,1)) > 0)
	    hand++;
	  else
	    hand--;

	if ( dat(i+1,j,k,1)*dat(i+1,j+1,k,1) < 0 )
	  if( (dat(i+1,j,k,1)*dat(i+1,j+1,k,0)-dat(i+1,j,k,0)*dat(i+1,j+1,k,1)) > 0)
	    hand++;
	  else
	    hand--;
	
	if ( dat(i+1,j+1,k,1)*dat(i,j+1,k,1) < 0 )
	  if( (dat(i+1,j+1,k,1)*dat(i,j+1,k,0)-dat(i+1,j+1,k,0)*dat(i,j+1,k,1)) > 0)
	    hand++;
	  else
	    hand--;
	
	if ( dat(i,j+1,k,1)*dat(i,j,k,1) < 0 )
	  if( (dat(i,j+1,k,1)*dat(i,j,k,0)-dat(i,j+1,k,0)*dat(i,j,k,1)) > 0)
	    hand++;
	  else
	    hand--;
	
	if ( (hand == 2) || (hand == -2) )
	  der(i,j,k,0) += 1.0;
	hand = 0;

	  if ( dat(i,j,k,1)*dat(i+1,j,k,1) < 0 )
	    if( (dat(i,j,k,1)*dat(i+1,j,k,0)-dat(i,j,k,0)*dat(i+1,j,k,1)) > 0)
	      hand++;
	    else
	      hand--;

	if ( dat(i+1,j,k,1)*dat(i+1,j,k+1,1) < 0 )
	  if( (dat(i+1,j,k,1)*dat(i+1,j,k+1,0)-dat(i+1,j,k,0)*dat(i+1,j,k+1,1)) > 0)
	    hand++;
	  else
	    hand--;

	if ( dat(i+1,j,k+1,1)*dat(i,j,k+1,1) < 0 )
	  if( (dat(i+1,j,k+1,1)*dat(i,j,k+1,0)-dat(i+1,j,k+1,0)*dat(i,j,k+1,1)) > 0)
	    hand++;
	  else
	    hand--;

	if ( dat(i,j,k+1,1)*dat(i,j,k,1) < 0 )
	  if( (dat(i,j,k+1,1)*dat(i,j,k,0)-dat(i,j,k+1,0)*dat(i,j,k,1)) > 0)
	    hand++;
	  else
	    hand--;
	
	if ( (hand == 2) || (hand == -2) )
	  der(i,j,k,0) += 1.0;
	hand = 0;

	if ( dat(i,j,k,1)*dat(i,j+1,k,1) < 0 )
	  if( (dat(i,j,k,1)*dat(i,j+1,k,0)-dat(i,j,k,0)*dat(i,j+1,k,1)) > 0)
	    hand++;
	  else
	    hand--;

	if ( dat(i,j+1,k,1)*dat(i,j+1,k+1,1) < 0 )
	  if( (dat(i,j+1,k,1)*dat(i,j+1,k+1,0)-dat(i,j+1,k,0)*dat(i,j+1,k+1,1)) > 0)
	    hand++;
	  else
	    hand--;

	if ( dat(i,j+1,k+1,1)*dat(i,j,k+1,1) < 0 )
	  if( (dat(i,j+1,k+1,1)*dat(i,j,k+1,0)-dat(i,j+1,k+1,0)*dat(i,j,k+1,1)) > 0)
	    hand++;
	  else
	    hand--;

	if ( dat(i,j,k+1,1)*dat(i,j,k,1) < 0 )
	  if( (dat(i,j,k+1,1)*dat(i,j,k,0)-dat(i,j,k+1,0)*dat(i,j,k,1)) > 0)
	    hand++;
	  else
	    hand--;

	if ( (hand == 2) || (hand == -2) )
	  der(i,j,k,0) += 1.0;
 	hand = 0;

      });
    }


    void deraxvel(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
		  const FArrayBox& datfab, const Geometry& geomdata,
		  Real time, const int* /*bcrec*/, int level)
    {
     
      auto const dat = datfab.array();
      auto const der = derfab.array();

      // Here dat contains (AxRe,AxIm,AxvRe,AxvIm)

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
	
	Real vA;
	Real da[10] = {0.0};
	bool skipLUT = false;
	bool LUTcorr = true;
	Real pre = 0.5*Nyx::sigmalut;
	da[0] = dat(i,j,k,Nyx::AxRe);
	da[1] = dat(i,j,k,Nyx::AxIm);
	da[2] = dat(i,j,k,Nyx::AxvRe);
	da[3] = dat(i,j,k,Nyx::AxvIm);
	if (LUTcorr) {
	  da[4] = pre*(dat(i+1,j,k,Nyx::AxRe)-dat(i-1,j,k,Nyx::AxRe));
	  if (std::abs(da[0]) > std::abs(5*da[4]))
	    { skipLUT = true; goto jmp;}
	  da[7] = pre*(dat(i+1,j,k,Nyx::AxIm)-dat(i-1,j,k,Nyx::AxIm));
	  if (std::abs(da[1]) > std::abs(5*da[7]))
	    { skipLUT = true; goto jmp;}
	  da[5] = pre*(dat(i,j+1,k,Nyx::AxRe)-dat(i,j-1,k,Nyx::AxRe));
	  if (std::abs(da[0]) > std::abs(5*da[5]))
	    { skipLUT = true; goto jmp;}
	  da[8] = pre*(dat(i,j+1,k,Nyx::AxIm)-dat(i,j-1,k,Nyx::AxIm));
	  if (std::abs(da[1]) > std::abs(5*da[8]))
	    { skipLUT = true; goto jmp;}
	  da[6] = pre*(dat(i,j,k+1,Nyx::AxRe)-dat(i,j,k-1,Nyx::AxRe));
	  if (std::abs(da[0]) > std::abs(5*da[6]))
	    { skipLUT = true; goto jmp;}
	  da[9] = pre*(dat(i,j,k+1,Nyx::AxIm)-dat(i,j,k-1,Nyx::AxIm));
	  if (std::abs(da[1]) > std::abs(5*da[9]))
	    { skipLUT = true; goto jmp;}
	  stringcorre(da,vA);
	}
      jmp:
	if (!LUTcorr || skipLUT)
	  vA = (da[0]*da[3]-da[1]*da[2])/(da[0]*da[0]+da[1]*da[1]);

	der(i,j,k,0) = vA;

      });
    }


  inline void stringcorre( Real *da, Real &re)
  {
    Real ma  = std::sqrt(da[4]*da[4]+da[5]*da[5]+da[6]*da[6]);
    Real mb  = std::sqrt(da[7]*da[7]+da[8]*da[8]+da[9]*da[9]);
    bool flip = false;
    if (ma < mb)
      flip=true;
    Real e   = flip ? ma/mb : mb/ma;
    Real cb  = (da[4]*da[7]+da[5]*da[8]+da[6]*da[9])/(ma*mb);
    Real sb  = std::sqrt(1-cb*cb);
    Real A   = flip ? da[1]/mb : da[0]/ma;
    Real B   = flip ? da[0]/mb : da[1]/ma;
    Real s1  = 1;
    Real s2  = 1;
    if (A < 0){ A = -A; B = -B; s1 = -s1; s2 = -s2;};
    if (B < 0){ B = -B; cb = -cb; s2 = -s2;};
    if ( (A <= 5) && (B <= 5) && ( e <= 1 ) && (cb <= 1) && (cb >= -1))
      {
	size_t iA = A*(N_A-1)/5;
	size_t iB = B*(N_B-1)/5;
	size_t iE = e*(N_E-1);
	size_t iC = (cb+1)*(N_C-1)/2;
	size_t i0 = iC + N_C*iE + N_CE*iB + N_CEB*iA;
	Real dA  = A*(N_A-1)/5-iA;
	Real dB  = B*(N_B-1)/5-iB;
	Real dE  = e*(N_E-1)-iE;
	Real dC  = (cb+1)*(N_C-1)/2-iC;
	Real dA1  = 1-dA;
	Real dB1  = 1-dB;
	Real dE1  = 1-dE;
	Real dC1  = 1-dC;
	Real L000 = (Real) (LUTc[i0]*dC1 + LUTc[i0+1]*dC);
	Real L001 = (Real) (LUTc[i0+N_C]*dC1 + LUTc[i0+N_C+1]*dC);
	Real L010 = (Real) (LUTc[i0+N_E]*dC1 + LUTc[i0+N_E+1]*dC);
	Real L011 = (Real) (LUTc[i0+N_E+N_C]*dC1 + LUTc[i0+N_E+N_C+1]*dC);
	Real L100 = (Real) (LUTc[i0+N_B]*dC1 + LUTc[i0+N_B+1]*dC);
	Real L101 = (Real) (LUTc[i0+N_B+N_C]*dC1 + LUTc[i0+N_B+N_C+1]*dC);
	Real L110 = (Real) (LUTc[i0+N_B+N_E]*dC1 + LUTc[i0+N_B+N_E+1]*dC);
	Real L111 = (Real) (LUTc[i0+N_B+N_E+N_C]*dC1 + LUTc[i0+N_B+N_E+N_C+1]*dC);
	L000 = L000*dE1 + L001*dE;
	L010 = L010*dE1 + L011*dE;
	L100 = L100*dE1 + L101*dE;
	L110 = L110*dE1 + L111*dE;
	L000 = L000*dB1 + L010*dB;
	L100 = L100*dB1 + L110*dB;
	Real Lc = s1*(L000*dA1 + L100*dA);
	L000 = (Real) (LUTs[i0]*dC1             + LUTs[i0+1]*dC);
	L001 = (Real) (LUTs[i0+N_C]*dC1         + LUTs[i0+N_C+1]*dC);
	L010 = (Real) (LUTs[i0+N_E]*dC1         + LUTs[i0+N_E+1]*dC);
	L011 = (Real) (LUTs[i0+N_E+N_C]*dC1     + LUTs[i0+N_E+N_C+1]*dC);
	L100 = (Real) (LUTs[i0+N_B]*dC1         + LUTs[i0+N_B+1]*dC);
	L101 = (Real) (LUTs[i0+N_B+N_C]*dC1     + LUTs[i0+N_B+N_C+1]*dC);
	L110 = (Real) (LUTs[i0+N_B+N_E]*dC1     + LUTs[i0+N_B+N_E+1]*dC);
	L111 = (Real) (LUTs[i0+N_B+N_E+N_C]*dC1 + LUTs[i0+N_B+N_E+N_C+1]*dC);
	L000 = L000*dE1 + L001*dE;
	L010 = L010*dE1 + L011*dE;
	L100 = L100*dE1 + L101*dE;
	L110 = L110*dE1 + L111*dE;
	L000 = L000*dB1 + L010*dB;
	L100 = L100*dB1 + L110*dB;
	Real Ls = s2*(L000*dA1 + L100*dA);
	re = flip ? ((Real) Ls*da[3]-Lc*da[2])/mb : ((Real) Lc*da[3]-Ls*da[2])/ma;
      }
    else
      re = (da[0]*da[3]-da[1]*da[2])/(da[0]*da[0]+da[1]*da[1]);
  }

    void deraxgrad(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
		   const FArrayBox& datfab, const Geometry& geomdata,
		   Real time, const int* /*bcrec*/, int level)
    {
     
      auto const dat = datfab.array();
      auto const der = derfab.array();

      // Here dat contains (AxRe,AxIm)

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
	Real invdeltasq = 1.0;
	Real axregrad= compute_laplace_operator(dat,i,j,k,Nyx::AxRe,invdeltasq,Nyx::neighbors);
	Real aximgrad= compute_laplace_operator(dat,i,j,k,Nyx::AxIm,invdeltasq,Nyx::neighbors);
	der(i,j,k,0) = std::max(axregrad,aximgrad);
      });
    }

    void deraxgraddens(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
		       const FArrayBox& datfab, const Geometry& geomdata,
		       Real time, const int* /*bcrec*/, int level)
    {
     
      auto const dat = datfab.array();
      auto const der = derfab.array();

      const Box& densbox = grow(bx,Nyx::neighbors);
      FArrayBox densfab(densbox, 1);
      auto const dens = densfab.array();

      auto const dx = geomdata.CellSizeArray();
      const Real invdeltasq = std::pow(time,3)/dx[0];

      // Here dat contains (AxRe,AxIm)

      amrex::ParallelFor(densbox,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
	dens(i,j,k,0) = dat(i,j,k,Nyx::AxRe)*dat(i,j,k,Nyx::AxRe)+dat(i,j,k,Nyx::AxIm)*dat(i,j,k,Nyx::AxIm);
      });


      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
	der(i,j,k,0) = compute_laplace_operator(dens,i,j,k,0,invdeltasq,Nyx::neighbors);
      });

    }

    void deraxphas(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
		   const FArrayBox& datfab, const Geometry& geomdata,
		   Real time, const int* /*bcrec*/, int level)
    {
     
      auto const dat = datfab.array();
      auto const der = derfab.array();

      // Here dat contains (AxRe,AxIm)

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
	der(i,j,k,0) = std::atan2(dat(i,j,k,Nyx::AxIm),dat(i,j,k,Nyx::AxRe));
      });
    }

    void deraxdens(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
		   const FArrayBox& datfab, const Geometry& geomdata,
		   Real time, const int* /*bcrec*/, int level)
    {
     
      auto const dat = datfab.array();
      auto const der = derfab.array();

      // Here dat contains (AxRe,AxIm)

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
	der(i,j,k,0) = dat(i,j,k,Nyx::AxRe)*dat(i,j,k,Nyx::AxRe)+dat(i,j,k,Nyx::AxIm)*dat(i,j,k,Nyx::AxIm);
      });
    }

    void deraxener(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
		   const FArrayBox& datfab, const Geometry& geomdata,
		   Real time, const int* /*bcrec*/, int level)
    {
     
      // Here dat contains (AxRe,AxIm,AxvRe,AxvIm)

      auto const dat = datfab.array();
      auto const der = derfab.array();

      auto const dx = geomdata.CellSizeArray();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
	der(i,j,k,0) = 0.5*(pow(dat(i,j,k,Nyx::AxRe)*dat(i,j,k,Nyx::AxvIm)-dat(i,j,k,Nyx::AxIm)*dat(i,j,k,Nyx::AxvRe),2)
	  +0.25/dx[0]/dx[0]*pow(dat(i,j,k,Nyx::AxRe)*(dat(i+1,j,k,Nyx::AxIm)-dat(i-1,j,k,Nyx::AxIm))-dat(i,j,k,Nyx::AxIm)*(dat(i+1,j,k,Nyx::AxRe)-dat(i-1,j,k,Nyx::AxRe)),2)
	  +0.25/dx[1]/dx[1]*pow(dat(i,j,k,Nyx::AxRe)*(dat(i,j+1,k,Nyx::AxIm)-dat(i,j-1,k,Nyx::AxIm))-dat(i,j,k,Nyx::AxIm)*(dat(i,j+1,k,Nyx::AxRe)-dat(i,j-1,k,Nyx::AxRe)),2)
	  +0.25/dx[2]/dx[2]*pow(dat(i,j,k,Nyx::AxRe)*(dat(i,j,k+1,Nyx::AxIm)-dat(i,j,k-1,Nyx::AxIm))-dat(i,j,k,Nyx::AxIm)*(dat(i,j,k+1,Nyx::AxRe)-dat(i,j,k-1,Nyx::AxRe)),2))
	  /(dat(i,j,k,Nyx::AxRe)*dat(i,j,k,Nyx::AxRe)+dat(i,j,k,Nyx::AxIm)*dat(i,j,k,Nyx::AxIm))
	  -Nyx::string_massterm*pow(time/Nyx::string_time1,Nyx::string_n+2)/Nyx::string_time1/Nyx::string_time1*dat(i,j,k,Nyx::AxIm)/
	  std::sqrt(dat(i,j,k,Nyx::AxRe)*dat(i,j,k,Nyx::AxRe)+dat(i,j,k,Nyx::AxIm)*dat(i,j,k,Nyx::AxIm));
      });
    }


#endif

#ifdef __cplusplus
}
#endif
