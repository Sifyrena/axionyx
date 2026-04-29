#include <AMReX_REAL.H>

#include <Nyx.H>

using namespace amrex;

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef FDM
    void derfdmlohner(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
		      const FArrayBox& datfab, const Geometry& geomdata,
		      Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      // Here dat contains (AxDens, AxRe, AxIm)

      //DO TO: fix these parameters!
      Real eps = 0.2;
      Real mindens = Nyx::meandens;

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

	Real t,b;
	der(i,j,k,0) = 0.0;

	if ( (dat(i,j,k,0)>mindens) && (dat(i,j,k,1)!=0.0) ){
	  
	  t =             pow(dat(i+1,j  ,k  ,1)-dat(i  ,j  ,k  ,1)-dat(i  ,j  ,k  ,1)+dat(i-1,j  ,k  ,1),2) + 
	                  pow(dat(i  ,j+1,k  ,1)-dat(i  ,j  ,k  ,1)-dat(i  ,j  ,k  ,1)+dat(i  ,j-1,k  ,1),2) + 
	                  pow(dat(i  ,j  ,k+1,1)-dat(i  ,j  ,k  ,1)-dat(i  ,j  ,k  ,1)+dat(i  ,j  ,k-1,1),2) + 
	            0.125*pow(dat(i+1,j+1,k  ,1)-dat(i+1,j-1,k  ,1)-dat(i-1,j+1,k  ,1)+dat(i-1,j-1,k  ,1),2) + 
	            0.125*pow(dat(i+1,j  ,k+1,1)-dat(i+1,j  ,k-1,1)-dat(i-1,j  ,k+1,1)+dat(i-1,j  ,k-1,1),2) + 
	            0.125*pow(dat(i  ,j+1,k+1,1)-dat(i  ,j+1,k-1,1)-dat(i  ,j-1,k+1,1)+dat(i  ,j-1,k-1,1),2);
	    
	   b =            pow(std::abs(dat(i+1,j  ,k  ,1) -         dat(i  ,j  ,k  ,1))+std::abs(dat(i  ,j  ,k  ,1) -         dat(i-1,j  ,k  ,1))     + 
			 eps*(std::abs(dat(i+1,j  ,k  ,1))+std::abs(dat(i  ,j  ,k  ,1))+std::abs(dat(i  ,j  ,k  ,1))+std::abs(dat(i-1,j  ,k  ,1))),2) + 
	                  pow(std::abs(dat(i  ,j+1,k  ,1) -         dat(i  ,j  ,k  ,1))+std::abs(dat(i  ,j  ,k  ,1) -         dat(i  ,j-1,k  ,1))     + 
			 eps*(std::abs(dat(i  ,j+1,k  ,1))+std::abs(dat(i  ,j  ,k  ,1))+std::abs(dat(i  ,j  ,k  ,1))+std::abs(dat(i  ,j-1,k  ,1))),2) + 
	                  pow(std::abs(dat(i  ,j  ,k+1,1) -         dat(i  ,j  ,k  ,1))+std::abs(dat(i  ,j  ,k  ,1) -         dat(i  ,j  ,k-1,1))     + 
	                 eps*(std::abs(dat(i  ,j  ,k+1,1))+std::abs(dat(i  ,j  ,k  ,1))+std::abs(dat(i  ,j  ,k  ,1))+std::abs(dat(i  ,j  ,k-1,1))),2) + 
	             pow(0.5*(std::abs(dat(i+1,j+1,k  ,1) -         dat(i-1,j+1,k  ,1))+std::abs(dat(i+1,j-1,k  ,1) -         dat(i-1,j-1,k  ,1)))    + 
	                 eps*(std::abs(dat(i+1,j+1,k  ,1))+std::abs(dat(i-1,j+1,k  ,1))+std::abs(dat(i+1,j-1,k  ,1))+std::abs(dat(i-1,j-1,k  ,1))),2) + 
	             pow(0.5*(std::abs(dat(i+1,j+1,k  ,1) -         dat(i+1,j-1,k  ,1))+std::abs(dat(i-1,j+1,k  ,1) -         dat(i-1,j-1,k  ,1)))    + 
	                 eps*(std::abs(dat(i+1,j+1,k  ,1))+std::abs(dat(i+1,j-1,k  ,1))+std::abs(dat(i-1,j+1,k  ,1))+std::abs(dat(i-1,j-1,k  ,1))),2) + 
	             pow(0.5*(std::abs(dat(i+1,j  ,k+1,1) -         dat(i-1,j  ,k+1,1))+std::abs(dat(i+1,j  ,k-1,1) -         dat(i-1,j  ,k-1,1)))    + 
	                 eps*(std::abs(dat(i+1,j  ,k+1,1))+std::abs(dat(i-1,j  ,k+1,1))+std::abs(dat(i+1,j  ,k-1,1))+std::abs(dat(i-1,j  ,k-1,1))),2) + 
	             pow(0.5*(std::abs(dat(i+1,j  ,k+1,1) -         dat(i+1,j  ,k-1,1))+std::abs(dat(i-1,j  ,k+1,1) -         dat(i-1,j  ,k-1,1)))    + 
	                 eps*(std::abs(dat(i+1,j  ,k+1,1))+std::abs(dat(i+1,j  ,k-1,1))+std::abs(dat(i-1,j  ,k+1,1))+std::abs(dat(i-1,j  ,k-1,1))),2) + 
	             pow(0.5*(std::abs(dat(i  ,j+1,k+1,1) -         dat(i  ,j-1,k+1,1))+std::abs(dat(i  ,j+1,k-1,1) -         dat(i  ,j-1,k-1,1)))    + 
	                 eps*(std::abs(dat(i  ,j+1,k+1,1))+std::abs(dat(i  ,j-1,k+1,1))+std::abs(dat(i  ,j+1,k-1,1))+std::abs(dat(i  ,j-1,k-1,1))),2) + 
	             pow(0.5*(std::abs(dat(i  ,j+1,k+1,1) -         dat(i  ,j+1,k-1,1))+std::abs(dat(i  ,j-1,k+1,1) -         dat(i  ,j-1,k-1,1)))    + 
	                 eps*(std::abs(dat(i  ,j+1,k+1,1))+std::abs(dat(i  ,j+1,k-1,1))+std::abs(dat(i  ,j-1,k+1,1))+std::abs(dat(i  ,j-1,k-1,1))),2);
	    
	    der(i,j,k,0) = sqrt(t/b);
	    
	      }
	if ( (dat(i,j,k,0)>mindens) && (dat(i,j,k,2)!=0.0) ){

	  t =           pow(dat(i+1,j  ,k  ,2)-dat(i  ,j  ,k  ,2)-dat(i  ,j  ,k  ,2)+dat(i-1,j  ,k  ,2),2) + 
	                pow(dat(i  ,j+1,k  ,2)-dat(i  ,j  ,k  ,2)-dat(i  ,j  ,k  ,2)+dat(i  ,j-1,k  ,2),2) + 
	                pow(dat(i  ,j  ,k+1,2)-dat(i  ,j  ,k  ,2)-dat(i  ,j  ,k  ,2)+dat(i  ,j  ,k-1,2),2) + 
	          0.125*pow(dat(i+1,j+1,k  ,2)-dat(i+1,j-1,k  ,2)-dat(i-1,j+1,k  ,2)+dat(i-1,j-1,k  ,2),2) + 
	          0.125*pow(dat(i+1,j  ,k+1,2)-dat(i+1,j  ,k-1,2)-dat(i-1,j  ,k+1,2)+dat(i-1,j  ,k-1,2),2) + 
	          0.125*pow(dat(i  ,j+1,k+1,2)-dat(i  ,j+1,k-1,2)-dat(i  ,j-1,k+1,2)+dat(i  ,j-1,k-1,2),2);
	    
	  b =           pow(std::abs(dat(i+1,j  ,k  ,2) -         dat(i  ,j  ,k  ,2))+std::abs(dat(i  ,j  ,k  ,2) -         dat(i-1,j  ,k  ,2))     + 
		       eps*(std::abs(dat(i+1,j  ,k  ,2))+std::abs(dat(i  ,j  ,k  ,2))+std::abs(dat(i  ,j  ,k  ,2))+std::abs(dat(i-1,j  ,k  ,2))),2) + 
	                pow(std::abs(dat(i  ,j+1,k  ,2) -         dat(i  ,j  ,k  ,2))+std::abs(dat(i  ,j  ,k  ,2) -         dat(i  ,j-1,k  ,2))     + 
	               eps*(std::abs(dat(i  ,j+1,k  ,2))+std::abs(dat(i  ,j  ,k  ,2))+std::abs(dat(i  ,j  ,k  ,2))+std::abs(dat(i  ,j-1,k  ,2))),2) + 
	                pow(std::abs(dat(i  ,j  ,k+1,2) -         dat(i  ,j  ,k  ,2))+std::abs(dat(i  ,j  ,k  ,2) -         dat(i  ,j  ,k-1,2))     + 
	               eps*(std::abs(dat(i  ,j  ,k+1,2))+std::abs(dat(i  ,j  ,k  ,2))+std::abs(dat(i  ,j  ,k  ,2))+std::abs(dat(i  ,j  ,k-1,2))),2) + 
	           pow(0.5*(std::abs(dat(i+1,j+1,k  ,2) -         dat(i-1,j+1,k  ,2))+std::abs(dat(i+1,j-1,k  ,2) -         dat(i-1,j-1,k  ,2)))    + 
	               eps*(std::abs(dat(i+1,j+1,k  ,2))+std::abs(dat(i-1,j+1,k  ,2))+std::abs(dat(i+1,j-1,k  ,2))+std::abs(dat(i-1,j-1,k  ,2))),2) + 
	           pow(0.5*(std::abs(dat(i+1,j+1,k  ,2) -         dat(i+1,j-1,k  ,2))+std::abs(dat(i-1,j+1,k  ,2) -         dat(i-1,j-1,k  ,2)))    + 
		       eps*(std::abs(dat(i+1,j+1,k  ,2))+std::abs(dat(i+1,j-1,k  ,2))+std::abs(dat(i-1,j+1,k  ,2))+std::abs(dat(i-1,j-1,k  ,2))),2) + 
	           pow(0.5*(std::abs(dat(i+1,j  ,k+1,2) -         dat(i-1,j  ,k+1,2))+std::abs(dat(i+1,j  ,k-1,2) -         dat(i-1,j  ,k-1,2)))    + 
		       eps*(std::abs(dat(i+1,j  ,k+1,2))+std::abs(dat(i-1,j  ,k+1,2))+std::abs(dat(i+1,j  ,k-1,2))+std::abs(dat(i-1,j  ,k-1,2))),2) + 
	           pow(0.5*(std::abs(dat(i+1,j  ,k+1,2) -         dat(i+1,j  ,k-1,2))+std::abs(dat(i-1,j  ,k+1,2) -         dat(i-1,j  ,k-1,2)))    + 
		       eps*(std::abs(dat(i+1,j  ,k+1,2))+std::abs(dat(i+1,j  ,k-1,2))+std::abs(dat(i-1,j  ,k+1,2))+std::abs(dat(i-1,j  ,k-1,2))),2) + 
	           pow(0.5*(std::abs(dat(i  ,j+1,k+1,2) -         dat(i  ,j-1,k+1,2))+std::abs(dat(i  ,j+1,k-1,2) -         dat(i  ,j-1,k-1,2)))    + 
		       eps*(std::abs(dat(i  ,j+1,k+1,2))+std::abs(dat(i  ,j-1,k+1,2))+std::abs(dat(i  ,j+1,k-1,2))+std::abs(dat(i  ,j-1,k-1,2))),2) + 
	           pow(0.5*(std::abs(dat(i  ,j+1,k+1,2) -         dat(i  ,j+1,k-1,2))+std::abs(dat(i  ,j-1,k+1,2) -         dat(i  ,j-1,k-1,2)))    + 
		       eps*(std::abs(dat(i  ,j+1,k+1,2))+std::abs(dat(i  ,j+1,k-1,2))+std::abs(dat(i  ,j-1,k+1,2))+std::abs(dat(i  ,j-1,k-1,2))),2);
	    
	    der(i,j,k,0) = max(der(i,j,k,0),sqrt(t/b));
	    }


      });
    }

    void derfdmepot(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
		    const FArrayBox& datfab, const Geometry& geomdata,
		    Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {
      
      auto const dat = datfab.array();
      auto const der = derfab.array();

      // Here dat contains (AxDens, Phi_Grav)
      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
#ifdef GRAVITY
	der(i,j,k,0) = -dat(i,j,k,0)*dat(i,j,k,1)/2.0;
#else
	der(i,j,k,0) = 0.0;
#endif
      });
    }

    void derfdmekin(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                    const FArrayBox& datfab, const Geometry& geomdata,
                    Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {
        auto const dat = datfab.array();
        auto const der = derfab.array();
    
        auto const dx = geomdata.CellSizeArray();
    
        const Real inv12  = 1.0_rt / 12.0_rt;
        const Real prefac = 0.5_rt * Nyx::hbaroverm * Nyx::hbaroverm;
    
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            auto sq = [] AMREX_GPU_DEVICE (Real x) -> Real { return x * x; };
    
            Real dAxRe_dx =
                (-dat(i+2,j,k,0) + 8.0_rt*dat(i+1,j,k,0)
                 -8.0_rt*dat(i-1,j,k,0) + dat(i-2,j,k,0))
                * inv12 / dx[0];
    
            Real dAxIm_dx =
                (-dat(i+2,j,k,1) + 8.0_rt*dat(i+1,j,k,1)
                 -8.0_rt*dat(i-1,j,k,1) + dat(i-2,j,k,1))
                * inv12 / dx[0];
    
            Real dAxRe_dy =
                (-dat(i,j+2,k,0) + 8.0_rt*dat(i,j+1,k,0)
                 -8.0_rt*dat(i,j-1,k,0) + dat(i,j-2,k,0))
                * inv12 / dx[1];
    
            Real dAxIm_dy =
                (-dat(i,j+2,k,1) + 8.0_rt*dat(i,j+1,k,1)
                 -8.0_rt*dat(i,j-1,k,1) + dat(i,j-2,k,1))
                * inv12 / dx[1];
    
            Real dAxRe_dz =
                (-dat(i,j,k+2,0) + 8.0_rt*dat(i,j,k+1,0)
                 -8.0_rt*dat(i,j,k-1,0) + dat(i,j,k-2,0))
                * inv12 / dx[2];
    
            Real dAxIm_dz =
                (-dat(i,j,k+2,1) + 8.0_rt*dat(i,j,k+1,1)
                 -8.0_rt*dat(i,j,k-1,1) + dat(i,j,k-2,1))
                * inv12 / dx[2];
    
            der(i,j,k,dcomp) = prefac * (
                  sq(dAxRe_dx) + sq(dAxIm_dx)
                + sq(dAxRe_dy) + sq(dAxIm_dy)
                + sq(dAxRe_dz) + sq(dAxIm_dz)
            );
        });
    }

          
    void derfdmekinrho(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
		       const FArrayBox& datfab, const Geometry& geomdata,
		       Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {
      
      auto const dat = datfab.array();
      auto const der = derfab.array();

      auto const dx = geomdata.CellSizeArray();

      // Here dat contains (AxDens)

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
		der(i,j,k,0) = Nyx::hbaroverm*Nyx::hbaroverm / 2.0 *
			( pow((-std::sqrt(std::abs(dat(i+2,j,k,0)))+8.0*std::sqrt(std::abs(dat(i+1,j,k,0)))-8.0*std::sqrt(std::abs(dat(i-1,j,k,0)))+std::sqrt(std::abs(dat(i-2,j,k,0))))/(12.0*dx[0]),2)
		     +pow((-std::sqrt(std::abs(dat(i,j+2,k,0)))+8.0*std::sqrt(std::abs(dat(i,j+1,k,0)))-8.0*std::sqrt(std::abs(dat(i,j-1,k,0)))+std::sqrt(std::abs(dat(i,j-2,k,0))))/(12.0*dx[1]),2)
		     +pow((-std::sqrt(std::abs(dat(i,j,k+2,0)))+8.0*std::sqrt(std::abs(dat(i,j,k+1,0)))-8.0*std::sqrt(std::abs(dat(i,j,k-1,0)))+std::sqrt(std::abs(dat(i,j,k-2,0))))/(12.0*dx[2]),2) );
      });
    }

    void derfdmekinv(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
		     const FArrayBox& datfab, const Geometry& geomdata,
		     Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {
      
      auto const dat = datfab.array();
      auto const der = derfab.array();

      auto const dx = geomdata.CellSizeArray();
      Real diff[]={0.0,0.0,0.0};

      // Here dat contains (AxDens, AxPhas)

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
		//Have to consider the cases where the phase jumps from close to -pi to pi and vice versa
		Real diff0 = min( std::abs(dat(i+1,j,k,1)-dat(i-1,j,k,1)) , std::abs(dat(i+1,j,k,1)-dat(i-1,j,k,1)-2.0*M_PI) , std::abs(dat(i+1,j,k,1)-dat(i-1,j,k,1)+2.0*M_PI) );
		Real diff1 = min( std::abs(dat(i,j+1,k,1)-dat(i,j-1,k,1)) , std::abs(dat(i,j+1,k,1)-dat(i,j-1,k,1)-2.0*M_PI) , std::abs(dat(i,j+1,k,1)-dat(i,j-1,k,1)+2.0*M_PI) );
		Real diff2 = min( std::abs(dat(i,j,k+1,1)-dat(i,j,k-1,1)) , std::abs(dat(i,j,k+1,1)-dat(i,j,k-1,1)-2.0*M_PI) , std::abs(dat(i,j,k+1,1)-dat(i,j,k-1,1)+2.0*M_PI) );

		der(i,j,k,0) = Nyx::hbaroverm*Nyx::hbaroverm / 2.0 * dat(i,j,k,0) *
		( pow( diff0/(2.0*dx[0]),2 )
		 +pow( diff1/(2.0*dx[1]),2 )
		 +pow( diff2/(2.0*dx[2]),2 ));
      });
    }

    void derfdmangmomx(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
		       const FArrayBox& datfab, const Geometry& geomdata,
		       Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {
      
      auto const dat = datfab.array();
      auto const der = derfab.array();

      auto const dx = geomdata.CellSizeArray();
      auto const plo = geomdata.ProbLoArray();

      // Here dat contains (AxRe, AxIm)

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
	der(i,j,k,0) = Nyx::hbaroverm*( -dat(i,j,k,0)*(dat(i,j,k+1,1)-dat(i,j,k-1,1))/(2.0*dx[2])*((j+0.5)*dx[1]+plo[1])
					+dat(i,j,k,1)*(dat(i,j,k+1,0)-dat(i,j,k-1,0))/(2.0*dx[2])*((j+0.5)*dx[1]+plo[1])
					+dat(i,j,k,0)*(dat(i,j+1,k,1)-dat(i,j-1,k,1))/(2.0*dx[1])*((k+0.5)*dx[2]+plo[2])
					-dat(i,j,k,1)*(dat(i,j+1,k,0)-dat(i,j-1,k,0))/(2.0*dx[1])*((k+0.5)*dx[2]+plo[2]) );
      });
    }

    void derfdmangmomy(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
		       const FArrayBox& datfab, const Geometry& geomdata,
		       Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {
      
      auto const dat = datfab.array();
      auto const der = derfab.array();

      auto const dx = geomdata.CellSizeArray();
      auto const plo = geomdata.ProbLoArray();

      // Here dat contains (AxRe, AxIm)

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
	der(i,j,k,0) = Nyx::hbaroverm*( -dat(i,j,k,0)*(dat(i+1,j,k,1)-dat(i-1,j,k,1))/(2.0*dx[0])*((k+0.5)*dx[2]+plo[2])
					+dat(i,j,k,1)*(dat(i+1,j,k,0)-dat(i-1,j,k,0))/(2.0*dx[0])*((k+0.5)*dx[2]+plo[2])
					+dat(i,j,k,0)*(dat(i,j,k+1,1)-dat(i,j,k-1,1))/(2.0*dx[2])*((i+0.5)*dx[0]+plo[0])
					-dat(i,j,k,1)*(dat(i,j,k+1,0)-dat(i,j,k-1,0))/(2.0*dx[2])*((i+0.5)*dx[0]+plo[0]) );

      });
    }

    void derfdmangmomz(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
		       const FArrayBox& datfab, const Geometry& geomdata,
		       Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {
      
      auto const dat = datfab.array();
      auto const der = derfab.array();

      auto const dx = geomdata.CellSizeArray();
      auto const plo = geomdata.ProbLoArray();

      // Here dat contains (AxRe, AxIm)

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
	der(i,j,k,0) = Nyx::hbaroverm*( -dat(i,j,k,0)*(dat(i,j+1,k,1)-dat(i,j-1,k,1))/(2.0*dx[1])*((i+0.5)*dx[0]+plo[0])
					+dat(i,j,k,1)*(dat(i,j+1,k,0)-dat(i,j-1,k,0))/(2.0*dx[1])*((i+0.5)*dx[0]+plo[0])
					+dat(i,j,k,0)*(dat(i+1,j,k,1)-dat(i-1,j,k,1))/(2.0*dx[0])*((j+0.5)*dx[1]+plo[1])
					-dat(i,j,k,1)*(dat(i+1,j,k,0)-dat(i-1,j,k,0))/(2.0*dx[0])*((j+0.5)*dx[1]+plo[1]) );
      });
    }
#endif

#ifdef __cplusplus
}
#endif
