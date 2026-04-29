#include "Nyx.H"
#include "Prob.H"

using namespace amrex;

#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <fstream>
#include <vector>

void prob_param_special_fill(GpuArray<Real,max_prob_param>& prob_param)
{}

void prob_errtags_default(Vector<AMRErrorTag>& errtags)
{
  AMRErrorTagInfo info;
  errtags.push_back(AMRErrorTag(1,AMRErrorTag::GREATER,"AxDens",info));
}

#include "Nyx.H"
#include "Prob.H"

using namespace amrex;

AMREX_GPU_DEVICE AMREX_FORCE_INLINE

// The following function is modified to directly read from an external H5 file.
void prob_initdata_state(const int i,
                         const int j,
                         const int k,
                         Array4<Real> const& state,
                         Array4<Real> const& axion,
                         GeometryData const& geomdata,
                         const GpuArray<Real,max_prob_param>& prob_param)
{
  Real hubl = Nyx::comoving_h;
  Real meandens = 2.775e+11*pow(hubl,2)*Nyx::comoving_OmM;

  amrex::Real h_species = prob_param[h_species_comp];

  const Real* prob_lo = geomdata.ProbLo();
  const Real* prob_hi = geomdata.ProbHi();
  const Real* dx = geomdata.CellSize();

  // Centres
  Real posx = (prob_hi[0]+prob_lo[0])/2.0;
  Real posy = (prob_hi[1]+prob_lo[1])/2.0;
  Real posz = (prob_hi[2]+prob_lo[2])/2.0;

}

void prob_initdata(const int i,
                   const int j,
                   const int k,
                   amrex::Array4<amrex::Real> const& state,
                   amrex::Array4<amrex::Real> const& diag_eos,
                   amrex::Array4<amrex::Real> const& axion,
                   amrex::GeometryData const& geomdata,
                   const amrex::GpuArray<amrex::Real,max_prob_param>& prob_param)
{

  ;//prob_initdata_state(i, j ,k, state, axion, geomdata, prob_param);

}

void prob_initdata_on_box(
    const Box& bx,
    Array4<amrex::Real> const& axion,
    GeometryData const& geomdata,
    const GpuArray<Real, max_prob_param>& prob_param)
{
   ; // FW NOTICE: GO!
  }
  
void prob_initdata_FDM_on_box(
    const Box& bx,
    Array4<amrex::Real> const& axion,
    GeometryData const& geomdata,
    const GpuArray<Real, max_prob_param>& prob_param,
    const amrex::Vector<double>& axionRebuffer,
    const amrex::Vector<double>& axionImbuffer)
{
    // FW NOTICE: GO!
    /*
   amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Box& domain = geomdata.Domain();

    // Get the number of cells in each dimension for the full domain
    const int nx = domain.length(0); 

    
      int idx = (j * (nx) * (nx) + i * (nx) + k);
      // Initialise Axion Wavefunction
      //axion(i,j,k, Nyx::AxRe) = axionRebuffer[idx];
      //axion(i,j,k, Nyx::AxIm) = axionImbuffer[idx];

      axion(i,j,k,Nyx::AxPhas) = 0;//std::atan2(axion(i,j,k,Nyx::AxIm),axion(i,j,k,Nyx::AxRe));
      
      axion(i,j,k,Nyx::AxDens) = axion(i,j,k,Nyx::AxRe)*axion(i,j,k,Nyx::AxRe) + axion(i,j,k,Nyx::AxIm)*axion(i,j,k,Nyx::AxIm);

    }); 
    */
}

void prob_initdata_state_on_box(const Box& bx,
                                Array4<amrex::Real> const& axion,
                                GeometryData const& geomdata,
                                const GpuArray<Real,max_prob_param>& prob_param)
{
 //load_hdf5(state, bx);
}
