#include "Nyx.H"
#include "Prob.H"

using namespace amrex;

// ABC: scaffolding for the next campaign (USE_HYDRO=TRUE, NO_HYDRO=FALSE).
// Every function below is a compliant but empty stub, matching the
// interface Prob.H declares under !defined(NO_HYDRO) with AXIONYX=TRUE.
// None of them initialise any state yet.

void prob_param_special_fill(GpuArray<Real,max_prob_param>& prob_param)
{}

void prob_errtags_default(Vector<AMRErrorTag>& errtags)
{}

void prob_initdata(const int i,
                    const int j,
                    const int k,
                    Array4<Real> const& state,
                    Array4<Real> const& diag_eos,
                    Array4<Real> const& axion,
                    GeometryData const& geomdata,
                    const GpuArray<Real,max_prob_param>& prob_param)
{}

void prob_initdata_state_on_box(const Box& bx,
                                 Array4<Real> const& axion,
                                 GeometryData const& geomdata,
                                 const GpuArray<Real,max_prob_param>& prob_param)
{}

void prob_initdata_state(const int i,
                          const int j,
                          const int k,
                          Array4<Real> const& state,
                          Array4<Real> const& axion,
                          GeometryData const& geomdata,
                          const GpuArray<Real,max_prob_param>& prob_param)
{}

void prob_initdata_on_box(const Box& bx,
                           Array4<Real> const& state,
                           Array4<Real> const& diag_eos,
                           Array4<Real> const& axion,
                           GeometryData const& geomdata,
                           const GpuArray<Real,max_prob_param>& prob_param)
{}

void prob_initdata_state_on_box(const Box& bx,
                                 Array4<Real> const& state,
                                 Array4<Real> const& diag_eos,
                                 Array4<Real> const& axion,
                                 GeometryData const& geomdata,
                                 const GpuArray<Real,max_prob_param>& prob_param)
{}

void prob_initdata_Basic_on_box(const Box& bx,
                                 Array4<Real> const& state,
                                 Array4<Real> const& diag_eos,
                                 Array4<Real> const& axion,
                                 GeometryData const& geomdata,
                                 const GpuArray<Real,max_prob_param>& prob_param)
{}

void prob_initdata_Fluid_on_box(const Box& bx,
                                 Array4<Real> const& state,
                                 Array4<Real> const& diag_eos,
                                 Array4<Real> const& axion,
                                 GeometryData const& geomdata,
                                 const GpuArray<Real,max_prob_param>& prob_param,
                                 const amrex::Vector<double>& axionRebuffer,
                                 const amrex::Vector<double>& axionImbuffer)
{}
