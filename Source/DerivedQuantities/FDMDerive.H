#ifndef _FDMDerive_H_
#define _FDMDerive_H_
#include <AMReX_BLFort.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef FDM
  void derfdmlohner
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real /*time*/, const int* /*bcrec*/, int level);

  void derfdmepot
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real /*time*/, const int* /*bcrec*/, int level);

  void derfdmekin
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real /*time*/, const int* /*bcrec*/, int level);

  void derfdmekinrho
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real /*time*/, const int* /*bcrec*/, int level);

  void derfdmekinv
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real /*time*/, const int* /*bcrec*/, int level);

  void derfdmangmomx
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real /*time*/, const int* /*bcrec*/, int level);

  void derfdmangmomy
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real /*time*/, const int* /*bcrec*/, int level);

  void derfdmangmomz
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real /*time*/, const int* /*bcrec*/, int level);
  
#endif

#ifdef __cplusplus
}
#endif

#endif
