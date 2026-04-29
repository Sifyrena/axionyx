#ifndef _AxComplexDerive_H_
#define _AxComplexDerive_H_
#include <AMReX_BLFort.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef AXCOMPLEX
  void derstring
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real time, const int* /*bcrec*/, int level);
  void deraxvel
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real time, const int* /*bcrec*/, int level);
  void deraxgrad
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real time, const int* /*bcrec*/, int level);
  void deraxgraddens
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real time, const int* /*bcrec*/, int level);
  void deraxphas
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real time, const int* /*bcrec*/, int level);
  void deraxdens
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real time, const int* /*bcrec*/, int level);
  void deraxener
    (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
     const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
     amrex::Real time, const int* /*bcrec*/, int level);
#endif

#ifdef __cplusplus
}
#endif

#endif
