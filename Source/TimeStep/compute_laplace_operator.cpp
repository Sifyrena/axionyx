#include "Nyx.H"
#include <AMReX_MultiFab.H>

using namespace amrex;

Real compute_laplace_operator (const Array4<const Real>& arr, int i, int j, int k, int comp, Real invdeltasq, int neighbors)
{
  switch(neighbors){
  case 1:
  return (  6.0*arr(i+1,j  ,k  ,comp)
	   +6.0*arr(i-1,j  ,k  ,comp)
	   +6.0*arr(i  ,j+1,k  ,comp)
	   +6.0*arr(i  ,j-1,k  ,comp)
	   +6.0*arr(i  ,j  ,k+1,comp)
	   +6.0*arr(i  ,j  ,k-1,comp)
	   +3.0*arr(i+1,j+1,k  ,comp)
	   +3.0*arr(i+1,j-1,k  ,comp)
	   +3.0*arr(i-1,j+1,k  ,comp)
	   +3.0*arr(i-1,j-1,k  ,comp)
	   +3.0*arr(i+1,j  ,k+1,comp)
	   +3.0*arr(i+1,j  ,k-1,comp)
	   +3.0*arr(i-1,j  ,k+1,comp)
	   +3.0*arr(i-1,j  ,k-1,comp)
	   +3.0*arr(i  ,j+1,k+1,comp)
	   +3.0*arr(i  ,j+1,k-1,comp)
	   +3.0*arr(i  ,j-1,k+1,comp)
	   +3.0*arr(i  ,j-1,k-1,comp)
	   +2.0*arr(i+1,j+1,k+1,comp)
	   +2.0*arr(i+1,j+1,k-1,comp)
	   +2.0*arr(i+1,j-1,k+1,comp)
	   +2.0*arr(i+1,j-1,k-1,comp)
	   +2.0*arr(i-1,j+1,k+1,comp)
	   +2.0*arr(i-1,j+1,k-1,comp)
	   +2.0*arr(i-1,j-1,k+1,comp)
	   +2.0*arr(i-1,j-1,k-1,comp)
	  -88.0*arr(i  ,j  ,k  ,comp))*invdeltasq/26.0;
  case 2:
  return (-     arr(i-2,j  ,k  ,comp)
	  +16.0*arr(i-1,j  ,k  ,comp)
	  +16.0*arr(i+1,j  ,k  ,comp)
	  -     arr(i+2,j  ,k  ,comp)
	  -     arr(i  ,j-2,k  ,comp)
	  +16.0*arr(i  ,j-1,k  ,comp)
	  +16.0*arr(i  ,j+1,k  ,comp)
	  -     arr(i  ,j+2,k  ,comp)
	  -     arr(i  ,j  ,k-2,comp)
	  +16.0*arr(i  ,j  ,k-1,comp)
	  +16.0*arr(i  ,j  ,k+1,comp)
	  -     arr(i  ,j  ,k+2,comp)
	  -90.0*arr(i  ,j  ,k  ,comp))*invdeltasq/12.0;
  default:
	  amrex::Error("Need one or two neighbors for Laplace operator");
	  return 0.0;
  }	  
}
