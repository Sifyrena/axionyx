#ifndef _compute_laplace_operator_H_
#define _compute_laplace_operator_H_

#include <AMReX_FArrayBox.H>

using namespace amrex;

Real compute_laplace_operator (const Array4<const Real>& arr, int i, int j, int k, int comp, Real invdeltasq, int neighbors);

#endif
