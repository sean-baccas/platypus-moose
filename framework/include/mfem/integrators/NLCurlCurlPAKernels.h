//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "libmesh/ignore_warnings.h"
#include "mfem.hpp"
#include "mfem/fem/integ/bilininteg_hcurl_kernels.hpp"
#include "libmesh/restore_warnings.h"
#include "ScalarLaws.h"

/*
All this header does is declare the functions that will hold
the calls to mfem kernels.

Note that we deliberately do not use the #ifdef MOOSE_MFEM_ENABLED macro
here. We don't include any moose headers (eg MooseConfig.h) so we don't
include the definition for this macro.
*/

namespace Moose::MFEM
{

// Evaluates the reference-element curl of an ND field on hexahedra at the tensor-product
// quadrature points, using sum factorization. This is the curl stage of
// mfem::internal::PACurlCurlApply3D. The physical curl is J c / det(J).
// The params are:
// const int D1D - the number of closed 1D basis functions (ND order + 1)
// const int Q1D - the number of quadpoints in each dimension
// const int NE - the number of elements
// const Array<real_t> & bo - open 1D basis values at the quadpoints (DofToQuadOpen B)
// const Array<real_t> & bc - closed 1D basis values at the quadpoints (DofToQuad B)
// const Array<real_t> & gc - closed 1D basis derivatives at the quadpoints (DofToQuad G)
// const Vector & x - the field as an E-vector in ElementDofOrdering::LEXICOGRAPHIC, which is
//    what NonlinearForm partial assembly passes to AddMultPA and AssembleGradPA
// Vector & curl - the reference curl, laid out as (3, NQ, NE) with the quadpoint index
//    q = qx + Q1D * (qy + Q1D * qz). This matches the layout of a vector QuadratureFunction and
//    the quadpoint ordering of GeometricFactors. Must be sized 3 * Q1D^3 * NE by the caller
void NLCurlCurlReferenceCurl3D(const int D1D,
                               const int Q1D,
                               const int NE,
                               const mfem::Array<mfem::real_t> & bo,
                               const mfem::Array<mfem::real_t> & bc,
                               const mfem::Array<mfem::real_t> & gc,
                               const mfem::Vector & x,
                               mfem::Vector & curl);

// Maps the reference curl from NLCurlCurlReferenceCurl3D to the physical curl in place, using
// (curl u)(F(x)) = J curl_hat(u_hat)(x) / det(J) (Monk, Finite Element Methods for Maxwell's
// Equations, eq. 3.76).
// The params are:
// const int Q1D - the number of quadpoints in each dimension
// const int NE - the number of elements
// const Vector & j - the jacobian matrices from GeometricFactors, laid out as (NQ, 3, 3, NE)
// Vector & curl - on entry the reference curl, on exit the physical curl, laid out as
//    (3, NQ, NE)
void
NLCurlCurlPhysicalCurl3D(const int Q1D, const int NE, const mfem::Vector & j, mfem::Vector & curl);

// The params are:
// const PowerLaw law - this struct has an eval() method to evaluate k().
//    Eventually, we plan to make this more general and support e.g. tabulated function
// const int Q1D - the number of quadpoints in each dimension
// (we skip coeff_dim, as we expect the functions to be scalar functions)
// const int NE - the number of elements
// const Array<real_t> & w, the quadrature weights - get these from the integration rule
// const Vector &j - the jacobian matrices in the form of a flat vector. Its length is
//    9 * number of elements * number of quadpoints, i.e. one 3x3 matrix for each quadpoint
//    in each element
// const Vector& c - this is the curl, evaluated at each qpoint
// Vector & op - this is the diagonal operator
void NLCurlCurlGradPASetup(const Law law,
                           const int Q1D,
                           const int NE,
                           const mfem::Array<mfem::real_t> & w,
                           const mfem::Vector & j,
                           const mfem::Vector & c,
                           mfem::Vector & op,
                           mfem::real_t curlu_zero_tol);

// Same, but creates the op used for AddMultPA
// const PowerLaw law - this struct has an eval() method to evaluate k().
//    Eventually, we plan to make this more general and support e.g. tabulated function
// const int Q1D - the number of quadpoints in each dimension
// (we skip coeff_dim, as we expect the functions to be scalar functions)
// const int NE - the number of elements
// const Array<real_t> & w, the quadrature weights - get these from the integration rule
// const Vector &j - the jacobian matrices in the form of a flat vector. Its length is
//    9 * number of elements * number of quadpoints, i.e. one 3x3 matrix for each quadpoint
//    in each element
// const Vector& c - this is the curl, evaluated at each qpoint
// Vector & op - this is the diagonal operator
void NLCurlCurlPASetup(const Law law,
                       const int Q1D,
                       const int NE,
                       const mfem::Array<mfem::real_t> & w,
                       const mfem::Vector & j,
                       const mfem::Vector & c,
                       mfem::Vector & op);

} // namespace Moose::MFEM
