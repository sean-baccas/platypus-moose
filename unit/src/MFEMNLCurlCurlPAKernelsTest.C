//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "gtest/gtest.h"
#include "NLCurlCurlPAKernels.h"
#include "libmesh/int_range.h"

namespace
{
/// Smooth field whose curl has a nonzero, spatially varying component in every direction
void
field(const mfem::Vector & p, mfem::Vector & u)
{
  u(0) = std::sin(p(1)) * p(2) * p(2);
  u(1) = std::cos(p(2)) * p(0);
  u(2) = std::exp(p(0)) * p(1);
}

/// Warps the unit cube so that the element Jacobians are non-affine and vary within each element
void
warp(const mfem::Vector & x, mfem::Vector & y)
{
  y = x;
  y(0) += 0.1 * std::sin(M_PI * x(1)) * x(2);
  y(1) += 0.1 * std::sin(M_PI * x(2)) * x(0);
  y(2) += 0.1 * std::sin(M_PI * x(0)) * x(1);
}
}

/**
 * The sum-factorized curl kernels must reproduce the element-by-element host evaluation of
 * mfem::CurlGridFunctionCoefficient at every quadrature point.
 */
TEST(MFEMNLCurlCurlPAKernelsTest, PhysicalCurlMatchesCurlGridFunctionCoefficient)
{
  for (const int order : {1, 2, 3})
  {
    auto mesh = mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::HEXAHEDRON);
    // Second-order geometry, so that warp produces curved elements
    mesh.SetCurvature(2);
    mesh.Transform(warp);

    mfem::ND_FECollection fec(order, 3);
    mfem::FiniteElementSpace fes(&mesh, &fec);
    mfem::GridFunction u(&fes);
    mfem::VectorFunctionCoefficient u_coef(3, field);
    u.ProjectCoefficient(u_coef);

    // The rule NLCurlCurlIntegrator::PreAssemblySetup uses when none is given
    const auto & el = dynamic_cast<const mfem::VectorTensorFiniteElement &>(*fes.GetTypicalFE());
    const auto & ir = mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder());
    const auto & maps_c = el.GetDofToQuad(ir, mfem::DofToQuad::TENSOR);
    const auto & maps_o = el.GetDofToQuadOpen(ir, mfem::DofToQuad::TENSOR);
    const auto * geom = mesh.GetGeometricFactors(ir, mfem::GeometricFactors::JACOBIANS);
    const int ne = mesh.GetNE();

    // The E-vector that NonlinearForm partial assembly passes to the integrator
    const auto * restriction = fes.GetElementRestriction(mfem::ElementDofOrdering::LEXICOGRAPHIC);
    mfem::Vector ue(restriction->Height());
    restriction->Mult(u, ue);

    mfem::Vector curl(3 * ir.GetNPoints() * ne);
    Moose::MFEM::NLCurlCurlReferenceCurl3D(
        maps_c.ndof, maps_c.nqpt, ne, maps_o.B, maps_c.B, maps_c.G, ue, curl);
    Moose::MFEM::NLCurlCurlPhysicalCurl3D(maps_c.nqpt, ne, geom->J, curl);

    mfem::QuadratureSpace qspace(mesh, ir);
    mfem::QuadratureFunction expected(qspace, 3);
    mfem::CurlGridFunctionCoefficient curl_coef(&u);
    curl_coef.Project(expected);

    // Both sides evaluate the same discrete field and differ only in summation order, so they
    // should agree to a modest multiple of machine epsilon relative to the curl's magnitude
    const mfem::real_t tol = 1e-12 * expected.Normlinf();
    curl.HostRead();
    for (const auto i : libMesh::make_range(curl.Size()))
      EXPECT_NEAR(curl(i), expected(i), tol) << "order " << order << ", entry " << i;
  }
}

#endif
