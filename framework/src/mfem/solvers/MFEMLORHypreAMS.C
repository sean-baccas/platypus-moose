//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMLORHypreAMS.h"
#include "MFEMProblem.h"
#include "MFEMScalarLawCoefficient.h"

#include "libmesh/ignore_warnings.h"
#include "mfem/fem/lor/lor_batched.hpp"
#include "libmesh/restore_warnings.h"

registerMooseObject("MooseApp", MFEMLORHypreAMS);

InputParameters
MFEMLORHypreAMS::validParams()
{
  InputParameters params = Moose::MFEM::LORLinearSolverBase<mfem::HypreAMS>::validParams();
  params.addClassDescription(
      "Hypre auxiliary-space Maxwell preconditioner on the low-order-refined discretization of "
      "curl(alpha curl u) + beta u, rebuilt for every operator it preconditions.");
  params.addParam<MFEMScalarCoefficientName>(
      "alpha_coefficient",
      "1.",
      "Name of scalar coefficient used in the curl-curl component of the surrogate form. If it is "
      "the k_coefficient of an MFEMNLCurlCurlKernel and is a scalar law such as "
      "MFEMPowerLawFunction, it is evaluated on the device at the current Newton iterate.");
  params.addParam<MFEMScalarCoefficientName>(
      "beta_coefficient",
      "1.",
      "Name of scalar coefficient used in the mass component of the surrogate form.");
  params.addParam<int>("print_level", 0, "Set the solver verbosity.");
  // Always an LOR solver
  params.setParameters("low_order_refined", true);
  params.suppressParameter<bool>("low_order_refined");
  return params;
}

MFEMLORHypreAMS::MFEMLORHypreAMS(const InputParameters & parameters)
  : Moose::MFEM::LORLinearSolverBase<mfem::HypreAMS>(parameters),
    _alpha_coef(getScalarCoefficient("alpha_coefficient")),
    _beta_coef(getScalarCoefficient("beta_coefficient")),
    _alpha_law_coef(dynamic_cast<MFEMScalarLawCoefficient *>(&_alpha_coef))
{
  ConstructSolver();
}

void
MFEMLORHypreAMS::ConstructSolver()
{
  _solver = std::make_unique<Rebuilder>(*this);
}

void
MFEMLORHypreAMS::UpdateEquationSystemContext()
{
  if (isParamSetByUser("preconditioner"))
    paramError("preconditioner", "LOR solver cannot take a preconditioner");

  SetupLOR();
  auto & fespace = *_a->ParFESpace();

  // Release in reverse order of construction, since each member references the one before it
  _lor_ams.reset();
  _surrogate.reset();
  _k_coef.reset();
  _k_qf.reset();
  _k_qspace.reset();

  mfem::Coefficient * alpha = &_alpha_coef;
  if (_alpha_law_coef)
  {
    // Batched LOR assembly takes a QuadratureFunctionCoefficient by reference rather than
    // evaluating it, after checking only the mesh and the order of the rule; the rule must also
    // have the same points, so it is the one batched LOR uses
    _k_ir = mfem::GetCollocatedIntRule(fespace);
    _k_qspace = std::make_unique<mfem::QuadratureSpace>(*fespace.GetMesh(), _k_ir);
    _k_qf = std::make_unique<mfem::QuadratureFunction>(*_k_qspace);
    // Placeholder for solves that never assemble a gradient
    *_k_qf = 1.0;
    _k_coef = std::make_unique<mfem::QuadratureFunctionCoefficient>(*_k_qf);
    _alpha_law_coef->setKTarget(_k_qf.get());
    alpha = _k_coef.get();
  }

  _surrogate = std::make_unique<mfem::ParBilinearForm>(&fespace);
  _surrogate->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*alpha));
  _surrogate->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(_beta_coef));

  // Otherwise mfem silently falls back to host assembly on a refined mesh
  if (!mfem::BatchedLORAssembly::FormIsSupported(*_surrogate))
    mooseError("MFEMLORHypreAMS requires batched LOR assembly, which needs a tensor-product "
               "(hexahedral) mesh.");

  _stale = true;
}

void
MFEMLORHypreAMS::SetSolverParameters(mfem::HypreAMS & solver)
{
  solver.iterative_mode = getParam<bool>("use_initial_guess");
  solver.SetPrintLevel(getParam<int>("print_level"));
}

void
MFEMLORHypreAMS::BuildSolver()
{
  // Batched LOR assembly, gradient matrix and vertex coordinates. hypre performs the AMS setup
  // itself on the first apply, so that cost is timed in "apply".
  TIME_SECTION("setup", 3);
  // Free the previous hypre setup before allocating the next
  _lor_ams.reset();
  _lor_ams = std::make_unique<mfem::LORSolver<mfem::HypreAMS>>(*_surrogate, _ess_tdofs);
  SetSolverParameters(_lor_ams->GetSolver());
  _stale = false;
  // Kernel launches are asynchronous; wait for them so their time lands in this section
  MFEM_DEVICE_SYNC;
}

void
MFEMLORHypreAMS::ApplySolver(const mfem::Vector & x, mfem::Vector & y)
{
  TIME_SECTION("apply", 3);
  _lor_ams->Mult(x, y);
  MFEM_DEVICE_SYNC;
}

void
MFEMLORHypreAMS::Rebuilder::SetOperator(const mfem::Operator & op)
{
  height = op.Height();
  width = op.Width();
  // Build on the next Mult rather than here: the operator set before the Newton solve is never
  // applied, and every gradient assembly runs before its SetOperator, so alpha is current then
  _owner._stale = true;
}

void
MFEMLORHypreAMS::Rebuilder::Mult(const mfem::Vector & x, mfem::Vector & y) const
{
  if (_owner._stale)
    _owner.BuildSolver();
  _owner.ApplySolver(x, y);
}

#endif

