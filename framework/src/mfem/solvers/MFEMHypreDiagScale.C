//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMHypreDiagScale.h"
#include "MFEMProblem.h"

registerMooseObject("MooseApp", MFEMHypreDiagScale);

InputParameters
MFEMHypreDiagScale::validParams()
{
  InputParameters params = MFEMSolverBase::validParams();
  params.addClassDescription("Hypre DiagScale solver and preconditioner for the iterative solution "
                             "of MFEM equation systems.");
  params.addParam<double>("l_tol", 1e-5, "Set the relative tolerance.");
  params.addParam<int>("l_max_its", 10000, "Set the maximum number of iterations.");
  params.addParam<int>("print_level", 2, "Set the solver verbosity.");
  return params;
}

MFEMHypreDiagScale::MFEMHypreDiagScale(const InputParameters & parameters)
  : MFEMSolverBase(parameters)
{
  constructSolver(parameters);
}

void
MFEMHypreDiagScale::constructSolver(const InputParameters & /*parameters*/)
{
  auto solver = std::make_unique<mfem::HypreDiagScale>();
  _solver = std::move(solver);
}

void
MFEMHypreDiagScale::updateSolver(mfem::ParBilinearForm & /*a*/, mfem::Array<int> & /*tdofs*/)
{
  if (_lor)
    mooseError("LOR solver not configured yet");
}
#endif
