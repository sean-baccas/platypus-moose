//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#pragma once
#include "MFEMSolverBase.h"
#include "libmesh/ignore_warnings.h"
#include <mfem.hpp>
#include "libmesh/restore_warnings.h"
#include <memory>

/**
 * Wrapper for mfem::HypreBoomerAMG solver.
 */
 class MFEMMinResSolver : public MFEMSolverBase
 {
public:
  static InputParameters validParams();

  MFEMMinResSolver(const InputParameters &);

  /// Returns a shared pointer to the instance of the Solver derived-class.
  void updateSolver(mfem::ParBilinearForm & a, mfem::Array<int> & tdofs) override;

protected:
  void constructSolver(const InputParameters & parameters) override;
};

#endif
