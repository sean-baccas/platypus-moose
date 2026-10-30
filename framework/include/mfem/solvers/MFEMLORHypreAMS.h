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

#include "MFEMLORLinearSolverBase.h"

class MFEMScalarLawCoefficient;

/**
 * Hypre auxiliary-space Maxwell preconditioner on the low-order-refined discretization of the
 * surrogate form CurlCurl(alpha) + VectorFEMass(beta), rebuilt for every operator it is given.
 * The surrogate form is independent of the equation system's bilinear form, so this can
 * precondition a Jacobian whose curl-curl term belongs to the nonlinear form. When alpha is an
 * MFEMScalarLawCoefficient, NLCurlCurlIntegrator::AssembleGradPA evaluates it on the device and
 * the LOR system is assembled on the device.
 */
class MFEMLORHypreAMS : public Moose::MFEM::LORLinearSolverBase<mfem::HypreAMS>
{
public:
  static InputParameters validParams();

  MFEMLORHypreAMS(const InputParameters &);

  void ConstructSolver() override;

protected:
  void UpdateEquationSystemContext() override;
  void SetSolverParameters(mfem::HypreAMS & solver) override;

private:
  // Builds the LOR HypreAMS solver from the surrogate form
  void BuildSolver();
  // Applies the solver built by BuildSolver()
  void ApplySolver(const mfem::Vector & x, mfem::Vector & y);

  /**
   * Marks the solver stale on every SetOperator call and rebuilds it on the next Mult. Within a
   * Newton solve the Krylov solver passes each new Jacobian to its preconditioner only through
   * mfem::Solver::SetOperator. The operator is used only for its size.
   */
  class Rebuilder : public mfem::Solver
  {
  public:
    explicit Rebuilder(MFEMLORHypreAMS & owner) : _owner(owner) {}
    void SetOperator(const mfem::Operator & op) override;
    void Mult(const mfem::Vector & x, mfem::Vector & y) const override;

  private:
    MFEMLORHypreAMS & _owner;
  };

  mfem::Coefficient & _alpha_coef;
  mfem::Coefficient & _beta_coef;
  /// _alpha_coef if it is a law the device can evaluate, otherwise nullptr
  MFEMScalarLawCoefficient * const _alpha_law_coef;

  /// Gauss-Lobatto rule collocated with the LOR vertices, at which batched LOR assembly reads
  /// alpha. Referenced by _k_qspace.
  mfem::IntegrationRule _k_ir;
  std::unique_ptr<mfem::QuadratureSpace> _k_qspace;
  /// k at the points of _k_ir, filled on the device by NLCurlCurlIntegrator::AssembleGradPA
  std::unique_ptr<mfem::QuadratureFunction> _k_qf;
  std::unique_ptr<mfem::QuadratureFunctionCoefficient> _k_coef;
  /// Surrogate form CurlCurl(alpha) + VectorFEMass(beta); never assembled at high order
  std::unique_ptr<mfem::ParBilinearForm> _surrogate;
  /// Whether the operator changed since the solver was last built
  bool _stale = true;
  std::unique_ptr<mfem::LORSolver<mfem::HypreAMS>> _lor_ams;
};

#endif

