#ifdef MOOSE_MFEM_ENABLED

#pragma once

#include "Function.h"
#include "MFEMProblem.h"
#include "ScalarLaws.h"
// #include <variant>

/**
 * Base class for MOOSE objects declaring a scalar law k(s) that MFEM integrators can evaluate in
 * device kernels. Declares a host coefficient of the same name evaluating k at the coefficient
 * named by 'argument', for consumers that need an mfem::Coefficient.
 */
class MFEMScalarLawFunctionBase : public Function
{
public:
  // This should be a variant, but for now we just point it at the alias
  // statement in ScalarLaws.h
  using Law = Moose::MFEM::Law;

  static InputParameters validParams();

  MFEMScalarLawFunctionBase(const InputParameters & parameters);

  void initialSetup() override;

  /// The law, with any pointers it holds valid on the device if @a on_device, otherwise on the
  /// host, until the next call that changes the law such as extendRange
  virtual Law getLaw(bool on_device) const = 0;
  /// Make the law accurate for arguments up to @a s; a no-op for laws valid for all s
  virtual void extendRange(mfem::real_t /*s*/) {}
  /// Name of the argument coefficient, so consumers can check it is the one they evaluate
  const MFEMScalarCoefficientName & argumentName() const { return _argument_name; }

protected:
  /// reference to the MFEMProblem instance
  MFEMProblem & _mfem_problem;
  /// Name of the coefficient s the law is evaluated at
  const MFEMScalarCoefficientName & _argument_name;
};

#endif