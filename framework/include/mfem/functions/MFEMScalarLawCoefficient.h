#ifdef MOOSE_MFEM_ENABLED

#pragma once

#include "libmesh/ignore_warnings.h"
#include "mfem.hpp"
#include "libmesh/restore_warnings.h"

class MFEMScalarLawFunctionBase;

/**
 * Host coefficient evaluating a scalar law k(s), or s dk/ds, at the value of an argument
 * coefficient s. Consumers that can evaluate the law in device kernels retrieve it through
 * lawFunction() instead.
 */
class MFEMScalarLawCoefficient : public mfem::Coefficient
{
public:
  /// Which quantity Eval returns
  enum class Quantity
  {
    K,
    S_DK_DS
  };

  MFEMScalarLawCoefficient(const MFEMScalarLawFunctionBase & law_function, Quantity quantity);

  /// Set the coefficient s the law is evaluated at
  void SetArgument(mfem::Coefficient & argument) { _argument = &argument; }
  /// The MOOSE object owning the law
  const MFEMScalarLawFunctionBase & lawFunction() const { return _law_function; }

  /// Set the quadrature function that NLCurlCurlIntegrator::AssembleGradPA fills, on the device,
  /// with k at the current state at the points of its quadrature space. Not owned.
  void setKTarget(mfem::QuadratureFunction * k_target) { _k_target = k_target; }
  /// The quadrature function set by setKTarget, or nullptr if no consumer requested one
  mfem::QuadratureFunction * kTarget() const { return _k_target; }

  void SetTime(mfem::real_t t) override;
  mfem::real_t Eval(mfem::ElementTransformation & T, const mfem::IntegrationPoint & ip) override;

private:
  const MFEMScalarLawFunctionBase & _law_function;
  const Quantity _quantity;
  /// Argument s, set in MFEMScalarLawFunctionBase::initialSetup
  mfem::Coefficient * _argument = nullptr;

  /// Quadrature function filled with k by NLCurlCurlIntegrator::AssembleGradPA. Not owned.
  mfem::QuadratureFunction * _k_target = nullptr;
};


#endif
