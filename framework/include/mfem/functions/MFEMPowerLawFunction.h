#ifdef MOOSE_MFEM_ENABLED

#pragma once

#include "MFEMScalarLawFunctionBase.h"

/**
 * Scalar law k(s) = a + b s^n.
 */
class MFEMPowerLawFunction : public MFEMScalarLawFunctionBase
{
public:
  static InputParameters validParams();

  MFEMPowerLawFunction(const InputParameters & parameters);

  Law getLaw(bool /*on_device*/) const override { return Moose::MFEM::PowerLaw{_a, _b, _n}; }

private:
  const Real _a;
  const Real _b;
  const Real _n;
};

#endif
