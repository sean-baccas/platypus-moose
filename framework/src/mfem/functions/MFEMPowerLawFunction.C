#ifdef MOOSE_MFEM_ENABLED

#include "MFEMPowerLawFunction.h"

registerMooseObject("MooseApp", MFEMPowerLawFunction);

InputParameters
MFEMPowerLawFunction::validParams()
{
  InputParameters params = MFEMScalarLawFunctionBase::validParams();
  params.addClassDescription("Scalar law $k(s) = a + b s^n$, evaluable in MFEM device kernels.");
  // a > 0, b >= 0 and n >= 0 make k positive and s k(s) strictly increasing for all s >= 0,
  // which keeps the nonlinear curl-curl problem well posed
  params.addRequiredRangeCheckedParam<Real>("a", "a > 0", "Value of $k$ at $s = 0$.");
  params.addRequiredRangeCheckedParam<Real>("b", "b >= 0", "Coefficient of the power term.");
  params.addRequiredRangeCheckedParam<Real>("n", "n >= 0", "Exponent of the power term.");
  return params;
}

MFEMPowerLawFunction::MFEMPowerLawFunction(const InputParameters & parameters)
  : MFEMScalarLawFunctionBase(parameters),
    _a(getParam<Real>("a")),
    _b(getParam<Real>("b")),
    _n(getParam<Real>("n"))
{
}

#endif
