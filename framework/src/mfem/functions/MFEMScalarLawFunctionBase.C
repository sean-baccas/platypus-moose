#ifdef MOOSE_MFEM_ENABLED

#include "MFEMScalarLawFunctionBase.h"
#include "MFEMScalarLawCoefficient.h"

InputParameters
MFEMScalarLawFunctionBase::validParams()
{
  InputParameters params = Function::validParams();
  params.addRequiredParam<MFEMScalarCoefficientName>(
      "argument",
      "Scalar coefficient s the law is evaluated at by host consumers, e.g. "
      "'<variable>_curl_mag'.");
  return params;
}

MFEMScalarLawFunctionBase::MFEMScalarLawFunctionBase(const InputParameters & parameters)
  : Function(parameters),
    _mfem_problem(
        cast_ref<MFEMProblem &>(*parameters.getCheckedPointerParam<SubProblem *>("_subproblem"))),
    _argument_name(getParam<MFEMScalarCoefficientName>("argument"))
{
  _mfem_problem.getCoefficients().declareScalar<MFEMScalarLawCoefficient>(
      name(), *this, MFEMScalarLawCoefficient::Quantity::K);
}

void
MFEMScalarLawFunctionBase::initialSetup()
{
  // Variable coefficients such as '<variable>_curl_mag' may be declared after this object
  auto & coef = cast_ref<MFEMScalarLawCoefficient &>(
      _mfem_problem.getCoefficients().getScalarCoefficient(name()));
  coef.SetArgument(_mfem_problem.getCoefficients().getScalarCoefficient(_argument_name));
}

#endif