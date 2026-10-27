#ifdef MOOSE_MFEM_ENABLED

#include "MFEMScalarLawCoefficient.h"
#include "MFEMScalarLawFunctionBase.h"

MFEMScalarLawCoefficient::MFEMScalarLawCoefficient(const MFEMScalarLawFunctionBase & law_function,
                                                   const Quantity quantity)
  : _law_function(law_function), _quantity(quantity)
{
}

void
MFEMScalarLawCoefficient::SetTime(const mfem::real_t t)
{
  mfem::Coefficient::SetTime(t);
  if (_argument)
    _argument->SetTime(t);
}

mfem::real_t
MFEMScalarLawCoefficient::Eval(mfem::ElementTransformation & T, const mfem::IntegrationPoint & ip)
{
  mooseAssert(_argument, "The argument of '" + _law_function.name() + "' has not been set.");
  const mfem::real_t s = _argument->Eval(T, ip);
  mfem::real_t k, s_dk_ds;
  // std::visit([&](const auto & law) { law.eval(s, k, s_dk_ds); },
  //            _law_function.getLaw(/*on_device=*/false));
  _law_function.getLaw(/*on_device=*/false).eval(s, k, s_dk_ds);
  return _quantity == Quantity::K ? k : s_dk_ds;
}

#endif
