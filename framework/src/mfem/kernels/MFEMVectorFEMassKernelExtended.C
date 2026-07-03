//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMVectorFEMassKernelExtended.h"
#include "MFEMProblem.h"

registerMooseObject("MooseApp", MFEMVectorFEMassKernelExtended);

// lazily define new class, inheriting from the old one
class VectorFEMassIntegratorExtended : public mfem::VectorFEMassIntegrator
{
public:
  VectorFEMassIntegratorExtended(mfem::Coefficient & coef) : mfem::VectorFEMassIntegrator(coef) {}

  // next, implement the AssembleGradPA method.
  // remember, the BilinearFormIntegrators DO inherit from NonlinearFormIntegrator
  void AssembleGradPA(const mfem::Vector &, const mfem::FiniteElementSpace &fes) override {
    mfem::VectorFEMassIntegrator::AssemblePA(fes);
  }

  void AddMultGradPA(const mfem::Vector &x, mfem::Vector &y) const override {
    mfem::VectorFEMassIntegrator::AddMultPA(x, y);
  }

  void AssembleGradDiagonalPA(mfem::Vector &diag) const override {
    auto * self = const_cast<VectorFEMassIntegratorExtended *>(this);
    self->mfem::VectorFEMassIntegrator::AssembleDiagonalPA(diag);
  }
};


InputParameters
MFEMVectorFEMassKernelExtended::validParams()
{
  InputParameters params = MFEMKernel::validParams();
  params.addClassDescription("Adds the domain integrator to an MFEM problem for the bilinear form "
                             "$(k \\vec u, \\vec v)_\\Omega$ "
                             "arising from the weak form of the mass operator "
                             "$k \\vec u$.");
  params.addParam<MFEMScalarCoefficientName>(
      "coefficient", "1.", "Name of property k to multiply the integrator by");
  return params;
}

MFEMVectorFEMassKernelExtended::MFEMVectorFEMassKernelExtended(const InputParameters & parameters)
  : MFEMKernel(parameters), _coef(getScalarCoefficient("coefficient"))
// FIXME: The MFEM bilinear form can also handle vector and matrix
// coefficients, so ideally we'd handle all three too.
{
}

mfem::BilinearFormIntegrator *
MFEMVectorFEMassKernelExtended::createBFIntegrator()
{
  return new VectorFEMassIntegratorExtended(_coef);
}

#endif
