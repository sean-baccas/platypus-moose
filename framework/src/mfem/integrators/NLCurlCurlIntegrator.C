//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*

#ifdef MOOSE_MFEM_ENABLED

#include "NLCurlCurlIntegrator.h"
#include "NLCurlCurlPAKernels.h"
#include "MFEMScalarLawFunctionBase.h"
#include "MFEMScalarLawCoefficient.h"

namespace Moose::MFEM
{
NLCurlCurlJacMatrixCoefficient::NLCurlCurlJacMatrixCoefficient(mfem::Coefficient & k,
                                                               mfem::Coefficient & curlu_dk_dcurlu,
                                                               mfem::VectorCoefficient & curlu_vec,
                                                               mfem::real_t curlu_zero_tol)
  : mfem::SymmetricMatrixCoefficient(curlu_vec.GetVDim()),
    _k_coef(k),
    _curlu_dk_dcurlu_coef(curlu_dk_dcurlu),
    _curlu_zero_tol(curlu_zero_tol),
    _curlu_hat_coef(curlu_vec, _curlu_zero_tol)
{
}

void
NLCurlCurlJacMatrixCoefficient::SetTime(mfem::real_t t)
{
  SymmetricMatrixCoefficient::SetTime(t);
  _k_coef.SetTime(t);
  _curlu_dk_dcurlu_coef.SetTime(t);
  _curlu_hat_coef.SetTime(t);
}

void
NLCurlCurlJacMatrixCoefficient::Eval(mfem::DenseSymmetricMatrix & K,
                                     mfem::ElementTransformation & T,
                                     const mfem::IntegrationPoint & ip)
{
  const int dim = GetSize();
  mfem::Vector curlu_hat(dim);

  _curlu_hat_coef.Eval(curlu_hat, T, ip);
  const mfem::real_t k = _k_coef.Eval(T, ip);
  const mfem::real_t curlu_dk_dcurlu = _curlu_dk_dcurlu_coef.Eval(T, ip);

  // K is sized by the caller (ProjectSymmetric points it at the quadrature data)
  for (const auto i : make_range(dim))
    for (const auto j : make_range(i, dim))
      K(i, j) = (i == j ? k : 0.0) + curlu_dk_dcurlu * curlu_hat(i) * curlu_hat(j);
}

NLCurlCurlIntegrator::NLCurlCurlIntegrator(mfem::Coefficient & k,
                                           mfem::Coefficient & curlu_dk_dcurlu,
                                           mfem::VectorCoefficient & curlu,
                                           mfem::real_t curlu_zero_tol,
                                           const mfem::IntegrationRule * ir)
  : mfem::NonlinearFormIntegrator(ir), // so the PA path sees the rule
    _curlcurl_res_integ(k, ir),
    _curlcurl_jac_matrix_coef(k, curlu_dk_dcurlu, curlu, curlu_zero_tol),
    _curlcurl_jac_integ(_curlcurl_jac_matrix_coef, ir),
    _k_coef(k),
    _curlu_dk_dcurlu_coef(curlu_dk_dcurlu),
    _curlu_coef(curlu),
    _curlu_zero_tol(curlu_zero_tol),
    _law_coef(dynamic_cast<const MFEMScalarLawCoefficient*>(&k))
{
}

void
NLCurlCurlIntegrator::AssembleElementVector(const mfem::FiniteElement & el,
                                            mfem::ElementTransformation & Tr,
                                            const mfem::Vector & elfun,
                                            mfem::Vector & elvect)
{
  _curlcurl_res_integ.AssembleElementVector(el, Tr, elfun, elvect);
}

void
NLCurlCurlIntegrator::AssembleElementGrad(const mfem::FiniteElement & el,
                                          mfem::ElementTransformation & Tr,
                                          const mfem::Vector & elfun,
                                          mfem::DenseMatrix & elmat)
{
  _curlcurl_jac_integ.AssembleElementGrad(el, Tr, elfun, elmat);
}

void
NLCurlCurlIntegrator::AssembleGradPA(const mfem::Vector & x,
                                     const mfem::FiniteElementSpace & fes)
{
  PreAssemblySetup(fes);
  _pa_grad_data.SetSize(_ndata * _nq * _ne, mfem::Device::GetMemoryType());

  if (_law_coef)
  {
    ComputeCurlAtQuadPoints(x);
    Moose::MFEM::NLCurlCurlGradPASetup(_law_coef->lawFunction().getLaw(/*on_device=*/true),
                                       _quad1D,
                                       _ne,
                                       _qspace_ir->GetWeights(),
                                       _geom->J,
                                       _curl_qp,
                                       _pa_grad_data,
                                       _curlu_zero_tol);
    if (_law_coef->kTarget())
      ComputeKAtTarget(x, fes);
    return;
  }

  mfem::CoefficientVector jac_cv(*_qspace, mfem::CoefficientStorage::SYMMETRIC);
  jac_cv.Project(_curlcurl_jac_matrix_coef);

  mfem::internal::PACurlCurlSetup3D(
      _quad1D, jac_cv.GetVDim(), _ne, _qspace_ir->GetWeights(), _geom->J, jac_cv, _pa_grad_data);
}

void
NLCurlCurlIntegrator::AssemblePA(const mfem::FiniteElementSpace & fes)
{
  // We delay projecting k onto qpoints until AddMultPA. By then, we
  // have access to the current state x as an E-vector
  PreAssemblySetup(fes);

  // If there is a ScalarCoefficientLaw, then we can evaluate k(s) and sk'(s)
  // on the device
  if (_law_coef)
    return;

  // original path: projec k onto the coefficient vector and then use the standard setup
  // from the mfem::internal namespace
  mfem::CoefficientVector k_cv(*_qspace, mfem::CoefficientStorage::FULL);
  k_cv.Project(_k_coef);
  _pa_res_data.SetSize(_ndata * _nq * _ne, mfem::Device::GetMemoryType());
  mfem::internal::PACurlCurlSetup3D(
      _quad1D, _k_coef.GetVDim(), _ne, _qspace_ir->GetWeights(), _geom->J, k_cv, _pa_res_data);
}

void
NLCurlCurlIntegrator::AddMultGradPA(const mfem::Vector & x, mfem::Vector & y) const
{
  mfem::CurlCurlIntegrator::ApplyPAKernels::Run(_dim,
                                                _dofs1D,
                                                _quad1D,
                                                _dofs1D,
                                                _quad1D,
                                                _symmetric,
                                                _ne,
                                                _mapsO->B,
                                                _mapsC->B,
                                                _mapsO->Bt,
                                                _mapsC->Bt,
                                                _mapsC->G,
                                                _mapsC->Gt,
                                                _pa_grad_data,
                                                x,
                                                y,
                                                false);
}

void
NLCurlCurlIntegrator::AddMultPA(const mfem::Vector & x, mfem::Vector & y) const
{
  // x is the state as an E-vector, since NonlinearForm::Mult is the only caller
  if (_law_coef)
  {
    ComputeCurlAtQuadPoints(x);
    _pa_res_data.SetSize(_ndata * _nq * _ne, mfem::Device::GetMemoryType());
    Moose::MFEM::NLCurlCurlPASetup(_law_coef->lawFunction().getLaw(/*on_device=*/true),
                                   _quad1D,
                                   _ne,
                                   _qspace_ir->GetWeights(),
                                   _geom->J,
                                   _curl_qp,
                                   _pa_res_data);
  }

  mfem::CurlCurlIntegrator::ApplyPAKernels::Run(_dim,
                                                _dofs1D,
                                                _quad1D,
                                                _dofs1D,
                                                _quad1D,
                                                _symmetric,
                                                _ne,
                                                _mapsO->B,
                                                _mapsC->B,
                                                _mapsO->Bt,
                                                _mapsC->Bt,
                                                _mapsC->G,
                                                _mapsC->Gt,
                                                _pa_res_data,
                                                x,
                                                y,
                                                false);
}

void
NLCurlCurlIntegrator::AssembleGradDiagonalPA(mfem::Vector & diag) const
{
  mfem::CurlCurlIntegrator::DiagonalPAKernels::Run(_dim,
                                                   _dofs1D,
                                                   _quad1D,
                                                   _dofs1D,
                                                   _quad1D,
                                                   _symmetric,
                                                   _ne,
                                                   _mapsO->B,
                                                   _mapsC->B,
                                                   _mapsO->G,
                                                   _mapsC->G,
                                                   _pa_grad_data,
                                                   diag);
}

// This gets called by both AssemblePA and AssembleGradPA, since they
// both need to store the transformation jacobians.
void
NLCurlCurlIntegrator::PreAssemblySetup(const mfem::FiniteElementSpace & fes)
{
  const mfem::FiniteElement * fel = fes.GetTypicalFE();
  mfem::Mesh * mesh = fes.GetMesh();

  // crucial check to see if it casts into VTFE
  const mfem::VectorTensorFiniteElement * el =
      dynamic_cast<const mfem::VectorTensorFiniteElement *>(fel);

  mooseAssert(el, "Only VectorTensorFiniteElement is supported!");
  mooseAssert(el->GetDerivType() == mfem::FiniteElement::CURL, "Unknown kernel type");

  // Pick integration rule: either the one that was passed into the constructor,
  // or the default found in mfem::CurlCurlIntegrator::AssembleElementMatrix. This
  // is so that we can make sure it matches the legacy version.
  const mfem::IntegrationRule * rule =
      IntRule ? IntRule : &mfem::IntRules.Get(el->GetGeomType(), 2 * el->GetOrder());
  mooseAssert(el->GetDim() == 3, "Following methods are only implemented in 3D");

  _nq = rule->GetNPoints();
  _dim = mesh->Dimension();
  mooseAssert(_dim == 3, "Following methods are only implemented in 3D");

  _ne = fes.GetNE();
  _geom = mesh->GetGeometricFactors(*rule, mfem::GeometricFactors::JACOBIANS);
  _mapsC = &el->GetDofToQuad(*rule, mfem::DofToQuad::TENSOR);
  _mapsO = &el->GetDofToQuadOpen(*rule, mfem::DofToQuad::TENSOR);
  _dofs1D = _mapsC->ndof;
  _quad1D = _mapsC->nqpt;

  mooseAssert(_dofs1D == _mapsO->ndof + 1 && _quad1D == _mapsO->nqpt,
              "Must be one fewer open basis function per _dim than closed per dim.");

  if (!_qspace || _qspace_mesh != mesh || _qspace_mesh_sequence != mesh->GetSequence() ||
      _qspace_ir != rule)
  {
    _qspace = std::make_unique<mfem::QuadratureSpace>(*mesh, *rule);
    _qspace_mesh = mesh;
    _qspace_mesh_sequence = mesh->GetSequence();
    _qspace_ir = rule;
  }
}

void NLCurlCurlIntegrator::ComputeCurlAtQuadPoints(const mfem::Vector & x) const
{
  _curl_qp.SetSize(3 * _nq * _ne, mfem::Device::GetMemoryType());
  Moose::MFEM::NLCurlCurlReferenceCurl3D(
      _dofs1D, _quad1D, _ne, _mapsO->B, _mapsC->B, _mapsC->G, x, _curl_qp);
  Moose::MFEM::NLCurlCurlPhysicalCurl3D(_quad1D, _ne, _geom->J, _curl_qp);
}

void
NLCurlCurlIntegrator::ComputeKAtTarget(const mfem::Vector & x,
                                       const mfem::FiniteElementSpace & fes)
{
  mfem::QuadratureFunction & k = *_law_coef->kTarget();
  mooseAssert(k.GetSpace()->GetMesh() == fes.GetMesh(),
              "The k target must be defined on the mesh of the integrator's FE space.");

  // The maps and geometric factors are cached by mfem against the address of the rule, which the
  // target's quadrature space keeps alive
  const mfem::IntegrationRule & ir = k.GetSpace()->GetIntRule(0);
  const auto * el = static_cast<const mfem::VectorTensorFiniteElement *>(fes.GetTypicalFE());
  const mfem::DofToQuad & maps_c = el->GetDofToQuad(ir, mfem::DofToQuad::TENSOR);
  const mfem::DofToQuad & maps_o = el->GetDofToQuadOpen(ir, mfem::DofToQuad::TENSOR);
  const mfem::GeometricFactors * geom =
      fes.GetMesh()->GetGeometricFactors(ir, mfem::GeometricFactors::JACOBIANS);
  const int quad1D = maps_c.nqpt;

  _k_target_curl.SetSize(3 * k.Size(), mfem::Device::GetMemoryType());
  Moose::MFEM::NLCurlCurlReferenceCurl3D(
      _dofs1D, quad1D, _ne, maps_o.B, maps_c.B, maps_c.G, x, _k_target_curl);
  Moose::MFEM::NLCurlCurlPhysicalCurl3D(quad1D, _ne, geom->J, _k_target_curl);
  Moose::MFEM::NLCurlCurlLawK(_law_coef->lawFunction().getLaw(/*on_device=*/true),
                              k.Size(),
                              _k_target_curl,
                              k);
}
}

#endif
