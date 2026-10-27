//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*

#include "NLCurlCurlPAKernels.h"

namespace Moose::MFEM
{

// TODO: write a unit test that checks this at multiple orders
// on a non-affine mesh
//
// Note: x is supposed to be lemenet vector
void
NLCurlCurlReferenceCurl3D(const int D1D,
                          const int Q1D,
                          const int NE,
                          const mfem::Array<mfem::real_t> & bo,
                          const mfem::Array<mfem::real_t> & bc,
                          const mfem::Array<mfem::real_t> & gc,
                          const mfem::Vector & x,
                          mfem::Vector & curl)
{
  MFEM_VERIFY(D1D <= mfem::DeviceDofQuadLimits::Get().HCURL_MAX_D1D, "D1D > HCURL_MAX_D1D");
  MFEM_VERIFY(Q1D <= mfem::DeviceDofQuadLimits::Get().HCURL_MAX_Q1D, "Q1D > HCURL_MAX_Q1D");

  auto Bo = mfem::Reshape(bo.Read(), Q1D, D1D - 1);
  auto Bc = mfem::Reshape(bc.Read(), Q1D, D1D);
  auto Gc = mfem::Reshape(gc.Read(), Q1D, D1D);
  auto X = mfem::Reshape(x.Read(), 3 * (D1D - 1) * D1D * D1D, NE);
  auto C = mfem::Reshape(curl.Write(), 3, Q1D, Q1D, Q1D, NE);

  // One thread per element. Each component of u has one open direction (D1D - 1 basis
  // functions) and two closed directions (D1D); the components are stored one after another in
  // the E-vector, which osc offsets into.
  mfem::forall(NE,
               [=] MFEM_HOST_DEVICE(int e)
               {
                 constexpr int MQ1D = mfem::DofQuadLimits::HCURL_MAX_Q1D;

                 for (int qz = 0; qz < Q1D; ++qz)
                   for (int qy = 0; qy < Q1D; ++qy)
                     for (int qx = 0; qx < Q1D; ++qx)
                       for (int c = 0; c < 3; ++c)
                         C(c, qx, qy, qz, e) = 0.0;

                 int osc = 0;

                 {
                   // x component: the reference curl is [0, (u_0)_{x_2}, -(u_0)_{x_1}]
                   const int D1Dz = D1D;
                   const int D1Dy = D1D;
                   const int D1Dx = D1D - 1;

                   for (int dz = 0; dz < D1Dz; ++dz)
                   {
                     mfem::real_t gradXY[MQ1D][MQ1D][2];
                     for (int qy = 0; qy < Q1D; ++qy)
                       for (int qx = 0; qx < Q1D; ++qx)
                         gradXY[qy][qx][0] = gradXY[qy][qx][1] = 0.0;

                     for (int dy = 0; dy < D1Dy; ++dy)
                     {
                       mfem::real_t massX[MQ1D];
                       for (int qx = 0; qx < Q1D; ++qx)
                         massX[qx] = 0.0;

                       for (int dx = 0; dx < D1Dx; ++dx)
                       {
                         const mfem::real_t t = X(dx + ((dy + (dz * D1Dy)) * D1Dx) + osc, e);
                         for (int qx = 0; qx < Q1D; ++qx)
                           massX[qx] += t * Bo(qx, dx);
                       }

                       for (int qy = 0; qy < Q1D; ++qy)
                       {
                         const mfem::real_t wy = Bc(qy, dy);
                         const mfem::real_t wDy = Gc(qy, dy);
                         for (int qx = 0; qx < Q1D; ++qx)
                         {
                           gradXY[qy][qx][0] += massX[qx] * wDy;
                           gradXY[qy][qx][1] += massX[qx] * wy;
                         }
                       }
                     }

                     for (int qz = 0; qz < Q1D; ++qz)
                     {
                       const mfem::real_t wz = Bc(qz, dz);
                       const mfem::real_t wDz = Gc(qz, dz);
                       for (int qy = 0; qy < Q1D; ++qy)
                         for (int qx = 0; qx < Q1D; ++qx)
                         {
                           C(1, qx, qy, qz, e) += gradXY[qy][qx][1] * wDz; // (u_0)_{x_2}
                           C(2, qx, qy, qz, e) -= gradXY[qy][qx][0] * wz;  // -(u_0)_{x_1}
                         }
                     }
                   }

                   osc += D1Dx * D1Dy * D1Dz;
                 }

                 {
                   // y component: the reference curl is [-(u_1)_{x_2}, 0, (u_1)_{x_0}]
                   const int D1Dz = D1D;
                   const int D1Dy = D1D - 1;
                   const int D1Dx = D1D;

                   for (int dz = 0; dz < D1Dz; ++dz)
                   {
                     mfem::real_t gradXY[MQ1D][MQ1D][2];
                     for (int qy = 0; qy < Q1D; ++qy)
                       for (int qx = 0; qx < Q1D; ++qx)
                         gradXY[qy][qx][0] = gradXY[qy][qx][1] = 0.0;

                     for (int dx = 0; dx < D1Dx; ++dx)
                     {
                       mfem::real_t massY[MQ1D];
                       for (int qy = 0; qy < Q1D; ++qy)
                         massY[qy] = 0.0;

                       for (int dy = 0; dy < D1Dy; ++dy)
                       {
                         const mfem::real_t t = X(dx + ((dy + (dz * D1Dy)) * D1Dx) + osc, e);
                         for (int qy = 0; qy < Q1D; ++qy)
                           massY[qy] += t * Bo(qy, dy);
                       }

                       for (int qx = 0; qx < Q1D; ++qx)
                       {
                         const mfem::real_t wx = Bc(qx, dx);
                         const mfem::real_t wDx = Gc(qx, dx);
                         for (int qy = 0; qy < Q1D; ++qy)
                         {
                           gradXY[qy][qx][0] += wDx * massY[qy];
                           gradXY[qy][qx][1] += wx * massY[qy];
                         }
                       }
                     }

                     for (int qz = 0; qz < Q1D; ++qz)
                     {
                       const mfem::real_t wz = Bc(qz, dz);
                       const mfem::real_t wDz = Gc(qz, dz);
                       for (int qy = 0; qy < Q1D; ++qy)
                         for (int qx = 0; qx < Q1D; ++qx)
                         {
                           C(0, qx, qy, qz, e) -= gradXY[qy][qx][1] * wDz; // -(u_1)_{x_2}
                           C(2, qx, qy, qz, e) += gradXY[qy][qx][0] * wz;  // (u_1)_{x_0}
                         }
                     }
                   }

                   osc += D1Dx * D1Dy * D1Dz;
                 }

                 {
                   // z component: the reference curl is [(u_2)_{x_1}, -(u_2)_{x_0}, 0]
                   const int D1Dz = D1D - 1;
                   const int D1Dy = D1D;
                   const int D1Dx = D1D;

                   for (int dx = 0; dx < D1Dx; ++dx)
                   {
                     mfem::real_t gradYZ[MQ1D][MQ1D][2];
                     for (int qz = 0; qz < Q1D; ++qz)
                       for (int qy = 0; qy < Q1D; ++qy)
                         gradYZ[qz][qy][0] = gradYZ[qz][qy][1] = 0.0;

                     for (int dy = 0; dy < D1Dy; ++dy)
                     {
                       mfem::real_t massZ[MQ1D];
                       for (int qz = 0; qz < Q1D; ++qz)
                         massZ[qz] = 0.0;

                       for (int dz = 0; dz < D1Dz; ++dz)
                       {
                         const mfem::real_t t = X(dx + ((dy + (dz * D1Dy)) * D1Dx) + osc, e);
                         for (int qz = 0; qz < Q1D; ++qz)
                           massZ[qz] += t * Bo(qz, dz);
                       }

                       for (int qy = 0; qy < Q1D; ++qy)
                       {
                         const mfem::real_t wy = Bc(qy, dy);
                         const mfem::real_t wDy = Gc(qy, dy);
                         for (int qz = 0; qz < Q1D; ++qz)
                         {
                           gradYZ[qz][qy][0] += massZ[qz] * wy;
                           gradYZ[qz][qy][1] += massZ[qz] * wDy;
                         }
                       }
                     }

                     for (int qx = 0; qx < Q1D; ++qx)
                     {
                       const mfem::real_t wx = Bc(qx, dx);
                       const mfem::real_t wDx = Gc(qx, dx);
                       for (int qy = 0; qy < Q1D; ++qy)
                         for (int qz = 0; qz < Q1D; ++qz)
                         {
                           C(0, qx, qy, qz, e) += gradYZ[qz][qy][1] * wx;  // (u_2)_{x_1}
                           C(1, qx, qy, qz, e) -= gradYZ[qz][qy][0] * wDx; // -(u_2)_{x_0}
                         }
                     }
                   }
                 }
               });
}

void
NLCurlCurlPhysicalCurl3D(const int Q1D, const int NE, const mfem::Vector & j, mfem::Vector & curl)
{
  const int NQ = Q1D * Q1D * Q1D;
  auto J = mfem::Reshape(j.Read(), NQ, 3, 3, NE);
  auto C = mfem::Reshape(curl.ReadWrite(), 3, NQ, NE);

  mfem::forall(NE * NQ,
               [=] MFEM_HOST_DEVICE(int idx)
               {
                 const int q = idx % NQ;
                 const int e = idx / NQ;
                 const mfem::real_t J11 = J(q, 0, 0, e);
                 const mfem::real_t J21 = J(q, 1, 0, e);
                 const mfem::real_t J31 = J(q, 2, 0, e);
                 const mfem::real_t J12 = J(q, 0, 1, e);
                 const mfem::real_t J22 = J(q, 1, 1, e);
                 const mfem::real_t J32 = J(q, 2, 1, e);
                 const mfem::real_t J13 = J(q, 0, 2, e);
                 const mfem::real_t J23 = J(q, 1, 2, e);
                 const mfem::real_t J33 = J(q, 2, 2, e);
                 const mfem::real_t inv_detJ =
                     1.0 / (J11 * (J22 * J33 - J32 * J23) - J21 * (J12 * J33 - J32 * J13) +
                            J31 * (J12 * J23 - J22 * J13));

                 const mfem::real_t c1 = C(0, q, e);
                 const mfem::real_t c2 = C(1, q, e);
                 const mfem::real_t c3 = C(2, q, e);
                 C(0, q, e) = inv_detJ * (J11 * c1 + J12 * c2 + J13 * c3);
                 C(1, q, e) = inv_detJ * (J21 * c1 + J22 * c2 + J23 * c3);
                 C(2, q, e) = inv_detJ * (J31 * c1 + J32 * c2 + J33 * c3);
               });
}

void
NLCurlCurlGradPASetup(const int Q1D,
                      const int NE,
                      const mfem::Array<mfem::real_t> & w,
                      const mfem::Vector & j,
                      const mfem::Vector & c,
                      mfem::Vector & k_coeff,
                      mfem::Vector & dk_coeff,
                      mfem::Vector & op,
                      mfem::real_t curlu_zero_tol)
{

  // number of quadpoints per element total
  const int NQ = Q1D * Q1D * Q1D;

  // next, turn all the important stuff into tensors
  auto W = w.Read();
  auto J = mfem::Reshape(j.Read(), NQ, 3, 3, NE);
  // hardcoding 1 here, since k expected to be scalar function
  auto K = mfem::Reshape(k_coeff.Read(), 1, NQ, NE);
  auto DK = mfem::Reshape(dk_coeff.Read(), 1, NQ, NE);

  // the Curl CoefficientVector
  auto C = mfem::Reshape(c.Read(), 3, NQ, NE);

  // finally, our operator
  auto Diag = mfem::Reshape(op.Write(), NQ, 6, NE);

  mfem::forall(NE * NQ,
               [=] MFEM_HOST_DEVICE(int idx)
               {
                 const int q = idx % NQ;
                 const int e = idx / NQ;
                 const mfem::real_t J11 = J(q, 0, 0, e);
                 const mfem::real_t J21 = J(q, 1, 0, e);
                 const mfem::real_t J31 = J(q, 2, 0, e);
                 const mfem::real_t J12 = J(q, 0, 1, e);
                 const mfem::real_t J22 = J(q, 1, 1, e);
                 const mfem::real_t J32 = J(q, 2, 1, e);
                 const mfem::real_t J13 = J(q, 0, 2, e);
                 const mfem::real_t J23 = J(q, 1, 2, e);
                 const mfem::real_t J33 = J(q, 2, 2, e);
                 const mfem::real_t detJ = J11 * (J22 * J33 - J32 * J23) -
                                           J21 * (J12 * J33 - J32 * J13) +
                                           J31 * (J12 * J23 - J22 * J13);
                 const mfem::real_t beta = W[q] / detJ;

                 // next, compute alpha. Simply read off the value of DK from this
                 // element/qpoint
                 const mfem::real_t alpha = DK(0, q, e);

                 // ditto for k
                 const mfem::real_t k = K(0, q, e);

                 // load the c_i - elements of the curl vector
                 mfem::real_t c1 = C(0, q, e);
                 mfem::real_t c2 = C(1, q, e);
                 mfem::real_t c3 = C(2, q, e);

                 // Unit curl vector
                 const mfem::real_t curl_mag = std::sqrt(c1 * c1 + c2 * c2 + c3 * c3);
                 const mfem::real_t inv_curl_mag = curl_mag > curlu_zero_tol ? 1.0 / curl_mag : 0.0;
                 c1 *= inv_curl_mag;
                 c2 *= inv_curl_mag;
                 c3 *= inv_curl_mag;

                 // next compute the g_i
                 const mfem::real_t g1 = J11 * c1 + J21 * c2 + J31 * c3;
                 const mfem::real_t g2 = J12 * c1 + J22 * c2 + J32 * c3;
                 const mfem::real_t g3 = J13 * c1 + J23 * c2 + J33 * c3;

                 // next the M_i. TODO - we can read these off the AddMultPA array
                 const mfem::real_t M11 = J11 * J11 + J21 * J21 + J31 * J31;
                 const mfem::real_t M22 = J12 * J12 + J22 * J22 + J32 * J32;
                 const mfem::real_t M33 = J13 * J13 + J23 * J23 + J33 * J33;
                 const mfem::real_t M12 = J11 * J12 + J21 * J22 + J31 * J32;
                 const mfem::real_t M13 = J11 * J13 + J21 * J23 + J31 * J33;
                 const mfem::real_t M23 = J12 * J13 + J22 * J23 + J32 * J33;

                 // Finally, let's write the 6 elements to the Diag tensor
                 // Do the upper triangle, in row major order
                 Diag(q, 0, e) = beta * (k * M11 + alpha * g1 * g1); // D11
                 Diag(q, 1, e) = beta * (k * M12 + alpha * g1 * g2); // D12
                 Diag(q, 2, e) = beta * (k * M13 + alpha * g1 * g3); // D13
                 Diag(q, 3, e) = beta * (k * M22 + alpha * g2 * g2); // D22
                 Diag(q, 4, e) = beta * (k * M23 + alpha * g2 * g3); // D23
                 Diag(q, 5, e) = beta * (k * M33 + alpha * g3 * g3); // D33
               });
}

void
NLCurlCurlPASetup(const int Q1D,
                  const int NE,
                  const mfem::Array<mfem::real_t> & w,
                  const mfem::Vector & j,
                  mfem::Vector & k_coeff,
                  mfem::Vector & op)
{
  // number of quadpoints per element total
  const int NQ = Q1D * Q1D * Q1D;

  // next, turn all the important stuff into tensors
  auto W = w.Read();
  auto J = mfem::Reshape(j.Read(), NQ, 3, 3, NE);
  // hardcoding 1 here, since k expected to be scalar function
  auto K = mfem::Reshape(k_coeff.Read(), 1, NQ, NE);

  // finally, our operator
  auto Diag = mfem::Reshape(op.Write(), NQ, 6, NE);

  // for each element
  for (int e = 0; e < NE; e++)
  {
    // for each qpoint
    for (int q = 0; q < NQ; q++)
    {
      const mfem::real_t J11 = J(q, 0, 0, e);
      const mfem::real_t J21 = J(q, 1, 0, e);
      const mfem::real_t J31 = J(q, 2, 0, e);
      const mfem::real_t J12 = J(q, 0, 1, e);
      const mfem::real_t J22 = J(q, 1, 1, e);
      const mfem::real_t J32 = J(q, 2, 1, e);
      const mfem::real_t J13 = J(q, 0, 2, e);
      const mfem::real_t J23 = J(q, 1, 2, e);
      const mfem::real_t J33 = J(q, 2, 2, e);
      const mfem::real_t detJ = J11 * (J22 * J33 - J32 * J23) - J21 * (J12 * J33 - J32 * J13) +
                                J31 * (J12 * J23 - J22 * J13);
      const mfem::real_t beta = W[q] / detJ;

      const mfem::real_t M11 = J11 * J11 + J21 * J21 + J31 * J31;
      const mfem::real_t M22 = J12 * J12 + J22 * J22 + J32 * J32;
      const mfem::real_t M33 = J13 * J13 + J23 * J23 + J33 * J33;
      const mfem::real_t M12 = J11 * J12 + J21 * J22 + J31 * J32;
      const mfem::real_t M13 = J11 * J13 + J21 * J23 + J31 * J33;
      const mfem::real_t M23 = J12 * J13 + J22 * J23 + J32 * J33;

      const mfem::real_t k = K(0, q, e);

      // finally, write into the diagonal operator
      Diag(q, 0, e) = beta * k * M11; // D11
      Diag(q, 1, e) = beta * k * M12; // D12
      Diag(q, 2, e) = beta * k * M13; // D13
      Diag(q, 3, e) = beta * k * M22; // D22
      Diag(q, 4, e) = beta * k * M23; // D23
      Diag(q, 5, e) = beta * k * M33; // D33
    }
  }
}
} // namespace Moose::MFEM
