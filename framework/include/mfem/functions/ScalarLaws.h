#ifdef MOOSE_MFEM_ENABLED

#pragma once

#include "libmesh/ignore_warnings.h"
#include "mfem.hpp"
#include "libmesh/restore_warnings.h"

namespace Moose::MFEM
{
/**
 * k(s) = a + b s^n. Plain data, copied by value into device kernels.
 */
struct PowerLaw
{
  mfem::real_t a;
  mfem::real_t b;
  mfem::real_t n;

  /// Evaluate k(s) and s dk/ds
  MFEM_HOST_DEVICE inline void
  eval(const mfem::real_t s, mfem::real_t & k, mfem::real_t & s_dk_ds) const;
};

// We are writing it this way to avoid virtual functions in device code.
MFEM_HOST_DEVICE inline void
PowerLaw::eval(const mfem::real_t s, mfem::real_t & k, mfem::real_t & s_dk_ds) const
{
  // pow(0, 0) == 1, so n = 0 gives the constant a + b
  const mfem::real_t b_sn = b * pow(s, n);
  k = a + b_sn;
  s_dk_ds = n * b_sn;
}

} // namespace Moose::MFEM

#endif
