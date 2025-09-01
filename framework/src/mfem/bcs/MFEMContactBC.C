//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMContactBC.h"

registerMooseObject("MooseApp", MFEMContactBC);

InputParameters
MFEMContactBC::validParams()
{
  return MFEMEssentialBC::validParams();
}

MFEMContactBC::MFEMContactBC(const InputParameters & parameters)
  : MFEMEssentialBC(parameters)
{
}

void
MFEMContactBC::ApplyBC(mfem::GridFunction & /*gridfunc*/)
{
}

#endif
