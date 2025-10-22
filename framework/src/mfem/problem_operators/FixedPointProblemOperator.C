
#ifdef MOOSE_MFEM_ENABLED

#include "FixedPointProblemOperator.h"

// Does a weighted += into the first from the second
static inline void Combine(mfem::BlockVector& a_vec, mfem::BlockVector& b_vec, float alpha)
{
  auto* a_ptr = a_vec.ReadWrite(false);
  auto* b_ptr = b_vec.ReadWrite(false);
  auto  size  = a_vec.Size();

  // write everything into old_vec
  for (int i=0; i<size; i++)
  {
    // could give it the option to offload this one
    a_ptr[i] = alpha*a_ptr[i] + (1.0-alpha)*b_ptr[i];
  }
}

namespace Moose::MFEM
{
void
FixedPointProblemOperator::SetGridFunctions()
{
  _test_var_names = GetEquationSystem()->TestVarNames();
  _trial_var_names = GetEquationSystem()->TrialVarNames();
  ProblemOperator::SetGridFunctions();

  // Make sure the sizes are set correctly here
  _true_x_old.Update(_block_true_offsets);

  // we also need to work out the block in which the displacement gridfunction data sits
  if (_problem.mesh().shouldDisplace())
  {
    _use_displacement = true;

    // we know this will return something
    std::string displacement_var = _problem.mesh().getMeshDisplacementVariable().value();

    // look it up in the list of gridfunctions
    // int displacement_index = 0;
    auto iter = _problem_data.gridfunctions.begin();
    auto end_iter = _problem_data.gridfunctions.end();

    for (; iter!=end_iter; iter++, _displacement_index++)
    {
      // if true, we've found which block we expect to find
      // the displacement variable
      if (iter->first==displacement_var) break;
    }

  }
}

void
FixedPointProblemOperator::Init(mfem::BlockVector & X)
{
  ProblemOperator::Init(X);

  GetEquationSystem()->BuildEquationSystem();
}

void
FixedPointProblemOperator::Solve()
{
  _problem.initMeshDisplacement();
  
  // we don't have to set this stuff to zero for correctness
  _true_x_old = 0.0;

  // make sure the old displacement is initialised
  for (int i=0; i<_max_iters; i++)
  {
    _problem.displaceMesh();

    GetEquationSystem()->BuildEquationSystem();
    GetEquationSystem()->BuildJacobian(_true_x, _true_rhs);
  
    if (_problem_data.jacobian_solver->isLOR() && _equation_system->_test_var_names.size() > 1)
      mooseError("LOR solve is only supported for single-variable systems");
  
    _problem_data.jacobian_solver->updateSolver(
        *_equation_system->_blfs.Get(_equation_system->_test_var_names.at(0)),
        _equation_system->_ess_tdof_lists.at(0));
  
    _problem_data.nonlinear_solver->SetSolver(_problem_data.jacobian_solver->getSolver());
    _problem_data.nonlinear_solver->SetOperator(*GetEquationSystem());
    _problem_data.nonlinear_solver->Mult(_true_rhs, _true_x);

    if (_use_displacement)
    {
      // we need to update the "old_displacement" in the mfem mesh
      
      _true_x.GetBlock(_displacement_index).SyncAliasMemory(_true_x);
      _problem.mesh().getOldDisplacement() = _true_x_old.GetBlock(_displacement_index);
    }

    Combine(_true_x, _true_x_old, _alpha);
  
    // GetEquationSystem()->RecoverFEMSolution(_true_x_old, _problem_data.gridfunctions);
    GetEquationSystem()->RecoverFEMSolution(_true_x, _problem_data.gridfunctions);

    if (!ContinueToIterate()) {
      std::cout << "Exiting after " << (i+1) << " iterations\n";
      break;
    }

    // stash the old ones
    _true_x_old = _true_x;
  }
}

bool
FixedPointProblemOperator::ContinueToIterate() const {
  // sum the differences
  mfem::real_t* old_x = _true_x_old.GetData();
  mfem::real_t* new_x = _true_x.GetData();

  mfem::real_t sum = 0;
  // This loop is brutally slow
  for (int i=0; i<_true_x.Size(); i++)
  {
    sum += (new_x[i] - old_x[i]) * (new_x[i] - old_x[i]);
  }

  std::cout << "Square difference is " << sum << "\n";

  return sum > _tol;
}

} // namespace Moose::MFEM

#endif
