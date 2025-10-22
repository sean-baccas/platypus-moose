#ifdef MOOSE_MFEM_ENABLED

#pragma once
#include "ProblemOperator.h"
#include "EquationSystemInterface.h"


namespace Moose::MFEM
{
class FixedPointProblemOperator : public ProblemOperator, public EquationSystemInterface
{
public:
  FixedPointProblemOperator(MFEMProblem& problem)
  : ProblemOperator(problem), _equation_system(_problem_data.eqn_system) {}

  void SetGridFunctions() override;
  void Init(mfem::BlockVector & X) override;
  virtual void Solve() override;

  [[nodiscard]] Moose::MFEM::EquationSystem * GetEquationSystem() const override
  {
    if (!_equation_system)
    {
      MFEM_ABORT("No equation system has been added to ProblemOperator.");
    }

    return _equation_system.get();
  }

  // Check that the _true_x current is not changing anymore
  bool ContinueToIterate() const;


private:
  std::shared_ptr<Moose::MFEM::EquationSystem> _equation_system{nullptr};

  // Weighting for previous iteration vs current iteration
  const float _alpha = 1.0;
  const float _tol   = 1e-3;
  const int   _max_iters = 20;

  // index within 
  int _displacement_index = 0;

  bool _use_displacement = false;

  // KEEP SOME COPY OF THE GRIDFUNCTIONS HERE!!
  mfem::BlockVector _true_x_old;
};


} // namespace Moose::MFEM






#endif
