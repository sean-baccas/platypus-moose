//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMSteady.h"
#include "MFEMProblem.h"
#include "EquationSystemProblemOperator.h"

#include "axom/slic.hpp"
#include "tribol/interface/tribol.hpp"
#include "tribol/interface/mfem_tribol.hpp"

registerMooseObject("MooseApp", MFEMSteady);

InputParameters
MFEMSteady::validParams()
{
  InputParameters params = MFEMProblemSolve::validParams();
  params += Executioner::validParams();
  params.addClassDescription("Executioner for steady state MFEM problems.");
  params.addParam<Real>("time", 0.0, "System time");
  return params;
}

MFEMSteady::MFEMSteady(const InputParameters & params)
  : Executioner(params),
    _mfem_problem(dynamic_cast<MFEMProblem &>(feProblem())),
    _mfem_problem_data(_mfem_problem.getProblemData()),
    _mfem_problem_solve(*this, getProblemOperators()),
    _system_time(getParam<Real>("time")),
    _time_step(_mfem_problem.timeStep()),
    _time(_mfem_problem.time()),
    _last_solve_converged(false)
{
  _time = _system_time;
  // If no ProblemOperators have been added by the user, add a default
  if (getProblemOperators().empty())
  {
    _mfem_problem_data.eqn_system = std::make_shared<Moose::MFEM::EquationSystem>();
    auto problem_operator =
        std::make_shared<Moose::MFEM::EquationSystemProblemOperator>(_mfem_problem);
    addProblemOperator(std::move(problem_operator));
  }
}

void
MFEMSteady::init()
{
  _mfem_problem.execute(EXEC_PRE_MULTIAPP_SETUP);
  _mfem_problem.initialSetup();

  initTribol();

  // Set up initial conditions
  _mfem_problem_data.eqn_system->Init(
      _mfem_problem_data.gridfunctions,
      _mfem_problem_data.fespaces,
      getParam<MooseEnum>("assembly_level").getEnum<mfem::AssemblyLevel>());

  for (const auto & problem_operator : getProblemOperators())
  {
    problem_operator->SetGridFunctions();
    problem_operator->Init(_mfem_problem_data.f);
  }
}

void
MFEMSteady::execute()
{
  if (_app.isRecovering())
  {
    _console << "\nCannot recover steady solves!\nExiting...\n" << std::endl;
    return;
  }

  _time_step = 0;
  _time = _time_step;
  _mfem_problem.outputStep(EXEC_INITIAL);
  _time = _system_time;

  preExecute();

  _mfem_problem.advanceState();

  // first step in any steady state solve is always 1 (preserving backwards compatibility)
  _time_step = 1;
  _mfem_problem.timestepSetup();

  _last_solve_converged = _mfem_problem_solve.solve();

  _mfem_problem.computeIndicators();
  _mfem_problem.computeMarkers();

  // need to keep _time in sync with _time_step to get correct output
  _time = _time_step;
  _mfem_problem.outputStep(EXEC_TIMESTEP_END);
  _time = _system_time;

  {
    TIME_SECTION("final", 1, "Executing Final Objects")
    _mfem_problem.execMultiApps(EXEC_FINAL);
    _mfem_problem.finalizeMultiApps();
    _mfem_problem.postExecute();
    _mfem_problem.execute(EXEC_FINAL);
    _time = _time_step;
    _mfem_problem.outputStep(EXEC_FINAL);
    _time = _system_time;
  }

  postExecute();
}

void
MFEMSteady::initTribol()
{
  // silence warning
  axom::slic::initialize();

  // Create a Tribol coupling scheme: defines contact surfaces and enforcement
  int coupling_scheme_id = 0;
  
  const int dimensions = _mfem_problem.mesh().dimension();
  tribol::initialize(dimensions, MPI_COMM_WORLD);

  // While there is a single mfem ParMesh for this problem, Tribol
  // defines a mortar and a nonmortar contact mesh, each with a unique mesh ID.
  // The Tribol mesh IDs for each contact surface are defined here.
  int mesh1_id = 0;
  int mesh2_id = 1;

  // take a reference to the pmesh
  mfem::ParMesh& pmesh = *(_mfem_problem.getProblemData().pmesh);

  // copied from the test case!
  std::set<int>  mortar_attrs    = {4};
  std::set<int>  nonmortar_attrs = {5};

  tribol::registerMfemCouplingScheme(
    coupling_scheme_id, mesh1_id, mesh2_id,
    pmesh, _mfem_problem.getCoords(), mortar_attrs, nonmortar_attrs,
    tribol::SURFACE_TO_SURFACE,
    tribol::NO_CASE,
    tribol::SINGLE_MORTAR,
    tribol::FRICTIONLESS,
    tribol::LAGRANGE_MULTIPLIER,
    tribol::BINNING_GRID
  );

  // Set Tribol options for Lagrange multiplier enforcement
  tribol::setLagrangeMultiplierOptions(
    coupling_scheme_id,
    tribol::ImplicitEvalMode::MORTAR_RESIDUAL_JACOBIAN
  );
 
  // #4: Update contact mesh decomposition so the on-rank Tribol meshes
  // coincide with the current configuration of the mesh. This must be called
  // before tribol::update().
  tribol::updateMfemParallelDecomposition();
  
  // #5: Update contact gaps, forces, and tangent stiffness contributions
  int cycle = 1;   // pseudo cycle
  mfem::real_t t = 1.0;  // pseudo time
  mfem::real_t dt = 1.0; // pseudo dt
  tribol::update(cycle, t, dt);
}

#endif
