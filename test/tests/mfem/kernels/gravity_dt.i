[Mesh]
  type = MFEMMesh
  file = ../mesh/beam_coarse.e
  displacement = "displacement"
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
  [H1FESpace]
    type = MFEMVectorFESpace
    fec_type = H1
    fec_order = FIRST
    ordering = "vdim"
  []
[]

[Variables]
  [displacement]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[BCs]
  [dirichlet]
    type = MFEMVectorDirichletBC
    variable = displacement
    boundary = '1'
  []
  [pull_down]
    type = MFEMVectorBoundaryIntegratedBC
    variable = displacement
    boundary = '2'
    vector_coefficient = '10.0 0.0 0.0'
  []
[]

[FunctorMaterials]
  [Rigidium]
    type = MFEMGenericFunctorMaterial
    prop_names = 'lambda mu'
    prop_values = '0 50'
  []
[]

[Kernels]
  [diff]
    type = MFEMLinearElasticityKernel
    variable = displacement
    lambda = lambda
    mu = mu
  []
[]

[Preconditioner]
  [superlu]
    type = MFEMHypreBoomerAMG
    fespace = H1FESpace
    l_max_its = 20
    print_level = 2
  []
[]

[Solver]
  type = MFEMSuperLU
[]

[Executioner]
  type = MFEMSteady
  device = "cpu"
  assembly_level = legacy
[]


[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/Gravity_dt
    vtk_format = ASCII
  []
[]
 