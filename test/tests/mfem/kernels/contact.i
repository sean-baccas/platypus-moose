[Mesh]
  type = MFEMMesh
  file = ../mesh/two-hex.mesh
  dim = 3
  uniform_refine = 2
  displacement = 'displacement'
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
  [H1FESpace]
    type = MFEMVectorFESpace
    fec_type = H1
    fec_order = FIRST
    ordering = NODES
  []
[]

[Variables]
  [displacement]
    type = MFEMVariable
    fespace = H1FESpace
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

[FunctorMaterials]
  [Rigidium]
    type = MFEMGenericFunctorMaterial
    prop_names = 'lambda mu'
    prop_values = '50.0 50.0'
  []
[]

[BCs]
  [interface]
    type = MFEMContactBC
    variable = displacement
  []
  [dirichlet]
    type = MFEMVectorDirichletBC
    variable = displacement
    boundary = "1 2 3 6"
    vector_coefficient = '0.0 0.0 0.0'
  []
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[Solver]
  type = MFEMMinResSolver
  preconditioner = hyprediag
  rel_tol   = 1e-12
  l_max_its = 1500  
[]

[Preconditioner]
  [hyprediag]
    type = MFEMHypreDiagScale
  []
[]

[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/Contact
    vtk_format = ASCII
  []
[]


