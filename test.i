[Mesh]
    [mesh]
      type = FileMeshGenerator
      file = '3D.exo'
    []
    
    [left]
        type = ParsedGenerateSideset
        combinatorial_geometry = 'x=0'
        new_sideset_name = 'left'
        input = 'mesh'
      []
      [right]
        type = ParsedGenerateSideset
        combinatorial_geometry = 'x=100'
        new_sideset_name = 'right'
        input = 'left'
      []
[]

  [GlobalParams]
    PorousFlowDictator = dictator
  []
  
  [Variables]
    [porepressure]
      initial_condition = 20E6
    []
  []
  
  [BCs]
    [left]
      type = DirichletBC
      variable = porepressure
      value = 20E6
      boundary = 'left'
    []
    [right]
        type = DirichletBC
        variable = porepressure
        value = 15E6
        boundary = 'right'
      []
  []

  
  [FluidProperties]
    [the_simple_fluid]
      type = SimpleFluidProperties
      thermal_expansion = 2E-4
      bulk_modulus = 2E9
      viscosity = 1E-3
      density0 = 1000
      cv = 4000.0
      cp = 4000.0
    []
  []
  
  [PorousFlowUnsaturated]
    porepressure = porepressure
    coupling_type = Hydro
    gravity = '0 0 0'
    fp = the_simple_fluid
  []
  
  [Materials]
    [porosity]
      type = PorousFlowPorosityConst # only the initial value of this is ever used
      porosity = 0.1
    []
    [biot_modulus]
      type = PorousFlowConstantBiotModulus
      solid_bulk_compliance = 1E-10
      fluid_bulk_modulus = 2E9
    []
    [permeability]
      type = PorousFlowPermeabilityConst
      permeability = '1E-12 0 0   0 1E-12 0   0 0 1E-12'
    []
    [thermal_expansion]
      type = PorousFlowConstantThermalExpansionCoefficient
      fluid_coefficient = 5E-6
      drained_coefficient = 2E-4
    []
    [thermal_conductivity]
      type = PorousFlowThermalConductivityIdeal
      dry_thermal_conductivity = '1 0 0  0 1 0  0 0 1'
    []
    [rock_heat]
      type = PorousFlowMatrixInternalEnergy
      density = 2500.0
      specific_heat_capacity = 1200.0
    []
  []
  
  [Preconditioning]
    active = basic
    [basic]
      type = SMP
      full = true
      petsc_options = '-ksp_diagonal_scale -ksp_diagonal_scale_fix'
      petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
      petsc_options_value = ' asm      lu           NONZERO                   2'
    []
    [preferred_but_might_not_be_installed]
      type = SMP
      full = true
      petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
      petsc_options_value = ' lu       mumps'
    []
  []
  
  [Executioner]
    type = Steady
    solve_type = Newton

  []
  
  [Outputs]
  file_base= '3D'
    exodus = true
  []
