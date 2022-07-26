[GlobalParams]
  var_name_base = {var_name_base}
  op_num = {op_num}
  #en_ratio = 1
[]

[Mesh]
  type = GeneratedMesh
  dim = {dim}
  nx = {nx}
  ny = {ny}
  nz = {nz}
  xmin = {min_0}
  xmax = {max_0}
  ymin = {min_1}
  ymax = {max_1}
  zmin = {min_2}
  zmax = {max_2}
  uniform_refine = {h_level}
  elem_type = QUAD4
[]

[Variables]
  [./c]
    #scaling = 10
  [../]
  [./w]
  [../]
  [./PolycrystalVariables]
  [../]
[]

[AuxVariables]
  [./bnds]
  [../]
  [./total_en]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./unique_grains]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./var_indices]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./centroids]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./cres]
    type = SplitCHParsed
    variable = c
    kappa_name = kappa_c
    w = w
    f_name = F
    args = '{gr_args_list}'
  [../]
  [./wres]
    type = SplitCHWRes
    variable = w
    mob_name = D
  [../]
  [./time]
    type = CoupledTimeDerivative
    variable = w
    v = c
  [../]
  [./PolycrystalSinteringKernel]
    c = c
    consider_rigidbodymotion = false
  [../]
[]

[AuxKernels]
  [./bnds]
    type = BndsCalcAux
    variable = bnds
    v = '{gr_args_list}'
  [../]
  [./Total_en]
    type = TotalFreeEnergy
    variable = total_en
    kappa_names = 'kappa_c {gr_kappa_list}'
    interfacial_vars = 'c {gr_args_list}'
  [../]
[]

[Materials]
  [./free_energy]
    type = SinteringFreeEnergy
    block = 0
    c = c
    v = '{gr_args_list}'
    derivative_order = 2
    #outputs = console
  [../]
  [./CH_mat]
    type = PFDiffusionGrowth
    block = 0
    rho = c
    v = '{gr_args_list}'
    outputs = console
    Dvol = {Dvol}
    Dvap = {Dvap}
    Dsurf = {Dsurf}
    Dgb = {Dgb}
  [../]
  [./constant_mat]
    type = GenericConstantMaterial
    block = 0
    prop_names = '  A    B   L   kappa_op kappa_c '
    prop_values = '{energy_A} {energy_B} {energy_L} {energy_kappa_op} {energy_kappa_c}'
  [../]
[]

[Postprocessors]
  [./elem_c]
    type = ElementIntegralVariablePostprocessor
    variable = c
  [../]
  [./elem_bnds]
    type = ElementIntegralVariablePostprocessor
    variable = bnds
  [../]
  [./total_energy]
    type = ElementIntegralVariablePostprocessor
    variable = total_en
  [../]
  [./free_en]
    type = ElementIntegralMaterialProperty
    mat_prop = F
  [../]
  [./dofs]
    type = NumDOFs
    system = 'NL'
  [../]
  [./tstep]
    type = TimestepSize
  [../]
  [./int_area]
    type = InterfaceAreaPostprocessor
    variable = c
  [../]
  [./gb_area]
    type = GrainBoundaryArea
  [../]
  [./neck]
    type = NeckAreaPostprocessor
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    coupled_groups = 'c,w c,{op_coupling_list} '
  [../]
[]

[Executioner]
  # Preconditioned JFNK (default)
  type = Transient
  scheme = BDF1
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -ksp_grmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  petsc_options_value = 'asm         31   preonly   lu      1'
  l_max_its = 20
  nl_max_its = 20
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-08
  l_tol = 1e-04
  end_time = {t_end}
  #dt = 0.01
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    growth_factor = 1.2
  [../]
[]

[Adaptivity]
  marker = bound_adapt
  max_h_level = {max_h_level}
  [./Indicators]
    [./error]
      type = GradientJumpIndicator
      variable = bnds
    [../]
  [../]
  [./Markers]
    [./bound_adapt]
      type = ValueRangeMarker
      lower_bound = 0.01
      upper_bound = 0.99
      variable = bnds
    [../]
  [../]
[]

[Outputs]
  print_linear_residuals = true
  csv = true
  gnuplot = true
  print_perf_log = true
  [./console]
    type = Console
    perf_log = true
  [../]
  [./exodus]
    type = Exodus
    elemental_as_nodal = true
  [../]
[]

[Debug]
  show_var_residual_norms = true
[]

[UserObjects]
  [./grain_center]
    type = GrainTracker
    outputs = none
    compute_var_to_feature_map = true
    execute_on = 'initial timestep_begin'
  [../]
[]

