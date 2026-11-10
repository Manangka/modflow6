bmi_sources = files(
    'bmi.f90',
    'mf6bmi.f90',
    'mf6bmiError.f90',
    'mf6bmiGrid.f90',
    'mf6bmiUtil.f90',
    'mf6xmi.F90',
)

message('MODFLOW 6 shared library name:', 'lib' + buildname)

# Only the public C API, i.e. all bind(C, name=...) entities in the sources
# above, is exported from the shared library. The export list is generated in
# the format of the platform's linker.
export_defines = []
if is_parallel_build
  export_defines += ['--define', '__WITH_MPI__']
endif
if build_machine.system() == 'windows'
  export_format = 'def'
elif build_machine.system() == 'darwin'
  export_format = 'exports'
else
  export_format = 'map'
endif
export_list = custom_target(
  'libmf6_exports',
  input: bmi_sources,
  output: 'libmf6.' + export_format,
  command: [
    import('python').find_installation(),
    files('generate_exports.py'),
    '--format', export_format,
    '--output', '@OUTPUT@',
    export_defines,
    '@INPUT@',
  ],
)

bmi_link_args = []
bmi_kwargs = {}
if export_format == 'def'
  bmi_kwargs += {'vs_module_defs': export_list}
elif export_format == 'exports'
  bmi_link_args += '-Wl,-exported_symbols_list,' + export_list.full_path()
else
  bmi_link_args += '-Wl,--version-script=' + export_list.full_path()
endif

library(
  buildname,
  bmi_sources,
  link_with: mf6core,
  name_prefix: 'lib',
  link_args: bmi_link_args,
  link_depends: export_list,
  install: true,
  kwargs: bmi_kwargs,
)
