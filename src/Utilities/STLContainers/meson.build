fypp_files = [
  'STLVecInt.fypp',
  'STLVecDbl.fypp',
]

generated_sources = []

foreach f : fypp_files
  generated_sources += custom_target(
    f.replace('.fypp', '') + '_gen',
    input  : f,
    output : f.replace('.fypp', '.f90'),
    command: [
      fypp,
      '-I' + meson.project_source_root() / 'src' / 'include',
      '@INPUT@',
      '@OUTPUT@'
    ]
  )
endforeach
