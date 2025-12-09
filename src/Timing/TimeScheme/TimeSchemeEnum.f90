module TimeSchemeEnumModule
  use KindModule, only: I4B

  implicit none

  ! Time scheme codes
  integer(I4B), parameter :: TIME_SCHEME_EULER = 0
  integer(I4B), parameter :: TIME_SCHEME_BDF2 = 1

end module TimeSchemeEnumModule
