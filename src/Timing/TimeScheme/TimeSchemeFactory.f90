module TimeSchemeFactoryModule

  use KindModule, only: I4B, DP
  use SimModule, only: store_error
  use TdisModule, only: delt_buffer, kper, kstp

  use TimeSchemeEnumModule
  use TimeSchemeInterfaceModule, only: TimeSchemeInterface
  use ImplicitEulerSchemeModule, only: ImplicitEulerSchemeType
  use BDF2SchemeModule, only: BDF2SchemeType

  implicit none
  private
  public :: create_time_scheme

contains

  function create_time_scheme(time_type) result(time_scheme)
    ! -- result
    class(TimeSchemeInterface), pointer :: time_scheme !< allocated concrete time scheme or null()
    ! -- dummy
    integer(I4B), intent(in) :: time_type !< enumerator from `TimeSchemeEnumModule`
    !
    ! -- Allocate time scheme instance
    select case (time_type)
    case (TIME_SCHEME_EULER)
      allocate (time_scheme, source=ImplicitEulerSchemeType(delt_buffer))
    case (TIME_SCHEME_BDF2)
      allocate (time_scheme, source=BDF2SchemeType(delt_buffer, kstp, kper))
    case default
      call store_error("Unknown time scheme", terminate=.TRUE.)
    end select

  end function create_time_scheme

end module TimeSchemeFactoryModule
