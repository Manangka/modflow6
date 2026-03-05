module BDF2SchemeModule
  use TimeSchemeInterfaceModule, only: TimeSchemeInterface
  use KindModule, only: DP, I4B, LGP
  use CircularBufferModule, only: CircularBufferType
  use ImplicitEulerSchemeModule, only: ImplicitEulerSchemeType

  implicit none
  private

  public :: BDF2SchemeType

  type, extends(TimeSchemeInterface) :: BDF2SchemeType
    private
    class(CircularBufferType), pointer :: delt_buffer => null()
    integer(I4B), pointer :: kstp => null()
    integer(I4B), pointer :: kper => null()
    integer(I4B) :: num_steps = 2 !< number of sub-steps in the time step
  contains
    procedure :: get_num_steps
    procedure :: get_time_iteration_steps
    procedure :: get_weight
  end type BDF2SchemeType

  interface BDF2SchemeType
    module procedure constructor
  end interface BDF2SchemeType

contains
  function constructor(delt_buffer, kstp, kper) Result(scheme)
    type(BDF2SchemeType) :: scheme
    ! -- dummy
    type(CircularBufferType), intent(in), target :: delt_buffer
    integer(I4B), pointer, intent(in) :: kstp
    integer(I4B), pointer, intent(in) :: kper

    scheme%kstp => kstp
    scheme%kper => kper
    scheme%delt_buffer => delt_buffer

  end function constructor

  function get_num_steps(this) result(num_steps)
    class(BDF2SchemeType), intent(in) :: this
    integer(I4B) :: num_steps

    num_steps = this%num_steps
  end function get_num_steps

  function get_time_iteration_steps(this) result(steps)
    ! -- dummy
    class(BDF2SchemeType), intent(in) :: this
    integer(I4B) :: steps
    ! -- local
    logical :: first
    type(ImplicitEulerSchemeType), allocatable :: euler_scheme

    first = this%kstp == 1 .and. this%kper == 1

    if (first) then
      euler_scheme = ImplicitEulerSchemeType(this%delt_buffer)
      steps = euler_scheme%get_time_iteration_steps()
      return
    else
      steps = 2
      return
    end if

  end function get_time_iteration_steps

  function get_weight(this, n) result(weight)
    use SimModule, only: store_error
    ! -- dummy
    class(BDF2SchemeType), intent(in) :: this
    real(DP) :: weight
    integer(I4B), intent(in) :: n
    ! -- local
    logical :: first
    real(DP) :: r
    real(DP), pointer, dimension(:) :: delt, delt_prev
    type(ImplicitEulerSchemeType), allocatable :: euler_scheme

    first = this%kstp == 1 .and. this%kper == 1
    delt => this%delt_buffer%rget(1)

    if (first) then
      euler_scheme = ImplicitEulerSchemeType(this%delt_buffer)
      weight = euler_scheme%get_weight(n)
      return
    else
      delt_prev => this%delt_buffer%rget(2)
      r = delt(1) / delt_prev(1)
      if (n == 1) then
        weight = (1.0_dp + 2.0_dp * r) / (1.0_dp + r)
      elseif (n == 2) then
        weight = -(1.0_dp + r)
      elseif (n == 3) then
        weight = r**2 / (1.0_dp + r)
      else
        call store_error("Weight calculation error", terminate=.TRUE.)
      end if
    end if

    weight = weight / delt(1)

  end function get_weight

end module BDF2SchemeModule
