module ImplicitEulerSchemeModule
  use TimeSchemeInterfaceModule, only: TimeSchemeInterface
  use KindModule, only: DP, I4B, LGP
  use CircularBufferModule, only: CircularBufferType

  implicit none
  private

  public :: ImplicitEulerSchemeType

  type, extends(TimeSchemeInterface) :: ImplicitEulerSchemeType
    private
    class(CircularBufferType), pointer :: delt_buffer => null()
    integer(I4B) :: num_steps = 1 !< number of sub-steps in the time step
  contains
    procedure :: update_delt
    procedure :: get_num_steps
    procedure :: get_time_iteration_steps
    procedure :: get_weight
    final :: destructor
  end type ImplicitEulerSchemeType

  interface ImplicitEulerSchemeType
    module procedure constructor
  end interface ImplicitEulerSchemeType

contains

  function constructor() Result(scheme)
    type(ImplicitEulerSchemeType) :: scheme
    ! -- dummy
    ! -- local
    allocate( scheme%delt_buffer, source=CircularBufferType(scheme%num_steps, 1, 'DELT_BUFFER', 'TDIS'))

  end function constructor

  subroutine destructor(this)
    ! -- dummy
    type(ImplicitEulerSchemeType), intent(inout) :: this
    deallocate (this%delt_buffer)

  end subroutine destructor

  subroutine update_delt(this, new_delt)
    class(ImplicitEulerSchemeType), intent(inout) :: this
    real(DP), intent(in) :: new_delt

    call this%delt_buffer%add([new_delt])
  end subroutine update_delt

  function get_num_steps(this) result(num_steps)
    class(ImplicitEulerSchemeType), intent(in) :: this
    integer(I4B) :: num_steps

    num_steps = this%num_steps
  end function get_num_steps

  function get_time_iteration_steps(this) result(steps)
    ! -- dummy
    class(ImplicitEulerSchemeType), intent(in) :: this
    integer(I4B) :: steps
    ! -- local
    steps = 1
    return

  end function get_time_iteration_steps

  function get_weight(this, n) result(weight)
    use SimModule, only: store_error
    ! -- dummy
    class(ImplicitEulerSchemeType), intent(in) :: this
    real(DP) :: weight
    integer(I4B), intent(in) :: n
    ! -- local
    real(DP), pointer :: delt(:)

    if (n == 1) then
      weight = 1.0_dp
    elseif (n == 2) then
      weight = -1.0_dp
    else
      call store_error("Weight calculation error", terminate=.TRUE.)
    end if

    delt => this%delt_buffer%rget(1)
    weight = weight / delt(1)

  end function get_weight

end module ImplicitEulerSchemeModule
