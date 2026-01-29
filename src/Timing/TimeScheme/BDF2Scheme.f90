module BDF2SchemeModule
  use TimeSchemeInterfaceModule, only: TimeSchemeInterface
  use KindModule, only: DP, I4B, LGP
  use CircularBufferModule, only: CircularBufferType

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
    procedure :: pop_delt
    procedure :: update_delt
    procedure :: get_num_steps
    procedure :: get_time_iteration_steps
    procedure :: get_weight
    final :: destructor
  end type BDF2SchemeType

  interface BDF2SchemeType
    module procedure constructor
  end interface BDF2SchemeType

contains
  function constructor(kstp, kper) Result(scheme)
    type(BDF2SchemeType) :: scheme
    ! -- dummy
    integer(I4B), pointer, intent(in) :: kstp
    integer(I4B), pointer, intent(in) :: kper

    scheme%kstp => kstp
    scheme%kper => kper
    allocate (scheme%delt_buffer, source= &
              CircularBufferType(scheme%num_steps, 1, 'DELT_BUFFER', 'TDIS'))

  end function constructor

  subroutine destructor(this)
    ! -- dummy
    type(BDF2SchemeType), intent(inout) :: this

    deallocate (this%delt_buffer)

  end subroutine destructor

  subroutine pop_delt(this)
    class(BDF2SchemeType), intent(inout) :: this

    call this%delt_buffer%pop()
  end subroutine pop_delt

  subroutine update_delt(this, new_delt)
    class(BDF2SchemeType), intent(inout) :: this
    real(DP), intent(in) :: new_delt

    call this%delt_buffer%add([new_delt])
  end subroutine update_delt

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

    first = this%kstp == 1 .and. this%kper == 1

    if (first) then
      steps = 1
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

    first = this%kstp == 1 .and. this%kper == 1
    delt => this%delt_buffer%rget(1)

    if (first) then
      if (n == 1) then
        weight = 1.0_dp
      elseif (n == 2) then
        weight = -1.0_dp
      else
        call store_error("Weight calculation error", terminate=.TRUE.)
      end if
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
