module TimeSchemeInterfaceModule
  use KindModule, only: I4B, DP

  implicit none
  private

  public :: TimeSchemeInterface

  type, abstract :: TimeSchemeInterface
  contains
    procedure(pop_delt), deferred :: pop_delt
    procedure(update_delt), deferred :: update_delt
    procedure(get_num_steps), deferred :: get_num_steps
    procedure(get_time_iteration_steps), deferred :: get_time_iteration_steps
    procedure(get_weight), deferred :: get_weight

  end type TimeSchemeInterface

  abstract interface
    subroutine pop_delt(this)
      import :: TimeSchemeInterface
      class(TimeSchemeInterface), intent(inout) :: this
    end subroutine pop_delt
  end interface

  abstract interface
    subroutine update_delt(this, new_delt)
      import :: TimeSchemeInterface
      import :: DP
      class(TimeSchemeInterface), intent(inout) :: this
      real(DP), intent(in) :: new_delt
    end subroutine update_delt
  end interface

  abstract interface
    function get_num_steps(this) result(num_steps)
      import :: TimeSchemeInterface
      import :: I4B
      class(TimeSchemeInterface), intent(in) :: this
      integer(I4B) :: num_steps
    end function get_num_steps
  end interface

  abstract interface
    function get_time_iteration_steps(this) result(steps)
      import :: TimeSchemeInterface
      import :: I4B
      class(TimeSchemeInterface), intent(in) :: this
      integer(I4B) :: steps
    end function get_time_iteration_steps
  end interface

  abstract interface
    function get_weight(this, n) result(weight)
      import :: TimeSchemeInterface
      import :: I4B, DP
      class(TimeSchemeInterface), intent(in) :: this
      integer(I4B), intent(in) :: n
      real(DP) :: weight
    end function get_weight
  end interface

end module TimeSchemeInterfaceModule
