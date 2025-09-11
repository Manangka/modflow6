module IMEXVectorDecoratorModule
  use VectorBaseModule
  use KindModule, only: I4B, DP, LGP
  use SimModule, only: store_error

  implicit none
  private
  
  type, public, extends(VectorBaseType) :: IMEXVectorType
    class(VectorBaseType), pointer :: vector => null()
    real(DP) :: coef
  contains
    procedure :: create_mm
    procedure :: destroy
    procedure :: get_array
    procedure :: get_ownership_range
    procedure :: get_size
    procedure :: get_value_local
    procedure :: zero_entries
    procedure :: set_value_local
    procedure :: add_value_local
    procedure :: axpy
    procedure :: norm2
    procedure :: print
  end type

contains

  subroutine create_mm(this, n, name, mem_path)
    class(IMEXVectorType) :: this
    integer(I4B) :: n
    character(len=*) :: name
    character(len=*) :: mem_path
    call store_error('Program error: create_mm not implemented.', terminate=.true.)
  end subroutine

  subroutine destroy(this)
    class(IMEXVectorType) :: this
    call store_error('Program error: destroy not implemented.', terminate=.true.)
  end subroutine

  function get_array(this) result(array)
    class(IMEXVectorType) :: this
    real(DP), dimension(:), pointer, contiguous :: array
    call store_error('Program error: get_array not implemented.', terminate=.true.)
    array => null()
  end function

  subroutine get_ownership_range(this, start, end)
    class(IMEXVectorType) :: this
    integer(I4B) :: start, end
    call store_error('Program error: get_ownership_range not implemented.', terminate=.true.)
    start = 0
    end = 0
  end subroutine

  function get_size(this) result(size)
    class(IMEXVectorType) :: this
    integer(I4B) :: size
    call store_error('Program error: get_size not implemented.', terminate=.true.)
    size = 0
  end function

  function get_value_local(this, idx) result(val)
    class(IMEXVectorType) :: this
    integer(I4B) :: idx
    real(DP) :: val
    call store_error('Program error: get_value_local not implemented.', terminate=.true.)
    val = 0.0_DP
  end function

  subroutine zero_entries(this)
    class(IMEXVectorType) :: this
    call store_error('Program error: zero_entries not implemented.', terminate=.true.)
  end subroutine

  subroutine set_value_local(this, idx, val)
    class(IMEXVectorType) :: this
    integer(I4B) :: idx
    real(DP) :: val
    call store_error('Program error: set_value_local not implemented.', terminate=.true.)
  end subroutine

  subroutine add_value_local(this, idx, val)
    class(IMEXVectorType) :: this
    integer(I4B) :: idx
    real(DP) :: val
    call this%vector%add_value_local(idx, val * this%coef)
  end subroutine

  subroutine axpy(this, alpha, vec_x)
    class(IMEXVectorType) :: this
    real(DP) :: alpha
    class(VectorBaseType), pointer :: vec_x
    call store_error('Program error: axpy not implemented.', terminate=.true.)
  end subroutine

  function norm2(this) result(n2)
    class(IMEXVectorType) :: this
    real(DP) :: n2
    call store_error('Program error: norm2 not implemented.', terminate=.true.)
    n2 = 0.0_DP
  end function

  subroutine print(this)
    class(IMEXVectorType) :: this
    call store_error('Program error: print not implemented.', terminate=.true.)
  end subroutine

end module IMEXVectorDecoratorModule
