module IMEXMatrixDecoratorModule
  use MatrixBaseModule
  use KindModule, only: I4B, DP
  use VectorBaseModule
  use SparseModule, only: sparsematrix
  use SimModule, only: store_error

  implicit none
  private
  
  type, public, extends(MatrixBaseType) :: IMEXMatrixType
    class(MatrixBaseType), pointer :: matrix => null()
    integer(I4B), dimension(:), pointer :: ja => null() !< column indices array from CSR format
    integer(I4B), dimension(:), pointer :: ia => null() !< row pointers
    real(DP), dimension(:), pointer, contiguous :: x => null()
    real(DP), dimension(:), pointer, contiguous :: rhs => null() !right-hand side vector

    real(DP) :: coef
    logical :: move_to_rhs = .true.
  contains
    procedure :: init
    procedure :: destroy
    procedure :: create_vec_mm
    procedure :: create_vec
    procedure :: get_value_pos
    procedure :: get_diag_value
    procedure :: set_diag_value
    procedure :: set_value_pos
    procedure :: add_value_pos
    procedure :: add_diag_value
    procedure :: zero_entries
    procedure :: zero_row_offdiag
    procedure :: get_first_col_pos
    procedure :: get_last_col_pos
    procedure :: get_column
    procedure :: get_position
    procedure :: get_position_diag
    procedure :: get_aij
    procedure :: get_row_offset
    procedure :: multiply
  end type

contains

  subroutine init(this, sparse, mem_path)
    class(IMEXMatrixType) :: this
    type(sparsematrix) :: sparse
    character(len=*) :: mem_path
    call store_error('Program error: init not implemented.', terminate=.true.)
  end subroutine

  subroutine destroy(this)
    class(IMEXMatrixType) :: this
    call store_error('Program error: destroy not implemented.', terminate=.true.)
  end subroutine

  function create_vec_mm(this, n, name, mem_path) result(vec)
    class(IMEXMatrixType) :: this
    integer(I4B) :: n
    character(len=*) :: name
    character(len=*) :: mem_path
    class(VectorBaseType), pointer :: vec
    call store_error('Program error: create_vec_mm not implemented.', terminate=.true.)
    vec => null()
  end function

  function create_vec(this, n) result(vec)
    class(IMEXMatrixType) :: this
    integer(I4B) :: n
    class(VectorBaseType), pointer :: vec
    call store_error('Program error: create_vec not implemented.', terminate=.true.)
    vec => null()
  end function

  function get_value_pos(this, ipos) result(value)
    class(IMEXMatrixType) :: this
    integer(I4B) :: ipos
    real(DP) :: value
    call store_error('Program error: get_value_pos not implemented.', terminate=.true.)
    value = 0.0_DP
  end function

  function get_diag_value(this, irow) result(diag_value)
    class(IMEXMatrixType) :: this
    integer(I4B) :: irow
    real(DP) :: diag_value
    call store_error('Program error: get_diag_value not implemented.', terminate=.true.)
    diag_value = 0.0_DP
  end function

  subroutine set_diag_value(this, irow, diag_value)
    class(IMEXMatrixType) :: this
    integer(I4B) :: irow
    real(DP) :: diag_value
    call store_error('Program error: set_diag_value not implemented.', terminate=.true.)
  end subroutine

  subroutine set_value_pos(this, ipos, value)
    class(IMEXMatrixType) :: this
    integer(I4B) :: ipos
    real(DP) :: value
    call store_error('Program error: set_value_pos not implemented.', terminate=.true.)
  end subroutine

  subroutine add_value_pos(this, ipos, value)
    class(IMEXMatrixType) :: this
    integer(I4B) :: ipos
    real(DP) :: value
    integer(I4B) :: icol, irow, i
    integer(I4B), dimension(:), pointer, contiguous :: ia
    integer(I4B), dimension(:), pointer, contiguous :: ja
    real(DP), dimension(:), pointer, contiguous :: amat

    call this%matrix%get_aij(ia, ja, amat)
    icol = ja(ipos)

    irow = 0
    do i = 1, size(ia) - 1
        if (ipos >= ia(i) .and. ipos < ia(i+1)) then
            irow = i
            exit
        end if
    end do

    if (this%move_to_rhs) then
        this%rhs(irow) = this%rhs(irow) - value * this%coef * this%x(icol)
    else
      call this%matrix%add_value_pos(ipos, value * this%coef)
    end if

  end subroutine

  subroutine add_diag_value(this, irow, value)
    class(IMEXMatrixType) :: this
    integer(I4B) :: irow
    real(DP) :: value
    call store_error('Program error: add_diag_value not implemented.', terminate=.true.)
  end subroutine

  subroutine zero_entries(this)
    class(IMEXMatrixType) :: this
    call store_error('Program error: zero_entries not implemented.', terminate=.true.)
  end subroutine

  subroutine zero_row_offdiag(this, irow)
    class(IMEXMatrixType) :: this
    integer(I4B) :: irow
    call store_error('Program error: zero_row_offdiag not implemented.', terminate=.true.)
  end subroutine

  function get_first_col_pos(this, irow) result(first_col_pos)
    class(IMEXMatrixType) :: this
    integer(I4B) :: irow
    integer(I4B) :: first_col_pos
    call store_error('Program error: get_first_col_pos not implemented.', terminate=.true.)
    first_col_pos = 0
  end function

  function get_last_col_pos(this, irow) result(last_col_pos)
    class(IMEXMatrixType) :: this
    integer(I4B) :: irow
    integer(I4B) :: last_col_pos
    call store_error('Program error: get_last_col_pos not implemented.', terminate=.true.)
    last_col_pos = 0
  end function

  function get_column(this, ipos) result(icol)
    class(IMEXMatrixType) :: this
    integer(I4B) :: ipos
    integer(I4B) :: icol
    call store_error('Program error: get_column not implemented.', terminate=.true.)
    icol = 0
  end function

  function get_position(this, irow, icol) result(ipos)
    class(IMEXMatrixType) :: this
    integer(I4B) :: irow
    integer(I4B) :: icol
    integer(I4B) :: ipos
    call store_error('Program error: get_position not implemented.', terminate=.true.)
    ipos = 0
  end function

  function get_position_diag(this, irow) result(ipos_diag)
    class(IMEXMatrixType) :: this
    integer(I4B) :: irow
    integer(I4B) :: ipos_diag
    call store_error('Program error: get_position_diag not implemented.', terminate=.true.)
    ipos_diag = 0
  end function

  subroutine get_aij(this, ia, ja, amat)
    class(IMEXMatrixType) :: this
    integer(I4B), dimension(:), pointer, contiguous :: ia
    integer(I4B), dimension(:), pointer, contiguous :: ja
    real(DP), dimension(:), pointer, contiguous :: amat
    call store_error('Program error: get_aij not implemented.', terminate=.true.)
  end subroutine

  function get_row_offset(this) result(offset)
    class(IMEXMatrixType) :: this
    integer(I4B) :: offset
    call store_error('Program error: get_row_offset not implemented.', terminate=.true.)
    offset = 0
  end function

  subroutine multiply(this, vec_x, vec_y)
    class(IMEXMatrixType) :: this
    class(VectorBaseType), pointer :: vec_x
    class(VectorBaseType), pointer :: vec_y
    call store_error('Program error: multiply not implemented.', terminate=.true.)
  end subroutine

end module IMEXMatrixDecoratorModule