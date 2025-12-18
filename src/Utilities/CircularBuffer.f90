module CircularBufferModule

  use KindModule, only: DP, I4B
  use MemoryManagerModule, only: mem_allocate, mem_deallocate

  implicit none
  private
  public :: CircularBufferType

  type :: CircularBufferType
    integer :: capacity !< Maximum number of elements in the buffer
    integer :: size = 0 !< Current number of elements in the buffer
    integer :: head_index = 0 !< Index of the head element
    real(DP), pointer, contiguous :: data(:, :) !< The buffer array

  contains
    procedure :: add
    procedure :: get
    procedure :: rget
    ! procedure :: is_full
    ! procedure :: is_empty
    final :: destructor
  end type CircularBufferType

  interface CircularBufferType
    module procedure constructor
  end interface CircularBufferType

contains
  function constructor(capacity, neq, name, mem_path) Result(buffer)
    type(CircularBufferType) :: buffer
    ! -- dummy
    integer(I4B), intent(in) :: capacity ! -- maximum number of elements
    integer(I4B), intent(in) :: neq ! -- number of equations (size of each element)
    character(len=*), intent(in) :: name !< variable name
    character(len=*), intent(in) :: mem_path !< path where variable is stored
    ! -- local

    call mem_allocate(buffer%data, neq, capacity, name, mem_path)
    buffer%capacity = capacity

  end function constructor

  subroutine destructor(this)
    ! -- dummy
    type(CircularBufferType), intent(inout) :: this

    call mem_deallocate(this%data)
  end subroutine destructor

  subroutine add(this, element)
    ! -- dummy
    class(CircularBufferType), target :: this
    real(DP), dimension(:), intent(in) :: element
    ! -- local
    integer :: insert_index

    if (this%size < this%capacity) then
      this%size = this%size + 1
    end if

    insert_index = mod(this%head_index, this%capacity) + 1
    this%data(:, insert_index) = element
    this%head_index = insert_index

  end subroutine add

  function get(this, n) result(element)
    ! -- dummy
    class(CircularBufferType), target :: this
    integer(I4B), intent(in) :: n !< 1-based index of the element to retrieve (1 = oldest)
    ! -- return
    real(DP), pointer, dimension(:) :: element
    ! -- local
    integer :: index

    if (n < 1 .or. n > this%size) then
      error stop "Index out of bounds in CircularBufferType%get"
    end if

    index = mod(this%head_index - this%size + n - 1 + this%capacity, &
                this%capacity) + 1
    element => this%data(:, index)

  end function get

  function rget(this, n) result(element)
    ! -- dummy
    class(CircularBufferType), target :: this
    integer(I4B), intent(in) :: n !< 1-based index of the element to retrieve (1 = newest)
    ! -- return
    real(DP), pointer, dimension(:) :: element
    ! -- local
    integer :: index
    if (n < 1 .or. n > this%size) then
      error stop "Index out of bounds in CircularBufferType%rget"
    end if

    index = mod(this%head_index - n + this%capacity, this%capacity) + 1
    element => this%data(:, index)

  end function rget

end module CircularBufferModule
