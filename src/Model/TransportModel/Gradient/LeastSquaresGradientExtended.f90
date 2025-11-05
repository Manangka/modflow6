module LeastSquaresGradientExtendedModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE, DSAME

  Use IGradient
  use BaseDisModule, only: DisBaseType
  use BoundaryFacesModule, only: BoundaryFacesType
  use LinearAlgebraUtilsModule, only: eye
  use PseudoInverseModule, only: pinv
  use DisUtilsModule, only: node_distance

  implicit none
  private

  public :: LeastSquaresGradientExtendedType

  type Array1D
    integer(I4B), dimension(:), allocatable :: data
  end type Array1D

  type Array2D
    real(DP), dimension(:, :), allocatable :: data
  end type Array2D

  type, extends(IGradientType) :: LeastSquaresGradientExtendedType
    private
    class(DisBaseType), pointer :: dis
    real(DP), dimension(:), pointer :: phi
    type(Array2D), allocatable, dimension(:) :: R ! Gradient reconstruction matrix
    type(BoundaryFacesType), allocatable :: boundary_faces

    type(Array1D), allocatable, dimension(:) :: connected_cells
  contains
    procedure :: get
    procedure :: set_field

    procedure, private :: compute_cell_gradient
    procedure, private :: create_gradient_reconstruction_matrix
    procedure, private :: find_connected_cells
    procedure, private :: has_relevant_boundaries
  end type LeastSquaresGradientExtendedType

  interface LeastSquaresGradientExtendedType
    module procedure Constructor
  end interface LeastSquaresGradientExtendedType

contains
  function constructor(dis) Result(gradient)
    ! --dummy
    class(DisBaseType), pointer, intent(in) :: dis
    !-- return
    type(LeastSquaresGradientExtendedType) :: gradient
    ! -- local
    integer(I4B) :: n, nodes

    gradient%dis => dis
    nodes = dis%nodes

    ! -- Create boundary Cells
    gradient%boundary_faces = BoundaryFacesType(dis)

    ! -- Find connected cells (1st and 2nd degree neighbors)
    allocate (gradient%connected_cells(dis%nodes))
    do n = 1, nodes
      gradient%connected_cells(n)%data = gradient%find_connected_cells(n)
    end do

    ! -- Compute the gradient operator
    allocate (gradient%R(dis%nodes))
    do n = 1, nodes
      gradient%R(n)%data = gradient%create_gradient_reconstruction_matrix(n)
    end do

  end function constructor

  function has_relevant_boundaries(this, n) result(has_boundary)
    ! -- dummy
    class(LeastSquaresGradientExtendedType) :: this
    integer(I4B), intent(in) :: n
    logical :: has_boundary
    ! -- local
    integer(I4B) :: ipos, ipos2, m
    integer(I4B) :: number_boundaries
    integer(I4B) :: isympos, ihc
    real(DP), dimension(3) :: boundary_normal, connection_normal
    logical :: in_same_plane

    ! Count boundaries that are in the same plane as connections
    number_boundaries = this%boundary_faces%ia(n + 1) - this%boundary_faces%ia(n)
    do ipos = this%boundary_faces%ia(n), this%boundary_faces%ia(n + 1) - 1
      boundary_normal = this%boundary_faces%get_normal(ipos)

      in_same_plane = .false.
      do ipos2 = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
        m = this%dis%con%ja(ipos2)
        isympos = this%dis%con%jas(ipos2)
        ihc = this%dis%con%ihc(isympos)
        call this%dis%connection_normal( &
          n, m, ihc, &
          connection_normal(1), connection_normal(2), connection_normal(3), &
          ipos2)

        if (abs(dot_product(boundary_normal, connection_normal)) > DSAME) then
          in_same_plane = .true.
          exit
        end if
      end do

      if (.not. in_same_plane) then
        number_boundaries = number_boundaries - 1
      end if
    end do

    has_boundary = number_boundaries > 0

  end function has_relevant_boundaries

  function find_connected_cells(this, n) result(res)
    use PtrHashTableModule, only: PtrHashTableType
    use IteratorModule, only: IteratorType
    ! -- dummy
    class(LeastSquaresGradientExtendedType) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), allocatable :: res(:)
    ! -- local
    integer(I4B) :: ipos, ipos2, m, mm, local_pos
    character(1) :: c
    class(*), pointer :: m_ptr, mm_ptr
    integer(I4B), pointer :: i_ptr
    type(PtrHashTableType) :: connected_cells
    class(IteratorType), allocatable :: itr
    logical :: has_boundary

    has_boundary = this%has_relevant_boundaries(n)

    ! Find 1st and 2nd degree connected cells (Direct neighbours and neighbours of neighbours)
    do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
      m = this%dis%con%ja(ipos)

      m_ptr => this%dis%con%ja(ipos)
      c = char(this%dis%con%ja(ipos))
      if (.not. connected_cells%contains(c)) then
        call connected_cells%add(c, m_ptr)
      end if

      if (.not. has_boundary) cycle

      do ipos2 = this%dis%con%ia(m) + 1, this%dis%con%ia(m + 1) - 1
        mm = this%dis%con%ja(ipos2)

        mm_ptr => this%dis%con%ja(ipos2)
        c = char(this%dis%con%ja(ipos2))
        if (.not. connected_cells%contains(c) .and. mm /= n) then
          call connected_cells%add(c, mm_ptr)
        end if
      end do
    end do

    ! Store the results in an array
    allocate (res(connected_cells%count()))
    local_pos = 1
    allocate (itr, source=connected_cells%iterator())
    do while (itr%has_next())
      call itr%next()

      select type (val => itr%value())
      type is (integer(I4B))
        i_ptr => val
      end select
      res(local_pos) = i_ptr
      local_pos = local_pos + 1
    end do

  end function find_connected_cells

  function create_gradient_reconstruction_matrix(this, n) result(R)
    ! -- dummy
    class(LeastSquaresGradientExtendedType) :: this
    integer(I4B), intent(in) :: n
    real(DP), dimension(:, :), allocatable :: R
    ! -- local
    integer(I4B) :: ipos, m
    real(DP) :: length
    real(DP), dimension(3) :: dnm
    real(DP), dimension(:, :), allocatable :: d
    real(DP), dimension(:, :), allocatable :: d_trans
    real(DP), dimension(:, :), allocatable :: grad_scale
    real(DP), dimension(:, :), allocatable :: W
    real(DP), dimension(3, 3) :: g
    real(DP), dimension(3, 3) :: g_inv
    integer(I4B), allocatable :: connected_cells(:)
    integer(I4B) :: num_connected_cells

    connected_cells = this%connected_cells(n)%data
    num_connected_cells = size(connected_cells)

    allocate (d(num_connected_cells, 3))
    allocate (d_trans(3, num_connected_cells))
    allocate (R(3, num_connected_cells))
    allocate (grad_scale(num_connected_cells, num_connected_cells))
    allocate (W(num_connected_cells, num_connected_cells))

    grad_scale = 0
    d = 0
    W = eye(num_connected_cells)

    ! Assemble the distance matrix
    do ipos = 1, num_connected_cells
      m = connected_cells(ipos)

      dnm = node_distance(this%dis, n, m)
      length = norm2(dnm)

      d(ipos, :) = dnm / length
      grad_scale(ipos, ipos) = 1.0_dp / length
    end do

    d_trans = transpose(d)

    ! Compute the G and inverse G matrices
    g = matmul(d_trans, matmul(W, d))
    g_inv = pinv(g)

    ! Compute the gradient operator
    R = matmul(matmul(matmul(g_inv, d_trans), W), grad_scale)

  end function create_gradient_reconstruction_matrix

  function get(this, n) result(grad_c)
    ! -- dummy
    class(LeastSquaresGradientExtendedType), target :: this
    integer(I4B), intent(in) :: n
    !-- return
    real(DP), dimension(3) :: grad_c

    grad_c = this%compute_cell_gradient(n)
  end function get

  subroutine set_field(this, phi)
    ! -- dummy
    class(LeastSquaresGradientExtendedType), target :: this
    real(DP), dimension(:), pointer, intent(in) :: phi

    this%phi => phi
  end subroutine set_field

  function compute_cell_gradient(this, n) result(grad_c)
    ! -- dummy
    class(LeastSquaresGradientExtendedType), target :: this
    integer(I4B), intent(in) :: n
    !-- return
    real(DP), dimension(3) :: grad_c
    ! -- local
    real(DP), dimension(:, :), pointer :: R
    integer(I4B) :: ipos

    integer(I4B) :: m
    real(DP), dimension(:), allocatable :: dc

    integer(I4B), allocatable :: connected_cells(:)
    integer(I4B) :: num_connected_cells

    connected_cells = this%connected_cells(n)%data
    num_connected_cells = size(connected_cells)

    ! Assemble the concentration difference vector
    allocate (dc(num_connected_cells))
    do ipos = 1, num_connected_cells
      m = connected_cells(ipos)
      dc(ipos) = this%phi(m) - this%phi(n)
    end do

    ! Compute the cells gradient
    R => this%R(n)%data
    grad_c = matmul(R, dc)

  end function compute_cell_gradient

end module
