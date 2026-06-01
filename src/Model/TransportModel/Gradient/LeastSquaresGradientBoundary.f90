module LeastSquaresGradientBoundaryModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE

  Use IGradient
  use BaseDisModule, only: DisBaseType
  use BoundaryFacesModule, only: BoundaryFacesType
  use PseudoInverseModule, only: pinv
  use DisUtilsModule, only: number_connected_faces, number_faces, node_distance

  implicit none
  private

  public :: LeastSquaresGradientBoundaryType

  type Array2D
    real(DP), dimension(:, :), allocatable :: data
  end type Array2D

  type, extends(IGradientType) :: LeastSquaresGradientBoundaryType
    class(DisBaseType), pointer :: dis
    real(DP), dimension(:), pointer :: phi
    type(Array2D), allocatable, dimension(:) :: R
    type(BoundaryFacesType), allocatable :: boundary_faces
  contains
    procedure :: get
    procedure :: set_field

    procedure, private :: compute_cell_gradient
    procedure, private :: create_gradient_reconstruction_matrix
  end type LeastSquaresGradientBoundaryType

  interface LeastSquaresGradientBoundaryType
    module procedure Constructor
  end interface LeastSquaresGradientBoundaryType

contains
  function constructor(dis) Result(gradient)
    ! --dummy
    class(DisBaseType), pointer, intent(in) :: dis
    !-- return
    type(LeastSquaresGradientBoundaryType) :: gradient
    ! -- local
    integer(I4B) :: n, nodes

    gradient%dis => dis
    nodes = dis%nodes

    ! -- Create boundary Cells
    gradient%boundary_faces = BoundaryFacesType(dis)

    ! -- Compute the gradient operator
    nodes = dis%nodes
    allocate (gradient%R(dis%nodes))
    do n = 1, nodes
      gradient%R(n)%data = gradient%create_gradient_reconstruction_matrix(n)
    end do
  end function constructor

  function create_gradient_reconstruction_matrix(this, n) result(R)
    ! -- dummy
    class(LeastSquaresGradientBoundaryType) :: this
    integer(I4B), intent(in) :: n
    real(DP), dimension(:, :), allocatable :: R
    ! -- local
    integer(I4B) :: number_connections
    integer(I4B) :: ipos, local_pos, m
    real(DP) :: length
    real(DP), dimension(3) :: dnm
    real(DP), dimension(:, :), allocatable :: d
    real(DP), dimension(:, :), allocatable :: d_trans
    real(DP), dimension(:, :), allocatable :: inverse_distance
    integer(I4B) :: number_boundaries, number_sides

    number_connections = number_connected_faces(this%dis, n)
    number_boundaries = this%boundary_faces%ia(n + 1) - this%boundary_faces%ia(n)
    number_sides = number_connections + number_boundaries

    allocate (d(number_sides, 3))
    allocate (d_trans(3, number_sides))
    allocate (R(3, number_sides))
    allocate (inverse_distance(number_sides, number_connections))

    inverse_distance = 0
    d = 0

    ! Assemble the distance matrix
    ! Handle the internal connections
    local_pos = 1
    do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
      m = this%dis%con%ja(ipos)

      dnm = node_distance(this%dis, n, m)
      length = norm2(dnm)

      d(local_pos, :) = dnm / length
      inverse_distance(local_pos, local_pos) = 1.0_dp / length

      local_pos = local_pos + 1
    end do

    ! Handle the boundary cells
    do ipos = this%boundary_faces%ia(n), this%boundary_faces%ia(n + 1) - 1
      d(local_pos, :) = this%boundary_faces%get_normal(ipos)
      local_pos = local_pos + 1
    end do

    d_trans = transpose(d)

    ! Compute the gradient reconstructions matrix
    R = matmul(pinv(d), inverse_distance)

  end function create_gradient_reconstruction_matrix

  function get(this, n) result(grad_c)
    ! -- dummy
    class(LeastSquaresGradientBoundaryType), target :: this
    integer(I4B), intent(in) :: n
    !-- return
    real(DP), dimension(3) :: grad_c

    grad_c = this%compute_cell_gradient(n)
  end function get

  subroutine set_field(this, phi)
    ! -- dummy
    class(LeastSquaresGradientBoundaryType), target :: this
    real(DP), dimension(:), pointer, intent(in) :: phi

    this%phi => phi
  end subroutine set_field

  function compute_cell_gradient(this, n) result(grad_c)
    ! -- return
    real(DP), dimension(3) :: grad_c
    ! -- dummy
    class(LeastSquaresGradientBoundaryType), target :: this
    integer(I4B), intent(in) :: n
    ! -- local
    real(DP), dimension(:, :), pointer :: R
    integer(I4B) :: ipos, local_pos
    integer(I4B) :: number_connections

    integer(I4B) :: m
    real(DP), dimension(:), allocatable :: dc

    ! Assemble the concentration difference vector
    number_connections = number_connected_faces(this%dis, n)
    allocate (dc(number_connections))
    local_pos = 1
    do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
      m = this%dis%con%ja(ipos)
      dc(local_pos) = this%phi(m) - this%phi(n)
      local_pos = local_pos + 1
    end do

    ! Compute the cells gradient
    R => this%R(n)%data
    grad_c = matmul(R, dc)

  end function compute_cell_gradient

end module LeastSquaresGradientBoundaryModule
