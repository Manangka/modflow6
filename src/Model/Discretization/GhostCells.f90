module GhostCellsModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE, DZERO, DSAME

  use BaseDisModule, only: DisBaseType
  use DisuModule, only: DisuType

  implicit none
  private
  public :: GhostCellsType
  public :: get_number_sides

  type GhostCellsType
    real(DP), dimension(:), allocatable :: xc
    real(DP), dimension(:), allocatable :: yc
    real(DP), dimension(:), allocatable :: zc
    integer(I4B), dimension(:), allocatable :: ia
  contains
    procedure, private :: create_ghost_cells
  end type GhostCellsType

  interface GhostCellsType
      module procedure constructor
  end interface GhostCellsType

contains

  function constructor(dis) result(ghost_cells)
  ! -- return
  type(GhostCellsType) :: ghost_cells
  ! -- dummy
  class(DisBaseType), intent(in) :: dis

  select type(dis)
  class is (DisuType)
    print*, "DisU not implemented"
  class default
    call create_ghost_cells(ghost_cells, dis)
  end select

  end function constructor


  subroutine create_ghost_cells(this, dis)
    ! -- dummy
    class(GhostCellsType) :: this
    class(DisBaseType), intent(in) :: dis
    ! -- local
    integer(I4B) :: nodes
    integer(I4B) :: n, ipos
    integer(I4B) :: number_ghostcells, number_local_ghostcells
    real(DP) , dimension(:, :, :), allocatable :: ghost_boundary_edges
    integer(I4B) :: ghost_cell_id, size_ghost_boundary_edges
    real(DP), dimension(:), allocatable :: v1, v2, ghost_centroid
    real(DP) :: node_xc, node_yc, node_zc

    nodes = dis%nodes
    number_ghostcells = 0

    ! Number of ghostcells
    do n = 1, nodes
      number_local_ghostcells = ghostcell_count(dis, n)
      number_ghostcells = number_ghostcells + number_local_ghostcells
    end do

    ! Allocate memory
    allocate(this%xc(number_ghostcells))
    allocate(this%yc(number_ghostcells))
    allocate(this%zc(number_ghostcells))
    allocate(this%ia(nodes + 1))

    ! Create ghostcells
    this%ia(1) = 1
    ghost_cell_id = 1
    do n = 1, nodes
      ! Increment ia for the new ghostcells
      number_local_ghostcells = ghostcell_count(dis, n)
      this%ia(n + 1) = this%ia(n) + number_local_ghostcells

      ! Get the edges belonging to the ghost cells
      ghost_boundary_edges = find_ghost_cell_boundary_edge(dis, n)

      ! Reflect the centroid of the node over the ghost boundary edges
      size_ghost_boundary_edges = size(ghost_boundary_edges, dim=1)
      do ipos = 1, size_ghost_boundary_edges
        v1 = ghost_boundary_edges(ipos, 1, :)
        v2 = ghost_boundary_edges(ipos, 2, :)

        node_xc = dis%xc(n)
        node_yc = dis%yc(n)
        node_zc = (dis%top(n) + dis%bot(n)) / 2.0_dp

        ghost_centroid = reflect_point(v1, v2, [node_xc, node_yc, node_zc])

        this%xc(ghost_cell_id) = ghost_centroid(1)
        this%yc(ghost_cell_id) = ghost_centroid(2)
        this%zc(ghost_cell_id) = ghost_centroid(3)

        ghost_cell_id = ghost_cell_id + 1
      end do
    end do
  end subroutine create_ghost_cells

  function get_number_sides(dis, n) result(number_sides)
    ! -- dummy
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: n
    integer(I4B) :: number_sides
    ! -- local
    real(DP), dimension(:,:), allocatable :: polyverts

    call dis%get_polyverts(n, polyverts)
    number_sides = size(polyverts, 2) + 2 ! We add 2 for the top and bottom
  end function get_number_sides


  function ghostcell_count(dis, n) result(number_local_ghostcells)
    ! -- dummy
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: n
    integer(I4B) :: number_local_ghostcells
    ! -- local
    integer(I4B) :: number_sides, number_connections
    real(DP), dimension(:,:), allocatable :: polyverts

    number_sides = get_number_sides(dis, n)
    number_connections = number_connected_nodes(dis, n)

    number_local_ghostcells = number_sides - number_connections

  end function ghostcell_count


  function number_connected_nodes(dis, n) result(number_connections)
    ! -- dummy
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: n
    integer(I4B) :: number_connections

    number_connections = dis%con%ia(n + 1) - dis%con%ia(n) - 1
  end function


  function  find_ghost_cell_boundary_edge(dis, n) result(boundary_edges)
    ! -- dummy
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: n
    real(DP), allocatable :: boundary_edges(:,:,:)
    ! -- local
    real(DP), dimension(:,:), allocatable :: polyverts
    integer(I4B) :: m, num_sides, iedge, ipos
    logical, dimension(:), allocatable :: is_ghost_boundary
    real(DP) :: x_dir, y_dir, z_dir, prod
    real(DP), dimension(3) :: edge_dir, cross
    integer(I4B) :: num_ghostcells
    integer(I4B) :: isympos, ihc
    real(DP) :: node_z
   

    call dis%get_polyverts(n, polyverts)
    num_sides = size(polyverts, 2)

    ! We start of by flagging all boundaries as ghost cells.
    ! Later we will check if they are connected to a real cell
    allocate(is_ghost_boundary(num_sides + 2))
    is_ghost_boundary = .true.
    node_z = (dis%top(n) + dis%bot(n)) / 2.0_dp
  
    ! Check to which edge the connection belongs.
    ! When we find the edge we unflag it as a ghost cell
    do ipos = dis%con%ia(n) + 1, dis%con%ia(n + 1) - 1
      m = dis%con%ja(ipos)
      call dis%connection_normal(n, m, 1,x_dir, y_dir, z_dir, ipos)

      ! Determine connection direction. 
      ! Horizontal and vertical connections are treated differently
      isympos = dis%con%jas(ipos)
      ihc = dis%con%ihc(isympos)

      if (ihc == 1) then
        ! Check if the horizontal connection is a ghost cell
        ! We first find the edge to which the connection belongs
        ! Then we take the cross product of the edge direction and the z unit vector
        ! This will give a vector that is perpendicular to the edge
        ! If the dot product of this vector and the connection direction is 1
        ! the connection is not a ghost cell
        do iedge = 1, num_sides
          edge_dir(1:2) = polyverts(:, mod(iedge, num_sides) + 1) - polyverts(:,iedge)
          edge_dir(3) = DZERO
          edge_dir = edge_dir / norm2(edge_dir)
          cross = cross_product(edge_dir, [0.0_dp, 0.0_dp, -1.0_dp])

          prod = dot_product(cross, [x_dir, y_dir, z_dir])
          if (abs(prod - DONE) < DSAME) then
            is_ghost_boundary(iedge) = .false.
            exit
          end if
        end do
      else
        ! Check if the vertical connection is a ghost cell
        ! By definition the last two edges are the bottom and top
        if (abs(dis%bot(n) - dis%top(m)) < DSAME) then
          is_ghost_boundary(num_sides + 1) = .false.
        end if

        if (abs(dis%top(n) - dis%bot(m)) < DSAME) then
          is_ghost_boundary(num_sides + 2) = .false.
        end if

      end if
    end do

    num_ghostcells = count(is_ghost_boundary)
    allocate(boundary_edges(num_ghostcells, 2, 3))

    ! Collect all the ghost cell edges
    ipos = 1
    boundary_edges = 0
    do iedge = 1, num_sides
      if (is_ghost_boundary(iedge)) then
        boundary_edges(ipos, 1, 1:2) = polyverts(:, iedge)
        boundary_edges(ipos, 2, 1:2) = polyverts(:, mod(iedge, num_sides) + 1)
        boundary_edges(ipos, 1:2, 3) = node_z
        ipos = ipos + 1
      end if
    end do

    ! The top and the bottom are not really edges but we treat them as such
    ! We only need a line on the top/bottom to reflect the centroid over
    if (is_ghost_boundary(num_sides + 1 )) then 
      boundary_edges(ipos, 1, :) = [dis%xc(n),  dis%yc(n), dis%bot(n)]
      boundary_edges(ipos, 2, :) = [dis%xc(n) + 1.0_dp,  dis%yc(n), dis%bot(n)]
      ipos = ipos + 1
    end if

    if (is_ghost_boundary(num_sides +2  )) then 
      boundary_edges(ipos, 1, :) = [dis%xc(n), dis%yc(n), dis%top(n)]
      boundary_edges(ipos, 2, :) = [dis%xc(n) + 1.0_dp, dis%yc(n), dis%top(n)]
      ipos = ipos + 1
    end if

  end function find_ghost_cell_boundary_edge
  

  function cross_product(a, b) result(c)
    ! -- return
    real(DP), dimension(3) :: c
    ! -- dummy
    real(DP), dimension(3), intent(in) :: a
    real(DP), dimension(3), intent(in) :: b

    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)
  end function cross_product


  function reflect_point(v1, v2, point) result(reflected_point)
    ! -- return
    real(DP), dimension(3) :: reflected_point
    ! -- dummy
    real(DP), dimension(3), intent(in) :: v1
    real(DP), dimension(3), intent(in) :: v2
    real(DP), dimension(3), intent(in) :: point
    ! -- local
    real(DP) , dimension(3) :: v12, p, b
    real(DP) :: proj, costheta

    v12 = v2 - v1
    p = point - v1
    proj = dot_product(p, v12)
    costheta = norm2(v12) * norm2(p)
    b = p - proj * v12 / norm2(v12)
    reflected_point = point - 2.0_dp * b

  end function reflect_point

end module GhostCellsModule