module TspAdvModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE, DZERO, DHALF, DTWO, DNODATA, DPREC, &
                             LINELENGTH
  use NumericalPackageModule, only: NumericalPackageType
  use BaseDisModule, only: DisBaseType
  use TspFmiModule, only: TspFmiType
  use TspAdvOptionsModule, only: TspAdvOptionsType
  use MatrixBaseModule
  use ForsytheMalcolmMoler

  implicit none
  private
  public :: TspAdvType
  public :: adv_cr

  type, extends(NumericalPackageType) :: TspAdvType

    integer(I4B), pointer :: iadvwt => null() !< advection scheme (0 up, 1 central, 2 tvd)
    real(DP), pointer :: ats_percel => null() !< user-specified fractional number of cells advection can move a particle during one time step
    integer(I4B), dimension(:), pointer, contiguous :: ibound => null() !< pointer to model ibound
    type(TspFmiType), pointer :: fmi => null() !< pointer to fmi object
    real(DP), pointer :: eqnsclfac => null() !< governing equation scale factor; =1. for solute; =rhow*cpw for energy

  contains

    procedure :: adv_df
    procedure :: adv_ar
    procedure :: adv_dt
    procedure :: adv_fc
    procedure :: adv_cq
    procedure :: adv_da

    procedure :: allocate_scalars
    procedure, private :: read_options
    procedure, private :: advqtvd
    procedure, private :: advtvd_bd
    procedure, private :: advqtvd_experimental
    procedure, private :: compute_cell_gradient
    procedure, private :: compute_cell_gradient_2dorder
    procedure, private :: node_distance
    procedure :: adv_weight
    procedure :: advtvd
    procedure :: limiter

  end type TspAdvType

contains

  !> @ brief Create a new ADV object
  !!
  !!  Create a new ADV package
  !<
  subroutine adv_cr(advobj, name_model, inunit, iout, fmi, eqnsclfac)
    ! -- dummy
    type(TspAdvType), pointer :: advobj
    character(len=*), intent(in) :: name_model
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    type(TspFmiType), intent(in), target :: fmi
    real(DP), intent(in), pointer :: eqnsclfac !< governing equation scale factor
    !
    ! -- Create the object
    allocate (advobj)
    !
    ! -- create name and memory path
    call advobj%set_names(1, name_model, 'ADV', 'ADV')
    !
    ! -- Allocate scalars
    call advobj%allocate_scalars()
    !
    ! -- Set variables
    advobj%inunit = inunit
    advobj%iout = iout
    advobj%fmi => fmi
    advobj%eqnsclfac => eqnsclfac
  end subroutine adv_cr

  !> @brief Define ADV object
  !!
  !! Define the ADV package
  !<
  subroutine adv_df(this, adv_options)
    ! -- dummy
    class(TspAdvType) :: this
    type(TspAdvOptionsType), optional, intent(in) :: adv_options !< the optional options, for when not constructing from file
    ! -- local
    character(len=*), parameter :: fmtadv = &
      "(1x,/1x,'ADV-- ADVECTION PACKAGE, VERSION 1, 8/25/2017', &
      &' INPUT READ FROM UNIT ', i0, //)"
    !
    ! -- Read or set advection options
    if (.not. present(adv_options)) then
      !
      ! -- Initialize block parser (adv has no define, so it's
      ! not done until here)
      call this%parser%Initialize(this%inunit, this%iout)
      !
      ! --print a message identifying the advection package.
      write (this%iout, fmtadv) this%inunit
      !
      ! --read options from file
      call this%read_options()
    else
      !
      ! --set options from input arg
      this%iadvwt = adv_options%iAdvScheme
    end if
  end subroutine adv_df

  !> @brief Allocate and read method for package
  !!
  !!  Method to allocate and read static data for the ADV package.
  !<
  subroutine adv_ar(this, dis, ibound)
    ! -- modules
    ! -- dummy
    class(TspAdvType) :: this
    class(DisBaseType), pointer, intent(in) :: dis
    integer(I4B), dimension(:), pointer, contiguous, intent(in) :: ibound
    ! -- local
    ! -- formats
    !
    ! -- adv pointers to arguments that were passed in
    this%dis => dis
    this%ibound => ibound
  end subroutine adv_ar

  !> @brief  Calculate maximum time step length
  !!
  !!  Return the largest time step that meets stability constraints
  !<
  subroutine adv_dt(this, dtmax, msg, thetam)
    ! dummy
    class(TspAdvType) :: this !< this instance
    real(DP), intent(out) :: dtmax !< maximum allowable dt subject to stability constraint
    character(len=*), intent(inout) :: msg !< package/cell dt constraint message
    real(DP), dimension(:), intent(in) :: thetam !< porosity
    ! local
    integer(I4B) :: n
    integer(I4B) :: m
    integer(I4B) :: ipos
    integer(I4B) :: nrmax
    character(len=LINELENGTH) :: cellstr
    real(DP) :: dt
    real(DP) :: flowmax
    real(DP) :: flowsumpos
    real(DP) :: flowsumneg
    real(DP) :: flownm
    real(DP) :: cell_volume
    dtmax = DNODATA
    nrmax = 0
    msg = ''

    ! If ats_percel not specified by user, then return without making
    ! the courant time step calculation
    if (this%ats_percel == DNODATA) then
      return
    end if

    ! Calculate time step lengths based on stability constraint for each cell
    ! and store the smallest one
    do n = 1, this%dis%nodes
      if (this%ibound(n) == 0) cycle
      flowsumneg = DZERO
      flowsumpos = DZERO
      do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
        if (this%dis%con%mask(ipos) == 0) cycle
        m = this%dis%con%ja(ipos)
        if (this%ibound(m) == 0) cycle
        flownm = this%fmi%gwfflowja(ipos)
        if (flownm < DZERO) then
          flowsumneg = flowsumneg - flownm
        else
          flowsumpos = flowsumpos + flownm
        end if
      end do
      flowmax = max(flowsumneg, flowsumpos)
      if (flowmax < DPREC) cycle
      cell_volume = this%dis%get_cell_volume(n, this%dis%top(n))
      dt = cell_volume * this%fmi%gwfsat(n) * thetam(n) / flowmax
      dt = dt * this%ats_percel
      if (dt < dtmax) then
        dtmax = dt
        nrmax = n
      end if
    end do
    if (nrmax > 0) then
      call this%dis%noder_to_string(nrmax, cellstr)
      write (msg, *) adjustl(trim(this%memoryPath))//'-'//trim(cellstr)
    end if
  end subroutine adv_dt

  !> @brief  Fill coefficient method for ADV package
  !!
  !!  Method to calculate coefficients and fill amat and rhs.
  !<
  subroutine adv_fc(this, nodes, matrix_sln, idxglo, cnew, rhs)
    use TdisModule, only: kstp, kper, delt
    ! -- modules
    ! -- dummy
    class(TspAdvType) :: this
    integer, intent(in) :: nodes
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), intent(in), dimension(:) :: idxglo
    real(DP), intent(in), dimension(:) :: cnew
    real(DP), dimension(:), intent(inout) :: rhs
    real(DP), allocatable, dimension(:) :: cnew2, rhs2, rhs2_old
    ! -- local
    integer(I4B) :: n, m, idiag, ipos
    real(DP) :: omega, qnm, q, dt
    integer(I4B) time_idx
    real(DP), dimension(3):: rk3_K_weights =(/0.0_dp, 0.5_dp, 2.0_dp /)
    real(DP), dimension(3):: rk3_rhs_weights =(/1.0_dp/6.0_dp, 4.0_dp/6.0_dp, 1.0_dp/6.0_dp /)
    real(DP), dimension(3):: rk3_rhs_old_weights =(/0.0_dp, 0.0_dp, -1.0_dp  /)

    dt = delt
    allocate(rhs2(nodes))
    rhs2_old = rhs2
    do time_idx = 1, 3
      cnew2 = cnew + delt * (rk3_K_weights(time_idx) * rhs2 + rk3_rhs_old_weights(time_idx) * rhs2_old)
      rhs2_old = rhs2
      rhs2 = 0
      !
      ! -- Calculate advection terms and add to solution rhs and hcof.  qnm
      !    is the volumetric flow rate and has dimensions of L^/T.
      do n = 1, nodes
        if (this%ibound(n) == 0) cycle
        idiag = this%dis%con%ia(n)
        do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
          if (this%dis%con%mask(ipos) == 0) cycle
          m = this%dis%con%ja(ipos)
          if (m <= n) cycle
          if (this%ibound(m) == 0) cycle
          qnm = this%fmi%gwfflowja(ipos) * this%eqnsclfac
          omega = this%adv_weight(this%iadvwt, ipos, n, m, qnm)
          ! call matrix_sln%add_value_pos(idxglo(ipos), qnm * (DONE - omega))
          ! call matrix_sln%add_value_pos(idxglo(idiag), qnm * omega)
          if (qnm > 0) then
            q = qnm * cnew2(m)
          else
            q = qnm * cnew2(n)
          end if
          rhs2(n) = rhs2(n) - q
          rhs2(m) = rhs2(m) + q

        end do
      end do
      !
      ! -- TVD
      if (this%iadvwt >= 2) then
        do n = 1, nodes
          if (this%ibound(n) == 0) cycle
          call this%advtvd(n, cnew2, rhs2)
        end do
      end if
      
      ! ! Heuns method
      ! cnew2 = cnew + delt * rhs2
      ! rhs = rhs + 0.5_dp * rhs2

      ! Runge-Kutta 3th
      rhs = rhs + rk3_rhs_weights(time_idx) * rhs2
    end do
  end subroutine adv_fc

  !> @brief  Calculate TVD
  !!
  !! Use explicit scheme to calculate the advective component of transport.
  !! TVD is an acronym for Total-Variation Diminishing
  !<
  subroutine advtvd(this, n, cnew, rhs)
    ! -- modules
    ! -- dummy
    class(TspAdvType) :: this
    integer(I4B), intent(in) :: n
    real(DP), dimension(:), intent(in) :: cnew
    real(DP), dimension(:), intent(inout) :: rhs
    ! -- local
    real(DP) :: qtvd
    integer(I4B) :: m, ipos
    !
    ! -- Loop through each n connection.  This will
    do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
      if (this%dis%con%mask(ipos) == 0) cycle
      m = this%dis%con%ja(ipos)
      if (m > n .and. this%ibound(m) /= 0) then
        ! qtvd = this%advqtvd(n, m, ipos, cnew)
        qtvd = this%advqtvd_experimental(n, m, ipos, cnew)
        rhs(n) = rhs(n) - qtvd
        rhs(m) = rhs(m) + qtvd
      end if
    end do
  end subroutine advtvd

  !> @brief  Calculate TVD
  !!
  !! Use explicit scheme to calculate the advective component of transport.
  !! TVD is an acronym for Total-Variation Diminishing
  !<
  function advqtvd(this, n, m, iposnm, cnew) result(qtvd)
    ! -- modules
    use ConstantsModule, only: DPREC
    ! -- return
    real(DP) :: qtvd
    ! -- dummy
    class(TspAdvType) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: m
    integer(I4B), intent(in) :: iposnm
    real(DP), dimension(:), intent(in) :: cnew
    ! -- local
    integer(I4B) :: ipos, isympos, iup, idn, i2up, j
    real(DP) :: qnm, qmax, qupj, elupdn, elup2up
    real(DP) :: smooth, cdiff, alimiter
    integer(I4B) :: ihc
    integer(I4B) :: ihc2
    real(DP) :: nx, ny, nz
    real(DP) :: nx2, ny2, nz2
    !
    ! -- initialize
    qtvd = DZERO
    !
    ! -- Find upstream node
    isympos = this%dis%con%jas(iposnm)
    qnm = this%fmi%gwfflowja(iposnm)
    if (qnm > DZERO) then
      ! -- positive flow into n means m is upstream
      iup = m
      idn = n
    else
      iup = n
      idn = m
    end if
    elupdn = this%dis%con%cl1(isympos) + this%dis%con%cl2(isympos)
    !
    ! -- Find connection direction
    ihc = this%dis%con%ihc(isympos)
    call this%dis%connection_normal(n, m, ihc, nx, ny, nz, ipos)
    !
    ! -- Find second node upstream to iup
    i2up = 0
    qmax = DZERO
    do ipos = this%dis%con%ia(iup) + 1, this%dis%con%ia(iup + 1) - 1
      j = this%dis%con%ja(ipos)
      isympos = this%dis%con%jas(ipos)
      ihc2 = this%dis%con%ihc(isympos)
      call this%dis%connection_normal(iup, j, ihc, nx2, ny2, nz2, ipos)

      if (((abs(nx) - abs(nx2)) > 1e-5) .or. (abs(ny) - abs(ny2) > 1e-5)) cycle

      if (this%ibound(j) == 0) cycle
      qupj = this%fmi%gwfflowja(ipos)

      if (qupj > qmax) then
        qmax = qupj
        i2up = j
        elup2up = this%dis%con%cl1(isympos) + this%dis%con%cl2(isympos)
      end if
    end do
    !
    ! -- Calculate flux limiting term
    if (i2up > 0) then
      smooth = DZERO
      cdiff = ABS(cnew(iup) - cnew(i2up)) ! adjusted code
      ! cdiff = ABS(cnew(idn) - cnew(iup)) ! original code
      if (cdiff > DPREC) then
        smooth = (cnew(idn) - cnew(iup)) / elupdn * &
                 elup2up / (cnew(iup) - cnew(i2up)) ! adjusted code
        ! smooth = (cnew(iup) - cnew(i2up)) / elup2up * &
        !          elupdn / (cnew(idn) - cnew(iup)) ! original code
      end if

      if (smooth > DZERO) then
        alimiter = this%limiter(smooth)

        qtvd = DHALF * alimiter * qnm * (cnew(iup) - cnew(i2up)) ! adjusted code
        ! qtvd = DHALF * alimiter * qnm * (cnew(idn) - cnew(iup)) ! original code
        qtvd = qtvd * this%eqnsclfac
      end if
    end if
  end function advqtvd

  function limiter(this, r, dplus, dmin) result(theta)
    ! -- return
    real(DP) :: theta ! limited slope
    ! -- dummy
    class(TspAdvType) :: this
    real(DP) :: r ! ratio of successive gradients
    real(DP), OPTIONAL :: dplus, dmin
    ! -- local
    real(DP) :: sign_dplus, H3, H3L, eta, alpha, cellsize, eps, beta, gamma, H3L_2

    select case (this%iadvwt)
    case (2) ! van Leer
      ! alimiter = DTWO * smooth / (DONE + smooth)
      theta = max(0.0_dp, min((r + dabs(r)) / (1.0_dp + dabs(r)), 2.0_dp))
    case (3) ! Koren
        theta = max(0.0_dp, min(2.0_dp * r, 1.0_dp / 3.0_dp + 2.0_dp / 3.0_dp * r, 2.0_dp))
    case (4) ! Superbee
      theta = max(0.0_dp, min(2.0_dp * r, 1.0_dp), min(r, 2.0_dp))
    case (5) ! van Albada
      theta = max(0.0_dp, (r * r + r) / (r * r + 1.0_dp))
    case (6) ! Koren modified
        ! theta = max(0.0_dp, min(4.0_dp * r * r + r, 1.0_dp/3.0_dp + 2.0_dp/3.0_dp * r, 2.0_dp))
      if (.not. present(dplus)) then
        theta = DZERO
      else  !-- 3th order
        sign_dplus = sign(1.0_dp, dplus)
        eps = 1e-6

        alpha = 2.0_dp
        cellsize = 17.5438596491228

        eta = sqrt(dmin**2 + dplus**2) / ((sqrt(5.0_dp/2.0_dp) * alpha * cellsize))
        H3 = (2.0_dp + r) / 3.0_dp
        H3L = sign_dplus * max(0.0_dp, min(sign_dplus * H3, &
          max(-sign_dplus * dmin, min(2.0_dp * sign_dplus * dmin, sign_dplus * H3, 1.5_dp * DABS(dplus)))))

        alpha = 1.0_dp
        beta = 2.0_dp
        gamma = 1.5_dp
        H3L_2 = max(0.0_dp, min(H3, max(-alpha * r, min(beta * r, H3, gamma))))
        if (eta < 0.01_dp - eps) then
          theta = H3
        else if (eta > 0.01_dp + eps) then
          theta = H3L_2
        else
          theta = 0.5_dp * ((1.0_dp - (eta - 1.0_dp)/eps) * H3 + (1.0_dp - (eta + 1.0_dp)/eps) * H3L_2)
        end if
        theta = H3L_2
        ! theta = theta / dplus
      end if
    CASE DEFAULT
      theta = DZERO
    end select

  end function

  function advqtvd_experimental(this, n, m, iposnm, cnew) result(qtvd)
    ! -- return
    real(DP) :: qtvd
    ! -- dummy
    class(TspAdvType) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: m
    integer(I4B), intent(in) :: iposnm
    real(DP), dimension(:), intent(in) :: cnew
    ! -- local
    integer(I4B) :: iup, idn, isympos
    real(DP) :: qnm, cup2
    real(DP), dimension(3) :: grad_c, dnm, grad1_c, grad2_c, grad2xyz_c, dxyz, dnu
    real(DP) :: smooth, alimiter, beta, smooth2, smooth3
    real(DP) :: cl1, cl2, rel_dist, c_virtual
    real(DP) :: dplus, dmin
    real(DP) :: dxplus, dxmin_half, dxplus_half
    integer(I4B) :: nnodes, number_connections
    real(DP), allocatable :: polyverts(:, :)

    !
    ! -- initialize
    qtvd = DZERO
    !
    ! -- Find upstream node
    isympos = this%dis%con%jas(iposnm)
    qnm = this%fmi%gwfflowja(iposnm)
    if (qnm > DZERO) then
      ! -- positive flow into n means m is upstream
      iup = m
      idn = n

      cl1 = this%dis%con%cl2(isympos)
      cl2 = this%dis%con%cl1(isympos)
    else
      iup = n
      idn = m

      cl1 = this%dis%con%cl1(isympos)
      cl2 = this%dis%con%cl2(isympos)
    end if
    !
    ! -- Return if straddled cells have same value
    if (abs(cnew(idn) - cnew(iup)) < 1e-8_dp) return
    !
    ! -- Return if upstream cell is a boundary
    ! call this%fmi%dis%get_polyverts(iup, polyverts)
    ! nnodes = size(polyverts, dim = 2)
    ! number_connections = this%dis%con%ia(iup + 1) - this%dis%con%ia(iup) - 1
    ! if (number_connections < nnodes) return

    ! call this%fmi%dis%get_polyverts(idn, polyverts)
    ! nnodes = size(polyverts, dim = 2)
    ! number_connections = this%dis%con%ia(idn + 1) - this%dis%con%ia(idn) - 1
    ! if (number_connections < nnodes) return

    !
    ! -- Compute cell concentration gradient
    call this%compute_cell_gradient(iup, cnew, grad_c)
    call this%compute_cell_gradient_2dorder(iup, cnew, grad1_c, grad2_c, grad2xyz_c)
    !
    ! -- Compute smoothness factor
    dnm = this%node_distance(iup, idn)
    smooth = 2.0_dp * (dot_product(grad_c, dnm)) / (cnew(idn) - cnew(iup)) - 1.0_dp
    !
    ! -- Correct smoothness factor to prevent negative concentration
    c_virtual = cnew(iup) - smooth * (cnew(idn) - cnew(iup))
    if (c_virtual <= DPREC) then
      smooth = cnew(iup) / (cnew(idn) - cnew(iup))
    end if

    dnu = -dnm 
    dxyz(1) = dnu(1) * dnu(2)
    dxyz(2) = dnu(2) * dnu(3)
    dxyz(3) = dnu(3) * dnu(1)
    cup2 = max(0.0_dp, cnew(iup) + dot_product(grad1_c, dnu) + dot_product(grad2_c, 0.5_dp * dnu**2) &
      + dot_product(grad2xyz_c, dxyz))
    smooth2 = (cnew(iup) - cup2)/(cnew(idn) - cnew(iup))
  
    c_virtual = cnew(iup) - smooth2 * (cnew(idn) - cnew(iup))
    if (c_virtual <= DPREC) then
      smooth2 = cnew(iup) / (cnew(idn) - cnew(iup))
    end if
    !
    ! -- Compute relative distance to face
    rel_dist = cl1 / (cl1 + cl2)
    !
    ! -- Compute limiter
    dplus = cnew(idn) - cnew(iup)
    dmin = smooth * dplus
    dxplus = cl2 * 2.0
    dxmin = cl2 * 2.0
    dxplus_half = cl1 + cl2
    dxmin_half = cl1 + cl2
    
    alimiter = this%limiter(smooth, dxmin_half * dplus / dxplus_half, dxplus_half * dmin / dxmin_half)
    ! -- Compute limited flux
    qtvd = 0.5 * alimiter * qnm * dplus
    qtvd = qtvd * this%eqnsclfac

  end function advqtvd_experimental

  function node_distance(this, n, m) result(d)
    ! -- return
    real(DP), dimension(3) :: d
    ! -- dummy
    class(TspAdvType) :: this
    integer(I4B), intent(in) :: n, m
    ! -- local
    real(DP) :: x_dir, y_dir, z_dir, length

    call this%dis%connection_vector(n, m, .true., 1.0_dp, 1.0_dp, 1, x_dir, y_dir, z_dir, length)
    d(1) = x_dir * length
    d(2) = y_dir * length
    d(3) = z_dir * length

  end function node_distance

  subroutine compute_cell_gradient_2dorder(this, n, cnew, grad_c, grad2_c, grad2xyz_c)
    ! -- dummy
    class(TspAdvType) :: this
    integer(I4B), intent(in) :: n
    real(DP), dimension(:), intent(in) :: cnew
    real(DP), dimension(3), intent(out) :: grad_c
    real(DP), dimension(3), intent(out) :: grad2_c
    real(DP), dimension(3), intent(out) :: grad2xyz_c
    ! -- local
    integer(I4B) :: ipos, local_pos, m
    integer(I4B) :: number_connections

    real(DP), dimension(3) :: dnm
    real(DP), dimension(:, :), allocatable :: d, d_trans
    real(DP), dimension(:, :), allocatable :: g, g_inv
    real(DP), dimension(:, :), allocatable :: grad_op
    real(DP), dimension(9) :: grad

    real(DP), dimension(:), allocatable :: dc

    number_connections = this%dis%con%ia(n + 1) - this%dis%con%ia(n) - 1
    if (number_connections == 1) then
      ! If a cell only has 1 neigbour compute the gradient using finite difference
      ! This case can happen if a triangle element is located in a cornor of a square domain
      ! with two sides being domain boundaries
      ipos = this%dis%con%ia(n) + 1
      m = this%dis%con%ja(ipos)
      dnm = this%node_distance(n, m)

      grad_c(1) = (cnew(m) - cnew(n)) / dnm(1)
      grad_c(2) = (cnew(m) - cnew(n)) / dnm(2)
      grad_c(3) = (cnew(m) - cnew(n)) / dnm(3)

      grad2_c(1) = 0
      grad2_c(2) = 0
      grad2_c(3) = 0

      grad2xyz_c(1) = 0
      grad2xyz_c(2) = 0
      grad2xyz_c(3) = 0
      return
    end if

    allocate (d(number_connections, 9))

    local_pos = 1
    do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
      m = this%dis%con%ja(ipos)
      dnm = this%node_distance(n, m)

      d(local_pos, 1) = dnm(1)
      d(local_pos, 2) = dnm(2)
      d(local_pos, 3) = dnm(3)
      d(local_pos, 4) = 0.5_dp*dnm(1)**2
      d(local_pos, 5) = 0.5_dp*dnm(2)**2
      d(local_pos, 6) = 0.5_dp*dnm(3)**2
      d(local_pos, 7) = dnm(1)*dnm(2)
      d(local_pos, 8) = dnm(2)*dnm(3)
      d(local_pos, 9) = dnm(3)*dnm(1)

      local_pos = local_pos + 1
    end do

    d_trans = transpose(d)
    g = matmul(d_trans, d)
    g_inv = pinv(g)
    grad_op = matmul(g_inv, d_trans)

    ! Assemble the concentration difference vector
    allocate (dc(number_connections))
    local_pos = 1
    do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
      m = this%dis%con%ja(ipos)
      dc(local_pos) = cnew(m) - cnew(n)
      local_pos = local_pos + 1
    end do

    grad = matmul(grad_op, dc)
    grad_c(1) = grad(1)
    grad_c(2) = grad(2)
    grad_c(3) = grad(3)

    grad2_c(1) = grad(4)
    grad2_c(2) = grad(5)
    grad2_c(3) = grad(6)

    grad2xyz_c(1) = grad(7)
    grad2xyz_c(2) = grad(8)
    grad2xyz_c(3) = grad(9)

  end subroutine

  subroutine compute_cell_gradient(this, n, cnew, grad_c)
    ! -- dummy
    class(TspAdvType) :: this
    integer(I4B), intent(in) :: n
    real(DP), dimension(:), intent(in) :: cnew
    real(DP), dimension(3), intent(out) :: grad_c
    ! -- local
    integer(I4B) :: ipos, local_pos
    integer(I4B) :: number_connections
    real(DP), dimension(:, :), allocatable :: d
    real(DP), dimension(:, :), allocatable :: d_trans
    real(DP), dimension(:, :), allocatable :: W2
    real(DP), dimension(:, :), allocatable :: grad_op
    real(DP), dimension(3, 3) :: g
    real(DP), dimension(3, 3) :: g_inv
    integer(I4B) :: m
    real(DP), dimension(3) :: dnm

    real(DP), dimension(:), allocatable :: dc

    number_connections = this%dis%con%ia(n + 1) - this%dis%con%ia(n) - 1
    if (number_connections == 1) then
      ! If a cell only has 1 neigbour compute the gradient using finite difference
      ! This case can happen if a triangle element is located in a cornor of a square domain
      ! with two sides being domain boundaries
      ipos = this%dis%con%ia(n) + 1
      m = this%dis%con%ja(ipos)
      dnm = this%node_distance(n, m)

      grad_c(1) = (cnew(m) - cnew(n)) / dnm(1)
      grad_c(2) = (cnew(m) - cnew(n)) / dnm(2)
      grad_c(3) = (cnew(m) - cnew(n)) / dnm(3)
      return
    end if

    allocate (d(number_connections, 3))
    allocate (d_trans(3, number_connections))
    allocate (grad_op(3, number_connections))
    allocate (W2(number_connections, number_connections))

    ! Assemble the distance and transposed distance matrices
    W2 = 0
    local_pos = 1
    do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
      m = this%dis%con%ja(ipos)
      dnm = this%node_distance(n, m)

      d(local_pos, 1) = dnm(1)
      d(local_pos, 2) = dnm(2)
      d(local_pos, 3) = dnm(3)

      d_trans(1, local_pos) = d(local_pos, 1)
      d_trans(2, local_pos) = d(local_pos, 2)
      d_trans(3, local_pos) = d(local_pos, 3)

      W2(local_pos, local_pos) = 1.0_dp / (dnm(1)**2 + dnm(2)**2 + dnm(3)**2)

      local_pos = local_pos + 1
    end do

    ! Compute the G and inverse G matrices
    g = matmul(d_trans, matmul(W2, d))
    g_inv = pinv(g)

    ! Compute the gradient operator
    grad_op = matmul(g_inv, matmul(d_trans, W2))

    ! Assemble the concentration difference vector
    allocate (dc(number_connections))
    local_pos = 1
    do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
      m = this%dis%con%ja(ipos)
      dc(local_pos) = cnew(m) - cnew(n)
      local_pos = local_pos + 1
    end do

    ! Compute the cells gradient
    grad_c = matmul(grad_op, dc)

  end subroutine compute_cell_gradient

  function pinv(A) result(B)
    real(DP), intent(in) :: A(:, :) !! Matrix
    real(DP) :: B(SIZE(A, DIM=2), SIZE(A, DIM=1)) !! Inverse matrix

    integer(I4B) :: pos, ierr
    real(DP), dimension(SIZE(A, DIM=1), SIZE(A, DIM=1)) :: U
    real(DP), dimension(SIZE(A, DIM=2), SIZE(A, DIM=2)) :: V
    real(DP), dimension(SIZE(A, DIM=2)) :: sigma
    real(DP), dimension(SIZE(A, DIM=2), SIZE(A, DIM=1)) :: sigma_plus

    CALL SVD(A, sigma, .TRUE., U, .TRUE., V, ierr)

    sigma_plus = 0
    do pos = 1, min(SIZE(A, DIM=1), SIZE(A, DIM=2))
      if (DABS(sigma(pos)) > 2.0_dp * DPREC) then
        sigma_plus(pos, pos) = 1.0_dp / sigma(pos)
      end if
    end do

    B = matmul(V, matmul(sigma_plus, transpose(U)))

  end function

  !> @brief Calculate advection contribution to flowja
  !<
  subroutine adv_cq(this, cnew, flowja)
    ! -- modules
    use TdisModule, only: kstp, kper, delt
    ! -- dummy
    class(TspAdvType) :: this
    real(DP), intent(in), dimension(:) :: cnew
    real(DP), intent(inout), dimension(:) :: flowja
    ! -- local
    integer(I4B) :: nodes
    integer(I4B) :: n, m, idiag, ipos
    real(DP) :: omega, qnm, q
    integer(I4B) time_idx
    real(DP), allocatable, dimension(:) :: cnew2, rhs2, rhs2_old
    real(DP), dimension(3):: rk3_K_weights =(/0.0_dp, 0.5_dp, 2.0_dp /)
    real(DP), dimension(3):: rk3_rhs_weights =(/1.0_dp/6.0_dp, 4.0_dp/6.0_dp, 1.0_dp/6.0_dp /)
    real(DP), dimension(3):: rk3_rhs_old_weights =(/0.0_dp, 0.0_dp, -1.0_dp  /)
    !
    ! -- Calculate advection and add to flowja. qnm is the volumetric flow
    !    rate and has dimensions of L^/T.
    allocate(rhs2(size(flowja)))
    rhs2_old = rhs2
    do time_idx = 1, 3
      cnew2 = cnew + delt * (rk3_K_weights(time_idx) * rhs2 + rk3_rhs_old_weights(time_idx) * rhs2_old)
      rhs2_old = rhs2
      rhs2 = 0
      
      nodes = this%dis%nodes
      do n = 1, nodes
        if (this%ibound(n) == 0) cycle
        idiag = this%dis%con%ia(n)
        do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
          m = this%dis%con%ja(ipos)
          if (m <= n) cycle
          if (this%ibound(m) == 0) cycle
          qnm = this%fmi%gwfflowja(ipos) * this%eqnsclfac
          omega = this%adv_weight(this%iadvwt, ipos, n, m, qnm)
          if (qnm > 0) then
            q = qnm * cnew2(m)
          else
            q = qnm * cnew2(n)
          end if
          rhs2(idiag) = rhs2(idiag) - q
          rhs2(ipos) = rhs2(ipos) + q
        end do
      end do
      !
      ! -- TVD
      if (this%iadvwt >= 2) call this%advtvd_bd(cnew2, rhs2)
      flowja = flowja + rk3_rhs_weights(time_idx) * rhs2
    end do
  end subroutine adv_cq

  !> @brief Add TVD contribution to flowja
  subroutine advtvd_bd(this, cnew, flowja)
    ! -- modules
    ! -- dummy
    class(TspAdvType) :: this
    real(DP), dimension(:), intent(in) :: cnew
    real(DP), dimension(:), intent(inout) :: flowja
    ! -- local
    real(DP) :: qtvd, qnm
    integer(I4B) :: nodes, n, m, ipos
    !
    nodes = this%dis%nodes
    do n = 1, nodes
      if (this%ibound(n) == 0) cycle
      do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
        m = this%dis%con%ja(ipos)
        if (this%ibound(m) /= 0) then
          qnm = this%fmi%gwfflowja(ipos)
          ! qtvd = this%advqtvd(n, m, ipos, cnew)
          qtvd = this%advqtvd_experimental(n, m, ipos, cnew)
          flowja(ipos) = flowja(ipos) + qtvd
        end if
      end do
    end do
  end subroutine advtvd_bd

  !> @brief Deallocate memory
  !<
  subroutine adv_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(TspAdvType) :: this
    !
    ! -- Deallocate arrays if package was active
    if (this%inunit > 0) then
    end if
    !
    ! -- nullify pointers
    this%ibound => null()
    !
    ! -- Scalars
    call mem_deallocate(this%iadvwt)
    call mem_deallocate(this%ats_percel)
    !
    ! -- deallocate parent
    call this%NumericalPackageType%da()
  end subroutine adv_da

  !> @brief Allocate scalars specific to the streamflow energy transport (SFE)
  !! package.
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr
    ! -- dummy
    class(TspAdvType) :: this
    ! -- local
    !
    ! -- allocate scalars in NumericalPackageType
    call this%NumericalPackageType%allocate_scalars()
    !
    ! -- Allocate
    call mem_allocate(this%iadvwt, 'IADVWT', this%memoryPath)
    call mem_allocate(this%ats_percel, 'ATS_PERCEL', this%memoryPath)
    !
    ! -- Initialize
    this%iadvwt = 0
    this%ats_percel = DNODATA
    !
    ! -- Advection creates an asymmetric coefficient matrix
    this%iasym = 1
  end subroutine allocate_scalars

  !> @brief Read options
  !!
  !! Read the options block
  !<
  subroutine read_options(this)
    ! -- modules
    use ConstantsModule, only: LINELENGTH
    use SimModule, only: store_error
    ! -- dummy
    class(TspAdvType) :: this
    ! -- local
    character(len=LINELENGTH) :: errmsg, keyword
    integer(I4B) :: ierr
    logical :: isfound, endOfBlock
    ! -- formats
    character(len=*), parameter :: fmtiadvwt = &
      &"(4x,'ADVECTION WEIGHTING SCHEME HAS BEEN SET TO: ', a)"
    !
    ! -- get options block
    call this%parser%GetBlock('OPTIONS', isfound, ierr, blockRequired=.false., &
                              supportOpenClose=.true.)
    !
    ! -- parse options block if detected
    if (isfound) then
      write (this%iout, '(1x,a)') 'PROCESSING ADVECTION OPTIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('SCHEME')
          call this%parser%GetStringCaps(keyword)
          select case (keyword)
          case ('UPSTREAM')
            this%iadvwt = 0
            write (this%iout, fmtiadvwt) 'UPSTREAM'
          case ('CENTRAL')
            this%iadvwt = 1
            write (this%iout, fmtiadvwt) 'CENTRAL'
          case ('TVD')
            this%iadvwt = 2
            write (this%iout, fmtiadvwt) 'TVD'
          case ('KOREN')
            this%iadvwt = 3
            write (this%iout, fmtiadvwt) 'KOREN'
          case ('SUPERBEE')
            this%iadvwt = 4
            write (this%iout, fmtiadvwt) 'SUPERBEE'
          case ('ALBADA')
            this%iadvwt = 5
            write (this%iout, fmtiadvwt) 'ALBADA'
          case ('KORENMOD')
            this%iadvwt = 6
            write (this%iout, fmtiadvwt) 'KORENMOD'
          case default
            write (errmsg, '(a, a)') &
              'Unknown scheme: ', trim(keyword)
            call store_error(errmsg)
            write (errmsg, '(a, a)') &
              'Scheme must be "UPSTREAM", "CENTRAL" or "TVD"'
            call store_error(errmsg)
            call this%parser%StoreErrorUnit()
          end select
        case ('ATS_PERCEL')
          this%ats_percel = this%parser%GetDouble()
          if (this%ats_percel == DZERO) this%ats_percel = DNODATA
          write (this%iout, '(4x,a,1pg15.6)') &
            'User-specified fractional cell distance for adaptive time &
            &steps: ', this%ats_percel
        case default
          write (errmsg, '(a,a)') 'Unknown ADVECTION option: ', &
            trim(keyword)
          call store_error(errmsg, terminate=.TRUE.)
        end select
      end do
      write (this%iout, '(1x,a)') 'END OF ADVECTION OPTIONS'
    end if
  end subroutine read_options

  !> @ brief Advection weight
  !!
  !! Calculate the advection weight
  !<
  function adv_weight(this, iadvwt, ipos, n, m, qnm) result(omega)
    ! -- return
    real(DP) :: omega
    ! -- dummy
    class(TspAdvType) :: this
    integer, intent(in) :: iadvwt
    integer, intent(in) :: ipos
    integer, intent(in) :: n
    integer, intent(in) :: m
    real(DP), intent(in) :: qnm
    ! -- local
    real(DP) :: lnm, lmn

    select case (iadvwt)
    case (1)
      ! -- calculate weight based on distances between nodes and the shared
      !    face of the connection
      if (this%dis%con%ihc(this%dis%con%jas(ipos)) == 0) then
        ! -- vertical connection; assume cell is fully saturated
        lnm = DHALF * (this%dis%top(n) - this%dis%bot(n))
        lmn = DHALF * (this%dis%top(m) - this%dis%bot(m))
      else
        ! -- horizontal connection
        lnm = this%dis%con%cl1(this%dis%con%jas(ipos))
        lmn = this%dis%con%cl2(this%dis%con%jas(ipos))
      end if
      omega = lmn / (lnm + lmn)
    case (0, 2:)
      ! -- use upstream weighting for upstream and tvd schemes
      if (qnm > DZERO) then
        omega = DZERO
      else
        omega = DONE
      end if
    end select
  end function adv_weight

end module TspAdvModule
