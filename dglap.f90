! ==============================================================================
! dglap.f90
!
! by Adam Freese <maxwellsdemon137@gmail.com> copyright 2017
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! This file is part of dglap90.
!
! dglap90 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dglap90 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with dglap90.  If not, see <http://www.gnu.org/licenses/>.
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! This file contains the modules:
!   - dglap_splitting_lo
!     :: Leading order splitting functions
!   - dglap_splitting_nlo
!     :: NLO splitting functions (TODO)
!   - dglap
!     :: Routines to perform DGLAP evolution
!
! Brute-force DGLAP evolution code, similar to the evolution code by Kumano
! and Miyama. There are a few differences I think justify this program's
! separate existence, however.
!
! The key differences from Kumano's and Miyama's code are:
!   - The x mesh changes from log to linear to reverse-log, so accurate
!     evolution obtains with a smaller grid. This is also necessary for
!     accurate evolution at high x.
!   - The maximum x value can be made greater than 1, which is needed for
!     studying PDFs in short range correlations.
!   - The interface is simpler.
!
! To get to the DGLAP module, search for xyzDGLAP
!
! The NLO splitting functions are not yet implemented.

! ==============================================================================
! //////////////////////////////////////////////////////////////////////////////
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! LO splitting function module follows [ xyzLO ]
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! //////////////////////////////////////////////////////////////////////////////
! ==============================================================================

module dglap_splitting_lo
  ! Module containing leading order splitting functions for DGLAP evolution
  use kinds, only : dp
  use qcd
  use runningQCD
  implicit none
  private

  public :: P_qq_0, P_qG_0, P_Gq_0, P_GG_0

contains

  ! ============================================================================
  ! ============================================================================

  function P_qq_0(x, pcode) result(res)
    ! for pcode,
    !  'r' :: regular part
    !  'p' :: plus prescription part (excluding (1-x) denominator)
    !  'd' :: delta part
    implicit none
    ! I/O
    real(dp), intent(in) :: x
    character(len=1), intent(in) :: pcode
    real(dp) :: res
    ! Switch for pcode
    select case(pcode)
    case('r')
      res = 0.0_dp
    case('p')
      res = CF * (1.0_dp + x**2)
    case('d')
      res = CF * 1.5_dp
    case default
      write(*,*) "Invalid pcode in splitting function"
      stop
    end select
    return
  end function P_qq_0

  ! ============================================================================
  ! ============================================================================

  function P_qG_0(x, pcode) result(res)
    ! for pcode,
    !  'r' :: regular part
    !  'p' :: plus prescription part (excluding (1-x) denominator)
    !  'd' :: delta part
    implicit none
    ! I/O
    real(dp), intent(in) :: x
    character(len=1), intent(in) :: pcode
    real(dp) :: res
    ! Switch for pcode
    select case(pcode)
    case('r')
      ! The factor of 2.0_dp*get_nFlavors() comes from the use of this in
      ! singlet evolution; we actually generate 2*nFl quarks instead of one.
      res = TF * ( x**2 + (1.0_dp - x)**2 )  * 2.0_dp*get_nFlavors()
    case('p')
      res = 0.0_dp
    case('d')
      res = 0.0_dp
    case default
      write(*,*) "Invalid pcode in splitting function"
      stop
    end select
    return
  end function P_qG_0

  ! ============================================================================
  ! ============================================================================

  function P_Gq_0(x, pcode) result(res)
    ! for pcode,
    !  'r' :: regular part
    !  'p' :: plus prescription part (excluding (1-x) denominator)
    !  'd' :: delta part
    implicit none
    ! I/O
    real(dp), intent(in) :: x
    character(len=1), intent(in) :: pcode
    real(dp) :: res
    ! Switch for pcode
    select case(pcode)
    case('r')
      res = CF * ( 1.0_dp + (1.0_dp - x)**2) / x
    case('p')
      res = 0.0_dp
    case('d')
      res = 0.0_dp
    case default
      write(*,*) "Invalid pcode in splitting function"
      stop
    end select
    return
  end function P_Gq_0

  ! ============================================================================
  ! ============================================================================

  function P_GG_0(x, pcode) result(res)
    ! for pcode,
    !  'r' :: regular part
    !  'p' :: plus prescription part (excluding (1-x) denominator)
    !  'd' :: delta part
    implicit none
    ! I/O
    real(dp), intent(in) :: x
    character(len=1), intent(in) :: pcode
    real(dp) :: res
    ! Switch for pcode
    select case(pcode)
    case('r')
      res = 2.0_dp * CA * ( (1.0_dp-x)/x + x*(1.0_dp-x) )
    case('p')
      res = 2.0_dp * CA * x
    case('d')
      res = 11.0_dp/6.0_dp * CA - 2.0_dp/3.0_dp * TF*get_nFlavors()
    case default
      write(*,*) "Invalid pcode in splitting function"
      stop
    end select
    return
  end function P_GG_0

  ! ============================================================================
  ! ============================================================================

end module dglap_splitting_lo

! ==============================================================================
! //////////////////////////////////////////////////////////////////////////////
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! NLO splitting function module follows [ xyzNLO ]
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! //////////////////////////////////////////////////////////////////////////////
! ==============================================================================

module dglap_splitting_nlo
  ! Module to contain next-to-eading-order splitting functions.
  ! TODO
  use kinds, only : dp
  use qcd
  use runningQCD
  implicit none
  private

!  public ::

contains

  ! ============================================================================
  ! ============================================================================



  ! ============================================================================
  ! ============================================================================

end module dglap_splitting_nlo

! ==============================================================================
! //////////////////////////////////////////////////////////////////////////////
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! DGLAP module follows [ xyzDGLAP ]
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! //////////////////////////////////////////////////////////////////////////////
! ==============================================================================

module dglap
  use kinds, only : dp
  use pointspace
  use qcd
  use runningQCD
  use dglap_splitting_lo
  use interpolation, only : interpolate_1D
  use dataframes

  implicit none
  private

  ! Some constants
  integer, parameter :: nFlavors = 6
  real(dp), parameter :: xlow  = 0.1_dp
  real(dp), parameter :: xhigh = 0.8_dp

  ! Evolution parameters that can be set by user
  integer :: Nx
  real(dp) :: xmin, xmax
  real(dp) :: Q2step

  ! Current Q2 value encoded by these
  real(dp) :: Q2, t

  ! Allocatable arrays to contain PDF info
  real(dp), dimension(:,:), allocatable :: pdf_mesh
  real(dp), dimension(:), allocatable :: x_mesh
  real(dp), dimension(:,:), allocatable :: q_minus, TNS
  real(dp), dimension(:), allocatable :: Sigma, Gluons

  ! Parton names for reading/writing CSV files
  character(len=2), dimension(-nFlavors:nFlavors), parameter :: &
      parton_names = &
      [ 'tb', 'bb', 'sb', 'cb', 'db', 'ub', 'g ', &
        'u ', 'd ', 's ', 'c ', 'b ', 't ' ]

  ! Some string manipulation for debugging
  character(len=120), parameter :: &
      debug_rrr = "(A20,ES20.10,ES20.10,ES20.10)", &
      debug_rrl = "(A20,ES20.10,ES20.10,L20)"

  ! Public methods
  public :: dglap_init, dglap_cleanup, &
      dglap_set_pdf, dglap_eval_pdf, &
      dglap_read_csv, dglap_write_csv, &
      dglap_evolve

  ! TODO :: file reading/writing, NLO

contains

  ! ============================================================================
  ! ////////////////////////////////////////////////////////////////////////////
  ! Public methods
  ! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  ! ============================================================================

  subroutine dglap_init( opt_Nx, opt_xmin, opt_xmax, opt_Q2step )
    ! Initializes the DGLAP module.
    ! Using this is optional: it's called by dglap_set_pdf anyway.
    ! However, the user can set non-default values for evolution parameters
    ! by calling this method.
    implicit none

    ! I/O
    integer, intent(in), optional :: opt_Nx
    real(dp), intent(in), optional :: opt_xmin, opt_xmax, opt_Q2step

    ! Default values for evolution parameters
    Nx     = 300
    xmin   = 1e-3_dp
    xmax   = 1.0_dp
    Q2step = 1.2_dp

    ! Process optional input
    if( present(opt_Nx) )     Nx     = opt_Nx
    if( present(opt_xmin) )   xmin   = opt_xmin
    if( present(opt_xmax) )   xmax   = opt_xmax
    if( present(opt_Q2step) ) Q2step = opt_Q2step

    ! In case user is changing options, call cleanup first
    call dglap_cleanup()

    ! Now allocate arrays based on parameters
    if( .not.allocated(x_mesh) )   allocate( x_mesh(Nx) )
    if( .not.allocated(pdf_mesh) ) allocate( pdf_mesh(-nFlavors:nFlavors, Nx) )

    ! And call the routine to initialize the x mesh
    call init_x_mesh()

    ! Let the RunningQCD module know we're doing LO, so alphaQCD is set right.
    call set_QCD_order(1)

    ! Be done
    return
  end subroutine dglap_init

  ! ============================================================================
  ! ============================================================================

  subroutine dglap_cleanup()
    ! Call this method when done with evolution, to clean up used memory.
    implicit none

    ! Deallocate memory that should be allocated if module is intialized.
    if( allocated(x_mesh) )   deallocate( x_mesh )
    if( allocated(pdf_mesh) ) deallocate( pdf_mesh )

    ! Just in case, deallocate anything else allocated that shouldn't be.
    if( allocated(q_minus) )  deallocate( q_minus )
    if( allocated(TNS) )      deallocate( TNS )
    if( allocated(Gluons) )   deallocate( Gluons )
    if( allocated(Sigma) )    deallocate( Sigma )

    return
  end subroutine dglap_cleanup

  ! ============================================================================
  ! ============================================================================

  subroutine dglap_read_csv( csv_file, Q2_value )
    ! Reads a CSV file into the PDF grid
    ! The Q2 value is not expected to be contained in the CSV file,
    ! so must be supplied by the user.
    implicit none
    ! I/O
    character(len=*), intent(in) :: csv_file
    real(dp), intent(in) :: Q2_value

    ! Intermediates
    type(dataframe) :: df_pdfs
    integer :: iFl

    ! First, call dglap_init to initialize some defaults; then free allocations.
    call dglap_init()
    call dglap_cleanup()

    ! Set the Q2 value to the one supplied by the user
    call set_Q2( Q2_value )
    Q2 = Q2_value

    ! Fill the dataframe by reading from file
    call df_pdfs%read_csv( csv_file )

    ! Read x mesh from CSV file. Fortran allocates automagically.
    ! Then set internal x-related constants based on what we read.
    x_mesh = df_pdfs%get('x')
    Nx = size( x_mesh )
    xmin = minval( x_mesh )
    xmax = maxval( x_mesh )

    ! Allocate the PDF grid
    allocate( pdf_mesh(-nFlavors:nFlavors, size(x_mesh)) )

    ! Read PDF arrays from CSV file.
    do iFl=-nFlavors, nFlavors, 1
      pdf_mesh(iFl,:) = df_pdfs%get( parton_names(iFl), failsafe=nil )
    end do

    ! TODO : add public dealloc method to dataframe type

    return
  end subroutine dglap_read_csv

  ! ============================================================================
  ! ============================================================================

  subroutine dglap_write_csv( csv_file )
    ! Writes the PDF grid to a CSV file
    implicit none
    ! I/O
    character(len=*), intent(in) :: csv_file

    ! Intermediates
    type(dataframe) :: df_pdfs
    integer :: iFl

    ! Write the internal x mesh and PDF mesh to a dataframe
    call df_pdfs%add_col( 'x', x_mesh )

    do iFl=-nFlavors, nFlavors, 1
      call df_pdfs%add_col( parton_names(iFl), pdf_mesh(iFl,:) )
    end do

    ! Write dataframe to file
    call df_pdfs%write_csv( csv_file )

    ! TODO : add public dealloc method to dataframe type

    return
  end subroutine dglap_write_csv

  ! ============================================================================
  ! ============================================================================

  subroutine dglap_set_pdf( pdf_function, Q2_value )
    ! Uses a given function of a single variable, x, to set the PDF grids.
    ! The function should return an array indexed from -6 to 6.
    ! Q2 is also passed, and is used internally to keep track of alphaQCD
    ! and the effective number of flavors.
    implicit none
    ! I/O
    real(dp), intent(in) :: Q2_value

    ! Interface for the pdf function
    interface
      function pdf_function(x) result(f)
        use kinds, only : dp
        implicit none
        real(dp), intent(in) :: x
        real(dp), dimension(-6:6) :: f
      end function pdf_function
    end interface

    ! Iterator
    integer i

    ! If user called without using dglap_init() first, we do so ourselves
    if( .not. allocated(x_mesh) ) call dglap_init()

    ! Read function evaluations into the PDF grid
    do i=1, Nx, 1
      pdf_mesh(:,i) = pdf_function( x_mesh(i) )
    end do

    ! Internalize given Q2 value, and set it in the runningQCD module
    Q2 = Q2_value
    call set_Q2(Q2)

    return

  end subroutine dglap_set_pdf

  ! ============================================================================
  ! ============================================================================

  function dglap_eval_pdf(x, iFl) result(res)
    ! Evaluates pdf flavor iFl at current Q2 value and given x value
    implicit none
    ! I/O
    real(dp), intent(in) :: x
    integer, intent(in) :: iFl
    real(dp) :: res
    ! Use interpolation
    res = interpolate_1D( x_mesh, pdf_mesh(iFl,:), x )
    return
  end function dglap_eval_pdf

  ! ============================================================================
  ! ============================================================================

  subroutine dglap_evolve(newQ2)
    implicit none

    ! I/O
    real(dp), intent(in) :: newQ2

    ! Intermediates
    real(dp) :: Q2_temp, Q2_mul
    real(dp), dimension(3), parameter :: hqm2 = [ cMass2, bMass2, tMass2 ]
    integer :: i, evsign

    ! Q2_mul is determined by whether we're evolving up or down
    if(newQ2 > Q2) then
      Q2_mul = Q2step
      evsign = +1
    else if(newQ2 < Q2) then
      Q2_mul = Q2step**(-1)
      evsign = -1
    else
      ! Why are we evolving in this case?
      return
    end if

    ! Prepare the decoupled singlet / non-singlet arrays that evolve
    call brief_arrays()

    ! Loop
    do
      ! Prepare a temporary Q2 value
      Q2_temp = Q2 * Q2_mul

      ! Special care is taken close to quark mass thresholds
      do i=1, 3, 1
        if( compare(hqm2(i), Q2) .and. compare(Q2_temp, hqm2(i)) ) then
          ! Evolve just short of the threshold
          call evolve_step( hqm2(i) - evsign*tiny(1.0_dp) )
          ! Set Q2 just past the threshold, so set Lambda and effective nFl
          call set_Q2( hqm2(i) + evsign*tiny(1.0_dp) )
        end if
      end do

      ! We can leave loop if we reach the final target
      if( compare(Q2_temp, newQ2) ) then
        call evolve_step(newQ2)
        exit
      end if

      ! If we're still here, now evolve to the temporary Q2 value
      call evolve_step(Q2_temp)

    end do

    ! Debrief the evolved arrays to get back our pdf mesh
    call debrief_arrays()

    return
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    contains
      function compare( init_Q2, target_Q2 ) result(tof)
        ! Lambda function which compares an initial and target Q**2.
        ! If the direction of evolution is upward,
        ! then we get True when init_Q2 is at least target_Q2.
        ! If the direction is downward,
        ! we get True when init_Q2 is at most target_Q2.
        real(dp), intent(in) :: init_Q2, target_Q2
        logical :: tof
        if( Q2_mul > 1.0_dp) then
          tof = ( init_Q2 >= target_Q2 )
        else
          tof = ( init_Q2 <= target_Q2 )
        end if
        return
      end function compare
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  end subroutine dglap_evolve

  ! ============================================================================
  ! ////////////////////////////////////////////////////////////////////////////
  ! Private (module-internal) methods used for implementation
  ! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  ! ============================================================================

  subroutine init_x_mesh()
    ! Initializes an x array in three regions:
    !   1. Between xmin and xlow, logarithmic spacing
    !   2. Between xlow and xhigh, linear spacing
    !   3. Between xhigh and max, reverse log spacing
    implicit none
    ! Allocate the x_mesh (if not allocated already)
    if( .not. allocated(x_mesh) ) allocate( x_mesh(Nx) )
    ! Prepare each of the regions
    x_mesh(      1 :   Nx/3 ) =    logspace(xmin, xlow,  Nx/3)
    x_mesh(   Nx/3 : 2*Nx/3 ) =    linspace(xlow, xhigh, Nx/3+1)
    x_mesh( 2*Nx/3 :   Nx   ) = revlogspace(xhigh, xmax, Nx/3+1)
    return
  end subroutine init_x_mesh

  ! ============================================================================
  ! ============================================================================

  subroutine brief_arrays()
    ! Creates linear combinations of "physical" PDFs that (mostly) decouple
    ! in their evolution.
    implicit none

    ! q_plus is a temporary array since it's not used in the evolution
    real(dp), dimension(:,:), allocatable :: q_plus
    integer i, k

    ! Gluons
    allocate( Gluons(Nx) )
    Gluons(:) = pdf_mesh(0,:)

    ! q+ and q- are easiest
    allocate( q_plus(nFlavors, Nx) )
    allocate( q_minus(nFlavors, Nx) )
    forall(i=1:nFlavors)
        q_plus(i,:)  = pdf_mesh(i,:) + pdf_mesh(-i,:)
        q_minus(i,:) = pdf_mesh(i,:) - pdf_mesh(-i,:)
    end forall

    ! Sigma for singlet evolution
    allocate( Sigma(Nx) )
    forall(i=1:Nx)
        Sigma(i) = sum( q_plus(:,i) )
    end forall

    ! T, like q_minus, evovles by non-singlet evolution
    allocate( TNS(2:nFlavors, Nx) )
    forall(i=1:Nx)
        forall(k=2:nFlavors)
            TNS(k,i) = sum( q_plus(1:k,i) ) - k*q_plus(k,i)
        end forall
    end forall

    ! Garbage collection and return
    deallocate( q_plus )

    return

  end subroutine brief_arrays

  ! ============================================================================
  ! ============================================================================

  subroutine debrief_arrays()
    ! Inversion of brief_arrays; using q_minus, TNS, Sigma, and G,
    ! we get back the pdf_mesh. To be used after evolution.
    implicit none

    real(dp), dimension(:,:), allocatable :: q_plus
    integer i, k

    ! Gluon is easy part
    pdf_mesh(0,:) = Gluons

    ! Reduction sequence for unpacking TNS
    allocate( q_plus(nFlavors, Nx) )
    do k=nFlavors, 2, -1
      q_plus(k,:) = ( Sigma(:) - TNS(k,:) ) / k
      Sigma(:) = ( (k-1)*Sigma(:) + TNS(k,:) ) / k
    end do
    q_plus(1,:) = Sigma(:)

    ! With q_plus and q_minus, the rest is easy
    forall(i=1:nFlavors)
        pdf_mesh( i,:) = 0.5_dp*( q_plus(i,:) + q_minus(i,:) )
        pdf_mesh(-i,:) = 0.5_dp*( q_plus(i,:) - q_minus(i,:) )
    end forall

    ! This step is partnered with brief_arrays, so clean up its garbage too.
    deallocate( Gluons, TNS, q_plus, q_minus, Sigma )

    return

  end subroutine debrief_arrays

  ! ============================================================================
  ! ============================================================================

  subroutine evolve_step(newQ2)
    implicit none

    ! I/O
    real(dp), intent(in) :: newQ2

    ! Intermediates
    real(dp) :: dt
    integer :: i, k

    ! Spacing between t=log(Q2)
    dt = log(newQ2/Q2)

    ! Evolve the non-singlet arrays

    !$OMP PARALLEL DO
    do k=2, nFlavors, 1
      call evolve_NS( TNS(k,:), P_qq_0, dt )
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO
    do i=1, nFlavors, 1
      call evolve_NS( q_minus(i,:), P_qq_0, dt )
    end do
    !$OMP END PARALLEL DO

    ! Evolve the singlet arrays
    call evolve_singlet( Sigma, Gluons, P_qq_0, P_qG_0, P_Gq_0, P_GG_0, dt )

    ! Set the new Q2 value in the runningQCD module, and acknowledge it here
    call set_Q2(newQ2)
    Q2 = newQ2

    ! Be done
    return

  end subroutine evolve_step

  ! ============================================================================
  ! ============================================================================

  subroutine evolve_NS( NS_array, P, dt )
    implicit none

    ! I/O
    real(dp), dimension(:), intent(inout) :: NS_array
    real(dp), intent(in) :: dt

    ! Interface for splitting function
    interface
      function P(z, pcode) result(res)
        use kinds, only : dp
        implicit none
        real(dp), intent(in) :: z
        character(len=1), intent(in) :: pcode
        real(dp) :: res
      end function P
    end interface

    ! Intermediates
    real(dp), dimension(:), allocatable :: NS_increments
    integer :: i, icode
    character(len=1), dimension(3), parameter :: pcodes = [ 'd', 'r', 'p' ]

    ! Allocate array of increments, initalize them to zero
    allocate( NS_increments(Nx) )
    NS_increments = 0.0_dp

    do i=1, Nx, 1
      ! Add delta, plus, and regular parts
      do icode = 1, 3, 1
        NS_increments(i) = NS_increments(i) &
            + dglap_convolution( NS_array, P, pcodes(icode), i )
      end do
    end do

    ! Now, add the increments (times factors) to the NS array
    NS_increments = NS_increments * ( get_alphaQCD() / (2.0_dp*pi) ) * dt
    NS_array = NS_array + NS_increments

    ! Collect garbage
    deallocate( NS_increments )
    return

  end subroutine evolve_NS

  ! ============================================================================
  ! ============================================================================

  subroutine evolve_singlet( q_array, G_array, P_qq, P_qG, P_Gq, P_GG, dt )
    implicit none

    ! I/O
    real(dp), dimension(:), intent(inout) :: q_array, G_array
    real(dp), intent(in) :: dt

    ! Interfaces for splitting functions
    interface
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function P_qq(z, pcode) result(res)
        use kinds, only : dp
        implicit none
        real(dp), intent(in) :: z
        character(len=1), intent(in) :: pcode
        real(dp) :: res
      end function P_qq
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function P_qG(z, pcode) result(res)
        use kinds, only : dp
        implicit none
        real(dp), intent(in) :: z
        character(len=1), intent(in) :: pcode
        real(dp) :: res
      end function P_qG
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function P_Gq(z, pcode) result(res)
        use kinds, only : dp
        implicit none
        real(dp), intent(in) :: z
        character(len=1), intent(in) :: pcode
        real(dp) :: res
      end function P_Gq
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function P_GG(z, pcode) result(res)
        use kinds, only : dp
        implicit none
        real(dp), intent(in) :: z
        character(len=1), intent(in) :: pcode
        real(dp) :: res
      end function P_GG
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end interface

    ! Intermediates
    real(dp), dimension(:), allocatable :: q_increments, G_increments
    integer :: i, icode
    character(len=1), dimension(3), parameter :: pcodes = [ 'd', 'r', 'p' ]

    ! Allocate arrays of increments, initalize them to zero
    allocate( q_increments(Nx), G_increments(Nx) )
    q_increments = 0.0_dp
    G_increments = 0.0_dp

    do i=1, Nx, 1
      ! Add delta, plus, and regular parts
      do icode = 1, 3, 1
        ! qq
        q_increments(i) = q_increments(i) &
            + dglap_convolution( q_array, P_qq, pcodes(icode), i )
        ! qG
        q_increments(i) = q_increments(i) &
            + dglap_convolution( G_array, P_qG, pcodes(icode), i )
        ! Gq
        G_increments(i) = G_increments(i) &
            + dglap_convolution( q_array, P_Gq, pcodes(icode), i )
        ! GG
        G_increments(i) = G_increments(i) &
            + dglap_convolution( G_array, P_GG, pcodes(icode), i )
      end do
    end do

    ! Now, add the increments (times factors) to the arrays
    q_increments = q_increments * ( get_alphaQCD() / (2.0_dp*pi) ) * dt
    G_increments = G_increments * ( get_alphaQCD() / (2.0_dp*pi) ) * dt
    q_array = q_array + q_increments
    G_array = G_array + G_increments

    ! Collect garbage
    deallocate( q_increments, G_increments )
    return

  end subroutine evolve_singlet

  ! ============================================================================
  ! ============================================================================

  function dglap_convolution( f, P, pcode, imin ) result(res)
    ! Numerical approximation for the convolution integral:
    !   I(x0) = int(x0,xmax) dy/dy f(y) P(x/y)
    implicit none

    ! I/O
    real(dp), dimension(:), intent(in) :: f
    character(len=1), intent(in) :: pcode
    integer, intent(in) :: imin
    real(dp) :: res

    ! Interface for splitting function
    interface
      function P(z, pcode) result(res)
        use kinds, only : dp
        implicit none
        real(dp), intent(in) :: z
        character(len=1), intent(in) :: pcode
        real(dp) :: res
      end function P
    end interface

    ! Array of integrand evaluations
    real(dp), dimension(:), allocatable :: intd, z
    real(dp) :: logterm

    ! Lower limit of integral
    real(dp) :: xlim

    ! Iterator
    integer i

    ! Lower integration limit
    xlim = x_mesh(imin)

    ! If we're dealing with the delta part of the splitting function,
    ! then be quick about this.
    if( pcode=='d' ) then
      res = P(1.0_dp, pcode) * f(imin)
      return
    end if

    ! For regular and plus parts, allocate the needed memory
    allocate( intd(imin:Nx), z(imin:Nx) )

    ! z = x/y is really x(imin) / x(i)
    z(imin:Nx) = x_mesh(imin) / x_mesh(imin:Nx)

    ! Regular part is easy
    if( pcode=='r' ) then
      do i=imin, Nx, 1
        intd(i) = f(i) * P(z(i), pcode) / x_mesh(i)
      end do
    end if

    ! Plus prescription is a bit tricky
    if( pcode=='p' ) then
      intd(imin) = 0.0_dp
      do i=imin+1, Nx, 1
        intd(i) = ( f(i)*P(z(i), pcode) - z(i)*f(imin)*P(1.0_dp, pcode) ) &
            / ( ( 1.0_dp - z(i) ) * x_mesh(i) )
      end do
      ! There's also a log term
      if( imin < Nx ) then
        logterm = f(imin) * P(1.0_dp, pcode) * log( 1.0_dp - xlim/xmax )
      else
        logterm = 0.0_dp
      end if
    end if

    ! Start result at 0, and add trapezoidal integral approximations
    res = 0.0_dp
    do i=imin, Nx-1, 1
      res = res + int_trapezoid( x_mesh(i), x_mesh(i+1), intd(i), intd(i+1) )
    end do

    ! If plus prescription, there is also a logarithm to add
    if( pcode=='p' ) res = res + logterm

    ! Collect garbage
    deallocate( intd, z )
    return

  end function dglap_convolution

  ! ============================================================================
  ! ============================================================================

  function int_trapezoid( x0, x1, y0, y1 ) result(res)
    ! Integral approxmation, assuming function is m*x+b
    implicit none
    ! I/O
    real(dp), intent(in) :: x0, x1, y0, y1
    real(dp) :: res
    ! This can be done quickly
    res = 0.5_dp * ( y1 + y0 ) * ( x1 - x0 )
    return
  end function int_trapezoid

  ! ============================================================================
  ! ============================================================================

end module dglap
