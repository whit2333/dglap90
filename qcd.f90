! ==============================================================================
! qcd.f90
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
!   - QCD
!   :: Mathematical and physical constants.
!   - RunningQCD
!   :: QCD parameters that depend on Q**2
!
! Everything in these modules assumes the MS-bar scheme, and zero-mass variable
! flavor number scheme.

! ==============================================================================
! //////////////////////////////////////////////////////////////////////////////
! ==============================================================================

module QCD
  ! This module contains constants relevant to QCD, assuming the MS-bar scheme.
  use kinds, only : dp
  implicit none
  ! Casimir invariants
  real(dp), public, parameter :: CF = 4.0_dp/3.0_dp
  real(dp), public, parameter :: CA = 3.0_dp
  real(dp), public, parameter :: TF = 1.0_dp/2.0_dp
  ! Mathematical constants
  real(dp), public, parameter :: pi = acos(-1.0_dp)
  real(dp), public, parameter :: zeta2 = pi**2 / 6.0_dp
  real(dp), public, parameter :: zeta3 = 1.2020569031_dp
  real(dp), public, parameter :: zeta4 = pi**4 / 90.0_dp
  ! Current quark masses for effective quark number thresholds
  real(dp), private, parameter :: cMass = 1.29_dp
  real(dp), private, parameter :: bMass = 4.5_dp
  real(dp), private, parameter :: tMass = 172.44_dp
  real(dp), public, parameter :: cMass2 = cMass**2
  real(dp), public, parameter :: bMass2 = bMass**2
  real(dp), public, parameter :: tMass2 = tMass**2
  ! eps is a tiny quantity with which to avoid singularities
  real(dp), public, parameter :: eps = 0.000000001_dp
  ! Parton charges (only their squares are really needed)
  real(dp), private, parameter :: eMinus = -1.0_dp/3.0_dp
  real(dp), private, parameter :: ePlus = 2.0_dp/3.0_dp
  real(dp), private, parameter, dimension(-6:6) :: e = &
       & (/ -ePlus, -eMinus, -ePlus, -eMinus, -eMinus, -ePlus, 0.0_dp, &
       & ePlus, eMinus, eMinus, ePlus, eMinus, ePlus /)
  real(dp), public, parameter, dimension(-6:6) :: e2 = e**2
  ! Lambda values used for quark threshold matching
  ! From PDG 2016:
  !  Lambda(nFl = 3) :: 332 +/- 17 MeV
  !  Lambda(nFl = 4) :: 292 +/- 16 MeV
  !  Lambda(nFl = 5) :: 210 +/- 14 MeV
  !  Lambda(nFl = 6) ::  89 +/-  6 MeV
  real(dp), public, parameter, dimension(3:6) :: &
       & LambdaFl = [ 0.332_dp, 0.292_dp, 0.210_dp, 0.089_dp ]
end module QCD

! ==============================================================================
! //////////////////////////////////////////////////////////////////////////////
! ==============================================================================

module runningQCD
  ! A module for dealing with the running coupling of QCD.
  ! Assumes the MS-bar scheme.
  ! Q2 is set module-wide
  use kinds, only : dp
  use QCD
  implicit none
  private
  real(dp) :: Q2, tau, nFl, alphaQCD, Lambda, Lambda2
  integer :: nFlavors
  integer :: nOrder = 1 ! LO by default
  public :: &
      ! Mutators
      set_Q2, set_QCD_order, &
      ! Accessors
      get_Q2, get_alphaQCD, get_nFlavors, get_Lambda2
contains

  ! ============================================================================
  ! ****************************************************************************
  ! *** Public methods ***
  ! ****************************************************************************
  ! ============================================================================

  subroutine set_Q2(givenQ2)
    ! Mutator for Q2
    implicit none
    real(dp), intent(in) :: givenQ2
    Q2 = givenQ2
    call set_number_flavors()
    call set_lambda_QCD()
    tau = log(Q2/Lambda2)
    call set_alpha_QCD()
    return
  end subroutine set_Q2

  ! ============================================================================
  ! ============================================================================

  function get_Q2() result(res)
    ! Accessor for Q2
    implicit none
    ! I/O
    real(dp) :: res
    res = Q2
    return
  end function get_Q2

  ! ============================================================================
  ! ============================================================================

  function get_alphaQCD(Q2val) result(res)
    ! Accessor for alphaQCD
    implicit none
    ! I/O
    real(dp), intent(in), optional :: Q2val
    real(dp) :: res
    ! If the optional argument was passed, set Q2 using it.
    if( present(Q2val) ) call set_Q2(Q2val)
    ! Return the current internal value for alphaQCD
    res = alphaQCD
    return
  end function get_alphaQCD

  ! ============================================================================
  ! ============================================================================

  function get_nFlavors() result(res)
    ! Accessor for number of effective flavors
    implicit none
    ! I/O
    integer :: res
    res = nFlavors
    return
  end function get_nFlavors

  ! ============================================================================
  ! ============================================================================

  function get_Lambda2() result(res)
    ! Accessor for number of effective flavors
    implicit none
    ! I/O
    real(dp) :: res
    res = Lambda2
    return
  end function get_Lambda2

  ! ============================================================================
  ! ============================================================================

  subroutine set_QCD_order(given_order)
    ! Set the order: 1 for LO, 2 for NLO, etc.
    implicit none
    ! I/O
    integer, intent(in) :: given_order
    nOrder = given_order
    return
  end subroutine set_QCD_order

  ! ============================================================================
  ! ****************************************************************************
  ! *** Private (internal) methods ***
  ! ****************************************************************************
  ! ============================================================================

  subroutine set_number_flavors()
    ! Set the effective number of flavors with a zero-mass scheme
    implicit none
    if(Q2 .gt. tMass2) then
       nFlavors = 6
    else if(Q2 .gt. bMass2) then
       nFlavors = 5
    else if(Q2 .gt. cMass2) then
       nFlavors = 4
    else
       nFlavors = 3
    end if
    ! A real version of nFlavors to be plugged into floating-point arithmetic
    nFl = nFlavors
    return
  end subroutine set_number_flavors

  ! ============================================================================
  ! ============================================================================

  subroutine set_lambda_QCD()
    ! Set the effective value of Lambda_QCD
    implicit none
    Lambda = LambdaFl(nFlavors)
    Lambda2 = Lambda**2
    return
  end subroutine set_lambda_QCD

  ! ============================================================================
  ! ============================================================================

  subroutine set_alpha_QCD()
    ! Set the value of alphaQCD
    ! Alpha is determined in the MS-bar scheme according to
    !  [1] P. A. Baikov, et al., PRL 118, 082002 (2017)
    implicit none
    real(dp), dimension(0:4) :: alpha, beta
    ! (constant) matrices used to build beta coefficients
    real(dp), dimension(0:4,0:4) :: betaMatrix = reshape ( (/ &
        2.75_dp,    -0.166667_dp, 0.0_dp,       0.0_dp,       0.0_dp, &
        6.375_dp,   -0.791667_dp, 0.0_dp,       0.0_dp,       0.0_dp, &
        22.3203_dp, -4.36892_dp,  0.0940394_dp, 0.0_dp,       0.0_dp, &
        114.23_dp,  -27.1339_dp,  1.58238_dp,   0.0058567_dp, 0.0_dp, &
        524.56_dp,  -181.8_dp,    17.16_dp,    -0.2258_dp,   -0.0017993_dp &
        /), (/5,5/), order=(/2,1/) )
    ! Vector matrix with powers of the number of flavors
    real(dp), dimension(0:4) :: nFlMatrix
    integer i

    ! Fill up flavor array
    forall(i=0:4)
        nFlMatrix(i) = nFlavors**i
    end forall

    ! Fill up the beta values
    ! The convention used in [1] is missing factors of pi, so include them here.
    forall(i=0:4)
        beta(i) = sum( betaMatrix(i,:)*nFlMatrix(:) ) / pi**(i+1)
    end forall

    ! Get the alpha contributions from various orders
    ! NB these are only approximations anyway (except at LO).
    alpha(0) = 1.0_dp
    alpha(1) = -( beta(1)*log(tau) ) / ( beta(0)**2 * tau )
    alpha(2) = ( beta(1)**2*( log(tau)**2 - log(tau) - 1.0_dp ) &
        + beta(0)*beta(2) ) / ( beta(0)**4 * tau )
    alpha(3) = - ( beta(1)**3*( log(tau)**3 - 2.5_dp*log(tau)**2 &
        - 2.0_dp*log(tau) + 0.5_dp ) + 3.0_dp*beta(0)*beta(1)*beta(2)*log(tau) &
        - 0.5_dp*beta(0)**2*beta(3) ) / ( beta(0)**6 * tau**3 )
    alpha(4) = 0.0_dp ! TODO

    ! Sum up the contributions
    alphaQCD = sum( alpha(0:nOrder-1) )  / ( beta(0)*tau )

    return
  end subroutine set_alpha_QCD

  ! ============================================================================
  ! ============================================================================

end module RunningQCD
