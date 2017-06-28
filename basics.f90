! ==============================================================================
! basics.f90
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
! This file contains the following modules
!   - kinds
!     :: defines a double precision (dp) kind for reals
!   - pointspace
!     :: routines to generate arrays of spaced points
!   - interpolation
!     :: interpolation routines. Uses fourth-degree polynomial interpolation.

! ==============================================================================
! //////////////////////////////////////////////////////////////////////////////
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! Double precision module follows [ xyzDP ]
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! //////////////////////////////////////////////////////////////////////////////
! ==============================================================================

module kinds
  integer, parameter :: dp = kind(1d0)
end module kinds

! ==============================================================================
! //////////////////////////////////////////////////////////////////////////////
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! Module for creating spaced points follows [ xyzSPACE ]
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! //////////////////////////////////////////////////////////////////////////////
! ==============================================================================

module pointspace
  ! A module with routines for creating arrays of spaced points
  use kinds, only : dp
  implicit none

contains

  ! ============================================================================
  ! ============================================================================

  pure function linspace(startVal, endVal, numPoints) result(val_array)
    ! Creates an array of linearly-spaced values. numpy.linspace clone
    implicit none
    ! I/O
    real(dp), intent(in) :: startVal, endVal
    integer, intent(in) :: numPoints
    real(dp), dimension(numPoints) :: val_array
    ! Intermediate variables
    integer :: i
    real(dp) :: stepSize
    ! Step size
    stepSize = ( endVal - startVal ) / ( numPoints - 1 )
    ! Array creation
    do i=1, numPoints, 1
       val_array(i) = stepSize*(i-1) + startVal
    end do
    ! Be done
    return
  end function linspace

  ! ============================================================================
  ! ============================================================================

  function logspace(startVal, endVal, numPoints) result(val_array)
    ! Creates an array of logarithmically-spaced values.
    implicit none
    ! I/O
    real(dp), intent(in) :: startVal, endVal
    integer, intent(in) :: numPoints
    real(dp), dimension(numPoints) :: val_array
    ! Intermediate variables
    integer :: i
    real(dp) :: stepSize
    ! Step size
    stepSize = ( endVal / startVal )**( 1.0_dp/(numPoints-1) )
    ! Create the logarithmically spaced points
    do i=1, numPoints, 1
       val_array(i) = startVal * stepSize**(i-1)
    end do
  end function logspace

  ! ============================================================================
  ! ============================================================================

  function revlogspace(startVal, endVal, numPoints) result(val_array)
    ! Creates a space of reverse-logarithmically spaced values,
    ! which become ever-closer to the endpoint
    implicit none
    ! I/O
    real(dp), intent(in) :: startVal, endVal
    integer, intent(in) :: numPoints
    real(dp), dimension(numPoints) :: val_array
    ! Intermediate variables
    integer :: i
    real(dp) :: stepSize
    ! Step size
    stepSize = ( startVal / endVal)**( 1.0_dp/(numPoints-1) )
    ! Create reverse-logarithmically spaced points
    do i=1, numPoints, 1
      val_array(i) = startVal + endVal*( 1.0_dp - stepSize**(i-1) )
    end do
    return
  end function revlogspace

  ! ============================================================================
  ! ============================================================================

end module pointspace

! ==============================================================================
! //////////////////////////////////////////////////////////////////////////////
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! Module for interpolating follows [ xyzINTERPOL ]
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! //////////////////////////////////////////////////////////////////////////////
! ==============================================================================

module interpolation
  use kinds, only : dp
  implicit none
  private
  integer, parameter :: N = 4 ! order of interpolation
  public :: interpolate_1D, interpolate_2D, interpolate_3D, interpolate_4D

contains

  pure function interpolate_4D(xx, yy, zz, ww, ff, x, y, z, w) result(f)
    ! Inputs.
    !   xx is an ordered array of x values
    !   yy is an ordered array of y values
    !   zz is an ordered array of z values
    !   ww is an ordered array of w values
    !   ff is a grid of f values. f(i,j,k,l) = f(x(i),y(j),z(k),w(l))
    ! This code assumes dimension of is size(x) by size(y) by size(z) by size(w)
    implicit none

    ! I/O
    real(dp), intent(in) :: xx(:), yy(:), zz(:), ww(:), ff(:,:,:,:), x, y, z, w
    real(dp) :: f
    ! Intermediates
    integer Nx, Ny, Nz, Nw, ixLoc, iyLoc, izLoc, iwLoc
    integer ixStart, iyStart, izStart, iwStart, ixEnd, iyEnd, izEnd, iwEnd
    real(dp), dimension(N) :: fyzw, fxzw, fxyw, fxyz
    real(dp) :: fx, fy, fz, fw
    real(dp) :: sigma_x, sigma_y, sigma_z, sigma_w
    integer i

    ! Sizes, indices
    Nx = size(xx)
    Ny = size(yy)
    Nz = size(zz)
    Nw = size(ww)
    ixLoc = locate(xx, x)
    iyLoc = locate(yy, y)
    izLoc = locate(zz, z)
    iwLoc = locate(ww, w)
    ixStart = min( max(ixLoc-(N-1)/2, 1), Nx+1-N)
    ixEnd = ixStart + N - 1
    iyStart = min( max(iyLoc-(N-1)/2, 1), Ny+1-N)
    iyEnd = iyStart + N - 1
    izStart = min( max(izLoc-(N-1)/2, 1), Nz+1-N)
    izEnd = izStart + N - 1
    iwStart = min( max(iwLoc-(N-1)/2, 1), Nw+1-N)
    iwEnd = iwStart + N - 1

    ! Partially interpolated arrays of values
    forall (i=1:N)
       fyzw(i) = interpolate_3D(yy, zz, ww, ff(ixStart+i-1,:,:,:), y, z, w)
       fxzw(i) = interpolate_3D(xx, zz, ww, ff(:,iyStart+i-1,:,:), x, z, w)
       fxyw(i) = interpolate_3D(xx, yy, ww, ff(:,:,izStart+i-1,:), x, y, w)
       fxyz(i) = interpolate_3D(xx, yy, zz, ff(:,:,:,iwStart+i-1), x, y, z)
    end forall

    ! Interpolate the partially interpolated arrays, getting their errors
    call interpolate_1D_error( xx(ixStart:ixEnd), fyzw, x, fx, sigma_x )
    call interpolate_1D_error( yy(iyStart:iyEnd), fxzw, y, fy, sigma_y )
    call interpolate_1D_error( zz(izStart:izEnd), fxyw, z, fz, sigma_z )
    call interpolate_1D_error( ww(iwStart:iwEnd), fxyz, w, fw, sigma_w )

    ! Weigh by errors, or just use one particular order if the error is zero
    if (sigma_x .eq. 0.0_dp) then
       f = fx
    else if (sigma_y .eq. 0.0_dp) then
       f = fy
    else if (sigma_z .eq. 0.0_dp) then
       f = fz
    else if (sigma_w .eq. 0.0_dp) then
       f = fw
    else
       f = ( fx/sigma_x**2 + fy/sigma_y**2 + fz/sigma_z**2 + fw/sigma_w**2 ) &
            / ( sigma_x**(-2) + sigma_y**(-2) + sigma_z**(-2) + sigma_w**(-2) )
    end if

    return
  end function interpolate_4D

  ! ============================================================================
  ! ============================================================================

  pure function interpolate_3D(xx, yy, zz, ff, x, y, z) result(f)
    ! Inputs.
    !   xx is an ordered array of x values
    !   yy is an ordered array of y values
    !   zz is an ordered array of z values
    !   ff is a grid of f values. f(i,j,k) = f(x(i),y(j),z(k))
    ! This code assumes dimension of is size(x) by size(y) by size(z).
    implicit none

    ! I/O
    real(dp), intent(in) :: xx(:), yy(:), zz(:), ff(:,:,:), x, y, z
    real(dp) :: f
    ! Intermediates
    integer Nx, Ny, Nz, ixLoc, iyLoc, izLoc
    integer ixStart, iyStart, izStart, ixEnd, iyEnd, izEnd
    real(dp), dimension(N) :: fxy, fyz, fxz
    real(dp) :: fx, fy, fz
    real(dp) :: sigma_x, sigma_y, sigma_z
    integer i

    ! Sizes, indices
    Nx = size(xx)
    Ny = size(yy)
    Nz = size(zz)
    ixLoc = locate(xx, x)
    iyLoc = locate(yy, y)
    izLoc = locate(zz, z)
    ixStart = min( max(ixLoc-(N-1)/2, 1), Nx+1-N)
    ixEnd = ixStart + N -1
    iyStart = min( max(iyLoc-(N-1)/2, 1), Ny+1-N)
    iyEnd = iyStart + N -1
    izStart = min( max(izLoc-(N-1)/2, 1), Nz+1-N)
    izEnd = izStart + N -1

    ! Partially interpolated arrays of values
    forall (i=1:N)
       fxy(i) = interpolate_2D(xx, yy, ff(:,:,izStart+i-1), x, y)
       fyz(i) = interpolate_2D(yy, zz, ff(ixStart+i-1,:,:), y, z)
       fxz(i) = interpolate_2D(xx, zz, ff(:,iyStart+i-1,:), x, z)
    end forall

    ! Interpolate the partially interpolated arrays, getting their errors
    call interpolate_1D_error( xx(ixStart:ixEnd), fyz, x, fx, sigma_x )
    call interpolate_1D_error( yy(iyStart:iyEnd), fxz, y, fy, sigma_y )
    call interpolate_1D_error( zz(izStart:izEnd), fxy, z, fz, sigma_z )

    ! Weigh by errors, or just use one particular order if the error is zero
    if (sigma_x .eq. 0.0_dp) then
       f = fx
    else if (sigma_y .eq. 0.0_dp) then
       f = fy
    else if (sigma_z .eq. 0.0_dp) then
       f = fz
    else
       f = ( fx / sigma_x**2 + fy / sigma_y**2 + fz / sigma_z**2 ) &
            / ( sigma_x**(-2) + sigma_y**(-2) + sigma_z**(-2) )
    end if

    return
  end function interpolate_3D

  ! ============================================================================
  ! ============================================================================

  pure function interpolate_2D(xx, yy, ff, x, y) result(f)
    ! Inputs.
    !   xx is an ordered array of x values
    !   yy is an ordered array of y values
    !   ff is a grid of f values. f(i,j) = f(x(i),y(j))
    ! This code assumes dimension of is size(x) by size(y).
    implicit none

    ! I/O
    real(dp), intent(in) :: xx(:), yy(:), ff(:,:), x, y
    real(dp) :: f
    ! Intermediates
    integer Nx, Ny, ixLoc, iyLoc, ixStart, iyStart, ixEnd, iyEnd
    real(dp), dimension(N) :: fx, fy
    real(dp) :: f_xy, f_yx, sigma_xy, sigma_yx
    integer i

    ! Sizes, indices
    Nx = size(xx)
    Ny = size(yy)
    ixLoc = locate(xx, x)
    iyLoc = locate(yy, y)
    ixStart = min( max(ixLoc-(N-1)/2, 1), Nx+1-N)
    ixEnd = ixStart + N -1
    iyStart = min( max(iyLoc-(N-1)/2, 1), Ny+1-N)
    iyEnd = iyStart + N -1

    ! Partially interpolated arrays of values
    forall (i=1:N)
       fx(i) = interpolate_1D(yy, ff(ixStart+i-1,:), y)
       fy(i) = interpolate_1D(xx, ff(:,iyStart+i-1), x)
    end forall

    ! Interpolate the partially interpolated arrays, getting their errors
    call interpolate_1D_error( xx(ixStart:ixEnd), fx, x, f_xy, sigma_xy )
    call interpolate_1D_error( yy(iyStart:iyEnd), fy, y, f_yx, sigma_yx )

    ! Weigh by errors, or just use one particular order if the error is zero
    if (sigma_xy .eq. 0.0_dp) then
       f = f_xy
    else if (sigma_yx .eq. 0.0_dp) then
       f = f_yx
    else
       f = ( f_xy / sigma_xy**2 + f_yx / sigma_yx**2 ) &
            / ( sigma_xy**(-2) + sigma_yx**(-2) )
    end if

    return
  end function interpolate_2D

  ! ============================================================================
  ! ============================================================================

  pure function interpolate_1D(xx, ff, x) result(f)
    implicit none
    ! I/O
    real(dp), intent(in) :: x, xx(:), ff(:)
    real(dp) :: f
    ! Intermediate dummy variable
    real(dp) :: ferr
    call interpolate_1D_error(xx, ff, x, f, ferr)
    return
  end function interpolate_1D

  ! ============================================================================
  ! ============================================================================

  pure subroutine interpolate_1D_error(xx, ff, x, f, ferr)
    implicit none
    ! I/O
    real(dp), intent(in) :: x, xx(:), ff(:)
    real(dp), intent(out) :: f, ferr
    ! Intermediates
    integer iLoc, iStart, iEnd, Nx

    Nx = size(xx)
    iLoc = locate(xx, x)
    iStart = min( max(iLoc-(N-1)/2, 1), Nx+1-N )
    iEnd = iStart + N - 1

    call polint( xx(iStart:iEnd), ff(iStart:iEnd), x, f, ferr)
    return
  end subroutine interpolate_1D_error

  ! ============================================================================
  ! ============================================================================

  pure subroutine polint(xx, yy, x, y, dy)
    ! Polynomial interpolation based on Newton's divided differences algorithm
    implicit none
    ! I/O
    real(dp), dimension(:), intent(in) :: xx, yy
    real(dp), intent(in) :: x
    real(dp), intent(out) :: y, dy
    ! Intermediates
    integer n, i, j
    real(dp), dimension(size(xx)) :: d
    ! Get size of array
    n = size(xx)
    ! Fill the divided difference table
    d = yy
    do i = 1, n, 1
       do j = n, i+1, -1
          d(j) = ( d(j)-d(j-1) ) / ( xx(j)-xx(j-i) )
       end do
    end do
    ! Nested multiplication for the polynomials
    y = d(n)
    do i = n-1, 1, -1
       y = d(i) + (x-xx(i))*y
    end do
    ! Error estimate
    dy = d(minloc(abs(x-xx),1))
    return
  end subroutine polint

  ! ============================================================================
  ! ============================================================================

  pure function locate(xx,x) result(iLoc)
    ! Finds an index iLoc such that xx(iLoc) < x < xx(iLoc+1)
    ! Assumes the sequence is ordered in ascending order.
    implicit none
    ! I/O
    real(dp), intent(in) :: xx(:), x
    integer :: iLoc
    ! Intermediates
    integer :: N, iStart, iMid, iEnd
    ! First, determine the array size
    N = size(xx)
    ! If x is too small or too large, then return an endpoint
    if( x <= xx(1) ) then
       iLoc = 1
       return
    else if( x >= xx(n) ) then
       iLoc = n - 1
       return
    end if
    ! Define initial endpoints for our search
    iStart = 0
    iEnd = N+1
    ! Search until the start and midpoints come together
    do while( iEnd-iStart > 1 )
       ! Find the midpoint of current start and end
       iMid = (iStart + iEnd) / 2
       ! If x is above the midpoint, move the start there. Else move the end.
       if( x >= xx(iMid) ) then
          iStart = iMid
       else
          iEnd = iMid
       end if
    end do
    ! The desired location is now the narrowed start point
    iLoc = iStart
    return
  end function locate

  ! ============================================================================
  ! ============================================================================

end module interpolation
