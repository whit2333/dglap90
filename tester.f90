! ==============================================================================
! tester.f90
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
! This file contains a test program that evolves a CJ15 LO PDF from
! Q**2 = 5 GeV**2 to Q**2 = 1000 GeV**2.

program test
  use kinds, only : dp
  use dglap
  implicit none

  ! Read the low-Q2 file
  call dglap_read_csv('pdf_grids/cj15_Q2_5_GeV2.csv', 5.0_dp)

  ! Evolve it to high Q2
  call dglap_evolve( 1000.0_dp )

  ! Save the result to a new CSV file
  call dglap_write_csv('pdf_grids/evolved_Q2_1000_GeV2.csv')

  ! Free up memory and quit
  call dglap_cleanup()
  return
end program test
