
!!----------------------------------------------------------------------
!!This source file is free software: you can redistribute it and/or modify
!!it under the terms of the GNU General Public License as published by
!!the Free Software Foundation, either version 3 of the License, or
!!(at your option) any later version.
!!
!!This program is distributed in the hope that it will be useful,
!!but WITHOUT ANY WARRANTY; without even the implied warranty of
!!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!GNU General Public License for more details.
!!
!!You should have received a copy of the GNU General Public License
!!along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!----------------------------------------------------------------------


! SIMULATION STATE DATA
!
! This module contains a series of flags indicating the state
! of the simulation.  This state can then be queried.
!
! -------------------------------------------------------------

module simstate
  implicit none
  private
  ! Set to true to indicate that all charge interactions
  ! are to be ignored.
  logical, private :: nocharge=.false.

  ! Indicator of whether the program is simulating under
  ! bulk (periodic boundary conditions) or not
  logical, private :: in_bulk_simulation = .false.
  logical, private :: in_widom_trial = .false.

  ! Accessors and modifiers

  public :: start_bulk, start_cell, is_bulk, use_electrostatic, do_electrostatic
  public :: start_widom_trials, end_widom_trials, in_widom

  contains

  subroutine use_electrostatic(doit)
    implicit none
    logical, intent(in) :: doit
    nocharge = doit
  end subroutine use_electrostatic

  function do_electrostatic ()
    implicit none
    logical :: do_electrostatic
    do_electrostatic=.not.nocharge
  end function do_electrostatic

  subroutine start_bulk
    implicit none
    in_bulk_simulation = .true.
  end subroutine start_bulk

  subroutine start_cell
    implicit none
    in_bulk_simulation = .false.
  end subroutine start_cell

  subroutine start_widom_trials
    implicit none
    in_widom_trial = .true.
  end subroutine start_widom_trials

  subroutine end_widom_trials
    implicit none
    in_widom_trial = .false.
  end subroutine end_widom_trials

  logical function is_bulk()
    implicit none
    is_bulk = in_bulk_simulation
  end function is_bulk

  logical function in_widom()
    implicit none
    in_widom = in_widom_trial
  end function in_widom

end module simstate
