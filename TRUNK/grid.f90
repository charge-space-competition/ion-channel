
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


! -------------------------------------------------------------
! INITIAL CONFIGURATION HELPER
! ----------------------------
!
! The V17 version of the generate of the initial configuration 
! simply attempted to add particles to random positions until an
! acceptable position was found.  This type of generator is
! good for gas like configurations, but the chance of adding
! particles approachs 0 well before maximum packing is reached.
! The grid configuration helpers work by constructing a grid
! of points separated by the minimum approach distance of the
! largest particle with itself that fits within a particular
! shape (here either a cube or a cylinder).  The user then
! requests locations that are randomly selected from the grid
! until all grid-points are exhausted.
!
! ------------------------------------------------------------
! Location of generated positions (from gridnext)
!   cubegrid
!     x, y and z are generated starting from 1..width * length
!       where width*length =< cubelength - 2*length
!
!   tubegrid
!     x and y are generated within circle of (radius - length)
!     z is generated in interval +/- 1/2 * (tubelength - 2*length)
!       with the closest points to 0 being +/- 1/2 length and
!       half the points lying below and above 0.
!
! ------------------------------------------------------------
! ACCESS SUMMARY
!
! TYPE: cubegrid, tubegrid:
!        Cubic grids that fit within the particular shapes
!
! access via:
!
! METH: gridinit(grid,length [,radius], gridspacing) 
!        Initialise the cube or tube grids with the given lengths
!        (and for the tube, the tube radius)
! 
! METH: gridnext(grid, xout, yout, zout)
!        Get the next location from the grid. The grid will 'stop'
!        if there is no more grid locations
! 
! METH: int gridsize(grid)
!        Get the remaining number of points
!
! METH: int gridkill(grid)
!        Release any memory for this object
!
! METH: int gridisok(grid)
!        Can this grid be used?
!

module grid
use const
implicit none
private

type cubegrid
  ! inter-grid spacing
  double precision :: length
  ! maximum grid widths in a direction
  integer :: width
  ! Grid point selector array
  integer, dimension(:), allocatable :: selector
  ! The number of grid points still accessible.
  integer :: selector_size
end type cubegrid

type tubegrid
  ! inter-grid spacing
  double precision :: length
  ! maximum grid widths in x, y and z
  integer :: xmax, ymax, zmax
  ! The number of xelements at a particular y
  integer, dimension(:), allocatable :: width
  ! The number of elements in an xy slice
  integer :: gridslice
  ! Grid point selector array
  integer, dimension(:), allocatable :: selector
  ! The number of grid points still accessible.
  integer :: selector_size
end type tubegrid

public :: tubegrid, cubegrid

interface gridinit
  module procedure cubegridinit
  module procedure tubegridinit
end interface gridinit

interface gridnext
  module procedure cubegridnext
  module procedure tubegridnext
end interface gridnext

interface gridsize
  module procedure cubegridsize
  module procedure tubegridsize
end interface gridsize

interface gridisok
  module procedure cubegridisok
  module procedure tubegridisok
end interface gridisok

interface gridkill
  module procedure cubegridkill
  module procedure tubegridkill
end interface gridkill


public :: gridinit, gridnext, gridsize, gridisok, gridkill

contains


  ! METH: gridinit(grid,length,radius, gridspacing) 
  !        Initialise the cube grids with the given lengths.
  subroutine cubegridinit(self,cubelength,length)
    implicit none
    type (cubegrid), intent(inout) :: self  
    double precision, intent(in) :: cubelength, length
    integer :: xidx, stat
    if (dbc) then
      if (gridisok(self)) then
        stop "Call to gridinit on already initialised object"
      endif
    endif

    ! grid spacing
    self%length = length
    ! width : maximum x, y or z grid coordinate
    self%width = floor((cubelength - 2 * self%length)/self%length)
    ! make width even
    self%width = ibclr(self%width,0)

    ! total number of grid points
    self%selector_size = self%width**3

    ! set up selector grid
    allocate(self%selector(self%selector_size), STAT=stat)
    if (stat.ne.0) stop "Memory allocation of SELF%SELECTOR failed"
    do xidx=1,self%selector_size
      self%selector(xidx)=xidx
    enddo
    ! Randomise the selector entries
    call shuffle(self%selector, self%selector_size)
    if (dbc) then
      if (.not.gridisok(self)) then
        stop "Call to gridinit failed to initialise object"
      endif
    endif

  end subroutine cubegridinit

  ! gridnext(grid, xout, yout, zout)
  !        Get the next location from the grid. The grid will 'stop'
  !        if there is no more grid locations
  subroutine cubegridnext(self, xout, yout, zout)
    implicit none
    type (cubegrid), intent(inout) :: self
    double precision, intent(out) :: xout, yout, zout
    integer :: xindex, yindex, zindex
    integer :: widsq
    if (dbc) then
      if (.not.gridisok(self)) then
        stop "Call to gridnext with uninitialised object"
      endif
    endif
    widsq = self%width**2
    ! decrement grid index 
    self%selector_size = self%selector_size - 1
    if (self%selector_size.lt.1) then
      write(unit=fidlog,fmt=*)"Stopping after ", self%width**3, &
           " attempts exhausted all gridpoints"
      stop "All grid points used"
    endif

    ! Get grid index from selector matrix
    xindex = self%selector(self%selector_size)

    ! Generate zindex and yindex from xindex
    zindex = 1
    do while (xindex.gt.widsq)
      xindex = xindex - widsq
      zindex = zindex + 1
    enddo
    yindex = 1
    do while (xindex.gt.self%width)
      xindex = xindex - self%width
      yindex = yindex + 1
    enddo
    
    ! create particle at current grid point
    xout = self%length * xindex
    yout = self%length * yindex
    zout = self%length * zindex
  end subroutine cubegridnext

  ! int gridsize(grid)
  !        Get the remaining number of points
  function cubegridsize(self)
    implicit none
    type (cubegrid), intent(inout) :: self
    integer :: cubegridsize

    if (dbc) then
      if (.not.gridisok(self)) then
        stop "Call to gridnext with uninitialised object"
      endif
    endif

    cubegridsize=self%selector_size
  end function cubegridsize

  ! int gridisok(grid)
  !        Test if the grid object can be used
  function cubegridisok(self)
    implicit none
    type (cubegrid), intent(inout) :: self
    logical :: cubegridisok

    cubegridisok = allocated(self%selector).and.self%selector_size.gt.0
  end function cubegridisok

  ! gridkill(grid)
  !        Release any allocated memory
  subroutine cubegridkill(self)
    implicit none
    type (cubegrid), intent(inout) :: self
    integer :: stat
    self%selector_size = 0
    stat = 0
    if (allocated(self%selector)) deallocate(self%selector, STAT=stat)
    if (stat.ne.0) stop "Memory deallocation of failed"
  end subroutine cubegridkill


  ! METH: gridinit(grid,length,radius, gridspacing) 
  !        Initialise the tube grids with the given lengths.
  subroutine tubegridinit(self,tubelength,radius,length)
    implicit none
    type (tubegrid), intent(inout) :: self  
    double precision, intent(in) :: radius, tubelength, length
    integer :: xidx, yidx, error, tidx, idx, widx
    integer, dimension(:), allocatable :: top
    integer :: stat

    if (dbc) then
      if (gridisok(self)) then
        stop "Call to gridinit on already initialised object"
      endif
    endif

    ! grid spacing
    self%length = length
    ! radius : the radius of the tube (rl5)
    ! xmax : maximum x grid coordinate
    self%xmax = int(floor((radius - self%length)/ self%length))
    ! ymax : initially maximum xmax = ymax coordinate (xmax/sqrt(2))
    self%ymax = int(radius / (sqrt(2.D0) * self%length))

    ! tubelength : the length of the tube 2 * (zl4 - zl2)
    ! zmax : max z coord (tublength/length)
    self%zmax = floor((tubelength - 2 * self%length)/self%length)
    ! make zmax even
    self%zmax = ibclr(self%zmax,0)
      
    ! width : number of cubes along the xidx axis at a point on the yidx axis
    allocate(self%width(self%xmax), STAT=stat)
    if (stat.ne.0) stop "Memory allocation of SELF%WIDTH failed"
    self%width = 0
 
   ! xidx, yidx, z : grid coordinates
    xidx = self%xmax
    yidx = 1  ! Fortan counting
    
    ! top : number of cubes along the xidx axis at a point on the yidx axis
    !       above ymax 
    allocate(top(self%xmax), STAT=stat)
    if (stat.ne.0) stop "Memory allocation of SELF%XMAX failed"
    
    ! error: value that controls movement in the xidx direction. When
    !        error > 0 we increment xidx.  This only works up to ymax
    error = -xidx
    
    widx = 0
    tidx = 0
    do while (xidx.ge.yidx)
      widx = widx + 1  ! Fortran counting
      self%width(widx) =  xidx
      error = error + yidx
      yidx = yidx + 1
      error = error + yidx
      if (error.ge.0) then
        error = error - xidx
        xidx = xidx - 1
        error = error - xidx
        tidx = tidx + 1 ! Fortran counting
        top(tidx) = yidx - 1
      endif
    enddo
    ! Last entry in top will be equal to last in width for
    ! even cubes
    if (self%width(widx).eq.top(tidx)) then
      ! drop last entry in top
      tidx = tidx - 1
    endif
    if (tidx.ge.1) then
      do idx=tidx,1,-1
        widx = widx + 1  ! Fortran counting
        self%width(widx) = top(idx)
      enddo
    endif
    deallocate(top, STAT=stat)
    if (stat.ne.0) stop "Memory deallocation of failed"

    if (widx.gt.self%xmax) then
      stop "Array over run in cubicgrid"
    endif

    ! reset ymax to be the maximum idx in width
    self%ymax = widx

    ! ymax is now maximum index of width matrix

    ! number of gridpoints on the xy plane
    self%gridslice = 4*sum(self%width(1:self%ymax))
    ! total number of grid points
    self%selector_size = self%gridslice*self%zmax

    ! set up selector grid
    allocate(self%selector(self%selector_size), STAT=stat)
    if (stat.ne.0) stop "Memory allocation of SELF%SELECTOR failed"
    do xidx=1,self%selector_size
      self%selector(xidx)=xidx
    enddo
    ! Randomise the selector entries
    call shuffle(self%selector, self%selector_size)

    if (debug) then
      write(unit=fidlog,fmt=*)"TUBEGRID parameters:"
      write(unit=fidlog,fmt=*)"Tube Radius :",radius
      write(unit=fidlog,fmt=*)"Grid spacing:",self%length
      write(unit=fidlog,fmt=*)"XY-max      : +/-",self%xmax
      write(unit=fidlog,fmt=*)"Z-max       :    ",self%zmax


      write(unit=fidlog,fmt=*)"Grid points :    ",sum(self%width(1:self%ymax)) * 4 * self%zmax
      do idx=1,widx
        write(unit=fidlog,fmt=*)("  *",tidx=1,self%width(idx))
      enddo
      if (widx.ne.self%xmax) then
        write(unit=fidlog,fmt=*)"Odd Circle in inittubegrid:"
        write(unit=fidlog,fmt=*)" predicted max idx = ", self%xmax
        write(unit=fidlog,fmt=*)" found max idx     = ", widx
        if (widx.gt.self%xmax) then
          stop "Array over run in inittubegrid"
        endif
      endif
    endif

    if (dbc) then
      if (.not.gridisok(self)) then
        stop "Call to gridinit failed to initialise object"
      endif
    endif

  end subroutine tubegridinit

  ! gridnext(grid, xout, yout, zout)
  !        Get the next location from the grid. The grid will 'stop'
  !        if there is no more grid locations
  subroutine tubegridnext(self, xout, yout, zout)
    implicit none
    type (tubegrid), intent(inout) :: self
    double precision, intent(out) :: xout, yout, zout
    integer :: xindex, yindex, zindex, yincrement

    integer :: dbg_max

    if (dbc) then
      if (.not.gridisok(self)) then
        stop "Call to gridnext with uninitialised object"
      endif
    endif

    ! decrement grid index 
    self%selector_size = self%selector_size - 1
    if (self%selector_size.lt.1) then
      write(unit=fidlog,fmt=*)"Stopping after ", self%gridslice*self%zmax, &
           " attempts exhausted all gridpoints"
      stop "All grid points used"
    endif

    ! Get grid index from selector matrix
    xindex = self%selector(self%selector_size)

    ! Generate zindex from xindex
    zindex = 1
    do while (xindex.gt.self%gridslice)
      xindex = xindex - self%gridslice
      zindex = zindex + 1
    enddo
    
    ! adjust the xindex, yindex
    yindex = 1
    yincrement = 1 ! incremental value to move up and down the width array
    do while (xindex.gt.2*self%width(yindex))
      ! Detect end of an x row
      xindex = xindex - 2*self%width(yindex)
      yindex = yindex + yincrement
      if (yindex.gt.self%ymax) then
        ! Detect half way through the xy slice
        yincrement = -1
        yindex = self%ymax
      elseif (yindex.le.0) then
        ! Detect end of the xy slice
        write(unit=fidlog,fmt=*)"Went off the end of width array"
        write(unit=fidlog,fmt=*)"Particle # : ",self%selector_size
        write(unit=fidlog,fmt=*)"Grid Index : ",self%selector(self%selector_size)
        write(unit=fidlog,fmt=*)"Slice      : ",self%gridslice
        do dbg_max=1,self%ymax
           write(unit=fidlog,fmt=*)"Width[",dbg_max,"]  : ",self%width(dbg_max)
        enddo

        write(unit=fidlog,fmt=*)"X,Y,Z      : ",xindex," ",yindex," ",zindex
        stop "Went off the end of width array"
      endif
    enddo

    ! create particle at current grid point
    xout = self%length * dble(xindex - self%width(yindex))
    yout = self%length * dble(-1 * yincrement * yindex)
    zout = self%length * (dble(zindex - self%zmax/2) - 0.5)
  end subroutine tubegridnext

  ! int gridsize(grid)
  !        Get the remaining number of points
  function tubegridsize(self)
    implicit none
    type (tubegrid), intent(inout) :: self
    integer :: tubegridsize

    if (dbc) then
      if (.not.gridisok(self)) then
        stop "Call to gridnext with uninitialised object"
      endif
    endif

    tubegridsize=self%selector_size
  end function tubegridsize

  ! int gridisok(grid)
  !        Test if the grid object can be used
  function tubegridisok(self)
    implicit none
    type (tubegrid), intent(inout) :: self
    logical :: tubegridisok

    tubegridisok = allocated(self%selector).and.allocated(self%width).and.self%selector_size.gt.0
  end function tubegridisok

  ! gridkill(grid)
  !        Release any allocated memory
  subroutine tubegridkill(self)
    implicit none
    type (tubegrid), intent(inout) :: self
    integer :: stat
    self%selector_size = 0
    stat = 0
    if (allocated(self%selector)) deallocate(self%selector, STAT=stat)
    if (stat.ne.0) stop "Memory deallocation of failed"
    if (allocated(self%width)) deallocate(self%width, STAT=stat)
    if (stat.ne.0) stop "Memory deallocation of failed"
  end subroutine tubegridkill

end module grid
