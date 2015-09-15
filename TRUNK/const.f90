
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


!  ----------------------------------------------------------------------
!  GLOBAL CONSTANTS AND GENERIC ROUTINES
!
!  This file is where program constants and interfaces to external
!  routines needed in more than one module are to be placed.  It
!  is also the place for helper routines that are not specific
!  to the program.
!
!  Global string constants are contained in 'strngs' module.
!
!  ----------------------------------------------------------------------
!  ARRAY AND LOOP SIZES
!  
!  Public definitions of array sizes that must be known in
!  multiple modules.
!
!  ----------------------------------------------------------------------
!  SET INDICES
!
!  Fixed indices used across modules. These include named
!  constants for specific regions and for file ids.
!
!  ----------------------------------------------------------------------
!  PHYSICAL CONSTANTS
!  
!  The set of physical constants used in this program.  There
!  are standard and derived constants.  Some of the derived
!  constants may change depending on input, so the programmer
!  is required to verify they have their final value if they
!  need to use them during initialisation.
!
!  ----------------------------------------------------------------------
!  PROGRAM CONSTANTS
!
!  Constants that control the inclusion of debugging and 
!  design-by-contract code in the program.  
!
!  firun - Global string containing the current run number
!  fuuid - Global string containing a unique run identifier
!
!  ----------------------------------------------------------------------
!  INTERFACES TO EXTERNAL METHODS
!
!  This includes the _fuzzy_ double precision equality test
!  method 'dfeq' and the random number generator 'ranff'
!
!  ----------------------------------------------------------------------
!  HELPER FUNCTIONS
!  
!  This contains a set of generic routines that can be used
!  anywhere in the program.
! 
!   sqr - square of input
!   getnth - get index of 'nth' occurence of a value in an array.
!   getnof - use one array as indices to look for 'nth' occurence of a 
!            value in a second array.  
!   isort - sort an array of integers assuming sorted with one
!            new value appended
!   iswap - swap two integer variables
!   iranff - generate a random integer in a  specified range 
!   next2 - next higher power of 2
!
module const
  implicit none
  private
  ! ----------------------------------
  ! ARRAY AND LOOP SIZES
  ! ----------------------------------

  ! The kind to use for integer counters
  integer, public, parameter :: Cntr_K = selected_int_kind(18)

  ! maximum number of particles per specie
  integer, public, parameter :: nionmx=2048
  ! Maximum number of particles
  integer, public, parameter :: ntotmx=8192
  ! Maximum number of ICC patches
  integer, public, parameter :: npchmx=2048
  ! Maximum number of species
  integer, public, parameter :: nspcmx=16
  ! Maximum number of salts
  integer, public, parameter :: nsltmx=8
  ! Maximum nuber of particles per salt
  integer, public, parameter :: nnewmx=4
  ! Maximum nuber of regions
  integer, public, parameter :: nrgnmx=4

  ! Maximum number of histogram bins in z direction
  integer, public, parameter :: nzgmx=4096

  ! Default random seed 
  integer, public, parameter :: mag=12584210

  ! ----------------------------------
  ! SET INDICES
  ! ----------------------------------

  ! REGIONS
  !
  ! Note: regions 1-4 overlap and are concentric.

  ! Region in channel defined by +/-zlim
  integer, public, parameter :: izlim=1

  ! Region in channel defined by +/-zl1
  integer, public, parameter :: ifilt=2

  ! Region from channel end to end
  integer, public, parameter :: ichan=3

  ! Region outside channel
  integer, public, parameter :: ibulk=4

  ! FILE ID CONSTANTS
  !
  ! All files are accessed through file id numbers.  To avoid
  ! accidentally writing to the same file id in more than one
  ! place, the file ids are specified as parameters here.  The
  ! comments have the corresponding source and target file name.

  ! ALL: _stdout_
  integer, public, parameter :: fidlog=6

  ! channel.f90: 'channel.inp'
  ! [NOTE: nested input file ids are incremented from this value.]
  integer, public, parameter :: fidbas=33

  ! Log file for trial moves.
  integer, public, parameter :: fidmlg=30

  !  patch.f90: 'amx.XXX.dat'
  integer, public, parameter :: fidamx=12

  ! channel.f90: 'conf.XXX' 
  integer, public, parameter :: fidinp=13

  ! conf.f90: 'conf.XXX'
  integer, public, parameter :: fidcnf=14

  ! conf.f90: 'occvol.XXX.dat'
  integer, public, parameter :: fidvol=15

  ! accum.f90: 'vjz.XXX.dat'
  integer, public, parameter :: fidvjz=16

  ! accum.f90: 'o.XXX'
  integer, public, parameter :: fidooo=17

  ! accum.f90: 'gin-MM.XXX.dat'
  integer, public, parameter :: fidgin=18

  ! accum.f90: 'gREG-MM.XXX.dat'
  integer, public, parameter :: fidgrg=19

  ! accum.f90: 'occ.XXX'
  integer, public, parameter :: fidocc=20

  ! accum.f90: 'h.XXX'
  integer, public, parameter :: fidhpc=21

  ! accum.f90: 'aREG-MM.XXX'
  integer, public, parameter :: fidarg=22

  ! accum.g90: 'gz-MM.XXX'
  integer, public, parameter :: fidgzz=23

  ! accum.f90: 'rdf-MM-MM.XXX'
  integer, public, parameter :: fidrdf=24

  ! trial.f90: 'wid-MM.XXX.dat'
  integer, public, parameter :: fidwid=25

  ! patch.f90: 'patchgeom.XXX.dat'
  integer, public, parameter :: fidpch=26

  ! patch.f90: 'patchgeom.XXX.dat'
  integer, public, parameter :: fidtry=31

  ! ----------------------------------
  ! PHYSICAL CONSTANTS
  ! ----------------------------------
  ! 'pi'
  double precision, public, parameter :: pi=3.141592653589793D0
  ! Conversion from particles/per Ang**3 to Molar {to S.I.}
  ! [ (10-30)m3 / (N_av)mol * L / (10-3)m3 ] ~ [ (10-27/N_av) L/mol ]
  double precision, public, parameter :: tosi=1660.539276735512625080121D0

  ! 1 / kT for temperature of 300K (access using 'beta()')
  ! 2.41440919407021130D20 (/J)
  double precision, private :: b1o4x5=0
  ! Convert Amgstrom to meters (/m)
  double precision, public, parameter :: dsi=1.D-10
  ! Epsilon (unit-less)
  double precision, public, parameter :: eps0=8.8542D-12
  ! Avogardo's Number (not used in program, but listed as it is the
  ! value used in calculating 'tosi' (N)
  double precision, public, parameter :: avog=6.02214D23

  ! Boltzmann's constant (J/K)
  double precision, public, parameter :: boltz=1.3806D-23
  ! Charge of an electron (V)
  double precision, public, parameter :: echg=1.6021917D-19

  ! Charge factor (access using 'qstar()')
  double precision, private :: q7a0y5=0

  ! Temperature in Kelvin (no external access)
  ! 300.0D0
  double precision, private :: t8u5y5=0

  ! ----------------------------------
  ! PROGRAM CONSTANTS
  ! ----------------------------------
  ! Should debugging output be used?
  logical, public, parameter :: debug=.false.

  ! Should design-by-contract be used
  logical, public, parameter :: dbc=.true.
  ! What level of dbc?
  integer, public, parameter :: dbc_level=0

  ! What do the levels mean?  The design-by-contract philosphy is
  ! to apply a set of progressively more intrusive tests to help
  ! maintain a valid program state. In a completely tested program
  ! one could remove all of these tests.  In practice, leaving the
  ! 'dbc_require' tests can provide some protection against 
  ! invalid input with less performance impact than the others.

  ! CODE USAGE: Because fortran has no shortcut if statements
  ! it is safest to write multiple if clauses.  As all tests are
  ! for program parameters, this code should be completely 
  ! removed during 'dead code removal' optimisation.
  !
  ! if (dbc) then
  !   if (dbc_level.ge.dbc_index) then
  !     ! Checking for a valid index
  !   endif
  ! endif

  ! Exhaustively check for index validity
  integer, public, parameter :: dbc_index=6
  ! Exhaustively check for mathmatical domain violation
  integer, public, parameter :: dbc_valid=5
  ! Check an intemediate results is in expected range
  integer, public, parameter :: dbc_check=4
  ! Check for data-set validity. Objects and modules
  ! can provide an 'invariant' method that checks that
  ! the data-set internal state is consistent. If an 
  ! 'initialise' procedure is supplied then it is 
  ! legitimate to say an invariant only holds after 
  ! calling 'initialise'.
  integer, public, parameter :: dbc_invariant=3
  ! At the end of a procedure check the result is in
  ! the advertised domain.
  integer, public, parameter :: dbc_ensure=2
  ! At the beginning of a procedure check the arguments are
  ! in the advertised domain.
  integer, public, parameter :: dbc_require=1




  ! The run number as a three character filename extension
  character(3), public :: firun

  ! The unique identifier for a simulation
  character(32), public :: fuuid

  ! The input file version number
  integer, public, parameter :: filver = 1

  ! The input file version for file being read
  integer, public :: filcur

  ! The maximum input file version number the program understands
  integer, public, parameter :: fvermx = 1

  ! ----------------------------------
  ! EXTERNAL FUNCTION INTERFACES
  ! ----------------------------------
  interface
    pure elemental function dfeq(lhs,rhs)
      logical :: dfeq
      double precision, intent(in) :: lhs, rhs
    end function dfeq
    function ranff()
      double precision :: ranff
    end function ranff
  end interface

  interface swap
    module procedure swap
    module procedure iswap
  end interface swap

  public :: dfeq, ranff, iranff, sqr, tmptur, swap
  public :: beta, qstar, incnst, next2, next64, setrun, getnth, getnof, isort

  ! Error condition tests.
  public :: require_readable, require_writable

  public :: get_pi

contains
  double precision function get_pi()
    implicit none
    get_pi = pi
  end function get_pi

  ! ----------------------------------
  ! PHYSICAL CONSTANTS
  ! ----------------------------------
  ! Get the current beta (1/kT)
  ! (units 1/J)
  double precision function beta()
    implicit none
    if (dbc) then
      if (dbc_level.ge.dbc_require) then
        if (dfeq(b1o4x5,0.D0)) stop "Called function beta before it was initialised"
      end if
    endif
    beta=b1o4x5
  end function beta

  ! ----------------------------------
  ! Initialise program constants
  !
  ! This subroutine should be one of the first subprograms called.
  ! It should be called before calculating any random numbers or
  ! reading in the species.
  ! 
  ! + Set temperature
  ! + Set beta (1/kT) for a temperature other than the default 300K.
  ! + Set factor that converts specie valency to q
  !
  ! @param Temperature : The new temperature in Kelvin
  subroutine incnst(temperature)
    implicit none
    double precision, intent(in) :: temperature
    if (dbc) then
      if (dbc_level.ge.dbc_require) then
        if (temperature.lt.260) stop "Temperature is in Kelvin so T<260 makes no sense"
        if (temperature.gt.380) stop "Temperature is in Kelvin so T>380 makes no sense"
      endif
    endif
    t8u5y5=temperature
    b1o4x5=1 /(boltz*temperature)
    q7a0y5=echg*dsqrt(b1o4x5/(4*pi*eps0*dsi))
  end subroutine incnst

  ! ------------------------------------------------------------
  ! Check if the fid represents a readable file.
  !
  ! This routine calls 'stop' if the file is not readable
  subroutine require_readable(fid)
    implicit none
    integer, intent(in) :: fid
    logical :: isopn_
    integer :: iostt_
    character(16) :: isred_
    inquire(unit=fid,iostat=iostt_,action=isred_,opened=isopn_)
    if (0.ne.iostt_) then
      stop "Error: I/O error reading from unit ID"
    elseif (.not.isopn_) then
      stop "Error: Unit ID is unopen"
    elseif (.not.(isred_.eq.'READ'.or.isred_.eq.'READWRITE'.or.isred_.eq.'read'.or.isred_.eq.'readwrite')) then
      write(unit=fidlog,fmt=*)"Unit ID ",fid," is opened for :",trim(isred_)," need 'READ' or  'READWRITE'"
      stop "Error: Passed unit ID is unreadable"
    endif
  end subroutine require_readable

  ! ------------------------------------------------------------
  ! Check if the fid represents a writable file.
  !
  ! This routine calls 'stop' if the file is not writable
  subroutine require_writable(fid)
    implicit none
    integer, intent(in) :: fid
    logical :: isopn_
    integer :: iostt_
    character(16) :: isred_
    inquire(unit=fid,iostat=iostt_,action=isred_,opened=isopn_)
    if (0.ne.iostt_) then
      stop "Error: I/O error reading from unit ID"
    elseif (.not.isopn_) then
      stop "Error: Unit ID is unopen"
    elseif (.not.(isred_.eq.'WRITE'.or.isred_.eq.'READWRITE'.or.isred_.eq.'write'.or.isred_.eq.'readwrite')) then
      write(unit=fidlog,fmt=*)"Unit ID ",fid," is opened for :",trim(isred_)," need 'WRITE' or  'READWRITE'"
      stop "Error: Passed unit ID is unwritable"
    endif
  end subroutine require_writable

  ! ----------------------------------------
  ! Simulation run temperature
  double precision function tmptur()
    implicit none
    if (dbc) then
      if (dfeq(t8u5y5,0.D0)) stop "Called function tmptur before it was initialised"
    endif
    tmptur=t8u5y5
  end function tmptur

  ! ----------------------------------
  ! Factor to convert specie valency to q
  double precision function qstar()
    implicit none
    if (dbc) then
      if (dbc_level.ge.dbc_require) then
        if (dfeq(q7a0y5,0.D0)) stop "Called function qstar before it was initialised"
      endif
    endif
    qstar=q7a0y5
  end function qstar

  ! ----------------------------------
  ! PROGRAM CONSTANTS
  ! ----------------------------------
  ! Initialise firun and fuuid strings
  !
  ! This subroutine converts the input run number into a three
  ! digit string used to generate run specific output file names.
  ! 
  ! @param irun : The number of this run.
  subroutine setrun(irun)
    implicit none
    integer, intent(in) :: irun
    if (dbc) then
      if (dbc_level.ge.dbc_require) then
        if (irun.le.0) stop "Run number must be in range 1 to 999"
        if (irun.ge.1000) stop "Run number too large, must be less than 1000"
      endif
    endif
    write(firun,'(I0.3)')irun
    call uuidgn(fuuid)
  end subroutine setrun

  ! ----------------------------------
  ! HELPER FUNCTIONS
  ! ----------------------------------

  ! -------------------------------------------------------------
  ! Return position in arr(1:isize) of ith occurance of val 
  !
  ! Returns isize+1 if no element found
  integer function getnth(arr,val,i,isize)
    integer, intent (in) :: val,i,isize
    integer, dimension(:), intent(in) :: arr
    integer :: ii,j
    if (dbc) then
      if (dbc_level.ge.dbc_index) then
        if (i.lt.1) stop "Error: in search, 'nth' item must be greater than zero"
        if (i.gt.isize) stop "Error: number of values greater than array size"
        if (size(arr).lt.isize) stop "Error: array size smaller than specified"
      endif
    endif
    j=i
    getnth=1
    d7y5m5: do ii=1,isize
      getnth=ii
      if (arr(ii).eq.val) j=j-1
      if (j.eq.0) return
    enddo d7y5m5
    getnth=isize+1
  end function getnth

  ! -------------------------------------------------------------
  ! Return position in arr1(1:isize) of ith occurance of val in
  ! arr2(arr1(X)) 
  !
  ! Returns isize+1 if no element found
  integer function getnof(arr1,arr2,val,i,isize)
    integer, intent (in) :: val,i,isize
    integer, dimension(:), intent(in) :: arr1, arr2
    integer :: ii,j
    integer :: arr2_size
    if (dbc) then
      if (dbc_level.ge.dbc_index) then
        if (i.lt.1) stop "Error: in search, 'nth' item must be greater than zero"
        if (i.gt.isize) stop "Error: number of values greater than array size"
        if (size(arr1).lt.isize) stop "Error: array size smaller than specified"
      endif
    endif
    j=i
    getnof=1
    arr2_size = size(arr2)
    d7y5m5: do ii=1,isize
      getnof=ii
      if (dbc) then
        if (dbc_level.ge.dbc_index) then
          if (arr1(ii).lt.1) stop "Error: in getnof search, arr1(i) item must be greater than zero"
          if (arr1(ii).gt.arr2_size) stop "Error: in getnof search, arr1(i) greater than array2 size"
        endif
      endif
      if (arr2(arr1(ii)).eq.val) j=j-1
      if (j.eq.0) return
    enddo d7y5m5
    getnof=isize+1
  end function getnof

  ! -------------------------------------------------------------
  ! Sort an array of integers
  !
  ! @param ilist : array to sort
  ! @param isize : number of elements in array to sort
  ! @param lascnd : ascending if true, else descending
  !
  ! Assumes nearly sorted order with new element appended
  subroutine isort(ilist,isize,lascnd)
    implicit none
    integer, dimension(:), intent(inout) :: ilist
    integer, intent(in) :: isize
    logical, intent(in) :: lascnd
    integer i_
    if (dbc) then
      if (dbc_level.ge.dbc_index) then
        if (0.ge.isize) stop "Error: sort array of zero or less size"
        if (size(ilist).lt.isize) stop "Error: array is smaller than requested size"
      endif
    endif
    if (isize.lt.2) return
    i_=isize
    n7q3i2: do while (i_.gt.1)
      i_=i_-1
      u4k1g8: if (lascnd.neqv.(ilist(i_+1).gt.ilist(i_))) then
        call iswap(ilist(i_+1),ilist(i_))
        i_=isize
      endif u4k1g8
    enddo n7q3i2
  end subroutine isort

  ! ----------------------------------
  ! Swap two integers
  !
  pure elemental subroutine iswap(ii,jj)
    implicit none
    integer, intent(inout) :: ii,jj
    integer tmp
    tmp=ii
    ii=jj
    jj=tmp
  end subroutine iswap

  ! ----------------------------------
  ! Swap two doubles
  !
  pure elemental subroutine swap(ii,jj)
    implicit none
    double precision, intent(inout) :: ii,jj
    double precision tmp
    tmp=ii
    ii=jj
    jj=tmp
  end subroutine swap

  ! ----------------------------------
  ! Calculate a squared
  pure double precision function sqr(a)
    implicit none
    double precision, intent(in) :: a
    sqr=a*a
  end function sqr

  ! ----------------------------------
  ! Use ranff to generate an integer between 1 and a
  !
  ! This routine uses the 'ranff' external funtion, which uses
  ! the Mersenne Twister 19937 to generate random numbers.
  !
  ! @param a : maximum integer desired
  !
  ! @require a >= 1
  ! @ensure 1 <= iranff <= a
  integer function iranff(a)
    implicit none
    integer, intent(in) :: a
    integer :: a_
    logical :: isneg
    if (dbc) then
      if (dbc_level.ge.dbc_require) then
        if (a.eq.0) stop "Error: maximum value must not be 0"
      end if
    endif
    a7y9b1: if (a.lt.0) then
      a_=-a
      isneg=.true.
    else
      a_=a
      isneg=.false.
    endif a7y9b1
    u8f9k5: if (a_.eq.1) then
      iranff = 1
    else u8f9k5
      iranff = ceiling(ranff()*dble(a_))
    endif u8f9k5
    ! handle case where ranff ~= 0
    if (iranff.lt.1) iranff = 1
    ! handle case where ranff ~= 1
    if (iranff.gt.a_) iranff = a_
    if (isneg) iranff=-iranff
  end function iranff

  ! Nearest power of 2

  integer pure function next2(a)
    implicit none
    integer, intent(in) :: a
    next2=a-1
    next2 = ior(ishft(next2,-1),next2)
    next2 = ior(ishft(next2,-2),next2)
    next2 = ior(ishft(next2,-4),next2)
    next2 = ior(ishft(next2,-8),next2)
    next2 = ior(ishft(next2,-16),next2)
    next2=next2+1
  end function next2

  ! Nearest multiple of 64

  integer pure function next64(a)
    implicit none
    integer, intent(in) :: a
    integer :: mod_a
    mod_a = mod(a,64)
    next64=a
    if (mod_a.ne.0) next64=next64 - mod_a + 64
  end function next64


end module const
