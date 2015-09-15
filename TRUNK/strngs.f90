
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
!  GLOBAL TEXT CONSTANTS AND STRING METHODS
!
!  ----------------------------------------------------------------------
!  INPUT FILE STRINGS
!
!  Set of constant strings recognised in program input file
!
!  ----------------------------------------------------------------------
!  STRING METHODS
!
!  str - convert a value into a string (currently only logicals)
!  readnv - read a line of input and split into a name, value pair
!
module strngs
implicit none

interface str
  module procedure logstr
  module procedure fmtfln

  subroutine redflt(strval, a_array, a_size, a_cnt, a_istt)
    character(*), intent(in) :: strval
    double precision, dimension(*), intent(inout) :: a_array
    integer, intent(in) :: a_size
    integer, intent(out) :: a_cnt, a_istt
  end subroutine redflt
end interface str
  !  -------------------
  !  INPUT FILE STRINGS
  !  -------------------

  ! input file sections
  character(*), parameter :: fschnl='channel'
  character(*), parameter :: fsgeom='geom'
  character(*), parameter :: fsconf='conf'
  character(*), parameter :: fsptch='patch'
  character(*), parameter :: fsregn='region'
  character(*), parameter :: fssalt='salt'
  character(*), parameter :: fstry='trial'
  character(*), parameter :: fsaccu='accum'
  character(*), parameter :: fsspec='specie'
  character(*), parameter :: fssubs='subspecie'
  character(*), parameter :: fsfver='fileversion'

  ! include subfiles option name
  character(*), parameter :: fsincl='include'

  ! generic option names
  character(*), parameter :: fsend='end'
  character(*), parameter :: fsname='name'

  ! simulation option names
  character(*), parameter :: fsnstp='nstep'
  character(*), parameter :: fsnavr='naver'
  character(*), parameter :: fsnblk='nbulk'
  character(*), parameter :: fsnmcf='multiconf'
  character(*), parameter :: fschpt='usepot'
  character(*), parameter :: fsgrid='usegrid'
  character(*), parameter :: fsbulk='bulkonly'
  character(*), parameter :: fsupdt='updatef'
  character(*), parameter :: fstsi='kelvin'
  character(*), parameter :: fsnoch='nocharge'

  ! statistic option names
  character(*), parameter :: fscgin='calgin'
  character(*), parameter :: fscrdf='calrdf'
  character(*), parameter :: fsclac='calacc'
  character(*), parameter :: fsclmb='calmob'
  character(*), parameter :: fswidm='calwid'
  character(*), parameter :: fsiwid='iwidom'
  character(*), parameter :: fsdrg='drg'
  character(*), parameter :: fsisav='isave'
  character(*), parameter :: fsgzoc='zocc'

  ! trial module option names
  !      uses  fsrtmv='ratmov'
  character(*), parameter :: fsadd='add'
  character(*), parameter :: fsdrmi='drmaxin'
  character(*), parameter :: fsdrmo='drmaxout'
  character(*), parameter :: fsrtsl='ratslt'
  character(*), parameter :: fsrtid='ratind'
  character(*), parameter :: fsrtjp='ratjmp'
  character(*), parameter :: fskmob='mobk'

  character(*), parameter :: fsdzg='dzg'
  ! geometry option names
  character(*), parameter :: fsgzl1='zl1'
  character(*), parameter :: fsgzl4='zl4'
  character(*), parameter :: fsgrl1='rl1'
  character(*), parameter :: fsgrl4='rl4'
  character(*), parameter :: fsgrl5='rl5'
  character(*), parameter :: fsgrlv='rlvest'
  character(*), parameter :: fsgrlc='rlcurv'
  character(*), parameter :: fsgzlm='zlimit'
  character(*), parameter :: fsntrg='ntarg'
  character(*), parameter :: fsnsrt='oldreg'

  ! specie type parameter names
  character(*), parameter :: fschex='chex'
  character(*), parameter :: fsd='d'
  character(*), parameter :: fsn='n'
  character(*), parameter :: fsrtgr='ratgr'
  character(*), parameter :: fsrtex='ratexc'
  character(*), parameter :: fsrtmv='ratmov'
  character(*), parameter :: fsrtsp='ratspc'
  character(*), parameter :: fsrtrg='ratreg'
  character(*), parameter :: fsz='z'
  character(*), parameter :: fstype='type'
  ! (allowed specie type values)
  character(*), parameter :: fsmobl='mob'
  character(*), parameter :: fsflxd='flex'
  character(*), parameter :: fschon='chonly'
  character(*), parameter :: fsfree='free'

  ! salt parameter names
  character(*), parameter :: fsislt='cation'
  character(*), parameter :: fsctrg='ctarg'

  ! subspecie parameter names
  character(*), parameter :: fsenth='enthalpy'
  character(*), parameter :: fsentr='entropy'
  character(*), parameter :: fsgrnd='ground'
  character(*), parameter :: fsexct='excited'
  character(*), parameter :: fsrtsw='ratswap'

  ! patch module option names
  character(*), parameter :: fsdxf='dxf'
  character(*), parameter :: fsdxw='dxw'
  character(*), parameter :: fsnsub='nsub'
  character(*), parameter :: fsepsp='epspr'
  character(*), parameter :: fsepsw='epsw'

contains
  ! --------------------------------------------
  ! Give a text representation of a logical value
  !
  ! Converts a logical value into '.true.' or
  ! '.false.'
  pure subroutine logstr(a_bool, output)
    implicit none
    logical, intent(in) :: a_bool
    character(*), intent(inout) :: output
    character(*), parameter :: yes = ".true."
    character(*), parameter :: no = ".false."
    if (a_bool) then
      output=yes(1:max(1,min(len(output),len(yes))))
    else
      output=no(1:max(1,min(len(output),len(no))))
    endif
  end subroutine logstr

  ! --------------------------------------------
  ! READ NAME,VALUE PAIR
  !
  ! Reads line from unit=fid, ignoring blank lines
  ! and deleting comments beginning with "#". Comments
  ! may appear anywhere on line.
  !
  ! Splits the line into a name, value pair
  !
  ! limits: lines longer than 1024 will be truncated
  !       : long names/values will be truncated to length
  !         of sname and svalue supplied by the caller.
  !
  ! @param fid : Input file id
  ! @param sname : (OUT) the name
  ! @param svalue : (OUT) the value (may be empty)
  ! @param istat : any read error condition (-1 => EOF)
  subroutine readnv(fid,sname,svalue,istat)
    use const
    implicit none
    integer, intent(in) :: fid ! File number
    character(*), intent(out) :: sname,svalue
    integer, intent(out) :: istat

    ! LOCALS
    character(1024) :: line
    integer :: ipos ! position of divider char

    if (dbc) then
      ! check for valid 'fid'
      if (dbc_level.ge.dbc_require) call require_readable(fid)
    end if
    istat=0
    do
      read(unit=fid,fmt='(A)',iostat=istat)line
      if (istat.ne.0) return
      ! remove comments (iostat used to report error)
      line=decmnt(line)

      if (len_trim(line).ne.0) then
        ! have name/value (look for ' ' or '='
        ipos=scan(line,' =')
        if (dbc) then
          if (dbc_level.ge.dbc_check) then
            if (ipos.ge.0) stop "Result from string 'scan' less than 0"
            if (ipos.le.len(line)) stop "Result from string 'scan' greater than string length"
          end if
        end if
        if (ipos.eq.0) then
          !     have only name
          sname=line
          svalue=' '
        else
          !     name and value
          sname=line(:ipos-1)
          svalue=adjustl(line(ipos+1:))
        endif
        return
      endif
    enddo
  end subroutine readnv

  ! ------------------------------------------------------------
  ! Method to pretty format floating point fields.
  !
  ! This method tries to select the best way to represent the number in
  ! the given width. Firstly one character is always reserved for the
  ! minus sign.  Numbers that can be nicely represented without exponent
  ! are formatted in fixed point, otherwise the numbers are formatted
  ! using an exponent.  Any decimal precision is reduced to the minimum
  ! representable by the type.  The formatted number is then checked for
  ! unnecessary zeros at the end of number, which are removed.
  !
  ! Examples
  ! -11.234 : Only as many decimal places as are non-zero are shown
  ! -0.3456 : Exponent is not used for numbers close to one
  ! -3.456E-001
  subroutine fmtfln(num, output)
    implicit none
    double precision, intent(in) :: num
    character(len=*), intent(inout) :: output
    interface 
      subroutine fmtfl2(val,str,sz)
        double precision, intent(in) :: val
        character(*), intent(inout) :: str
        integer, intent(in) :: sz
      end subroutine fmtfl2
    end interface
    call fmtfl2(num,output,len(output))
  end subroutine fmtfln
    

  subroutine fmtflt(num, output)
    use const
    implicit none
    double precision, intent(in) :: num
    character(len=*), intent(inout) :: output
    integer :: wid
    character(2) :: frmt1, frmt2, frmt3 ! format components
    character(10) :: frmt   ! Constructed format 
    double precision :: val ! Absolute value of num
    integer :: zeros   ! number of useless zeros
    integer :: precsn  ! number of decimal places
    integer :: decpt   ! index of decimal point
    integer :: expnt   ! index of beginning of exponent or string length
    integer :: i       ! loop index

    wid=len(output)
    if (dbc) then
      if (dbc_level.ge.dbc_require) then
        if (wid.lt.3) stop "Calling fmtflt with a string shorter than 3 makes no sense"
      end if
    end if
    ! write(*,*)wid
    write(frmt2,fmt='(I2)')wid

    ! ----
    ! treat zero as a special case
    if (dfeq(num,0.D0)) then
      output="0.00"
    else
      ! ----
      ! non-zero
      val=abs(num)
      ! Work in loop as added check to avoid
      ! '*'s
      do while(.true.)
        precsn=wid
        ! ----
        ! Choose best format letter and precision
        if (1.D0.gt.val) then
          ! have decimal
          if (val.ge.0.001D0) then
            ! decimal close to one
            frmt1='F'
            precsn=wid-3
          else
            ! decimal not close to one
            frmt1='ES'
            if (1.D-100.lt.val) then
              precsn=wid-7
            else
              precsn=wid-8
            endif
          endif
        else if (10D0**(wid-1).gt.val) then
          ! have number that can be expressed as simple decimal
          frmt1='F'
          precsn=wid-max(1,ceiling(log10(val))-2)
          if (num.ne.val) precsn=precsn-1
        else
          ! decimal not close to one
          frmt1='ES'
          if (1.D100.gt.val) then
            precsn=wid-7
          else
            precsn=wid-8
          endif
        endif
        if (precsn.lt.0) precsn = 0
        ! ----
        ! reduce precision to maximum useful
        do 
          if (1.D0.ne.(1.D0+(10D0**(-precsn)))) exit
          precsn = precsn-1
        enddo
        write(frmt3,fmt='(I2)')precsn

        ! ----
        ! concatenate format parts and generate initial result
        frmt='('//frmt1//trim(adjustl(frmt2))//'.'//trim(adjustl(frmt3))//')'
        write(output,fmt=trim(frmt))num
        !
        ! Check for '*'s
        if (scan(output,'*').eq.0) then
          ! No '*'s, exit loop
          exit
        end if
        ! reduce wid
        wid = wid - 1
        if (wid.lt.3) then
          write(fidlog,*)"WARNING: Failed call to fmtflt(",num,",S[",len(output),"])"
          ! stop "Width in fmtflt reached 2 without valid output"
          return
        end if
      end do

      ! ----
      ! try to eliminate zeros
      if (frmt1.ne.'F') then
        ! have exponent
        decpt=scan(output,'.')
        expnt=scan(output,'ED')-1
        ! check for lack of exponent
        if (expnt.le.0) expnt=len_trim(output)
      else
        ! no exponent
        decpt=scan(output,'.')
        expnt=len_trim(output)
      endif
      if (expnt.gt.(decpt+4).and.decpt.ne.0) then
        zeros=0
        ! use decpt+2 so we have at least two decimal places
        do i=expnt-1,decpt+3,-1
          if (output(i:i).eq.'0') then
            zeros=zeros+1
          else
            exit
          endif
        enddo
        if (output(expnt:expnt).eq.'0'.or.zeros.gt.2) then
          ! for E we need to shift exponent
          if (frmt1.ne.'F') then
            output(expnt-zeros:)=output(expnt+1:)
          else
            output(expnt-zeros:expnt)=' '
          endif
        endif
      endif
    endif
    ! ----
    ! Make final result right justified
    frmt1='A'
    frmt='('//frmt1//trim(frmt2)//')'
    write(output,fmt=frmt)trim(output)
  end subroutine fmtflt

  ! ----------------------------------------------------------------------
  ! Convert a string of numbers into a double precision array
  !
  ! NOTE: this method assumes that the only spacing characters in strval
  ! will be spaces or tabs
  !
  subroutine redMYfltSS(strval,a_arry,a_size,a_cnt,a_istt)
    use const
    implicit none
    character(len=*), intent(in) :: strval ! the line of numbers 
    double precision, dimension(*), intent(inout) :: a_arry ! the target array
    integer, intent(in) :: a_size ! the maximum number of values to read
    integer, intent(out) :: a_cnt ! the actual number read
    integer, intent(out) :: a_istt ! iostat code of last read, /= 0 means error
    integer :: left,right,end ! substr limits
    integer :: idx ! array index
    double precision :: val
    character(32) :: strsub ! sub-string of strval
    character(2) :: spcset
    spcset(1:1) = " "
    spcset(2:2) = achar(9)
    right = 1
    end = len_trim(strval)
    !if (dbc) then
    !  if (dbc_level.ge.dbc_require) then
    !    if (size(a_arry).lt.a_size) stop "In redflt: Given array length is larger than the array"
    !  end if
    !end if
    do idx=1,a_size
      ! set left to next non-blank character
      left = right - 1 + verify(strval(right:),spcset)
      if (left.eq.0.or.left.ge.end) then
        ! end of string
        return
      endif

      ! set right to first blank after left
      right = scan(strval(left:),spcset)
      if (right.eq.0.or.right.gt.end - left + 2) then
        right = end + 1
        strsub=strval(left:right)
      else
        right = left - 1 + right
        strsub=strval(left:right - 1)
      endif

      ! read value
      a_istt = 0
      read(strsub,fmt='(D20.13)',iostat=a_istt)val
      if (a_istt.ne.0) then
        write(unit=fidlog,fmt=*)"ERROR: problem converting """,val,""" in a floating point number."
        stop "Error converting string into a real number"
      endif
      ! update array if not error
      a_arry(idx) = val
      a_cnt = idx
    enddo
  end subroutine redMYfltSS


! ----------------------------------------------------------------------
! Dequote a string
!
! Attempt to remove comments from a string
pure function decmnt(string)
  implicit none
  character(len=*), intent(in) :: string
  character(len=0), parameter :: nullstr = ""
  character(len=len_trim(string)) :: decmnt
  integer :: ipos ! position of '#' char
  ! remove leading whitespace
  decmnt=adjustl(string)

  ! Check for comment lines
  ipos=index(decmnt,'#')
  if (ipos.eq.1) decmnt=nullstr

  ! Remove end-of-line comments 
  if (ipos.ne.0) decmnt=decmnt(:ipos-1)
end function decmnt

! ----------------------------------------------------------------------
! Dequote a string
!
! Attempt to remove leading and trailing quotes from a string.
pure function dequote(string)
  implicit none
  character(len=*), intent(in) :: string
  character(len=len_trim(string)) :: dequote
  character(*), parameter :: quotes = """'"
  dequote=trim(adjustl(string))
  if (1.eq.scan(dequote,quotes)) then
    if (dequote(1:1).eq.dequote(len_trim(dequote):len_trim(dequote))) then
      dequote=dequote(2:len_trim(dequote)-1)
    endif
  endif
end function dequote

end module strngs
