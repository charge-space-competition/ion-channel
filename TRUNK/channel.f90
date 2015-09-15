
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

!  Version GRS400 17
!
!    This is a major change from version 16.
!    Features:
!    - particles stored contiguously
!    - "no" code duplication of core routines
!    - use of fortran 90/95 array operations
!    - more uniform specie handling
!    - accept run number and random seed on command-line
!    - active data set is not changed until commit operation, all
!      methods operate on /trial/ before commit.
!
! EXTERNAL ROUTINES
!
! Uses dgetrs,dgetrf LAPACK routines for manipulating ICC matrix
!   - any LAPACK library should do!
!
! ranff: External function that creates double precision random 
!   numbers using the Mersenne Twister with 19937 parameters.
! sranff: Seed subroutine for the random number generator
!
! dfeq: (double precision fuzzy equals) Floating point numbers that
!   are intuitively equal are rarely equal in practise. This is due
!   to the limited precision that results in an accumulation of 
!   numerical errors. The method provided in cutil.cpp says two
!   floating point numbers are equal if they only vary in the last
!   five binary bits.
!
! DATA INPUT FILE
!
! The program now uses a single data input file. This file consists
!   of labelled sections and labelled data elements.  The data can now
!   be input in any order so requires a two stage read process. Each
!   module should implement a 'rd{section}' subroutine for the sections
!   it is responsible for. This method may be called multiple times
!   depending on the input file and section type.  A second method
!   'rf{mod}' is called after the input file has been read.  This method
!   is where input file derived data should be calculated.  The order
!   these second methods are called is fixed, so data that has inter-
!   module dependencies can be finalised there.
! 
! INPUT FINALISATION ORDER
!
! * MODULES THAT ARE FINALISED BEFORE INPUT
!
! 0: STRNGS, VERS
!
! * MODULES THAT FINALISE IN FIRST STAGE
!
! 0: CONST, CHANNEL
!
! * MODULES REQUIRING SECOND STAGE
!
! 1: SPECIE
! 2: GEOM
! 3: CONF
! 4: PATCH
! 5: TRIAL
! 6: ACCUM 
!
! New code can explicitly assume that all modules above them in
! this table are in their finalised state in the second stage
! of finalisation.  Therefore any new code written for these modules
! must fulfil the condition that the object is in its final valid
! state at the end of the second stage.
!
! EXCEPTIONS to this rule must be clearly documented _and_ any 
! access to such a feature must occur with a check for validity
! and generating a program-stopping error if feature has not been
! finalised.
!
! -------------------------------------------------------------
! Methods:
! cmdlin
!    process the command line.
! readin
!    manage process of reading in the input file
! runsim, runblk
!    perform the simulation.
! zero
!    reset statistical data; calls zeroXX subroutines of 
!    individual modules.
!
! -------------------------------------------------------------
! Data related to running the sumulations.
!
!     ninner -> Inner loop size of the MC iteration.
!     nstep -> the total number of iterations
!     naver -> the number of iterations to allow the system to equilibrate
!
program channel
use accum
use const
use conf
use geom
use patch
use spec
use trial
implicit none
  ! The size of the MC inner loop
  integer, parameter :: ninner=1024
  ! Number of iterations for: total, equilibration and bulk
  integer :: nstep,naver,nblk
  ! Signal handler's integer
  integer :: onsigterm
  ! Whether to just do the bulk simulation (default = false)
  logical :: bulk_only
  ! Integer value for update method (accum%malas1, accum%malas2, accum%accept, accum%lamperski)
  integer :: bulk_update_method
  integer :: bulk_update_variant
  logical :: bulk_clamp
  ! Title for this run
  character(1024) :: runame
  ! usechm = whether to use chemical potentials from input or perform
  !        simulation of bulk to generate them (default true)
  logical :: usechm

  ! Optional over-ride of the ntarget value from module geom.
  integer :: ntarg_opt

  ! Make the default value below 0 so that it will not override
  ! any value set in the geometry.
  ntarg_opt = -1
  nstep     = 0
  naver     = 0
  nblk      = 0
  onsigterm = 0
  bulk_only = .false.
  usechm    = .true.
  bulk_update_method = 0
  bulk_update_variant = 0
  bulk_clamp=.false.

  ! Output information about the program and host
  call prinfo(fidlog)

  ! Process command line and read program input
  call cmdlin
  write(unit=fidlog,fmt='("# UUID ",A32)')fuuid
  write(unit=fidlog,fmt='(72("-"))')
  write(unit=fidlog,fmt='(4/)')

  ! read input file

  call readin

  if (bulk_only) then
    ! Compute chem potentials only
    call genrbk
    call runblk
  else
    !   Create and save an initial configuration
    call genrcf("dat/conf."//firun//".chk")

    ! Use input chem potentials or generate them
    q9g1g4: if (.not.usechm) then
      ! Create a bulk initial configuration
      write(unit=fidlog,fmt='(72("-"))')
      write(unit=fidlog,fmt=*)"Computing chemical potentials."
      write(unit=fidlog,fmt='("Bulk simulation volume and length: ",F12.2,F12.2)')volblk(),lenblk()
      call genrbk
      call runblk
      ! reset statistics
      call zero
      ! read back the created configuration
      call readcf("dat/conf."//firun//".chk")
    endif q9g1g4

    call wrtinp("dat/initial."//firun//".inp")

    ! Run the simulation
    call runsim("dat","channel")
  endif
  contains

  ! -------------------------------------------------------------
  ! Process the command line
  subroutine cmdlin
  implicit none

  external sranff ! random seed
  ! LOCALS
  character(256) arun ! command line argument
  integer irun       ! for run number
  integer istat      ! check for bad conversions

  ! read job number
  irun=1
  u5d6d5: if (command_argument_count().ge.1) then
    call get_command_argument(1,arun)
    read(arun,fmt=*,iostat=istat)irun
    if (istat.ne.0) then
      call help(fidlog)
      stop 1
    endif
  endif u5d6d5

  ! --- make firun, character chain for index of run
  call setrun (irun)

  ! if a second value is present use it as random seed
  ! else use first value
  v3i2s8: if (command_argument_count().ge.2) then
    call get_command_argument(2,arun)
    read(arun,fmt=*,iostat=istat)irun
    if (istat.ne.0) then
      call help(fidlog)
      stop 1
    endif
    ! initialise the random number generator
    call sranff(irun)
  else v3i2s8
    ! initialise the random number generator
    call sranff(mag+irun)
  endif v3i2s8

  end subroutine cmdlin


  ! ----------------------------------------------------------------------
  ! Write out the 'channel' section parameters
  !
  subroutine ecchnl(fid)
  use strngs
  use conf
  use simstate
  implicit none
  integer, intent(in) :: fid
  character(20) :: str_out

  write(unit=fid,fmt='(A)')fschnl
  write(unit=fid,fmt='(A,1X,I7)')fsnstp,nstep
  write(unit=fid,fmt='(A,1X,I7)')fsnavr,naver
  if (nblk.ne.0) then
     write(unit=fid,fmt='(A,1X,I7)')fsnblk,nblk
  end if
  call str(tmptur(), str_out)
  write(unit=fid,fmt='(A,1X,A)')fstsi,trim(adjustl(str_out))
  write(unit=fid,fmt='(A,1X,"""",A,"""")')fsname,trim(runame)
  call str(usechm, str_out)
  write(unit=fid,fmt='(A,1X,A)')fschpt,trim(str_out)
  call str(usegrid, str_out)
  write(unit=fid,fmt='(A,1X,A)')fsgrid,trim(str_out)
  call str(byiter, str_out)
  write(unit=fid,fmt='(A,1X,A)')fsnmcf,trim(str_out)
  call str(.not.do_electrostatic(), str_out)
  write(unit=fid,fmt='(A,1X,A)')fsnoch,trim(str_out)
  select case (bulk_update_method)
    case (lamperski)
      write(unit=fid,fmt='(A,1X,A)')fsupdt,"lamperski"
    case (malas1)
      write(unit=fid,fmt='(A,1X,A)')fsupdt,"malas1"
    case (malas2)
      write(unit=fid,fmt='(A,1X,A)')fsupdt,"malas2"
    case (accept)
      select case (bulk_update_variant)
      case (0)
        write(unit=fid,fmt='(A,1X,A)')fsupdt,"accept"
      case default
        write(unit=fid,fmt='(A,1X,A,I2)')fsupdt,"accept",bulk_update_variant
      end select
  end select
  call str(bulk_clamp, str_out)
  write(unit=fid,fmt='(A,1X,A)')"bulk_clamp",trim(str_out)
  if (ntarg_opt.ne.0) then
    write(unit=fid,fmt='(A,1X,I7)')fsntrg,ntarg_opt
  endif
  call str(byiter, str_out)
  write(unit=fid,fmt='(A,1X,A)')fsnmcf,trim(str_out)
  call str(bulk_only, str_out)
  write(unit=fid,fmt='(A,1X,A)')fsbulk,trim(str_out)
  call str(dotlog, str_out)
  write(unit=fid,fmt='(A,1X,A)')"trial_log",trim(str_out)
  write(unit=fid,fmt='(A)')fsend
  write(unit=fid,fmt=*)

  end subroutine ecchnl

  ! ------------------------------------------------------------
  ! Print global help information 
  subroutine help(fid)
    implicit none
    integer, intent(in) :: fid
    if (dbc) then
      if (dbc_level.ge.dbc_require) then
        call require_writable(fid)
      endif
    endif
    write(unit=fid,fmt=*)"Usage information"
    write(unit=fid,fmt=*)
    write(unit=fid,fmt=*)"<program name> [run number [random seed]]"
    write(unit=fid,fmt=*)
    write(unit=fid,fmt='(A," : ",A)')"run number","[optional, number, default=1 (range 1 to 999)] the run number."
    write(unit=fid,fmt='(9X,A)')             "The run number is used as part of all output file names and so is the"
    write(unit=fid,fmt='(9X,A)')             "way the program allows you to distinguish between simulations. Running"
    write(unit=fid,fmt='(9X,A)')             "two simulations with the same run number in a single folder will cause"
    write(unit=fid,fmt='(9X,A)')             "the results of the two simulations to over-write each other, and should"
    write(unit=fid,fmt='(9X,A)')             "be avoided. The run number is converted to a three letter label, with"
    write(unit=fid,fmt='(9X,A)')             "'0's added at the front, when used in file names, for example run "
    write(unit=fid,fmt='(9X,A)')             "number 4 becomes '004' and 111 becomes '111'."
    write(unit=fid,fmt='(9X,A)')             "Additionally, the program can use run numbers to choose from a set of "
    write(unit=fid,fmt='(9X,A)')             "input files. When input files with names like 'channel.004.inp', where"
    write(unit=fid,fmt='(9X,A)')             "4 is the run-number, exist the program selects any matching input file"
    write(unit=fid,fmt='(9X,A)')             "allowing multiple simulations to be run in the same folder. Alternately,"
    write(unit=fid,fmt='(9X,A)')             "if only 'channel.inp' exists then the same simulation will be safely run "
    write(unit=fid,fmt='(9X,A)')             "multiple times (in the same folder) with different run numbers."
    write(unit=fid,fmt='(A," : ",A)')"random seed","[optional, number] random seed used to initialise the internal"
    write(unit=fid,fmt='(9X,A)')             "random number generator. Simulations with the same seed value and input"
    write(unit=fid,fmt='(9X,A)')             "will have identical results. If not given the 'run number' is used as"
    write(unit=fid,fmt='(9X,A)')             "the seed. (Note the ionch/utilities folder contains the 'irand32' program "
    write(unit=fid,fmt='(9X,A)')             "that can be used to generate good seed values)"

    write(unit=fid,fmt=*)
    write(unit=fid,fmt=*)"Input file information"
    write(unit=fid,fmt=*)
    call help_channel(fid)
    write(unit=fid,fmt=*)
    call help_specie(fid)
    write(unit=fid,fmt=*)
    call help_salt(fid)
    write(unit=fid,fmt=*)
    call help_subspecie(fid)
    write(unit=fid,fmt=*)
  end subroutine help

  ! ------------------------------------------------------------
  ! Print channel help information 
  subroutine help_channel(fid)
    use strngs
    implicit none
    integer, intent(in) :: fid
    if (dbc) then
      if (dbc_level.ge.dbc_require) then
        call require_writable(fid)
      endif
    endif
    write(unit=fidlog,fmt=*)"channel input section definition"
    write(unit=fid,fmt='(A6," : ",A)')fsnstp,"[required, number] target number of simulation steps"
    write(unit=fid,fmt='(A6," : ",A)')fsnavr,"[required, number] number of thermalisation steps"
    write(unit=fid,fmt='(A6," : ",A)')fsnblk,"[optional, number] number of bulk simulation steps (default="""//fsnavr//""")"
    write(unit=fid,fmt='(A6," : ",A)')fstsi,"[optional, number, default=300] simulation temperature in Kelvin"
    write(unit=fid,fmt='(A6," : ",A)')fsname,"[optional, text] title for the simulation"
    write(unit=fid,fmt='(A9," : ",A)')fsnmcf,"[optional, boolean, default=false] if true the program will write"
    write(unit=fid,fmt='(9X,A)')             "step specific 'conf' output files during each checkpoint."
    write(unit=fid,fmt='(A6," : ",A)')fsntrg,"[optional, number] if a value is given for this option here it will"
    write(unit=fid,fmt='(9X,A)')             "over-ride any value specified in the 'geom' section.  When you have"
    write(unit=fid,fmt='(9X,A)')             "the 'geom' input section in a separate file, this allows adjustment"
    write(unit=fid,fmt='(9X,A)')             "of the particle number without needing multiple geometry files."
    write(unit=fid,fmt='(A6," : ",A)')fschpt,"[optional, boolean, default=true] if false the program will perform"
    write(unit=fid,fmt='(9X,A)')             "a bulk simulation to determine the chemical potentials before the"
    write(unit=fid,fmt='(9X,A)')             "main simulation. When true the chemical potentials from the input"
    write(unit=fid,fmt='(9X,A)')             "are used."
    write(unit=fid,fmt='(A6," : ",A)')fsnoch,"[optional, boolean, default=false] if true the program will perform"
    write(unit=fid,fmt='(9X,A)')             "the simulation ignoring all charges and charge-charge interactions."
    write(unit=fid,fmt='(A7," : ",A)')fsgrid,"[optional, boolean, default=false] if true the program will use an "
    write(unit=fid,fmt='(9X,A)')             "alternate initial configuration generation method that works by placing"
    write(unit=fid,fmt='(9X,A)')             "particles on a grid. This method allows much higher initial "
    write(unit=fid,fmt='(9X,A)')             "concentrations to be used than the normal method."
    write(unit=fid,fmt='(A9," : ",A)')fsbulk,"[optional, boolean, default=false] if true the program will perform"
    write(unit=fid,fmt='(9X,A)')             "only the bulk simulation for calculating the chemical excess."
    write(unit=fid,fmt='(A9," : ",A)')fsupdt,"[optional, (malas1|malas2|accept|lamperski)] if only calculating the chemical"
    write(unit=fid,fmt='(9X,A)')             "excess, this is the update method to use."

  end subroutine help_channel

  ! ----------------------------------------------------------------------
  ! print program/system information
  !
  ! It is important to keep accurate records of the machine,
  ! compiler, libraries and program version that went to create a
  ! result.  This routine gathers and prints this information.
  !
  ! Some of this information is provided in the vers module.  The
  ! content of the vers module is created by the build system each
  ! time the program is built and contains information about how
  ! the program was compiled and any static version information.
  !
  ! The 'libver' routine gives either static or dynamic version
  ! information depending on the library used.
  !
  ! Information about the machine running the program is performed
  ! by the 'machine_info' method (in the 'cutil.cpp' file).  This
  ! prints the information directly to stdout and requires the 
  ! 'flush(fidlog)' to ensure correct ordering of output.
  !
  subroutine prinfo(fid)
  use vers
  implicit none
  integer, intent(in) :: fid
  character(32) :: date_
  integer :: idx
  double precision :: dummy
  character(2048) :: buffer_
  integer :: buffersize_
  intrinsic :: get_command_argument

    if (dbc) then
      if (dbc_level.ge.dbc_require) then
        call require_writable(fid)
      endif
    endif
  idx = 0
  dummy = 0
  write(unit=fid,fmt=*)"PROGRAM:  ",compnm
  ! run date information
  date_(1:32) = ' '
  call date_and_time(date_(1:8), date_(10:19), date_(21:25))
  write(unit=fid,fmt=*)"RUN DATE: ",trim(date_)
  write(unit=fid,fmt='(72("-"))')
  ! program compilation information
  write(unit=fid,fmt=*)"COMPILED VERSION:  ",compvr
  write(unit=fid,fmt=*)"COMPILED DATE:  ",compdt
  write(unit=fid,fmt=*)"COMPILER INFORMATION"
  write(unit=fid,fmt=*)"FORTRAN COMPILER: ",cmpflr
  write(unit=fid,fmt=*)"FORTRAN TARGET: ",cmpftg
  write(unit=fid,fmt=*)"FORTRAN FLAGS: ",cmpflg
  write(unit=fid,fmt=*)"CPP COMPILER: ",cmpclr
  write(unit=fid,fmt=*)"CPP TARGET: ",cmpctg
  write(unit=fid,fmt=*)"CPP FLAGS: ",cmpclg
  write(unit=fid,fmt=*)"-- SIZE LIMITS --"
  write(unit=fid,fmt='(A20,":",I6)')"SPECIES",nspcmx
  write(unit=fid,fmt='(A20,":",I6)')"FREE SPECIES",noccmx
  write(unit=fid,fmt='(A20,":",I6)')"SALTS",nsltmx
  write(unit=fid,fmt='(A20,":",I6)')"SALT COMPONENTS",nnewmx
  write(unit=fid,fmt='(A20,":",I6)')"TOTAL ION PARTICLES",ntotmx
  write(unit=fid,fmt='(A20,":",I6)')"STRUCTURAL IONS",nionmx
  write(unit=fid,fmt='(A20,":",I6)')"REGIONS",nrgnmx
  write(unit=fid,fmt='(A20,":",I6)')"RADIAL HISTOGRAM BINS",nrgmx
  write(unit=fid,fmt='(A20,":",I6)')"AXIAL HISTOGRAM BINS",nzgmx
  write(unit=fid,fmt='(A20,":",I6)')"INNER LOOP COUNT",ninner
  write(unit=fid,fmt='(A20,":",I6)')"ICC GRID-POINTS",npchmx

  write(unit=fid,fmt='(72("-"))')
  write(unit=fid,fmt=*)"LAPACK LIBRARY: ",trim(libver())
  ! host and environment information
  buffersize_ = 2048
  call machine_info(buffer_,buffersize_,onsigterm)
  if (buffersize_.gt.0) then
    write(unit=fid,fmt='(A)',advance='NO')buffer_(:buffersize_)
  endif

  write(unit=fid,fmt=*)"MACHINE INTEGER SIZE",bit_size(idx)
  write(unit=fid,fmt=*)"REAL NUMBERS:"
  write(unit=fid,fmt=*)" EPSILON: ",epsilon(dummy)
  write(unit=fid,fmt=*)" HUGE: ",huge(dummy)
  write(unit=fid,fmt=*)" PRECISION: ",precision(dummy)
  write(unit=fid,fmt='(72("-"))')
  ! Command line arguments
  write(unit=fid,fmt=*)"COMMAND LINE ARGUMENTS"
  do idx=1,command_argument_count()
    call get_command_argument(idx,date_)
    select case (idx)
    case (1)
      write(unit=fid,fmt='(1X,A12,": ",A)')"RUN NUMBER",trim(date_)
    case (2)
      write(unit=fid,fmt='(1X,A12,": ",A)')"RANDOM SEED",trim(date_)
    case default
      write(unit=fid,fmt='(1X,"ARGUMENT[",I2,"]: ")')idx,trim(date_)
    end select
  enddo
  write(unit=fid,fmt='(72("-"))')
  end subroutine prinfo

  ! -------------------------------------------------------------
  ! Read in the key simulation parameters
  !
  subroutine rdchnl(fid,sname,svalue,istat)
  use strngs
  use conf
  use simstate
  implicit none
  integer, intent(in) :: fid
  character(len=*), intent(in) :: sname,svalue
  integer, intent(out) :: istat
  ! default temperature
  double precision :: tsi=300.0D0
  logical, dimension(2) :: mask_
  character(32) :: nme_
  character(1024) :: val_
  logical :: nocharge
  runame='ion channel'

  if (dbc) then
    if (sname.ne.fschnl) stop "Error: incorrect section name"
    if (len_trim(svalue).ne.0) call d1s8f2("Error: section does not take any parameters")
  endif
  mask_=.false.
  bulk_clamp=.false.
  nocharge=.false.
  g0s3p2: do
    val_ = " "
    call readnv(fid,nme_,val_,istat)
    ! exit on bad read or section 'end'
    if (istat.ne.0) return
    if (nme_.eq.fsend) exit g0s3p2
    ! looking for nstep,naver,ratmv, ratslt, ratind
    j4c2w5: select case (nme_)
      case (fsnstp) j4c2w5
        read(val_,*)nstep
        mask_(1)=.true.
      case (fsnavr) j4c2w5
        read(val_,*)naver
        mask_(2)=.true.
      case (fsnblk) j4c2w5
        read(val_,*)nblk
      case (fstsi) j4c2w5
        read(val_,'(D20.13)')tsi
       case (fsntrg) j4c2w5
          read(val_,*)ntarg_opt
      case (fsname) j4c2w5
        runame=val_
        ! remove matching leading and trailing quotes
        runame=dequote(runame)
      case (fschpt) j4c2w5
        read(val_,*)usechm
      case (fsnoch) j4c2w5
        read(val_,*)nocharge
      case (fsnmcf) j4c2w5
        read(val_,*)byiter
      case (fsgrid) j4c2w5
        read(val_,*)usegrid
      case (fsbulk) j4c2w5
        read(val_,*)bulk_only
      case ("trial_log") j4c2w5
        read(val_,*)dotlog
      case ("bulk_clamp") j4c2w5
        read(val_,*)bulk_clamp
      case (fsupdt) j4c2w5
        bulk_update_variant = 0
        select case (val_(1:6))
        case ("lamper")
          bulk_update_method = lamperski
        case ("malas1")
          bulk_update_method = malas1
        case ("malas2")
          bulk_update_method = malas2
        case ("accept")
          bulk_update_method = accept
          if (len(val_).gt.6) then
            read(val_(7:),*)bulk_update_variant
          else
            bulk_update_variant = 1
          endif
        end select
      case default j4c2w5
        call d1s8f2("Name "//nme_//" is not valid in simulation parameter (channel) section")
    end select j4c2w5
  enddo g0s3p2
  w6k6z5: if (.not.all(mask_)) then
    call d1s8f2("Not all required tags were present")
  endif w6k6z5

  call use_electrostatic(nocharge)
  call incnst(tsi)
  
  end subroutine rdchnl

  subroutine d1s8f2(msg)
  use strngs
  implicit none
  character(len=*), intent(in) :: msg
    write(unit=fidlog,fmt=*)"## Bad simulation parameter (channel) section in input. ##"
    write(unit=fidlog,fmt=*)msg
    write(unit=fidlog,fmt=*)"## Bad simulation parameter (channel) section in input. ##"
    call help_channel(fidlog)
    stop 1
  end subroutine d1s8f2

  ! ----------------------------------------------------------------
  ! Read the system input file
  !
  ! This is the main input routine for the Monte Carlo simulation.
  ! It reads in the data necessary for the simulation phase and
  ! initialises the simulation data structures.
  !
  ! The input file 'channel.inp' is made up of a series of labelled
  ! subsections. The module specific subsections of the input file
  ! may appear in any order and the content of these subsections
  ! may also appear in any order. To allow sections and input data
  ! to appear in any order the read-in process is split into two
  ! phases.
  !
  ! PHASE ONE: read input file
  !
  ! In the first phase this method calls a module specific routine
  ! for each subsection (based on the subsection label).  Within
  ! a subsection these routines must handle getting input in any
  ! order, optional input and basic input validity.  In general
  ! non-optional input sets a flag that is tested when the subsection
  ! ends.  If the input values can be validated at this time they
  ! should be.  Also, any initialisation that is independent of other
  ! modules may also occur.
  !
  ! The input file can include other sub-files.  The ability to
  ! include a sub-file is only possible at the top-most level
  ! outside a subsection.  However, sub-files may include other
  ! sub-files but recursive inclusion is undefined.
  !
  ! Duplicate subsections are not detected specifically, with later
  ! subsections triggering a repeat call to the module specific input
  ! method.  For sections that are expected to appear multiple
  ! times (eg particle and salt definitions) repeat calls to the
  ! module method must provide well-defined outcomes.  For sections
  ! that are expected to appear once (eg the geometry definition)
  ! the results of a repeat call is allowed to be undefined.
  !
  ! PHASE TWO: initialisation
  !
  ! In the second phase an initialisation routine in each module
  ! is called in a specific order.  This order allows a module to
  ! be initialised after the modules it depends on.  These routines
  ! should finalise initialisation for their module, perhaps logging
  ! the consequences of the input data.  The routines are also
  ! required to write to the log an input file subsection that
  ! contains the input data (and any defaults).  Input data that
  ! is normalised should be shown as the normalised value, otherwise
  ! the logged values should be the same as the input file.
  !
  ! NOTE: it is important that the logged input data subsections
  ! can be used to perform an exact repeat of the current calculation.
  !
  ! INPUT ERROR HANDLING:
  !
  ! In phase one, if any subroutine detects syntatic invalid or
  ! missing input it should print an error message detailing the
  ! problem and what correct input should be.  In phase two, errors
  ! are more likely to be range errors.  These should also be
  ! reported, along with suggestions for making the input valid.
  !

  subroutine readin
  use accum
  use conf
  use geom
  use patch
  use spec
  use strngs
  use simstate
  implicit none

  integer :: fid
  integer :: istat ! IO error
  logical, dimension(7) :: mask_
  character(32) :: nme_
  character(1024) :: val_
  integer :: ispec
  logical :: isexst ! does file exist
  character(32) :: isrdbl ! is file readable
  
  write(unit=fidlog,fmt=*)"Checking for """//fschnl//"."//firun//".inp"" as input file."
  ! PHASE ONE
  fid=fidbas
  ! try openning 'channel.XXX.inp' before 'channel.inp'
  inquire(file=fschnl//"."//firun//".inp", iostat=istat, exist=isexst, read=isrdbl)
  if (isexst) then
    if (istat.eq.0) then
      istat=1
      if (isrdbl.ne.'NO') then
        open(fid,file=fschnl//"."//firun//".inp", iostat=istat, action="read")
        if (istat.eq.0) then
          write(unit=fidlog,fmt=*)"Using """//fschnl//"."//firun//".inp"" as input file."
        endif
      endif
    endif
  else
    istat=1
  endif
  if (istat.ne.0) then
    write(unit=fidlog,fmt=*)"File """//fschnl//"."//firun//".inp"" not found, checking for """//fschnl//".inp""."
    inquire(file=fschnl//".inp", iostat=istat, exist=isexst, read=isrdbl)
    if (isexst) then
      if (istat.eq.0) then
        istat=1
        if (isrdbl.ne.'NO') then
          open(fid,file=fschnl//".inp", iostat=istat, action="read")
          if (istat.eq.0) then
            write(unit=fidlog,fmt=*)"Using """//fschnl//".inp"" as input file."
          endif
        endif
      endif
    else
      istat=1
    endif
  endif
  if (istat.ne.0) then
    write(unit=fidlog,fmt=*)"Neither """//fschnl//"."//firun//&
         &".inp"" nor """//fschnl//".inp"" name a readable file in the current directory."
    stop "No input file"
  endif
  mask_=.false.
  m4h5o1: do
    val_ = " "
    call readnv(fid,nme_,val_,istat)
    ! exit loop on EOF
    if (istat.lt.0) then
      ! End of file
      close(fid)
      if (fid.le.fidbas) exit m4h5o1
      fid=fid-1
      cycle m4h5o1
    endif
    ! error on bad read
    if (istat.ne.0) call c7x0s5("File read error occured", mask_)
    ! looking for section names:
    !  fschnl -> 'channel'
    !  fsgeom -> 'geom'
    !  fsconf -> 'conf'
    !  fsptch -> 'patch'
    !  fsregn -> 'region'
    !  fssalt -> 'salt'
    !  fstry  -> 'trial'
    !  fsaccu -> 'accum'
    !  fsspec -> 'specie'
    !  fssubs -> 'subspecie'
    !  fsincl -> 'include' a subfile.
    s0i4x7: select case (nme_)
      case (fsincl) s0i4x7
        fid=fid+1
        inquire(file=trim(val_), iostat=istat, exist=isexst, read=isrdbl)
        if (isexst) then
          if (istat.eq.0) then
            istat=1
            if (isrdbl.ne.'NO') then
              open(fid,file=trim(val_),action="read",iostat=istat)
              if (istat.eq.0) then
                write(unit=fidlog,fmt=*)"Reading input file: """//trim(val_)//"""."
              else
                call c7x0s5("Error occurred opening input file """//trim(val_)//""".", mask_)
              endif
            else
              call c7x0s5("Included input file """//trim(val_)//""" is not readable.", mask_)
            endif
          else
            call c7x0s5("Error occurred looking up input file """//trim(val_)//""".", mask_)
          endif
        else
          call c7x0s5("Included input file """//trim(val_)//""" was not found.", mask_)
        endif
      case (fschnl) s0i4x7
        call rdchnl(fid,nme_,val_,istat)
        mask_(1)=(istat.eq.0)
      case (fsgeom) s0i4x7
        call rdgeom(fid,nme_,val_,istat)
        mask_(2)=(istat.eq.0)
      case (fsptch) s0i4x7
        call rdptch(fid,nme_,val_,istat)
        mask_(3)=(istat.eq.0)
      case (fssubs) s0i4x7
        call rdsubs(fid,nme_,val_,istat)
        ! TODO remove sharing of mask with salt
        mask_(4)=(istat.eq.0)
      case (fssalt) s0i4x7
        call rdsalt(fid,nme_,val_,istat)
        mask_(4)=(istat.eq.0)
      case (fsspec) s0i4x7
        call rdspec(fid,nme_,val_,istat)
        mask_(5)=(istat.eq.0)
      case (fstry) s0i4x7
        call rdtral(fid,nme_,val_,istat)
        mask_(6)=(istat.eq.0)
      case (fsaccu) s0i4x7
        call rdaccu(fid,nme_,val_,istat)
        mask_(7)=(istat.eq.0)
      case (fsfver) s0i4x7
        filcur=0
        read(val_,*)filcur
        if (fvermx.lt.filcur) then
          call c7x0s5("File version "//val_//" is too recent for this program version", mask_)
        endif
        if (0.ge.filcur) then
          call c7x0s5("File version "//val_//" is invalid", mask_)
        endif
      case default s0i4x7
        call c7x0s5("Name "//nme_//" is not valid in simulation parameter (channel) section", mask_)
    end select s0i4x7
    ! check iostat was not set in a subroutine
    ! (this is error here for EOF too!)
    if (istat.ne.0) call c7x0s5("File read error occured", mask_)
  enddo m4h5o1
  w6k6z5: if (.not.all(mask_)) then
    call c7x0s5("Not all required sections were present", mask_)
  endif w6k6z5

  do while (fid.ge.fidbas)
    close(unit=fid)
    fid=fid-1
  enddo

  ! PHASE TWO
  ! input file successfully read, call finalisation methods

  call rfspec(rl(1), zlimit())
  call rfgeom(ntarg_opt)
  call rfconf
  if (do_electrostatic()) call rfptch
  call rftral
  call rfaccu

  ! System initialised, print input data summary
  write(unit=fidlog,fmt='(3/)')
  write(unit=fidlog,fmt='(72("-"))')
  write(unit=fidlog,fmt=*)"SIMULATION SUMMARY"
  write(unit=fidlog,fmt='(72("-"))')
  ! ---------------------------------------------------------------
  write(unit=fidlog,fmt=*)
  write(unit=fidlog,fmt=*)" Simulation parameters in SI units"
  write(unit=fidlog,fmt=*)
  write(unit=fidlog,fmt='(1X,A27,F10.2)')" T [K] = ",1/(beta()*boltz)
  ! print reduced charge factor
  write(unit=fidlog,fmt='(1X,A27,F10.5)')"reduced charge factor = ",qstar()
  write(unit=fidlog,fmt=*)
  write(unit=fidlog,fmt='(1X,A27,F10.2)')" filter radius = ", rl(1)
  write(unit=fidlog,fmt='(1X,A27,F10.2)')" vestibule outer radius = ", rl(2)
  write(unit=fidlog,fmt='(1X,A27,F10.2)')" flat wall outer radius = ", rl(3)
  write(unit=fidlog,fmt='(1X,A27,F10.2)')" outer cylinder radius = ", rl(4)
  write(unit=fidlog,fmt='(1X,A27,F10.5)')" cell cylinder radius = ", rl(5)
  write(unit=fidlog,fmt=*)
  write(unit=fidlog,fmt='(1X,A27,F10.2)')" filter length = ", zl(1)
  write(unit=fidlog,fmt='(1X,A27,F10.2)')" channel length = ", zl(2)
  write(unit=fidlog,fmt='(1X,A27,F10.2)')" outer cylinder length = ", zl(3)
  write(unit=fidlog,fmt='(1X,A27,F10.5)')" cell length = ", zl(4)
  write(unit=fidlog,fmt=*)
  write(unit=fidlog,fmt='(1X,A27,F10.2)')" vestibule radius = ", rlvest()
  write(unit=fidlog,fmt='(1X,A27,F10.2)')" curvature radius = ", rlcurv()
  write(unit=fidlog,fmt=*)
  write(unit=fidlog,fmt='(1X,A27,F10.2)')" eps solvent = ", epsw
  write(unit=fidlog,fmt='(1X,A27,F10.2)')" eps protein = ", epspr
  write(unit=fidlog,fmt='(52("-"))')
  write(unit=fidlog,fmt='(1X,A6,1X,A8,3X,A7,3X,A7)')"Specie","N(input)","q [e]","d [A]"
  write(unit=fidlog,fmt='(52("-"))')
  o7a4o6: do ispec=1,nspec()
    write(unit=fidlog,fmt='(1X,A6,5X,I4,3X,F7.2,3X,F7.2)')fspc(ispec),ni(ispec),xz(ispec),xri(ispec)*2
  enddo o7a4o6
  write(unit=fidlog,fmt=*)
  write(unit=fidlog,fmt='(72("-"))')
  call ecchnl(fidlog)
  write(unit=fidlog,fmt='(72("-"))')
  write(unit=fidlog,fmt=*)"Initialisation of simulation complete"
  write(unit=fidlog,fmt='(72("-"))')
  write(unit=fidlog,fmt=*)

  end subroutine readin

  ! ERROR REPORT routine for 'readin'
  subroutine c7x0s5(msg,mask)
  use strngs
  implicit none
  character(len=*), intent(in) :: msg
  logical, dimension(:), intent(in) :: mask
  integer :: idx ! loop index
  integer :: siz
  character(128) :: ok_sec ! string for completed section names
    write(unit=fidlog,fmt=*)"A problem occurred reading the input file."
    write(unit=fidlog,fmt=*)msg
    write(unit=fidlog,fmt=*)"Required sections:"
    write(unit=fidlog,fmt='(1X,A,6(",",1X,A))')fschnl,fsgeom,fsptch,fssalt,fsspec,fstry,fsaccu
    write(unit=fidlog,fmt=*)"Optional parameter: ",fsfver
    write(unit=fidlog,fmt=*)
    siz=1
    ok_sec(1:128)=' '
    do idx=1,7
      if (mask(idx)) then
        select case (idx)
        case (1)
          ok_sec(siz:siz+len_trim(fschnl))=trim(fschnl)
          siz=siz+len_trim(fschnl)+1
        case (2)
          ok_sec(siz:siz+len_trim(fsgeom))=trim(fsgeom)
          siz=siz+len_trim(fsgeom)+1
        case (3)
          ok_sec(siz:siz+len_trim(fsptch))=trim(fsptch)
          siz=siz+len_trim(fsptch)+1
        case (4)
          ok_sec(siz:siz+len_trim(fssalt))=trim(fssalt)
          siz=siz+len_trim(fssalt)+1
        case (5)
          ok_sec(siz:siz+len_trim(fsspec))=trim(fsspec)
          siz=siz+len_trim(fsspec)+1
        case (6)
          ok_sec(siz:siz+len_trim(fstry))=trim(fstry)
          siz=siz+len_trim(fstry)+1
        case (7)
          ok_sec(siz:siz+len_trim(fsaccu))=trim(fsaccu)
          siz=siz+len_trim(fsaccu)
        end select
      endif
    enddo
    write(unit=fidlog,fmt=*)"Successfully read sections: ",trim(ok_sec)
    stop 1
  end subroutine c7x0s5

  ! -------------------------------------------------------------
  ! Perform the simulation of ion channel
  !
  subroutine runsim(dirname,filename)
  use accum
  use conf
  use spec
  use trial
  use patch
  use simstate
  implicit none
  character(len=*),intent(in) :: dirname,filename
  ! LOCALS
  integer :: i,j,istep
  logical :: iswidm
  double precision :: dgfree
  integer :: mixtme

  call start_cell

  ! Initialise patch vectors c and h
  if (do_electrostatic().and..not.is_homogeneous()) then
    call genrch
  endif

  write(unit=fidlog,fmt=*)' START GCMC for ',trim(runame),' (run ',firun,')'
  write(unit=fidlog,fmt=*)
  ! turn off Widom sampling in Equilibration
  iswidm=calwid
  calwid=.false. 
  !  EQUILIBRATION
  !
  ! Markov-chain Monte Carlo has mixing time scaling with n log(n) where
  ! n is the degrees of freedom.  In our case the degrees of freedom is
  ! equal to the number of particles by three and the number of ICC tiles.
  ! I use particles * 3 because they have 3D movement and 1 * tiles as they
  ! only vary charge.

  dgfree = 3 * ntot() + npatch
  mixtme = max(1,nint(dgfree * log(dgfree) / dble(ninner)))
  write(unit=fidlog,fmt=*)"Degrees of freedom factor   : ",nint(dgfree)
  write(unit=fidlog,fmt=*)"Scale factor for mixing time: ",mixtme
  write(unit=fidlog,fmt=*)"Trials per simulation step  : ",mixtme*ninner
  write(unit=fidlog,fmt=*)

  if (naver.ne.0) then
    write(unit=fidlog,fmt='(72("-"))')
    write(unit=fidlog,fmt=*)" THERMALISATION STAGE 1: Constrained GC"
    m2j0e4: do istep=1,naver
      do i=1,mixtme*ninner
        if (onsigterm.ne.0) goto 100
        call trypeq
        call thermal_accmlt
      end do
    ! NO-LOCAL RDF  call print_rdf
    enddo m2j0e4
    call concentration_report(fidlog)
    call wrtinp(dirname//"/equil."//firun//".inp")
    call writcf(dirname//"/equil."//firun//".chk",1)
    write(unit=fidlog,fmt='(72("-"))')
    write(unit=fidlog,fmt=*)" THERMALISATION STAGE 2: Full GC"
    m2j0e3: do istep=1,naver
      do i=1,mixtme*ninner
        if (onsigterm.ne.0) goto 100
        call tryeq
        call thermal_accmlt
      end do
    ! NO-LOCAL RDF  call print_rdf
    enddo m2j0e3
    call concentration_report(fidlog)
    ! Output statistics even during thermalization
    call saves(0,mixtme*ninner)

    call wrtinp(dirname//"/equil."//firun//".inp")
    call writcf(dirname//"/equil."//firun//".chk",1)
    call zero
  endif

  ! restore calwid
  calwid=iswidm

  if (nstep.eq.0) stop "No simulation steps requested"

  ! check mobile ions are now in valid positions.
  if (have_localized()) then
    call chkmob
  endif

  ! SIMULATION
  write(unit=fidlog,fmt=*)" START SIMULATION"
  j = 0
  r3a4i2: do istep=0,nstep,isave
    do j=1,isave
      f3y5z7: do i=1,mixtme*ninner
        if (onsigterm.ne.0) goto 100
        call tryeq
        call accmlt
      enddo f3y5z7
      call hist
  ! NO-LOCAL RDF    call print_rdf
    enddo
    write(unit=fidlog,fmt=*)" MAIN SIMULATION"
    call saves(istep+isave,mixtme*ninner)
    call wrtinp(dirname//"/"//filename//"."//firun//".inp")
    call writcf(dirname//"/"//filename//"."//firun//".chk",istep+isave)
  enddo r3a4i2

100  if (onsigterm.ne.0) then
    call wrtinp(dirname//"/"//filename//"."//firun//".inp.TERM")
    call writcf(dirname//"/"//filename//"."//firun//".chk.TERM",istep+isave)
!    call saves(istep+j,mixtme*ninner)
  endif
    
  write(unit=fidlog,fmt=*)" End of simulation"

  end subroutine runsim

  ! -------------------------------------------------------------
  ! Perform the simulation of bulk solution
  !
  subroutine runblk
    use accum
    use conf
    use trial
    use spec
    use simstate
    use rnstat
    implicit none
    logical, parameter :: isbulk=.true.
    ! chempi sum variable (chptav) and counter (chpcnt) for calculating average
    ! chemical potentials.  The average potential is used to update the chemical
    ! potentials at the end of the second and third estimation cycles.
    double precision, dimension(nspcmx) :: chdelta
    type (runstat), dimension(nspcmx) :: chptav
    integer :: ispec
    ! Mix time.
    double precision :: dgfree
    double precision, dimension(2) :: tolerance
    integer :: mixtme, nbulk
    ! LOCALS
    integer :: i,j,istep
    do ispec=idxcl(),nspec()
      call rs_init(chptav(ispec))
    end do

    call start_bulk
    call zero
    write(unit=fidlog,fmt='(72("-"))')
    write(unit=fidlog,fmt=*)"References for deriving chemical potentials:"
    write(unit=fidlog,fmt=*)"Attila Malasics, Dirk Gillespie and Dezso"" Boda ""Simulating" 
    write(unit=fidlog,fmt=*)"prescribed particle densities in the grand canonical"
    write(unit=fidlog,fmt=*)"ensemble using iterative algorithms"", The Journal of Chemical "
    write(unit=fidlog,fmt=*)"Physics, 2008, 128, 124102"
    write(unit=fidlog,fmt=*)
    write(unit=fidlog,fmt=*)"S. Lamperski, Mol. Simul. 33, 1193, 2007"

    ! We calculate the mix time as N.log(N) where N is the degrees of freedom. In the bulk simulation
    ! this is three times the number of particles, and for the channel it includes the number of patches.
    ! This number is then divided by 'ninner' and rounded to the nearest integer. This number is then
    ! multiplied by ninnner in inner loops; resulting in an integral multiple of ninner for the inner  
    ! loops.
    dgfree = 3 * ntot()
    if (.not.isbulk) then
       dgfree  = dgfree + npatch
    endif
    mixtme = max(1,nint(dgfree * log(dgfree) / dble(ninner)))

    write(unit=fidlog,fmt=*)"BULK: Degrees of freedom factor   : ",nint(dgfree)
    write(unit=fidlog,fmt=*)"BULK: Scale factor for mixing time: ",mixtme
    write(unit=fidlog,fmt=*)"BULK: Trials per simulation step  : ",mixtme*ninner
    write(unit=fidlog,fmt=*)"BULK: Length of PBC Cube          : ",lenblk()
    write(unit=fidlog,fmt=*)

    if (nblk.eq.0) then
       nbulk = min(100,naver)
    else
       nbulk = nblk
    endif

    chdelta=0
    if (bulk_update_method.eq.0) bulk_update_method = malas1
    write(unit=fidlog,fmt=*)
    write(unit=fidlog,fmt=*)" Start chemical excess in bulk calculation"
    ! We do ten "accept" runs to rapidly approach the chemical excesses
    ! when we are performing the estimation as part of a simulation.
    if (.not.bulk_only) then
       do istep=1,10
          do j=1,isave
             do i=1,ninner*mixtme
                if (onsigterm.ne.0) goto 100
                call tryblk(.true.)
                call bulk_accmlt(accept)
             enddo
             call reset_wasdel
          end do
    ! NO-LOCAL RDF      call print_rdf
          call iterat(accept,0)
       end do
    end if
    do istep=1,nbulk
       do ispec=idxcl(),nspec()
          chdelta(ispec) = chexi(ispec)
       enddo
       do j=1,isave
          do i=1,ninner*mixtme
             if (onsigterm.ne.0) goto 100
             call tryblk(bulk_clamp)
             call bulk_accmlt(bulk_update_method)
          enddo
          if (bulk_update_method.eq.accept) call reset_wasdel
       end do
    ! NO-LOCAL RDF   call print_rdf
       tolerance = 0.D0
       do ispec=idxcl(),nspec()
          tolerance(2) = tolerance(2) + abs(bulk_concentration(ispec) - ctargi(ispec))
       enddo
       call iterat(bulk_update_method,bulk_update_variant)
       do ispec=idxcl(),nspec()
          if (istep.gt.nbulk/2) then
             call rs_push(chptav(ispec), chexi(ispec))
          end if
          tolerance(1) = tolerance(1) + abs(chdelta(ispec) - chexi(ispec))
       enddo
       write(unit=fidlog,fmt='(72("-"))')
       write(unit=fidlog,fmt=*)"CHEMICAL EXCESS CHANGE : ", tolerance(1)
       write(unit=fidlog,fmt=*)"Concentration Error : ", tolerance(2)
    enddo

100 if (rs_count(chptav(idxcl())).gt.1) then
       do ispec=idxcl(),nspec()
         call chexi_set(ispec, rs_mean(chptav(ispec)))
       end do
       write(unit=fidlog,fmt='(72("-"))')
       write(unit=fidlog,fmt=*)"FINAL ESTIMATE of CHEMICAL POTENTIALS"
       write(unit=fidlog,fmt='(72("-"))')
       call output_specie_potentials(fidlog)
       write(unit=fidlog,fmt='(72("-"))')
       call output_salt_potentials(fidlog)
       write(unit=fidlog,fmt='(72("-"))')
    endif
    if (onsigterm.ne.0) then
       write(unit=fidlog,fmt=*)" End of chemical excess in bulk calculation after receiving SIGTERM"
       stop "End of simulation after receiving SIGTERM"
    endif
    ! we do not rescale chptav to sum here as it is the last time
    ! we use it.
    call zero
    write(unit=fidlog,fmt=*)" End of chemical excess in bulk calculation"
    write(unit=fidlog,fmt=*)

  end subroutine runblk

  subroutine wrtinp(filename)
  use strngs
  implicit none
  character(len=*), intent(in) :: filename

  open(unit=fidinp,file=filename,action="write")
    write(unit=fidinp,fmt='(A,1X,I3)')fsfver,fvermx
    write(unit=fidinp,fmt='("# UUID ",A32)')fuuid
    write(unit=fidinp,fmt=*)
    ! writing sections in order:
    !  'channel'
    call ecchnl(fidinp)
    !  'geom'
    call ecgeom(fidinp)
    !  'patch'
    call ecptch(fidinp)
    !  'trial'
    call ectral(fidinp)
    !  'accum'
    call ecaccu(fidinp)
    !  'salt'
    call ecsalt(fidinp)
    !  'subs'
    call ecsubs(fidinp)
    !  'specie' (via conf module)
    call wrconf(fidinp)
  close(unit=fidinp)
  end subroutine wrtinp

  ! -------------------------------------------------------------
  ! Ensure statistics data is (re)set to zero
  !
  subroutine zero
  use accum
  use trial
    call zeroac
    call zeroav
    call zerotr
  end subroutine zero

end program channel

