
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


! TRIAL DATA
!
! This module contains the main work methods of the MC and GCMC
! step.  The module proposes a new configuration and tests this
! new configuration for acceptance.  Only if the trial configuration
! is accepted does the module update the main configuration.
!
! NOTE: the previous version (version 16) updated part of the main
! configuration during a trial and wound back the change if a move
! was not accepted.  Looking at the code showed that both schemes
! required the same amount of data transfer, however modifying
! the main configuration only when the trial was accepted made
! the 'trial' and 'conf' modules more independent.  The modify
! -on-accept scheme was therefore implemented.  As a side benefit
! implementing test methods that gather data by attempting to
! insert or delete particles to calculate the energy change but
! never commit the change are straightforward to implement.
!
! -------------------------------------------------------------
! STANDARD TRIAL MOVE SEQUENCE
!
! This module attempts to use a similar sequence of operations
! for all move attempts.  In particular there is only one method
! to calculate the particle-particle interaction and one
! method to calculate the particle-patch interaction which are
! only accessible to the 'moveeq' method.
!
! The driver routines for GCMC steps are the 'tryXXX' methods.
! These methods all follow the same basic stages of making a
! move attempt:
!
! (1) Select the move type using 'selXXX' method
!     - set istate,iregnw,ispcnw,igcnw using variation of
!       selXXX method
! (2) Generate move, create or destroy using a 'chgXXX' method.
!     - set some of rxnw,rynw,rznw,r2nw,indxnw. These methods
!       are responsible for checking the new positions do not
!       overlap the system geometry and to sets the volume
!       parameter needed when adding or deleting particles.
! (3) Calculate the change in energy using 'moveeq'
!     - uses the cergit and pergit method for calculating energy
!       change.
!     - Decides to accept or reject change
!     - If accepted, call 'commit' to update 'conf' module
!
! -------------------------------------------------------------
!   hnw,dcnw -> trial h and change in c vectors 
!   ripnw -> particle-patch distances for new position
!   riinw -> inter-particle distances for new position
!   rxnw,rynw,rznw,r2nw -> x,y,z,r of changing particles. If delete
!     this contains the old value, if move or add it is the new values
!   riicrn -> inter-particle distances between the particles being added
!   amove,ajump,ajin,ajout -> particle move statistics
!   acreat,adest,acrea1,adest1 -> particle add/remove statistics
!   indxnw -> global index list of change
!   ispcnw -> specie list of change
!   nchgnw -> the number of particles that are changing
!   igcnw -> the salt that is changing
!   iregnw -> where change is happening
!   istate - change state flags
!     CHANGE TYPE
!     istate[0] = 0,+/-1,+/-2
!      0 = move (no change in nactv)
!      +/-1 = add or delete individual ion
!      +/-2 = add or delete salt
!     MOVETYPE
!     istate[1] = 1,2,3,4
!      1 = sperical move
!      2 = jump (in bulk or non-structural, in channel for 
!           structural)
!      3 = jump into channel
!      4 = jump out of channel
!
! ----------------------------------------------------------------------
! MODULE PUBLIC OPERATIONS AND ATTRIBUTES
!
! accrat
!    Print the GCMC move acceptance data if calacc is .true.
!
! accget
!    Get the log of the ratio of create attempts over delete attempts.
!
! avergi
!    Print information about the Metropolis acceptance values
!
! avgget
!    Get difference between create and destroy attempts divided
!    by the total number of attempts.
!
! calwid  ('calwid' input variable of 'accum' section)
!    Boolean value indicating if data for the Widom insertion method
!    for calculating the Chem. Pot. should be collected
!
! ectral
!    Echo the attributes of the module in input file format.
!
! nwdtry  ('iwidom' input variable of 'accum' section) 
!    The Widom method used in this code can make use of the single
!    particle insertion moves of the GCMC.  A number of insertion trials
!    will therefore be done for "free".  However, this can lead to 
!    different sampling across the simulation.  This variable sets the
!    increment in the number of trials that must exist at the beginning 
!    of each 'save' step (minimum trials = nwdtry * save number).  If 
!    too few trials exist in a bin, then insertion attempts are made until
!    this minimum is reached.
!
!    As the minimum number is nwdtry * save number it is correlated with
!    the number of steps in each save interval ('isave' input variable).
!
!    * TODO *  CHANGE TO NUMBER OF TRIALS IN A FIXED NUMBER OF STEPS
!    TO ELIMINATE CORRELATION WITH 'isave'
!    ie  minimum = { nwdtry / steps } * step-count
!
! rdtral
!    subroutine for reading 'trial' section of input file
!
! rftral
!    subroutine to finialise initialisation of module in 'readin'
!    process
!
! savewd
!    Subroutine for publishing the Widom insertion trial results.  It is
!    this method that optional performs insertion trials to reach the 
!    minimum as well as writing a summary in the log file and in the
!    output file.
!
! tryblk
!    Driver routine for a single trial of the bulk simulation with
!    periodic boundary conditions.
!
! tryeq
!    Driver routine for a single trial of the membrane/ion channel
!    simulation within a simulation box
!
! trypeq
!    Driver routine for a single trial in the equilibraion phase of
!    the ion-channel simulation
!
! trywd
!    Driver routine for a Widom insertion trial.  Unlike other 'try*'
!    routines this performs multiple trials as needed to reach minimum
!    in each histogram bin
!
! zeroav
!    Reset the acceptance ratio data separately from other
!    statistical data.
!
! zerotr
!    Reset statistical data recorded during a simulation.  
!
! ----------------------------------------------------------------------
module trial
  use const
  use rnstat
  use spec, only : sphr_l, jump_l, jmpin_, jmpou_, swap_l
  implicit none
  private
  ! logical, private, parameter :: dbctrial=.true.
  ! logical, private, parameter :: debugtrial=.true.
  ! ISTATE values
  ! Constants defining the type of move
  integer, private, parameter :: move_l=0
  integer, private, parameter :: add1_l=1
  integer, private, parameter :: addslt=2
  integer, private, parameter :: rem1_l=-1
  integer, private, parameter :: remslt=-2

  ! Enumeration of Chem. Pot. Update types (see subroutine 'accum::iterat' for
  ! more details)
  integer, parameter, public :: malas1=1, malas2=2, accept=3, lamperski=4


  ! Constant multiplier for mobile ion constraint. Note this
  ! is the maximum value that the penalty energy can contribute
  ! to the potential energy as "mobpen" scales Hooke's law
  ! potential energy by rsr^2 (k=1). As mobchk constrains ion 
  ! within rsr, mobpen will give value between [0:1].
  double precision, private :: k_mobl = 0.0D0

  ! Variables for estimation of the chemical potential
  logical, private, dimension(nspcmx) :: wasdel=.false.
  double precision, private, dimension(:), allocatable :: hnw,dcnw
  double precision, private, dimension(:,:), allocatable :: ripnw
  double precision, private, dimension(:,:), allocatable :: riinw
  double precision, private, dimension(nnewmx) :: rxnw,rynw,rznw,r2nw
  double precision, private, dimension(nnewmx,nnewmx) :: riicrn
  double precision, private, dimension(:), allocatable :: dh

  double precision, private, dimension(2,nspcmx) :: amove,ajump,ajin,ajout
  double precision, private, dimension(nspcmx) :: aovrlp
  double precision, private, dimension(2,nrgnmx,nsltmx) :: acreat,adest
  double precision, private, dimension(2,nrgnmx,nspcmx) :: acrea1,adest1

  ! Widom insertion method uses the 'gz' z-axial bins to
  ! use the average densities from there
  type widm_t
    ! mean field (not scaled by qstar())
    double precision, allocatable, dimension(:,:) :: umfi
    ! interparticle coulomb only interaction
    double precision, allocatable, dimension(:,:) :: uchsi
    ! total interaction
    double precision, allocatable, dimension(:,:) :: uchs
    ! exponential interparticle coulomb only interaction
    double precision, allocatable, dimension(:,:) :: euchsi
    ! exponential total interaction
    double precision, allocatable, dimension(:,:) :: euchs
    ! number of tests
    !   - attempts
    integer (kind=Cntr_K), allocatable, dimension(:,:) :: trials
    !   - successes
    integer (kind=Cntr_K), allocatable, dimension(:,:) :: succss
    ! Total number of attempts
    integer (kind=Cntr_K) :: trys
  end type widm_t

  type (widm_t), private, pointer :: widobj

  ! Running total of change in energy (average and variance)
  type (runstat), private :: adeltaf
  ! Running total of change in energy per move type
  type (runstat), private, dimension(nspcmx) :: amoveng, achgeng, adeleng

  ! Rates for different MC trials: movement, salt add/del, indiv. add/del, swap
  ! This array contains the probabilities of the top level trial types so after 
  ! normalization sum(rates_) == 1
  ! INVARIANT sum(rates_) == 1
  double precision, private, dimension(4) :: rates_=0.D0

  ! Update radii inside and outside the channel
  double precision, private :: drmax(2)=0.D0
  integer, private, dimension(nnewmx) :: indxnw,ispcnw
  ! NOTE: 
  ! The value of istate(1) (move_l, add1_l, addslt, rem1_l, remslt)
  ! is also the net change in particle number.
  !
  ! The value of istate(2) is the subtype change iff istate(1)=0
  integer, private, dimension(2) :: istate
  integer, private :: nchgnw,igcnw,iregnw
  ! Proportionality factors (formerly nvterm and chmpot)
  double precision, private :: propscal, propexpo

  ! only check salt/indiv add/del if true
  logical, private :: dosalt=.false., doindv=.false.,doswap=.false.

  ! Store information on the average RDF
  type (runstat), dimension(nspcmx,nspcmx), private :: rdf_est

  ! ----------------------------------------------------------------------
  !
  !  PUBLIC ROUTINES AND ATTRIBUTES

  ! Are we using Widom calculation?
  logical, public :: calwid =.false.

  ! the target number of widom trials per save
  integer, public :: nwdtry =0

  ! report chemical potentials
  public :: savewd

  ! Print statistical data
  public :: accrat, print_rdf

  ! Simulation-quality monitoring statistic
  public :: avergeng, avergi, accget, zeroav

  ! Simulation update driver routine
  public :: tryeq, trypeq, tryblk, trywd

  ! Data read and initialisation
  public :: ectral, rdtral, rftral, zerotr

  ! Data logging
  public :: logmov
  logical, public :: dotlog=.false.
contains
 ! --------------------------------------------------
  ! Write out acceptance ratios
  !
  ! The information output by this method can be used to evaluate
  ! the performance of the parameters and methods used in the GCMC step.
  subroutine accrat(fid)
    use conf
    use geom
    use spec
    use simstate
    implicit none
    integer, intent(in):: fid
    double precision, dimension(8) :: mratio
    integer :: igc,ireg,ispec
    n5x3j7: if (.not.dfeq(ratmv(),0.D0)) then
      if (is_bulk()) then
        write(unit=fid,fmt='(36("-"))')
        write(unit=fid,fmt='(3X,"|     MOVE ","|     JUMP ")')
        write(unit=fid,fmt='(35("-"))')
        d3y8z7: do ispec=idxcl(),nspec()
          mratio = 0
          if (amove(1,ispec).ne.0) mratio(1)=amove(2,ispec)/amove(1,ispec)
          if (ajump(1,ispec).ne.0) mratio(2)=ajump(2,ispec)/ajump(1,ispec)
          write(unit=fid,fmt='(A2," ",2("|",2X,F8.5))')fspc(ispec),mratio(1),mratio(2)
        enddo d3y8z7
      else
        write(unit=fid,fmt='(72("-"))')
        write(unit=fid,fmt='(3X,"|     MOVE ","|     JUMP "&
             &,"|       IN ","|      OUT "&
             &,"|    ACC IN ","|   ACC OUT ")')
        write(unit=fid,fmt='(72("-"))')
        d3y8z6: do ispec=1,nspec()
          mratio = 0
          o9z4l2: if (.not.isfree(ispec)) then
            if (amove(1,ispec).ne.0) mratio(1)=amove(2,ispec)/amove(1,ispec)
            if (ajump(1,ispec).ne.0) mratio(2)=ajump(2,ispec)/ajump(1,ispec)
            write(unit=fid,fmt='(A2," ",2("|",2X,F8.5))')fspc(ispec),mratio(1),mratio(2)
          else o9z4l2
            if (amove(1,ispec).ne.0) mratio(1)=amove(2,ispec)/amove(1,ispec)
            if (ajump(1,ispec).ne.0) mratio(2)=ajump(2,ispec)/ajump(1,ispec)
            if (ajin(1,ispec).ne.0) mratio(3)=ajin(2,ispec)/ajin(1,ispec)
            if (ajout(1,ispec).ne.0) mratio(4)=ajout(2,ispec)/ajout(1,ispec)               
  
            write(unit=fid,fmt='(A2," ",4("|",2X,F8.5),2("|",I11))')&
                 &        fspc(ispec),mratio(1),mratio(2),mratio(3),mratio(4)&
                 &        ,nint(ajin(2,ispec)),nint(ajout(2,ispec))
          endif o9z4l2
        enddo d3y8z6
      endif
    endif n5x3j7

    d0o2i7: if (.not.dfeq(ratslt(),0.D0)) then
      if (is_bulk()) then
        write(unit=fid,fmt='(36("-"))')
        write(unit=fid,fmt='(1X,A4,1X,2(1X,A7))')'SALT','CREATE','DESTROY'
        write(unit=fid,fmt='(36("-"))')
        a3k1v5: do igc=1,nsalt()
          ispec=isalt(igc)
          mratio=0
          if (acreat(1,ibulk,igc).ne.0) mratio(1)=acreat(2,ibulk,igc)/acreat(1,ibulk,igc)
          if (adest(1,ibulk,igc).ne.0) mratio(2)=adest(2,ibulk,igc)/adest(1,ibulk,igc)
          write(unit=fid,fmt='(1X,A4,1X,2(1X,F7.4))')fsalt(igc),mratio(1),mratio(2)
        enddo a3k1v5
      else
        write(unit=fid,fmt='(72("-"))')
        write(unit=fid,fmt='(1X,A4,1X,2(1X,A30,1X))')'SALT','CREATE','DESTROY'
        write(unit=fid,fmt='(6X,8(1X,A7))')(freg(ireg),ireg=izlim,ibulk),(freg(ireg),ireg=izlim,ibulk)
        write(unit=fid,fmt='(72("-"))')
        a3k1v6: do igc=1,nsalt()
          ispec=isalt(igc)
          mratio=0
          z8c8z6: do ireg=izlim,ibulk
            k6p6s6: if (.not.dfeq(ratreg(ireg,ispec),0.D0)) then
              if (acreat(1,ireg,igc).ne.0) mratio(ireg)=acreat(2,ireg,igc)/acreat(1,ireg,igc)
              if (adest(1,ireg,igc).ne.0) mratio(nrgnmx+ireg)=adest(2,ireg,igc)/adest(1,ireg,igc)
            endif k6p6s6
          enddo z8c8z6
          write(unit=fid,fmt='(1X,A4,1X,8(1X,F7.4))')fsalt(igc),(mratio(ireg),ireg=1,2*nrgnmx)
        enddo a3k1v6
      endif
    endif d0o2i7

    x5t6m6: if (.not.dfeq(ratind(),0.0D0)) then
      if (is_bulk()) then
        write(unit=fid,fmt='(36("-"))')
        write(unit=fid,fmt='(1X,A4,1X,2(1X,A7))')'ION','CREATE','DESTROY'
        write(unit=fid,fmt='(36("-"))')
        x4f6u1: do ispec=idxcl(),nspec()
          mratio=0
          if (acrea1(1,ibulk,ispec).ne.0) mratio(1)=acrea1(2,ibulk,ispec)/acrea1(1,ibulk,ispec)
          if (adest1(1,ibulk,ispec).ne.0) mratio(2)=adest1(2,ibulk,ispec)/adest1(1,ibulk,ispec)
          write(unit=fid,fmt='(1X,A4,1X,2(1X,F7.4))')fspc(ispec),mratio(1),mratio(2)
        enddo x4f6u1
      else
        write(unit=fid,fmt='(72("-"))')
        write(unit=fid,fmt='(1X,A4,1X,2(1X,A30,1X))')'ION','CREATE','DESTROY'
        write(unit=fid,fmt='(6X,8(1X,A7))')(freg(ireg),ireg=izlim,ibulk),(freg(ireg),ireg=izlim,ibulk)
        write(unit=fid,fmt='(72("-"))')
        x4f6u6: do ispec=idxcl(),nspec()
          mratio=0
          x8c8b6: do ireg=izlim,ibulk
            k7c6k8: if (.not.dfeq(ratreg(ireg,ispec),0.D0))then
              if (acrea1(1,ireg,ispec).ne.0) mratio(ireg)=acrea1(2,ireg,ispec)/acrea1(1,ireg,ispec)
              if (adest1(1,ireg,ispec).ne.0) mratio(nrgnmx+ireg)=adest1(2,ireg,ispec)/adest1(1,ireg,ispec)
            endif k7c6k8
          enddo x8c8b6
          write(unit=fid,fmt='(1X,A4,1X,8(1X,F7.4))')fspc(ispec),(mratio(ireg),ireg=1,2*nrgnmx)
        enddo x4f6u6
      endif
    endif x5t6m6
    if (is_bulk()) then
      write(unit=fid,fmt='(36("-"))')
      ! summary of MOVE trials
      mratio(1)=sum(amove(1,idxcl():nspec()))
      mratio(2)=sum(amove(2,idxcl():nspec()))
      mratio(3)=sum(ajump(1,idxcl():nspec()))
      mratio(4)=sum(ajump(2,idxcl():nspec()))
      if (dfeq(mratio(1),0.D0)) then
        mratio(2)=0.D0
      else
        mratio(2)=mratio(2)/mratio(1)
      endif
      if (dfeq(mratio(3),0.D0)) then
        mratio(4)=0.D0
      else
        mratio(4)=mratio(4)/mratio(3)
      endif
      write(unit=fid,fmt='(3(1X,A9,1X,"|"))')"MOVES","MOVE","JUMP"
      write(unit=fid,fmt='((1X,A9,1X,"|"),2(1X,F9.6,1X,"|"))')"ACC",mratio(2),mratio(4)
      mratio(2)=mratio(1)+mratio(3)
      if (dfeq(mratio(2),0.D0)) mratio(2)=1.D0
      write(unit=fid,fmt='((1X,A9,1X,"|"),2(1X,F9.6,1X,"|"))')"TRIALS",mratio(1)/mratio(2),mratio(3)/mratio(2)
      ! summary of STEP type trials
      write(unit=fid,fmt='(4(1X,A9,1X,"|"))')"ITER TYPE","MOVE","INDIV","SALT"
      mratio(1)=mratio(2)
      mratio(2)=sum(amove(2,idxcl():nspec()))+sum(ajump(2,idxcl():nspec()))
      mratio(3)=sum(acrea1(1,4,idxcl():nspec()))+sum(adest1(1,4,idxcl():nspec()))
      mratio(4)=sum(acrea1(2,4,idxcl():nspec()))+sum(adest1(2,4,idxcl():nspec()))
      mratio(5)=sum(acreat(1,4,1:nsalt()))+sum(adest(1,4,1:nsalt()))
      mratio(6)=sum(acreat(2,4,1:nsalt()))+sum(adest(2,4,1:nsalt()))
      do ispec=1,5,2
        if (dfeq(mratio(ispec),0.D0)) then
          mratio(ispec+1)=0.D0
        else
          mratio(ispec+1)=mratio(ispec+1)/mratio(ispec)
        endif
      enddo
      write(unit=fid,fmt='((1X,A9,1X,"|"),3(1X,F9.6,1X,"|"))')"ACC",mratio(2),mratio(4),mratio(6)
      if (dfeq(ratslt(),0.D0).and.dfeq(ratind(),0.D0)) return
      mratio(2)=sum(mratio(1:5:2))
      if (dfeq(mratio(2),0.D0)) mratio(2)=1
      write(unit=fid,fmt='((1X,A9,1X,"|"),3(1X,F9.6,1X,"|"))')"TRIALS",mratio(1)/mratio(2),mratio(3)/mratio(2),mratio(5)/mratio(2)
      write(unit=fid,fmt='(72("-"))')
    else
      write(unit=fid,fmt='(72("-"))')
      ! summary of MOVE trials
      mratio(1)=sum(amove(1,1:nspec()))
      mratio(2)=sum(amove(2,1:nspec()))
      mratio(3)=sum(ajump(1,1:nspec()))
      mratio(4)=sum(ajump(2,1:nspec()))
      mratio(5)=sum(ajin(1,1:nspec()))
      mratio(6)=sum(ajin(2,1:nspec()))
      mratio(7)=sum(ajout(1,1:nspec()))
      mratio(8)=sum(ajout(2,1:nspec()))
      do ispec=1,7,2
        if (dfeq(mratio(ispec),0.D0)) then
          mratio(ispec+1)=0.D0
        else
          mratio(ispec+1)=mratio(ispec+1)/mratio(ispec)
        endif
      enddo
      write(unit=fid,fmt='(5(1X,A9,1X,"|"),1X,A10,"|")')"MOVES","MOVE","JUMP","JUMP-IN","JUMP-OUT","SUM INOUT"
      write(unit=fid,fmt='((1X,A9,1X,"|"),5(1X,F9.6,1X,"|"))')"ACC",mratio(2),mratio(4),mratio(6),mratio(8)&
           &     ,(mratio(6)*mratio(5)+mratio(8)*mratio(7))/(mratio(5)+mratio(7)+0.0001)
      mratio(2)=sum(mratio(1:7:2))
      if (dfeq(mratio(2),0.D0)) mratio(2)=1.D0
      write(unit=fid,fmt='((1X,A9,1X,"|"),5(1X,F9.6,1X,"|"))')"TRIALS",mratio(1)/mratio(2),mratio(3)/mratio(2)&
           &   ,mratio(5)/mratio(2),mratio(7)/mratio(2),(mratio(7)+mratio(5))/mratio(2)
      ! summary of STEP type trials
      write(unit=fid,fmt='(4(1X,A9,1X,"|"))')"ITER TYPE","MOVE","INDIV","SALT"
      mratio(1)=mratio(2)
      mratio(2)=sum(amove(2,1:nspec()))+sum(ajump(2,1:nspec()))+sum(ajin(2,1:nspec()))+sum(ajout(2,1:nspec()))
      mratio(3)=sum(acrea1(1,1:4,1:nspec()))+sum(adest1(1,1:4,1:nspec()))
      mratio(4)=sum(acrea1(2,1:4,1:nspec()))+sum(adest1(2,1:4,1:nspec()))
      mratio(5)=sum(acreat(1,1:4,1:nsalt()))+sum(adest(1,1:4,1:nsalt()))
      mratio(6)=sum(acreat(2,1:4,1:nsalt()))+sum(adest(2,1:4,1:nsalt()))
      do ispec=1,5,2
        if (dfeq(mratio(ispec),0.D0)) then
          mratio(ispec+1)=0.D0
        else
          mratio(ispec+1)=mratio(ispec+1)/mratio(ispec)
        endif
      enddo
      write(unit=fid,fmt='((1X,A9,1X,"|"),3(1X,F9.6,1X,"|"))')"ACC",mratio(2),mratio(4),mratio(6)
      if (dfeq(ratslt(),0.D0).and.dfeq(ratind(),0.D0)) return
      mratio(2)=sum(mratio(1:5:2))
      if (dfeq(mratio(2),0.D0)) mratio(2)=1
      write(unit=fid,fmt='((1X,A9,1X,"|"),3(1X,F9.6,1X,"|"))')"TRIALS",mratio(1)/mratio(2),mratio(3)/mratio(2),mratio(5)/mratio(2)
      ! summary of ADD/DEL trials by region
      mratio(1)=sum(acrea1(1,1,1:nspec()))+sum(adest1(1,1,1:nspec()))+sum(acreat(1,1,1:nsalt()))+sum(adest(1,1,1:nsalt()))
      mratio(3)=sum(acrea1(1,2,1:nspec()))+sum(adest1(1,2,1:nspec()))+sum(acreat(1,2,1:nsalt()))+sum(adest(1,2,1:nsalt()))
      mratio(5)=sum(acrea1(1,3,1:nspec()))+sum(adest1(1,3,1:nspec()))+sum(acreat(1,3,1:nsalt()))+sum(adest(1,3,1:nsalt()))
      mratio(7)=sum(acrea1(1,4,1:nspec()))+sum(adest1(1,4,1:nspec()))+sum(acreat(1,4,1:nsalt()))+sum(adest(1,4,1:nsalt()))
  
      mratio(2)=sum(acrea1(2,1,1:nspec()))+sum(adest1(2,1,1:nspec()))+sum(acreat(2,1,1:nsalt()))+sum(adest(2,1,1:nsalt()))
      mratio(4)=sum(acrea1(2,2,1:nspec()))+sum(adest1(2,2,1:nspec()))+sum(acreat(2,2,1:nsalt()))+sum(adest(2,2,1:nsalt()))
      mratio(6)=sum(acrea1(2,3,1:nspec()))+sum(adest1(2,3,1:nspec()))+sum(acreat(2,3,1:nsalt()))+sum(adest(2,3,1:nsalt()))
      mratio(8)=sum(acrea1(2,4,1:nspec()))+sum(adest1(2,4,1:nspec()))+sum(acreat(2,4,1:nsalt()))+sum(adest(2,4,1:nsalt()))
      do ispec=1,7,2
        if (dfeq(mratio(ispec),0.D0)) then
          mratio(ispec+1)=0.D0
        else
          mratio(ispec+1)=mratio(ispec+1)/mratio(ispec)
        endif
      enddo
      write(unit=fid,fmt='(5(1X,A9,1X,"|"))')"SALT REGN",freg(1),freg(2),freg(3),freg(4)
      write(unit=fid,fmt='((1X,A9,1X,"|"),4(1X,F9.6,1X,"|"))')"ACC",mratio(2),mratio(4),mratio(6),mratio(8)
      mratio(2)=sum(mratio(1:7:2))
      if (dfeq(mratio(2),0.D0)) then
        write(unit=fid,fmt='((1X,A9,1X,"|"),4(1X,F9.6,1X,"|"))')"TRIALS",0.D0,0.D0,0.D0,0.D0
      else
        write(unit=fid,fmt='((1X,A9,1X,"|"),4(1X,F9.6,1X,"|"))')"TRIALS",mratio(1)/mratio(2),mratio(3)/mratio(2)&
             &   ,mratio(5)/mratio(2),mratio(7)/mratio(2)
      endif
  
      write(unit=fid,fmt='(72("-"))')
    endif
  end subroutine accrat

  ! --------------------------------------------------
  ! Calculate statistics for the energy factors used to determine
  ! acceptance for addition, move and creation steps.
  !
  ! The information output here can be used to investigate the
  ! energetics of add/deletion of particles, something that may
  ! give an indication of problems with the chemical potentials.
  subroutine avergi(ispec)
    use spec
    use strngs
    implicit none
    integer, intent(in) :: ispec
    character(15) :: mavstr, cavstr, davstr
    character(15) :: mvarstr, cvarstr, dvarstr
    call str (rs_mean(achgeng(ispec)), cavstr)
    call str (rs_variance(achgeng(ispec)), cvarstr)
    call str (rs_mean(adeleng(ispec)), davstr)
    call str (rs_variance(adeleng(ispec)), dvarstr)
    call str (rs_mean(amoveng(ispec)), mavstr)
    call str (rs_variance(amoveng(ispec)), mvarstr)
    write(unit=fidlog,fmt='(1X,A6," | ",A15,": ",A15)')fspc(ispec),mavstr,mvarstr
    write(unit=fidlog,fmt='(1X,A6," | ",A15,": ",A15)')" ",cavstr,cvarstr
    write(unit=fidlog,fmt='(1X,A6," + ",A15,": ",A15)')" ",davstr,dvarstr
  end subroutine avergi

  ! --------------------------------------------------
  ! Calculate statistics for the running energe.
  subroutine avergeng
    use spec
    use strngs
    implicit none
    double precision :: df_mean
    character(10) :: df_mean_str, df_var_str, df_total_str
    df_mean = rs_mean(adeltaf)
    call str(df_mean, df_mean_str)
    call str(df_mean*dble(rs_count(adeltaf)), df_total_str)
    call str(rs_variance(adeltaf), df_var_str)
    write(unit=fidlog,fmt='(1X,"DeltaE ",3(" | ",A6,": ",A10))')"TOTAL",&
         df_total_str,"MEAN",df_mean_str,"VAR.",df_var_str
  end subroutine avergeng

  ! --------------------------------------------------
  ! Get log of ratio try_create/try_destroy
  !
  ! Get the ratio of create/destroy attempts for free
  ! species.
  !
  ! One of the equilibration methods limits the create/destroy
  ! attempts by repeating only create or only destroy attempts
  ! when a GCMC move is selected until the attempt succeeds, at
  ! which time the reverse operation is repeated.  This ensures
  ! that the number of particles of each specie will remain close
  ! to the starting value.  If the chemical potential is exactly
  ! correct then the number of attempts to create a particle should
  ! be equal to the attempts to destroy one.  Thus the value
  ! returned here can be used to self-consistently adjust the
  ! chemical potentials.
  subroutine accget(ispec,mu_a,mu_d,p_a,p_d,accv)
    use spec
    use const
    use rnstat
    implicit none
    integer, intent(in) :: ispec
    ! add/delete energies
    double precision, intent(out) :: mu_a, mu_d
    ! Probabilities
    double precision, intent(out) :: p_a, p_d
    ! Accessible volume fraction
    double precision, intent(out) :: accv

    double precision, save, dimension(nspcmx) :: noverlap=0
    double precision, save, dimension(nspcmx) :: ncret_denom=0, ndest_denom=0
    double precision, save, dimension(nspcmx) :: ncret_numer=0, ndest_numer=0
    ! Probabilty numerators and denominators
    double precision :: nover, cret_denom, dest_denom, cret_numer, dest_numer
    if (dbc) then
      if ((ispec.lt.idxcl().or.ispec.gt.nspec())) stop "Error: specie index out of range"
    endif
    nover = noverlap(ispec)
    noverlap(ispec)=aovrlp(ispec)
    nover = noverlap(ispec) - nover

    cret_denom = ncret_denom(ispec)
    ncret_denom(ispec) = sum(acrea1(1,1:4,ispec))
    cret_denom = ncret_denom(ispec) - cret_denom

    dest_denom = ndest_denom(ispec)
    ndest_denom(ispec) = sum(adest1(1,1:4,ispec))
    dest_denom = ndest_denom(ispec) - dest_denom

    cret_numer = ncret_numer(ispec)
    ncret_numer(ispec) = sum(acrea1(2,1:4,ispec))
    cret_numer = ncret_numer(ispec) - cret_numer

    dest_numer = ndest_numer(ispec)
    ndest_numer(ispec) = sum(adest1(2,1:4,ispec))
    dest_numer = ndest_numer(ispec) - dest_numer

    ! Average energy for addition
    ! avg = -log(exp(-deltu-conc-mu)/conc)
    !
    mu_a = rs_mean(achgeng(ispec))
    call rs_reset(achgeng(ispec))
 
    ! avg = -log(exp(-deltu+conc+mu)*propscal)
    !
    mu_d = rs_mean(adeleng(ispec))
    call rs_reset(adeleng(ispec))
   
    if (dfeq(cret_denom,0.D0).or.dfeq(dest_denom,0.D0)) then
      ! stop "No trials!"
      ! for very low concentration species on of these may
      ! be zero so we just return zero which will make no changes
      p_a=0.D0
      p_d=0.D0
      accv=1.0D0
    else
      ! Calculate probabilities
      p_a = cret_numer/cret_denom
      p_d = dest_numer/dest_denom
      write(unit=fidlog,fmt=*)"#########       P_dest          ##########",p_d
      write(unit=fidlog,fmt=*)"#########       P_cret          ##########",p_a
      if (dfeq(p_a,0.D0).or.dfeq(p_d,0.D0)) then
        ! stop "No trials!"
        ! for very low concentration species on of these may
        ! be zero so we just return zero which will make no changes
        p_a=dest_denom
        p_d=cret_denom
      endif
      ! Accessible volume fraction is
      !  (attempts - overlap)/attempts
      accv = (cret_denom - nover)/cret_denom
    endif
  end subroutine accget

  ! --------------------------------------------------
  ! Save energy data
  !
  ! Update the internal data sets that monitor the energy
  ! changes (printed in 'avergi')
  subroutine avgset(bfcr)
    implicit none
    double precision, intent(in) :: bfcr
    integer :: i
    double precision :: energy
    ! save the energy changes for the ensemble
    energy = -dlog(bfcr)
    call rs_push(adeltaf,energy)
    ! save changes per move type
    if (isamov()) then
      call rs_push(amoveng(ispcnw(1)),energy)
    else if (isadd()) then
      do i=1,nchgnw
        call rs_push(achgeng(ispcnw(i)),energy)
      enddo
    else if (isdele()) then
      do i=1,nchgnw
        call rs_push(adeleng(ispcnw(i)),energy)
      enddo
    endif
  end subroutine avgset

  ! ----------------------------------------------------------------
  ! Make a change for a particle in the bulk simulation
  subroutine chgblk(ovrlap)
    use conf
    use geom
    use patch
    use spec
    implicit none
    logical, intent(out) :: ovrlap
    ! Local variables
    double precision :: drmaxi ! update radii
    integer :: ispec   ! specie data
    integer :: ii      ! particle indices
    logical :: update  ! true if r2 needs to change
    ovrlap=.false.
    if (debug) write(unit=fidlog,fmt=*)"! chgblk called with state: ",istate
    ! ADD particles INIT
    if (isadd()) then
      if (nactv+1.ge.confsz()) stop 'maximum number of atoms in box'
      ispec=ispcnw(1)
      propscal=volblk()/(ni(ispec)+1)
      ! Set Chemical potentials (non-zero for add/delete)
      if (nchgnw.ne.1) stop 'salt GC moves not implemented for bulk sims'
      if (issalt()) stop 'salt GC moves not implemented for bulk sims'
      !  propexpo = chemps(igcnw)
      !else
        propexpo = chempi(ispcnw(1))
      !endif
      call cubmov(rxnw(1), rynw(1), rznw(1), r2nw(1), lenblk())
      !! EPS: can assume all particles are in water phase in bulk
    endif
    ! DELETE particles INIT
    if (isdele()) then
      if (nactv-1.le.0) stop 'no more particles in box!'
      ispec=ispcnw(1)
      if (ni(ispec).lt.1) then
        ! no more of a specie in box
        ovrlap=.true.
        return
      endif
      propscal=(ni(ispec))/volblk()
      ! Set Chemical potentials (non-zero for add/delete)
      if (nchgnw.ne.1) stop 'salt GC moves not implemented for bulk sims'
      if (issalt()) stop 'salt GC moves not implemented for bulk sims'
      ! then
      !  propexpo = -chemps(igcnw)
      ! else
        propexpo = -chempi(ispcnw(1))
      ! endif
      indxnw(1)=getnth(ispcbk,ispec,iranff(ni(ispec)),nactv)
      rxnw(1) = rx(indxnw(1))
      rynw(1) = ry(indxnw(1))
      rznw(1) = rz(indxnw(1))
      r2nw(1) = r2(indxnw(1))
      !! EPS: can assume all particles are in water phase in bulk
    endif
    ! MOVE particles INIT
    if (isamov()) then
      ! get ion
      ispec=ispcnw(1)
      ii = indxnw(1)
      if (ismove()) then
        drmaxi=drmax(2)*xri(ispec)*2
        call cubmov(rxnw(1), rynw(1), rznw(1), r2nw(1), rx(ii), ry(ii), rz(ii), -drmaxi, drmaxi)
        update=.false.
        if (rxnw(1).gt.lenblk()) then
          rxnw(1)=rxnw(1)-lenblk()
          update=.true.
        endif
        if (rxnw(1).lt.0.D0) then
          rxnw(1)=lenblk()+rxnw(1)
          update=.true.
        endif
        if (rynw(1).gt.lenblk()) then
          rynw(1)=rynw(1)-lenblk()
          update=.true.
        endif
        if (rynw(1).lt.0.D0) then
          rynw(1)=lenblk()+rynw(1)
          update=.true.
        endif
        if (update) r2nw(1) = dsqrt(sqr(rxnw(1))+sqr(rynw(1)))
        if (rznw(1).gt.lenblk()) rznw(1)=rznw(1)-lenblk()
        if (rznw(1).lt.0.D0) rznw(1)=rznw(1)+lenblk()
        propscal=1.D0
        propexpo=0.D0
      else if (isjump()) then
        call cubmov(rxnw(1), rynw(1), rznw(1), r2nw(1), lenblk())
        propscal=1.D0
        propexpo=0.D0
      else if (isswap()) then
        rxnw(1) = rx(indxnw(1))
        rynw(1) = ry(indxnw(1))
        rznw(1) = rz(indxnw(1))
        r2nw(1) = r2(indxnw(1))
        propexpo = chempi(ispcnw(2)) - chempi(ispec)
        propscal = ni(ispec)/(ni(ispcnw(2))+1)
      endif
      !! EPS: can assume all particles are in water phase in bulk
    endif
  end subroutine chgblk

  ! ----------------------------------------------------------------
  ! Make a change for a particle in the channel simulation
  !
  !  The propscal is the number-volume component of the Metropolis
  !  acceptance probability.
  !  
  !  * (addition) is V/N+1, for salts we use (Va/(Na+1) . (Vb/(Nb+1))^z)
  !
  !  * (deletion) is N/V, for salts we use (Na/Va . (Nb/Vb)^z)
  !
  !      where z is the charge of the anion 'a', 'b' is chloride
  !
  subroutine chgsim(ovrlap)
    use conf
    use geom
    use patch
    use spec
    implicit none
    logical, intent(out) :: ovrlap

    ! Local variables
    double precision drmaxi ! update radii
    double precision epsi ! eps value
    integer   gzbin ! gz bin for Widom trials
    integer   icount   ! salt data
    integer   ispec   ! specie data
    integer ii      ! particle indices

    ! Update stats and set up move state
    ovrlap=.false.
    if (debug) write(unit=fidlog,fmt=*)"! chgsim called with state: ",istate
    propscal=1.0D0

    ! ADD particles INIT
    if (isadd()) then
      if (nactv+nchgnw.ge.confsz()) stop 'maximum number of atoms in box'
      ! Set Chemical potentials (non-zero for add/delete)
      if (issalt()) then
        propexpo = chemps(igcnw)
      else
        propexpo = chempi(ispcnw(1))
      endif
      do icount=1,nchgnw
        ispec=ispcnw(icount)
        ! Use icount /= 1 to select 'chloride' in salts and not in individual
        if (icount.ne.1) then
          ! Volumes used for Cl in salt is always the entire bulk
          ! because we always add salt Cl into the entire box.
          propscal=propscal*vreg(ibulk,idxcl())/(ni(idxcl())+icount-1)
        else
          ! ispec != counter ion
          if (iregnw.eq.ibulk) then
            propscal=propscal*vreg(iregnw,ispec)/(ni(ispec)+1)
          else
            propscal=propscal*vreg(iregnw,ispec)/(nin(iregnw,ispec)+1)
          endif
        endif
      enddo

      do icount=1,nchgnw
        ispec=ispcnw(icount)
        if (icount.ne.1) then
          call jmpmov(rznw(icount), r2nw(icount), -zreg(ibulk,ispec), zreg(ibulk,ispec), rreg(ibulk,ispec))
        else
          call jmpmov(rznw(icount), r2nw(icount), -zreg(iregnw,ispec), zreg(iregnw,ispec), rreg(iregnw,ispec))
        endif
        ovrlap=.true.
        call wall(ispec,rznw(icount),r2nw(icount),ovrlap)
        if (ovrlap) return
      enddo

      do icount=1,nchgnw
        ispec=ispcnw(icount)
        epsi=0
        call jmpfin(rxnw(icount),rynw(icount),r2nw(icount))
      enddo

      !   Update data for predicting the chemical potential
      if (calwid) then
        if (isadd()) then
          if (.not.issalt()) then
            ! get bin and increment counters
            ispec=ispcnw(1)
            gzbin=gz_bin(rznw(1))
            widobj%trys=widobj%trys+1
            widobj%trials(gzbin,ispec)=widobj%trials(gzbin,ispec)+1
          endif
        endif
      endif
    endif
    ! DELETE particles INIT
    if (isdele()) then
      if (nactv-nchgnw.le.0) stop 'no more particles in box!'
      ! Set Chemical potentials (non-zero for add/delete)
      if (issalt()) then
        propexpo = -chemps(igcnw)
      else
        propexpo = -chempi(ispcnw(1))
      endif
      do icount=1,nchgnw
        ispec=ispcnw(icount)
        if (ni(ispec).lt.1) then
          ! no more of a specie in box
          ovrlap=.true.
          return
        endif
        ! Use icount /= 1 to select 'chloride' in salts and not in individual
        if (icount.ne.1) then
          ! Volumes used for Cl in salt is always the entire bulk
          ! because we always add/delete salt Cl into the entire box.
          propscal=propscal*ni(idxcl())/vreg(ibulk,idxcl())
          if (ni(idxcl()).lt.1) then
            ! no more chloride in box
            ovrlap=.true.
            return
          endif
          indxnw(icount)=getnth(ispcbk,idxcl(),iranff(ni(idxcl())),nactv)
        else
          if (iregnw.eq.ibulk) then
            propscal=propscal*ni(ispec)/vreg(iregnw,ispec)
            indxnw(icount)=getnth(ispcbk,ispec,iranff(ni(ispec)),nactv)
          else
            ! Having no ions of a specie in the channel
            ! is not an error.
            if (nin(iregnw,ispec).lt.1) then
              ovrlap=.true.
              return
            endif
            propscal=propscal*nin(iregnw,ispec)/vreg(iregnw,ispec)
            ! indreg gets particles by region
            indxnw(icount)=indreg(iranff(nin(iregnw,ispec)),iregnw,ispec)
          endif
        endif
        if (dbc) then
          if ((indxnw(icount).le.0.or.indxnw(icount).gt.nactv)) &
               & stop "Error: [delete] could not find nth particle (nin)"
          if (ispcbk(indxnw(icount)).ne.ispec) &
               & stop "Error: found particle is not of right specie (ni)"
        endif
      enddo
      ! Check if we have duplicate anions if nchgnw > 2
      if (nchgnw.gt.2.and.indxnw(2).eq.indxnw(3)) then
        ovrlap=.true.
        return
      endif
      if (nchgnw.gt.3.and.(indxnw(2).eq.indxnw(4).or.indxnw(3).eq.indxnw(4))) then
        ovrlap=.true.
        return
      endif
      do icount=1,nchgnw
        rxnw(icount) = rx(indxnw(icount))
        rynw(icount) = ry(indxnw(icount))
        rznw(icount) = rz(indxnw(icount))
        r2nw(icount) = r2(indxnw(icount))
     enddo
    endif
    ! MOVE particles INIT
    if (isamov()) then
      ! choose ion
      ispec=ispcnw(1)
      ii = indxnw(1)
      if (ismove()) then
        if (abs(rz(ii)).le.zreg(ichan,ispec)) then
          drmaxi=drmax(1)*xri(ispec)
        else
          drmaxi=drmax(2)*xri(ispec)
        endif
        call cubmov(rxnw(1), rynw(1), rznw(1), r2nw(1), rx(ii), ry(ii), rz(ii), -drmaxi, drmaxi)
        propscal=1.D0

        ! handle mobile structural ions 
        if (localized(ispec)) then
          call mobchk(ii, rxnw(1), rynw(1), rznw(1), ovrlap)
          ! if relaxmob is active, check previous position.
          if (ovrlap.and.relaxmob) then
            call mobchk(ii, rx(ii), ry(ii), rz(ii), ovrlap)
            ovrlap = .not.ovrlap
          endif
          if (ovrlap) return
        endif
        ! Set Chemical potentials (non-zero for add/delete)
        propexpo = 0.D0
      else if (isjump()) then
        if (chonly(ispec)) then
          call jmpmov(rznw(1), r2nw(1), -zreg(izlim,ispec), zreg(izlim,ispec), rreg(izlim,ispec))
        else
          call jmpmov(rznw(1), r2nw(1), -zreg(ibulk,ispec), zreg(ibulk,ispec), rreg(ibulk,ispec))
        endif
        propscal=1.D0
        ! Set Chemical potentials (non-zero for add/delete)
        propexpo = 0.D0
      else if (isjin()) then
        call jmpmov(rznw(1), r2nw(1), -zreg(ichan,ispec), zreg(ichan,ispec), rreg(ichan,ispec))
        propscal=vin(ispec)*(ni(ispec)-nin(ichan,ispec))/(vout(ispec)*(nin(ichan,ispec)+1))
        ! Set Chemical potentials (non-zero for add/delete)
        propexpo = 0.D0
      else if (isjout()) then
        call jmpmov(rznw(1), r2nw(1), -zreg(ibulk,ispec), zreg(ibulk,ispec),&
             & -zreg(ichan,ispec), zreg(ichan,ispec), rreg(ibulk,ispec))
        propscal=vout(ispec)*nin(ichan,ispec)/(vin(ispec)*(ni(ispec)-nin(ichan,ispec)+1))
        ! Set Chemical potentials (non-zero for add/delete)
        propexpo = 0.D0
      else if (isswap()) then
        propexpo = chempi(ispcnw(2)) - chempi(ispec)
        if (iregnw.eq.ibulk) then
          propscal=propscal*ni(ispec)/vreg(iregnw,ispec)
          propscal=propscal*vreg(iregnw,ispcnw(2))/(ni(ispcnw(2))+1)
        else
          propscal=propscal*nin(iregnw,ispec)/vreg(iregnw,ispec)
          propscal=propscal*vreg(iregnw,ispcnw(2))/(nin(iregnw,ispcnw(2))+1)
        endif
        rxnw(1) = rx(ii)
        rynw(1) = ry(ii)
        rznw(1) = rz(ii)
        r2nw(1) = r2(ii)
        ! get swap specie from second entry
        if (xri(ispec).gt.xri(ispcnw(2))) then
          ! new specie is small then old specie, no overlap possible
          ! so skip wall test
          return
        endif
        ! set ispec to ispcnw(2) for wall test
        ispec = ispcnw(2)
      endif

      ovrlap=.true.
      call wall(ispec, rznw(1), r2nw(1), ovrlap)
      if (.not.ovrlap) then
        if (isajmp()) then
          call jmpfin(rxnw(1), rynw(1), r2nw(1))
        end if
      end if
    end if
  end subroutine chgsim

  ! -------------------------------------------------------------
  ! Commit changes after a MCGC move has been accepted
  !
  subroutine commit
    use const
    use conf
    use geom
    use patch
    use spec
    use simstate
    implicit none
    ! LOCALS
    logical :: isbulk
    double precision :: qi
    integer :: ichg  ! changed particle indexs
    integer :: ii,jj ! particle indices
    integer :: jspec
    isbulk = is_bulk()

    g5t3r5: if (isadd().or.isamov()) then
      ! ADD or MOVE
      w0w3t1: do ichg=1,nchgnw
        ! set region before updating position
        c0d1e8: if (isadd()) then
          ! Get new II from conf
          indxnw(ichg)=idxget(ispcnw(ichg))
          ii=indxnw(ichg)
          !--------------------
          ! ADD
          if (.not.isbulk) call addreg(ispcnw(ichg),ii,rznw(ichg))
          ! Iff ichg not 1 then we need to merge riicrn into riinw
          if (ichg.ne.1) then
            do jj=1,ichg-1
              riinw(indxnw(jj),ichg)=riicrn(jj,ichg)
            enddo
          endif
          charge = charge + xz(ispcnw(ichg))
        else c0d1e8
          ! Get II from indxnw
          ii=indxnw(ichg)
          if (isswap()) then
            if (.not.isbulk) call delreg(ispcnw(1),ii,rznw(ichg))
            if (.not.isbulk) call addreg(ispcnw(2),ii,rznw(ichg))
            ispcbk(ii)=ispcnw(2)
          else
            !--------------------
            ! MOVE
            if (.not.isbulk) call setreg(ispcnw(ichg),ii,rz(ii),rznw(ichg))
          endif
        endif c0d1e8
        ! Do not need to update position for swap
        if (.not.isswap()) then
          ! update position
          call setpos(ii,rxnw(ichg),rynw(ichg),rznw(ichg),r2nw(ichg))
          ! update rqqii matrix
          qi=xq(ispcnw(ichg))
          do jj=1,nactv
            if (jj.eq.ii) cycle
            jspec=ispcbk(jj)
            if (jspec.eq.0) cycle
            if (dbc) then
              if (dfeq(riinw(jj,ichg),0.D0)) stop "Rii_nw was zero."
            endif
            if (.not.(dfeq(xq(ispcnw(ichg)),0.D0).or.dfeq(xq(jspec),0.D0))) then
              rqqii(jj,ii)=xq(jspec)*qi/riinw(jj,ichg)
            else
              rqqii(jj,ii)=riinw(jj,ichg)
            endif
            rqqii(ii,jj)=rqqii(jj,ii)
          end do
          rqqii(ii,ii)=0
          ! update rip matrix
          c3z8a6: if (.not.isbulk.and..not.is_homogeneous()) then
            rip(1:npatch,ii)=ripnw(1:npatch,ichg)
            iprip(1:npatch,ii)=parea(1:npatch)/rip(1:npatch,ii)
          endif c3z8a6
        endif
      enddo w0w3t1
    else g5t3r5
      !--------------------
      ! DELETE
      s1y4m4: do ichg=1,nchgnw
        if (.not.isbulk) call delreg(ispcnw(ichg),indxnw(ichg),rznw(ichg))
        call idxrel(ispcnw(ichg),indxnw(ichg))
        charge = charge - xz(ispcnw(ichg))
        ! idxrel sets ISPCBK(?) to zero to indicate invalid particle
      enddo s1y4m4
    endif g5t3r5

    r1b2a3: if (.not.isbulk.and.do_electrostatic().and..not.is_homogeneous()) then
      h(1:npatch)=hnw(1:npatch)
      c(1:npatch)=c(1:npatch)+dcnw(1:npatch)
    endif r1b2a3

  end subroutine commit

  ! --------------------------------------------------------------
  ! Echo input file content
  !
  ! This writes the interpreted and default values for
  ! module parameters that may have been read from the 
  ! input file.  It writes this information in the same
  ! format as an input file.
  subroutine ectral(fid)
    use strngs
    implicit none
    integer, intent(in) :: fid
    character(20) :: outstr
    write(unit=fid,fmt='(A)')fstry
    call str (drmax(1), outstr)
    write(unit=fid,fmt='(A,1X,A)')fsdrmi,trim(adjustl(outstr))
    call str (drmax(2), outstr)
    write(unit=fid,fmt='(A,1X,A)')fsdrmo,trim(adjustl(outstr))
    call str (rates_(1), outstr)
    write(unit=fid,fmt='(A,1X,A)')fsrtmv,trim(adjustl(outstr))
    call str (rates_(2), outstr)
    write(unit=fid,fmt='(A,1X,A)')fsrtsl,trim(adjustl(outstr))
    call str (rates_(3), outstr)
    write(unit=fid,fmt='(A,1X,A)')fsrtid,trim(adjustl(outstr))
    call str (rates_(4), outstr)
    write(unit=fid,fmt='(A,1X,A)')fsrtsw,trim(adjustl(outstr))
    call str (k_mobl, outstr)
    write(unit=fid,fmt='(A,1X,A)')fskmob,trim(adjustl(outstr))
    write(unit=fid,fmt='(A)')fsend
    write(unit=fid,fmt=*)
  end subroutine ectral

  ! ----------------------------------------------------------------------
  ! - QUERY functions


  ! --------------------------------------------------
  ! Does the 'move' involve a salt?
  logical function issalt()
    implicit none
    issalt=istate(1).eq.addslt.or.istate(1).eq.remslt
  end function issalt

  ! --------------------------------------------------
  ! Does the 'move' involve a specie swap?
  logical function isswap()
    implicit none
    isswap=istate(1).eq.0.and.istate(2).eq.swap_l
  end function isswap

  ! --------------------------------------------------
  ! Does move involve adding particles (either indiv or as salt)?
  logical function isadd()
    implicit none
    isadd=istate(1).gt.0
  end function isadd

  ! --------------------------------------------------
  ! Does move involve deleting particles (either indiv or as salt)?
  logical function isdele()
    implicit none
    isdele=istate(1).lt.0
  end function isdele

  ! --------------------------------------------------
  ! Is any one of the move types (not add or delete)
  logical function isamov()
    implicit none
    isamov=istate(1).eq.0
  end function isamov

  ! --------------------------------------------------
  ! Is the move a small displacement move
  logical function ismove()
    implicit none
    ismove=istate(1).eq.0.and.istate(2).eq.sphr_l
  end function ismove

  ! --------------------------------------------------
  ! is any one of the jump types (in region, into or
  ! out of channel)
  logical function isajmp()
    implicit none
    isajmp=istate(1).eq.0.and.istate(2).ge.jump_l
  end function isajmp

  ! --------------------------------------------------
  ! Is the move a jump within a region
  logical function isjump()
    implicit none
    isjump=istate(1).eq.0.and.istate(2).eq.jump_l
  end function isjump

  ! --------------------------------------------------
  ! Is move a jump into the channel
  logical function isjin()
    implicit none
    isjin=istate(1).eq.0.and.istate(2).eq.jmpin_
  end function isjin

  ! --------------------------------------------------
  ! Is move a jump out of the channel
  logical function isjout()
    implicit none
    isjout=istate(1).eq.0.and.istate(2).ge.jmpou_
  end function isjout

  ! --------------------------------------------------
  ! Are any internal variables on the current move
  ! NaNs? (only runs if  dbc=.true.)
  subroutine isnanstate()
    use spec
    implicit none
#include "require.h"

    integer :: idx
    if (dbc) then
      if (isnan(propscal)) then
        stop "The proportionality scale is NaN"
      elseif (isnan(propexpo)) then
        stop "The proportionality exponent is NaN"
      else
        do idx=1,nchgnw
          if (isnan(rxnw(idx))) then
            stop "A proposed x coord is NaN"
          else if (isnan(rynw(idx))) then
            stop "A proposed y coord is NaN"
          else if (isnan(rznw(idx))) then
            stop "A proposed z coord is NaN"
          else if (isnan(r2nw(idx))) then
            stop "A proposed radial value is NaN"
          else if (1.gt.ispcnw(idx).or.nspec().lt.ispcnw(idx)) then
            stop "A proposed specie value is out of range"
          endif
        enddo
      endif
    endif
  end subroutine isnanstate

  ! --------------------------------------------------
  ! Log move
  !
  ! Log a trial move (need)
  ! success id type oldx,y,z newx,y,z indxnw,ispcnw
  subroutine logmov(success)
    use conf
    logical, intent(in) :: success
    integer :: istat, indx, idx
    logical :: isopn
    character(16) :: isred
    istat=0
    inquire(unit=fidmlg,iostat=istat,action=isred,opened=isopn)
    if (istat.ne.0) then
      write(*,*)"Error: I/O error reading from trial log ID"
      stop 1
    endif
    if (.not.isopn) then
      open(unit=fidmlg,iostat=istat,file="dat/trial."//firun)
      if (istat.ne.0) then
        write(unit=fidlog,fmt=*)"Error: opening trial log file"
        stop 1
      endif
      inquire(unit=fidmlg,iostat=istat,action=isred,opened=isopn)
      if (istat.ne.0) then
        write(unit=fidlog,fmt=*)"Error: I/O error reading from trial log ID"
        stop 1
      elseif (.not.isopn) then
        stop "Error: Unit ID is unopen"
      endif
    endif
    if (.not.(isred.eq.'READ'.or.isred.eq.'READWRITE'.or.isred.eq.'read'.or.isred.eq.'readwrite')) then
      write(unit=fidlog,fmt=*)"Unit ID ",fidmlg," is opened for :",trim(isred)," need 'READ' or  'READWRITE'"
      stop 1
    endif
    do idx=1,nchgnw
      indx=indxnw(idx)
      write(unit=fidmlg,fmt="(L1,' ',I2,' ',2(I1,' '),I4,6(' ',F7.2))",iostat=istat)&
          &success,istate(1),istate(2),ispcnw(idx),indxnw(idx),rxnw(idx),rynw(idx),rznw(idx),rx(indx),ry(indx),rz(indx)
      if (istat.ne.0) then
        write(*,*)"Error writing trial log file"
        stop 1
      endif
    enddo
  end subroutine logmov

   ! ----------------------------------------------------------------
  ! UNIFIED MOVE
  ! 
  ! This method is called after setting up the move using one
  ! of the 'selXXX' methods then 'chgXXX' methods.  This
  ! method then calculates the energy change and converts this
  ! into a acceptance probability.  If the move is accepted it
  ! calls the 'commit' method to update the 'conf' module.
  !
  ! NOTE: this method always exits without calling commit 
  subroutine moveeq(riifnc,docmmt)
    use conf
    use geom
    use patch
    use spec
    use simstate
    implicit none
    ! Should we perform commit operation if move accepted?
    logical, intent(in) :: docmmt

#include "require.h"
    ! Function for calculating radial distance between two points
    interface
      double precision function riifnc(x1,y1,z1,x2,y2,z2)
        double precision, intent(in) :: x1,y1,z1,x2,y2,z2
      end function riifnc
    end interface

    ! Local variables
    double precision :: uii,uip,deltu,bfcr ! energy terms
    double precision :: umfi ! Mean field (only if calwid)
    double precision :: umob, mobtmp ! potential from mobile ion
    double precision :: proptmp ! check chemical potential in dbc mode
    logical :: ovrlap
    double precision :: rij, xq1
    integer          :: idx

    ! Verify state of move.
    if (dbc) call isnanstate()

    ! COULOMB and ICC
    ovrlap=.true.
    uip=0.D0
    umfi=0.D0
    umob=0.D0
    ! swap can ignore cergit as electrostatic uii must be 0
    if (.not.isswap().or..not.do_electrostatic()) then
      call cergit(uii,umfi,ovrlap,riifnc)
    else
      ! with swap there is no energy change but we need to check
      ! distances if the new particle is larger than the old one
      if (xri(ispcnw(2)).gt.xri(ispcnw(1)).or..not.do_electrostatic()) then
        xq1 = xq(ispcnw(2))
        !$omp parallel do private(rij) firstprivate(xq1) reduction(.or.:ovrlap)
        do idx=1,nactv
          if (indxnw(1).ne.idx.and.ispcbk(idx).ne.0) then
            rij = xq1 * xq(ispcbk(idx)) / rqqii(indxnw(1), idx)
            if (rij.le.dd_get(ispcnw(2), ispcbk(idx))) ovrlap=.true.
          end if
        end do
        !$omp end parallel do
      end if
    endif
    if (ovrlap) then
      if (is_bulk().and.isadd()) then
        aovrlp(ispcnw(1))=aovrlp(ispcnw(1))+1
      endif
      return
    endif
    ! As mobile ions can not be part of a salt, this test
    ! will only be true when moving a mobile ion.
    if (localized(ispcnw(1))) then
      ! old penalty
      call mobpen(indxnw(1),rx(indxnw(1)),ry(indxnw(1)),rz(indxnw(1)),mobtmp)
      umob = mobtmp
      ! new penalty
      call mobpen(indxnw(1),rxnw(1),rynw(1),rznw(1),mobtmp)
      umob = k_mobl * (mobtmp - umob)
    endif

    ! swap can ignore pergit as uip must be 0
    if (.not.isswap().and..not.is_bulk().and..not.is_homogeneous().and.do_electrostatic()) then
      call pergit(uip)
    endif

    ! Test Chemical potentials (non-zero for add/delete)
    if (dbc) then
      if (.not.isamov()) then
        if (issalt()) then
          proptmp = chemps(igcnw)
        else
          proptmp = chempi(ispcnw(1))
        endif
        if (isdele()) then
          proptmp = -proptmp
        endif
        if (.not.dfeq(proptmp,propexpo)) then
          write(unit=fidlog,fmt=*)"Check chemical potential ", proptmp,  &
               " is not equal to provided potential ",propexpo
          write(unit=fidlog,fmt=*)"State[", istate, "], IGC [",igcnw,"] SPC [",ispcnw(1),"]"
          stop "Chemical potential not correctly set"
        endif
      endif
    endif

    deltu=0
    if (debug) write(unit=fidlog,fmt=*)"POT: UII=[",uii,"] UIP=[",uip,"] UMOB=[",umob,"]"
    if (dbc) then
      if (isnan(uii)) then
        stop "The electrostatic potential energy component Uii is nan"
      elseif (isnan(uip)) then
        stop "The induced-charge potential energy component Uip is nan"
      elseif (isnan(uii)) then
        stop "The localisation potential energy component Umob is nan"
      endif
    endif
    if (do_electrostatic()) then
      deltu=uii+uip
    endif
    deltu=deltu+umob-propexpo

    bfcr = dexp(-deltu)*propscal

    if (bfcr.ge.ranff()) then
      ! TRIAL ACCEPTED

      ! Update statistics for a successful attempt
      call on_cmt

      ! Record acceptance probability value for accepted moves
      call avgset(bfcr)

      ! Update data for predicting the chemical potential only for
      ! accepted trials.
      if (calwid) then
        if (isadd().and..not.issalt()) then
          ! (ignore umob as mobile ions will never be part of widom trial)
          call widom_trial (rznw(1),ispcnw(1),umfi,uii+uip,uii)
        endif
      endif

      ! Call log
      if (dotlog) call logmov(.true.)

      ! COMMIT changes
      if (docmmt) call commit
    else
      ! Call log
      if (dotlog) call logmov(.false.)
    endif

  contains

    ! ----------------------------------------------------------------
    ! UNIFIED COLOUMB (UII) ENERGY ROUTINE
    !
    ! If ovrlap is true uii is undefined
    !
    ! Calculate the change in coulomb energy.
    !
    ! Dezso"" Boda, Dirk Gillespie, Wolfgang Nonner, Douglas Henderson, and Bob Eisenberg
    ! "Computing induced charges in inhomogeneous dielectric media:
    ! Application in a Monte Carlo simulation of complex ionic systems"
    ! PHYSICAL REVIEW E 69, 046702 2004 
    ! left half of eqn(32) with eqn(33)
    !
    ! Reference for canonical moves (radial,jump) and preference sampling
    ! (jump-in,jump-out):
    ! Dezso"" Boda, Douglas Henderson, and David Busath ""Monte Carlo
    ! study of the selectivity of calcium channels: improved
    ! geometrical model"", Molecular Physics, 2002, 110(14) 2361-2368
    ! 
    !   radial,jump: P = min[ 1, exp(-\delta U / kT) ]
    !
    !   preference sampling for jump-in/jump-out:  
    !   P = min [1, V_end/V_start * exp(-\delta U / kT) ]
    ! 
    ! Reference for grand-canonical insertion and deletion
    ! Dezso"" Boda, Mo'nika' Valisko', Bob Eisenberg, Wolfgang Nonner,
    ! Douglas Henderson, and Dirk Gillespie ""The effect of protein
    ! dielectric coefficient on the ionic selectivity of a calcium channel""
    ! The Journal of Chemical Physics, 2006, 125, 034901
    !
    !  P = min [ 1,  (N_Cation,start! N_Anion,start!) exp( (\nu B - \delta U)/kT ) ]
    !          [      (N_Cation,end! N_Anion,end!)                                 ]
    !  with
    !        B =\mu_Ca + 2\mu_Cl + ln(V_Ca/\lambda^3_Ca) + 2 ln (V_Cl / \lambda^3_Cl)
    !  and
    !        \lamda_i = h/(2\pi.k.T.m_i)^(1/2)
    !
    !  coded as (addition)
    !      V_cation*V_anion.../(N_Cation+1).(N_Anion+1)... exp (\mubar-\delta U)/kT
    !  (deletion)
    !      (N_Cation).(N_Anion)/V_cation*V_anion... exp (\mubar-\delta U)/kT
    !
    ! SEC1: addition of salt we account for interaction between added particles
    ! 
    !  sum(salt_i,salt_j:j > 1) q_i * q_j * {(1/eps_i + 1/eps_j)/2} / (rnw_ij)
    ! 
    ! SEC2: for moved and new particles we add interaction between particle and all
    !  other particles
    !
    !  sum(chg_i,j: j /=i)  q_i * q_j * {(1/eps_i + 1/eps_j)/2} / (rnw_ij)
    !
    ! SEC3: for moved or deleted particles subtract the interaction with all
    !  other particles (in original position 'ii')
    !
    !  sum(chg_ii, j: j /= i)  q_i * q_j * {(1/eps_i + 1/eps_j)/2} / (r_ij)
    !
    ! TOTAL: SEC1 + SEC2 - SEC3
    subroutine cergit(uii,umfi,ovrlap,riifnc)
      use conf
      use spec
      implicit none
      double precision, intent(out) :: uii ! Coulomb energy
      double precision, intent(out) :: umfi ! Mean field (only if calwid)
      logical, intent(out) :: ovrlap ! if particles overlap
      ! Function for calculating radial distance between two points
      interface
        function riifnc(x1,y1,z1,x2,y2,z2)
          double precision riifnc
          double precision, intent(in) :: x1,y1,z1,x2,y2,z2
        end function riifnc
      end interface

      ! LOCALS
      double precision :: uiinw,uiiold,uiitmp
      integer          :: i  ! generic loop index
      integer          :: icount ! multiple change loop counter
      integer          :: ispec
      uiinw=0
      uiiold=0
      ovrlap=.false.
      umfi=0
      l7z0d8: do icount=1,nchgnw
        ispec=ispcnw(icount)

        ! SEC1: Energy change from new particle + new particle
        d7o0k9: if (issalt().and.isadd().and.icount.ne.nchgnw) then
          uiitmp=0.D0
          riicrn(icount,icount)=0.D0
          a7d3j8: do i=icount+1,nchgnw
            riicrn(i,icount)=riifnc(rxnw(i),rynw(i),rznw(i),rxnw(icount),rynw(icount),rznw(icount))
            ! Check if riicrn indicates overlap
            j5l8k3: if (riicrn(i,icount).le.dd_get(ispec, ispcnw(i))) then
              ovrlap=.true.
              return
            else j5l8k3
              riicrn(icount,i)=riicrn(i,icount)
              if (.not.(dfeq(xq(ispec),0.D0).or.dfeq(xq(ispcnw(i)),0.D0))) then
                uiitmp =uiitmp + xq(ispcnw(i))*(2/epsw) &
                     &                    / riicrn(i,icount)
              endif
            endif j5l8k3
          enddo a7d3j8

          if (.not.dfeq(xq(ispec),0.D0)) then
            ! apply constant factor
            uiinw=uiinw + uiitmp * xq(ispec)/2
          endif
        endif d7o0k9

        ! SEC2: Energy change from new particle/position + existing particles
        x5l5t7: if (isamov().or.isadd()) then
          uiitmp=0.D0
          !$omp parallel do shared(riinw), firstprivate(nactv,icount), &
          !$omp & reduction(+:uiitmp) reduction(.or.:ovrlap)
          y8u5z9: do i=1,nactv
            if (ovrlap) cycle y8u5z9
            z9a7f3: if (indxnw(icount).ne.i.and.ispcbk(i).ne.0) then
              riinw(i,icount)=riifnc(rx(i),ry(i),rz(i),rxnw(icount),rynw(icount),rznw(icount))
              ! Check if riinw indicates overlap
              y6s6y7: if (riinw(i,icount).le.dd_get(ispec, ispcbk(i))) then
                ovrlap=.true.
              endif y6s6y7
              if (.not.(dfeq(xq(ispec),0.D0).or.dfeq(xq(ispcbk(i)),0.D0))) then
                uiitmp =uiitmp + xq(ispcbk(i))*(2/epsw)/riinw(i,icount)
                ! calculate mean field
                if (calwid) then
                  if (isadd().and..not.issalt()) umfi=umfi + xq(ispec)/(epsw*riinw(i,icount))
                endif
              endif
            else z9a7f3
              riinw(i,icount)=0.D0
            endif z9a7f3
          enddo y8u5z9
          !$omp end parallel do
          if (ovrlap) return
          ! apply constant factor
          if (.not.dfeq(xq(ispec),0.D0)) then
            uiinw=uiinw + uiitmp * xq(ispec)/2
          endif
        endif x5l5t7

        ! SEC3: Energy from particle in old position
        a0e6a5: if (isamov().or.isdele()) then
          if (.not.dfeq(xq(ispec),0.D0)) then
            uiitmp=0.D0
            !$omp parallel do reduction(+:uiitmp) firstprivate(nactv,icount)
            p8l3k4: do i=1,nactv
              e9g9k7: if (indxnw(icount).ne.i.and.ispcbk(i).ne.0) then
                if (.not.dfeq(xq(ispcbk(i)),0.D0)) then
                  uiitmp =uiitmp + rqqii(i,indxnw(icount))*(2/epsw)
                endif
              endif e9g9k7
            enddo p8l3k4
            !$omp end parallel do
            uiiold=uiiold+uiitmp/2
          endif
        endif a0e6a5
      enddo l7z0d8
      uii=uiinw-uiiold
    end subroutine cergit

    ! ----------------------------------------------------------------
    ! Calculate ICC part of energy
    !
    ! Dezso"" Boda, Dirk Gillespie, Wolfgang Nonner, Douglas Henderson, and Bob Eisenberg
    ! "Computing induced charges in inhomogeneous dielectric media:
    ! Application in a Monte Carlo simulation of complex ionic systems"
    ! PHYSICAL REVIEW E 69, 046702 2004 
    ! right half of eqn(32) with eqn(33), rewritten as eqn(35) for move of particle i:
    !
    !   \D W_i = z(i).e/(8.pi) sum(p){area(p) [h_nw(p)/r_ipnw]-[h_old(p)/r_ipold]} +
    !         sum(j; j/=i){z(j).e/(8.pi) sum(p){area(p) [h_nw(p)-h_old(p)]/r_jpold]}}
    !
    !
    ! SEC1 : {all} change in energy caused by change in 'h' for particles that
    !   have not moved (lower part of eqn above)
    !
    !      sum(i,a) = q(i).area(a).hnw(a)-q(i).area(a).hold(a) ]/{2 * r(i,a)}
    !
    ! SEC2 : {add,move} contrib for moved/new particle from new 'h'
    !   
    !      sum(i,a) = q(i).area(a).hnw(a) ]/{2 * r(i_new,a)}
    !
    ! SEC3 : {move,del} contrib from old position/del particle from new 'h'
    !
    !      sum(i,a) = q(i).area(a).hold(a) ]/{2 * r(i_old,a)}
    !
    !  total = SEC1 + SEC2 - SEC3
    subroutine pergit(uip)
      use conf
      use patch
      use spec
      implicit none
      double precision, intent(out) :: uip
      double precision :: ddot
      external :: ddot
      ! locals
      integer :: i, ii, ipch
      double precision :: uiptmp
      double precision :: xqspc

      call calch()

      dh(1:npatch) = h(1:npatch) - hnw(1:npatch)
      uip    = 0.D0
      uiptmp = 0.D0
      ! Calculate \delta uip for all original position particles + patches
      ! based on the change in 'h'
      do ii=1,nactv
        if (ispcbk(ii).eq.0) cycle
        xqspc=xq(ispcbk(ii))
        if (.not.dfeq(xqspc,0.D0)) then
          !     -- add for originals
          uiptmp = uiptmp + xqspc*0.5D0*ddot(npatch,dh,1,iprip(1,ii),1)
          !! xqspc*dh(ipch)*0.5D0*iprip(ipch,ii))
        endif
      enddo
      uip = uip + uiptmp

      ! Adjust for change
      do i=1,nchgnw
        xqspc=xq(ispcnw(i))
        if (.not.dfeq(xqspc,0.D0)) then

          if (isadd().or.isamov()) then
            uiptmp=0.D0
            !   -- new configuration: from rnew calculated in calch ---
            do ipch=1,npatch
              uiptmp = uiptmp + xqspc*hnw(ipch)*parea(ipch)/(2*ripnw(ipch,i))
            enddo
            uip=uip + uiptmp
          endif
          if (isamov().or.isdele()) then
            ii=indxnw(i)
            uiptmp=0.D0
            !     tmp1 = 0
            !     daxpy(h, dh, tmp1)
            !     tmp2 = 0
            !     gemv (1/rip, tmp1, tmp2 )
            !     sum = 0
            !     gemv (parea, tmp2, sum)
            !     uiptmp = xqspc/2 * sum 

            !   -- remove old configuration: lookup table --
            do ipch=1,npatch
              !       USE (h + dh) here as we need to remove value added in the first loop.
              uiptmp = uiptmp - xqspc*(h(ipch)+dh(ipch))*0.5D0*iprip(ipch,ii)
            enddo
            uip=uip + uiptmp
          endif
        endif
      enddo
    end subroutine pergit

    ! ------------------------------------------------------------------
    !
    ! Calculate 'c' vector, then back substitute into 'A^-1' to generate 'h'
    !
    ! Dezso"" Boda, Dirk Gillespie, Wolfgang Nonner, Douglas Henderson, and Bob Eisenberg
    ! "Computing induced charges in inhomogeneous dielectric media:
    ! Application in a Monte Carlo simulation of complex ionic systems"
    ! PHYSICAL REVIEW E 69, 046702 2004 
    ! with eqn(31) and eqn(28)
    !
    !  c(a) = sum(i, p : all p) = deps(p) * [q(i)/{ eps(i) * 4pi * r_ip^2 }] * [dot(V_ip, V_pnorm) / rip]
    !
    !   where [q(i) / {eps(i) * 4pi * r_ip^2}] * [dot(V_ip, V_pnorm)/rip] is the screened electric field from i 
    !         deps(p) is function of change of eps at boundary point p {(eps_in - eps_out}/{eps_in+eps_out)}
    !         V_ip/r_ip is unit vector in direction i to p and V_pnorm is unit vector normal to surface at p
    !
    subroutine calch
      use conf
      use patch
      use spec
      implicit none

      ! LOCALS
      double precision  alfa,ri             ! particle data
      double precision rxki,ryki,rzki,rki,rij ! particle-patch data
      integer ipch       ! patch index
      integer icount

      ! ------- CALCULATE VECTOR C -----------------------------
      dcnw=0.D0
      do icount=1,nchgnw
        ri=xri(ispcnw(icount))
        if (.not.dfeq(xq(ispcnw(icount)),0.D0)) then
          alfa=xq(ispcnw(icount))/(4*pi)
          ! ----- NEW ----------------------------------------------
          if (isadd().or.isamov()) then
            u8z1n8: do ipch=1,npatch
              rxki = prx(ipch)-rxnw(icount)
              ryki = pry(ipch)-rynw(icount)
              rzki = prz(ipch)-rznw(icount)
              rij = dsqrt(sqr(rxki)+sqr(ryki)+sqr(rzki))
              rki  = rxki*pux(ipch)+ryki*puy(ipch)+rzki*puz(ipch)
              dcnw(ipch) = dcnw(ipch) - deps(ipch)*alfa*rki/(epsw*rij*rij*rij)
              ! Calculate and store distances between patches and particle ii
              ripnw(ipch,icount)=rij
              ! (stop if overlap)
              if (dbc) then
                if (dbc_level.ge.dbc_check) then
                  if (rij.le.ri) ovrlap = .true.
                end if
              end if
            enddo u8z1n8
            if (dbc) then
              if (dbc_level.ge.dbc_check) then
                if (ovrlap) stop "Overlap between article and patch should be impossible"
              end if
            end if
          endif
          ! ----- OLD ----------------------------------------------
          if (isamov().or.isdele()) then
            do ipch=1,npatch
              rij = rip(ipch,indxnw(icount))
              if (dbc) then
                if (dbc_level.ge.dbc_check) then
                  if (dfeq(rij,0.D0)) then
                    ovrlap = .true.
                    write(unit=fidlog,fmt=*)"Error: R_i_p is zero (patch, element, nactv):",&
                         &ipch,indxnw(icount),nactv,isamov(),isdele()
                  endif
                endif
              endif
              rki=(prx(ipch)-rx(indxnw(icount)))*pux(ipch) + (pry(ipch)-ry(indxnw(icount)))*puy(ipch)&
                   &  + (prz(ipch)-rz(indxnw(icount)))*puz(ipch)
              dcnw(ipch) = dcnw(ipch) + deps(ipch)*alfa*rki/(epsw*rij*rij*rij)
            enddo
            if (dbc) then
              if (dbc_level.ge.dbc_check) then
                if (ovrlap) then
                  stop 1
                endif
              endif
            endif
          endif
        end if
      enddo
      do ipch = 1,npatch
        hnw(ipch) = c(ipch) - dcnw(ipch)
      enddo
    
      ! ------- CALCULATE VECTOR H -----------------------------

      call baksub(hnw)

    end subroutine calch

  end subroutine moveeq

  ! --------------------------------------------------
  ! ratmv, ratslt and ratind
  ! Get relative rates of particle moves and salt or individual 
  ! insertions. After 'rftral', ratmv + ratslt + ratind =~= 1

  ! --------------------------------------------------
  ! Relative rate of movement
  double precision function ratmv()
    implicit none
    ratmv=rates_(1)
  end function ratmv

  ! --------------------------------------------------
  ! Relative rate of salt add/del
  double precision function ratslt()
    implicit none
    ratslt=rates_(2)
  end function ratslt

  ! --------------------------------------------------
  ! Relative rate of indiv. ion add/del
  double precision function ratind()
    implicit none
    ratind=rates_(3)
  end function ratind

  ! --------------------------------------------------
  ! Relative rate of swap type moves
  double precision function ratswap()
    implicit none
    ratswap=rates_(4)
  end function ratswap

  ! -------------------------------------------------------------
  ! Read in the trial input section
  !
  subroutine rdtral(fid,sname,svalue,istat)
    use strngs
    implicit none
    integer, intent(in) :: fid
    character(len=*), intent(in) :: sname,svalue
    integer, intent(out) :: istat
    logical, dimension(5) :: mask_
    character(32) :: nme_
    character(1024) :: val_
    rates_ = 0
    if (dbc) then
      if (sname.ne.fstry) stop "Error: incorrect section name"
      if (len_trim(svalue).ne.0) stop "Error: section does not any parameters"
    endif
    mask_=.false.
    c1w8h3: do
      val_ = " "
      call readnv(fid,nme_,val_,istat)
      ! exit on bad read or section 'end'
      if (istat.ne.0) return
      if (nme_.eq.fsend) exit c1w8h3
      ! looking for drmaxin,drmaxout,ratmv, ratslt, ratind
      b4i6i7: select case (nme_)
      case (fsdrmi) b4i6i7
        read(val_,'(D20.13)')drmax(1)
        mask_(1)=.true.
      case (fsdrmo) b4i6i7
        read(val_,'(D20.13)')drmax(2)
        mask_(2)=.true.
      case (fsrtmv) b4i6i7
        read(val_,'(D20.13)')rates_(1)
        rates_(1) = abs(rates_(1))
        mask_(3)=.true.
      case (fsrtsl) b4i6i7
        read(val_,'(D20.13)')rates_(2)
        rates_(2) = abs(rates_(2))
        mask_(4)=.true.
      case (fsrtid) b4i6i7
        read(val_,'(D20.13)')rates_(3)
        rates_(3) = abs(rates_(3))
        mask_(5)=.true.
      case (fsrtsw) b4i6i7
        read(val_,'(D20.13)')rates_(4)
        rates_(4) = abs(rates_(4))
      case (fskmob) b4i6i7
        read(val_,'(D20.13)')k_mobl
      case default b4i6i7
        call s9b8m1("Name "//nme_//" is not valid in simulation parameter (channel) section")
      end select b4i6i7
    enddo c1w8h3
    t4p3z0: if (.not.all(mask_)) then
      call s9b8m1("Not all required tags were present")
    endif t4p3z0

  contains

    subroutine s9b8m1(msg)
      implicit none
      character(len=*), intent(in) :: msg
      write(unit=fidlog,fmt=*)"Bad simulation parameter (trial) section in input."
      write(unit=fidlog,fmt=*)msg
      write(unit=fidlog,fmt=*)"Required tags are:"
      write(unit=fidlog,fmt=*)fsdrmi,fsdrmo,fsrtsl,fsrtid,fsrtmv
      write(unit=fidlog,fmt=*)"Optional tags are:"
      write(unit=fidlog,fmt=*)fskmob,fsrtsw,fsrtjp
      stop 1
    end subroutine s9b8m1

  end subroutine rdtral

  ! --------------------------------------------------
  ! Initialise trial module after reading the input file
  subroutine rftral
    use spec
    implicit none
    double precision :: s_
    write(unit=fidlog,fmt='(72("-"))')
    write(unit=fidlog,fmt=*)"Key references"
    write(unit=fidlog,fmt='(72("-"))')
    write(unit=fidlog,fmt=*)"Reference for original MC"
    write(unit=fidlog,fmt=*)"D Boda and D D Busath and D Henderson and S Sokolowski. ""Monte"
    write(unit=fidlog,fmt=*)"Carlo Simulations of the Mechanism of Channel Selectivity: The"
    write(unit=fidlog,fmt=*)"competition between Volume Exclusion and Charge Neutrality. "
    write(unit=fidlog,fmt=*)"J Phys Chem B, 2000, 104, 8903"
    write(unit=fidlog,fmt=*)
    write(unit=fidlog,fmt=*)"Reference for toroid geometry"
    write(unit=fidlog,fmt=*)"DezsH{o} Boda and Douglas Henderson and David D Busath, ""Monte"
    write(unit=fidlog,fmt=*)"Carlo study of the selectivity of calcium channels: improved"
    write(unit=fidlog,fmt=*)"geometrical model"", Molecular Physics, 2002, 100, 2361-2368"
    write(unit=fidlog,fmt=*)
    write(unit=fidlog,fmt=*)"Reference for canonical moves (radial,jump) and preference sampling"
    write(unit=fidlog,fmt=*)"(jump-in,jump-out):"
    write(unit=fidlog,fmt=*)"Dezso"" Boda, Douglas Henderson, and David Busath ""Monte Carlo" 
    write(unit=fidlog,fmt=*)"study of the selectivity of calcium channels: improved"
    write(unit=fidlog,fmt=*)"geometrical model"", Molecular Physics, 2002, 110(14) 2361-2368"
    write(unit=fidlog,fmt=*)
    write(unit=fidlog,fmt=*)"Reference for grand-canonical insertion and deletion"
    write(unit=fidlog,fmt=*)"Dezso"" Boda, Mo'nika' Valisko', Bob Eisenberg, Wolfgang Nonner,"
    write(unit=fidlog,fmt=*)"Douglas Henderson, and Dirk Gillespie ""The effect of protein"
    write(unit=fidlog,fmt=*)"dielectric coefficient on the ionic selectivity of a calcium channel"""
    write(unit=fidlog,fmt=*)"The Journal of Chemical Physics, 2006, 125, 034901"
    write(unit=fidlog,fmt=*)
    write(unit=fidlog,fmt=*)"References for induced charge:"
    write(unit=fidlog,fmt=*)"R. Allen and J.-P. Hansen and S. Melchionna ""Electrostatic potential"
    write(unit=fidlog,fmt=*)"inside ionic solutions confined by dielectrics: a variational approach"""
    write(unit=fidlog,fmt=*)"Phys Chem Chem Physics, 2001, 3, 4177-4186"
    write(unit=fidlog,fmt=*)"  this implementation based on:"
    write(unit=fidlog,fmt=*)"Dezso"" Boda, Dirk Gillespie, Wolfgang Nonner, Douglas Henderson, "
    write(unit=fidlog,fmt=*)"and Bob Eisenberg ""Computing induced charges in inhomogeneous "
    write(unit=fidlog,fmt=*)"dielectric media: Application in a Monte Carlo simulation of "
    write(unit=fidlog,fmt=*)"complex ionic systems"" Physical Review E, 2004, 69, 046702"
    write(unit=fidlog,fmt=*)
    if (calwid) then
      write(unit=fidlog,fmt=*)"References for estimation of chemical potential."
      write(unit=fidlog,fmt=*)"Widom, B, ""Some Topics in the Theory of Fluids"", J. Chem. Phys.,"
      write(unit=fidlog,fmt=*)"1963, 39(11), 2808-2812."
    endif
    write(unit=fidlog,fmt='(72("-"))')
    write(unit=fidlog,fmt=*)"Summary of the simulation step parameters (input section 'trial')"

    if (.not.allocated(hnw)) call d8v1p4

    ! normalise rates
    s_ = sum(rates_)
    if (dfeq(s_,0.D0)) stop "Error: sum of update type rates is zero."
    rates_=rates_/s_
    write(unit=fidlog,fmt='(72("-"))')
    write(unit=fidlog,fmt='(1X,4(A12,1X,"|"))')"Update type","move","salt add","indiv add"
    write(unit=fidlog,fmt='(1X,(A12,1X,"|"),4(F12.6,1X,"|"))')"probability",ratmv(),ratslt(),ratind()

    dosalt=.not.dfeq(ratslt(),0.D0).and.nsalt().ge.1
    doswap=.not.dfeq(ratswap(),0.D0).and.nsubs().ge.1
    doindv=.not.dfeq(ratind(),0.D0)

    write(unit=fidlog,fmt='(72("-"))')
    write(unit=fidlog,fmt=*)"Interpreted simulation step input"
    write(unit=fidlog,fmt='(72("-"))')
    call ectral(fidlog)
    write(unit=fidlog,fmt='(72("-"))')
  contains

    ! --------------------------------------------------------------
    ! Allocate internal data arrays
    subroutine d8v1p4
      use conf
      use geom
      use spec
      use patch
      use simstate
      implicit none
      integer :: i, j

      ! Coloumb interaction attribute
      allocate(riinw(confsz(),nnewmx))
      riinw=0
      do i=1,nspcmx
        call rs_init(achgeng(i))
        call rs_init(adeleng(i))
        call rs_init(amoveng(i))
        do j=1,nspcmx
          call rs_init(rdf_est(j,i))
        end do
      end do
      call rs_init(adeltaf)
      ! Widom sampling attributes
      if (calwid) then
        allocate(widobj)
        allocate(widobj%umfi(gz_max,nspec()))
        allocate(widobj%uchsi(gz_max,nspec()))
        allocate(widobj%uchs(gz_max,nspec()))
        allocate(widobj%trials(gz_max,nspec()))
        allocate(widobj%euchsi(gz_max,nspec()))
        allocate(widobj%euchs(gz_max,nspec()))
        allocate(widobj%succss(gz_max,nspec()))
        widobj%umfi(:,:) = 0.D0
        widobj%uchsi(:,:) = 0.D0
        widobj%uchs(:,:) = 0.D0
        widobj%euchsi(:,:) = 0.D0
        widobj%euchs(:,:) = 0.D0
        widobj%trials(:,:) = 0
        widobj%succss(:,:) = 0
        widobj%trys = 0
      endif

      ! Induced Charge Computation interaction attributes
      if (.not.is_homogeneous().and.do_electrostatic()) then
        allocate(hnw(npchmx),dcnw(npchmx),dh(npatch),ripnw(npchmx,nnewmx))
        dh = 0
        hnw = 0
        dcnw = 0
        ripnw = 0
      endif
    end subroutine d8v1p4

  end subroutine rftral

  ! --------------------------------------------------
  ! Indicate change type for create/destroy
  !
  ! @param itype [add1_l|addslt|rem1_l|remslt]
  ! @param ireg  -> region for add/del
  ! @param ispec -> specie _or_ salt to add/del
  subroutine selchg(itype,ireg,ispec)
    use spec
    implicit none
    integer, intent(in) :: itype, ireg,ispec
    if (itype.eq.0) stop "Error: A add/del change requires a type"
    if (ireg.lt.1.or.ireg.gt.4) stop "Error: Invalid region index"
    if (itype.lt.remslt.or.itype.gt.addslt) stop "Error: Invalid change type"
    istate(1)=itype
    istate(2)=0
    ispcnw=0
    iregnw=ireg
    indxnw=0
    i0i9p7: select case (itype)
    case(add1_l) i0i9p7
      nchgnw=1
      ispcnw(1)=ispec
      igcnw=0
      acrea1(1,ireg,ispec)=acrea1(1,ireg,ispec)+1
    case(rem1_l) i0i9p7
      nchgnw=1
      ispcnw(1)=ispec
      igcnw=0
      adest1(1,ireg,ispec)=adest1(1,ireg,ispec)+1
    case(addslt) i0i9p7
      acreat(1,ireg,ispec)=acreat(1,ireg,ispec)+1
      igcnw=ispec
      ispcnw(1)=isalt(igcnw)
      nchgnw=nint(xz(ispcnw(1)))+1
      ispcnw(2:nchgnw)=idxcl()
    case(remslt) i0i9p7
      adest(1,ireg,ispec)=adest(1,ireg,ispec)+1
      igcnw=ispec
      ispcnw(1)=isalt(igcnw)
      nchgnw=nint(xz(ispcnw(1)))+1
      ispcnw(2:nchgnw)=idxcl()
    end select i0i9p7
  end subroutine selchg

  ! --------------------------------------------------
  ! Indicate change type for movement
  !
  ! @param itype [sphr_l|jump_l|jmpin_|jmpin_]
  ! @param ispec -> specie moving
  !
  ! @pre itype in [sphr_l|jump_l|jmpin_|jmpin_]
  ! @pre ni(ispec) /= 0
  ! @post indxnw(1) /=0 && > nactv
  subroutine selmov(itype,ispec)
    use geom
    use spec
    use conf
    implicit none
    integer, intent(in) :: itype, ispec
    if (itype.lt.sphr_l.or.itype.gt.jmpou_) stop "Invalid move type"
    if (dbc) then
      if (ni(ispec).lt.0) stop "Error: [selmov] ni(ispec) == 0"
    endif
    istate(1)=0
    istate(2)=itype
    nchgnw=1
    ispcnw=0
    iregnw=4
    ispcnw(1)=ispec
    indxnw=0
    indxnw(1)=getnth(ispcbk,ispec,iranff(ni(ispec)),nactv)
    if (dbc) then
      if (indxnw(1).gt.nactv) stop "Error: [move] could not find nth particle (ni)"
      if (indxnw(1).lt.1) stop "Error: [move] nth particle (ni) was 0"
    endif
    igcnw=0
    y5c6c8: select case (itype)
    case(sphr_l) y5c6c8
      amove(1,ispec)=amove(1,ispec)+1
    case(jump_l) y5c6c8
      ajump(1,ispec)=ajump(1,ispec)+1
    case(jmpin_) y5c6c8
      u1o0a8: if (inchan(rz(indxnw(1)))) then
        ajout(1,ispec)=ajout(1,ispec)+1
        istate(2)=jmpou_
      else u1o0a8
        ajin(1,ispec)=ajin(1,ispec)+1
      endif u1o0a8
    case(jmpou_) y5c6c8
      u1o0a9: if (inchan(rz(indxnw(1)))) then
        ajout(1,ispec)=ajout(1,ispec)+1
      else u1o0a9
        ajin(1,ispec)=ajin(1,ispec)+1
        istate(2)=jmpin_
      endif u1o0a9
    end select y5c6c8
  end subroutine selmov

  ! --------------------------------------------------
  ! Indicate change type for subspecie swap
  !
  ! @param itype [swap_l]
  ! @param isubs -> subspecie to swap
  subroutine selswp(isubs)
    use spec
    use conf
    use geom
    implicit none
    integer, intent(in) :: isubs
    integer :: ispec1, ispec2
    ! select particle from either subspecie
    ispec1 = subspecie_index(isubs,1)
    ispec2 = subspecie_index(isubs,2)
    indxnw(1)=iranff(ni(ispec1) + ni(ispec2))
    ! if in second specie, swap ispec1 and ispec2
    if (indxnw(1).gt.ni(ispec1)) then
      indxnw(1) = indxnw(1) - ni(ispec1)
      call swap(ispec1, ispec2)
    endif
    indxnw = 0
    indxnw(1) = getnth(ispcbk,ispec1,indxnw(1),nactv)
    ispcnw = 0
    ispcnw(1) = ispec1
    ispcnw(2) = ispec2
    ! set istate(1) to zero to indicate no net particle number change
    istate(1) = move_l
    istate(2) = swap_l
    ! get region from the current particle.
    call inregn(rz(indxnw(1)),ispec1,iregnw)
    nchgnw = 1
    igcnw = 0
    adest1(1,iregnw,ispec1)=adest1(1,iregnw,ispec1)+1
    acrea1(1,iregnw,ispec2)=acrea1(1,iregnw,ispec2)+1
  end subroutine selswp

  ! --------------------------------------------------
  ! indicate that the change was accepted
  !
  ! Internal routine for managing statistics.
  subroutine on_cmt
  use conf
  implicit none
  ! Drift test variables
  integer :: ichg ! change counter
  integer, save :: cutoff = 0
  integer, save :: x0 = 0
  double precision, save :: sum_ = 0.D0
  double precision, save :: cnt_ = 0.D0
  double precision, save :: off_ = 0.D0
  if (.not.isamov()) then
    ! initialise x0
    if (0.eq.x0) call setx0
    cnt_ = cnt_ + 1.D0
    sum_ = sum_ + (ntot() - x0)**3

    if (cutoff.lt.abs(sum_/cnt_)) then
      if (cnt_.gt.off_) then
        write(unit=fidlog,fmt=*)'WARNING: Significant skewness in population number: '&
          ,sum_/(cnt_*sqrt(dble(x0)))
        write(unit=fidlog,fmt='(1X,A,1X,I4,1X,I4,1X,A,1X,I12)')&
          'WARNING: Original and current population:',x0,ntot(),'Steps:',int(cnt_)
        ! stop 'Skew in population number'
        off_=cnt_ + 1000.D0
        ! reset x0 to ntot()
        x0 = ntot()
      endif
    endif
  endif
  m6n6q3: select case (istate(1))
  case(addslt) m6n6q3
    acreat(2,iregnw,igcnw)=acreat(2,iregnw,igcnw)+1
  case(add1_l) m6n6q3
    acrea1(2,iregnw,ispcnw(1))=acrea1(2,iregnw,ispcnw(1))+1
  case(move_l) m6n6q3
    s2a7a2: select case (istate(2))
    case(sphr_l) s2a7a2
      amove(2,ispcnw(1))=amove(2,ispcnw(1))+1
    case(jump_l) s2a7a2
      ajump(2,ispcnw(1))=ajump(2,ispcnw(1))+1
    case(jmpin_) s2a7a2
      ajin(2,ispcnw(1))=ajin(2,ispcnw(1))+1
    case(jmpou_) s2a7a2
      ajout(2,ispcnw(1))=ajout(2,ispcnw(1))+1
    case(swap_l) s2a7a2
      adest1(2,iregnw,ispcnw(1))=adest1(2,iregnw,ispcnw(1))+1
      acrea1(2,iregnw,ispcnw(2))=acrea1(2,iregnw,ispcnw(2))+1
    end select s2a7a2
  case(rem1_l) m6n6q3
    adest1(2,iregnw,ispcnw(1))=adest1(2,iregnw,ispcnw(1))+1
  case(remslt) m6n6q3
    adest(2,iregnw,igcnw)=adest(2,iregnw,igcnw)+1
  end select m6n6q3
 
  s1y4m4: do ichg=1,nchgnw
    call on_cmt_rdf(ispcnw(ichg),ispcbk,riinw(1,ichg),nactv)
  enddo s1y4m4

  contains
  subroutine setx0
  use spec
  use geom
  implicit none
    integer :: igc

    do igc=1,nsalt()
      x0=x0+nint(ctargs(igc)*vtotal(isalt(igc))/tosi)
      x0=x0+nint(ctargs(igc)*vtotal(idxcl())*xz(isalt(igc))/tosi)
    enddo
    cutoff = x0*x0
  end subroutine setx0
  
  end subroutine on_cmt

  ! ----------------------------------------------------------------------
  ! Use statistics of RDF to estimate convergence to equilibrium
  !
  ! update rdf_est
  subroutine on_cmt_rdf(jspec, ispecs, dist, sz)
    use rnstat
    implicit none
    integer, intent(in) :: jspec, sz
    integer, dimension(sz), intent(in) :: ispecs
    double precision, dimension(sz), intent(in) :: dist
    integer :: idx
    do idx=1,sz
      if (ispecs(idx).ne.0) then
        call rs_push(rdf_est(ispecs(idx),jspec), dist(idx))
      endif
    end do
  end subroutine on_cmt_rdf

  ! ----------------------------------------------------------------------
  ! Print statistics of RDF used estimate convergence to equilibrium
  !
  ! print rdf_est
  subroutine print_rdf
    use rnstat
    use spec
    implicit none
    integer :: idx, jdx
    write(unit=fidlog,fmt="(17(A10))")"ION",(fspc(idx),idx=1,nspec())
    do idx=1,nspec()
      write(unit=fidlog,fmt="(A10,16(F10.4))")" MEAN:"//fspc(idx),(rs_mean(rdf_est(jdx,idx)),jdx=1,nspec())
      write(unit=fidlog,fmt="(A10,16(F10.4))")"  VAR:"//fspc(idx),(rs_variance(rdf_est(jdx,idx)),jdx=1,nspec())
    end do
    write(unit=fidlog,fmt="(70('-'))")
  end subroutine print_rdf

  ! ----------------------------------------------------------------------
  ! Attempt a move in simulation of bulk
  !
  ! ----------------------------------------------------------------------
  ! Attempt a simple move in equilibration phase of simulation of bulk
  !
  subroutine tryblk(doclamp)
    use conf
    use geom
    use spec
    implicit none
    logical, intent(in) :: doclamp

    double precision :: rndnum
    integer :: ispec
    logical :: ovrlap

    ! The particle count before an attempted add/remove
    integer :: n_orig

    rndnum = ranff()

    if (ranff().lt.0.5D0) then
      ispec = select_gc_specie(rndnum)
      !! Handle 'accept' method's special limits on add/rem
      if (doclamp) then
        ! 'accept' method alternate add/rem
        if (wasdel(ispec)) then
          call selchg(add1_l,4,ispec)
        else
          call selchg(rem1_l,4,ispec)
        endif
        n_orig = ni(ispec)
      else
        ! Other estimation methods
        if (ranff().lt.0.5D0) then
          call selchg(add1_l,4,ispec)
        else
          call selchg(rem1_l,4,ispec)
        endif
      endif
      ovrlap=.false.
      call chgblk(ovrlap)
      if (.not.ovrlap) call moveeq(dispbc,.true.)

      !! Handle 'accept' method's special limits on add/rem
      if (doclamp) then
        ! check if particle number changed
        if (n_orig.gt.ni(ispec)) then
          wasdel(ispec)=.true.
        elseif (n_orig.lt.ni(ispec)) then
          wasdel(ispec)=.false.
        endif
      endif
    else
      rndnum=ranff()
      ispec = select_mv_specie(rndnum, .true.)

      if (ni(ispec).eq.0) return
      !   choose kind of move
      if (ranff().gt.ratmov(ispec)) then
        call selmov(sphr_l,ispec)
      else
        call selmov(jump_l,ispec)
      endif
      ovrlap = .false.
      call chgblk(ovrlap)
      if (.not.ovrlap) call moveeq(dispbc,.true.)
    endif
  end subroutine tryblk

  ! ----------------------------------------------------------------
  ! Attempt an MC/GC move in simulation of channel
  !
  ! This is the routine for the main part of a GCMC simulation
  subroutine tryeq
    use conf
    use geom
    use spec
    implicit none

    ! LOCALS
    double precision :: rndnum  ! random number
    integer :: igc,ireg,ispec ! selected salt,region or specie and general index
    logical :: ovrlap

    rndnum = ranff()
    ovrlap = .false.

    k7q3s9: if (dosalt.and.rndnum.lt.ratslt()) then

      !   try to make salt GC
      !   choose salt
      rndnum=ranff()
      ! start with igc = nsalt() so to avoid invalid value at end.
      igc = select_salt(rndnum)

      ispec = isalt(igc)
      !   choose region
      rndnum = ranff()
      ! set ireg to ibulk so region will always be valid
      ireg = select_region(ispec, rndnum)

      !   choose GC step
      o4o8s1: if (ranff().lt.0.5D0) then
        call selchg(addslt,ireg,igc)
        call chgsim(ovrlap)
        if (ovrlap) return
        call moveeq(disbox,.true.)
      else o4o8s1
        call selchg(remslt,ireg,igc)
        call chgsim(ovrlap)
        if (ovrlap) return
        call moveeq(disbox,.true.)
      endif o4o8s1

    elseif (doindv.and.rndnum.lt.(ratslt()+ratind())) then k7q3s9

      !   try to make individual GC
      rndnum=ranff()
      !   choose ion species
      ispec = select_gc_specie(rndnum)
      !   choose region
      rndnum=ranff()
      ! set ireg to ibulk so region will always be valid
      ireg = select_region(ispec, rndnum)

      !   choose GC step
      w7p4i1: if (ranff().lt.0.5D0) then
        call selchg(add1_l,ireg,ispec)
        call chgsim(ovrlap)
        if (ovrlap) return
        call moveeq(disbox,.true.)
      else w7p4i1
        call selchg(rem1_l,ireg,ispec) 
        call chgsim(ovrlap)
        if (ovrlap) return
        call moveeq(disbox,.true.)
      endif w7p4i1

    elseif (doswap.and.rndnum.lt.(ratslt()+ratind()+ratswap())) then k7q3s9

      !   try to make swap move
      rndnum=ranff()
      !   choose subspecies
      ispec = select_subs(rndnum)
      
      call selswp(ispec)
      call chgsim(ovrlap)
      if (ovrlap) return
      call moveeq(disbox,.true.)

    else k7q3s9

      !   try to make an ion move
      !   choose ion species
      rndnum=ranff()
      ispec = select_mv_specie(rndnum, .false.)

      if (ni(ispec).eq.0) return
      !   choose kind of move

      call selmov(select_move(ispec, ranff()),ispec)
      call chgsim(ovrlap)
      if (ovrlap) return
      call moveeq(disbox,.true.)
    endif k7q3s9

  end subroutine tryeq

  ! ----------------------------------------------------------------------
  ! Attempt a simple move in pre-equilibration phase of simulation of channel
  !
  ! In this version of the trial moves we make all moves types equally
  ! likely for all species.
  subroutine trypeq
    use conf
    use geom
    use patch
    use spec
    implicit none

    logical :: ovrlap
    ! variable to force add/del operations to alternate
    integer :: ispec, i
    ! select ispec from randomly chosen particle

    do while (.true.)
      ! Randomly select any specie
      ispec = ispcbk(max(1,ceiling (ranff() * nactv)))

      ! capture limiting case
      if (0.eq.ispec) cycle

      i = ni(ispec)

      !! BREAKS WHEN REMOVING DO LOOP
      if (ranff().lt.0.5D0) then
        ! add/remove (only free species)
        if (.not.isfree(ispec)) cycle

        if (wasdel(ispec)) then
          call selchg(add1_l,ibulk,ispec)
        else
          call selchg(rem1_l,ibulk,ispec)
        endif
      else
        if (ranff().lt.0.75D0) then
          if (localized(ispec)) then
            call selmov(sphr_l,ispec)
          elseif (ranff().lt.0.33D0) then
            call selmov(jump_l,ispec)
          else
            call selmov(jmpin_,ispec)
          endif
        else
          call selmov(sphr_l,ispec)
        endif
      endif
      call chgsim(ovrlap)
      if (ovrlap) return
      call moveeq(disbox,.true.)
      ! check if particle number changed
      if (i.ne.ni(ispec)) wasdel(ispec)=.not.wasdel(ispec)
      return
    end do

  end subroutine trypeq

  ! ----------------------------------------------------------------
  ! Attempt GCMC single insertions in a simulation of channel
  ! as part of widom method for calculating chem. pot.
  subroutine trywd(asteps)
  use conf
  use geom
  use spec
  use simstate
  implicit none
  integer, intent(in) :: asteps ! iteration steps

  ! LOCALS
  integer :: ispec ! specie index
  integer :: ireghi ! region index
  logical :: ovrlap
  integer :: ntotal ! total number of particles
  integer :: bin ! gz bin index
  integer :: ntrgt ! target number of attempted individual insertions
  double precision :: lo, hi, mid ! start, end and midpoint of a gz bin
  double precision :: rbox ! allowed max radius
  if (dbc) then
    if (asteps.eq.0) stop "Attempt to perform Widom sampling on step zero"
  endif
  if (0.eq.nwdtry) return
  call start_widom_trials

  ! calculate current trial target
  ntrgt=nwdtry*asteps/1000
  if (0.eq.ntrgt) return

  ! reset the internals as per selchg
  istate(1)=add1_l
  istate(2)=0
  ispcnw=0
  iregnw=0
  ireghi=0
  indxnw=0
  nchgnw=1
  igcnw=0
  ntotal=ntot()

  ! THE propscal is not used in the recording of data for the Widom
  ! inserton method, so is arbitrarily set to one.
  propscal=1.D0

  ! The idea here is to normalise the sampling for each bin, we work through the 
  ! bins for each specie updating them only if they have insufficient trials at
  ! this point in the calculation.
  do ispec=idxcl(),nspec()
    ! Use selchg to set up state
    call selchg(add1_l,ibulk,ispec)
    ! set chemical potential
    propexpo=chempi(ispec)
    ! back track acrea1 update from selchg
    acrea1(1,ibulk,ispec)=acrea1(1,ibulk,ispec)-1
    do bin=1,gz_max
      if (debug) then
         write(unit=fidlog,fmt=*)"! bin[",bin,"] specie[",fspc(ispec),"] trials[", &
           & widobj%trials(bin,ispec),"]:[",ntrgt,"]"
      endif
      if (widobj%trials(bin,ispec).ge.ntrgt) cycle

      lo  = gz_lo(bin)
      mid = gz_mid(bin)
      hi  = gz_hi(bin)
      ! Find regions at top and bottom of interval
      call inregn(lo,ispec,iregnw)
      call inregn(hi,ispec,ireghi)
      ! Reset iregnw to choose outer of two possible region
      ! ASSUME: higher region encompasses all regions with lower numbers
      iregnw = max(ireghi, iregnw)
      iregnw = max(izlim, iregnw)
      if (dbc) then
        if (iregnw.lt.1.or.iregnw.gt.4) stop "Invalid region"
      endif

      rbox = rreg(iregnw,ispec)
      if (iregnw.eq.ibulk) then
        propscal=vreg(iregnw,ispec)/(ni(ispec)+1)
      else
        propscal=vreg(iregnw,ispec)/(nin(iregnw,ispec)+1)
      endif

      do while (widobj%trials(bin,ispec).lt.ntrgt)
        ! set up particle
        call jmpmov(rznw(1),r2nw(1),lo,hi,rbox)
        ! register attempt
        widobj%trials(bin,ispec)=widobj%trials(bin,ispec)+1
        widobj%trys=widobj%trys+1
        acrea1(1,iregnw,ispec)=acrea1(1,iregnw,ispec)+1

        ! check particle in allowed geometry
        ovrlap=.true.
        call wall(ispec,rznw(1),r2nw(1),ovrlap)
        if (ovrlap) cycle ! repeat if overlap in wall

        call jmpfin(rxnw(1),rynw(1),r2nw(1))

        ! run trial
        call moveeq(disbox,.false.)
      enddo
    enddo
  enddo
  if (ntotal.ne.ntot()) then
    write(unit=fidlog,fmt=*)"PROGRAM ERROR: Total number of particles changed during Widom sampling for Chem. Pot.."
    stop "Change in particle number in Widom routine"
  endif
  call end_widom_trials
  end subroutine trywd

  ! -------------------------------------------------- 
  !  Record the number of successful trial adds in a gz bin
  !
  !  pre: calwid
  subroutine widom_trial(rznew,ispec,umfi,uchs,uchsi)
  use spec
  use geom
  implicit none
  double precision, intent(in) :: rznew
  integer, intent(in) :: ispec
  double precision, intent(in) :: umfi, uchs, uchsi
  integer :: gzbin_
  if (dbc) then
    if (ispec.lt.1.or.ispec.gt.nspec()) stop "specie index out of range"
  endif

  gzbin_=gz_bin(rznew)

  widobj%umfi(gzbin_,ispec)=widobj%umfi(gzbin_,ispec)+umfi
  widobj%uchs(gzbin_,ispec)=widobj%uchs(gzbin_,ispec)+uchs
  widobj%uchsi(gzbin_,ispec)=widobj%uchsi(gzbin_,ispec)+uchsi
  widobj%euchs(gzbin_,ispec)=widobj%euchs(gzbin_,ispec)+exp(-1.0*uchs)
  widobj%euchsi(gzbin_,ispec)=widobj%euchsi(gzbin_,ispec)+exp(-1.0*uchsi)
  widobj%succss(gzbin_,ispec)=widobj%succss(gzbin_,ispec)+1

  end subroutine widom_trial

  ! -------------------------------------------------- 
  ! write results of the widom calculation
  subroutine savewd
    use spec
    use conf
    use geom
    use strngs
    implicit none
    ! (NOTE: in units of RT)
    !
    ! Free energy = ln ( mean(Euchs) ) = ln (mean([])/[]o) + mu_ex
    !  where for some bin:
    !  * the standard concentration []o == 1 
    !  * Euchs = euchs + (trials - success) : to add 1{=e^0} for each failed trial
    !  * mean(Euchs) = Euchs/trials
    !  * mean([]) = mean(gz) / gzvol = gzaver / gzvol
    !
    ! Therefore for bin number ibin and specie ispec
    ! gzaver(ibin,ispec,ninbin)
    ! mean_Euchs = (euchs(ibin,ispec)+(trials(ibin,ispec)-succss(ibin,ispec)))/trials(ibin,ispec)
    ! mean_conc  = ninbin/(tosi*gzvol(ibin))
    ! mu_ex = log(mean_Euchs) - log(mean_conc)
    !
    !
    interface
      subroutine gzaver(ibin,ispec,naver)
        integer, intent(in) :: ibin, ispec
        double precision, intent(out) :: naver
      end subroutine gzaver
    end interface

    integer :: ispec,idx
    double precision :: ninbin    ! average particle count in a bin

    double precision :: mu        ! chem pot of current bin
    double precision :: mu_ex     ! excess chem pot of current bin
    double precision :: mid       ! absolute value of mid-point
    double precision :: mu_sum, mu_blk, minblk  ! sum (total and bulk)
    double precision :: mu_cnt, mu_cbk  ! sample count (total and bulk)
    double precision :: ratio1, ratio2
    write(unit=fidlog,fmt='(72("-"))')
    write(unit=fidlog,fmt=*)"Chemical potential estimated using insertion sampling (Widom). First"
    write(unit=fidlog,fmt=*)"value is estimated from subpart of bulk region, second value is over"
    write(unit=fidlog,fmt=*)"entire simulation (last value is the minimum value in bulk region)."
    write(unit=fidlog,fmt='(1X,A3,"   ",A16,"   ",A16,"   ",A16)')"ION","CHEM. POT. BULK",&
         &"CHEM. POT. AVER.","BULK MIN"
    write(unit=fidlog,fmt='(72("-"))')
    do ispec=idxcl(),nspec()
      mu_sum=0
      mu_blk=0
      mu_cnt=0
      mu_cbk=0
      minblk=0
      open(unit=fidtry,file='res/wid-'//fspc(ispec)//'.'//firun//'.dat')
      write(unit=fidtry,fmt=*)'# UUID ',fuuid
      write(unit=fidtry,fmt='(A)')'# label zcoord    mu   mu_ex  trials  umfi  uchs  uchsi  euchs  euchsi vol   succss'
      write(unit=fidtry,fmt='(A)')'# unit  ang       ENG  ENG    count   ENG   ENG   ENG    PROB   PROB   ang3  count'

      do idx=1,gz_max
        ! term used to adjust Chem. Pot. for ideal gas [V/N+1]
        ! skip bins with no trial
        if (widobj%trials(idx,ispec).ne.0) then
          ! average concentration
          call gzaver(idx,ispec,ninbin)
          ! calculate mean_Euchs first

          mid=abs(gz_mid(idx))
          ! We add widobj%trials(idx,ispec)-widobj%succss(idx,ispec) as these are the number
          ! of trials that failed metropolis. They contribute 0 to the uchs but 1 to euchs
          ! 
          mu = widobj%euchs(idx,ispec) + dble(widobj%trials(idx,ispec) - widobj%succss(idx,ispec))
          mu = mu / dble(widobj%trials(idx,ispec))
          ! now chem pot
          if (mu.le.0.D0) then
            write(fidlog,*)"WARNING: skipping mu value less than 0 [",mu,"] for specie ",fspc(ispec)," at GZ index ",idx
            mu = 0.D0
            mu_ex = 0.D0
          else
            mu=-log(mu)
            ! now chem excess
            if (dfeq(ninbin,0.D0)) then
              mu_ex=0.D0
            else
              mu_ex=mu-log(tosi*ninbin/gz_vol(idx,ispec))
            endif
            mu_sum=mu_sum+mu
            mu_cnt=mu_cnt+1
            if (mid.ge.zbulk1.and.mid.le.zbulk2) then
              minblk=min(minblk,mu)
              mu_blk=mu_blk+mu
              mu_cbk=mu_cbk+1
            endif
          endif
          write(unit=fidtry,fmt='(1X,3(G16.4),I9,6(G16.4),I9)')gz_mid(idx),mu,mu_ex            &
               & ,widobj%trials(idx,ispec),widobj%umfi(idx,ispec)&
               & ,widobj%uchs(idx,ispec),widobj%uchsi(idx,ispec) &
               & ,widobj%euchs(idx,ispec),widobj%euchsi(idx,ispec)&
               & ,gz_vol(idx,ispec),widobj%succss(idx,ispec)
        endif
      enddo
      close(fidtry)
      ratio1=0
      if (.not.dfeq(mu_cbk,0.D0)) ratio1 = mu_blk / mu_cbk
      ratio2=0
      if (.not.dfeq(mu_cnt,0.D0)) ratio2 = mu_sum / mu_cnt
      write(unit=fidlog,fmt='(1X,A3,"   ",F16.4,"   ",F16.4,"   ",F16.4)')fspc(ispec),ratio1,ratio2,minblk
    enddo
    write(unit=fidlog,fmt='(72("-"))')
  end subroutine savewd
 
  ! -------------------------------------------------- 
  ! Reset the energy factor statistic accumulators
  subroutine zeroav
  implicit none
  integer :: i
  do i=1,nspcmx
    call rs_reset(achgeng(i))
    call rs_reset(adeleng(i))
    call rs_reset(amoveng(i))
  end do
  call rs_reset(adeltaf)
  end subroutine zeroav

  ! ----------------------------------------------------
  ! reset the statistics counters
  subroutine zerotr
  implicit none

  amove=0
  ajump=0
  ajin=0
  ajout=0
  acreat=0
  adest=0
  acrea1=0
  adest1=0
  if (calwid) then
     call zerowd
  endif   
  end subroutine zerotr

  ! ----------------------------------------------------
  ! reset the widom energy accumulators
  subroutine zerowd
    implicit none
    if (associated(widobj)) then
      widobj%umfi = 0.D0
      widobj%uchsi = 0.D0
      widobj%uchs = 0.D0
      widobj%euchsi = 0.D0
      widobj%euchs = 0.D0
      widobj%trials = 0
      widobj%succss = 0
      widobj%trys = 0
    endif
  end subroutine zerowd

end module trial
