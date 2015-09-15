
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
! PARTICLE POSITION DATA
! ----------------------
!
! [gindx=global index; Gindx=limited global index; lindx=specie
! local index; ispec=specie index; ireg=region index]
!
! Analysis of access patterns in the channel16 program suggested
! merging or removing a number of arrays and changing of access.
! The major change has been to make all the particles consecutive
! in the arrays.  This allows rapid traversal of particle data
! using a 'do gindx=1,nactv' loop.
!
! When a particle is deleted its index is added to 'idelst',
! where it is reused for the next addition, and the entry in
! ispcbk is set to zero.  Thus there are unused entries in the
! main arrays for deleted particles that are tested for using
! 'ispcbk(gindx).ne.0'.
!
! Having a compact data array and deletion list requires only
! a little programmer management.
!
! (1) For particle index checking the following relation holds:
!   ntot = nactv - ndel
! where ntot is the number of particles, nactv the maximum
! particle index and ndel the number of 'inactive' entries (also
! the length of the idelst).
!
! This means that you access the particle arrays between 1 and
! 'nactv' not 1 and 'ntot'.
!
! (2) The standard particle traversal is now:
!
! do ii=1,nactv
!   ispec=ispcbk(ii)
!   if (ispec.eq.0) cycle ! VERY IMPORTANT !
!   ...
!
! compared to 
!
! do ispec=1,nspec()
!   do i=1,n(ispec)
!     ii=indspc(i)
!     ...
!
! Not only is the new loop linear access to ispcbk (and rx etc)
! but it is also a single loop. Using this loop now meant 'indspc'
! and 'indreg' were infrequently called.  The code for managing
! these arrays was removed and they are now provided as functions
! that finds the nth occurence of the specie in ispcbk - hence
! it is now a lot more costly!
!
! Conversely, the new loop structure requires many more access to
! the specie data.  However, this data set is comparatively tiny
! and should always be fast.
!
! Other key changes are:
!
!   rij(ii,jj) was also mainly used as 1/2*rij(ii,jj).  Additionally
!   it was used with q(i) and q(j).  Thus the data was transformed
!   to:
!     q(i)*q(j)/(2*rij(ii,jj)) --> rqqii(ii,jj).
!
! While the main position data is stored in the rxyo7r array,
! this is private.  So access is now through functions rx, ry,
! rz, r2 which appear like the previous array accesses.
!
! The change also reduced the complexity of adding or removing
! a particle.  The cost of managing the map of particles to
! region and specie was reduced to just a map of particles to
! region.
!
! ------------------------------------------------------------
! ACCESS SUMMARY
!
! ARRY: rx,ry,rz,r2(gindx) -> matrix of x, y, z coords and r (=sqrt(x*x+y*y))
! access via:
! METH: setpos(gindx,x,y,z,r) !!NOTE this does not manage entries
!      in other arrays such as indreg, ispcbk!!
!
! ARRY: rqqii(gindx) -> matrix of qi * qj / inter-particle distance
!                 or -> inter-particle distance {qi = 0 or qj = 0}
!       NOTE: when one of the particles has zero charge, rqqii is the 
!             simple inter-particle distance 
!
! ARRY: ispcbk(gindx) -> particle specie
! ARRY: ni(ispec) -> number of a specie
!
! METH: ntot() -> total number of particles (ntot = nactv - ndel)
!
! ARRY: idellst(??) = deletion list (sorted)
! INT : nactv = number of active array entries
! INT : ndel = number of deleted particle entries
!
! METH: indreg(lindx,ireg,ispec) function to map specie/region index to global index
!
! METH: indspc(lindx,ispec) = function to map specie index to global index
!
! TODO: The mechanics of managing the data for a particle position
! are too exposed.  Currently, other modules must manually manage
! the changes to region occupancy via addreg/setreg/delreg.  This
! should all be hidden within this module.
!
module conf
use const
implicit none
private
  ! The allocated size of the conf arrays (<= ntotmx)
  ! access via confsz()
  integer, private :: ntotsz = ntotmx

  ! rqqii -> q_i*q_j/(2 * inter-particle distance)
  double precision, public, allocatable, dimension(:,:) :: rqqii

  ! rxyo7r(gindx) -> particle position and radial from z axis {sqrt(x^2+y^2)}
  ! ACCESS:
  ! rx,ry,rz -> x,y,z coordinates
  ! r2 -> radial from z-axis
  double precision, public, allocatable, dimension(:) :: rx,ry,rz,r2

  ! charge -> the charge of the system
  double precision, public :: charge

  ! indreg -> map of particle to region used by indreg function
  integer, private, allocatable, dimension(:,:) :: p4c6h2

  ! ispcbk -> specie of particle
  integer, public, allocatable, dimension(:) :: ispcbk

  ! idelst -> list of currently inactive indices (of rx,ry etc)
  integer, private, allocatable, dimension(:) :: idelst

  ! ni -> number of particles in a specie
  integer, public, dimension(nspcmx) :: ni

  ! nreg -> number of particles in a particular region
  integer, public, dimension(nrgnmx) :: nreg

  ! nin -> number of particle specie is in a particular region
  integer, public, dimension(nrgnmx,nspcmx) :: nin

  ! array sizes
  ! nactv -> maximum in-use index in rx,ry etc
  ! ndel -> number of deleted particles in idelst
  integer, public :: nactv,ndel

  ! Set to true to append iteration number to end of 
  ! filename in writcf.  This will make each save
  ! unique. (Currently this is compile only option)
  logical, public :: byiter=.false.

  ! Indicate if the grid or random placement initial
  ! configuration sub-algorithm is to be used.
  logical, public :: usegrid=.false.

  ! Set to true to indicate that a more relaxed 
  ! check of mobile ion positions is made during
  ! the thermalisation phase
  logical, public :: relaxmob=.false.

  public :: lookup, wrconf, rfconf, genrcf, genrbk, readcf, writcf
  public :: addreg, delreg, setreg, idxget, idxrel, confsz
  public :: indspc, indreg, chkmob, ntot, volslb
  public :: setpos
  ! verification functions
  public :: check_specie_count
contains
  ! --------------------------------------------------
  ! check that the number of particles of a particular
  ! specie found from 'ispcbk' match 'ni' 
  ! 
  subroutine check_specie_count(msg)
    use spec
    implicit none
    character(len=*), intent(in) :: msg
    integer, dimension(0:nspcmx) :: local_ni
    integer :: idx
    logical :: fail
    fail = .false.
    local_ni = 0
    do idx = 1,nactv
      local_ni(ispcbk(idx)) = local_ni(ispcbk(idx)) + 1
    enddo
    do idx = 1,nspec()
      if (local_ni(idx).ne.ni(idx)) fail = .true.
    enddo
    do idx = nspec()+1,nspcmx
      if (local_ni(idx).ne.0) fail = .true.
    enddo
    if (fail) then 
      write(*,*)msg
      write(*,fmt="(1X,A13,' = ',I4)")"empty places",local_ni(0)
      do idx = 1,nspec()
        write(*,fmt="(1X,A2,1X,A10,' = ',I4,' expected ',I4)")fspc(idx),'found',local_ni(idx),ni(idx)
      enddo
      do idx = nspec()+1,nspcmx
        write(*,fmt="(1X,A13,' = ',I4,I4,' expected 0')")'non-specie',local_ni(idx)
      enddo
      stop "Mismatch between content of 'ispcbk' and 'ni'"
    endif
  end subroutine check_specie_count


  ! Access functions
  ! ------------------------------------------------------------------
  ! Get the maximum possible number of particles
  pure integer function confsz()
  implicit none
  confsz=ntotsz
  end function confsz

  ! ------------------------------------------------------------------
  ! Set the x,y,z and r coordinates of particle i
  !
  ! NOTE: This method _only_ sets the coordinates, it does not manage
  ! any of the other arrays required to manage a particle
  subroutine setpos(ii,x,y,z,r)
  implicit none
  integer, intent(in) :: ii
  double precision, intent(in) :: x,y,z,r
  if (dbc) then
    if (dbc_level.ge.dbc_index) then
      if (ii.gt.nactv) stop "Index out of range of active particle set"
      if (ii.lt.1) stop "Index out of range of active particle set"
    end if
  endif
  rx(ii)=x
  ry(ii)=y
  rz(ii)=z
  r2(ii)=r
  end subroutine setpos

  ! ------------------------------------------------------------------
  ! Get the global index of a particle given a per-specie index
  ! and specie
  !
  ! This method finds the nth occurence of ispec in ispcbk. This
  ! means that it is reasonable fast, but not as fast as a direct
  ! array lookup. DO NOT USE THIS METHOD TO LOOP THROUGH ALL
  ! PARTICLES OF A SPECIE, if that is required use something like:
  !
  ! do ii=1,nactv
  !   if (ispcbk(ii).eq.ispec) then
  !     ...
  !
  ! NOTE: in earlier code scanning through per specie was generally
  ! only done when gathering data for all species.  Obviously, such
  ! code can be transformed into a single loop over particles with
  ! collecting data for all species.
  integer function indspc(lindx,ispec)
  implicit none
  integer, intent(in) :: lindx, ispec
  if (dbc) then
    if (dbc_level.ge.dbc_index) then
      if (lindx.gt.ni(ispec)) stop "Index out of range of active particle set"
      if (lindx.lt.1) stop "Index out of range of active particle set"
    end if
  endif
  indspc=getnth(ispcbk,ispec,lindx,nactv)
  end function indspc

  ! ------------------------------------------------------------------
  ! Get the global index of a particle given a per-specie index,
  ! region and specie
  !
  ! This module uses array 'p4c6h2' to map particles of any specie
  ! into a region.  Thus to find the nth particle of a specie in
  ! a region it looks for the nth entry in 'p4c6h2' with the desired
  ! specie type in 'ispcbk'. DO NOT USE THIS METHOD TO LOOP THROUGH
  ! ALL PARTICLES OF A SPECIE IN A REGION, if that is required use
  ! something like:
  !
  ! do ii=1,nactv
  !   if (ispcbk(ii).eq.ispec) then
  !     if (rz(ii).gt...) then
  !       ...
  !
  ! NOTE: in earlier code scanning through per specie nd region
  ! was generally only done when gathering data for all species.
  ! Obviously, such code can be transformed into a single loop
  ! over particles with collecting data for all species and region.
  integer function indreg(lindx,ireg,ispec)
  implicit none
  integer, intent(in) :: lindx, ireg, ispec
  integer :: i
  ! getnof returns index in p4c6h2
  if (dbc) then
    if (dbc_level.ge.dbc_index) then
      if (ireg.lt.1.or.ireg.gt.ibulk) stop "Region index out of range"
      if (lindx.gt.nin(ireg,ispec)) stop "Index out of range of active particle set"
      if (lindx.lt.1) stop "Index out of range of active particle set"
    end if
  endif
  indreg=getnof(p4c6h2(:,ireg),ispcbk,ispec,lindx,nreg(ireg))
  if (debug) then
    if (indreg.gt.nreg(ireg)) then
      write(unit=fidlog,fmt=*)"! R[",ireg,"]S[",ispec,"]NIN[",nin(ireg,ispec),"]INDX[" &
                  & ,lindx,"]NREG[",nreg(ireg),"]NACTV[",nactv,"]: "
      do i=1,nreg(ireg)
        write(unit=fidlog,fmt='(1X,"!  indx, spec",I4,1X,I4)')p4c6h2(i,ireg),ispcbk(p4c6h2(i,ireg))
      enddo
      stop 1
    endif
  endif
  ! convert p4c6h2 index into global index
  indreg=p4c6h2(indreg,ireg)
  if (dbc) then
    if (dbc_level.ge.dbc_ensure) then
      if (indreg.lt.1.or.indreg.gt.nactv) stop "Index out of range of active particle set"
      if (ispcbk(indreg).ne.ispec) stop "Specie of particle does not math request"
    end if
  endif
  end function indreg

  ! ------------------------------------------------------------------
  ! Get the total number of particles
  !
  pure integer function ntot()
  implicit none
  ntot=nactv-ndel
  end function ntot

  ! ------------------------------------------------------------------
  ! MANAGE indreg: add particle to a region
  !
  ! Add a particle to a region index if rz is in region's z scope
  !
  ! Avoid checking radial as we assume particle is in a valid position
  !
  ! @param ispec   particle specie
  ! @param i       particle's _global_ index
  ! @param rz      particle's xz coordinate

  subroutine addreg(ispec,i,rz)
  use spec
  use geom
  implicit none
  integer, intent(in) :: ispec,i
  double precision, intent(in) :: rz

  ! LOCALS
  integer ireg
  if (dbc) then
    if (dbc_level.ge.dbc_index) then
      if (ispec.lt.1.or.ispec.gt.nspec()) stop "Error: specie index is out of range"
      if (i.lt.1.or.i.gt.nactv) stop "Error: particle index out of range"
    endif
    if (dbc_level.ge.dbc_require) then
      if (ispcbk(i).ne.ispec) stop "Error: specie does not match index"
    endif
  endif
  k0t2t5: do ireg=1,3
    if ((abs(rz).le.zreg(ireg,ispec))) &
      & call ce1c59(ireg,ispec,i)
  enddo k0t2t5

  end subroutine addreg

  ! ------------------------------------------------------------------
  ! Private subroutine to actually add elements to a region
  !
  ! @param ireg : region
  ! @param ispec : specie
  ! @param ii : global index
  subroutine ce1c59(ireg,ispec,ii)
  implicit none
  integer, intent(in) :: ireg,ispec,ii

  nreg(ireg)=nreg(ireg)+1
  nin(ireg,ispec) = nin(ireg,ispec)+1
  p4c6h2(nreg(ireg),ireg) = ii
  end subroutine ce1c59

  ! ------------------------------------------------------------------
  ! MANAGE indreg: remove a particle from a region
  !
  ! Remove a particle from a region index if rz it was in region's z scope
  !
  ! Avoid checking radial as we assume particle is in a valid position
  !
  ! @param ispec   particle specie
  ! @param ii      particle's _global_ index
  ! @param rz      particle's xz coordinate

  subroutine delreg(ispec,ii,rz)
  use spec
  use geom
  implicit none
  integer, intent(in) :: ispec,ii
  double precision, intent(in) :: rz

  ! LOCALS
  integer ireg

  if (dbc) then
    if (dbc_level.ge.dbc_index) then
      if (ispec.lt.1.or.ispec.gt.nspec()) stop "Error: specie index is out of range"
      if (ii.lt.1.or.ii.gt.nactv) stop "Error: particle index out of range"
    endif
    if (dbc_level.ge.dbc_require) then
      if (ispcbk(ii).ne.ispec) then
        write(unit=fidlog,fmt=*)"Error in delreg: specie arg ",ispec," does not match ispcbk: ",ispcbk(ii)
        stop "Error: specie does not match index"
      endif
    endif
  endif
  j2l6c1: do ireg=1,3
    if ((abs(rz).le.zreg(ireg,ispec))) &
      & call ed73b5(ireg,ispec,ii)
  enddo j2l6c1

  end subroutine delreg
 
  ! ------------------------------------------------------------------
  ! Private subroutine to actually delete elements from a region
  !
  ! @param ireg : region
  ! @param ispec : specie
  ! @param ii : global index

  subroutine ed73b5(ireg,ispec,ii)
  implicit none
  integer, intent(in) :: ireg,ispec,ii

  ! LOCALS
  integer j1,j,i1

  j1=ii
  y6d3i7: do j=nreg(ireg),1,-1
    i1=p4c6h2(j,ireg)
    p4c6h2(j,ireg)=j1
    o1v9k6: if (ii.eq.i1) then
      nin(ireg,ispec)=nin(ireg,ispec)-1
      nreg(ireg)=nreg(ireg)-1
      return
    endif o1v9k6
    j1=i1
  enddo y6d3i7
  stop "Error deleting particle from region: possibly corrupt indices."
  end subroutine ed73b5
 
  ! ------------------------------------------------------------------
  ! GENERATE AN INITIAL MC CONFIGURATION FOR SIMULATION OF BULK
  ! 
  ! This method sets the program up to do a simulation of bulk 
  ! solution with periodic boundary conditions.  At the time it
  ! is called the system will have already started to generate
  ! the channel simulation.  The channel simulation data needs
  ! to be saved before calling this method as the conf module
  ! is reset.
  !
  ! The bulk simulation is performed in a cube with periodic
  ! boundary conditions.  The structural ions are ignored and a
  ! system containing just the free salts is generated.  Position
  ! for ions is random in all three dimensions and a test is made
  ! to ensure no overlap (including overlap occuring across a
  ! boundary.
  !
  ! This method generates entries for the added ions in the rqqii
  ! array as it progresses so the 'lookup' method does not need
  ! to be called after this method ends for the new ions.
  subroutine genrbk
  use grid
  use geom
  use spec
  implicit none

  ! externals
  double precision, external :: epsblk
  ! Locals
  double precision :: rzinw_,r2inw_,rxinw_,ryinw_ ! trial data
  integer :: ispec_, igc ! Specie and salt indices 
  integer :: ii_      ! particle global indices
  integer :: lidx_,itrys_ ! counters
  logical :: ovrlap
  double precision :: ri,rii
  integer :: jspec,jj
  integer :: nchg ! charge as an integer
  integer :: ntrg_ ! target number of particles
  double precision :: rmolar ! ratio of particles per molar fraction
  double precision :: length ! grid spacing
  integer :: particlecount ! Number of particles to create
 
  type (cubegrid) :: genrgrid ! possible grid for create particels

  if (dbc) then
    if (dbc_level.ge.dbc_check) then
      if (.not.allocated(rqqii)) stop "Error: conf%genrbk called before initialising conf module"
    endif
  endif

  write(unit=fidlog,fmt='(72("-"))')
  write(unit=fidlog,fmt=*)" Making initial simulation of bulk configuration"
  write(unit=fidlog,fmt='(72("-"))')

  ! -- reset conf
  if (nactv > 0) then
    p4c6h2=0
    nreg=0
    nin=0
    ispcbk=0
    ni=0
    nactv=0
    ndel=0
    rqqii=0  
  endif

  ! ------- Calculate particle numbers --------------------
  rmolar=volblk()/tosi
  do igc=1,nsalt()
    ni(isalt(igc))=ni(isalt(igc))+nint(ctargs(igc)*rmolar)
    ni(idxcl())=ni(idxcl())+nint(ctargs(igc)*rmolar*xz(isalt(igc)))
  enddo

  ! Double check for electroneutrality due to possible rounding
  charge=0
  m0u9b7: do ispec_=idxcl(),nspec()
    charge=charge+xz(ispec_)*ni(ispec_)
  enddo m0u9b7
  ! -----recompute N(Cl) if cell is not neutral ------------
  if (.not.dfeq(charge,0.D0)) then
    nchg=nint(charge)
    write(unit=fidlog,fmt='(1X,A,I4)')"Total cell charge = ",nchg
    write(unit=fidlog,fmt='(1X,"Cell is not neutral, recomputing (Cl) number")')
    write(unit=fidlog,fmt=*)"Old N(Cl) = ",ni(idxcl())
    if (nint(xz(idxcl())).eq.-1) then
      ni(idxcl())=nchg+ni(idxcl())
    else
      stop "Input error: valency of specie (Cl) should be -1!"
    endif
    write(unit=fidlog,fmt=*)"New N(Cl) = ",ni(idxcl())
    write(unit=fidlog,fmt='(72("-"))')
  endif

  ! Reset charge
  charge=0
  do ispec_=idxcl(),nspec()
    charge=charge+xz(ispec_)*ni(ispec_)
  enddo

  ovrlap=.false.

  particlecount = sum(ni(idxcl():nspec()))
  s6q4r8: if (usegrid) then
    write(unit=fidlog,fmt=*)"Using grid for generating initial conformation"
    ! ------------------------
    ! USE GRID POSITION FINDER
    ! ------------------------

    ! length : the minimum interparticle distance (2 * max(xri(i:1,nspec())))
    length = 0.D0
    do jj=1,nspec()
      length = max(length, xri(jj))
    enddo
    length = length * 2.D0

    ! tubelength : 2 * (zl(4) - zl(2))
    call gridinit(genrgrid, lenblk(), length)

    ! do we have enough gridpoints?
    if (gridsize(genrgrid).le.particlecount) then
       write(unit=fidlog,fmt=*)"Grid size (",gridsize(genrgrid),&
            &") less than number of particles (",particlecount,")"
       stop "Insufficient positions in grid"
    endif
    if (gridsize(genrgrid)/2.le.particlecount) then
      write(unit=fidlog,fmt=*)"WARNING: Grid size is less than twice the number of particles"
    endif
  endif s6q4r8
  write(unit=fidlog,fmt=*)"! Generating ",particlecount," initial particles"

  do ispec_=idxcl(),nspec()
    ntrg_=ni(ispec_)
    write(unit=fidlog,fmt=*)"! Generating ",ntrg_," particles for ",fspc(ispec_)
    do lidx_=1,ntrg_
      itrys_=0 ! counter to stop loop repeating ad infinitum
      ii_=nactv+1
      ovrlap=.true.
      do while (ovrlap)
        ovrlap=.false.
        itrys_=itrys_+1
        if (itrys_.gt.1000) &
          & stop "Unable to generate initial guess in reasonable time"
        if (usegrid) then
          call gridnext(genrgrid, rxinw_, ryinw_, rzinw_)
          r2inw_ = sqrt(sqr(rxinw_) + sqr(ryinw_))
        else
          call cubmov(rxinw_, ryinw_, rzinw_, r2inw_, lenblk())
        endif
        ! Check for overlap with existing particles and build rqqii matrix
        ri=xri(ispec_)
        do jj=1,nactv
          jspec=ispcbk(jj)
          if (jspec.eq.0) cycle
          rii=dispbc(rxinw_, ryinw_, rzinw_, rx(jj), ry(jj), rz(jj))
          if (rii.lt.ri+xri(jspec)) then
            ovrlap=.true.
            exit
          endif
          if (dfeq(xq(ispec_),0.D0).or.dfeq(xq(jspec),0.D0)) then
            rqqii(jj,ii_)=rii
          else
            rqqii(jj,ii_)=xq(ispec_)*xq(jspec)/(2 * rii)
          endif
          rqqii(ii_,jj)=rqqii(jj,ii_)
        enddo
      enddo
      nactv=nactv+1
      rx(ii_)=rxinw_
      ry(ii_)=ryinw_
      rz(ii_)=rzinw_
      r2(ii_)=r2inw_
      ispcbk(ii_)=ispec_
    enddo
  enddo

  if (usegrid) then
    call gridkill(genrgrid)
  endif
  write(unit=fidlog,fmt=*)" Initial configuration made"
  end subroutine genrbk

  ! ------------------------------------------------------------------
  ! GENERATE AN INITIAL MC CONFIGURATION
  ! 
  ! This generates the initial simulation system.  At entry into
  ! this method any structural ions should already have been added
  ! to the system.
  !
  ! Particles are added one-by-one checking for overlap after each
  ! addition, which causes the particle to be rejected.  To avoid
  ! this loop from continuing ad nauseum the number of possible
  ! insertion attempts is limited.
  !
  ! The trial position can be anywhere within the simulation
  ! boundary.  Only particle-particle and particle-system overlap
  ! is used to reject a trial during particle addition; no account
  ! is taken of the energy so it is not even calculated.
  !
  ! This method generates entries for the added ions in the rqqii
  ! array as it progresses so the 'lookup' method does not need
  ! to be called after this method ends for the new ions.  When
  ! structural ions are present the 'lookup' method must be called
  ! before this method.
  !
  ! ENSURE:
  ! * verify that any predefined particle counts fit to the requested 
  !   particle number. (accept N +/- 2*sqrt(N)
  ! ** error if too _many_ predefined particles [ > N + tol ]
  ! ** add particles if too _few_ predifined particles [ < N + tol ]
  ! ** _no_ change if predefined number is within tolerance
  subroutine genrcf(filename)
    use grid
    use geom
    use spec
    implicit none
    character(len=*), intent(in) :: filename

    ! Locals
    double precision :: rzinw_,r2inw_,rxinw_,ryinw_ ! trial data
    integer :: ispec_,igc ! Specie and salt indices 
    integer :: ii_        ! particle global indices
    integer :: ntrg_,lidx_,itrys_ ! counters
    logical :: ovrlap
    double precision :: ri,dda,riasq, ncut
    integer :: jspec,jj
    double precision :: charge ! sum charge
    integer :: nchg ! charge as an integer
    double precision :: zreg_l,zreg_r,rreg_ ! cached region variables

    ! Number of loc_ni to create for each specie
    integer :: particlecount
    integer, dimension(nspcmx) :: loc_ni    ! cache particle numbers
    ! inter-grid spacing
    double precision :: length

    type (tubegrid) :: genrgrid

    if (dbc) then
      if (dbc_level.ge.dbc_check) then
        if (.not.allocated(rqqii)) stop "Error: conf%genrcf called before initialising conf module"
      endif
    endif

    ! count the number of loc_ni to create
    loc_ni = 0
    ! ------- Calculate particle numbers --------------------
    ! (Do not handle any "structural" ions that are the equiv.
    ! specie to a free ion)
    ! Use [spc] = n_spc * tosi / V_tot(spc)
    do igc=1,nsalt()
      loc_ni(isalt(igc))=loc_ni(isalt(igc))+nint(ctargs(igc)*vtotal(isalt(igc))/tosi)
      loc_ni(idxcl())=loc_ni(idxcl())+nint(ctargs(igc)*vtotal(idxcl())*xz(isalt(igc))/tosi)
    enddo

    ! Double check for electroneutrality due to structural ions
    ! and rounding
    charge=0
    do ispec_=1,nspec()
      charge=charge+xz(ispec_)*loc_ni(ispec_)
    enddo
    ! -----recompute N(Cl) if cell is not neutral ------------
    if (.not.dfeq(charge,0.D0)) then
      nchg=nint(charge)
      write(unit=fidlog,fmt=*)"Total cell charge = ",nchg
      write(unit=fidlog,fmt=*)"Cell is not neutral, recomputing chloride (Cl) number"
      write(unit=fidlog,fmt=*)"Old N(Cl) = ",loc_ni(idxcl())
      if (nint(xz(idxcl())).eq.-1) then
        loc_ni(idxcl())=nchg+loc_ni(idxcl())
      else
        stop "Input error: Valency of specie (Cl) should be -1!"
      endif
      write(unit=fidlog,fmt=*)"New N(Cl) = ",loc_ni(idxcl())
      write(unit=fidlog,fmt='(72("-"))')
    endif
    ! ------- Check nionmx is not exceeded -----
    write(unit=fidlog,fmt=*)"Target particle numbers"
    particlecount = 0
    do ispec_=1,nspec()
       if (dbc) then
          if (dbc_level.ge.dbc_check) then
             if (loc_ni(ispec_).lt.0) stop "Error: A negative number of particles for a specie was required"
          endif
       endif
       if (ispec_.ge.idxcl()) then
          particlecount = particlecount + loc_ni(ispec_)
          write(unit=fidlog,fmt='(1X,A2,6X,I6)')fspc(ispec_),loc_ni(ispec_)
       endif
    enddo
    write(unit=fidlog,fmt='(1X,"TOTAL",3X,I6)')particlecount

    ! Reset charge
    charge=0
    do ispec_=1,nspec()
      charge=charge+xz(ispec_)*ni(ispec_)
    enddo
    ovrlap=.false.
    write(unit=fidlog,fmt=*)" Making initial configuration"

    d1c0m3: do ispec_=idxcl(),nspec()
      ! At this point we have the target particle numbers in loc_ni. Now
      ! compare them to existing particle count.
      if (ni(ispec_).ne.0) then
        ! have existing particles (minimum cutoff of 10)
        ncut = max(10,nint(2*sqrt(dble(loc_ni(ispec_)))))
        if (ni(ispec_).gt.loc_ni(ispec_)+ncut) then
          ! error, too many particles defined
          write(fidlog,*)"Specie ",fspc(ispec_)," has ",ni(ispec_)," predifined particles."
          write(fidlog,*)"This is greater than the tolerance (",ncut,") from ",loc_ni(ispec_),"."
          stop "Too many particles in initial configuration"
          loc_ni(ispec_) = 0
        else if (ni(ispec_).gt.loc_ni(ispec_)-ncut) then
          ! ok, predefined particles within tolerance
          !
          ! NOTE: we do not change the particle count if the number of
          ! predefined particles is within tolerance
          loc_ni(ispec_) = 0
          cycle
        else
          ! else too few loc_ni -> calculate how many to generate.
          loc_ni(ispec_) = loc_ni(ispec_) - ni(ispec_)
        endif
      endif
    enddo d1c0m3

    particlecount = sum(loc_ni)

    write(unit=fidlog,fmt=*)"! Generating ",particlecount," initial particles"
    s6q4r9: if (usegrid) then
      write(unit=fidlog,fmt=*)"Using grid for generating initial conformation"
      ! ------------------------
      ! USE GRID POSITION FINDER
      ! ------------------------

      ! length : the minimum interparticle distance (2 * max(xri(i:1,nspec())))
      length = 0.D0
      do jj=1,nspec()
        length = max(length, xri(jj))
      enddo
      length = length * 2.D0

      ! radius : the radius of the tube (rl5)
      ! tubelength : 2 * (zl(4) - zl(2))
      call gridinit(genrgrid, 2 * (zl(4) - zl(2)), rl(5), length)

      ! do we have enough gridpoints?
      if (gridsize(genrgrid).le.particlecount) then
         write(unit=fidlog,fmt=*)"Grid size (",gridsize(genrgrid),&
              &") less than number of particles (",particlecount,")"
         stop "Insufficient positions in grid"
      endif
      if (gridsize(genrgrid)/2.le.particlecount) then
        write(unit=fidlog,fmt=*)"WARNING: Grid size is less than twice the number of particles"
      endif
    endif s6q4r9

    ! ------------------------------------
    ! LOOP OVER NUMBER OF PARTICLES TO ADD
    ! ------------------------------------      
    g6x9b6: do lidx_=particlecount,1,-1

      ! randomly select a specie by selecting
      ! a random particle number (jspec) between 1 and lidx
      ! then finding which specie this belongs to.
      ispec_ = 0
      jspec = ceiling(lidx_ * ranff())
      j0v7p7: do jj = idxcl(), nspec()
        jspec = jspec - loc_ni(jj)
        if (jspec.le.0) then
          ispec_ = jj
          exit j0v7p7
        endif
      enddo j0v7p7
      if (ispec_.eq.0) then
        if (dbc) then
          write(unit=fidlog,fmt=*)"Problem selecting specie:"
          write(unit=fidlog,fmt=*)"Particle count:",lidx_
          write(unit=fidlog,fmt=*)"Chooser:",jspec,jj
        endif
        stop "Problem selecting specie:"
      endif

      d5j5w2: if (.not.usegrid) then
        ! --------------------------
        ! USE RANDOM POSITION FINDER
        ! --------------------------
        ! Cache insertion regions for this specie
        zreg_l=-zreg(4,ispec_)
        zreg_r=zreg(4,ispec_)
        rreg_=rreg(4,ispec_)
        itrys_=0 ! counter to stop loop repeating ad infinitum
      endif d5j5w2

      ii_=nactv+1
      ovrlap=.true.
      j8m8d9: do while (ovrlap)
        ovrlap = .false.
        c9d3k8: if (usegrid) then
          ! ------------------------
          ! USE GRID POSITION FINDER
          ! ------------------------
          call gridnext(genrgrid, rxinw_, ryinw_, rzinw_)

          ! adjust z to be outside zl(2)
          if (0.D0.lt.rzinw_) then
            rzinw_ = rzinw_ + zl(2)
          else
            rzinw_ = rzinw_ - zl(2)
          endif

          r2inw_=sqrt(sqr(rxinw_) + sqr(ryinw_))

          ! NO wall test as grid should ensure that would pass
          if (dbc) then
            if (dbc_level.ge.dbc_check) then
              call wall(ispec_,rzinw_,r2inw_,ovrlap)
              if (ovrlap) then
                write(unit=fidlog,fmt=*)"particle on grid overlaps wall"
                write(unit=fidlog,fmt=*)"position:   ",rxinw_, ryinw_, rzinw_, r2inw_
                stop "particle on grid overlaps wall"
              endif
            endif
          endif
        else c9d3k8
          ! --------------------------
          ! USE RANDOM POSITION FINDER
          ! --------------------------
          itrys_=itrys_+1
          if (itrys_.gt.1000) then
            write(unit=fidlog,fmt=*)'Stopping after ',itrys_,' attempts, created only '&
                 &,lidx_,' of ',ntrg_,' particles for ',fspc(ispec_)
            stop "Unable to generate initial guess in reasonable time"
          endif
          call jmpmov(rzinw_, r2inw_, zreg_l,zreg_r,rreg_)
          call wall(ispec_,rzinw_,r2inw_,ovrlap)
        endif c9d3k8
        h3d0v2: if (.not.ovrlap) then
          call jmpfin(rxinw_, ryinw_, r2inw_)
          ! Check for overlap with existing particles
          ri=xri(ispec_)
          i5c2k9: do jj=1,nactv
            jspec=ispcbk(jj)
            if (jspec.eq.0) cycle i5c2k9
            dda=sqr(ri+xri(jspec))
            riasq=sqr(rxinw_-rx(jj))+sqr(ryinw_-ry(jj))+sqr(rzinw_-rz(jj))
            d2k6t6: if (riasq.lt.dda) then
              ovrlap=.true.
              exit i5c2k9
            endif d2k6t6
            if (dfeq(xq(ispec_),0.D0).or.dfeq(xq(jspec),0.D0)) then
              ! handle uncharged species specially
              rqqii(jj,ii_)=dsqrt(riasq)
            else
              rqqii(jj,ii_)=xq(ispec_)*xq(jspec)/(2 * dsqrt(riasq))
            endif
            rqqii(ii_,jj)=rqqii(jj,ii_)
          enddo i5c2k9
        endif h3d0v2
      enddo j8m8d9
      nactv=nactv+1
      rx(ii_)=rxinw_
      ry(ii_)=ryinw_
      rz(ii_)=rzinw_
      r2(ii_)=r2inw_
      ispcbk(ii_)=ispec_
      call addreg(ispec_,ii_,rzinw_)
      loc_ni(ispec_) = loc_ni(ispec_) - 1
      ni(ispec_)=ni(ispec_)+1
      !     write(unit=fidlog,fmt=*)ii_,ispec_,ni,rxinw_,ryinw_,rzinw_
    enddo g6x9b6

    if (usegrid) then
      call gridkill(genrgrid)
    endif
    write(unit=fidlog,fmt=*)" Initial configuration made"
    call writcf(filename,0)

  end subroutine genrcf

  ! --------------------------------------------------
  ! Get an available global index
  !
  ! As we now use a compact array and a deletion list the next
  ! available index for particle insertion must be obtained from
  ! checking the deletion list as well as the array's active size.
  ! This method performs that function and returns the lowest
  ! numbered global index where a new particle can be defined.  It
  ! updates the deletion list if it was not empty. It also updates
  ! the 'ispcbk' and 'ni' data sets.
  !
  ! @param ispec : specie of new particle
  integer function idxget(ispec)
  use spec
  implicit none
  integer, intent(in) :: ispec
  if (dbc) then
    if(ispec.lt.1.or.ispec.gt.nspec()) stop "Specie index out of range"
  endif

  t9t5p8: if (ndel.eq.0) then
    nactv=nactv+1
    if (nactv.gt.ntotsz) stop "Too many particles being created." 
    idxget=nactv
  else t9t5p8
    idxget=idelst(ndel)
    ndel=ndel-1
  endif t9t5p8
  ispcbk(idxget)=ispec
  ni(ispec)=ni(ispec)+1

  end function idxget

  ! --------------------------------------------------
  ! release a global index
  !
  ! As we now use a compact array and a deletion list we 
  ! need to manage deleting a particle.  This method
  ! handles the addition of the unused index into the
  ! deletion list, which we keep in sorted order.  The
  ! 'ispcbk' and 'ni' data sets are also updated here.
  !
  ! @param ispec : specie of new particle
  ! @param ii : index of old particle
  subroutine idxrel(ispec,ii)
  use spec
  implicit none
  integer, intent(in) :: ii, ispec
  if (dbc) then
      if (ispec.lt.1.or.ispec.gt.nspec()) stop "Specie index out of range"
      if (ii.lt.1.or.ii.gt.nactv) stop "Particle index out of range"
      if (ispec.ne.ispcbk(ii)) stop "Specie of index does not match argument"
  endif
  ! Handle case where particle is last one in array
  if (ii.eq.nactv) then
    nactv=nactv-1
  else
    ndel=ndel+1
    if (ndel.gt.nionmx) stop "Too many particles being destroyed"
    idelst(ndel)=ii
    if (ndel.gt.1) call isort(idelst,ndel,.false.)
  endif
  ! SETTING ISPCBK to zero to indicate invalid particle
  ispcbk(ii)=0
  ni(ispec)=ni(ispec)-1
  end subroutine idxrel

  ! ------------------------------------------------------------------
  ! Generate lookup data for particles
  !
  ! This method generates the rqqii array for particles whose
  ! position information has been read in.  This method is best
  ! called directly after the ions have been defined.  For example
  ! if structural ions are present in the input then this method
  ! should be called before a generator method is called to add
  ! the rest of the particles.  Similarly this method is called
  ! within 'readcf' after all the particles have been read in.
  !
  ! ENSURE:
  ! * NOTE: calling rfgeom before rfconf/lookup is required. 
  ! * all particles defined so far do not overlap the wall OR other
  !   particles.
  ! * 'mobile' particles are within their allowed spheres
  subroutine lookup
  use spec
  use geom
  implicit none

  ! LOCALS
  double precision rxi,ryi,rzi  ! particle coords
  double precision rijsq,ri,ddi ! inter-particle or patch data
  integer ii,jj       ! dummy vars or particle indices
  integer ispec,jspec ! particle indices
  logical :: ovrlap
  ovrlap = .false.

  !  -------------------------
  !  Handle particle positions
  !  -------------------------
  e7v0j8: do ii=1,nactv
    ispec=ispcbk(ii)
    if (ispec.eq.0) cycle e7v0j8
    rxi = rx(ii)
    ryi = ry(ii)
    rzi = rz(ii)
    ri  = xri(ispec)
    rqqii(ii,ii) = 0.D0

    if (localized(ispec)) then
      ! if particle is 'mobile' type, check its distance from its origin
      ! Here we only report the error.
      call mobchk(ii, rxi, ryi, rzi, ovrlap)
      a1s8q1: if (ovrlap) then
        write(fidlog,*)"WARNING: particle ",ii," (",fspc(ispec),") is ",&
                 &sqrt(sqr(rxi-rsx(ii)) + sqr(ryi-rsy(ii)) + sqr(rzi-rsz(ii))),&
                 &" from its origin, beyond allowed value of ",sqrt(rsr(ii))
        if (.not.relaxmob) then
          relaxmob=.true.
          write(fidlog,*)"WARNING: Removing mobile/flexible cut-off during thermalisation"
        endif
      endif a1s8q1
      ovrlap = .false.
    endif

    ! Determine distance to all particles higher in array than us
    e1j1j4: do jj=ii+1,nactv
      jspec=ispcbk(jj)
      if (jspec.eq.0) cycle e1j1j4
      ! Calculate min distance squared
      ddi=sqr(ri+xri(jspec))
      ! Calculate distance
      rijsq=sqr(rxi-rx(jj))+sqr(ryi-ry(jj))+sqr(rzi-rz(jj))
      a1s8q0: if (rijsq.lt.ddi) then
        write(fidlog,*)"Particle ",ii," (",fspc(ispec),") only ",sqrt(rijsq),"(",sqrt(ddi),&
                 &") from ",jj," (",fspc(jspec),")"
        stop "Overlap in initial configuration"
      endif a1s8q0
      if (dfeq(xq(ispec),0.D0).or.dfeq(xq(jspec),0.D0)) then
        rqqii(jj,ii)=sqrt(rijsq)
      else
        rqqii(jj,ii)=xq(ispec)*xq(jspec)/(2 * sqrt(rijsq))
      endif
      rqqii(ii,jj)=rqqii(jj,ii)
    enddo e1j1j4
  enddo e7v0j8
  end subroutine lookup

  ! ------------------------------------------------------------------
  ! Check the position of mobile/flexible particles at end of thermalisation
  !
  subroutine chkmob
  use spec
  implicit none

  ! LOCALS
  integer ii    ! dummy vars or particle indices
  integer ispec ! specie indices
  logical :: ovrlap
  relaxmob=.false.
  c4w8m7: do ii=1,nactv
    ispec=ispcbk(ii)
    ! exit at first 'free' particle
    if (isfree (ispec)) return
    ! cycle if particle is not mobile
    if (localized(ispec)) then
      ! if particle is 'mobile' or 'flexible' type, check its distance from its origin
      call mobchk(ii, rx(ii), ry(ii), rz(ii), ovrlap)
      d6p9s4: if (ovrlap) then
        write(fidlog,*)"Particle ",ii," (",fspc(ispec),") is ",&
             &sqrt(sqr(rx(ii)-rsx(ii)) + sqr(ry(ii)-rsy(ii)) + sqr(rz(ii)-rsz(ii))),&
             &" from its origin, beyond allowed value of ",sqrt(rsr(ii))
        stop "Bad mobile/flexible ion position"
      endif d6p9s4
    endif
  enddo c4w8m7
  end subroutine chkmob

  ! -------------------------------------------------------------
  ! READ CONF CHECKPOINT (saved using 'writcf')
  ! 
  ! readcf and writcf save particle specific data, all
  ! other information must be read in some other way.
  ! In particular details of the simulation system are
  ! not saved so changes to the system geometry will
  ! invalidate the contents.
  ! 
  subroutine readcf(filename)
  use spec
  implicit none
  character(len=*), intent(in) :: filename

  ! LOCALS
  integer ii,ispec,ispx   ! indices
  integer imap(nspcmx) ! map old specie index to new index
  character(2) fspnam
  
  ! Check arrays are allocated
  if (dbc) then
    if (dbc_level.ge.dbc_check) then
      if (.not.allocated(rqqii)) stop "Error: conf%readcf called before initialising conf module"
    endif
  endif
  ! reset variables
  call x2a4c2

  p7u8e0: do ispec=1,nspcmx
    imap(ispec)=ispec
  enddo p7u8e0

  open(unit=fidcnf,file=filename,action="read")
    read(unit=fidcnf,fmt='(I5)') nactv

  !   Read specie indices and create mapping
    u2p2o1: do ispx=1,nspec()
      read(unit=fidcnf,fmt=*) ii,fspnam
      m8u5e8: do ispec=1,nspcmx
        a4n2z5: if (fspnam.eq.fspc(ispec)) then
          imap(ii)=ispec
          exit m8u5e8
        endif a4n2z5
      enddo m8u5e8
    enddo u2p2o1
    read(unit=fidcnf,fmt=*)
  !   Read particles
    w1s2h4: do ii=1,nactv
      read(unit=fidcnf,fmt=*)rx(ii),ry(ii),rz(ii),ispcbk(ii)
      ispcbk(ii)=imap(ispcbk(ii))
      ni(ispcbk(ii))=ni(ispcbk(ii))+1
      r2(ii)=dsqrt(sqr(rx(ii)) + sqr(ry(ii)))
      call addreg(ispcbk(ii),ii,rz(ii))
    enddo w1s2h4
  close(unit=fidcnf)

  ! Reset charge
  charge=0
  do ispec=1,nspec()
    charge=charge+xz(ispec)*ni(ispec)
  enddo

  ! Fill up index arrays and lookup tables
  call lookup
  contains

  subroutine x2a4c2
    ! Reset the conf variables
    rqqii = 0
    rx = 0
    ry = 0
    rz = 0
    r2 = 0

    p4c6h2 = 0
    ispcbk = 0
    idelst = 0
    charge=0
    ni=0
    nreg=0
    nin=0
    nactv=0
    ndel=0
  end subroutine x2a4c2

  end subroutine readcf

  ! ------------------------------------------------------------------
  ! FINALISE INITIALISATION
  !
  ! This method is called to finalise the process of initialising
  ! the conf module when an input file is being read.  This method
  ! is one of the input finalisation method that is called in a
  ! particular order in 'channel%readin' method.  To ensure correct
  ! initialisation the order these methods are called is critical.
  !
  ! No consideration is given to the possibility this method is
  ! called more than once in a program run.
  !
  ! The key activities are:
  !   - check, allocate and zero in-use arrays
  !   - copy particle position information present in input file from
  !        the 'spec' module.
  !   - tell the 'spec' module that the input position information
  !        it has is no longer needed.
  !   - use 'lookup' to generate inter-particle information (for 
  !        particles copied above.)
  !   - compute the initial volume profile for structural ions in the
  !        filter.
  !
  subroutine rfconf
  use spec
  use geom
  implicit none
  double precision :: x_,y_,z_
  integer :: ispec ! loop indices
  integer :: idx ! structural ion local index
  integer :: ii ! global index
  integer :: stat ! indicate if memory allocation succeeded
  ntotsz=next64(nionmx + ntargt()*2) ! count of all structural ions < nionmx
  if (.not.allocated(rqqii)) then
    allocate(rqqii(ntotsz,ntotsz), STAT=stat)
    if (stat.ne.0) stop "Memory allocation of RQQII failed"
    rqqii=0
    allocate(rx(ntotsz),ry(ntotsz),rz(ntotsz),r2(ntotsz), STAT=stat)
    if (stat.ne.0) stop "Memory allocation of failed"
    rx=0
    ry=0
    rz=0
    r2=0
    allocate(ispcbk(ntotsz), STAT=stat)
    if (stat.ne.0) stop "Memory allocation of ISPCBK failed"
    ispcbk=0
    allocate(idelst(nionmx), STAT=stat)
    if (stat.ne.0) stop "Memory allocation of IDELST failed"
    idelst=0
    allocate(p4c6h2(ntotsz,nrgnmx), STAT=stat)
    if (stat.ne.0) stop "Memory allocation of P4C6H2 failed"
    p4c6h2=0
    nreg=0
    nin=0
    ni=0
    nactv=0
    ndel=0
  endif

  ! copy any predifined positions from spec module
  write(unit=fidlog,fmt='(72("-"))')
  write(unit=fidlog,fmt=*)"Initial position of structural ions (*mobile update radii and centre)"
  write(unit=fidlog,fmt='(72("-"))')
  write(unit=fidlog,fmt='(A3,2X,A3,7(1X,A8))')"Spc","Idx","X","Y","Z","Radii","Xc","Yc","Zc"
  write(unit=fidlog,fmt='(72("-"))')
  do ispec=1,nspec()
    if (input_count(ispec).gt.0) then
      do idx=1,input_count(ispec)
        call struks(ispec,idx,x_,y_,z_)
        nactv=nactv+1
        rx(nactv)=x_
        ry(nactv)=y_
        rz(nactv)=z_
        r2(nactv)=dsqrt(sqr(x_)+sqr(y_))
        if (localized(ispec)) then
          write(unit=fidlog,fmt='(A3,2X,I3,7(1X,F8.2))')fspc(ispec),nactv,rx(nactv),&
               ry(nactv),rz(nactv),sqrt(rsr(ispec,idx)),rsx(ispec,idx),rsy(ispec,idx),rsz(ispec,idx)
        else
          write(unit=fidlog,fmt='(A3,2X,I3,4(1X,F8.2))')fspc(ispec),nactv,rx(nactv),ry(nactv),rz(nactv)
        endif
        ispcbk(nactv)=ispec
        call addreg(ispec,nactv,z_)
        ni(ispec)=ni(ispec)+1
      enddo
      ! Tell spec module it no longer needs to store xyz for this specie
      call deltmp(ispec)
    end if
  enddo
  write(unit=fidlog,fmt='(72("-"))')
  write(unit=fidlog,fmt=*)"(*)The given radii is the distance a 'mobile'/'flexible' type particle can"
  write(unit=fidlog,fmt=*)"move from the Xc,Yc,Zc position (may be the same as the starting position)."
  write(unit=fidlog,fmt='(72("-"))')

  ! Fill up index arrays and lookup tables for predifined particles
  call lookup

  do ii=1,nactv
    ispec=ispcbk(ii)
    if (ispec.eq.0) continue
    if (isfree(ispec)) exit
    if (localized(ispec)) then
      write(unit=fidlog,fmt='(A3,2X,I3,7(1X,F8.2))')fspc(ispec),ii,rx(ii),&
             ry(ii),rz(ii),sqrt(rsr(ii)),rsx(ii),rsy(ii),rsz(ii)
    else
      write(unit=fidlog,fmt='(A3,2X,I3,4(1X,F8.2))')fspc(ispec),ii,rx(ii),ry(ii),rz(ii)
    endif
  enddo

  ! Calculate volume ratios 
  call volslb(fidlog,10,0)

  end subroutine rfconf

  ! --------------------------------------------------------------
  ! OCCLUDED VOLUME
  ! 
  ! Compute volume fractions for spheres in the filter region.
  !
  ! This model works by a combination of restricted volume and
  ! high charge density in the filter region.  Previous versions
  ! of the program did not limit the structural ions beyond
  ! restricting them to the filter region.  This version introduces
  ! the ability to reduce the mobility of structural ions.  This
  ! makes knowledge about the  profile of free space in the filter
  ! interesting.  This method calculates the percent volume occupied
  ! by structural ions in a set of intervals in filter.
  !
  ! This information is saved into a file 'res/occvol.num.dat' and
  ! written to the input file descriptor.
  subroutine volslb(fid,nival,istep)
  use geom
  use spec
  implicit none
  integer, intent(in) :: fid
  integer, intent(in) :: nival  ! number of intervals to use
  integer, intent(in) :: istep  ! current save step
  double precision :: left,right,width,sumsq,val
  integer :: idx ! interval index
  ! total filter region
  write(unit=fid,fmt=*)
  write(unit=fid,fmt='(A)')'Percentage of filter region volume occupied by structural particles'
  write(unit=fid,fmt='(72("-"))')
  write(unit=fid,fmt='(2(A9),A11)')"Left","Right","Occupied"
  write(unit=fid,fmt='(30("-"))')
  if (istep.eq.0)  then
    open(unit=fidvol,file="res/occvol."//firun//".dat")
    ! WRITE FILE HEADER
    write(unit=fidvol,fmt='("# UUID ",A32)')fuuid
    write(unit=fidvol,fmt='(A)')'# label step  5   15  25  35  45  55  65  75  85  95  aver  var'
    write(unit=fidvol,fmt='(A)')'# units count %   %   %   %   %   %   %   %   %   %   %     %'
  else
    open(unit=fidvol,file="res/occvol."//firun//".dat",position="APPEND")
  endif 
  write(unit=fidvol,fmt='(1X,I8)',advance="NO")istep
  sumsq=0.D0
  right=0.D0
  val=0.D0
  ! slices
  left=-zlimit()
  width=(2*zlimit()/nival)
  do idx=1,nival
    ! left=-zocc+idx*width
    right=left+width
    val=sphvol(left,right,.true.)/(pi*sqr(rl(1))*(width))
    sumsq=sumsq+val*val
    write(unit=fid,fmt='(2X,F7.2,2X,F7.2,2X,F7.2," %")')left,right,100*val
    write(unit=fidvol,fmt='(2X,F13.6)',advance="NO")val
    left=right
  enddo

  ! total
  left=-zlimit()
  right=zlimit()
  width=2*zlimit()
  val=sphvol(left,right,.true.)/(pi*sqr(rl(1))*(width))
  write(unit=fid,fmt='(30("-"))')
  write(unit=fid,fmt='("  [",F7.2,":",F7.2,"]",1X,F7.2," %")')left,right,100*val
  write(unit=fidvol,fmt='(2X,F13.6)',advance="NO")val

  ! compute sample standard deviation
  val=sqrt((sumsq-nival*sqr(val))/(nival-1))
  write(unit=fid,fmt='(1X,A18,1X,F7.2," %")')"Standard Deviation",val*100
  write(unit=fidvol,fmt='(2X,F13.6)')val
  write(unit=fid,fmt='(72("-"))')
  close(unit=fidvol)
  contains

  ! -------
  ! calculate the volume contribution of spheres within z-interval
  !
  ! This internal method checks all structural ions and calculates
  ! the volume it occupies within an interval of the filter.  This
  ! includes calculating the part of the volume inside the interval
  ! for those spheres that overlap the interval boundaries.
  double precision function sphvol(left,right,ignorefree)
  use geom
  ! assumes that spheres do not overlap with
  ! each other and with curved cylinder walls
  implicit none
  double precision, intent(in) :: left,right
  logical, intent(in) :: ignorefree

  integer :: gidx ! particle index
  integer :: ispec ! specie index
  double precision :: zcentr,spcrad,partvl,leftmx,rigtmx,h

  sphvol=0.D0

  ! foreach particle
  do gidx=1,nactv
    ispec=ispcbk(gidx) ! sphere specie

    ! ignore unsused particles
    if (ispec.eq.0) cycle

    ! ignore free species?
    if (ignorefree.and.isfree(ispec)) cycle
!    if (ignorefree.and.flexible(ispec)) cycle

    zcentr=rz(gidx)    ! sphere centre 
    spcrad=xri(ispec)  ! sphere radii
    leftmx=zcentr-spcrad
    rigtmx=zcentr+spcrad
    partvl=0.D0
    h=0.D0
 
    if ((rigtmx.gt.left).and.(leftmx.lt.right)) then
      ! sphere in region, add total volume
      partvl=4*pi*(spcrad**3)/3 ! volume of sphere in cylinder

      ! check for partial overlap
      if (leftmx.lt.left) then
        ! overlap to left
        if (zcentr.gt.left) then
          ! overlap from inside: delete extra volume
          h=left-leftmx
          partvl=partvl-pi*sqr(h)*(3*spcrad-h)/3
        else
          ! overlap from outside
          h=rigtmx-left
          partvl=pi*sqr(h)*(3*spcrad-h)/3
        endif
      endif
      ! note: we use separate if for right-side overlap so we
      ! handle case where a sphere overlaps both boundaries.
      if (rigtmx.gt.right) then
        ! overlap to right
        if (zcentr.lt.right) then
          ! overlap from inside: delete extra volume
          h=rigtmx-right
          partvl=partvl-pi*sqr(h)*(3*spcrad-h)/3
        else
          ! overlap from outside
          h=right-leftmx
          partvl=pi*sqr(h)*(3*spcrad-h)/3
        endif
      endif
    endif
  
    sphvol=sphvol+partvl
  enddo
  endfunction sphvol

  end subroutine volslb

  ! ------------------------------------------------------------------
  ! MANAGE indreg: move a particle between regions.
  !
  ! Update a particle's region settings after moving.
  !
  ! This method only updates settings where a particle has moved into
  ! or out of a region.
  !
  ! @param ispec : Specie index
  ! @param ii : global index
  ! @param rzold : z-coord of old position
  ! @param rznew : z-coord of new position
  subroutine setreg(ispec,ii,rzold,rznew)
  use spec
  use geom
  implicit none
  integer, intent(in) :: ispec,ii
  double precision, intent(in) :: rzold,rznew
  ! LOCALS
  integer ireg 
  if (dbc) then
    if (dbc_level.ge.dbc_require) then
      if (ispec.lt.1.or.ispec.gt.nspec()) stop "invalid specie index"
      if (ii.lt.1.or.ii.gt.nactv) stop "invalid particle index"
      if (ispec.ne.ispcbk(ii)) stop "invalid particle index"
    endif
  endif
  z8h8b0: do ireg=1,3
    if ((abs(rzold).gt.zreg(ireg,ispec)).and. &
      &  (abs(rznew).le.zreg(ireg,ispec))) &
      & call ce1c59(ireg,ispec,ii)

    if ((abs(rznew).gt.zreg(ireg,ispec)).and. &
      & (abs(rzold).le.zreg(ireg,ispec))) &
      & call ed73b5(ireg,ispec,ii)
  enddo z8h8b0
  end subroutine setreg

  ! -------------------------------------------------------------
  ! Write configuration in input file format
  !
  ! This method separates the particles into individual species.
  ! Once separated, it passes the array of particles for each 
  ! specie to spec@ecspec to write out as per the input file.  
  !
  subroutine wrconf(fid)
  use spec
  implicit none
  integer, intent(in) :: fid
 
  type savecf
    integer :: sz_
    double precision, dimension(:,:), allocatable :: xyzs_
  end type savecf

  ! LOCALS
  type(savecf), dimension(:), allocatable :: specs_
  integer :: ispec ! specie
  integer :: ii   ! Particle indices
  integer :: stat ! indicate status of memory allocation.deallocation

  if (nspec().le.0) stop "Number of species is zero"
  allocate(specs_(nspec()), STAT=stat)
  if (stat.ne.0) stop "Memory allocation of SPECS_ failed"
  do ispec=1,nspec()
    if (ni(ispec).ne.0) then
      if (localized(ispec)) then
        allocate(specs_(ispec)%xyzs_(7,ni(ispec)), STAT=stat)
        if (stat.ne.0) stop "Memory allocation of SPECS_%XYZS failed"
      else
        allocate(specs_(ispec)%xyzs_(3,ni(ispec)), STAT=stat)
        if (stat.ne.0) stop "Memory allocation of SPEC_%XYZS failed"
      endif
    else
      write(*,*)"Specie ",fspc(ispec)," has no particles."
      allocate(specs_(ispec)%xyzs_(0,0), STAT=stat)
      if (stat.ne.0) stop "Memory allocation of SPECS_%XYZS(0) failed"
    endif
    specs_(ispec)%sz_=0
  enddo
  ! get xyz (etc) data for each specie
  do ii=1,nactv
    ispec=ispcbk(ii)
    if (0.ne.ispec) then
      specs_(ispec)%sz_ = specs_(ispec)%sz_ + 1
      specs_(ispec)%xyzs_(1,specs_(ispec)%sz_) = rx(ii)
      specs_(ispec)%xyzs_(2,specs_(ispec)%sz_) = ry(ii)
      specs_(ispec)%xyzs_(3,specs_(ispec)%sz_) = rz(ii)
    endif
  enddo
  ! write out each specie
  do ispec=1,nspec()
    call ecspec(fid,ispec,specs_(ispec)%xyzs_,specs_(ispec)%sz_)
    deallocate(specs_(ispec)%xyzs_, STAT=stat)
    if (stat.ne.0) stop "Memory deallocation of SPECS_%XYZS failed"
  enddo
  deallocate(specs_, STAT=stat)
  if (stat.ne.0) stop "Memory deallocation of SPECS_ failed"

  end subroutine wrconf

  ! -------------------------------------------------------------
  ! CHECKPOINT CONF
  !
  ! readcf and writcf only save particle specific data, all other
  ! information must be read in some other way.
  !
  ! This method is designed for the purpose of generating a snapshot
  ! of the particle configuration.  It can be used to save then
  ! read a configuration within the same program run.  It can also
  ! be used to generate a set of particle configuration files that
  ! can be analysed in a separate step.
  subroutine writcf(filename,iter)
  use spec
  implicit none
  character(len=*), intent(in) :: filename
  integer, intent(in) :: iter
  character(8) :: stepnm

  ! LOCALS
  integer ii   ! Particle indices

  if (debug) write(unit=fidlog,fmt='(1("! ",32(1x,I1)))')ispcbk(1:nactv)
  if (byiter) then
    write(stepnm,'(".",I0.7)')iter
  else
    stepnm(1:8)=' '
  endif
  open(unit=fidcnf,file=trim(filename//stepnm),action="write")
    write(unit=fidcnf,fmt='(I5,4X,I8)') nactv-ndel,iter

  !   Write specie indices in case they change
    d5a4r1: do ii=1,nspec()
      write(unit=fidcnf,fmt='(I2,3X,A2)') ii,fspc(ii)
    enddo d5a4r1
    write(unit=fidcnf,fmt=*)

    k8g0b7: do ii=1,nactv
      if (0.ne.ispcbk(ii)) &
        & write(unit=fidcnf,fmt='(3(D24.16,1X),I2)')rx(ii),ry(ii),rz(ii),ispcbk(ii)
    enddo k8g0b7
  close(unit=fidcnf)

  end subroutine writcf

end module conf

