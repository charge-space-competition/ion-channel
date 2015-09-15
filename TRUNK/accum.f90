
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
! Statistics Module
! 
! -------------------------------------------------------------
! Methods:
! accmlt
!    accumulate statistical data that requires no calculation.
! hist
!    process MC system data into histogram bins.
! rdaccu
!    read 'accum' input file section
! rfaccu
!    finalise initialisation after reading the input file
! saves
!    output a digest of the statistics
! zeroac
!    reset the statistic counters
! -------------------------------------------------------------
! Constants
! nzgmx,nrgmx : histrogram bins in z and r directions
!
module accum
use const
use rnstat
implicit none
private
  ! maximum number of histogram bins in r direction
  integer, public, parameter :: nrgmx=2048

  ! maximum number of particles in filter region we count for
  ! generating the co-occupancy matrix 'aocc'
  integer, public, parameter :: noccmx=4

  ! Statistics relating to mobile ions (average: ~/athist)
  double precision, dimension(nionmx), private :: amobdl, amobvr, amobdx, amobdy, amobdz

  ! Sum of h vector from patch module (average: ah/ataccu)
  double precision, private, allocatable, dimension(:) :: ah

  ! co-occupancy matrix - record what sets of particles
  ! are in the filter region at the same time
  double precision, private, dimension(0:noccmx,0:noccmx,0:noccmx,0:noccmx,3) :: aocc

  ! Histogram of particles in each region
  !
  ! (merged 'an' into anin as anin(ibulk,*)
  type (hist2array), private :: anin

  ! Sum of per-specie particle number in bulk region (average: abulk/athist)
  type(runstat), private, dimension(nspcmx) :: abulk

  ! Sum of charges (average: acharge/ataccu)
  type (runstat), private :: acharge

  ! Count of calls to 'hist' routine; used to determine average
  ! of statistics collected in 'hist' routine
  integer (kind=Cntr_K), private :: athist

  ! Count of calls to 'accmlt' routine; used to determine
  ! average of statistics collected in 'accmlt' routine
  integer (kind=Cntr_K), private :: ataccu

  ! Width of bins (both axial and radial) in 2D 'gin' distribution
  ! and 3D 'gxyz' and 'grtz' distributions.
  double precision, private :: drg

  ! The axial outer-most point of the 2D 'gin' distribution
  double precision, private :: zinlft

  ! The radial outer-most point of the 2D 'gin distribution
  double precision, private :: rinup

  ! Width of histogram bins in the inter-particle radial distribution
  double precision, private :: drdf

  ! Width to use for the inner-most occupancy calculations (defaults to zlimit)
  double precision, private :: zocc=0.D0

  ! volumes for the 2D 'gin' histogram bins (independent of z and specie)
  double precision, private, allocatable, dimension(:) :: ginvol

  ! The 2D 'gin' distribution
  type (hist3array), private :: gin

  ! ------------------------------------------------------------------------
  ! 3D distributions

  ! The 3D 'gxyz' and 'grtz' distributions.
  type (hist4array), private :: gxyz, grtz

  ! Index of the particle to use as the theta zero
  integer, private :: grtz_zero

  ! Width of bins (both axial and radial) in 3D 'gxyz' and 'grtz' distributions.
  double precision, private :: d3_width

  ! The radial outer-most point of the 3D 'gxyz'/'grtz' distributions
  double precision, private :: d3_rmax, d3_zmax

  ! The number of radial bins in 3D 'dxyz'/'drtz' distributions
  integer, private :: d3_zgrid_size, d3_rgrid_size

  ! volumes for the 3D 'grtz' histogram bins (independent of z and specie)
  double precision, private, allocatable, dimension(:) :: d3_vol


  ! The counts (numerator) for the inter-particle radial distribution
  ! INDEX (ireg, jspec, ispec)
  type (hist1array), private, allocatable, dimension (:,:,:) :: rdfhist

  ! The z-axial 'gz' distribution data
  type (hist2array), private :: gzhist

  ! The number of in-use bins in the 'gz' distribution
  integer, private :: nrg

  ! The number of radial bins in 2D 'gin' distribution
  integer, private :: nrgr

  ! The number of axial bins in 2D 'gin' distribution
  integer, private :: nrgz

  ! The number of active bins in the inter-particle radial distribution
  integer, private :: nrdf

  ! The number of times we have saved a digest (ie called 'saves')
  integer, private :: ksub=0

  ! [INPUT] How often the statistic data is saved.
  integer, public :: isave

  ! [INPUT] calgin: compute 2D (z,r) profiles
  logical, private :: calgin=.false.

  ! [INPUT] calrdf: compute inter-particle radial distribution
  logical, private :: calrdf=.false.

  ! [INPUT] calacc: show move acceptance ratios
  logical, private :: calacc=.false.

  ! [INPUT] calmob: calculate information about 'mobile' ions
  logical, private :: calmob=.false.

  ! fluctuation variable sum(N_i * N_j) used in iterat
  double precision, dimension(nspcmx,nspcmx) :: ninj
  ! ----------------------------------------------------------------------
  !
  ! PUBLIC ROUTINES AND DATA
  !
  ! INDEPENDENT METHOD
  !   gzaver
  !
  ! ----------------------------------------------------------------------

  public :: accmlt, thermal_accmlt, bulk_accmlt, iterat, ecaccu, gz_av, hist
  public :: concentration_report, bulk_concentration, rdaccu, rfaccu
  public :: reset_wasdel, saves, zeroac

contains
  
  ! -------------------------------------------------------------
  ! Accumulate simple collective statistics regardless of simulation
  ! state.
  !
  subroutine base_accmlt
    use conf
    use simstate
    implicit none
    ataccu=ataccu+1
    if (do_electrostatic()) then
      call rs_push(acharge, charge)
    endif
  end subroutine base_accmlt

  ! -------------------------------------------------------------
  ! Accumulate simple collective statistics for main simulation
  !
  !  This routine collects filter occupancy data from 'nin' arrays.
  !  It also saves the 'h' vector and charge.
  subroutine accmlt
  use conf
  use patch
  use spec
  use simstate
  implicit none

  ! LOCALS
  integer ireg,ispec  ! loop counters for regions,species,patches
  integer i1,i2,i3,i4 ! histogram indices
  call thermal_accmlt

  h3c4v2: do ireg=izlim,ichan
    i1=min(nin(ireg,idxcl()),noccmx)
    i2=min(nin(ireg,idxcl()+1),noccmx)
    i3=min(nin(ireg,idxcl()+2),noccmx)
    i4=min(nin(ireg,idxcl()+3),noccmx)
    aocc(i4,i3,i2,i1,ireg)=aocc(i4,i3,i2,i1,ireg)+1
  enddo h3c4v2
  l0z8a5: do ispec=1,nspec()
    call hist_push (anin, izlim, ispec, nin(izlim,ispec))
    call hist_push (anin, ifilt, ispec, nin(ifilt,ispec))
    call hist_push (anin, ichan, ispec, nin(ichan,ispec))
    call hist_push (anin, ibulk, ispec, ni(ispec))
  enddo l0z8a5

  if (do_electrostatic().and..not.is_homogeneous()) ah(1:npatch)=ah(1:npatch)+h(1:npatch)
  call hist_end_push (anin)
  end subroutine accmlt

  ! --------------------------------------------------
  ! Accumulate statistics for thermalisation and main
  ! part of simulation.
  !
  !  This lesser accmlt is intended for use in the equilibration stage
  !  where the chemical potentials are being adjusted.
  subroutine thermal_accmlt
  use conf
  use geom
  use spec
  implicit none

  ! LOCALS
  integer :: ii     ! loop indices
  double precision :: zabs ! absolute z value
  integer, dimension(nspcmx) :: counter
  call base_accmlt

  counter = 0
  do ii=1,nactv
    if (ispcbk(ii).lt.idxcl()) cycle
    zabs=abs(rz(ii))
    if (zabs.ge.zbulk1.and.zabs.le.zbulk2) then
      if (r2(ii).le.rbulk) counter(ispcbk(ii))=counter(ispcbk(ii)) + 1
    endif
  enddo
  do ii=idxcl(),nspec()
    call rs_push(abulk(ii),dble(counter(ii)))
  enddo
  end subroutine thermal_accmlt

  ! ------------------------------------
  ! accumulate base statistics during bulk simulation
  !
  ! This is a reduced version of accmlt.  It only records the
  ! charge and total particle numbers.
  subroutine bulk_accmlt(cpmeth)
  use conf
  use spec
  use trial
  implicit none
  integer, intent(in) :: cpmeth

  integer :: ispec,jspec
  
  call base_accmlt

  do ispec=idxcl(),nspec()
    call rs_push(abulk(ispec),dble(ni(ispec)))
  enddo
  if (cpmeth.eq.malas2) then
    do ispec=idxcl(),nspec()
      do jspec=idxcl(),nspec()
        ninj(ispec,jspec) = ninj(ispec,jspec) + ni(ispec)*ni(jspec)
      enddo
    enddo
  endif

  end subroutine bulk_accmlt
  
  ! ------------------------------------
  ! The current concentration of a species in the bulk subregion
  !
  ! For bulk simulations this is the concentration in the entire
  ! simulation volume. For simulations of a cell this is the
  ! concentration in the 'bulk' smpling sub-region.
  double precision function bulk_concentration(ispec)
    use simstate
    use geom
    implicit none
    integer, intent(in) :: ispec
    double precision :: volume
    if (is_bulk()) then
      volume = volblk()
    else
      volume = vbulk
    end if
    bulk_concentration = rs_mean(abulk(ispec))*tosi/volume
  end function bulk_concentration

  ! ------------------------------------------------------------------
  ! Write out module parameters as per input file
  subroutine ecaccu(fid)
  use strngs
  use trial
  implicit none
  integer, intent(in) :: fid
  character(8) :: lglout
  character(20) :: fltout
  write(unit=fid,fmt='(A)')fsaccu
  call str(calgin, lglout)
  write(unit=fid,fmt='(A,1X,A)')fscgin,trim(lglout)
  call str(calrdf, lglout)
  write(unit=fid,fmt='(A,1X,A)')fscrdf,trim(lglout)
  call str(calacc, lglout)
  write(unit=fid,fmt='(A,1X,A)')fsclac,trim(lglout)
  call str(calmob, lglout)
  write(unit=fid,fmt='(A,1X,A)')fsclmb,trim(lglout)
  call str(calwid, lglout)
  write(unit=fid,fmt='(A,1X,A)')fswidm,trim(lglout)
  write(unit=fid,fmt='(A,1X,I6)')fsiwid,nwdtry
  call str(drg, fltout)
  write(unit=fid,fmt='(A,1X,A)')fsdrg,trim(adjustl(fltout))
  call str(zocc, fltout)
  write(unit=fid,fmt='(A,1X,A)')fsgzoc,trim(adjustl(fltout))
  write(unit=fid,fmt='(A,1X,I6)')fsisav,isave
  write(unit=fid,fmt='(A)')fsend
  write(unit=fid,fmt=*)
  end subroutine ecaccu

  ! ------------------------------------------------------------------
  ! Return the average number of particles in a 'gz' bin.
  !
  double precision function gz_av(ibin, ispec)
  implicit none
  integer, intent(in) :: ibin, ispec
  gz_av=hist_mean(gzhist,ibin,ispec)
  end function gz_av

  ! ------------------------------------------------------------------
  ! Accumulate simple statistics
  !
  !  This subroutine gathers statistics during a simulation run.  It is
  !  intended that this method be as fast as possible as it is the most
  !  often called statistics gathering method.  The method runs through
  !  the particle list once.
  !
  !  Collection of the z-axial distribution (gz) always occurs.  It uses
  !  gz_bin to get the axial bin number.
  !
  !  The 'calmob' set of statistics are related to motion of the mobile
  !  structural ions.  These ions are new in version 17 and this method
  !  is to quantify how well the code is behaving with respect to
  !  movement around the centre-point.
  !
  !  The 'calgin' set of statistics is a 2D distribution around the
  !  z-axis and r-radial.  Unlike the gz distribution, the data is only
  !  collected here in a subregion of the simulation space.  It
  !  determines the histogram bins by dividing the z and r dimensions by
  !  drg.
  !
  !  The 'calrdf' set of statistics is the inter-particle radial
  !  distribution around particles in region 1 or 2.  The bin number is
  !  determined from the radial distance divided by 'drdf'.  (Note: The
  !  maximum number of bins is < 15/drdf.)
  !
  subroutine hist
  use conf
  use geom
  use spec
  use trial
  implicit none

  ! LOCALS
  double precision rzi,r2i,zabs    ! xz and r positions of particles
  double precision rij        ! interparticle distances
  integer ireg,ispec,jspec ! loop counters for regions,species
  integer z_bin_,r_bin_,x_bin_,y_bin_,t_bin_    ! histogram indices
  integer ii,jj          ! particle indices (global,local)
  double precision :: dx, dy, dz, rsq ! Mobile ion statistics
  double precision :: grtz_theta ! zero angle to use
  athist=athist+1
  b0n8h5: if (calgin) then
    grtz_theta = atan2(rx(grtz_zero), ry(grtz_zero))
  endif b0n8h5

  t6j9p4: do ii=1,nactv
    ispec=ispcbk(ii)
    if (ispec.eq.0) cycle t6j9p4

    rzi=rz(ii)
    zabs=dabs(rzi) 

    ! Collect z-axial distribution
    z_bin_=gz_bin(rzi)
    call hist_push(gzhist,z_bin_,ispec,1)

    ! Collect statistics about the mobile ions
    if (calmob) then
      if (localized(ispec)) then
        dx=rx(ii)-rsx(ii)
        dy=ry(ii)-rsy(ii)
        dz=rzi-rsz(ii)
        rsq=sqr(dx)+sqr(dy)+sqr(dz)
        amobdl(ii)=amobdl(ii)+sqrt(rsq)
        amobvr(ii)=amobvr(ii)+rsq
        amobdx(ii)=amobdx(ii)+dx
        amobdy(ii)=amobdy(ii)+dy
        amobdz(ii)=amobdz(ii)+dz
      endif
    endif

    ! Collect z-azial & r-radial distrbution
    b0n8h1: if (calgin) then
      r2i=r2(ii)
      u0m0i7: if (zabs.lt.-zinlft.and.r2i.lt.rinup) then
        z_bin_=max(1,ceiling((rzi-zinlft)/drg))
        r_bin_=max(1,ceiling(r2i/drg))
        call hist_push (gin, r_bin_, z_bin_, nspec()+1, int(xz(ispec)))
        call hist_push (gin, r_bin_, z_bin_, ispec, 1)
      endif u0m0i7
      if (zabs.lt.d3_zmax.and.r2i.lt.d3_rmax) then
        z_bin_=max(1,ceiling((rzi + d3_zmax)/d3_width))
        x_bin_=max(1,ceiling((rx(ii) + d3_rmax)/d3_width))
        y_bin_=max(1,ceiling((ry(ii) + d3_rmax)/d3_width))
        call hist_push (gxyz, x_bin_, y_bin_, z_bin_, ispec, 1)
        t_bin_=theta_bin(rx(ii), ry(ii), grtz_theta)
        r_bin_=max(1,ceiling(r2i / d3_width))
        call hist_push (grtz, r_bin_, t_bin_, z_bin_, ispec, 1)
      endif
    endif b0n8h1

    ! inter-particle radial distribution in region 1 and 2
    m4w8g4: if (calrdf) then
      call inregn(rzi,ispec,ireg)
      if (ireg.le.ifilt) then
        m6l4o9: do jj=1,nactv
          jspec=ispcbk(jj)
          if (jspec.eq.0) cycle m6l4o9
          
          j4d4z9: if (jj.ne.ii) then
            if (dbc) then
              if (dfeq(rqqii(jj,ii),0.D0)) then
                write(unit=fidlog,fmt=*)"ERROR"
                write(unit=fidlog,fmt=*)"rqqii(jj=",jj,",ii=",ii,")=",rqqii(jj,ii)
                write(unit=fidlog,fmt=*)"ispec=",ispec,", jspec=",jspec
                write(unit=fidlog,fmt=*)"xq(ispec)=",xq(ispec),", xq(jspec)=",xq(jspec)
                stop "Error in rqqii"
              endif
            endif
            if (.not.(dfeq(xq(ispec),0.D0).or.dfeq(xq(jspec),0.D0))) then
              rij=xq(ispec)*xq(jspec)/rqqii(jj,ii)
            else
              rij=rqqii(jj,ii)
            endif
            m7a9v3: if (rij.le.15.0D0) then
              if (dbc) then
                if (rij.lt.0.D0) stop "Rij is negative"
              endif
              r_bin_=max(1,ceiling((rij-dd_get(ispec,jspec))/drdf))
              if (dbc) then
                if (r_bin_.le.0) then
                  write(unit=fidlog,fmt=*)"Distance between particles ",ii," (",fspc(ispec),") and ",jj," (",fspc(jspec),&
                               ") to small: ",rij," < ",dd_get(ispec,jspec)
                  stop "Particle overlap detected."
                endif
              endif
              ! ireg <= ifilt
              call hist_push(rdfhist(ireg,jspec,ispec), r_bin_, 1)
              call hist_end_push(rdfhist(ireg,jspec,ispec))
              if (ireg.eq.izlim) then ! assume region 1 is inside region 2
                call hist_push(rdfhist(ifilt,jspec,ispec), r_bin_, 1)
                call hist_end_push(rdfhist(ifilt,jspec,ispec))
              endif
            endif m7a9v3
          endif j4d4z9
        enddo m6l4o9
      endif
    endif m4w8g4

  enddo t6j9p4
  call hist_end_push(gzhist)
  call hist_end_push(gin)
  call hist_end_push(gxyz)
  call hist_end_push(grtz)
  end subroutine hist


  subroutine reset_wasdel
    use trial
    use spec
    use geom
    implicit none
    ! integer :: ispec
    !do ispec=idxcl(),nspec()
    !   if (bulk_concentration(ispec).lt.ctargi(ispec)) then
    !      wasdel(ispec)=.true.
    !   else
    !      wasdel(ispec)=.false.
    !   endif
    !end do
  end subroutine reset_wasdel

  ! ----------------------------------------------------------------------
  ! Iterative method for computing chemical potentials to use in Grand
  ! Canonical part of simulation.
  !
  ! In the pre-equilibration stage, the programmer can select a method to
  ! iteratively adjust the chemical potentials.  In the protein simulation
  ! cell, the data must be collected from the bulk sub-region.
  !
  ! (acept1) Iterative method for computing chemical potentials
  ! to use in Grand Canonical part of simulation, during the pre-equilibration
  ! phase of channel simulation when particle addition/deletion are controlled
  !
  ! (malas1 is method I from)
  ! Attila Malasics, Dirk Gillespie and Dezso¨ Boda ""Simulating 
  ! prescribed particle densities in the grand canonical
  ! ensemble using iterative algorithms"", The Journal of Chemical 
  ! Physics, 2008, 128, 124102
  !
  ! (malas2 is method II from)
  ! Attila Malasics, Dirk Gillespie and Dezso¨ Boda ""Simulating 
  ! prescribed particle densities in the grand canonical
  ! ensemble using iterative algorithms"", The Journal of Chemical 
  ! Physics, 2008, 128, 124102
  subroutine iterat(cpmeth,variant)
  use conf
  use strngs
  use spec
  use geom
  use patch
  use trial
  use simstate
  implicit none
  ! simstate%is_bulk : indicate if PBC Bulk simulation or Bounded Cell simulation
  ! The update method to use
  integer, intent(in) :: cpmeth, variant

  ! Constant for non-zero charge.
  double precision :: csloth
  ! Factor for non-zero charge 
  double precision :: chcons
  ! Non-zero charge correcton
  double precision :: chisp
  ! Average charge
  double precision :: avchg
  ! Volume ('volblk' for bulk sim and 'vbulk' for channel sim)
  double precision :: bulk_conc
  double precision :: volume
  integer :: ispec, ii, jspec
  ! current particle count in "bulk" sampling subregion of channel sim
  integer, dimension(nspcmx) :: ncurr
  double precision :: zabs
  ! Lamperski delta
  double precision, save :: dchex = 0.05
  ! Lamperski: indicate if in the last step the target concentration
  ! was above the measured concentration. Then when direction changes reduce
  ! dchex by 80%
  !LMPSKI_OPT logical, save :: lmp_up = .true.
  ! malasic 2 solver data
  integer :: dimA, lda, ldb, nrhs, info ! solver data
  integer, dimension(idxcl():nspec()) :: ipiv ! solver data
  double precision, dimension(idxcl():nspec()) :: rhs, navrg
  double precision, dimension(idxcl():nspec(),idxcl():nspec()) :: fluctn
  ! scratch variables for block-local temporaries
  double precision :: tmp, tmp2, ac6n, ac6nt
  double precision :: mu_a, mu_d, p_a, p_d
  double precision :: accvolfrac
  double precision, dimension(nspcmx), save :: accvolfracold = 0 

  character(15) :: methname
  external :: dgesv
  csloth = 2*(6*dlog(2+dsqrt(3.0d0))-pi)
  ipiv=0
  avchg=0
  chisp=0
  zabs=0
  athist=athist+1
  rhs=0
  navrg=0
  fluctn=0
  chcons=0.0D0
  write(unit=fidlog,fmt='(72("-"))')
  write(unit=fidlog,fmt='(" Chem. Pot. estimation iteration ",i6)')athist
  if (is_bulk()) then
    if (do_electrostatic()) then
      avchg=rs_mean(acharge)
      write(unit=fidlog,fmt=*)"Average system charge:",avchg
      chcons=csloth*(qstar()/epsw)**2/(8*lenblk())
    endif
    ncurr=ni
    volume=volblk()
  else
    ncurr=0
    ! not in bulk, need to manually get the elements in bulk???????
    do ii=1,nactv
      if (ispcbk(ii).lt.idxcl()) cycle
      zabs=abs(rz(ii))
      if (zabs.ge.zbulk1.and.zabs.le.zbulk2) then
        if (r2(ii).le.rbulk) ncurr(ispcbk(ii))=ncurr(ispcbk(ii)) + 1
      endif
    enddo
    if (do_electrostatic()) then
      chcons=csloth*(qstar()/epsw)**2/(8*(vbulk**(1.D0/3.D0)))
    endif
    volume=vbulk
  endif

  if (cpmeth.eq.malas2) then
    do ispec=idxcl(),nspec()
      navrg(ispec) = rs_mean(abulk(ispec))
      rhs(ispec)=ctargi(ispec)-(navrg(ispec)*tosi/volume)
    enddo
    do ispec=idxcl(),nspec()
      do jspec=idxcl(),nspec()
        fluctn(ispec,jspec) = ((ninj(ispec,jspec)/dble(ataccu))-(navrg(ispec)*navrg(jspec)))*beta()/volume
      enddo
    enddo
    if (debug) then
      write(unit=fidlog,fmt=*)"(rhs matrix)"
      write(unit=fidlog,fmt=*)(rhs(ispec),ispec=idxcl(),nspec())
      write(unit=fidlog,fmt='(72("-"))')
      write(unit=fidlog,fmt=*)"(Fluctn matrix)"
      do ispec=idxcl(),nspec()
        write(unit=fidlog,fmt=*)(fluctn(ispec,jspec),jspec=idxcl(),nspec())
      enddo
      write(unit=fidlog,fmt='(72("-"))')
    endif
    dimA=nspec()-idxcl()+1
    nrhs=1
    lda=dimA
    ldb=dimA
    info=0
    if (debug) then
      write(unit=fidlog,fmt=*)' dimA=',dimA,' nrhs=',nrhs,' lda=',lda,' ldb=',ldb,' info=',info
    endif
    call dgesv(dimA, nrhs, fluctn, lda, ipiv, rhs, ldb, info)
    if (debug) then
      write(unit=fidlog,fmt=*)' after call:'
      write(unit=fidlog,fmt=*)' dimA=',dimA,' nrhs=',nrhs,' lda=',lda,' ldb=',ldb,' info=',info
      write(unit=fidlog,fmt=*)"(rhs matrix)"
      write(unit=fidlog,fmt=*)(rhs(ispec),ispec=idxcl(),nspec())
    endif
  endif
  do ispec=idxcl(),nspec()
    write(unit=fidlog,fmt='(72("-"))')
    ! mu_i(n) = kT.log(rho^targ_i)+mu_ex(n)
    !
    ! mu^ex_i(n) = mubar^ex_i(n-1)
    ! mubar^ex_i(n-1)=mu_i(n) - kT log(rhobar_i(n))
    ! 
    ! chempi(n)*beta()={dlog(ctargi(ispec)) - dlog(tosi) + chexi(ispec)*beta()
    !                           - dlog(bulk_conc(ispec)) + dlog(tosi)}
    ! chempi(n)       ={dlog(ctargi(ispec)/bulk_conc(ispec))/beta() + chempi(ispec)}

    ! non-zero charge correction
    chisp=-xz(ispec)*avchg*chcons
    if (dfeq(rs_mean(abulk(ispec)),0.D0)) then
       write(unit=fidlog,fmt=*)"Specie ",fspc(ispec)," has no particles, ignoring."
       cycle
    endif

    bulk_conc = bulk_concentration(ispec)
    call accget (ispec,mu_a,mu_d,p_a,p_d,accvolfrac)

    select case (cpmeth)
    case (0) ! no update
    case (lamperski)
      methname = "Lamperski"
      if (bulk_conc.gt.ctargi(ispec)) then
        !LMPSKI_OPT if (lmp_up) dchex = dchex * 0.8
        !LMPSKI_OPT lmp_up=.false.
        call chexi_set (ispec, chexi(ispec) - dchex)
      else
        !LMPSKI_OPT if (.not.lmp_up) dchex = dchex * 0.8
        !LMPSKI_OPT lmp_up=.true.
        call chexi_set (ispec, chexi(ispec) + dchex)
      endif

    case (malas1)
      methname = "Malasics 1"
      ! change in chemical potential is ratio of target and current densities
      ! use average density when close to convergence
      ! concentration in current system
      call chexi_set (ispec, chempi(ispec) - log(bulk_conc/tosi) + chisp)
    case (malas2)
      methname = "Malasics 2"
      ! Only update if info was zero
      if (info.eq.0) then
        call chexi_set (ispec, chexi(ispec) + rhs(ispec) + chisp)
      else
        write(unit=fidlog,fmt=*) "** warning: failed matrix inversion **"
      endif
    case (accept)
      write(methname,'("Accept",I4)')variant
      !if (rs_mean(abulk(ispec)).lt.rs_count(abulk(ispec))) then
      !  write(unit=fidlog,fmt=*)"Specie ",fspc(ispec)," has too few particles to give reliable results, ignoring."
      !  cycle
      !endif
      ! Approximate density ratio with trial ratios
      ! log(try(+)/try(-))/2
      navrg(ispec)=rs_mean(abulk(ispec))
      select case (variant)
      case (13)
        write(unit=fidlog,fmt=*)"#########       Variant 13      ##########"
        ac6nt=ctargi(ispec)*volblk()/tosi
        ac6n=floor(ac6nt)
        tmp=(1.D0-erf((0.5D0-ac6nt)/sqrt(2.D0*ac6nt)))/2.D0
        write(*,*)"Multiply factor : N_t=",ac6nt," N=",ac6n," factor=",tmp
        call chexi_set (ispec, chexi(ispec) + dlog(tmp*p_d/p_a)/2.D0 + chisp)
        write(unit=fidlog,fmt=*)"##### ln(1-erf(1/2-N_t)/sqrt(2*N_t)) #####"
      case (12)
        write(unit=fidlog,fmt=*)"#########       Variant 12      ##########"
        ac6nt=ctargi(ispec)*volblk()/tosi
        ac6n=floor(ac6nt)
        if (ac6nt.gt.1.D0) then
          tmp=(p_d)
        else
          tmp=p_d*ac6nt
        endif
        write(*,*)"Multiply factor : N_t=",ac6nt," N=",ac6n," factor=",tmp
        call chexi_set (ispec, chexi(ispec) + dlog(tmp/p_a)/2.D0 + chisp)
        write(unit=fidlog,fmt=*)"##### ln(N_t<1?P_d*N_t:P_d) #####"
      case (11)
        write(unit=fidlog,fmt=*)"#########       Variant 11      ##########"
        ac6nt=ctargi(ispec)*volblk()/tosi
        ac6n=floor(ac6nt)
        tmp=(ac6nt+0.5D0)/(ac6n+0.5D0)
        write(*,*)"Multiply factor : N_t=",ac6nt," N=",ac6n," factor=",tmp
        call chexi_set (ispec, chexi(ispec) + dlog((p_d/p_a) / tmp)/2.D0 + chisp)
        write(unit=fidlog,fmt=*)"##### ln((P_d/P_a)/(N_t+0.5D0/N+0.5D0))/2 #####"
      case (10)
        write(unit=fidlog,fmt=*)"#########       Variant 10      ##########"
        ac6nt=ctargi(ispec)*volblk()/tosi
        ac6n=floor(ac6nt)
        if (ac6nt.gt.1.D0) then
          tmp=(ac6nt*ac6nt)/(ac6n*ac6n+ac6n)
        else
          tmp=6.D0*ac6nt
        endif
        write(*,*)"Multiply factor : N_t=",ac6nt," N=",ac6n," factor=",tmp
        call chexi_set (ispec, chexi(ispec) + dlog((p_d/p_a) / tmp)/2.D0 + chisp)
        write(unit=fidlog,fmt=*)"##### ln((P_d/P_a)/(N_t**2/N(N+1)))/2 #####"
      case (9)
      write(unit=fidlog,fmt=*)"#########       Variant 9      ##########"
      if (p_a.gt.p_d) then
        !LMPSKI_OPT ! if (lmp_up) dchex = dchex * 0.8
        !LMPSKI_OPT lmp_up=.false.
        call chexi_set (ispec, chexi(ispec) - dchex + chisp)
      else
        !LMPSKI_OPT ! if (.not.lmp_up) dchex = dchex * 0.8
        !LMPSKI_OPT lmp_up=.true.
        call chexi_set (ispec, chexi(ispec) + dchex + chisp)
      endif
      write(unit=fidlog,fmt=*)"#########  (P_d>P_a) ? +D : -D  ##########"
      case (8)
        ! Use adjustment for accessible volume fraction
        write(unit=fidlog,fmt=*)"#########       Variant 8      ##########"
        call chexi_set (ispec, chexi(ispec) + dlog(p_d/p_a*accvolfrac)/2.0D0 + chisp)

        !  - mu_a - mu_d
        write(unit=fidlog,fmt=*)"#########  ln(P_d/P_a*AccVol)/2  ##########"
        !  if (dchex.lt.1.0) then
        !    dchex=dchex+0.05
        !  endif
       case (7)
        ! update is ln(PA/PR) / (-2 - kTe(-mu)/PA)
        ! tmp=(2 + (exp(-chexi(ispec))/tmp))
        write(unit=fidlog,fmt=*)"#########       Variant 7      ##########"
        if (accvolfracold(ispec).eq.0.0D0) then
           accvolfracold(ispec)=accvolfrac
        endif
        call chexi_set (ispec, chexi(ispec) + dlog(p_d/p_a)/2.0D0 &
                + dlog(accvolfracold(ispec)/accvolfrac) + chisp)
        accvolfracold(ispec)=accvolfrac
        !  - mu_a - mu_d
        write(unit=fidlog,fmt=*)"######## ln(P_d*AV(i-1)/P_a*AV(i)) #######"
        !  if (dchex.lt.1.0) then
        !    dchex=dchex+0.05
        !  endif
       case (6)
        ! update is ln({Pd/Pa}*{1/(N + 1 - N_t)-1}) / 2
        ! tmp=(2 + (exp(-chexi(ispec))/tmp))
        write(unit=fidlog,fmt=*)"#########       Variant 6      ##########"
        ! number density as real
        ac6nt=ctargi(ispec)*volblk()/tosi
        ac6n=floor(ac6nt)
        tmp=p_a/(p_a + p_d)
        tmp2=p_d/(p_a+p_d)
        call chexi_set (ispec, chexi(ispec) + dlog(ac6nt/(tmp*(ac6n+1)+tmp2*ac6n))/2.D0 + chisp)
        write(unit=fidlog,fmt=*)"## ln(Nt/{Pa(N+1)/(Pa+Pd)+PdN/(Pa+Pd)})/2 ##"
        !  if (dchex.lt.1.0) then
        !    dchex=dchex+0.05
        !  endif
      case (5)
        write(unit=fidlog,fmt=*)"#########       Variant 5      ##########"
        ac6nt=ctargi(ispec)*volblk()/tosi
        ac6n=floor(ac6nt)
        if (ac6nt.gt.1.D0) then
          tmp=(ac6nt*ac6nt)/(ac6n*ac6n+ac6n)
        else
          tmp=6.D0*ac6nt
        endif
        write(*,*)"Multiply factor : N_t=",ac6nt," N=",ac6n," factor=",tmp
        call chexi_set (ispec, chexi(ispec) + dlog((p_d/p_a) * tmp)/2.D0 + chisp)
        write(unit=fidlog,fmt=*)"##### ln((P_d/P_a)*(N_t**2/N(N+1)))/2 #####"
      case (4)
        write(unit=fidlog,fmt=*)"#########       Variant 4      ##########"
        ac6nt=ctargi(ispec)*volblk()/tosi
        ac6n=floor(ac6nt)
        if (ac6nt.gt.1.D0) then
          tmp=(ac6nt*ac6nt)/(ac6n*ac6n+ac6n)
        else
          tmp=6.D0*ac6nt
        endif
        write(*,*)"Multiply factor : N_t=",ac6nt," N=",ac6n," factor=",tmp
        call chexi_set (ispec, chexi(ispec) + dlog((p_d/(p_a*accvolfrac)) * tmp)/2.D0 + chisp)
        write(unit=fidlog,fmt=*)"##### ln((P_d/P_a*p)*(N_t**2/N(N+1)))/2 #####"
      case (3)
        write(unit=fidlog,fmt=*)"#########       Variant 3      ##########"
        ac6nt=ctargi(ispec)*volblk()/tosi
        ac6n=floor(ac6nt)
        if (ac6nt.gt.1.D0) then
          tmp=(ac6nt*ac6nt)/(ac6n*ac6n+ac6n)
        else
          tmp=6.D0*ac6nt
        endif
        write(*,*)"Multiply factor : N_t=",ac6nt," N=",ac6n," factor=",tmp
        call chexi_set (ispec, chexi(ispec) + dlog((p_d*accvolfrac/(p_a)) * tmp)/2.D0 + chisp)
        write(unit=fidlog,fmt=*)"##### ln((P_d*p/P_a)*(N_t**2/N(N+1)))/2 #####"
      case (2)
        write(unit=fidlog,fmt=*)"#########       Variant 2      ##########"
        call chexi_set (ispec, chexi(ispec) + dlog(p_d*accvolfrac/p_a) /2.D0 + chisp)
        write(unit=fidlog,fmt=*)"##### ln(P_d*p/P_a)/2 #####"
      case default
        write(unit=fidlog,fmt=*)"##########      Variant 1      ##########"
        call chexi_set (ispec, chexi(ispec) + (dlog(p_d/p_a))/2.D0 + chisp)
        write(unit=fidlog,fmt=*)"##########    ln(P_d/P_a)/2    ##########"
      end select
    case default
      stop "Unknown chemical potential update method"
    end select

    write(unit=fidlog,fmt=*)"SPECIE: ",fspc(ispec)
    write(unit=fidlog,fmt='(1X,"SPC","|",5X,"METHOD",1X,"|",3X,"<[BULK]>",1X,& 
        &"|",5X,"[BULK]",1X,"| TOTAL|",3X,"CHEM POT",1X,"|EXCESS C.P.",1X)') 
    bulk_conc = bulk_concentration(ispec)
    write(unit=fidlog,fmt='(1X,A2,1X,"|",1X,A10,1X,"|",1X,F10.5,1X,"|",1X,F10.5,1X,& 
        &"|",1X,I5,"|",1X,F10.5,1X,"|",1X,F10.5,1X)')& 
        &fspc(ispec),methname,bulk_conc,tosi*ni(ispec)/vtotal(ispec),& 
        &ni(ispec),chempi(ispec),chexi(ispec)
    write(unit=fidlog,fmt='(1X,A25,4(1X,F10.6))')"MU_a, MU_d, P_a, P_d :",mu_a,mu_d,p_a,p_d
    call avergi(ispec)
  enddo
  write(unit=fidlog,fmt='(72("-"))')
  call avergeng
  
  write(unit=fidlog,fmt='(72("-"))')
  call accrat(fidlog)
  ! report chemcial potentials
  call chemical_potential_report(fidlog)

  do ispec=1,nspec()
    call rs_reset(abulk(ispec))
  end do
  if (do_electrostatic()) call rs_reset(acharge)
  ataccu=0
  call zeroav

  end subroutine iterat

  ! -----------------------------------------------------
  ! Read statistical data control information section
  !
  ! This routine is called to read-in the input data for the 'accum'
  ! module.  See 'channel%readin' method for more details of how
  ! system is initialised.
  !
  ! @param fid : input unit number
  ! @param sname : the name value that caused this function to be called
  ! @param svalue : the value associated with the name (may be empty string)
  !
  ! @pre sname=fsaccu and len(svalue)=0
  subroutine rdaccu(fid,sname,svalue,istat)
  use strngs
  use trial
  implicit none
  integer, intent(in) :: fid
  character(len=*), intent(in) :: sname,svalue
  integer, intent(out) :: istat
  logical, dimension(2) :: mask_
  character(32) :: nme_
  character(1024) :: val_

  if (dbc) then
    if (sname.ne.fsaccu) stop "Error: incorrect section name"
    if (len_trim(svalue).ne.0) stop "Error: section does not any parameters"
  endif
  mask_=.false.
  !
  ! INPUT SECTION
  !
  ! accum
  ! drg REAL
  ! calgin (.true.|.false.)
  ! calacc (.true.|.false.)
  ! calrdf (.true.|.false.)
  ! calmob (.true.|.false.)
  ! calwid (.true.|.false.)
  ! iwidom INT
  ! isave INT
  ! end
  !
  t7k4b6: do
    val_ = " "
    call readnv(fid, nme_, val_, istat)
    ! return on bad read
    if (istat.ne.0) return
    ! exit loop on section 'end'
    if (nme_.eq.fsend) exit t7k4b6
    ! looking for calgin, calrdf, etc
    u2h7j1: select case (nme_)
      case (fscgin) u2h7j1
        read(val_,*)calgin
      case (fscrdf) u2h7j1
        read(val_,*)calrdf
      case (fsclac) u2h7j1
        read(val_,*)calacc
      case (fsclmb) u2h7j1
        read(val_,*)calmob
      case (fswidm) u2h7j1
        read(val_,*)calwid
      case (fsiwid) u2h7j1
        read(val_,*)nwdtry
      case (fsgzoc) u2h7j1
        read(val_,*)zocc
      case (fsdrg) u2h7j1
        read(val_,*)drg
        mask_(1)=.true.
      case (fsisav) u2h7j1
        read(val_,*)isave
        mask_(2)=.true.
      case default u2h7j1
        call v7t9f4("Name "//nme_//" is not valid in statistic (accum) section")
    end select u2h7j1
  enddo t7k4b6
  d0q8g6: if (.not.all(mask_)) then
    call v7t9f4("Not all required tags were present.")
  endif d0q8g6
 
  if (calgin) then
     if (dfeq(drg,0.D0)) call v7t9f4("When collecting 3D statistics "//fsdrg//" option may not be zero.")
  endif

  contains

  ! print brief error message on incorrect input section
  subroutine v7t9f4(msg)
  use strngs
  implicit none
  character(len=*), intent(in) :: msg

    write(unit=fidlog,fmt=*)"Bad statistic (accum) section in input:"
    write(unit=fidlog,fmt=*)msg
    write(unit=fidlog,fmt=*)"Required tags are:"
    write(unit=fidlog,fmt=*)fsdrg," and ",fsisav
    write(unit=fidlog,fmt=*)"These optional tags control what results are obtained (default to false):"
    write(unit=fidlog,fmt=*)fscgin,", ",fscrdf,", ",fsclac,", ",fsclmb," and ",fswidm
    write(unit=fidlog,fmt=*)"This optional tag is the width of inner occupancy region (defaults to zlimit):"
    write(unit=fidlog,fmt=*)fsgzoc
    write(unit=fidlog,fmt=*)"If ",fswidm," is .true. then ",fsiwid," indicates the minimum number of"
    write(unit=fidlog,fmt=*)"test particle insertions to be used."
   
    stop 1
  end subroutine v7t9f4

  end subroutine rdaccu

  ! ---------------------------------------------------------------
  ! Initialise data arrays and calculate dimensions/regions for 
  ! statistics collection
  !
  ! Key functions:
  !  - allocate the arrays that are needed
  !  - on 'calgin' initialise the 'gin' distribution and geometry parameters.
  !  - initialise the 'gzhist' distribution (call geometry def. in trial.mod)
  !  - on 'calrdf' initialise the 'rdf' distribution
  !  - echo interpreted content of input file
  subroutine rfaccu
  use conf
  use geom
  use spec
  use trial
  implicit none

  write(unit=fidlog,fmt='(72("-"))')
  write(unit=fidlog,fmt=*)"Statistical data accumulation parameters"
  write(unit=fidlog,fmt='(72("-"))')

  if (dfeq(zocc,0.D0)) then
    zocc=zlimit()
  else
    if (zocc.gt.zl(1)) then
      write(unit=fidlog,fmt=*)"Limit of occupancy region 1 is outside central cylinder (zocc > zl1)"
      stop 1
    endif
  endif

  ! Geometry for rdf
  if (calrdf) then
    drdf=drg
    nrdf=int(15.D0/drdf)+1
  endif

  write(unit=fidlog,fmt='(72("-"))')
  write(unit=fidlog,fmt=*)"Interpreted statistic/accumulation parameters from input"
  write(unit=fidlog,fmt='(72("-"))')
  call ecaccu(fidlog)
  write(unit=fidlog,fmt='(72("-"))')

  ! initialise gin 2D distribution
  if (calgin) call v8x4m0

  ! allocate arrays now we know there minimum sizes
  call t7n1g5

  contains

  ! --------------------------------------------------
  ! Ensure allocation of arrays

  subroutine t7n1g5
  use spec
  use geom
  use simstate
  implicit none
  ! Locals
  double precision :: factor    ! temporaries
  integer          :: i,j        ! index
  integer :: zsize_, rsize_, nsize_ ! computed array sizes
  if (.not.allocated(ah)) then
     zsize_ = next64(gz_max)
     nsize_ = next2(nspec())
     call hist_init(gzhist,zsize_,nsize_)
     call hist_init(anin,nrgnmx,nsize_)
     if (do_electrostatic()) then
       call rs_init(acharge)
       allocate(ah(npchmx))
       ah=0
     endif
  endif
  if (calrdf) then
    rsize_ = next64(nrdf)
    nsize_ = next2(nspec())
    if (.not.allocated(rdfhist)) allocate(rdfhist(izlim:ifilt,nsize_,nsize_))
    do i=1,nspec()
      do j=1,nspec()
        call hist_init(rdfhist(izlim,j,i), rsize_)
        call hist_init(rdfhist(ifilt,j,i), rsize_)
      enddo
    enddo
  endif
  if (calgin) then
    zsize_ = next64(nrgz)
    rsize_ = next64(nrgr)
    if (.not.allocated(ginvol)) allocate(ginvol(rsize_))
    nsize_ = next2(nspec()+1)
    call hist_init(gin, rsize_,zsize_,nsize_)

    zsize_ = next64(d3_zgrid_size)
    rsize_ = next64(d3_rgrid_size)
    nsize_ = next2(nspec())
    if (.not.allocated(d3_vol)) allocate(d3_vol(rsize_))
    call hist_init(gxyz,rsize_,rsize_,zsize_,nsize_)
    call hist_init(grtz,rsize_/2,rsize_,zsize_,nsize_)

    ! --------------------------------------------------
    ! Compute volumes of slabs
    factor=drg**3*pi
    b8g5l1: do j=1,nrgr
      ! -- formula --
      ! ginvol = drg * 3 * pi *((j*drg)**2 - ((j-1)*drg)**2)
      ginvol(j) = factor*(2*j-1)
    enddo b8g5l1
    factor=d3_width**3*pi/dble(d3_rgrid_size)
    l1b8g5: do j=1,d3_rgrid_size
      ! -- formula --
      ! d3_vol = d3_width * 3 * pi *((j*d3_width)**2 - ((j-1)*d3_width)**2) 
      !          -----------------------------------------------------------
      !                             d3_rgrid_size
      d3_vol(j) = factor*(2*j-1)
    enddo l1b8g5
  endif
  end subroutine t7n1g5

  ! --------------------------------------------------
  ! Geometry for gin 2D distributions
  subroutine v8x4m0
  use geom
  use trial
  use strngs
  implicit none

  if (dbc) then
    if (.not.calgin) stop "ERROR: accum%v8xm0 called when calgin was false"
  endif
  grtz_zero=1
  rinup  = rl(4)+10
  zinlft = zl(2)+10
  nrgr   = max(1,nint(rinup/drg+1))
  rinup  = nrgr*drg
  nrgz   = max(1,nint(zinlft/drg+1))
  zinlft = -nrgz*drg
  nrgz   = 2*nrgz
  write(unit=fidlog,fmt=*)"Parameters for 2D histogram of Z-axis and Radius"
  write(unit=fidlog,fmt='(" Histogram bins in radial dimension = ",I8," width = ",F8.4)')nrgr,drg
  write(unit=fidlog,fmt='(" Histogram bins in z-axial dimension= ",I8," width = ",F8.4)')nrgz,drg

  d3_zmax      = zl(1) + min(rl(1),rlvest())
  d3_rmax      = rl(1) + min(rl(1),rlvest())
  ! want total size of 3D hist to equal the same as 2D hist
  ! Increment d3_width until we get a reasonable size
  d3_width     = drg
  d3_zgrid_size= 2*max(1,nint(d3_zmax/d3_width))
  d3_rgrid_size= 2*max(1,nint(d3_rmax/d3_width))
  do while(d3_zgrid_size*d3_rgrid_size*d3_rgrid_size.gt.2*nrgr*nrgz)
    d3_width     = d3_width+drg
    d3_zgrid_size= 2*max(1,nint(d3_zmax/d3_width))
    d3_rgrid_size= 2*max(1,nint(d3_rmax/d3_width))
  enddo
  d3_zmax      = d3_width * d3_zgrid_size/2
  d3_rmax      = d3_width * d3_rgrid_size/2

  write(unit=fidlog,fmt=*)"Parameters for 3D histograms of X,Y,Z and R,theta,Z axes"
  write(unit=fidlog,fmt='(" Histogram bins in x,y,theta,r dimensions = ",I8," width = ",F8.4)')d3_rgrid_size,d3_width
  write(unit=fidlog,fmt='(" Histogram bins in z-axial dimension      = ",I8," width = ",F8.4)')d3_zgrid_size,d3_width
  write(unit=fidlog,fmt='(72("-"))')
  v8i1x6: if (nrgr.gt.nrgz) then
    nrg=nrgr
  else v8i1x6
    nrg=nrgz
  endif v8i1x6
  v8b7t5: if (nrg.gt.nrgmx) then
    write(unit=fidlog,fmt=*) "The number of requested bins exceed program limit"
    write(unit=fidlog,fmt='(" request = ",i6,"   limit = ",i6)')nrg,nrgmx
    write(unit=fidlog,fmt=*) "SUGGESTED SOLUTION: increase '"//fsdrg//"' option in input file"
    stop 'Program array size limit exceeded'
  endif v8b7t5

  end subroutine v8x4m0

  end subroutine rfaccu

  ! ------------------------------------------------------------------
  ! Save digest of statistics
  !
  ! This is the most important method of the accumulation module
  ! as it is responsible for creating more in-depth statistics
  ! from the simulation.  Some of the data analysis is delegated
  ! to other modules (particularly when the data can be kept private
  ! to the other module.)
  !
  !
  !
  subroutine saves(astep,loopsz)
  use const
  use conf
  use geom
  use patch
  use spec
  use trial
  use simstate
!$ use omp_lib
  implicit none
  integer, intent(in) :: astep,loopsz

  ! LOCALS
  double precision z_mid,r_mid, vol_ ! histogram dimensions
  double precision x_mid,y_mid,t_mid,xyz_vol ! histogram dimensions
  double precision anz     ! average particle count
  double precision ddi ! inter-particle distances
  double precision rlow,rhigh,ri ! histogram bin data
  double precision vshell        ! histogram bin data
  double precision rdfi1,rdfi2   ! radial dist. counts
  double precision rdf1vs,rdf2vs ! dist. density
  double precision azcl          ! occupancy sums
  double precision :: zocspc  ! per-specie occupancy limit
  integer igc,ireg,ispec,ipch ! salt,region, specie, patch indices
  integer g_,h_,k_,l_        ! histogram bin indices
  integer kend,lend      ! histogram bin ends
  integer strgn1,strgn2,strgn3,enrgn1,enrgn2,enrgn3 ! histogram bin coord for conductance
  integer i,j,jspec      ! particle indices and second specie
  integer :: fid_     ! loop unique File Ids
  double precision, dimension(nzgmx,nspcmx) :: cz
  double precision, dimension(nrgnmx+2,nspcmx), save :: anaca=0
  double precision, dimension(3,nspcmx) :: cnduct
  logical, dimension(3) :: usecd
  double precision :: cin

  cz=0
  cin=0
  cnduct=0
  kend=0
  lend=0

  if (dbc) then
    if (ataccu.eq.0) stop 'Attempt to save results before accumulating any data'
  endif

  write(unit=fidlog,fmt='(3/)')
  write(unit=fidlog,fmt='(72("-"))')
  write(unit=fidlog,fmt=*)" BEGIN SAVES FOR STEP ",astep
  write(unit=fidlog,fmt='(72("-"))')
  open(unit=fidooo,file="res/o."//firun)
  write(unit=fidooo,fmt='(" UUID    = ",A32)')fuuid
  write(unit=fidooo,fmt='(" STEP    = ",I10)')astep
  write(unit=fidooo,fmt='(" INNER   = ",I10)')loopsz
  write(unit=fidooo,fmt='(" RUN     = ",A10)')firun
  write(unit=fidooo,fmt='(" TEMP    = ",F10.3)')tmptur()
  write(unit=fidooo,fmt=*)
  write(unit=fidooo,fmt='(" zl4 [A]     = ",f10.5)')zl(4)
  write(unit=fidooo,fmt='(" rl5 [A]     = ",f10.5)')rl(5)
  write(unit=fidooo,fmt='(" rl1 [A]     = ",f10.5)')rl(1)
  if (do_electrostatic()) then
    call patchsaves(fidooo,astep)
  else
    write(unit=fidooo,fmt='(A)')"** No electrostatic potential energy included. **"
    write(unit=fidlog,fmt='(A)')"** No electrostatic potential energy included. **"
  endif

  write(unit=fidooo,fmt=*)
  q0o9y2: do igc=1,nsalt()
    write(unit=fidooo,fmt='(a4,": ctarg = ",f11.8,", chemp = ",f10.5)')&
     &       fsalt(igc),ctargs(igc),chemps(igc)
  enddo q0o9y2
  q0o9y3: do ispec=1,nspec()
    write(unit=fidooo,fmt='(a4,": ctarg = ",f11.8,", chemp = ",f10.5,", excess = ",f10.6)')&
     &       fspc(ispec),ctargi(ispec),chempi(ispec),chexi(ispec)
  enddo q0o9y3
   write(unit=fidooo,fmt=*)
  write(unit=fidooo,fmt='(12x,"   N   ","  q [e]  ","  d [A]  ")')
  a0g9f9: do ispec=1,nspec()
  write(unit=fidooo,fmt='(a2," ion: ",1x,i7,2x,f5.2,2x,2x,f8.4)')&
     &      fspc(ispec),ni(ispec),xz(ispec),xri(ispec)*2
  enddo a0g9f9
  
  if (do_electrostatic()) then
    write(unit=fidlog,fmt=*)"Average system charge:",rs_mean(acharge)
  endif
  ksub=ksub+1
  ! ----- Calculate and write out concentration profiles -------

  ! ----- Write out occupancies and stuff ------------------
  call u9c6u7
  call concentration_report(fidlog)
  call concentration_report(fidooo)

  if (calwid) then
    call trywd(astep)
    call savewd
  endif
  write(unit=fidlog,fmt='(72("-"))')
 
  ! ----- Acceptance ratios --------------------------------
  write(unit=fidlog,fmt=*)"STATISTICS OF METROPOLIS ACCEPTANCE FACTOR"
  write(unit=fidlog,fmt='(1X,A6,3(" | ",A18))')"Specie","Move (aver:var)","Create (aver:var)","Destroy (aver:var)"
  write(unit=fidlog,fmt='(72("-"))')
  do ispec=idxcl(),nspec()
    call avergi(ispec)
  enddo
  call avergeng
  write(unit=fidlog,fmt='(72("-"))')

  ! Reset accumulators used in avergi and avergeng
  call zeroav

  ! Report trial success data
  l9m3d8: if (calacc) then
    call accrat(fidlog)
    call accrat(fidooo)
  endif l9m3d8

  close(unit=fidooo)

  ! ----- Occupancies --------------------------------------
  open(unit=fidocc,file="res/occ."//firun)
  write(unit=fidocc,fmt='("# UUID ",A32)')fuuid
  write(unit=fidocc,fmt='(A)')'# label  g  h  k  l  count'
  write(unit=fidocc,fmt='(A)')'# units  -  -  -  -  count'

    if (nfree().ge.3) kend=4 ! need k loop for 3 species
    if (nfree().ge.4) lend=4 ! need l loop for 4 species
    z5p5s2: do ireg=izlim,ichan
      m9l5i8: do g_=0,4
      e4y8n5: do h_=0,4
      h5j7b7: do k_=0,kend
      w4d5u5: do l_=0,lend
         write(unit=fidocc,fmt='(4(i1,1x),f15.10)')g_,h_,k_,l_,aocc(l_,k_,h_,g_,ireg)/dble(ataccu)
      enddo w4d5u5
      enddo h5j7b7
      enddo e4y8n5
      enddo m9l5i8
      write(unit=fidocc,fmt=*)
    enddo z5p5s2
    write(unit=fidocc,fmt='(72("#"))')
    write(unit=fidocc,fmt='(A)')'# label  h  k  l  runsum'
    write(unit=fidocc,fmt='(A)')'# units  -  -  -  count'

    y9s3i3: do ireg=izlim,ichan
      e7g0x5: do h_=0,4
      g2f4m1: do k_=0,kend
      v1m0m2: do l_=0,lend
        azcl=0
        d0h0t7: do g_=0,4
           azcl=azcl+aocc(l_,k_,h_,g_,ireg)/dble(ataccu)
        enddo d0h0t7
        write(unit=fidocc,fmt='(3(i1,1x),f15.10)')h_,k_,l_,azcl
      enddo v1m0m2
      enddo g2f4m1
      enddo e7g0x5
      write(unit=fidocc,fmt=*)
    enddo y9s3i3
  close(unit=fidocc)

  ! 'h' vector
  if (do_electrostatic().and..not.is_homogeneous()) then
    open(unit=fidhpc,file="res/h."//firun)
    write(unit=fidhpc,fmt='("# UUID ",A32)')fuuid
    write(unit=fidhpc,fmt='(A)')'# label zcoord   h'
    write(unit=fidhpc,fmt='(A)')'# units ang      ENG'

    y8c6i5: do ipch=1,npatch
      write(unit=fidhpc,fmt=*)ipch,ah(ipch)/dble(ataccu)
    enddo y8c6i5
    close(unit=fidhpc)
  endif
  ! write out details of volumes occupied by the ions
  if (idxcl().ne.1) call volslb(fidlog,10,astep)

  ! ------------------------
  ! Check for histogram data
  ! ------------------------
  if (athist.eq.0) return

  ! write out statistics about the movement of 'mobile' structural ions
  if (have_localized()) call amobsv 

  k2l1f2: if (calgin) then
    ! 2D distribution 
    !$omp parallel do private(fid_,ispec,k_,l_,z_mid,r_mid,cin)
    r6o1e6: do ispec=1,nspec()
      if (debug) then
        !$omp critical
        write(unit=fidlog,fmt=*)"! gin profile for ",fspc(ispec)
        !$omp end critical
      end if
      fid_ = fidgin
!$    fid_ = fid_ + omp_get_thread_num()
      open(unit=fid_,file="res/grz-"//fspc(ispec)//"."//firun)
      write(unit=fid_,fmt='("# UUID ",A32)')fuuid
      write(unit=fid_,fmt='(A)')'# label zcoord  rcoord  conc  vol'
      write(unit=fid_,fmt='(A)')'# unit  ang     ang     molar ang3'
      c4v9i5: do k_=1,nrgz
        z_mid=zinlft+(k_-1)*drg+drg/2
        e6r0w8: do l_=1,nrgr
          vol_ = ginvol(l_)
          cin=tosi*hist_mean(gin, l_,k_,ispec)/vol_
          r_mid=(l_-1)*drg+drg/2
          write(unit=fid_,fmt=*)z_mid,r_mid,cin,vol_
        enddo e6r0w8
        write(unit=fid_,fmt=*)
      enddo c4v9i5
      close(unit=fid_)
    enddo r6o1e6
    ! $omp end parallel do

    ! 2D distribution of charge
    if (debug) write(unit=fidlog,fmt=*)"! gin profile for charge"
    open(unit=fidgin,file="res/grz-charge."//firun)
    write(unit=fidgin,fmt='("# UUID ",A32)')fuuid
    write(unit=fidgin,fmt='(A)')'# label zcoord  rcoord  charge/dens  vol'
    write(unit=fidgin,fmt='(A)')'# unit  ang     ang     molar ang3'
    c4v9i6: do k_=1,nrgz
      z_mid=zinlft+(k_-1)*drg+drg/2
      e6r0w6: do l_=1,nrgr
        cin=tosi*hist_mean(gin, l_,k_,nspec()+1)/ginvol(l_)
        r_mid=(l_-1)*drg+drg/2
        write(unit=fidgin,fmt='(5 (4x,F16.10))')z_mid,r_mid,cin,ginvol(l_)
      enddo e6r0w6
      write(unit=fidgin,fmt=*)
    enddo c4v9i6
    close(unit=fidgin)

    ! 3D distrbutions
    xyz_vol=drg**3
    r6o1e4: do ispec=1,nspec()
      if (debug) write(unit=fidlog,fmt=*)"! gxyz profile for ",fspc(ispec)
      open(unit=fidgin,file="res/gxyz-"//fspc(ispec)//"."//firun)
      write(unit=fidgin,fmt='("# UUID ",A32)')fuuid
      write(unit=fidgin,fmt='(A)')'# label xcoord  ycoord  zcoord conc  vol'
      write(unit=fidgin,fmt='(A)')'# unit  ang     ang     ang    molar ang3'
      c4v9i4: do g_=1,d3_zgrid_size ! zcoord
        !! z-range :: -d3_zmax +  d3_width/2 : (d3_zgrid_size-1)*d3_width + d3_width/2 - d3_zmax
        z_mid=(g_-1)*d3_width + d3_width/2 - d3_zmax
        e6r0w4: do k_=1,d3_rgrid_size ! ycoord
          !! y-range :: d3_width/2 - d3_rmax : (d3_rgrid_size/2 - 1)*d3_width + d3_width/2 - d3_rmax
          y_mid=(k_-1)*d3_width + d3_width/2 - d3_rmax
          e6r0w3: do l_=1,d3_rgrid_size  ! xcoord
            !! x-range :: d3_width/2 - d3_rmax : (d3_rgrid_size/2 - 1)*d3_width + d3_width/2 - d3_rmax
            x_mid=(l_-1)*d3_width + d3_width/2 - d3_rmax
            cin=tosi*hist_mean(gxyz, l_, k_, g_,ispec)/xyz_vol
            write(unit=fidgin,fmt='(5 (4x,F16.10))')x_mid,y_mid,z_mid,cin,xyz_vol
          enddo e6r0w3
          write(unit=fidgin,fmt=*)
        enddo e6r0w4
        write(unit=fidgin,fmt=*)
      enddo c4v9i4
      close(unit=fidgin)
    enddo r6o1e4

    e4r6o1: do ispec=1,nspec()
      if (debug) write(unit=fidlog,fmt=*)"! grtz profile for ",fspc(ispec)
      open(unit=fidgin,file="res/grtz-"//fspc(ispec)//"."//firun)
      write(unit=fidgin,fmt='("# UUID ",A32)')fuuid
      write(unit=fidgin,fmt='(A)')'# label rcoord  theta   zcoord conc  vol'
      write(unit=fidgin,fmt='(A)')'# unit  ang     radian  ang    molar ang3'
      i4c4v9: do g_=1,d3_zgrid_size ! zcoord
        !! z-range :: -d3_zmax +  d3_width/2 : (d3_zgrid_size-1)*d3_width + d3_width/2 - d3_zmax
        z_mid=(g_-1)*d3_width + d3_width/2 - d3_zmax
        w4e6r0: do k_=1,d3_rgrid_size ! thetacoord
          !! theta : -pi : pi - chunk
          t_mid=(k_-1)*2*pi/dble(d3_rgrid_size)-pi
          w3e6r0: do l_=1,d3_rgrid_size/2  ! rcoord
            !! r-range :: d3_width/2 : (d3_rgrid_size/2 - 1)*d3_width + d3_width/2
            r_mid=(l_-1)*d3_width + d3_width/2
            cin=tosi*hist_mean(grtz, l_, k_, g_,ispec)/d3_vol(l_)
            write(unit=fidgin,fmt='(5 (4x,F16.10))')r_mid,t_mid,z_mid,cin,d3_vol(l_)
          enddo w3e6r0
          write(unit=fidgin,fmt=*)
        enddo w4e6r0
        write(unit=fidgin,fmt=*)
      enddo i4c4v9
      close(unit=fidgin)
    enddo e4r6o1

    if (debug) write(unit=fidlog,fmt=*)"! gin profiles complete"
  endif k2l1f2

  k8a2o9: do ispec=1,nspec()
    cnduct(:,ispec)=0
    usecd=.true.

    ! FOUND BUG: Previously assumed xri > zocc
    zocspc=max(0.0D0,zocc-xri(ispec))

    strgn3=gz_bin(-zreg(ichan,ispec))
    strgn2=gz_bin(-zreg(ifilt,ispec))
    strgn1=gz_bin(-zocspc)
    enrgn1=gz_bin(zocspc)
    enrgn2=gz_bin(zreg(ifilt,ispec))
    enrgn3=gz_bin(zreg(ichan,ispec))
    ! FOUND BUG: check for valid intervals, 
    if (strgn3.gt.enrgn3) then
      usecd(ichan)=.false.
    end if
    if (strgn2.gt.enrgn2) then
      usecd(ifilt)=.false.
    end if
    if (strgn1.gt.enrgn1) then
      usecd(izlim)=.false.
    endif
    
    ! Output concentration profile and calculate slope conductance
    ! at the same time.
    open(unit=fidgzz,file="res/gz-"//fspc(ispec)//"."//firun)
      write(unit=fidgzz,fmt='("# UUID ",A32)')fuuid
      write(unit=fidgzz,fmt='(A)')'# label zcoord  conc   amount  vol   width'
      write(unit=fidgzz,fmt='(A)')'# units ang     molar  mole    ang3  ang'
      t2q4q1: do i=1,gz_max
        k8b4r5: if (hist_sum(gzhist,i,ispec).eq.0) then
          anz=0
          cz(i,ispec)=0

          if (i.ge.strgn3.and.i.le.enrgn3) then
            usecd(ichan)=.false.
            if (i.ge.strgn2.and.i.le.enrgn2) then
              usecd(ifilt)=.false.
              if (i.ge.strgn1.and.i.le.enrgn1) then
                usecd(izlim)=.false.
              endif
            endif
          endif
        else k8b4r5
          anz = hist_mean(gzhist,i,ispec)
          cz(i,ispec) = (tosi*anz)/gz_vol(i,ispec)
          if (i.ge.strgn3.and.i.le.enrgn3) then
            cnduct(ichan,ispec)=cnduct(ichan,ispec)+1.D0/anz
            if (i.ge.strgn2.and.i.le.enrgn2) then
              cnduct(ifilt,ispec)=cnduct(ifilt,ispec)+1.D0/anz
              if (i.ge.strgn1.and.i.le.enrgn1) then
                cnduct(izlim,ispec)=cnduct(izlim,ispec)+1.D0/anz
              endif
            endif
          endif
        endif k8b4r5
        write(unit=fidgzz,fmt='(5 (4x,F16.10))')gz_mid(i),cz(i,ispec)&
   &                   ,tosi*anz,gz_vol(i,ispec),gzwdth
      enddo t2q4q1
    close(unit=fidgzz)
    ! finalise conductances for specie
    do ireg=izlim,ichan
      if (usecd(ireg)) then
        if (dbc) then
          if(dbc_level.le.dbc_check) then
            if (dfeq(cnduct(ireg,ispec),0.D0)) then
              stop "Unexpected zero value"
            endif
          end if
        endif
        cnduct(ireg,ispec)=1.D0/cnduct(ireg,ispec)
      else
        cnduct(ireg,ispec)=0.D0
      endif
    enddo

    ! end ispec loop
  enddo k8a2o9

  ! ----- Output conductances -------------------------------

  u2j1y3: do ispec=idxcl(),nspec()
    l4i3o7: do ireg=izlim,ichan
      if (ksub.eq.1) then
        open(unit=fidgrg,file="res/s-"//freg(ireg)//"-"//fspc(ispec)//"."//firun)
          write(unit=fidgrg,fmt='("# UUID ",A32)')fuuid
          write(unit=fidgrg,fmt='(A)')'# label step  slope-conduct'
          write(unit=fidgrg,fmt='(A)')'# units count UNK'
      else
        open(unit=fidgrg,file="res/s-"//freg(ireg)//"-"//fspc(ispec)//"."//firun,position="APPEND")
      endif
          write(unit=fidgrg,fmt='(I5,4x,F16.12)')ksub,cnduct(ireg,ispec)
        close(unit=fidgrg)
    enddo l4i3o7
  enddo u2j1y3

  if (debug) write(unit=fidlog,fmt=*)"! gz profiles complete"

  ! Calculate radial distribution of ions near channel
  p0g3b9: if (calrdf) then
    i1n6j0: do ispec=1,nspec()
      w0g4c1: do jspec=1,nspec()

        ddi=dd_get(ispec,jspec)
        open(unit=fidrdf,file="res/rdf-"//fspc(ispec)//"-"//&
     &        fspc(jspec)//"."//firun)
        write(unit=fidrdf,fmt='("# UUID ",A32)')fuuid
        write(unit=fidrdf,fmt='(A)')'# label dist  nzlim  nfilt  dzlim   dfilt'
        write(unit=fidrdf,fmt='(A)')'# units ang   count  count  ang-3   ang-3'
        d4l0y1: do j=1,nrdf
          rlow=ddi+(j-1)*drdf
          rhigh=rlow+drdf
          ri=rlow+drdf/2
          vshell=4*pi*(rhigh**3-rlow**3)/3
          f8j0q0: if (hist_count(rdfhist(izlim,jspec,ispec)).ne.0) then
            rdfi1=hist_mean(rdfhist(izlim,jspec,ispec),j)
            rdf1vs=rdfi1/vshell
          else f8j0q0
            rdfi1=0
            rdf1vs=0
          endif f8j0q0
          r2w1i8: if (hist_count(rdfhist(ifilt,jspec,ispec)).ne.0) then
            rdfi2=hist_mean(rdfhist(ifilt,jspec,ispec),j)
            rdf2vs=rdfi2/vshell
          else r2w1i8
            rdfi2=0
            rdf2vs=0
          endif r2w1i8
          write(unit=fidrdf,fmt='(5 (4x,F16.10))')ri,rdfi1,rdfi2,rdf1vs,rdf2vs
        enddo d4l0y1
        close(unit=fidrdf)
      enddo w0g4c1
    enddo i1n6j0
  endif p0g3b9
  if (debug) write(unit=fidlog,fmt=*)"! rdf profiles complete"

  contains

  ! write mobile ion statistics
  subroutine amobsv
  double precision :: av, var
  write(unit=fidlog,fmt=*)"Statistics for mobile ions"
  write(unit=fidlog,fmt='(A7,5(A8))')"SPC","AVRG","VAR","DX","DY","DZ"
  j = 0 ! using j as particle count 
  do ispec=1,nspec()
    if (localized(ispec)) then
      do i=1,ni(ispec)
        av=amobdl(i+j)/dble(athist)
        var=(amobvr(i+j)/dble(athist)-sqr(av))
        write(unit=fidlog,fmt='(A2,"[",I3,"]",5(F8.2))')fspc(ispec),i,av,var,&
           & amobdx(i+j)/dble(athist),amobdy(i+j)/dble(athist),amobdz(i+j)/dble(athist)
      end do
    endif
    j = j + ni(ispec)
  enddo
  end subroutine amobsv

  ! Update output file with the number of ions in the filter 
  subroutine u9c6u7
  use spec
  implicit none
  double precision, dimension (nrgnmx) :: nreg_
  double precision :: asuba
  ! ----- Calculate bulk concentration and ---------------------
  ! ----- number of ions in the filter -------------------------

  n7n4e6: do ispec=idxcl(),nspec()
    nreg_(izlim) = hist_mean(anin, izlim, ispec)
    nreg_(ifilt) = hist_mean(anin, ifilt, ispec)
    nreg_(ichan) = hist_mean(anin, ichan, ispec)
    ! abulk used here to get data from the 'bulk' sampling sub-zone
    nreg_(ibulk) = rs_mean(abulk(ispec))

    ! Output population profiles for ions
    f8g3o6: do ireg=izlim,ibulk
       t2b3u5: if (ksub.gt.1) then
          open(unit=fidarg,file="res/a-"//freg(ireg)//"-"//fspc(ispec)//"."//firun,position="APPEND")
          asuba =ksub*nreg_(ireg) - (ksub-1)*anaca(ireg,ispec)
          anaca(ireg,ispec)=nreg_(ireg)
       else t2b3u5
          open(unit=fidarg,file="res/a-"//freg(ireg)//"-"//fspc(ispec)//"."//firun)
          write(unit=fidarg,fmt='("# UUID ",A32)')fuuid
          write(unit=fidarg,fmt='(A)')'# label step  aver_cml  aver_step'
          write(unit=fidarg,fmt='(A)')'# units count count     count'

          anaca(ireg,ispec)=nreg_(ireg)
          asuba=nreg_(ireg)
       endif t2b3u5
       write(unit=fidarg,fmt='(I5,2 (4X,F16.10))')ksub,anaca(ireg,ispec),asuba
       close(unit=fidarg)
    enddo f8g3o6
  enddo n7n4e6
  end subroutine u9c6u7
  end subroutine saves

  !--------------------------------------------------
  ! Report concentrations and chem potentials
  subroutine chemical_potential_report(fid)
    use strngs
    use spec
    implicit none
    integer, intent(in) :: fid
    double precision :: bulk_conc
    integer :: ispec, igc
    character(15) :: fltout1, fltout2, fltout3, fltout4
    write(unit=fid,fmt=*)"ESTIMATED error in CHEMICAL POTENTIALS based on current ensemble data"
    write(unit=fid,fmt='(72("-"))')
    write(unit=fid,fmt='(1X,A6,4(1X,A15))')"SPECIE","TARGET CONC.","CURRENT CONC.","CHEM. EXCESS","ERROR."
    do ispec=idxcl(),nspec()
      bulk_conc = bulk_concentration(ispec)
      call str (ctargi(ispec), fltout1)
      call str (chexi(ispec), fltout3)
      if (dfeq(bulk_conc,0.D0)) then
        ! Skip on zero bulk concentration
        write(unit=fid,fmt='(1X,A6,4(1X,A15))')fspc(ispec),fltout1,'0.0',fltout3,'-----'
      else
        call str (bulk_conc, fltout2)
        ! log(bulk_conc/ctargi(ispec)) is the energy difference
        ! required to change the concentration from current to target
        ! which is the error in the mu_ex ( dE = -log(Ct/Cc) )
        call str (log(bulk_conc/ctargi(ispec)), fltout4)
        write(unit=fid,fmt='(1X,A6,4(1X,A15))')fspc(ispec),fltout1,fltout2,fltout3,fltout4
      endif
    enddo
    write(unit=fid,fmt='(1X,A6,4(1X,A15))')"SALT","TARGET CONC.","CURRENT CONC.","CHEM. POT."
    do igc=1,nsalt()
      ispec  = isalt(igc)
      bulk_conc = bulk_concentration(ispec)
      call str (ctargs(igc), fltout1)
      call str (chemps(igc), fltout3)
      if (dfeq(bulk_conc,0.D0)) then
        ! Skip on zero bulk concentration
        write(unit=fid,fmt='(1X,A6,3(1X,A15))')fsalt(igc),fltout1,'0.0',fltout3
      else
        call str (bulk_conc, fltout2)
        write(unit=fid,fmt='(1X,A6,3(1X,A15))')fsalt(igc),fltout1,fltout2,fltout3
      endif
    enddo
  end subroutine chemical_potential_report

  !--------------------------------------------------
  ! report concentrations
  subroutine concentration_report(fid)
    use conf
    use geom
    use spec
    use strngs
    implicit none
    integer, intent(in) :: fid
    integer :: ispec
    double precision :: bulk_conc
    double precision, dimension(nspcmx) :: cfilt,cchan,ctot,czlim
    double precision, dimension(nrgnmx) :: charg_ ! region charge
    cfilt=0
    cchan=0
    ctot=0
    czlim=0
    charg_=0

    ! ----- Calculate bulk concentration and ---------------------
    write(unit=fid,fmt='(72("-"))')
    write(unit=fid,fmt='(10X,"|",5X,"[BULK]",1X,"|",3X,"CHANNEL",1X,"|",3X,"FILTER ",1X,"|",3X,"CENTER ",1X,"|",7X,"TOTAL")')
    ! ----- number of ions in the filter -------------------------

    u1k8f1: do ispec=1,nspec()
      r7f9c4: if (isfree(ispec).or.flexible(ispec)) then
        czlim(ispec)=hist_mean(anin,izlim,ispec)
        cfilt(ispec)=hist_mean(anin,ifilt,ispec)
        cchan(ispec)=hist_mean(anin,ichan,ispec)
        ctot(ispec) =hist_mean(anin,ibulk,ispec)
      else r7f9c4
        czlim(ispec)=ni(ispec)
        cfilt(ispec)=ni(ispec)
        cchan(ispec)=ni(ispec)
        ctot(ispec)=ni(ispec)
      endif r7f9c4
      bulk_conc = bulk_concentration(ispec)
      charg_(izlim)=charg_(izlim)+czlim(ispec)*xz(ispec)
      charg_(ifilt)=charg_(ifilt)+cfilt(ispec)*xz(ispec)
      charg_(ichan)=charg_(ichan)+cchan(ispec)*xz(ispec)
      charg_(ibulk)=charg_(ibulk)+xz(ispec)*bulk_conc

      b0f1x5: if(isfree(ispec)) then
        write(unit=fid,fmt='(1x,a2,"  ion: |",1X,F10.5,1X,3("|",1X,F9.3,1X),"|",F7.3,1X,I4)')&
             &      fspc(ispec),bulk_conc,cchan(ispec),cfilt(ispec),czlim(ispec),ctot(ispec),ni(ispec)
        write(unit=fid,fmt='(2x," Cell av|",1X,F10.5,1X)')tosi*ctot(ispec)/vtotal(ispec)
        write(unit=fid,fmt='(2x,"Cell cur|",1X,F10.5,1X)')tosi*ni(ispec)/vtotal(ispec)
      else b0f1x5
        write(unit=fid,fmt='(1x,a2,"  ion: |",12X,3("|",1X,F9.3,1X),"|",F7.3,1X,I4)')&
             &      fspc(ispec),cchan(ispec),cfilt(ispec),czlim(ispec),ctot(ispec),ni(ispec)
      endif b0f1x5
    enddo u1k8f1

    write(unit=fid,fmt='(72("-"))')
    write(unit=fid,fmt='("  Charge: |",1X,F10.5,1X,3("|",1X,F9.3,1X))')&
         &    charg_(ibulk)*vbulk/tosi,charg_(ichan),charg_(ifilt),charg_(izlim)

  end subroutine concentration_report

  ! ------------------------------------------------------------------
  ! Calculate theta bin of this angle 
  ! 
  function theta_bin(x, y, theta0)
  implicit none
  double precision, intent(in) :: x,y,theta0
  double precision :: theta
  integer :: theta_bin
  theta=atan2(x,y)-theta0
  if (theta.lt.0.D0) theta = theta + 2*pi
  ! bin numbering is -pi == 0, pi == 2*nrgr
  theta_bin=max(1,ceiling((theta+pi)*d3_rgrid_size/(2*pi)))
  ! Check is in bounds
  if (theta_bin.le.0) theta_bin=theta_bin + d3_rgrid_size
  if (theta_bin.gt.d3_rgrid_size) theta_bin=theta_bin - d3_rgrid_size
  if (dbc) then
    if (theta_bin.le.0) then
      write(unit=fidlog,fmt=*)"Calculated bin [",theta_bin,"] for angle=",theta," is less than 0"
      stop "Calculated bin is less than 0"
    elseif (theta_bin.gt.d3_rgrid_size) then
      write(unit=fidlog,fmt=*)"Calculated bin [",theta_bin,"] for angle=",theta," is greater than ",d3_rgrid_size
      stop "Calculated bin is greater than grid size"
    endif
  endif
  end function theta_bin

  ! ------------------------------------------------------------------
  ! reset the accumulators
  !
  subroutine zeroac
  use spec
  use patch
  use simstate
  implicit none
  integer :: i, j
  amobdl=0
  amobvr=0
  amobdx=0
  amobdy=0
  amobdz=0
  if (do_electrostatic()) then
    if (.not.is_homogeneous()) ah=0
    call rs_reset (acharge)
  endif
  aocc=0
  call hist_reset (anin)
  athist = 0
  ataccu = 0
  if (calgin) then
    call hist_reset(gin)
    call hist_reset(gxyz)
    call hist_reset(grtz)
    call hist_reset(gzhist)
  endif
  if (calrdf) then
    do i=1,nspec()
      do j=1,nspec()
        call hist_reset(rdfhist(izlim,j,i))
        call hist_reset(rdfhist(ifilt,j,i))
      enddo
      call rs_init(abulk(i))
    enddo
  endif
  ninj=0
  end subroutine zeroac

end module accum

subroutine gzaver(ibin,ispec,naver)
use accum
implicit none
integer, intent(in) :: ibin, ispec
double precision, intent(out) :: naver
naver=gz_av(ibin,ispec)
end subroutine gzaver
