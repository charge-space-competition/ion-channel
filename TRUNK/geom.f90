
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


!  SIMULATION GEOMETRIC PARAMETERS (non-volatile)
!
!  ASSUMPTIONS
!
!    * SHAPES (channel protein, membrane, simulation box)
!      * All shapes are symmetric with respect to z=0 plane
!      * All shapes are radially symmetric with respect to z-axis (radial direction)
!    * REGIONS (centre zone=1, filter=2, channel=3, bulk=4)
!      * regions 1 and 2 are cylinders with the same radii
!      * region 3 need not be a cylinder (ie can include vestibule arches)
!        * membrane is assumed to have identical z-extent as region 3
!      * region 4 is bulk region---the whole volume excluding membrane
!        * membrane extends to boundary in radial direction
!
!  INPUT SECTION
!
!    geom
!     zl1 1.23
!     zl4 1.23
!     rl1 1.23
!     rl4 1.23
!     rl5 1.23
!     rlvest 1.23
!     rlcurv 1.23
!     zlimit 1.23
!     ntarg 234
!    end
!                 ^
!                 ^ rl5
!                 |
!        <- zl3  >|
!        ---------|---------
!       /|        |  ^rl4   \
!      / |        |  |   ____\
!     |<-> = rlcurv  |  ^rl3  |
!     |<--> = rlvest |  |_____|
!     \   |       |  |  ^rl2  /
!      \  |       |  |  |    /
!       \ |       |  |  |   /
!        \|<-zl1->|  |  |  /
!         ----------------- 
!           ^ rl1 |  |  |
!      <-----zl2->|  |  |
!   ________|_____|_ |__|______ z-axis
!<<-------------->| zl4
!          ??<--->|(z=0)
!            zlimit
!
!     * zl1 = half-length of region 2
!     * zl4 = nominal half-length of simulation box
!     * zlimit = limit of structural ions 
!       ** (zlimit <= zl1)
!     * rl1 = radius of regions 1 and 2
!     * rl4 = outer radius of protein
!     * rl5 = nominal radius of simulation box
!     * rlvest= radius of vestibule arc
!     * rlcurv= radius of protein outer arc
!       ** (rlvest + rlcurv <= rl4 - rl1)
!       ** (rlcurv <= zl1 + rlvest)
!
!     * ntarg = number of 'free' particles in box
!       ** zl4 and rl5 will be adjusted when ntarg is set to give
!         bulk volume that achieves a specified target concentration 
!         with given number of particles.
!
!     Derived from above:
!     * half-length of region 3 is zl1 + rlvest
!     * maximum radius of region 3 is rl1 + rlvest
!
! LENGTH COMPARISONS
!
! When testing if a particle is too close to an object the comparison
! used is less-than or greater-than. Conversely to test if a particle
! is far-enough away from another object the less-equal and greater-
! equal comparisons are used.  For tests if a particle is within a 
! region use less-than comparisons.
!
! REGIONS AND VOLUMES
!
! To be totally clear how the above definitions translate into 
! the regions and volumes used in the program is explained here.
!
! VOLUMES
!
! Volumes used in calculations are always defined by the extents of
! the ranges of z and r. For example if z and r are in range 0-1
! then the cylinder volume is \pi r^2 z = \pi.  This is regardless 
! of any possible obstructions within this range that would invalidate
! a particle at any particular z and r.  This applies both when we are
! generating the z and r randomly in the case of MC moves or sampling
! the system where the volume used is deduced from the allowed ranges
! of z and r.
!
! REGIONS
!
! Region 1: (+/- zlimit)
!
! This is a region that is independent of the rest of the geometry.
! It is used to specify the region available to structural ions. It is
! an explicit assumption in the code that Region 1 is completely within
! Region 2.
!
! Region 2: (+/- zl(1))
!
! This is the region of the channel excluding the vestibules.
!
! Region 3: (+/- zl(2))
!
! This is the entire channel.
!
! Region 4: (+/- zl(4))
!
! This is the entire simulation cell.
!
! NOTE:
!
! There is one more 'region' refered to in the channel simulation,
! the "bulk-sampling" region. This should not be confused with
! Region 4 which is sometimes refered to as the 'bulk' region and
! has the index parameter 'ibulk'.  The "bulk-sampling" space
! is the space outside the membrane with all dimensions
! reduced so that it is reasonable to assume the
! effects of the simulation cell boundaries are
! small. Thus we assume the ensemble of particles within this space 
! will have properties similar to the ensemble in the
! bulk simulation (which has PBC). The variables zbulk1, zbulk2
! rbulk and vbulk are the geometric extents of this space.
!
! NOTE ON POSITION GENERATION
!
! In the event of a disjoint interval (eg volume for jump-out) using
! one comparison leads to asymmetric intervals. For example if the 
! random number is [-1:1] then splitting using < 0 gives intervals [-1:0) [0:1]. 
! Therefore two comparisons are required < 0 and > 0 with a new random
! number selected until one is found in either interval.  This gives
! intervals [-1:0) and (0:1]. Code like the following could be used.
! """
! do
!   rnd = ranff() ...
!   if (rnd.lt.0.D0.or.rnd.gt.0.D0) exit
! enddo
! """
!
! Generating two closed intervals is more problematic.
!
!
! REGION BOUNDARIES
!
! The region boundaries generally used are adjusted to improve sampling
! efficiency and smoothness across the geometric boundaries.
!
! Region 1: z [-zlimit:zlimit], r [0:rl(1)-r(spec)]
! Region 2: z [-zl(1)-r(spec):zl(1)+r(spec)], r [0:rl(1)-r(spec)]
! Region 3: z [-zl(2)-r(spec):zl(2)+r(spec)], r [0:rl(2)]
! Region 4: z [-zl(4)+r(spec):zl(4)-r(spec)], r [0:rl(4)-r(spec)]
!
! ZREG/RREG/VREG: These define insertion boundaries and volumes
! in terms of the particle center-point and are therefore specie
! dependent.  A particle is considered to be inserted into
! a region if the center-point is in these bounds.
!
! NOTE: Model is constructed with symmetry around 0.0
!
! ZLJIN == ZREG(3): The position of a particle center-point that
! allows any part of the particle to be inside the channel.
!
! ZLJOUT == ZREG(4)-ZREG(3): The width of the cell sub-part where
! no part of a particle can be inside the channel
!
! GZ STATISTICAL DISTRIBUTION
!
! * New properties
! ** The bin z-axial parameters are common for all species
! ** The radial and volume are specie dependent
! ** The volumes used are those accessible to specie center-points.
! ** Inclusion in a bin is dependent on half-open intervals, with open
! end towards ..TODO.. negative.
! 
! VGZ/VJIN: These are the slab volumes that a particle may
! occupy. The volumes are defined by the spatial extent of the
! particle center-point and are specie dependent.
!
! NEW VOLUME OF ROTATION FUNCTION(S)
!
! VOLROT(Z0,Z1,RC,ZC,R): This function determines the volume of
! rotation under the arc centered at RC,ZC with radius R between
! Z0 and Z1. This is limited to: Z0 and Z1 being in range [ZC:ZC+R],
! Z0 < Z1 and, R being in range (0:RC]. It calculates the volume
! using an analytical integration of the volume of rotation
! integral \pi \intg_Z0^Z1 (f(Z))^2 dZ where f(Z) is equation of
! arc in the z,r plane.
!
! VOLROT(RC,R): This function determines volume of rotation as
! for 5 parameter version over the range Z0=ZC-R,Z1=ZC+R. This is
! the volume under the toroid defined by the arc. Using the
! entire region under arc allows simplification of the integral making
! this calculation less expensive than the more general one.
!
! OTHER METHODS
!
! WALL: The centrepoint of a particle must be further than
! its radius from any geometric obstruction. Furthermore
! structural ions must have no part lieing outside the region 
! defined by +/- zlimit (ie center-point limited by zlimit-r(spec))
!
! VOLBLK: volume used in calculating bulk properties within
! the non-bulk simulation cell, where inclusion in a bulk sub-region is
! based on center-point being in bulk sub-cell.
!
! ----------------------------------------------------------------------
! PUBLIC ATTRIBUTES AND METHODS
!
!  cubmov
!     generate new x, y, z coordinates in a given box
!
!  defgeo(donut)
!     fill geomdf structure with the complete geometry
!     parameters needed for generating ICC patchs
!
!  ecgeom
!     write out the geometry options
!
!  gz_hi, gz_lo, gz_mid
!     specie independent z-coordinates of the histogram bins.
!
!  gz_max
!     the maximum bin number, range is 1 to gz_max.  Note that
!     the number of bins is determined from the maximum specie
!     specific extent so bins at either extreme will not be
!     accessible to all species.
!
!  gz_vol
!     specie specific volume of histogram bins.  This volume is
!     defined by the space that can be occupied by the center-points
!     of particles of a specie.
!
!  gzwdth
!     specie independent histogram bin width.
!
!  ingeom(radius,zcoord,rcoord,chonly,overlap) 
!     test if particle of given radius overlaps any walls
!
!  jmpfin
!     generate new x and y coordinates from a given r coordinate
!
!  jmpmov
!     generate new z and r coordinates in a given region
!
!  rbulk
!     radial dimension for the "bulk-sampling" space within the
!     channel simulation
!
!  rdgeom(fid,name,value,ios)
!     read in 'geom' section of input file
!
!  rfgeom
!     finalise geometry after all sections of the input file
!     have been read
!
!  rl(1-5)
!     protein and MC box radial dimensions
!
!  rlvest(), rlcurv()
!     arc radii for the inner and outer sections of protein
!
!  rreg
!     Specie specific r-dimension of a region
!
!  vbulk
!     volume defined by zbulk1, zbulk2 and rbulk; the "bulk-sampling" 
!     space within the channel simulation
!
!  vin,vout
!     specie specific volumes used in preference sampling
!     of moves into and out of channel. Defined by (zljin,
!      rreg(3)) and (zljout, rreg(4)).
!
!  volrot
!     calculate the volume under an arc using volume
!     of rotation formula.
!
!  vreg
!     Specie specific volume defined by zreg and rreg
!
!  vtotal
!     Specie specific total voume of the simulation.
!
!  wall(ispec,zcoord,rcoord,overlap)
!     test if particle of given specie overlaps any walls
!
!  zbulk1,zbulk2
!     inner and outer z-coords for the "bulk-sampling" space within the
!     channel simulation
!
!  zl(1-4)
!     protein and MC box z-dimensions
!
!  zljin,zljout
!     specie specific half-width of regions for preference sampling
!     of moves into and out of channel. zljin=zreg(3), zljout=
!     zreg(4)-zreg(3)
!
!  zlimit()
!     z-dimension of the central filter
!
!  zreg
!     Specie specific z-dimension of a region
!
! ----------------------------------------------------------------------
module geom
  use const
  implicit none
  private
  ! dimensions of protein donut and MC space
  ! r-coordinates of protein donut
  double precision, private, dimension(8) :: q4b3p4

  ! left z-coordinate of bulk
  double precision, public :: zbulk1 = 0.D0
  ! right z-coordinate of bulk
  double precision, public :: zbulk2 = 0.D0
  ! radial bulk
  double precision, public :: rbulk = 0.D0
  ! volume of bulk
  double precision, public :: vbulk = 0.D0
  ! length and volume for simulation of bulk
  double precision, private :: x8h5t7 = 0.D0, s9d1h4 = 0.D0
  ! have we read in the geometry?
  logical, private :: irdgeo = .false.
  ! whether to use the older insertion region definition or newer
  logical, private :: useold = .true.
  ! region names
  character(4), public, dimension(nrgnmx) :: freg

  ! target number of particles for given ionstr and volume
  integer, private :: ntrg = 0

  ! Region boundaries and volumes
  double precision, public, dimension(nrgnmx,nspcmx) :: zreg = 0.D0
  double precision, public, dimension(nrgnmx,nspcmx) :: rreg = 0.D0
  double precision, public, dimension(nrgnmx,nspcmx) :: vreg = 0.D0
  ! Preference move (jump-in/jump-out of channel) data
  double precision, public, dimension(nspcmx) :: zljin = 0.D0
  double precision, public, dimension(nspcmx) :: zljout = 0.D0
  double precision, public, dimension(nspcmx) :: vin = 0.D0
  double precision, public, dimension(nspcmx) :: vout = 0.D0
  double precision, public, dimension(nspcmx) :: vtotal = 0.D0

  public :: cubmov, jmpmov, jmpfin, volrot

  ! -------------------------------------------------------------- 
  ! PROTEIN DEFINITION USED BY PATCH MODULE
  !
  ! Type for transferring protein shape information needed by
  ! the patch module.  This module fills this type with detailed
  ! geometry data of the patches used in the 'patch' module. This
  ! information is only required during the initialisation of
  ! the 'patch' module so is only created temporarily when required.
  !
  ! The data in this object is obtained by passing an object of this
  ! type to the 'defgeo' subroutine
  !
  type geomdf
     double precision, dimension(8,300) :: dfia,tta1,tta2
     double precision, dimension(6,300) :: dfil,zzl1,zzl2,rrl1,rrl2
     integer, dimension(8,300) :: nfia
     integer, dimension(6,300) :: nfil
     integer, dimension(20) :: ilsign
     double precision, dimension(8) :: za0,ra0,ra,ta1,ta2,uasign
     double precision, dimension(6) :: ulsign,zl1,rl1,zl2,rl2,tgalfa,alfa
     integer, dimension(8) :: nta
     integer, dimension(6) :: nzl,nrl
     integer :: narch,nline,nwall,nwall2,nunif,ncent
  end type geomdf

  public :: geomdf

  ! -------------------------------------------
  ! Calculate volumes under an arc
  interface volrot
     module procedure vrot5, vrot2
  end interface volrot

  ! -------------------------------------------
  ! Random number in interval
  interface rndint
     module procedure rndit4, rndit2
  end interface rndint

  ! -------------------------------------------
  ! Jump move (sets x,r)
  !
  ! This interface generates a random x,r position information
  ! within an interval.  One method takes a radius and a two
  ! parameter interval.  The other takes a radius and four parameter
  ! interval, with the result in the interval of the first two
  ! parameters excluding the interval of the second two paramters.
  interface jmpmov
     module procedure jmpmv4, jmpmv2
  end interface jmpmov

  ! ------------------------------------------------------------------
  ! Cubic move (sets x,y,z,r)
  !
  ! Two interfaces exist, the first generates a position within a
  ! cube centred on the current position, the second generates a
  ! position in a cube centred on a second position.
  interface cubmov
     module procedure cubmv1, cubmv2
  end interface cubmov

  public :: lenblk, inchan, ingeom, inregn, dispbc, disbox, rl, zl, zlimit
  public :: ntargt, rlvest, rlcurv, volblk  

  public :: ecgeom, rdgeom, rfgeom, wall, defgeo

  ! DEBUG routines
  public :: dumpgdfn, initgdfn


  ! --------------------------------------------------------------
  ! GZ Z-AXIAL HISTOGRAM BIN PARAMETERS
  !

  ! [INPUT] The input target width for bins in the 'gz' distribution
  double precision, public :: gzwdth

  ! The lowest bin edge
  double precision, private :: gzlow_ = 0.D0

  ! The total number of bins
  integer, public :: gz_max = 0

  ! Half the number of bins
  integer, private :: nhbin

  ! The centre-point of each histogram bin in 'gz' statistical distributon
  ! (specie independent)
  double precision, public, allocatable, dimension(:) :: gz_mid

  ! Volumes for bins/slabs in the z-axial 'gz' statistical distribution
  ! (specie dependent)
  double precision, public, allocatable, dimension(:,:) :: gz_vol

  ! enquiry functions/subroutines
  public :: gz_lo,gz_hi,gz_bin

contains

  ! -------------------------------------------
  ! Random number in interval
  ! Generate a random number that is between outlft and outrht but
  ! not between inlft and inrht
  !
  !  1  --  3 - 4 -- 6
  !  width 6-1+3-4 = 4
  !  guess=ran*width+1
  !  if (guess>3) guess+(4-3)
  double precision function rndit4(outlft,outrht,inlft,inrht)
    implicit none
    double precision, intent(in) :: outlft,outrht,inlft,inrht
    double precision :: width
    width=outrht-outlft+inlft-inrht
    do ! reject random numbers on inlft boundary
       rndit4=width*ranff()+outlft
       if (dfeq(rndit4,inlft)) cycle
       exit
    enddo
    if (rndit4.gt.inlft) then
       rndit4=rndit4+inrht-inlft
    endif
  end function rndit4

  ! -------------------------------------------
  ! Random number in interval
  ! Generate a random number that is between left and right
  double precision function rndit2(left,right)
    implicit none
    double precision, intent(in) :: left, right
    rndit2=(right-left)*ranff()+left
  end function rndit2

  ! -------------------------------------------
  ! Distance in non-periodic system
  !
  ! Simple distance between two points.
  pure double precision function disbox(a,b,c,x,y,z)
    implicit none
    double precision, intent(in) :: a,b,c,x,y,z
    disbox=dsqrt(sqr(a-x)+sqr(b-y)+sqr(c-z))
  end function disbox

  ! ------------------------------------------------------------------
  ! Displacement in periodic cube.
  !
  ! Minimum distance of two points accounting for periodic boundary
  ! conditions.
  pure double precision function dispbc(a,b,c,x,y,z)
    implicit none
    double precision, intent(in) :: a,b,c,x,y,z
    dispbc=dsqrt(sqr(do_pbc(a-x))+sqr(do_pbc(b-y))+sqr(do_pbc(c-z)))
  end function dispbc

  ! ------------------------------------------------------------------
  ! Radial move (sets x,y,z,r)
  !
  ! Generate a random position in a cube that extends from 0 to
  ! maxmov in all directions.  This method is used to perform a
  ! jump move in the simulation of a bulk solution.
  subroutine cubmv1(x, y, z, r, maxmov)
    implicit none
    double precision, intent(out) :: x, y, z, r
    double precision, intent(in) :: maxmov
    x=ranff()*maxmov
    y=ranff()*maxmov
    z=ranff()*maxmov
    r=dsqrt(sqr(x)+sqr(y))
  end subroutine cubmv1

  ! ------------------------------------------------------------------
  ! Radial move (sets x,y,z,r)
  !
  ! Generate a random position in the cube defined by the point
  ! x0, y0 and z0 that extends in all directions by interval
  ! left, right.  This method is used to perform a move that is
  ! a small displacement from the original (x0,y0,z0) position.
  subroutine cubmv2(x, y, z, r, x0, y0, z0, left,right)
    implicit none
    double precision, intent(out) :: x, y, z, r
    double precision, intent(in) :: x0, y0, z0, left, right
    x=x0+rndint(left,right)
    y=y0+rndint(left,right)
    z=z0+rndint(left,right)
    r=dsqrt(sqr(x)+sqr(y))
  end subroutine cubmv2

  ! ------------------------------------------------------------------
  ! Jump move (sets z,r)
  !
  ! Generate random Z, R in cylinder defined by left, right and radius
  ! but excluding interval between inlft and inrght
  subroutine jmpmv4(z, r, left, right, inlft, inrght, radius)
    implicit none
    double precision, intent(out) :: z, r
    double precision, intent(in) :: left, right, inlft, inrght, radius
    z=rndint(left, right, inlft, inrght)
    r=dsqrt(ranff())*radius
  end subroutine jmpmv4

  ! ------------------------------------------------------------------
  ! Jump move (sets z,r) 
  !
  ! Generate random Z, R in cylinder defined by left, right and radius
  subroutine jmpmv2(z, r, left, right, radius)
    implicit none
    double precision, intent(out) :: z, r
    double precision, intent(in) :: left, right, radius
    z=rndint(left, right)
    r=dsqrt(ranff())*radius
  end subroutine jmpmv2

  ! ------------------------------------------------------------------
  ! Generate a random Y and X from R
  !
  ! Use this method to generate an x and y position randomly on
  ! a circle of radius R in the x-y plane.
  subroutine jmpfin(x, y, r)
    implicit none
    double precision, intent(in) :: r
    double precision, intent(out) :: x, y
    double precision :: phi
    phi = ranff() * 2*pi
    x=r*dcos(phi)
    y=r*dsin(phi)
  end subroutine jmpfin

  ! ------------------------------------------------------------------
  ! Effective distance in the bulk box with periodic boundary
  ! conditions
  pure double precision function do_pbc(a)
    implicit none
    double precision, intent(in) :: a
    do_pbc=a-dnint(a/lenblk())*lenblk()
  end function do_pbc

  ! ------------------------------------------------------------
  ! Determine the gz bin the item belongs in
  integer function gz_bin(rzi)
    implicit none
    double precision, intent(in) :: rzi
    ! solution
    gz_bin = max(1,ceiling((rzi - gzlow_)/gzwdth))
    ! ensure in range
    gz_bin = max(1, min(gz_bin, gz_max))
  end function gz_bin

  ! ------------------------------------------------------------
  ! Lowest z-coordinate of a bin
  double precision function gz_lo(idx)
    implicit none
    integer, intent(in) :: idx
    gz_lo=gz_mid(idx)-(gzwdth/2.D0)
  end function gz_lo

  ! ------------------------------------------------------------
  ! Highest z-coordinate of histogram bin
  double precision function gz_hi(idx)
    implicit none
    integer, intent(in) :: idx
    gz_hi=gz_mid(idx)+(gzwdth/2.D0)
  end function gz_hi

  ! ------------------------------------------------------------------
  ! Test if rznew is in the channel region.
  !
  ! @param rznew : z-coord of particle
  logical function inchan(rznew)
    implicit none
    double precision, intent(in) :: rznew
    inchan=(-zl(2).le.rznew.and.rznew.le.zl(2))
  end function inchan

  ! ------------------------------------------------------------
  ! Which region is a particle in?
  !
  ! Determine the lowest numbered region a particle is in. Note
  ! region 2 contains region 1 and region 3 contains region 2.
  ! This requires the programmer to account for the inclusions.
  ! For example membership of region 2 and 3 when ireg==1.
  ! 
  subroutine inregn(z,ispec,ireg)
    use spec
    implicit none
    double precision, intent(in) :: z
    integer, intent(in) :: ispec
    integer, intent(out) :: ireg
    integer :: idx
    double precision :: z_
    ! handle structural ions
    if (.not.isfree(ispec)) then
       ireg=izlim
       return
    endif
    ! now non-structural ions
    z_=abs(z)
    ireg = ibulk ! assume in bulk if not in channel
    do idx=izlim,ichan
       if (z_.le.zreg(idx,ispec)) then
          ireg=idx
          return
       endif
    enddo
  end subroutine inregn

  ! ------------------------------------------------------------------
  ! CHECK FOR PARTICLE - PROTEIN/MEMBRANE OVERLAP
  !
  ! Test if rznew and r2inew are a valid position for a particle
  ! with the given radius.  This method is a more flexible version
  ! of the 'wall' method.  The main difference is this method is
  ! independent of the 'spec' module and allows you to remove the
  ! restriction that structural ions may only exist in the filter
  ! region.  This restriction is valid in a standard simulation
  ! where the 'wall' method should always be used.  Alternatively,
  ! this method can be used to probe the geometry with spheres of
  ! any size or with structural ions anywhere in the system.
  !
  ! @param ri    : radius of particle
  ! @param rznew : z-coord of particle
  ! @param r2new : radial distance of particle
  ! @param chonly : if the particle can only exist in the channel
  ! @param ovrlap : (OUT) whether the particle overlaps a wall

  subroutine ingeom(ri,rznew,r2inew,ichnly,ovrlap)
    implicit none
    double precision, intent(in) :: ri,rznew,r2inew
    logical, intent(in) :: ichnly
    logical, intent(out) :: ovrlap

    ! Locals
    double precision drznew,rsq
    double precision zl10,zl20,zl40,rl50,rl4,rl10  ! Environment geom data
    if (dbc) then
       if (.not.irdgeo) stop "Error: geometry test called before defining geometry"
    endif
    ovrlap=.false.
    drznew=dabs(rznew)
    ! channel cylinder extreme in radial direction
    ! Point at which particle only need to worry about central cylinder
    rl10=rl(1)-ri
    ! z-axis extreme for structural ions
    zl10=zlimit()-ri
    i4l7f5: if (ichnly) then
       ovrlap = (r2inew.gt.rl10.or.drznew.gt.zl10)
       return
    else i4l7f5
       ! Point at which overlap with protein not possible
       zl20=zl(2)+ri
       if (drznew.ge.zl20) then
          ! in bulk
          ! check if outside simulation cell
          ! Extreme of bulk in r
          rl50=rl(5)-ri
          ! Extreme of bulk in z
          zl40=zl(4)-ri
          ovrlap = (r2inew.gt.rl50.or.drznew.gt.zl40)
          return
       endif
       ! BELOW can assume in channel
       if (r2inew.le.rl10) then
          ! inside minimum channel radius
          return
       endif
       ! BELOW can assume in channel but outside rl10
       if (drznew.lt.zl(1).or.r2inew.gt.rl(2)) then
          ! in inner cylinder but outside rl10 
          ! or in vestibule outside rl(2)
          ovrlap=.true.
          return
       endif
       ! BELOW can assume in vestibule with r2inew > rl10
       ! calculate (sqr) of arc radius (extended by specie radius)
       rl4=sqr(rlvest()+ri)
       ! calculate distance from particle to centerpoint
       rsq=sqr(rl(2)-r2inew)+sqr(drznew-zl(1))
       ovrlap = (rsq.lt.rl4)
    endif i4l7f5
  end subroutine ingeom

  ! --------------------------------------------------
  ! Length of simulation of bulk box
  pure double precision function lenblk()
    implicit none
    lenblk=x8h5t7
  end function lenblk

  ! --------------------------------------------------
  ! Target number of particles
  pure integer function ntargt()
    implicit none
    ntargt=ntrg
  end function ntargt

  ! --------------------------------------------------
  ! volume of simulation of bulk box
  pure double precision function volblk()
    implicit none
    volblk=s9d1h4
  end function volblk

  ! --------------------------------------------------
  ! Get a radial dimension of system
  double precision function rl(idx)
    implicit none
    integer, intent(in) :: idx
    select case (idx)
    case(1)
       rl=q4b3p4(3)
    case(2)
       rl=q4b3p4(3)+q4b3p4(6)
    case(3)
       rl=q4b3p4(4)-q4b3p4(7)
    case(4)
       rl=q4b3p4(4)
    case(5)
       rl=q4b3p4(5)
    case default
       stop "invalid index for array rl"
    end select
  end function rl

  ! Get a z-coord dimension of system
  double precision function zl(idx)
    implicit none
    integer, intent(in) :: idx
    select case (idx)
    case(1)
       zl=q4b3p4(1)
    case(2)
       zl=q4b3p4(1)+q4b3p4(6)
    case(3)
       zl=q4b3p4(1)+q4b3p4(6)-q4b3p4(7)
    case(4)
       zl=q4b3p4(2)
    case default
       stop "invalid index for array zl"
    end select
  end function zl

  pure double precision function rlvest()
    implicit none
    rlvest=q4b3p4(6)
  end function rlvest

  pure double precision function rlcurv()
    implicit none
    rlcurv=q4b3p4(7)
  end function rlcurv

  pure double precision function zlimit()
    implicit none
    zlimit=q4b3p4(8)
  end function zlimit

  ! -----------------------------------------------------
  ! Read geometry information section
  !
  ! @param fid : input unit number
  ! @param sname : the name value that caused this function to be called
  ! @param svalue : the value associated with the name (may be mepty string)
  !
  ! @pre sname=fsgeom
  subroutine rdgeom(fid,sname,svalue,istat)
    use strngs
    implicit none
    integer, intent(in) :: fid
    character(len=*), intent(in) :: sname,svalue
    integer, intent(out) :: istat
    logical, dimension(9) :: mask_
    character(32) :: nme_
    character(1024) :: val_

    freg=""
    if (dbc) then
       if (sname.ne.fsgeom) stop "Error: incorrect section name"
       if (len_trim(svalue).ne.0) stop "Error: section does not any parameters"
    endif
    mask_=.false.
    y2w9r4: do
       val_ = " "
       call readnv(fid,nme_,val_,istat)
       ! exit on bad read or section 'end'
       if (istat.ne.0) return
       if (nme_.eq.fsend) exit y2w9r4
       ! looking for rlvest,rlcurv,zl1,zl4,rl1,rl4,rl5,zlim
       j4c2w5: select case (nme_)
       case (fsgzl1) j4c2w5
          read(val_,'(D20.13)')q4b3p4(1)
          mask_(1)=.true.
       case (fsgzl4) j4c2w5
          read(val_,'(D20.13)')q4b3p4(2)
          mask_(2)=.true. ! sim cell 1/2 len
       case (fsgrl1) j4c2w5
          read(val_,'(D20.13)')q4b3p4(3)
          mask_(3)=.true.
       case (fsgrl4) j4c2w5
          read(val_,'(D20.13)')q4b3p4(4)
          mask_(4)=.true.
       case (fsgrl5) j4c2w5
          read(val_,'(D20.13)')q4b3p4(5)
          mask_(5)=.true. ! sim cell rad
       case (fsgrlv) j4c2w5
          read(val_,'(D20.13)')q4b3p4(6)
          mask_(6)=.true.
       case (fsgrlc) j4c2w5
          read(val_,'(D20.13)')q4b3p4(7)
          mask_(7)=.true.
       case (fsgzlm) j4c2w5
          read(val_,'(D20.13)')q4b3p4(8)
          mask_(8)=.true.
       case (fsntrg) j4c2w5
          read(val_,*)ntrg
          mask_(2)=.true.
       case (fsdzg) j4c2w5
          read(val_,*)gzwdth
          mask_(9)=.true.
       case (fsname) j4c2w5
          read(val_,*)freg
       case (fsnsrt) j4c2w5
          read(val_,*)useold
       case default j4c2w5
          call d1s8f2("Name "//nme_//" is not valid in geometry section, value ignored")
       end select j4c2w5
    enddo y2w9r4
    w6k6z5: if (.not.all(mask_)) then
       call d1s8f2("Not all required tags were present")
    endif w6k6z5
    if (.not.mask_(5)) then
      ! rl5 not set. Ok if ntrg is set
      mask_(5) = (ntrg.ne.0)
    endif
    ! Input requirements
    !   ** (zlimit <= zl1)
    !   ** (rlvest + rlcurv <= rl4 - rl1)
    !   ** (rlcurv <= zl1 + rlvest)
    if (zlimit().gt.zl(1)) stop "Invalid geometry input: zlimit is greater than zl1"
    if ((rlvest()+rlcurv()).gt.(rl(4)-rl(1))) stop "Invalid geometry input: rlvest + rlcurv greater than rl4 - rl1"
    if (rlcurv().gt.(zl(2))) stop "Invalid geometry input: rlcurv greater than zl1 + rlvest"
    if (ntrg.eq.0.and.(zl(4).le.zl(2).or.rl(5).le.rl(4))) &
         & stop "Invalid geometry input: channel protein extends outside simulation box boundary"
    irdgeo=.true.
  contains

    subroutine d1s8f2(msg)
      implicit none
      character(len=*), intent(in) :: msg
      write(unit=fidlog,fmt=*)"Bad geometry section in input."
      write(unit=fidlog,fmt=*)msg
      write(unit=fidlog,fmt=*)"Required tags are:"
      write(unit=fidlog,fmt='(1X,A6)')fsgzl1,fsgzl4,fsgrl1,fsgrl4,fsgrl5,fsgrlv,fsgrlc,fsgzlm,fsntrg,fsdzg
      write(unit=fidlog,fmt=*)"Optional tags: ",fsname," ",fsnsrt
      stop 1
    end subroutine d1s8f2

  end subroutine rdgeom

  ! --------------------------
  ! Finalise the read process for geometry module.
  subroutine rfgeom(ntarget_override)
    use spec
    use strngs
    implicit none
    integer, intent(in) :: ntarget_override
    integer :: ispec,igc,ireg ! loop indices
    double precision :: ri ! specie radius
    double precision :: ri_aver ! average free specie radius
    integer :: ri_free ! number of free specie radii
    double precision :: blksbv ! total volume - channel volume
    integer :: iv ! cation valency
    character(14) :: fltout ! buffer for writing formatted float

    ! -- Modify geometry parameters if requested --
    if (dbc) then
       if (.not.irdgeo) stop "Error: geometry initialisation called before defining geometry"
    endif

    write(unit=fidlog,fmt='(72("-"))')
    write(unit=fidlog,fmt=*)"Geometry data summary"
    write(unit=fidlog,fmt='(72("-"))')
    ! -- Need the average ionic specie diameter to expand
    !    the total volume by.
    ri_aver = 0
    ri_free = 0
    do ispec=1,nspec()
       if (isfree(ispec)) then
          ri_aver = ri_aver + xri(ispec)
          ri_free = ri_free + 1
       endif
    enddo
    if (ri_free.gt.0) ri_aver = ri_aver / dble(ri_free)

    ! Over ride ntrg if the overide > 0.
    if (ntarget_override.gt.0) ntrg=ntarget_override

    ! -- Compute simulation geometry parameters if required
    q0m9x6: if (ntrg.gt.0) then
       ! volume of bulk sim is volume from particle number over 
       ! concentration
       s9d1h4=dble(ntrg)*tosi/ionstr()
       ! subtract inner cylinder and vestibules
       blksbv=s9d1h4-(2*pi*sqr(rl(1))*zl(1))-volrot(rl(2),rlvest())


       ! for simulation with a channel we want sim box to be 
       ! as long as diameter (ie. zl4 - zl2 = rl5)
       ! (Here we neglect the channel volume when calculating
       ! channel simulation parameters from the bulk sim volume.)
       ! (Here we increase volume by the average free specie
       ! radius)
       q4b3p4(5)=((blksbv/(pi*2.D0))**(1.D0/3.D0)) + ri_aver
       q4b3p4(2)=zl(2)+q4b3p4(5)
       if (rl(5).le.rl(4)) then
          write(fidlog,*)"Error: Channel protein outer radius (",rl(4)&
               &,") extends outside the simulation box radius (",rl(5),")"
          write(fidlog,*)"calculated from salt concentrations and target particle number."
          write(fidlog,*)"SUGGESTED SOLUTION: increase "//fsntrg//" option in 'geom' input."
          stop "Computed simulation system geometry too small for protein."
       endif
    else q0m9x6
       ! calculate volume of bulk and then ntrg
       ! (add inner cylinder and vestibules)
       s9d1h4=sqr(rl(5))*pi*(zl(4)-zl(2))*2+(2*pi*sqr(rl(1))*zl(1))+volrot(rl(2),rlvest())
       ntrg=ionstr()*s9d1h4/tosi
       if (ntrg.lt.1) then
          write(fidlog,*)"Error: Channel protein occupies too much volume for requested salt concentrations."
          write(fidlog,*)"SUGGESTED SOLUTION: increase "//fsgrl5//" and/or "//fsgzl4//" option in 'geom' input."
          stop "The computed simulation particle number is less than one."
       endif
    endif q0m9x6
    ! for bulk sim we have a cube
    x8h5t7=s9d1h4**(1.D0/3.D0)

    ! ---- dimensions of region for sampling bulk conc. ------
    ! We set the bulk sampling region to be ~1/3 distance from
    ! the boundaries
    zbulk1 = zl(2) + 0.3D0*(zl(4)-zl(2))
    zbulk2 = zl(4) - 0.3D0*(zl(4)-zl(2))
    rbulk  = 2*rl(5)/3
    vbulk  = rbulk**2*pi*(zbulk2-zbulk1)*2

    ! -- Definition of trial insertion regions for each specie --
    w8u1y0: do ispec=1,nspec()
       ri=xri(ispec)
       if (useold) then
          ! Filter region insertion
          if (.not.isfree(ispec)) then
             ! Structural ions have all of particle in filter
             zreg(izlim,ispec)= zlimit()-ri
          else
             zreg(izlim,ispec)= zlimit()
          endif
          ! Centre region insertion
          zreg(ifilt,ispec)= zl(1)
          ! Total channel insertion
          zreg(ichan,ispec)= zl(2)
       else
          ! Filter region insertion (all of particle in filter)
          zreg(izlim,ispec)= zlimit()-ri
          ! Centre region insertion (any of particle in cylinder)
          zreg(ifilt,ispec)= zl(1)+ri
          ! Total channel insertion region (any part of particle in channel)
          zreg(ichan,ispec)= zl(2)+2*ri
       endif

       rreg(izlim,ispec) = rl(1)-ri 
       ! Centre region insertion
       rreg(ifilt,ispec) = rl(1)-ri
       ! Total channel insertion region
       rreg(ichan,ispec) = rl(2)

       ! Bulk insertion region (particle in simulation)
       zreg(ibulk,ispec) = zl(4)-ri
       rreg(ibulk,ispec) = rl(5)-ri

       w5e3o7: do ireg=izlim,ibulk
          vreg(ireg,ispec)=sqr(rreg(ireg,ispec))*pi*(2*zreg(ireg,ispec))
       enddo w5e3o7

       zljin(ispec)=zreg(ichan,ispec)
       vin(ispec)=vreg(ichan,ispec)
       zljout(ispec)=2*(zreg(ibulk,ispec) - zreg(ichan,ispec))
       vout(ispec)=sqr(rreg(ibulk,ispec))*pi*zljout(ispec)
       vtotal(ispec)=vout(ispec) + &
            & ( 2 * pi * sqr(rreg(ifilt,ispec)) * zreg(ifilt,ispec) ) + &
            & ( volrot (rl(2), rlvest() + ri) )

    enddo w8u1y0

    ! -- Wrote the regions for insertion and deletion to logfile
    write(unit=fidlog,fmt=*)" Regions for specie and salt insertion/deletion"
    write(unit=fidlog,fmt='(72("-"))')
    write(unit=fidlog,fmt='(5(1X,A10),1X,A14)')"Spec./salt","Region","Left Z","Right Z","Radius","Vol. Fact."
    write(unit=fidlog,fmt='(72("-"))')
    c2z2f1: do ispec=idxcl(),nspec()
       i0t7t7: do ireg=izlim,ibulk
         call str(vreg(ireg,ispec), fltout)
         write(unit=fidlog,fmt='(2(1X,A10),3(1X,F10.2),1X,A14)')fspc(ispec),freg(ireg),-zreg(ireg,ispec)&
               ,zreg(ireg,ispec),rreg(ireg,ispec),trim(fltout)
       enddo i0t7t7
       call str(vtotal(ispec), fltout)
       write(unit=fidlog,fmt='(2(1X,A10),3(11X),1X,A14)')fspc(ispec),"total",trim(fltout)
    enddo c2z2f1
    f0w4v3: do igc=1,nsalt()
       ispec=isalt(igc)
       iv=nint(xz(ispec))
       m1e3p1: do ireg=izlim,ibulk
          write(unit=fidlog,fmt='(2(1X,A10),3(1X,F10.2),1X,G14.4)')fsalt(igc),freg(ireg),-zreg(ireg,ispec)&
               ,zreg(ireg,ispec),rreg(ireg,ispec),vreg(ireg,ispec)*(vreg(ibulk,idxcl())**iv)
       enddo m1e3p1
    enddo f0w4v3
    write(unit=fidlog,fmt='(72("-"))')
    write(unit=fidlog,fmt=*)"Interpreted geometry input"
    write(unit=fidlog,fmt='(72("-"))')
    call ecgeom(fidlog)
    write(unit=fidlog,fmt='(72("-"))')

    ! initialise gz geometry
    call gz_def()

  contains

    ! ------------------------------------------------------------
    ! DEFINE Z-AXIAL HISTOGRAM BIN GEOMETRY
    !
    subroutine gz_def()
      use spec
      use strngs
      implicit none
      integer :: ispec ! specie loop counters
      integer :: ibin ! current bin number
      double precision :: ri ! specie radius
      double precision :: zup,zlow,vol      ! histogram bin dimension
      double precision, dimension(nspcmx) :: slcmax   ! Maximum zup value for slice in vestibule  
      gz_max=0
      do ispec=1,nspec()
         gzlow_=max(gzlow_,zreg(ibulk,ispec))
         slcmax(ispec)=zl(1)+rlvest()+xri(ispec) 
      enddo
      nhbin = max(1,ceiling(gzlow_/gzwdth))
      gz_max = nhbin * 2
      gzlow_ = -nhbin * gzwdth
      q3r5z7: if (gz_max.gt.nzgmx) then
         write(unit=fidlog,fmt=*)"Attempt to use more z-axial histograms than program limit"
         write(unit=fidlog,fmt='("Attempt = ",i6,"   Program limit = ",i6)')gz_max,nzgmx
         write(unit=fidlog,fmt=*)'SUGGESTED SOLUTION: increase '''//fsdzg//''' option in input file'
         stop 'Program array size limit exceeded'
      endif q3r5z7

      if (.not.allocated(gz_mid)) then
         ! 'gz' distribution attributes
         allocate(gz_mid(gz_max))
         allocate(gz_vol(gz_max,nspec()))
      endif
      gz_mid=0
      gz_vol=0

      ! ------- Geometry for gz distributions ------------------
      write(unit=fidlog,fmt=*)"Parameters for 1D histogram of Z-axis"
      write(unit=fidlog,fmt='(1X,A12,1X,F6.3)')"BIN WIDTH",gzwdth
      write(unit=fidlog,fmt='(1X,A12,1X,I6)')  "BIN COUNT",gz_max
      write(unit=fidlog,fmt='(72("-"))')

      ! This code assumes symmetric geometry around z = 0 and 
      ! calculates volumes starting from z = 0 to z = zmax and
      ! assigns the arrays symmetrically from the centre
      ! ie for volumes: gzvol[nhbin - ibin + 1] = gzvol[nhbin + ibin]
      do ibin=1,nhbin
         ! get bin bounds (as positive values)
         !
         !  i=1, zlow=abs(gzlow)-wid
         !  i : zlow = abs(mid(i)) - wd/2
         !    : zup = zlow + wd
         if (ibin.eq.1) then
            zlow = 0.D0
         else
            zlow = dble(ibin) * gzwdth 
         endif
 
         gz_mid(nhbin + ibin) = zlow + (gzwdth/2)
         gz_mid(nhbin - ibin + 1) = -gz_mid(nhbin + ibin)
         ! use ibin + 1 for one based counting

        if (ibin.eq.nhbin) then
            zup = abs(gzlow_)
         else
            zup = zlow + gzwdth
         endif
         ! calculate volumes
         do ispec = 1,nspec()
            if (zup.lt.zl(1)) then
               ! all in central cylinder
               vol = pi * sqr(rreg(ifilt,ispec)) * gzwdth
            else if (dfeq(zup,zl(1))) then
               ! calculate as all in central cylinder
               vol = pi * sqr(rreg(ifilt,ispec)) * gzwdth
            else if (zlow.gt.slcmax(ispec)) then
               ! all in bulk
               vol = pi * sqr(rreg(ibulk,ispec)) * gzwdth
            else if (dfeq(zlow,slcmax(ispec))) then
               ! calculate as all in bulk
               vol = pi * sqr(rreg(ibulk,ispec)) * gzwdth
            else
               ! somewhere around vestibule.
               ! Check if end-points are outside vestibule
               ri = xri(ispec)
               if (zlow.lt.zl(1)) then
                  ! part inner cylinder, part vol rotate
                  vol = (zl(1) - zlow)*pi*sqr(rreg(ifilt,ispec)) &
                       + volrot(zl(1), rl(2), rlvest() + ri, zl(1), zup)
               else if (zup.gt.slcmax(ispec)) then
                  ! part outer cylinder, part vol rotate
                  vol = (zup - slcmax(ispec))*pi*sqr(rreg(ibulk,ispec)) &
                      + volrot(zl(1), rl(2), rlvest() + ri, zlow, slcmax(ispec))
               else
                  ! all part vol rotate
                  vol = volrot(zl(1), rl(2), rlvest() + ri, zlow, zup)
               endif
            endif
            gz_vol(nhbin + ibin,ispec) = vol 
            gz_vol(nhbin - ibin + 1,ispec) = vol
         enddo
      enddo
    end subroutine gz_def

  end subroutine rfgeom

  ! ------------------------------------------------------------------
  ! Write out geometry input parameters as per input file.
  !
  subroutine ecgeom(fid)
    use strngs
    implicit none
    integer, intent(in) :: fid
    character(20) :: fltout

    write(unit=fid,fmt='(A)')fsgeom
    call str(zl(1), fltout)
    write(unit=fid,fmt='(A,1X,A)')fsgzl1,trim(adjustl(fltout))
    call str(zl(4), fltout)
    write(unit=fid,fmt='(A,1X,A)')fsgzl4,trim(adjustl(fltout))
    call str(rl(1), fltout)
    write(unit=fid,fmt='(A,1X,A)')fsgrl1,trim(adjustl(fltout))
    call str(rl(4), fltout)
    write(unit=fid,fmt='(A,1X,A)')fsgrl4,trim(adjustl(fltout))
    call str(rl(5), fltout)
    write(unit=fid,fmt='(A,1X,A)')fsgrl5,trim(adjustl(fltout))
    call str(rlvest(), fltout)
    write(unit=fid,fmt='(A,1X,A)')fsgrlv,trim(adjustl(fltout))
    call str(rlcurv(), fltout)
    write(unit=fid,fmt='(A,1X,A)')fsgrlc,trim(adjustl(fltout))
    call str(zlimit(), fltout)
    write(unit=fid,fmt='(A,1X,A)')fsgzlm,trim(adjustl(fltout))
    call str(gzwdth, fltout)
    write(unit=fid,fmt='(A,1X,A)')fsdzg,trim(adjustl(fltout))
    write(unit=fid,fmt='(A,1X,I4)')fsntrg,ntrg
    write(unit=fid,fmt='(A,4(1X,"""",A,""""))')fsname,freg(izlim:ibulk)
    call str(useold,fltout)
    write(unit=fid,fmt='(A,1X,A)')fsnsrt,trim(fltout)
    write(unit=fid,fmt='(A)')fsend
    write(unit=fid,fmt=*)

  end subroutine ecgeom


  ! ------------------------------------------------------------------
  ! Define the protein donut
  !
  ! Use the input geometry to define the information that the
  ! 'patch' module will use to generate a set of 'patches' that
  ! cover the protein surface.

  !  (e_w) (-zl2,tpi/2)(U+1)(zl2,tpi/2) (e_w)
  !  (U-1)       _-|-----L2------|-_    (U-1)
  !            _/  |   (rl4)     |  \_
  ! (-zl2,rl3)/__A3.rlc       rlc.A2__\(zl2,rl3,t0)
  ! (tpi,e_w) |                       | (e_w)
  !   (U-1)   W2(L4)   (e_pr)         W1(L3)  (U+1)
  ! (-zl2,rl2)|__A4.rlv       rlv.A1__| (zl2,rl2,t2pi)
  !   (tpi)    \_  |             |  _/
  ! (U-1,e_ch)   \_|_____L1______|_/ (U-1,e_ch)
  ! (-zl2,rl1,t3pi/2) (U-1,e_ch)  (zl2,rl1,t3pi/2)
  !
  subroutine defgeo(donut)
    implicit none
    type (geomdf), intent(inout) :: donut

    ! -----  Define doughnut geometry -----------------
    donut%narch=4
    donut%nline=2
    donut%nwall=2

    ! archs
    donut%za0(1)=  zl(1)
    donut%za0(2)=  zl(3)
    donut%za0(3)= -zl(3)
    donut%za0(4)= -zl(1)

    donut%ra0(1)=  rl(2)
    donut%ra0(2)=  rl(3)
    donut%ra0(3)=  rl(3)
    donut%ra0(4)=  rl(2)

    donut%ra(1)= rlvest()
    donut%ra(2)= rlcurv()
    donut%ra(3)= rlcurv()
    donut%ra(4)= rlvest()

    donut%ta1(1)= 3*pi/2
    donut%ta1(2)= 0
    donut%ta1(3)= pi/2
    donut%ta1(4)= pi

    donut%ta2(1)= 2*pi
    donut%ta2(2)= pi/2
    donut%ta2(3)= pi
    donut%ta2(4)= 3*pi/2

    donut%uasign(1)=-1
    donut%uasign(2)=-1
    donut%uasign(3)=-1
    donut%uasign(4)=-1

    ! lines
    donut%zl1(1)= -zl(1)
    donut%rl1(1)=  rl(1)
    donut%zl2(1)=  zl(1)
    donut%rl2(1)=  rl(1)
    donut%alfa(1)   = 0
    donut%tgalfa(1) = 0

    donut%zl1(2)= -zl(3)
    donut%rl1(2)=  rl(4)
    donut%zl2(2)=  zl(3)
    donut%rl2(2)=  rl(4)
    donut%alfa(2)   = 0
    donut%tgalfa(2) = 0

    donut%ulsign(1)=-1
    donut%ulsign(2)=1
    donut%ulsign(3)=1
    donut%ulsign(4)=-1

    ! walls
    donut%zl1(3)=  zl(2)
    donut%rl1(3)=  rl(2)
    donut%zl2(3)=  zl(2)
    donut%rl2(3)=  rl(3)

    donut%zl1(4)= -zl(2)
    donut%rl1(4)=  rl(2)
    donut%zl2(4)= -zl(2)
    donut%rl2(4)=  rl(3)

    donut%alfa(5)=pi/2
    donut%alfa(6)=pi/2

    call dumpgdfn(donut)
  end subroutine defgeo

  ! ------------------------------------------------------------------
  ! Initialise the donut geometry definition
  !
  ! Use the input geometry to define the information that the
  ! 'patch' module will use to generate a set of 'patches' that
  ! cover the protein surface.
  subroutine initgdfn(donut)
    implicit none
    type (geomdf), intent(inout) :: donut
    ! -----  zero doughnut geometry -----------------
    donut%dfia = 0
    donut%tta1 = 0
    donut%tta2 = 0
    donut%dfil = 0
    donut%zzl1 = 0
    donut%zzl2 = 0
    donut%rrl1 = 0
    donut%rrl2 = 0
    donut%nfia = 0
    donut%nfil = 0
    donut%ilsign = 0
    donut%za0 = 0
    donut%ra0 = 0
    donut%ra = 0
    donut%ta1 = 0
    donut%ta2 = 0
    donut%uasign = 0
    donut%ulsign = 0
    donut%zl1 = 0
    donut%rl1 = 0
    donut%zl2 = 0
    donut%rl2 = 0
    donut%tgalfa = 0
    donut%alfa = 0
    donut%nta = 0
    donut%nzl = 0
    donut%nrl = 0
    donut%narch=0
    donut%nline=0
    donut%nwall=0
    donut%nwall2=0
    donut%nunif=0
    donut%ncent=0
  end subroutine initgdfn

   ! ------------------------------------------------------------------
  ! Dump contents of the protein donut definition
  !
  ! Use the input geometry to define the information that the
  ! 'patch' module will use to generate a set of 'patches' that
  ! cover the protein surface.
  subroutine dumpgdfn(donut)
    implicit none
    type (geomdf), intent(inout) :: donut
    integer :: idx

    ! -----  Output doughnut geometry -----------------
    write(*,*)"ARCS ",donut%narch
    write(*,*)"LINE ",donut%nline
    write(*,*)"WALL ",donut%nwall
    do idx = 1,donut%narch
      write(*,*)"ARC : ",idx, donut%za0(idx), donut%ra0(idx), donut%ra(idx), &
           donut%ta1(idx), donut%ta2(idx), donut%uasign(idx)
    enddo
    if (donut%nline.ne.0) then
      do idx = 1,donut%nline
        write(*,*)"LINE : ",idx, donut%zl1(idx), donut%rl1(idx), donut%zl2(idx), &
             donut%rl2(idx), donut%alfa(idx), donut%tgalfa(idx), donut%ulsign(idx)
      enddo
    endif
    if (donut%nwall.ne.0) then
      do idx = 1+donut%nline,donut%nline+donut%nwall
        write(*,*)"WALL : ",idx, donut%zl1(idx), donut%rl1(idx), donut%zl2(idx), &
             donut%rl2(idx), donut%alfa(idx)
      enddo
    endif
  end subroutine dumpgdfn

  ! ------------------------------------------------------------
  ! Calculate the volume of rotation inside of 
  ! toroid
  !
  ! This performs the analytical integration of:
  !
  !  pi * \int^x1_x0  [ f(x) ]^2 dx
  !  
  !  where f(x) is  (y - yc)^2 = r^2 - (x - xc)^2
  !
  !   => pi * {(yc^2 - r^2)x + (x^3)/3 - yc(x\sqrt(r^2-x^2) + r^2\sin(x/r))}
  !    for x = -r to +r
  !   => pi * {2.yc^2.r - 2.r^3 + 2/3.r^3 - 2.yc.r^2.pi/4}
  !   => 2 * pi * {yc^2.r - yc.r^2.pi/4 - 2/3.r^3}
  !
  pure double precision function vrot2(arc_r, arc_rd)
    implicit none
    double precision, intent(in) :: arc_r, arc_rd

    vrot2 = 2.D0 * pi * arc_rd * ( sqr(arc_r) &
         - arc_r * arc_rd * pi / 4.D0 &
         + 2.D0 * sqr(arc_rd) / 3.D0 &
         )
  end function vrot2

  ! ------------------------------------------------------------
  ! Calculate the volume of rotation inside of toroid slice
  !
  ! \param arc_z Z coordinate of arc center-point
  ! \param arc_r R coordinate of arc center-point
  ! \param arc_rd Radius of arc
  ! \param z0 Z coordinate of start of slice
  ! \param z1 Z coordinate of end of slice (z1 > z0)
  ! 
  ! This performs the analytical integration of:
  !
  !  pi * \int_z0^z1 [ f(z) ]^2 dz
  !  
  !  where f(z) (in z,r plane) is (r - arc_r)^2 = arc_rd^2 - (z - arc_z)^2
  !
  !  => pi * {(arc_r^2 + arc_rd^2)z' - (z'^3)/3 - arc_r(z'\sqrt(arc_rd^2-z'^2)
  !                   + arc_rd^2\sin(z'/arc_rd))}
  !    where z' is (z - arc_z) and is defined for z' in [0:arc_rd]
  !
  ! pure 
  double precision function vrot5(arc_z, arc_r, arc_rd, z_0, z_1)
    implicit none
    double precision, intent(in) :: z_1, z_0, arc_z, arc_r, arc_rd
    double precision :: arc_r2, z1bar_, z0bar_
    if (dbc) then
      if (dfeq(z_0,z_1)) stop "Z0 can not equal Z1 in call to volrot"
      if (z_0.ge.z_1) stop "Z0 must be less than Z1 in call to volrot"
    endif
    arc_r2 = sqr(arc_rd)
    z1bar_ = z_1 - arc_z
    z0bar_ = z_0 - arc_z
    if (dfeq(z0bar_, 0.D0)) then
       ! Handle case where x0 = xc --> z0bar_ = 0
       vrot5 = pi * ( &
            (sqr(arc_r) + arc_r2) * (z_1 - z_0) &
            - ((z1bar_)**3)/3.0D0 &
            - arc_r * (z1bar_ * sqrt(arc_r2 - sqr(z1bar_)) + arc_r2 * asin(z1bar_/arc_rd)) &
            )
    else if (dfeq(z1bar_, arc_rd)) then
       ! handle case where x1 = xc + r --> z1bar_ = r)
       vrot5 = pi * ( &
            (sqr(arc_r) + arc_r2) * (z_1 - z_0) &
            + ((z0bar_)**3 - (z1bar_)**3)/3.0D0 &
            - arc_r * (arc_r2 * pi / 4.D0)             &
            + arc_r * (z0bar_ * sqrt(arc_r2 - sqr(z0bar_)) + arc_r2 * asin(z0bar_/arc_rd)) &
            )
    else
       vrot5 = pi * ( &
            (sqr(arc_r) + arc_r2) * (z_1 - z_0) &
            + ((z0bar_)**3 - (z1bar_)**3)/3.0D0 &
            - arc_r * (z1bar_ * sqrt(arc_r2 - sqr(z1bar_)) + arc_r2 * asin(z1bar_/arc_rd)) &
            + arc_r * (z0bar_ * sqrt(arc_r2 - sqr(z0bar_)) + arc_r2 * asin(z0bar_/arc_rd)) &
            )
    endif
  end function vrot5

  ! --------------------------------------------------------------
  ! CHECK FOR PARTICLE - PROTEIN/MEMBRANE OVERLAP
  ! 
  ! This method checks for particle overlap with any of the
  ! non-particle surfaces.  It is the general method for checking
  ! if a particle is in a valid position with respect to the
  ! geometry.
  subroutine wall(ispec,rznew,r2inew,ovrlap)
    use spec
    implicit none
    integer, intent(in) :: ispec
    double precision, intent(in) :: rznew,r2inew
    logical, intent(out) :: ovrlap
    call ingeom(xri(ispec),rznew,r2inew,(mobile(ispec).or.chonly(ispec)),ovrlap)
  end subroutine wall

end module geom

