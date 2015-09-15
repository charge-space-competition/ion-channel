
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


!  -------------------------------------------------------------
!  ICC data
!
!  @begin docfileformat
!   input file sections: Patch
!
!    patch
!    dxw [REAL] # approximate patch delta width
!    dxf [REAL] # approximate patch delta angle
!    nsub [INT] # integration grid array size parameter
!    epsw
!    epspr
!    end
!
!  @end docfileformat
!
!  baksub
!     convert 'c' vector to 'h' vector by back-substituting
!     from the decomposed 'A' matrix
!  matrix
!     build the initial 'A' matrix and decompose it.
!  rdptch
!     read patch section of input file
!  rfptch
!     initialise the patch system after reading the input file
!
!  Subroutine contained methods of note for patch module maintainers
!
!  readam/writam
!     read or write the inverted 'A' matrix
!  intwll,intarc,intlin
!     surface integration for creation of 'A' matrix
!
module patch
use const
implicit none
private
  ! filename
  character(*), parameter :: fnamx='dat/amx.' 
  ! Particle-patch distance matrix (rip)
  ! Patch area / Particle-patch distance matrix (iprip = parea/rip)
  double precision, public, dimension(:,:), allocatable :: rip,iprip
  ! 'A' matrix
  double precision, private, dimension(:,:), allocatable :: amx
  ! the allocated size of the main patch arrays (<= npchmx)
  integer, private :: npchsz
  ! Patch x,y,z,area,ux,uy,uz,deps date [8 dbls] coords
  ! X coordinate of a patch
  double precision, public, dimension(:), allocatable :: prx
  ! Y coordinate of a patch
  double precision, public, dimension(:), allocatable :: pry
  ! Z coordinate of a patch
  double precision, public, dimension(:), allocatable :: prz
  ! Surface area of a patch
  double precision, public, dimension(:), allocatable :: parea
  ! X dimension of normal vector to centre of patch 
  double precision, public, dimension(:), allocatable :: pux
  ! Y dimension of normal vector to centre of patch
  double precision, public, dimension(:), allocatable :: puy
  ! Z dimension of normal vector to centre of patch
  double precision, public, dimension(:), allocatable :: puz
  ! The effective dielectric constant on the outside of a patch
  double precision, public, dimension(:), allocatable :: deps
  ! H vector
  double precision, public, dimension(:), allocatable :: h
  ! c vector
  double precision, public, dimension(:), allocatable :: c
  ! back substitution index vector
  integer, public, dimension(:), allocatable :: indx
  ! water relative permittivity
  double precision, public :: epsw=80.0  !!INPUT!!
  ! protein relative permittivity
  double precision, public :: epspr=10.0  !!INPUT!!
  ! patch size factors
  double precision, private :: dxf=1.6   !!INPUT!!
  double precision, private :: dxw=1.6   !!INPUT!!
  ! number of patches
  integer, public :: npatch = 0
  ! patch integration factor
  integer, private :: nsub0=10   !!INPUT!!
  ! Has the inverted 'A' matrix been read or created?
  logical, private :: irdamx=.false.
  ! is the calculation using homogeneous permittivity?
  logical, private ::   homog=.false.

  public :: baksub, genrch, rdptch, rfptch, ecptch, patchsaves, is_homogeneous

contains
  ! --------------------------------
  ! Is the ICC being used?
  logical function is_homogeneous()
    use simstate
    implicit none
    is_homogeneous = homog
  end function is_homogeneous

  ! --------------------------------
  ! calculate a gaussbox over ALL surfaces, so parea must be composed
  ! of closed surfaces
  !

  double precision function gaussbox()
    implicit none

    gaussbox = sum(parea(1:npatch)*h(1:npatch)/deps(1:npatch))/qstar()
  end function gaussbox

  double precision function totalarea()
    implicit none
    
    totalarea = sum(parea(1:npatch))
  end function totalarea
    

  ! ------------------------
  ! Output gaussbox charge and area to fidooo
  !
  !
  ! astep : step number 1...
  !
  ! rchg = cumulative gauss sum/e0
  ! archg = cumulative abs(gauss)
  ! r2chg = cumulative gauss^2
  ! area = surface area for which gauss sum is calculated
  !
  ! chg = current gauss sum
  ! echg = avg rchg
  ! aechg = avg archg
  ! er2chg = avg r2chg
  !
  subroutine patchsaves(fid,astep)
    use strngs
    implicit none
    integer, intent(in) :: fid,astep
    double precision, save :: rchg = 0, archg = 0, r2chg = 0
    double precision, save :: area = 0
    double precision, save :: max_chg = 0
    double precision :: chg, echg, aechg, er2chg, abs_chg
    character(20) :: fltout

    if (area.eq.0) area = totalarea()
    
    chg = gaussbox()

    rchg = rchg + chg
    archg = archg + dabs(chg)
    r2chg = r2chg + chg**2

    echg = rchg/astep
    aechg = archg/astep
    er2chg = r2chg/astep
    
    abs_chg = dabs(chg)
    if (max_chg.lt.abs_chg) max_chg = abs_chg

    write(unit=fid, fmt='(" epspr       = ",f10.5)')  epspr
    call str(chg, fltout)
    write(unit=fid, fmt='(" gauss[a]    = ",A)')trim(adjustl(fltout))
    call str(echg, fltout)
    write(unit=fid, fmt=*)' <gauss>     = ',trim(adjustl(fltout))
    call str(aechg, fltout)
    write(unit=fid, fmt=*)' <|gauss|>   = ',trim(adjustl(fltout))
    call str(er2chg-echg**2, fltout)
    write(unit=fid, fmt=*)' Var(gauss)  = ',trim(adjustl(fltout))
    call str(er2chg-aechg**2, fltout)
    write(unit=fid, fmt=*)' Var(|gauss|)= ',trim(adjustl(fltout))
    call str(area, fltout)
    write(unit=fid, fmt=*)' area        = ',trim(adjustl(fltout))
    call str(max_chg, fltout)
    write(unit=fid, fmt=*)' max(|gauss|)= ', trim(adjustl(fltout))
  end subroutine patchsaves
 
  ! Access functions
  ! ------------------------------------------------------------------
  ! Calculate the induced charge per unit area of a patch using
  ! the ICC protocol
  !
  ! The initialisation phase of the ICC protocol generate a solution
  ! matrix for the set of simultaneous equations representing the
  ! patches.  To generate the induced charges on all the patches
  ! we perform a back substitution on the solution matrix using
  ! 'hmat'.  This process is performed by an external Lapack
  ! routine.
  subroutine baksub(hmat)
  implicit none
  integer :: n1,n2,info,nrhs,nx
  character(*), parameter :: trans="N"
  double precision, dimension(npchmx), intent(inout) :: hmat
  external dgetrs
  if (dbc) then
    if (homog) stop "Error: patch::baksub called when system is homogeneous"
    if (.not.irdamx) stop "Error: patch::baksub called with no inverted A matrix"
  endif
  info=0
  nrhs=1
  n1=npatch
  n2=npatch
  nx=npchsz
  ! Use a LAPACK method to back substitute the rhs (b) 
  ! from the solution from ludcmp.
  call dgetrs (trans, n1, nrhs, amx, nx, indx, hmat, n2, info);
  if (0.ne.info) stop "Matrix back-substitution failed."

  end subroutine baksub

  ! ------------------------------------------------------------------
  ! Generate the patch locations from the base geometry
  !
  ! @pre !homog
  subroutine defgrd(geodfn)
  use geom
  implicit none
  type (geomdf), intent(inout) :: geodfn
#include "require.h"

  ! LOCALS
  integer, dimension(5,npchmx) :: isurf
  double precision, dimension(8,300) :: dtta,ara,tac
  double precision, dimension(6,300) :: drl,arr,rlc,dzl,arl,zlc
  double precision, dimension(8) :: dta
  double precision, dimension(6) :: dll
  integer, dimension(5) :: nsurf
  integer ::  iii,is  ! Loop indices
  double precision, dimension(:,:), allocatable :: p__

  if (dbc) then
    if (homog) stop "Called patch method with homogeneous system"
  endif

  call defgeo(geodfn)

  call dumpgdfn(geodfn)

  write(unit=fidlog,fmt=*)" ICC: Calculating tiles on the protein surface"
  write(unit=fidlog,fmt='(72("-"))')
  ! ------ go over 1st line: filter -----------------

  npatch=0
  iii=0
  is=0
  nsurf=0
  allocate(p__(8,npchmx))

  call gofilt

  npatch=npatch+iii
  write(unit=fidlog,fmt=*)" NUMBER OF TILES IN FILTER = ",iii
  write(unit=fidlog,fmt='(50("-"))')

  call goline

  write(unit=fidlog,fmt=*)" NUMBER OF TILES IN OUTER LINE = ",iii-npatch
  write(unit=fidlog,fmt='(50("-"))')
  npatch=iii

  call goarch

  write(unit=fidlog,fmt=*)" NUMBER OF TILES IN ARCHS = ",iii-npatch
  write(unit=fidlog,fmt='(50("-"))')
  npatch=iii

  call gowall

  write(unit=fidlog,fmt=*)" NUMBER OF TILES IN WALLS = ",iii-npatch
  write(unit=fidlog,fmt='(50("-"))')
  write(unit=fidlog,fmt=*)" NUMBER OF TILES TOTAL    = ",iii
  write(unit=fidlog,fmt='(72("-"))')
  npatch=iii

  !----------------------------------------------------------------------
  ! set npchsz and allocate and copy data to patch main arrays
  npchsz=next64(npatch)

  ! If prx is already defined then this is second call to this method.
  if (.not.allocated(prx)) then
    allocate(prx(npchsz))
    allocate(pry(npchsz))
    allocate(prz(npchsz))
    allocate(parea(npchsz))
    allocate(pux(npchsz))
    allocate(puy(npchsz))
    allocate(puz(npchsz))
    allocate(deps(npchsz))
  endif

  prx(1:npatch)=p__(1,1:npatch)
  pry(1:npatch)=p__(2,1:npatch)
  prz(1:npatch)=p__(3,1:npatch)
  parea(1:npatch)=p__(4,1:npatch)
  pux(1:npatch)=p__(5,1:npatch)
  puy(1:npatch)=p__(6,1:npatch)
  puz(1:npatch)=p__(7,1:npatch)
  deps(1:npatch)=p__(8,1:npatch)

  deallocate(p__)

  contains

  ! For lines close and parallel to the z axis
  !
  subroutine gofilt
  implicit none
  integer ::  i,j,k  ! Loop indices
  double precision :: rc,fik,cir
  i=1    
  write(unit=fidlog,fmt='(1X,I1,A)')i,"th line"
  write(unit=fidlog,fmt='(A20,2(1X,F7.2))')"start (z,r) =",geodfn%zl1(i),geodfn%rl1(i)
  write(unit=fidlog,fmt='(A20,2(1X,F7.2))')"end (z,r) =",geodfn%zl2(i),geodfn%rl2(i)
  if (dbc) then
     if (.not.dfeq(geodfn%rl1(i),geodfn%rl2(i))) then
        stop "Line not parallel to z-axis"
     endif
  endif
  dll(i)=(geodfn%zl2(i)-geodfn%zl1(i))
  write(unit=fidlog,fmt='(A20,1X,F7.2)')"length (z) =",dll(i)

  ! --- loop over xz coordinate ---------------------
  geodfn%nzl(i)=max(10,int(dll(i)/dxf)+1)
  write(unit=fidlog,fmt='(A20,1X,I4)')"mesh size (z) =",geodfn%nzl(i)
  i3y5r5: do j=1,geodfn%nzl(i)
  !    write(unit=fidlog,fmt=*)j,"th ring"
    dzl(i,j)=(geodfn%zl2(i)-geodfn%zl1(i))/dble(geodfn%nzl(i))
    o3j9m0: if (j.eq.1) then
      geodfn%zzl1(i,j)=geodfn%zl1(i)
      geodfn%zzl2(i,j)=geodfn%zl1(i)+dzl(i,j)
    else o3j9m0
      geodfn%zzl1(i,j)=geodfn%zzl2(i,j-1)
      geodfn%zzl2(i,j)=geodfn%zzl1(i,j)+dzl(i,j)
    endif o3j9m0

    arl(i,j)=(geodfn%zzl2(i,j)-geodfn%zzl1(i,j))
    zlc(i,j)=(geodfn%zzl2(i,j)+geodfn%zzl1(i,j))/2

  ! --- loop over phi angle ------------------------
    rc=geodfn%rl1(i)
    cir=2*pi*rc
    geodfn%nfil(i,j)=max(16,int(cir/dxf)+1)
    geodfn%dfil(i,j)=2*pi/dble(geodfn%nfil(i,j))

    write(unit=fidlog,fmt='(A20,1X,I4)')"mesh size (r) =",geodfn%nfil(i,j)

    b3a4l6: do k=1,geodfn%nfil(i,j)
      iii=iii+1
      is=is+1
      nsurf(1)=nsurf(1)+1
      isurf(1,nsurf(1))=is
      fik=(k-0.5D0)*geodfn%dfil(i,j)
      p__(1,iii)=rc*dcos(fik)
      p__(2,iii)=rc*dsin(fik)
      p__(3,iii)=zlc(i,j)
      if (dbc) then
        if (isnan(p__(1,iii))) stop "px is NaN"
        if (isnan(p__(2,iii))) stop "py is NaN"
        if (isnan(p__(3,iii))) stop "pz is NaN"
      endif
      p__(5,iii)=geodfn%ulsign(i)*dcos(fik)
      p__(6,iii)=geodfn%ulsign(i)*dsin(fik)
      p__(7,iii)=-0.D0
      if (dbc) then
        if (isnan(p__(5,iii))) stop "ux is NaN"
        if (isnan(p__(6,iii))) stop "uy is NaN"
        if (isnan(p__(7,iii))) stop "uz is NaN"
      endif
      ! zlength = arl  rlength = 2.pi.r/nfil <==> dfil.r
      p__(4,iii)=arl(i,j)*geodfn%dfil(i,j)*rc
      if (dbc) then
        if (isnan(p__(4,iii))) stop "area is NaN"
        if (isnan(arl(i,j))) stop "arl is NaN"
      endif
      p__(8,iii)=2*(epspr-epsw)/(epspr+epsw)
    enddo b3a4l6
  enddo i3y5r5
  write(unit=fidlog,fmt="(50('-'))")
  end subroutine gofilt

  ! For lines distant and parallel to the z-axis
  subroutine goline
  implicit none
  integer ::  i,j,k  ! Loop indices
  double precision :: rc,dx,fik,cir
  ! ------ go over outer cylinder -------------------

  p6o0d0: do i=2,geodfn%nline
    write(unit=fidlog,fmt='(1X,I1,A)')i,"th line"
    write(unit=fidlog,fmt='(A20,1X,F7.2,1X,F7.2)')"start (z,r) =",geodfn%zl1(i),geodfn%rl1(i)
    write(unit=fidlog,fmt='(A20,1X,F7.2,1X,F7.2)')"end (z,r) =",geodfn%zl2(i),geodfn%rl2(i)
    dll(i)=(geodfn%zl2(i)-geodfn%zl1(i))
    if (dbc) then
       if (.not.dfeq(geodfn%rl1(i),geodfn%rl2(i))) then
          stop "Line not parallel to z-axis"
       endif
    endif
    write(unit=fidlog,fmt='(A20,1X,F7.2)')"length (z) =",dll(i)

  ! --- loop over z coordinate ---------------------
    dx=dxw*dsqrt(geodfn%rl2(i))
    geodfn%nzl(i)=max(4,int(dll(i)/dx)+1)
    write(unit=fidlog,fmt='(A20,1X,I4)')"mesh size (z) =",geodfn%nzl(i)
    p3o5e7: do j=1,geodfn%nzl(i)
  !      write(unit=fidlog,fmt=*)j,"th ring"
      dzl(i,j)=(geodfn%zl2(i)-geodfn%zl1(i))/dble(geodfn%nzl(i))
      e2g3j7: if (j.eq.1) then
        geodfn%zzl1(i,j)=geodfn%zl1(i)
        geodfn%zzl2(i,j)=geodfn%zl1(i)+dzl(i,j)
      else e2g3j7
        geodfn%zzl1(i,j)=geodfn%zzl2(i,j-1)
        geodfn%zzl2(i,j)=geodfn%zzl1(i,j)+dzl(i,j)
      endif e2g3j7
      arl(i,j)=(geodfn%zzl2(i,j)-geodfn%zzl1(i,j))
      zlc(i,j)=(geodfn%zzl2(i,j)+geodfn%zzl1(i,j))/2

  ! --- loop over phi angle ------------------------
      rc=geodfn%rl1(i)
      cir=2*pi*rc
      geodfn%nfil(i,j)=max(16,int(cir/dx)+1)
      geodfn%dfil(i,j)=2*pi/dble(geodfn%nfil(i,j))
      write(unit=fidlog,fmt='(A20,1X,I4)')"mesh size (r) =",geodfn%nfil(i,j)

      v4z2z1: do k=1,geodfn%nfil(i,j)
        iii=iii+1
        is=is+1
        nsurf(2)=nsurf(2)+1
        isurf(2,nsurf(2))=is
        fik=(k-0.5D0)*geodfn%dfil(i,j)
        p__(1,iii)=rc*dcos(fik)
        p__(2,iii)=rc*dsin(fik)
        p__(3,iii)=zlc(i,j)
        if (dbc) then
          if (isnan(p__(1,iii))) stop "px is NaN"
          if (isnan(p__(2,iii))) stop "py is NaN"
          if (isnan(p__(3,iii))) stop "pz is NaN"
        endif
        p__(5,iii)=geodfn%ulsign(i)*dcos(fik)
        p__(6,iii)=geodfn%ulsign(i)*dsin(fik)
        p__(7,iii)=-0.D0
        if (dbc) then
          if (isnan(p__(5,iii))) stop "ux is NaN"
          if (isnan(p__(6,iii))) stop "uy is NaN"
          if (isnan(p__(7,iii))) stop "uz is NaN"
        endif
        ! zlength = arl  rlength = 2.pi.r/nfil <==> dfil.r
        p__(4,iii)=arl(i,j)*geodfn%dfil(i,j)*rc
        if (dbc) then
          if (isnan(p__(4,iii))) stop "area is NaN"
          if (isnan(arl(i,j))) stop "arl is NaN"
        endif
        p__(8,iii)=2*(epspr-epsw)/(epspr+epsw)
      enddo v4z2z1
    enddo p3o5e7
    write(unit=fidlog,fmt="(50('-'))")
  enddo p6o0d0
  end subroutine goline

  !
  ! Jim Fonseca's method for calculating patches on an arch
  !
  subroutine goarch
  implicit none
  double precision :: dla
  integer ::  i,j,k  ! Loop indices
  double precision :: dx,rc,fik,szaml,cir
  c1q4v0: do i=1,geodfn%narch
  
    write(unit=fidlog,fmt='(1X,I1,A)')i,"th arch"
    write(unit=fidlog,fmt='(A20,1X,F7.2,1X,F7.2)')"center (z,r) =",geodfn%za0(i),geodfn%ra0(i)
    write(unit=fidlog,fmt='(A20,1X,F7.2)')"radius =",geodfn%ra(i)
    write(unit=fidlog,fmt='(A20,1X,F7.2,1X,F7.2)')"angles (begin,end) =",geodfn%ta1(i),geodfn%ta2(i)
    dta(i)=geodfn%ta2(i)-geodfn%ta1(i)
    write(unit=fidlog,fmt='(A20,1X,F7.2)')"arc =",dta(i)   
    dla=geodfn%ra(i)*dta(i)
    write(unit=fidlog,fmt='(A20,1X,F7.2)')"length (z) =",dla
    dx=dxw*dsqrt(geodfn%ra0(i))

  ! loop over theta angle 
  !theta is measured along channel axis
    geodfn%nta(i) = max(4,int(dla/dx)+1)
        
    write(unit=fidlog,fmt='(A20,1X,I4)')"mesh size (z=arc) =",geodfn%nta(i)
        
  !loop over arc        
    do j=1,geodfn%nta(i)
  !    write(unit=fidlog,fmt=*)j,"th ring"
  !each section of the arc has length dtta(arcnumber,ring)
      dtta(i,j)=dta(i)/dble(geodfn%nta(i))
      if (j.eq.1) then
  ! if it's the first ring, set to beginning theta of the arc keep
  ! in mind it seems that we are counting arcs from right to left,
  ! since we start at theta1, hopefully this doesn't make much of a
  ! difference (i.e. doesn't matter if go ta1 to ta2 or vice versa)
        geodfn%tta1(i,j)=geodfn%ta1(i)
        geodfn%tta2(i,j)=geodfn%ta1(i)+dtta(i,j)
      else
        geodfn%tta1(i,j)=geodfn%tta2(i,j-1)
        geodfn%tta2(i,j)=geodfn%tta1(i,j)+dtta(i,j)
      endif
  !find the center theta (which is not just the average of tta1,tta2)
      ara(i,j)=geodfn%ra0(i)*dtta(i,j) &
        -geodfn%ra(i)*(dcos(geodfn%tta2(i,j))-dcos(geodfn%tta1(i,j)))
      szaml=0.5*geodfn%ra0(i)*(geodfn%tta2(i,j)**2-geodfn%tta1(i,j)**2) &
        - geodfn%ra(i)*(geodfn%tta2(i,j)*dcos(geodfn%tta2(i,j))         &
               -geodfn%tta1(i,j)*dcos(geodfn%tta1(i,j)))+               &
        geodfn%ra(i)*(dsin(geodfn%tta2(i,j))-dsin(geodfn%tta1(i,j)))
      tac(i,j)=szaml/ara(i,j)
      rc=geodfn%ra0(i)+geodfn%ra(i)*dsin(tac(i,j))
  !circumference of this ring            
      cir=2*pi*rc
      geodfn%nfia(i,j)=max(16,int(cir/dx)+1)
      geodfn%dfia(i,j)=2*pi/dble(geodfn%nfia(i,j))
      write(unit=fidlog,fmt='(A20,3(1X,F7.2))')"theta (lo,mid,hi) =",geodfn%tta1(i,j),tac(i,j),geodfn%tta2(i,j)
      write(unit=fidlog,fmt='(A20,1X,I4)')"mesh size (r) =",geodfn%nfia(i,j)

  !loop over all the slices of phi            
      z9c4a8: do k=1,geodfn%nfia(i,j)
        iii=iii+1
        is=is+1
        y4h5p2: if (i.eq.1.or.i.eq.4) then
          nsurf(1)=nsurf(1)+1
          isurf(1,nsurf(1))=is
        else y4h5p2
          nsurf(2)=nsurf(2)+1
          isurf(2,nsurf(2))=is
        endif y4h5p2
        fik=(k-0.05D0)*geodfn%dfia(i,j)
        if (dbc) then
          if (isnan(fik)) stop "fik is NaN"
        endif
        p__(1,iii)=rc*dcos(fik)
        p__(2,iii)=rc*dsin(fik)
        p__(3,iii)=geodfn%za0(i)+geodfn%ra(i)*dcos(tac(i,j))
        if (dbc) then
          if (isnan(p__(1,iii))) stop "px is NaN"
          if (isnan(p__(2,iii))) stop "py is NaN"
          if (isnan(p__(3,iii))) stop "pz is NaN"
        endif
        p__(5,iii)=-geodfn%uasign(i)*dsin(tac(i,j))*dcos(fik)
        p__(6,iii)=-geodfn%uasign(i)*dsin(tac(i,j))*dsin(fik)
        p__(7,iii)=-geodfn%uasign(i)*dcos(tac(i,j))
        if (dbc) then
          if (isnan(p__(5,iii))) stop "ux is NaN"
          if (isnan(p__(6,iii))) stop "uy is NaN"
          if (isnan(p__(7,iii))) stop "uz is NaN"
        endif
        p__(4,iii)=geodfn%ra(i)*ara(i,j)*geodfn%dfia(i,j)
        if (dbc) then
          if (isnan(p__(4,iii))) stop "area is NaN"
          if (isnan(ara(i,j))) stop "ara is NaN"
        endif
        p__(8,iii)=2*(epspr-epsw)/(epspr+epsw)
      enddo z9c4a8
    enddo
    write(unit=fidlog,fmt="(50('-'))")
  enddo c1q4v0
  ! end of loop over archs 
  end subroutine goarch

  subroutine gowall
  implicit none
  integer ::  i,j,k  ! Loop indices
  double precision :: cir,ratio,rc,dfiu,fik,szaml
  integer :: nfiu
  
  ! ------ go over last walls of protein --------

  b0j3s2: do i=geodfn%nline+1,geodfn%nline+geodfn%nwall

    write(unit=fidlog,fmt='(1X,I1,A)')i,"th line/wall"
    write(unit=fidlog,fmt='(A20,2(1X,F7.2))')"start (z,r) =",geodfn%zl1(i),geodfn%rl1(i)
    write(unit=fidlog,fmt='(A20,2(1X,F7.2))')"end (z,r) =",geodfn%zl2(i),geodfn%rl2(i)
    dll(i)=geodfn%rl2(i)-geodfn%rl1(i)
    write(unit=fidlog,fmt='(A20,1X,F7.2)')"length (r) =",dll(i)

  ! --- loop over the rest with increasing tiles ---
    cir=2*pi*(geodfn%rl1(i)+dll(i)/2)
    nfiu=int(cir/dxf)+1
    dfiu=2*pi/dble(nfiu)
    ratio=(1+dfiu/2)/(1-dfiu/2)
    ! write(unit=fidlog,fmt=*)" Ratio old = ",ratio
    geodfn%nrl(i)=int(dlog(geodfn%rl2(i)/geodfn%rl1(i))/dlog(ratio))+1
    ! write(unit=fidlog,fmt=*)" N r for increasing tiles = ",geodfn%nrl(i)
    ratio=(geodfn%rl2(i)/geodfn%rl1(i))**(1/dble(geodfn%nrl(i)))
    ! write(unit=fidlog,fmt=*)" Ratio new = ",ratio
    write(unit=fidlog,fmt='(A20,1X,I4)')"mesh size (r)",geodfn%nrl(i)
    geodfn%rrl1(i,1)=geodfn%rl1(i)
    w8w7s2: do j=1,geodfn%nrl(i)
      geodfn%dfil(i,j)=dfiu
      geodfn%rrl2(i,j)=geodfn%rrl1(i,j)*ratio
      ! write(unit=fidlog,fmt=*)j,"th ring"
      drl(i,j)=geodfn%rrl2(i,j)-geodfn%rrl1(i,j)
      ! write(unit=fidlog,fmt=*)"   dr     = ",drl(i,j)
      arr(i,j)=0.5D0*(geodfn%rrl2(i,j)**2-geodfn%rrl1(i,j)**2)
      szaml=(1.D0/3.D0)*(geodfn%rrl2(i,j)**3-geodfn%rrl1(i,j)**3)
      rlc(i,j)=szaml/arr(i,j)

  ! --- loop over phi angle --------------------
      rc=rlc(i,j)
      geodfn%nfil(i,j)=nfiu
      geodfn%dfil(i,j)=dfiu
      write(unit=fidlog,fmt='(A20,3(1X,F7.2))')"r (lo,mid,hi) =",geodfn%rrl1(i,j),rlc(i,j),geodfn%rrl2(i,j)
      write(unit=fidlog,fmt='(A20,1X,I4)')"mesh size (r) =",geodfn%nfil(i,j)
      i6o4l5: do k=1,geodfn%nfil(i,j)
        iii=iii+1
        is=is+1
        nsurf(2)=nsurf(2)+1
        isurf(2,nsurf(2))=is
        fik=(k-0.5D0)*geodfn%dfil(i,j)
        p__(1,iii)=rc*dcos(fik)
        p__(2,iii)=rc*dsin(fik)
        p__(3,iii)=geodfn%zl1(i)
        p__(5,iii)=0
        p__(6,iii)=0
        p__(7,iii)=geodfn%ulsign(i)
        p__(4,iii)=arr(i,j)*geodfn%dfil(i,j)
        p__(8,iii)=2*(epspr-epsw)/(epspr+epsw)
      enddo i6o4l5
      geodfn%rrl1(i,j+1)=geodfn%rrl2(i,j)
    enddo w8w7s2
    write(unit=fidlog,fmt="(50('-'))")
  enddo b0j3s2

  end subroutine gowall

  end subroutine defgrd

  ! --------------------------------------------------------------
  ! Echo the input
  !
  ! Write the parameters that may be set in the input file in
  ! the same format as the input file.  Optional input will
  ! be written here as the default value.
  subroutine ecptch(fid)
  use strngs
  implicit none
  integer, intent(in) :: fid
  character(20) :: fltout
  write(unit=fid,fmt='(A)')fsptch
  call str(dxf, fltout)
  write(unit=fid,fmt='(A,1X,A)')fsdxf,trim(adjustl(fltout))
  call str(dxw, fltout)
  write(unit=fid,fmt='(A,1X,A)')fsdxw,trim(adjustl(fltout))
  write(unit=fid,fmt='(A,1X,I4)')fsnsub,nsub0
  call str(epsw, fltout)
  write(unit=fid,fmt='(A,1X,A)')fsepsw,trim(adjustl(fltout))
  call str(epspr, fltout)
  write(unit=fid,fmt='(A,1X,A)')fsepsp,trim(adjustl(fltout))
  write(unit=fid,fmt='(A)')fsend
  write(unit=fid,fmt=*)
  end subroutine ecptch

  !------------------------------------------------------
  !
  ! Generate the initial h and c vectors
  subroutine genrch
  use conf
  use spec
  use strngs
  implicit none

  ! LOCALS
  double precision rxi,ryi,rzi  ! particle coords
  double precision rxik,ryik,rzik,rzkj ! particle-patch data
  double precision chgt,areat   ! summary data for pacthes
  integer ipch  ! patch index
  integer ii    ! particle index
  integer ispec ! specie index
  character(20) :: fltout
 
  if (.not.irdamx) stop "patch%genrch called before processing 'A' matrix"

  c=0
  u3s0t5: do ii=1,nactv
    ispec=ispcbk(ii)
    if (ispec.eq.0) cycle u3s0t5
    rxi=rx(ii)
    ryi=ry(ii)
    rzi=rz(ii)
    l9c5v5: do ipch=1,npatch
      rxik=prx(ipch)-rxi
      ryik=pry(ipch)-ryi
      rzik=prz(ipch)-rzi
      rip(ipch,ii)=dsqrt(rxik*rxik+ryik*ryik+rzik*rzik)
      iprip(ipch,ii)=parea(ipch)/rip(ipch,ii)
      rzkj=rxik*pux(ipch)+ryik*puy(ipch)+rzik*puz(ipch)
      c(ipch)=c(ipch)-deps(ipch)*xq(ispec)*rzkj/(4*pi*epsw*(rip(ipch,ii)**3))
    enddo l9c5v5
  enddo u3s0t5

  forall (ipch=1:npatch) h(ipch)=c(ipch)

  ! ------- CALCULATE INITIAL H ------------------------

  call baksub(h)
  areat = totalarea()
  chgt = gaussbox()

  write(unit=fidlog,fmt='(72("-"))')
  write(unit=fidlog,fmt=*)"Initial ICC calculation:"
  call str(chgt, fltout)
  write(unit=fidlog,fmt=*) &
       ' Total initial induced charge / e : ', &
       trim(adjustl(fltout))
  call str(areat, fltout)
  write(unit=fidlog,fmt=*) &
       ' Total area / Ang^2 = ', &
       trim(adjustl(fltout))
  write(unit=fidlog,fmt='(72("-"))')
  end subroutine genrch

  ! ------------------------------------------------------------------
  ! Compute the A matrix
  !
  ! NOTE: The saved matrix does not adjust the i=j or multiply by deps

  subroutine matrix
  use geom
  implicit none
  type (geomdf), allocatable :: geodfn

  external dgetrf
  ! LOCALS
  double precision fi1,fi2    ! derived patch size params
  double precision  aij       ! integrator output
  ! LOCALS
  integer i,j,k ! loop indices
  integer nsub     ! patch number (input) parameters
  integer ipatch,jpatch              ! 'A' matrix indices
  integer np10  ! counter for percent-complete output
  character*(*) fmt101
  parameter (fmt101='(2X,I8)')
  ! LOCAL
  integer info,npch2

  if (dbc) then
    if (irdamx) stop "patch::matrix called when inverted matrix exists"
  endif
  ! if not aa then rdmtrx
  allocate(geodfn)

  ! -----  Define geometry ------------------------------
  call defgrd(geodfn)

  ! Npatch defined, allocate memory
  call f4s3s6
  ! ----  Fill matrix -----------------------------------

  write(unit=fidlog,fmt=*)" Filling matrix, patience"
  np10=npatch/10

  e2o3c9: do ipatch=1,npatch

    if (mod(ipatch,np10).eq.0) write(unit=fidlog,fmt='(i3," %")')10*ipatch/np10

    jpatch=0

    h8s9k7: do i=1,geodfn%nline
      x4x9c4: do j=1,geodfn%nzl(i)
        r5y4j7: do k=1,geodfn%nfil(i,j)
          jpatch=jpatch+1
          p8q0w3: if (ipatch.eq.jpatch) then
  !           if this patch is 'seeing' itself, we need lots of tiny patches
            nsub=5*nsub0
          else p8q0w3
            nsub=nsub0
          endif p8q0w3
  !         fi1/2 is start/end of distance (in radians) around the phi angle
          fi1=(k-1)*geodfn%dfil(i,j)
          fi2=k*geodfn%dfil(i,j)
          call intlin(nsub,geodfn%zl1(i),geodfn%rl1(i),geodfn%tgalfa(i),geodfn%zzl1(i,j),geodfn%zzl2(i,j),fi1,fi2,ipatch,aij)
          amx(ipatch,jpatch)=aij/(4*pi)
        enddo r5y4j7
      enddo x4x9c4
    enddo h8s9k7

    q1w3h2: do i=1,geodfn%narch
      h8l6c7: do j=1,geodfn%nta(i)
        f8g1o9: do k=1,geodfn%nfia(i,j)
          jpatch=jpatch+1
          k4v5z9: if (ipatch.eq.jpatch) then
  !           if this patch is 'seeing' itself, we need lots of tiny patches
            nsub=5*nsub0
          else k4v5z9
            nsub=nsub0
          endif k4v5z9
  !         fi1/2 is start/end of distance (in radians) around the phi angle
          fi1=(k-1)*geodfn%dfia(i,j)
          fi2=k*geodfn%dfia(i,j)
  !          write(unit=fidlog,fmt=*)ipatch,jpatch,nsub,geodfn%za0(i),geodfn%ra0(i),geodfn%ra(i),geodfn%tta1(i,j),
  ! :          geodfn%tta2(i,j),fi1,fi2,uasign(i)
          call intarc(nsub,geodfn%za0(i),geodfn%ra0(i),geodfn%ra(i),geodfn%tta1(i,j),geodfn%tta2(i,j),fi1,fi2,ipatch,aij)
          amx(ipatch,jpatch)=aij/(4*pi)
        enddo f8g1o9
      enddo h8l6c7
    enddo q1w3h2

    l3d7u8: do i=geodfn%nline+1,geodfn%nline+geodfn%nwall
      v0y7z9: do j=1,geodfn%nrl(i)
        a9e7o1: do k=1,geodfn%nfil(i,j)
          jpatch=jpatch+1
          x2t2t9: if (ipatch.eq.jpatch) then
  !           if this patch is 'seeing' itself, we need lots of tiny patches
            nsub=5*nsub0
          else x2t2t9
            nsub=nsub0
          endif x2t2t9
  !         fi1/2 is start/end of distance (in radians) around the phi angle
          fi1=(k-1)*geodfn%dfil(i,j)
          fi2=k*geodfn%dfil(i,j)
          call intwll(nsub,geodfn%zl1(i),geodfn%rrl1(i,j),geodfn%rrl2(i,j),fi1,fi2,ipatch,aij)
          amx(ipatch,jpatch)=aij/(4*pi)
        enddo a9e7o1
      enddo v0y7z9
    enddo l3d7u8
  enddo e2o3c9

  deallocate(geodfn)
  call dumpam( "amx0.dat" )

  ! -----------------------------------------------------
  ! Update A matrix with deps etc
  forall (i=1:npatch) amx(i,1:npatch)=deps(i)*amx(i,1:npatch)
  indx=0
  forall (i=1:npatch) amx(i,i)=amx(i,i)+1

  call dumpam( "amx1.dat" )

  ! -----------------------------------------------------
  ! Do LU decomposition of the A matrix
  ! call ludcmp
  write(unit=fidlog,fmt=*)"LU decomposing the 'A' matrix (using BLAS)"
  npch2=npchsz
  info=0
  call dgetrf(npatch, npchsz, amx, npch2, indx, info)
  if (debug) write(unit=fidlog,fmt=*)"! Matrix has been inverted, info = ",info
  if (0.ne.info) stop "Matrix inversion failed, no save."

  irdamx=.true.
  call writam

  contains

  ! -----------------------------------------------------
  ! Integrate a patch on an arc
  subroutine intarc(nsub,z0,r0,r,t1,t2,fi1,fi2,ii,aij)
  implicit none
  integer nsub,ii
  double precision z0,r0,r,t1,t2,fi1,fi2,aij
#include "require.h"

  ! LOCALS
  double precision ar,arel,cosfij,costtc,dfi,dfisub,dt,dtsub
  double precision fij,pxi,pxj,pxij,pyi,pyj,pyij,pzi,pzj,pzij
  double precision rij,rij3,rijsq,sinfij,sinttc,rc,szaml,tt1,tt2
  double precision ttc,uxi,uyi,uzi
  integer i,j ! loop indices

  ar=0
  pxi=prx(ii)
  pyi=pry(ii)
  pzi=prz(ii)
  uxi=pux(ii)
  uyi=puy(ii)
  uzi=puz(ii)

  dt=t2-t1
  dtsub=dt/dble(nsub)
  dfi=fi2-fi1
  dfisub=dfi/dble(nsub)

  aij=0
  if (dfeq(0.0D0,dfi)) stop "error in intarc, dfi = 0"
  ! --- double loop over subtiles ----------------------
  ! --- loop over theta angle --------------------------

  o3f5x0: do i=1,nsub
    r3o7r0: if (dfeq(0.0D0,dtsub)) then
      sinttc=0
      costtc=1
    else   r3o7r0
      tt1=t1+(i-1)*dtsub
      tt2=tt1+dtsub
      ar=r0*dtsub-r*(dcos(tt2)-dcos(tt1))
      szaml=0.5D0*r0*(tt2**2-tt1**2)-r*(tt2*dcos(tt2)-tt1*dcos(tt1))+r*(dsin(tt2)-dsin(tt1))
      ttc=szaml/ar
      if (dbc) then
        if (isnan(ttc)) stop "error in intarc, ttc is NaN"
      endif
      sinttc=dsin(ttc)
      costtc=dcos(ttc)
    endif r3o7r0
    rc=r0+r*sinttc
    arel=r*ar*dfisub
    if (dbc) then
      if (isnan(rc)) stop "error in intarc, rc is NaN"
      if (isnan(arel)) stop "error in intarc, arel is NaN"
    endif
  ! --- loop over phi angle ----------------------------
    pzj=z0+r*costtc

    y6p7q8: do j=1,nsub
      fij=fi1+(j-0.5D0)*dfisub
      cosfij=dcos(fij)
      sinfij=dsin(fij)
      pxj=rc*cosfij
      pyj=rc*sinfij
      pxij=pxi-pxj
      pyij=pyi-pyj
      pzij=pzi-pzj
      rijsq=pxij*pxij+pyij*pyij+pzij*pzij
      rij=dsqrt(rijsq)
      rij3=rijsq*rij
      aij=aij+(pxij*uxi+pyij*uyi+pzij*uzi)*arel/rij3
    enddo y6p7q8
  enddo o3f5x0
  if (dbc) then
    if (isnan(aij)) stop "error in intarc, aij is NaN"
  endif
  end subroutine intarc

  ! -----------------------------------------------------
  ! Integrate a patch on cylindrical surface

  subroutine intlin(nsub,z0,r0,talfa,z1,z2,fi1,fi2,ii,aij)
  implicit none
  integer nsub,ii
  double precision z0,r0,talfa,z1,z2,fi1,fi2,aij
#include "require.h"
  ! LOCALS
  double precision ar,arel,cosfij,dfi,dfisub,dz,dzsub
  double precision fij,pxi,pxj,pxij,pyi,pyj,pyij,pzi,pzj,pzij
  double precision rij,rij3,rijsq,sinfij,rc,szaml,zz1,zz2
  double precision zc,uxi,uyi,uzi,calfa
  integer i,j ! loop indices

  pxi=prx(ii)
  pyi=pry(ii)
  pzi=prz(ii)
  uxi=pux(ii)
  uyi=puy(ii)
  uzi=puz(ii)

  ! calfa is cos of alfa
  calfa=1/dsqrt(1+talfa**2)
  ! z2 and z1 are end and begin of this particular ring
  dz=z2-z1
  ! so split up each ring into pieces (10)
  dzsub=dz/dble(nsub)
  ! also around the phi angle
  dfi=fi2-fi1
  dfisub=dfi/dble(nsub)

  aij=0
  ! gonna sum up aij for each patch
  ! will sum over tiny patches (nsub*nsub tiny patches)

  ! --- double loop over subtiles ------
  ! --- loop over theta angle ----------
  ! so each patch is going to have nsub*nsub tiny patches
  i2a3m1: do i=1,nsub
      zz1=z1+(i-1)*dzsub
      zz2=zz1+dzsub
  ! get the centroid in this weird way i can't figure out
  ! r0 is radius if line at beginning (left) side
  ! z0 beginning (left) of line
  ! talfa is tangent of the slope of the line
  ! dzsub is spacing between each tiny patch
  ! zz2 is right side of tinypatch
  ! zz1 is left side of tinypatch
    ar=(r0-z0*talfa)*dzsub+talfa*0.5D0*(zz2**2-zz1**2)
    szaml=0.5D0*(r0-z0*talfa)*(zz2**2-zz1**2)+talfa*(1.D0/3.D0)*(zz2**3-zz1**3)
    zc=szaml/ar
    rc=r0-z0*talfa+zc*talfa
    arel=rc/calfa

  ! --- loop over phi angle ------------

    r3g2u7: do j=1,nsub
  ! loop over tiny patches around phi angle
  ! fij is center of tiny patch around phi angle
      fij=fi1+(j-0.5D0)*dfisub
      cosfij=dcos(fij)
      sinfij=dsin(fij)
  ! pxj,pyz, and pzj are centroid of tiny patch
      pxj=rc*cosfij
      pyj=rc*sinfij
      pzj=zc
  ! pxi,pyi,pzi is center of patch
      pxij=pxi-pxj
      pyij=pyi-pyj
      pzij=pzi-pzj
      rijsq=pxij*pxij+pyij*pyij+pzij*pzij
  ! rij is the distance from the patch centroid to the tiny patch
      rij=dsqrt(rijsq)
      rij3=rijsq*rij
      aij=aij+(pxij*uxi+pyij*uyi+pzij*uzi)*arel*dzsub*dfisub/rij3
    enddo r3g2u7
  enddo i2a3m1
  if (dbc) then
    if (isnan(aij)) stop "error in intlin, aij is NaN"
  endif
  end subroutine intlin

  ! -----------------------------------------------------
  ! Integrate patch on a disk surface
  subroutine intwll(nsub,z0,r1,r2,fi1,fi2,ii,aij)
  implicit none
  integer nsub,ii
  double precision z0,r1,r2,fi1,fi2,aij
#include "require.h"

  ! LOCALS
  double precision ar,arel,dfi,dfisub,dr,drsub
  double precision fij,pxi,pxj,pxij,pyi,pyj,pyij,pzi,pzj,pzij
  double precision rij,rij3,rijsq,rc,szaml,rr1,rr2,uxi,uyi,uzi
  integer i,j ! loop indices

  pxi=prx(ii)
  pyi=pry(ii)
  pzi=prz(ii)
  uxi=pux(ii)
  uyi=puy(ii)
  uzi=puz(ii)

  dr=r2-r1
  drsub=dr/dble(nsub)
  dfi=fi2-fi1
  dfisub=dfi/dble(nsub)

  aij=0
  ! --- double loop over subtiles ------
  ! --- loop over theta angle ----------

  q1e3h4: do i=1,nsub
    rr1=r1+(i-1)*drsub
    rr2=rr1+drsub
    ar=0.5D0*(rr2**2-rr1**2)
    szaml=(1.D0/3.D0)*(rr2**3-rr1**3)
    rc=szaml/ar
    arel=rc

  ! --- loop over phi angle ------------

    f1l0h6: do j=1,nsub
      fij=fi1+(j-0.5D0)*dfisub
      pxj=rc*dcos(fij)
      pyj=rc*dsin(fij)
      pzj=z0
      pxij=pxi-pxj
      pyij=pyi-pyj
      pzij=pzi-pzj
      rijsq=pxij*pxij+pyij*pyij+pzij*pzij
      rij=dsqrt(rijsq)
      rij3=rijsq*rij
      aij=aij+(pxij*uxi+pyij*uyi+pzij*uzi) *arel*drsub*dfisub/rij3
    enddo f1l0h6
  enddo q1e3h4
  if (dbc) then
    if (isnan(aij)) stop "error in intwll, aij is NaN"
  endif
  end subroutine intwll

  !----------------------------------------------------------------------
  ! save amx and indx
  !
  ! This is the counterpoint method to 'readam'.  It saves a digest
  ! of of the input parameters critical to defining the matrix.
  ! These are the protein geometry parameters, the patch integration
  ! grid parameters and the permittivity constants. Then saves the
  ! 'amx' matrix itself.
  subroutine writam
  use strngs
  use geom
  implicit none
  integer :: i,j ! loop indices
  ! check variables are:
  ! patch data: npatch, dxf,dxw,nsub0,
  ! geom data: zl,rl,rlvest(),rlcurv(),epsw,epspr
  write(unit=fidlog,fmt=*)"Saving inverted 'A' matrix"
  open(unit=fidamx,file=fnamx//firun//'.dat',action='write',status='replace',err=1012)
    write(unit=fidamx,fmt='(I11)')npatch,nsub0
    write(unit=fidamx,fmt='(I11)')nint(16384*zl(1)),nint(16384*zl(4))
    write(unit=fidamx,fmt='(I11)')nint(16384*rl(1)),nint(16384*rl(4)),nint(16384*rl(5))
    write(unit=fidamx,fmt='(I11)')nint(16384*rlvest()),nint(16384*rlcurv()),nint(16384*epsw) &
      & ,nint(16384*epspr)
    write(unit=fidamx,fmt=*)  
    x2t9z8: do i=1,npatch
      write(unit=fidamx,fmt='(I4)')indx(i)
      s0j2d6: do j=1,npatch
        write(unit=fidamx,fmt='(I4,1X,I4,1X,E26.18)')i,j,amx(i,j)
      enddo s0j2d6
      write(unit=fidamx,fmt=*)  
    enddo x2t9z8
  close(unit=fidamx)
  return

1012 continue
  write(unit=fidlog,fmt=*) "Unable to open file for writing the 'A' matrix, matrix not saved"
  end subroutine writam

  !----------------------------------------------------------------------
  ! save amx and indx
  !
  ! Write out amx and indx for debugging purposes.
  subroutine dumpam(fname)
  use strngs
  use geom
  implicit none
  integer :: i,j ! loop indices
  character(len=*), intent(in) :: fname
  ! check variables are:
  ! patch data: npatch, dxf,dxw,nsub0,
  ! geom data: zl,rl,rlvest(),rlcurv(),epsw,epspr
  write(unit=fidlog,fmt=*)"Dumping 'A' matrix to ",fname
  open(unit=fidamx,file=fnamx//fname,action='write',status='replace',err=1013)
    write(unit=fidamx,fmt='(I11)')npatch
    write(unit=fidamx,fmt=*)  
    x2t9z8: do i=1,npatch
      write(unit=fidamx,fmt='(I4)')indx(i)
      s0j2d6: do j=1,npatch
        write(unit=fidamx,fmt='(I4,1X,I4,1X,E26.18)')i,j,amx(i,j)
      enddo s0j2d6
      write(unit=fidamx,fmt=*)  
    enddo x2t9z8
  close(unit=fidamx)
  return

1013 continue
  write(unit=fidlog,fmt=*) "Unable to open file for writing the 'A' matrix, matrix not saved"
  end subroutine dumpam


  end subroutine matrix

  ! --------------------
  ! 'patch' parameters from input file
  subroutine rdptch(fid,sname,svalue,istat)
  use strngs
  use conf
  implicit none
  integer, intent(in) :: fid
  character(len=*), intent(in) :: sname,svalue
  integer, intent(out) :: istat
  
  character(32) :: nme_
  character(1024) :: val_
  
  if (dbc) then
    if (sname.ne.fsptch) stop "Error: incorrect section name"
    if (len_trim(svalue).ne.0) stop "Error: section does not any parameters"
  endif

  ! ALL PATCH INPUT VARIABLES HAVE DEFAULT VALUES !

  ! process section
  b7x5c9: do
    val_ = " "
    call readnv(fid,nme_,val_,istat)
    ! exit on bad read or section 'end'
    if (istat.ne.0) return
    if (nme_.eq.fsend) exit b7x5c9
    k9w1n4: select case (nme_)
      case (fsdxf) k9w1n4
        read(val_,'(D20.13)')dxf
      case (fsdxw) k9w1n4
        read(val_,'(D20.13)')dxw
      case (fsnsub) k9w1n4
        read(val_,*)nsub0
      case (fsepsw) k9w1n4
        read(val_,'(D20.13)')epsw
      case (fsepsp) k9w1n4
        read(val_,'(D20.13)')epspr
      case default k9w1n4
        call j2l1m4("Name "//nme_//" is not valid in patch section")
    end select k9w1n4
  enddo b7x5c9

  contains

  subroutine j2l1m4(msg)
  use strngs
  implicit none
  character(len=*), intent(in) :: msg

    write(unit=fidlog,fmt=*)"Bad patch section in input:"
    write(unit=fidlog,fmt=*)msg
    write(unit=fidlog,fmt=*)"Valid names are (all are optional):"
    write(unit=fidlog,fmt=*)fsdxf,", ",fsdxw,", ",fsnsub,", ",fsepsw," and ",fsepsp
    stop 1
  end subroutine j2l1m4

  end subroutine rdptch

  !----------------------------------------------------------------------
  ! Initialise the patch system
  !
  ! Allocate storage for the patch arrays iff not homog
  subroutine rfptch
  implicit none
  logical :: isok

  homog=dfeq(epspr,epsw)
  ! NOTE: solvent homog is only used if not homog. Simply setting
  ! it here may be more efficient as we don't need to wait for dfeq
  ! to return.
  if (.not.homog) then
    call readam(isok)
    if (.not.isok) call matrix
  endif

  write(unit=fidlog,fmt='(72("-"))')
  write(unit=fidlog,fmt=*)"Interpreted 'patch' input section"
  write(unit=fidlog,fmt='(72("-"))')
  call ecptch(fidlog)
  write(unit=fidlog,fmt='(72("-"))')

  call A42878()

  contains

  ! Write out patch information
  subroutine A42878()
  implicit none
  character(*), parameter :: fngpch='dat/ptchgeom.'
  integer :: istat, ii

  open(unit=fidpch,file=fngpch//firun//'.dat',action='write',iostat=istat)
  if (istat.ne.0) stop 'Unable to open file to save patch geometry.'
  write(unit=fidpch,fmt='("# UUID ",A32)')fuuid
  write(unit=fidpch,fmt='(A)')'# xcoord  ycoord  zcoord area xnorm ynorm znorm deps'
  write(unit=fidpch,fmt='(A)')'# ang     ang     ang    ang2 ang   ang   ang   1/permittivity'
  do ii=1,npatch
    write(unit=fidpch,fmt='(F16.10,7(2x,F16.10))')prx(ii),pry(ii),prz(ii),parea(ii),pux(ii),puy(ii),puz(ii),deps(ii)
  enddo
  close(unit=fidpch)
  end subroutine A42878


  ! -----------------------------------------------------
  ! Read solved matrix information (called by readin)
  !
  ! The AMX file format contains a digest of the input parameters
  ! critical to defining the matrix. These are the protein geometry
  ! parameters, the patch integration grid parameters and the
  ! permittivity constants. If these parameters do not match when
  ! reading in a file 'imatch' is set to false and the matrix is
  ! not read in. If the matrix is successfully read in then 'imatch'
  ! is set to true (as well as the 'irdmax' internal flag).
  !
  ! @param imatch : set to false if the matrix on disk
  !     was not generated with the current parameters
  subroutine readam(imatch)
    use strngs
    use geom
    implicit none
    logical, intent(out) :: imatch

    integer, dimension(4) :: ckzl,ckrl
    integer :: ckrlvs,ckrlcv,ckepsw,ckepsp
    integer cknpch,chnsb0
    integer ipch,jpch ! patch indices
    integer ii,jj     ! dummy vars or particle indices
    integer istat     ! see if there is a file to read 
    type (geomdf), allocatable :: geodfn
    logical :: isexst ! does file exist
    character(32) :: isrdbl ! is file readable
  
    imatch=.false.
    inquire(file=fnamx//firun//'.dat',iostat=istat, exist=isexst,read=isrdbl)
    if (isexst) then
      if (istat.eq.0) then
        if (isrdbl.ne.'NO') then
          open(fidamx,file=fnamx//firun//'.dat',iostat=istat,action="read")
          if (istat.eq.0) then
            write(unit=fidlog,fmt=*)"Checking """//fnamx//firun//".dat"" for previous A matrix."
          else
            close(unit=fidamx)
          endif
        endif
      endif
    else
      istat=1
    endif
    if (istat.ne.0) then
      return
    endif
    read(unit=fidamx,fmt='(I11)',iostat=istat)cknpch,chnsb0
    if (istat.ne.0) then
      close(unit=fidamx)
      return
    endif
    read(unit=fidamx,fmt='(I11)')ckzl(1),ckzl(2)
    read(unit=fidamx,fmt='(I11)')ckrl(1),ckrl(2),ckrl(3)
    read(unit=fidamx,fmt='(I11)')ckrlvs,ckrlcv,ckepsw,ckepsp
    read(unit=fidamx,fmt=*)
    if (npatch.eq.0) then
      allocate(geodfn)
      call initgdfn(geodfn)
      call dumpgdfn(geodfn)
      call defgrd(geodfn)
      call dumpgdfn(geodfn)
      deallocate(geodfn)
    endif
    imatch=.true.
    !   verify parameters
    call i0s5u3(cknpch,npatch,imatch)
    call i0s5u3(chnsb0,nsub0,imatch)
    call i0s5u3(ckzl(1),nint(16384*zl(1)),imatch)
    call i0s5u3(ckzl(2),nint(16384*zl(4)),imatch)
    call i0s5u3(ckrl(1),nint(16384*rl(1)),imatch)
    call i0s5u3(ckrl(2),nint(16384*rl(4)),imatch)
    call i0s5u3(ckrl(3),nint(16384*rl(5)),imatch)
    call i0s5u3(ckrlvs,nint(16384*rlvest()),imatch)
    call i0s5u3(ckrlcv,nint(16384*rlcurv()),imatch)
    call i0s5u3(ckepsw,nint(16384*epsw),imatch)
    call i0s5u3(ckepsp,nint(16384*epspr),imatch)
    if (.not.imatch) then
      write(unit=fidlog,fmt=*)"A matrix from file """//fnamx//firun//".dat"" does not match current parameters."
      close(unit=fidamx)
      return
    endif
    ! Npatch defined, allocate arrays
    call f4s3s6 
    write(unit=fidlog,fmt=*)" Reading in LU decomposed A and indx "
    d8x5q1: do ipch=1,npatch
      read(unit=fidamx,fmt=*)indx(ipch)
      g7e5l2: do jpch=1,npatch
        read(unit=fidamx,fmt=*)ii,jj,amx(ipch,jpch)
        if (ii.ne.ipch) stop "Index mismatch reading A matrix"
        if (jj.ne.jpch) stop "Index mismatch reading A matrix"
      enddo g7e5l2
      read(unit=fidamx,fmt=*)
    enddo d8x5q1
    close(unit=fidamx)
    irdamx=.true.
    imatch=.true.

  end subroutine readam

  ! check that integers ii and jj match
  subroutine i0s5u3(ii,jj,imatch)
  implicit none
  integer, intent(in) :: ii,jj
  logical, intent(inout) :: imatch
  double precision :: rii, rjj
  rii=ii
  rjj=jj
  if (ii.ne.jj) then
    write(unit=fidlog,fmt=*)"Mismatch between parameters for A matrix on disk and current parameters:"
    write(unit=fidlog,fmt=*)"Read [",rii/16384,"] Current [",rjj/16384,"]"
    imatch=.false.
  endif
  end subroutine i0s5u3

  end subroutine rfptch

  !----------------------------------------------------------------------
  ! Private method to allocate second half of patch system
  ! 
  subroutine f4s3s6
  use conf
  implicit none
  if (dbc) then
    if (homog) stop "Attempt to initialise patch when not in use"
  endif
  if (.not.allocated(rip)) then
    allocate(h(npchsz),amx(npchsz,npchsz),c(npchsz),rip(npchsz,confsz())&
         &,iprip(npchsz,confsz()),indx(npchsz))
    h=0
    amx=0
    c=0
    rip=0
    iprip=0
    indx=0
  endif
  end subroutine f4s3s6
  
end module patch

! have free function for solvent eps
double precision function epsblk()
use patch
implicit none
epsblk=epsw
end function epsblk
