
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
! SPECIE AND SALT DATA
! 
!
!  @begin docfileformat
!   input file sections: Specie
!
!    specie (mobile|flexible|free|chonly)
!    (type {mobile|flexible|free|chonly})
!    z {number} # specie charge
!    d {number} # specie diameter
!    {ctrg|ntrg} {count} # indicate particle count
!    n {count}
!     # position definitions
!    name {2 letter name}
!    (chex {number}) # excess chemical potential
!    ratspc {number} # chance to move this specie
!    ratmov {number} # chance to move/jump cf jumpin/out
!    ratexc {number} # chance to jumpin/jumpout (cf. move or jump)
!    ratgr {number} # chance to add/delete ion
!    end
!
!  NOTES:
!    (1) if 'type' is not present, defaults to 'free'
!    (2) if 'n' is given the next {count} lines contain xyz coords for
!        of the specific ions. Mobile and flexible ions also require
!        the r and optionally xyz (r=mobility radius, xyz=localisation 
!        pt). If the second set of xyz points is not included, the 
!        starting point is used for the localisation point.
!  @end docfileformat
!
!  @begin docfileformat
!   input file sections: Salt
!
!    salt
!    cation "Ka" (optional if name present)
!    name "KaCl" (optional if cation present)
!    chex 1.234
!    ctarg 0.1
!    ratgr 1.234
!    end
!
!  @end docfileformat
!
!  @begin docfileformat
!   input file sections: State
!
!    state
!    name "Ka"
!    state 1 "K1"
!    state 2 "K2"
!    enthalpy 1.234 # in kJ/mol
!    entropy 1.234  # in J/(mol.K)
!    end
!
!  @end docfileformat
! ######################################################################
! CHANGES for state:
! ##################
! SALT processing:
! A salt 'specie' can now be a superspecie. Impact on
!  - generating the initial concentration of the subspecies
!  - salt add/remove selecting subspecie
!  - salt data accumulation and reporting
!
! SPECIE processing:
! A specie can now be a subspecie. Impact on
!  - generating the initial concentration from salt
!
! ADDITIONS for state:
! ####################
! Superstate species require new output routines for
!  - Data accumulation and reporting
!   - Sum of subspecie histograms
!   - RDF ? Widom energy ?
! ######################################################################
!
!  Module functions [and corresponding attributes] or attributes
!
!    xri(ispec) [species(ispec)%radius_]         = particle radii
!    xz(ispec) [species(ispec)%valency_]         = particle charge/valency
!    xq(ispec) [species(ispec)%red_charge_]      = particle charge/valency * qstar
!    chempi(ispec)                               = specie chem. pot.
!    chexi(ispec) [species(ispec)%chem_excess_]  = specie chem. excess
!    chemps(idx)                                 = salt chem. pot.
!    chemps(idx)                                 = salt chem. excess
!    nspec = number of specie
!    nstr = number of structural specie
!    imax_mob_ = number of mobile structural specie
!    imax_flex_ = number of flexible structural specie
!    * number of 'chonly' structural species is nstr-(imax_mob_+imax_flex_)
!    idxcl = chloride anion index
!    isalt(idx) [salts(idx)%cation_index_]       = salt cation index
!    nsalt = number of salts
!    fsalt(idx) [salts(idx)%code_name_]          = salt name (4char)
!    fspc(ispec) [species(ispec)%code_name_]     = specie code name (2 char)
!    spcidx(name) = index of specie with name
!
! --------------------------------------------------------------
! NOTES ON SPECIE TYPES and INPUT
!
! ** 'chonly' type
!
! Particles of a channel-only specie are restricted to movement
! anywhere within the filter region.  They can move in increments
! from the current position or jump to anywhere in the filter.
! They can not be added or removed, not can jump into or out of
! the channel.
!
! ** 'flexible' type
!
! Particles of a flexible specie move within a sphere (as per mobile) 
! but they can exist outside zlimit.
!
! ** 'mobile' type
!
! Particles of a mobile specie are restricted in two way, firstly
! they must remain within the filter region and secondly they must
! each remain within a fixed radius of a defined point.  The only
! movement possible is small displacement moves.  The can not be
! added or deleted from the simulation.
!
! ** 'free' type
!
! Particles of a free ion specie may participate in any move
! type, be a component of a salt and be added or deleted from
! the system.  NOTE: only specie of this type may be a component
! of a salt and all components of all salts must be 'free' species 
!
! POSITION INPUT
!
! When the position of particles are defined in the input file
! they must all be non-overlapping valid positions. The program
! checks and exits if this is not true.  This restriction does
! not apply for 'mobile'/'flexible' ion localisation centre-points.
!
! SALT INPUT RESTRICTIONS
!
! Salts can only have a single cation
!
! A salt code name must contain a cation's code name as
! the first two letters. 

module spec
  use const
  implicit none
  private

  ! Enumeration of specie type
  integer, parameter, private :: nfre_l=0, nmob_l=1, nflx_l=2, ncho_l=3

  ! Constants defining subtype of move
  integer, public, parameter :: sphr_l=1
  integer, public, parameter :: jump_l=2
  integer, public, parameter :: jmpin_=3
  integer, public, parameter :: jmpou_=4
  integer, public, parameter :: swap_l=5

  ! Data structure use when collecting data from the input file.
  !
  ! The input file can have species in any order and we want them
  ! to be in a specific order.  Objects of this type are used to
  ! cache the specie input until after all input has been read.
  ! At that point the species are sorted into the desired order
  ! and the input data transferred into the module's dataset.
  type specie
    ! == specie data
    ! Specie radius (in Angstrom)
    double precision :: radius_
    ! Valency (in electrons)
    double precision :: valency_
    ! Reduced charge
    double precision :: red_charge_
    ! Chem. Ex. (in internal units of kT)
    double precision :: chem_excess_
    ! == container for any particle data  (in Angstrom) < x, y, z >
    double precision, dimension(:,:), pointer :: xyz_
    ! == container for any mobile ion data (in Angstrom) < r, x, y, z >
    double precision, dimension(:,:), pointer :: sxyzr_
    ! spctyp == nmob_l,nflx_p,nstr_l,nfre_l 
    integer :: type_ 
    ! Initial specie count
    integer :: input_count_
    ! specie code name
    character(2) :: code_name_
    ! MC rate parameters for this specie: 
    ! relative chance  for specie (relative probability)
    double precision :: rate_specie_
    ! chance to jump in/out or move/jump (relative probability)
    double precision :: rate_exchange_
    ! chance to move or jump (relative probability)
    double precision :: rate_move_  !!?? UNUSED ??!!
    ! chance to add/delete (relative probability)
    double precision :: rate_change_
    ! chance that an add/delete occurs in a region (relative probability)
    double precision, dimension(nrgnmx) :: rate_region_
    ! target salt concentrations (in SI Molar (mol/l)
    double precision :: target_concentration_
  end type specie

  type salt
    ! salt add/delete rate (relative probability)
    double precision :: rate_change_
    ! target salt concentrations (in SI Molar (mol/l)
    double precision :: target_concentration_
    ! index of cation specie for each salt
    integer :: cation_index_
    ! salt names
    character(4) :: code_name_
  end type salt

  type subspec
    ! super specie name
    character(2) :: code_name_
    ! The probability of swapping between these two sub
    ! species. (relative probability)
    double precision :: rate_swap_
    ! free energy terms (in S.I. units: kJ/mol)
    double precision :: enthalpy_
    ! free energy terms (in S.I. units: J/mol.K)
    double precision :: entropy_
    ! sub-specie indices, ground state specie is at index 1
    integer, dimension(2) :: sub_indices_
    ! sub-specie names, ground state specie is at index 1
    character(2), dimension(2) :: sub_names_
  end type subspec

  ! -----------------
  ! SPECIE PARAMETERS
  ! ----------
  type (specie), private, dimension(:), allocatable :: species
  ! DATA LIST SIZES AND INDICES
  ! nspec_ : number of specie
  ! idxcl_ : index of chloride specie (i<idxcl --> notfree: i>=idxcl --> free)
  ! nstr_  : number of structural ions
  ! imax_mob_ : maximum index of mobile ions
  ! imax_flex_ : maximum index of flexible ions
  ! NOTE because of sorting by type:
  !       number of mobile ions = imax_mob_ ( - 0 )
  !       number of flexible ions = imax_flex_ - imax_mob_
  !       number of channel-only ions = nstr_ - imax_flex_
  ! INVARIANTS nstr_ + 1 == idxcl_
  !            nspec_ > idxcl_
  ! (invariants arise because one ion must be chloride
  !  and there must be at least one anion)
  integer, private :: nspec_ = 0, idxcl_ = 0, nstr_ = 0, imax_mob_ = 0, imax_flex_ = 0

  ! ---------------
  ! SALT PARAMETERS
  !
  type (salt), private, dimension(:), allocatable :: salts
  ! Number of salt definitions
  integer, private :: nsalt_ = 0

  ! ---------------
  ! SUBSPECIE PARAMETERS
  !
  type (subspec), private, dimension(:), allocatable :: subspecies
  ! Number of subspecie definitions
  integer, private :: nsubs_ = 0

  double precision, private :: enlrg !!!

  ! ------------------------------------------------------------
  ! Public methods
  public :: nfree, have_localized, nspec, nsalt, nsubs, idxcl

  ! SPECIE ATTRIBUTE ACCESS
  public :: xri, xz, xq, chempi, chexi, deltmp, ratspc, ratexc, ratmov, ratgri
  public :: fspc, chexi_set, ratreg, mobchk, mobpen, input_count

  ! SALT ATTRIBUTE ACCESS
  public :: chemps, ratgrs, isalt, fsalt

  ! SUBSPECIE ATTRIBUTE ACCESS
  public :: subspecie_index

  ! random trial selection functions
  ! **select move may call ranff() 
  public :: select_salt, select_gc_specie, select_mv_specie, select_subs
  public :: select_region, select_move

  ! smallest allowable distance between two particles of given species
  public :: dd_get 

  ! read specie/salt sections of input file
  public :: rdspec, rdsalt, rdsubs, rfspec

  ! specie query functions
  public :: chonly, flexible, mobile, isfree, localized

  ! Get information about the specie initial particle configuration
  ! read in from the input file.
  public :: ctargi, ionstr, ctargs, struks, spcidx

  ! Output procedures
  public :: output_specie_potentials, output_salt_potentials
  public :: ecspec, ecsalt, ecsubs, help_salt, help_specie, help_subspecie

  interface rsx
    module procedure rsx1
    module procedure rsx2
  end interface rsx

  interface rsy
    module procedure rsy1
    module procedure rsy2
  end interface rsy

  interface rsz
    module procedure rsz1
    module procedure rsz2
  end interface rsz

  interface rsr
    module procedure rsr1
    module procedure rsr2
  end interface rsr

  public rsx, rsy, rsz, rsr
contains

  ! ------------------------------------------------------------
  ! Copy a salt
  subroutine assign_salt(lhs, rhs)
    implicit none
    type(salt), intent(inout) :: lhs
    type(salt), intent(inout) :: rhs
    lhs%rate_change_    = rhs%rate_change_
    lhs%target_concentration_ = rhs%target_concentration_
    lhs%cation_index_   = rhs%cation_index_
    lhs%code_name_      = rhs%code_name_
  end subroutine assign_salt

  ! ------------------------------------------------------------
  ! Copy a spece
  subroutine move_specie(lhs, rhs)
    implicit none
    type(specie), intent(inout) :: lhs
    type(specie), intent(inout) :: rhs
    lhs%radius_         = rhs%radius_         
    lhs%valency_        = rhs%valency_        
    lhs%red_charge_     = rhs%red_charge_     
    lhs%chem_excess_    = rhs%chem_excess_    
    lhs%type_           = rhs%type_           
    lhs%input_count_    = rhs%input_count_    
    lhs%code_name_      = rhs%code_name_      
    lhs%rate_specie_    = rhs%rate_specie_    
    lhs%rate_exchange_  = rhs%rate_exchange_  
    lhs%rate_move_      = rhs%rate_move_      
    lhs%rate_change_    = rhs%rate_change_    
    lhs%rate_region_    = rhs%rate_region_    
    lhs%target_concentration_ = rhs%target_concentration_ 
    lhs%xyz_            = rhs%xyz_ 
    lhs%sxyzr_          = rhs%sxyzr_ 
    nullify(rhs%xyz_)
    nullify(rhs%sxyzr_)
  end subroutine move_specie

  ! ------------------------------------------------------------
  ! Initialise a specie object
  !
  subroutine ctor__(self)
    implicit none
    type(specie), intent(inout) :: self
    self%radius_         = 0.D0
    self%valency_        = 0.D0
    self%red_charge_     = 0.D0
    self%chem_excess_    = 0.D0
    self%type_           = -1
    self%input_count_    = 0
    self%code_name_      = '  '
    self%rate_specie_    = 0.D0
    self%rate_exchange_  = 0.D0
    self%rate_move_      = 0.D0
    self%rate_change_    = 0.D0
    self%rate_region_    = 0.D0
    self%target_concentration_ = 0.D0
    nullify(self%xyz_)
    nullify(self%sxyzr_)

  end subroutine ctor__

  ! ------------------------------------------------------------
  ! initialise a salt
  subroutine ctor_salt(self)
    implicit none
    type(salt), intent(inout) :: self
    self%rate_change_    = 0.D0
    self%target_concentration_ = 0.D0
    self%cation_index_   = 0
    self%code_name_      = '    '
  end subroutine ctor_salt

  ! ------------------------------------------------------------
  ! initialise a subspecie
  subroutine ctor_subs(self)
    implicit none
    type(subspec), intent(inout) :: self
    self%code_name_    = '  '
    self%rate_swap_    = 0.D0
    self%enthalpy_     = 0.D0
    self%entropy_      = 0.D0
    self%sub_indices_  = 0
    self%sub_names_    = ' '

  end subroutine ctor_subs

  ! --------------------------------------------------
  ! Is this a structural channel-only specie?
  !
  ! Particles of a channel-only specie are restricted to movement
  ! anywhere within the filter region.  They can move in increments
  ! from the current position or jump to anywhere in the filter.
  ! They can not be added or removed, not can jump into or out of
  ! the channel.
  !
  ! chonly species are structural ions that are _not_ 'mobile' and
  ! _not_ 'flexible'.
  logical function chonly(ispec)
    implicit none
    integer, intent(in) :: ispec
    if (dbc) then
      if (dbc_level.ge.dbc_index) then
        if (ispec.lt.1.or.ispec.gt.nspec_) stop "Specie index out of range"
      endif
    endif
    chonly=(species(ispec)%type_.eq.ncho_l)
  end function chonly

  ! --------------------------------------------------
  ! Calculate a specie's chem. pot.
  double precision function chempi(ispec)
    implicit none
    integer, intent(in) :: ispec
    double precision :: rtargi
    if (dbc) then
      if (dbc_level.ge.dbc_require) then
        if (dfeq(species(ispec)%target_concentration_, 0.D0)) stop &
          & "Error: Want chemical potential but have zero concentration"
      endif
      if (dbc_level.ge.dbc_index) then
        if (ispec.lt.1.or.ispec.gt.nspec_) stop &
          & "Error: Specie index out of range"
      endif
    endif
    rtargi = species(ispec)%target_concentration_ / tosi
    chempi = species(ispec)%chem_excess_ + log(rtargi)
  end function chempi

  ! --------------------------------------------------
  ! Calculate a salt's chem. pot.
  double precision function chemps(igc)
    implicit none
    integer, intent(in) :: igc
    integer :: iv, isalt
    if (dbc) then
      if (dbc_level.ge.dbc_index) then
        if (igc.lt.1.or.igc.gt.nsalt_) stop &
          &"Error: Salt index out of range"
      endif
    endif
    isalt = salts(igc)%cation_index_
    iv = nint(species(isalt)%valency_)
    ! NOTE use of individual chemical potentials
    chemps = chempi(isalt) + iv * chempi(idxcl_)
  end function chemps

  ! --------------------------------------------------
  ! set specie chem. excess.
  !
  ! REQUIRE: idxcl <= ispec <= nspec AND ctargi(ispec) /= 0 
  subroutine chexi_set(ispec, val)
    implicit none
    integer, intent(in) :: ispec
    double precision, intent(in) :: val
    ! calculate ex. chem. pot. from  chem. pot.
    if (dbc) then
      if (dbc_level.ge.dbc_require) then
        if (ispec.lt.idxcl_.or.ispec.gt.nspec_) stop "Error: Specie index out of range"
      endif
    endif
    species(ispec)%chem_excess_ = val
  end subroutine chexi_set

  ! --------------------------------------------------
  ! Get a specie chem. excess
  double precision function chexi(ispec)
    implicit none
    integer, intent(in) :: ispec
    if (dbc) then
      if (dbc_level.ge.dbc_index) then
        if (ispec.lt.1.or.ispec.gt.nspec_) stop "Error: Specie index out of range"
      endif
    endif
    chexi=species(ispec)%chem_excess_
  end function chexi

  ! --------------------------------------------------
  ! Target concentration of a specie
  !
  ! The target concentration of specie is the sum of the partial
  ! target concentrations from all salts it is a component of.
  !
  ! @pre isfree(ispec)
  double precision function ctargi(ispec)
    implicit none
    integer, intent(in) :: ispec
    if (dbc) then
      if (dbc_level.ge.dbc_index) then
        if (ispec.lt.1.or.ispec.gt.nspec_) stop "Error: Specie index out of range"
      endif
    endif
    ctargi = species(ispec)%target_concentration_
  end function ctargi

  ! --------------------------------------------------
  ! Target concentration of a salt.
  !
  ! This is a fixed value from the input file giving the desired
  ! target concentration of a salt.  !
  double precision function ctargs(igc)
    implicit none
    integer, intent(in) :: igc
    if (dbc) then
      if (dbc_level.ge.dbc_index) then
        if (igc.lt.1.or.igc.gt.nsalt_) stop "Error: Specie index out of range"
      endif
    endif
    ctargs=salts(igc)%target_concentration_
  end function ctargs

  ! --------------------------------------------------
  ! Minimum distance between two particles of the given types
  !
  ! This is the sum of the individual radii for each specie
  ! 
  double precision function dd_get(ispec,jspec)
    implicit none
    integer, intent(in) :: ispec,jspec
    if (dbc) then
      if (dbc_level.ge.dbc_index) then
        if (ispec.lt.1.or.ispec.gt.nspec_) stop "Error: Specie index out of range"
        if (jspec.lt.1.or.jspec.gt.nspec_) stop "Error: Specie index out of range"
      endif
    endif
    dd_get=species(ispec)%radius_+species(jspec)%radius_
  end function dd_get

  ! --------------------------------------------------
  ! Delete unused data sets 
  !
  ! remove temporary data only used during reading of
  ! input file. 
  !
  ! Note: The input file allows particle positions to be specified
  ! in the specie section.  However, this information needs to be
  ! transferred to the 'conf' module.  Module inclusion ordering
  ! means that the specie module is initialised first and the conf
  ! module must pull the information from this module.  The 'conf'
  ! module therefore needs to call this module to indicate when
  ! it has finised using this data.  This protocol was favoured
  ! over a scheme that transferred data ownership of the data to
  ! the conf module.
  subroutine deltmp(ispec)
    implicit none
    integer, intent(in) :: ispec
    if (dbc) then
      if (dbc_level.ge.dbc_index) then
        if (ispec.lt.1.or.ispec.gt.nspec_) stop "Error: Specie index out of range"
      endif
    endif
    if (associated(species(ispec)%xyz_)) then
      deallocate(species(ispec)%xyz_)
      nullify(species(ispec)%xyz_)
    endif
  end subroutine deltmp

  ! --------------------------------------------------
  ! Echo individual spec section
  !
  ! Write out the interpreted specie data from the input
  ! file in the same format as an input file.  This includes program
  ! default values for optional input data and normalised rate
  ! values.
  subroutine ecspec(fid,ispec,xyzs,sz)
    use strngs
    implicit none
    integer, intent(in) :: fid, ispec, sz
    double precision, dimension(:,:), intent(in) :: xyzs
    ! locals
    integer :: idx, jdx ! loop counters
    character(20) :: fltout
    if (dbc) then
      if (dbc_level.ge.dbc_require) then
        call require_writable(fid)
      endif
      if (dbc_level.ge.dbc_index) then
        if (ispec.lt.1.or.ispec.gt.nspec_) stop "Error: Specie index out of range"
 !       if (sz.gt.0) then
 !TODO         if (size(xyzs).lt.sz) stop "Requested array size is larger than actual size"
 !       endif
      endif
    endif

    write(unit=fid,fmt='(A)')fsspec
    write(unit=fid,fmt='(A,1X,"""",A,"""")')fsname,species(ispec)%code_name_
    call str(species(ispec)%rate_specie_, fltout)
    write(unit=fid,fmt='(A,1X,A)')fsrtsp,trim(adjustl(fltout))
    call str(species(ispec)%rate_move_, fltout)
    write(unit=fid,fmt='(A,1X,A)')fsrtmv,trim(adjustl(fltout))
    call str(species(ispec)%rate_exchange_, fltout)
    write(unit=fid,fmt='(A,1X,A)')fsrtex,trim(adjustl(fltout))
    call str(species(ispec)%rate_change_, fltout)
    write(unit=fid,fmt='(A,1X,A)')fsrtgr,trim(adjustl(fltout))

    write(unit=fid,fmt='(A)',advance='NO')fsrtrg
    do idx=1,nrgnmx
      call str(species(ispec)%rate_region_(idx), fltout)
      write(unit=fid,fmt='(1X,A)',advance='NO') trim(adjustl(fltout))
    enddo
    write(unit=fid,fmt='("")')
    call str(species(ispec)%valency_, fltout)
    write(unit=fid,fmt='(A,1X,A)')fsz,trim(adjustl(fltout))
    call str(species(ispec)%radius_*2, fltout)
    write(unit=fid,fmt='(A,1X,A)')fsd,trim(adjustl(fltout))
    call str(species(ispec)%chem_excess_, fltout)
    write(unit=fid,fmt='(A,1X,A)')fschex,trim(adjustl(fltout))

    ! determine the specie type
    if (chonly(ispec)) then
      write(unit=fid,fmt='(A,1X,A)')fstype,fschon
    elseif (mobile(ispec)) then
      write(unit=fid,fmt='(A,1X,A)')fstype,fsmobl
    elseif (flexible(ispec)) then
      write(unit=fid,fmt='(A,1X,A)')fstype,fsflxd
    else
      write(unit=fid,fmt='(A,1X,A)')fstype,fsfree
    endif
    ! for non-free types print x y z (r) coordinates
    if (sz.ne.0) then
      write(unit=fid,fmt='(A,1X,I4)')fsn,sz
      if (mobile(ispec).or.flexible(ispec)) then
        do idx=1,sz
          do jdx=1,3
            call str(xyzs(jdx,idx), fltout)
            write(unit=fid,fmt='(1X,A)',advance='NO') trim(adjustl(fltout))
          enddo
          call str(sqrt(species(ispec)%sxyzr_(1,idx)), fltout)
          write(unit=fid,fmt='(1X,A)',advance='NO') trim(adjustl(fltout))
          do jdx=2,4
            call str(species(ispec)%sxyzr_(jdx,idx), fltout)
            write(unit=fid,fmt='(1X,A)',advance='NO') trim(adjustl(fltout))
          enddo
          write(unit=fid,fmt='("")')
        enddo
      else
        do idx=1,sz
          do jdx=1,3
            call str(xyzs(jdx,idx), fltout)
            write(unit=fid,fmt='(1X,A)',advance='NO') trim(adjustl(fltout))
          enddo
          write(unit=fid,fmt='("")')
        enddo
      endif
    endif
    write(unit=fid,fmt='(A)')fsend
    write(unit=fid,fmt=*)
  end subroutine ecspec

  ! --------------------------------------------------
  ! Write all the salt parameters
  !
  ! Write out the interpreted specie data from the input
  ! file in the same format as an input file.  This includes program
  ! default values for optional input data and normalised rate
  ! values.
  subroutine ecsalt(fid)
    use strngs
    implicit none
    integer, intent(in) :: fid
    integer :: igc
    character(20) :: fltout
    if (dbc) then
      if (dbc_level.ge.dbc_require) then
        call require_writable(fid)
      endif
    end if
    do igc=1,nsalt_
      write(unit=fid,fmt='(A)')fssalt
      write(unit=fid,fmt='(A,1X,"""",A,"""")')fsname,salts(igc)%code_name_
      call str(salts(igc)%rate_change_, fltout)
      write(unit=fid,fmt='(A,1X,A)')fsrtgr,trim(adjustl(fltout))
      call str(salts(igc)%target_concentration_, fltout)
      write(unit=fid,fmt='(A,1X,A)')fsctrg,trim(adjustl(fltout))
      write(unit=fid,fmt='(A)')fsend
      write(unit=fid,fmt=*)
    enddo
  end subroutine ecsalt

  ! --------------------------------------------------
  ! Write all the subs parameters
  !
  ! Write out the interpreted specie data from the input
  ! file in the same format as an input file.  This includes program
  ! default values for optional input data and normalised rate
  ! values.
  subroutine ecsubs(fid)
    use strngs
    implicit none
    integer, intent(in) :: fid
    integer :: igc
    character(20) :: fltout
    if (dbc) then
      if (dbc_level.ge.dbc_require) then
        call require_writable(fid)
      endif
    endif
    if (nsubs_.ge.1) then
      do igc=1,nsubs_
        write(unit=fid,fmt='(A)')fssubs
        write(unit=fid,fmt='(A,1X,"""",A,"""")')fsname,subspecies(igc)%code_name_
        call str(subspecies(igc)%rate_swap_, fltout)
        write(unit=fid,fmt='(A,1X,A)')fsrtsw,trim(adjustl(fltout))
        call str(subspecies(igc)%enthalpy_, fltout)
        write(unit=fid,fmt='(A,1X,A)')fsenth,trim(adjustl(fltout))
        call str(subspecies(igc)%entropy_, fltout)
        write(unit=fid,fmt='(A,1X,A)')fsentr,trim(adjustl(fltout))
        write(unit=fid,fmt='(A,1X,"""",A,"""")')fsgrnd,subspecies(igc)%sub_names_(1)
        write(unit=fid,fmt='(A,1X,"""",A,"""")')fsexct,subspecies(igc)%sub_names_(2)
        write(unit=fid,fmt='(A)')fsend
        write(unit=fid,fmt=*)
      enddo
    endif
  end subroutine ecsubs

  ! --------------------------------------------------
  ! specie code names
  character(2) function fspc(ispec)
    implicit none
    integer, intent(in) :: ispec
    if (dbc) then
      if (dbc_level.ge.dbc_index) then
        if (ispec.lt.1.or.ispec.gt.nspec_) stop "Error: Specie index out of range"
      endif
    endif
    fspc=species(ispec)%code_name_
  end function fspc

  ! --------------------------------------------------
  ! salt code names
  character(4) function fsalt(ispec)
    implicit none
    integer, intent(in) :: ispec
    if (dbc) then
      if (dbc_level.ge.dbc_index) then
        if (ispec.lt.1.or.ispec.gt.nsalt_) stop "Error: Salt index out of range"
      endif
    endif
    fsalt=salts(ispec)%code_name_
  end function fsalt

  ! --------------------------------------------------
  ! Is this a flexible structural specie type?
  !
  ! Particles of a flexible specie do move within a sphere, but can not
  ! be added or deleted.  Unlike the other structural ions they may exist
  ! outside zlimit.
  !
  ! flexible species are structural ions that are _not_ 'mobile' and
  ! _not_ 'chonly'
  logical function flexible(ispec)
    implicit none
    integer, intent(in) :: ispec
    if (dbc) then
      if (dbc_level.ge.dbc_index) then
        if (ispec.lt.1.or.ispec.gt.nspec_) stop "Error: Specie index out of range"
      endif
    endif
    flexible=(species(ispec)%type_.eq.nflx_l)
  end function flexible

  ! --------------------------------------------------
  ! The number of localised ions
  pure logical function have_localized()
    implicit none
    have_localized=(imax_mob_+imax_flex_).ne.0
  end function have_localized

  ! ------------------------------------------------------------
  ! Print a help message about the specie input
  !
  !
  subroutine help_salt(fid)
    use strngs
    implicit none
    integer, intent(in) :: fid
    if (dbc) then
      if (dbc_level.ge.dbc_require) then
        call require_writable(fid)
      endif
    endif
    write(unit=fidlog,fmt=*)"Salt input section definition:"
    write(unit=fid,fmt='(A6," : ",A)')fsctrg,"[required, number] target salt concentration"
    write(unit=fid,fmt='(A6," : ",A)')fsname,"*[required, four letters] salt code name"
    write(unit=fid,fmt='(A6," : ",A)')fsislt,"*[required, two letters] specie code name of the cation"
    write(unit=fid,fmt='(9X,5A)')             "* Only one of '",fsname,"' or '",fsislt,"' is required"
    write(unit=fid,fmt='(A6," : ",A)')fsrtgr,"[required, number] Probability this salt is used in a "
    write(unit=fid,fmt='(9X,A)')"grand-canonical trial"
  end subroutine help_salt

  ! ------------------------------------------------------------
  ! Print a help message about the specie input
  !
  !
  subroutine help_specie(fid)
    use strngs
    implicit none
    integer, intent(in) :: fid
    if (dbc) then
      if (dbc_level.ge.dbc_require) then
        call require_writable(fid)
      endif
    endif
    write(unit=fid,fmt=*)"Specie input section definition:"
    write(unit=fid,fmt='(A6," : ",A)')fsz,"[required, number] specie valency"
    write(unit=fid,fmt='(A6," : ",A)')fsd,"[required, number] specie diameter"
    write(unit=fid,fmt='(A6," : ",A)')fsname,"[required, two letters] specie code name"
    write(unit=fid,fmt='(A6," : ",3(A))')fstype,"[optional, default type is ",fsfree,&
         & "] specie model type"
    write(unit=fid,fmt='(9X,9(A))')"Argument to ",fstype," should be one of: """,fsmobl,&
         & """, """,fsflxd,""", """,fschon,""" or"
    write(unit=fid,fmt='(9X,3(A))')"""",fsfree,""" (only first 3 letters required)"
    write(unit=fid,fmt='(A6," : ",3(A))')fschex,"[required for ",fsfree,&
         & " type only, number] initial chemical excess"
    write(unit=fid,fmt='(A6," : ",9(A))')fsrtsp,"[optional for ",fsmobl,", ",fsflxd,&
         & " and ",fschon,", required for ",fsfree,", number]"
    write(unit=fid,fmt='(9X,A)')" Probability this specie is used in a move trial"
    write(unit=fid,fmt='(A6," : ",3(A))')fsrtex,"[",fsfree," type only, required, number] Once selected for a move, "
    write(unit=fid,fmt='(9X,A)')"probability of channel insertion (vs step) move"
    write(unit=fid,fmt='(A6," : ",3(A))')fsrtmv,"[",fsfree," type only, required, number] Once selected for a step move, "
    write(unit=fid,fmt='(9X,A)')"probability of a gas phase or liquid phase move"
    write(unit=fid,fmt='(A6," : ",3(A))')fsrtgr,"[",fsfree," type only, optional, number] Probability this specie is "
    write(unit=fid,fmt='(9X,A)')"used in a individual ion grand-canonical trial"
    write(unit=fid,fmt='(A6," : ",3(A))')fsrtrg,"[",fsfree," type only, required*, four numbers] Per-region probability "
    write(unit=fid,fmt='(9X,A)')"this specie is inserted into the particular region in a grand-"
    write(unit=fid,fmt='(9X,A)')"canonical trial. *Required if this specie is part of a salt or "
    write(unit=fid,fmt='(9X,3(A))')"if ",fsrtgr," is set"
    write(unit=fid,fmt='(A6," : ",5(A))')fsn,"[required for ",fsmobl," and ",fsflxd,&
         &", optional for other types, integer] "
    write(unit=fid,fmt='(9X,3(A))')"Initial position definition flag. The lines following """,&
         & fsn,""" contain "
    write(unit=fid,fmt='(9X,A)')"X number of lines of x,y,z etc information, where X is the argument "
    write(unit=fid,fmt='(9X,9(A))')"to """,fsn,""". This tag must appear after the ",fstype,&
         & " tag for ",fsmobl," and ",fsflxd," "
    write(unit=fid,fmt='(9X,5(A))')"species. ",fschon," or ",fsfree," require three numbers defining the initial"
    write(unit=fid,fmt='(9X,5(A))')"x,y,z position of the particle. ",fsmobl," and ",fsflxd,&
         &" additionally require "
    write(unit=fid,fmt='(9X,A)')"a fourth number with the localisation radius, and an optional three "
    write(unit=fid,fmt='(9X,A)')"numbers defining the x,y,z position of the localisation centre-point. "
    write(unit=fid,fmt='(9X,A)')"If these last three numbers are not included then the initial x,y,z "
    write(unit=fid,fmt='(9X,A)')"position is used as the localisation centre-point."

  end subroutine help_specie

  ! ------------------------------------------------------------
  ! Print a help message about the sub-specie input
  !
  !
  subroutine help_subspecie(fid)
    use strngs
    implicit none
    integer, intent(in) :: fid
    if (dbc) then
      if (dbc_level.ge.dbc_require) then
        call require_writable(fid)
      endif
    endif

    write(unit=fidlog,fmt=*)"Subspecie input section definition:"
    write(unit=fidlog,fmt='(A8,": ",A)')fsname,"The name of the sub-specie"
    write(unit=fidlog,fmt='(A8,": ",A)')fsrtsw,"The probability of an exchange between subspecies"
    write(unit=fidlog,fmt='(A8,": ",A)')fsenth,"The enthalpy of difference between subspecies (in SI units)"
    write(unit=fidlog,fmt='(A8,": ",A)')fsentr,"The entropy of difference between subspecies (in SI units)"
    write(unit=fidlog,fmt='(A8,": ",A)')fsgrnd,"The code name of the ground state specie definition that"
    write(unit=fidlog,fmt='(11X,A)')"forms part of the sub-specie pair."
    write(unit=fidlog,fmt='(A8,": ",A)')fsexct,"The code name of the excited state specie definition that"
    write(unit=fidlog,fmt='(11X,A)')"forms part of the sub-specie pair."
  end subroutine help_subspecie


  ! ------------------------------------------------------------
  ! The number of active species definitions
  pure integer function idxcl()
    implicit none
    idxcl = idxcl_
  end function idxcl

  ! --------------------------------------------------
  ! Target ionic strength.
  !
  ! This is the sum of the partial concentrations of every
  ! component of every salt.
  !
  pure double precision function ionstr()
    implicit none
    integer :: igc
    ionstr=0
    do igc=1,nsalt_
      ionstr=ionstr+salts(igc)%target_concentration_*(species(salts(igc)%cation_index_)%valency_+1)
    enddo
  end function ionstr

  ! --------------------------------------------------
  ! Number of particles defined in input for a species.
  !
  integer function input_count(ispec)
    implicit none
    integer, intent(in) :: ispec
    if (dbc) then
      if (dbc_level.ge.dbc_index) then
        if (ispec.lt.1.or.ispec.gt.nspec_) stop "Error: Specie index out of range"
      endif
    endif
    input_count = species(ispec)%input_count_
  end function input_count

  ! --------------------------------------------------
  ! Is this a non-structural (ie free) ion?
  !
  ! Particles of a free ion specie may participate in any move
  ! type, be a component of a salt and be added or deleted from
  ! the system.
  ! isfree species are _not_ structural ions
  logical function isfree(ispec)
    implicit none
    integer, intent(in) :: ispec
    if (dbc) then
      if (dbc_level.ge.dbc_index) then
        if (ispec.lt.1.or.ispec.gt.nspec_) stop "Error: Specie index out of range"
      endif
    endif
    isfree=(species(ispec)%type_.eq.nfre_l)
  end function isfree

  ! --------------------------------------------------
  ! Is this a member of a subspecie pair?
  !
  ! Particles of a free ion specie may participate in a
  ! subspecie pair.
  logical function issubs(ispec)
    implicit none
    integer, intent(in) :: ispec
    integer :: idx
    if (dbc) then
      if (dbc_level.ge.dbc_index) then
        if (ispec.lt.1.or.ispec.gt.nspec_) stop "Error: Specie index out of range"
      endif
    endif
    issubs=nsubs_.ne.0
    if (issubs) then
      issubs=(species(ispec)%type_.eq.nfre_l)
      if (issubs) then
        do idx=1,nsubs_
          if (ispec.eq.subspecies(idx)%sub_indices_(1)) return
          if (ispec.eq.subspecies(idx)%sub_indices_(2)) return
        enddo
        issubs=.false.
      endif
    endif
  end function issubs

  ! --------------------------------------------------
  ! index of cation specie for each salt
  integer function isalt(ispec)
    implicit none
    integer, intent(in) :: ispec
    if (dbc) then
      if (dbc_level.ge.dbc_index) then
        if (ispec.lt.1.or.ispec.gt.nsalt_) stop "Error: Salt index out of range"
      endif
    endif
    isalt=salts(ispec)%cation_index_
  end function isalt

  ! --------------------------------------------------
  ! Is this a localized structural specie type?
  !
  ! Particles of a localized specie move within a sphere, and can not
  ! be added or deleted. Localized species are structural ions that are
  ! either 'mobile' or 'flexible'
  logical function localized(ispec)
    implicit none
    integer, intent(in) :: ispec
    if (dbc) then
      if (dbc_level.ge.dbc_index) then
        if (ispec.lt.1.or.ispec.gt.nspec_) stop "Error: Specie index out of range"
      endif
    endif
    localized=(species(ispec)%type_.eq.nflx_l).or.(species(ispec)%type_.eq.nmob_l)
  end function localized

  ! --------------------------------------------------
  ! Is this a mobile structural ion specie type?
  !
  ! Particles of a mobile specie are restricted in two way, firstly
  ! they must remain within the filter region and secondly they must
  ! each remain within a fixed radius of a defined point.  The only
  ! movement possible is small displacement moves.  The can not be
  ! added or deleted from the simulation.
  !
  ! mobile species are structural ions that are _not_ 'chonly' and
  ! _not_ 'flexible'
  logical function mobile(ispec)
    implicit none
    integer, intent(in) :: ispec
    if (dbc) then
      if (dbc_level.ge.dbc_index) then
        if (ispec.lt.1.or.ispec.gt.nspec_) stop "Error: Specie index out of range"
      endif
    endif
    mobile=(species(ispec)%type_.eq.nmob_l)
  end function mobile

  ! ------------------------------------------------------------
  ! Check if mobile particle is within particle specific constraint
  !
  ! REQUIRE mobile(ispcbk(ii)) or flexible(ispcbk(ii))
  subroutine mobchk(ii, xnew, ynew, znew, ovrlap)
    implicit none
    integer, intent(in) :: ii
    double precision, intent(in) :: xnew, ynew, znew
    logical, intent(out) :: ovrlap
    integer :: ispec, jspec, idx
    double precision :: rsr, rsx, rsy, rsz
    idx = ii
    ispec = 0
    ovrlap = .false.
    do jspec=1,imax_mob_ + imax_flex_
      if (idx.le.species(jspec)%input_count_) then
        ispec = jspec
        exit
      endif
      idx = idx - species(jspec)%input_count_
    enddo
    if (ispec.eq.0) then
      stop "particle is not correct specie type in test of mobile ion position"
    endif
    if (dbc) then
      if (dbc_level.ge.dbc_check) then
        if (.not.(species(ispec)%type_.eq.nmob_l.or.species(ispec)%type_.eq.nflx_l)) then
          stop "particle is not correct specie type in test of mobile ion position"
        endif
      endif
    endif
    rsr = species(ispec)%sxyzr_(1,idx)
    rsx = species(ispec)%sxyzr_(2,idx)
    rsy = species(ispec)%sxyzr_(3,idx)
    rsz = species(ispec)%sxyzr_(4,idx)
    ovrlap = ((sqr(xnew-rsx) + sqr(ynew-rsy) + sqr(znew-rsz)).gt.(rsr))
  end subroutine mobchk

  ! ------------------------------------------------------------
  ! Calculate hooke's law energy penalty of displacement from 
  ! center-point normalised to [0:1].
  subroutine mobpen(ii, xnew, ynew, znew, penlty)
    implicit none
    integer, intent(in) :: ii
    double precision, intent(in) :: xnew, ynew, znew
    double precision, intent(out) :: penlty
    integer :: ispec, jspec, idx
    double precision :: rsr, rsx, rsy, rsz
    penlty = 0.D0
    idx = ii
    ispec = 0
    do jspec=1,imax_mob_ + imax_flex_
      if (idx.le.species(jspec)%input_count_) then
        ispec = jspec
        exit
      endif
      idx = idx - species(jspec)%input_count_
    enddo
    if (ispec.eq.0) then
      stop "particle is not correct specie type in test of mobile ion position"
    endif
    if (dbc) then
      if (dbc_level.ge.dbc_check) then
        if (.not.(species(ispec)%type_.eq.nmob_l.or.species(ispec)%type_.eq.nflx_l)) then
          stop "particle is not correct specie type in test of mobile ion position"
        endif
      endif
    endif
    rsr = species(ispec)%sxyzr_(1,idx)
    rsx = species(ispec)%sxyzr_(2,idx)
    rsy = species(ispec)%sxyzr_(3,idx)
    rsz = species(ispec)%sxyzr_(4,idx)
    penlty = (sqr(xnew-rsx) + sqr(ynew-rsy) + sqr(znew-rsz))/rsr
  end subroutine mobpen

  ! --------------------------------------------------
  ! The number of free ions
  pure integer function nfree()
    implicit none
    nfree=nspec_-nstr_
  end function nfree

  ! ------------------------------------------------------------
  ! The number of active species definitions
  pure integer function nspec()
    implicit none
    nspec = nspec_
  end function nspec

  ! ------------------------------------------------------------
  ! The number of salt definitions
  pure integer function nsalt()
    implicit none
    nsalt = nsalt_
  end function nsalt

  ! ------------------------------------------------------------
  ! The number of subspecie definitions
  pure integer function nsubs()
    implicit none
    nsubs = nsubs_
  end function nsubs

  ! --------------------------------------------------
  ! Output a table of the specie chemical potentials
  subroutine output_specie_potentials(thefid)
    use strngs
    implicit none
    integer, intent(in) :: thefid
    integer :: ispec
    character(15) :: fltout1, fltout2, fltout3
    if (dbc) then
      if (dbc_level.ge.dbc_require) then
        call require_writable(thefid)
      endif
    end if
    ! ---- DEFINITION OF VARIABLES FOR INDIVIDUAL GC -----------
    write(unit=thefid,fmt='(1X,A6,3(1X,A15))')"SPECIE","TARGET CONC.","CHEM. EXCESS","CHEM. POT."
    w9m6j8: do ispec=nstr_+1,nspec_
      call str (species(ispec)%chem_excess_, fltout1)
      call str (chempi(ispec), fltout2)
      call str (species(ispec)%target_concentration_, fltout3)
      write(unit=thefid,fmt='(1X,A6,3(1X,A15))')species(ispec)%code_name_,fltout3,fltout1,fltout2
    enddo w9m6j8
  end subroutine output_specie_potentials

  ! --------------------------------------------------
  ! Output a table of the salt chemical potentials
  subroutine output_salt_potentials(thefid)
    use strngs
    implicit none
    integer, intent(in) :: thefid
    integer :: igc
    character(15) :: fltout1, fltout2
    if (dbc) then
      if (dbc_level.ge.dbc_require) then
        call require_writable(thefid)
      endif
    endif
    ! -----salt chemical potentials-----------------------------
    write(unit=thefid,fmt='(1X,A6,2(1X,A15))')"SALT","TARGET CONC.","CHEM. POT."
    k6z8l3: do igc=1,nsalt_
      call str (salts(igc)%target_concentration_, fltout1)
      call str (chemps(igc), fltout2)
      write(unit=thefid,fmt='(1X,A6,2(1X,A15))')salts(igc)%code_name_,fltout1,fltout2
    enddo k6z8l3
  end subroutine output_salt_potentials

  ! ----------------------------------------------------
  ! TRIAL RATES
  ! ----------------------------------------------------

  
  ! ----------------------------------------------------
  ! Rate that a specie is selected for a move compared
  ! to other species. 
  double precision function ratspc(ispec)
    implicit none
    integer, intent(in) :: ispec
    if (dbc) then
      if (dbc_level.ge.dbc_index) then
        if (ispec.lt.1.or.ispec.gt.nspec_) stop "Error: Specie index out of range"
      endif
    endif
    ratspc=species(ispec)%rate_specie_
  end function ratspc

  ! ----------------------------------------------------
  ! Rate that a specie undergoes a simple move (jump or 
  ! displacement) or a enhanced move (jump-in or jump-out).
  ! This is used once a specie has been selected for a move
  double precision function ratexc(ispec)
    implicit none
    integer, intent(in) :: ispec
    if (dbc) then
      if (dbc_level.ge.dbc_index) then
        if (ispec.lt.1.or.ispec.gt.nspec_) stop "Error: Specie index out of range"
      endif
    endif
    ratexc=species(ispec)%rate_exchange_
  end function ratexc

  ! ----------------------------------------------------
  ! Rate that a specie undergoes a jump versus small
  ! displacement move. This is used once a specie has
  ! been selected for this type of move.
  double precision function ratmov(ispec)
    implicit none
    integer, intent(in) :: ispec
    if (dbc) then
      if (dbc_level.ge.dbc_index) then
        if (ispec.lt.1.or.ispec.gt.nspec_) stop "Error: Specie index out of range"
      endif
    endif
    ratmov=species(ispec)%rate_move_
  end function ratmov

  ! ----------------------------------------------------
  ! Rate that a specie undergoes an individual GC move
  ! compared to other species.
  double precision function ratgri(ispec)
    implicit none
    integer, intent(in) :: ispec
    if (dbc) then
      if (dbc_level.ge.dbc_index) then
        if (ispec.lt.1.or.ispec.gt.nspec_) stop "Error: Specie index out of range"
      endif
    endif
    ratgri=species(ispec)%rate_change_
  end function ratgri

  ! --------------------------------------------------
  ! Rate that a salt undergoes a GC add/delete move
  ! compared to other salts
  double precision function ratgrs(ispec)
    implicit none
    integer, intent(in) :: ispec
    if (dbc) then
      if (dbc_level.ge.dbc_index) then
        if (ispec.lt.1.or.ispec.gt.nsalt_) stop "Error: Salt index out of range"
      endif
    endif
    ratgrs=salts(ispec)%rate_change_
  end function ratgrs

  ! ----------------------------------------------------
  ! Rate that a particular region is selected for
  ! add/delete GC move.
  double precision function ratreg(irgn,ispec)
    implicit none
    integer, intent(in) :: irgn,ispec
    if (dbc) then
      if (dbc_level.ge.dbc_index) then
        if (irgn.lt.1.or.irgn.gt.ibulk) stop "Error: Region index out of range"
        if (ispec.lt.1.or.ispec.gt.nspec_) stop "Error: Specie index out of range"
      endif
    endif
    ratreg=species(ispec)%rate_region_(irgn)
  end function ratreg

  ! -----------------------------------------------------
  ! Read salt information section for a salt
  !
  ! @param fid : input unit number
  ! @param sname : the name value that caused this function to be called
  ! @param svalue : the value associated with the name (may be mepty string)
  !
  ! @pre sname=fssalt
  !
  subroutine rdsalt(fid,sname,svalue,istat)
    use strngs
    implicit none
    integer, intent(in) :: fid
    character(len=*), intent(in) :: sname,svalue
    integer, intent(out) :: istat
    logical, dimension(3) :: mask_
    character(32) :: nme_
    character(1024) :: val_
    if (dbc) then
      if (dbc_level.ge.dbc_require) then
        call require_readable(fid)
        if (sname.ne.fssalt) stop "Error: incorrect section name"
        if (len_trim(svalue).ne.0) stop "Error: section does not any parameters"
      endif
    endif
    mask_=.false.

    if (.not.allocated(salts)) then
      allocate(salts(nsltmx))
      nsalt_=0
    endif

    nsalt_=nsalt_+1
    call ctor_salt(salts(nsalt_))
    ! process section
    e2a4b7: do
      val_ = " "
      nme_ = " "
      call readnv(fid,nme_,val_,istat)
      ! exit on bad read or section 'end'
      if (istat.ne.0) return
      if (nme_.eq.fsend) exit e2a4b7
      ! looking for salts(name or cation)%code_name_,chex,ctarg,ratgrs
      w5k6s5: select case (nme_)
      case (fsislt) w5k6s5
        read(val_,*)salts(nsalt_)%code_name_
        salts(nsalt_)%code_name_(3:4)='Cl'
        mask_(1)=.true.
      case (fsname) w5k6s5
        read(val_,*)salts(nsalt_)%code_name_
        mask_(1)=.true.
      case (fsrtgr) w5k6s5
        read(val_,'(D20.13)')salts(nsalt_)%rate_change_
        mask_(2)=.true.
      case (fsctrg) w5k6s5
        read(val_,'(D20.13)')salts(nsalt_)%target_concentration_
        mask_(3)=.true.
      case default w5k6s5
        call d0s2r1("Name "//nme_//" is not valid in salt section")
      end select w5k6s5
    enddo e2a4b7

    d0q8g6: if (.not.all(mask_)) then
      call d0s2r1("Not all required tags were present.")
    endif d0q8g6

  contains

    subroutine d0s2r1(msg)
      use strngs
      implicit none
      character(len=*), intent(in) :: msg

      write(unit=fidlog,fmt=*)"## Bad salt section in input ##"
      write(unit=fidlog,fmt=*)msg
      write(unit=fidlog,fmt=*)"## Bad salt section in input ##"

      call help_salt(fidlog)
      stop 1
    end subroutine d0s2r1

  end subroutine rdsalt

  ! -----------------------------------------------------
  ! Read specie information section for one specie
  !
  ! A input section for a specie may contain:
  !
  !   specie (free|chonly|flexible|mobile)
  !       # alternate method of specifying specie type
  !     type  (free|chonly|flexible|mobile)
  !     z {REAL} # particle charge
  !     d {REAL} # particle diameter
  !       # Label for specie, must correspond to that used in the 
  !       # names of salt
  !     name "{2 character text}"
  !     chex {REAL} # Precomputed chemical potential
  !       # Parameters controlling a particle of this
  !       # specie's involvement in the MC move.
  !     ratspc {REAL}
  !     ratmov {REAL}
  !     ratexc {REAL}
  !     ratgr {REAL}
  !     ratreg {REAL} {REAL} {REAL} {REAL}
  !       # Particle locations - (only) optional for free type
  !       # species. The 'n' option must be followed by position
  !       # information for the given number of particles.
  !     n {INTEGER}
  !       # For 'free' and 'chonly' species require
  !       # only x, y, z coordinates for their starting position.
  !       #  X      Y      Z
  !       {REAL} {REAL} {REAL}
  !
  !       # For 'mobile' and 'flexible' species you must specify the starting
  !       # position as above and additionally specify
  !       # the allowed update radii.  This can be followed
  !       # by an additional set of x, y, z coordinates defining
  !       # the centre-point for the radius. If only one set
  !       # of coordinates is specified then they are used for
  !       # the starting and centre-points.  Note that there
  !       # is no restriction on using only 4 or only
  !       # 7 parameter definitions per input section.
  !       #
  !       #  X      Y      Z   Upd. R   X0     Y0     Z0 
  !       {REAL} {REAL} {REAL} {REAL}
  !       {REAL} {REAL} {REAL} {REAL} {REAL} {REAL} {REAL}
  !
  !     end
  !
  ! @param fid : input unit number
  ! @param sname : the name value that caused this function to be
  !   called
  ! @param svalue : the value associated with the name (may be
  !   empty string)
  !
  ! @pre sname=fsspec
  subroutine rdspec(fid,sname,svalue,istat)
    use strngs
    implicit none
    integer, intent(in) :: fid
    character(len=*), intent(in) :: sname,svalue
    integer, intent(out) :: istat
    logical, dimension(10) :: mask_
    character(32) :: nme_
    character(1024) :: val_
    character(3) :: type_
    integer, parameter :: notype=4
    integer, save :: allocsz=nspcmx
    type (specie), dimension(:), allocatable :: tmp_species
    integer :: idx, stat

    if (dbc) then
      if (dbc_level.ge.dbc_require) then
        call require_readable(fid)
        if (sname.ne.fsspec) stop "Error: incorrect section name"
      endif
    endif
    istat = 0
    mask_=.false.

    if (.not.allocated(species)) then
      allocate(species(allocsz), STAT=stat)
      if (stat.ne.0) stop "Memory allocation of SPECIES failed"
      nspec_=0
      imax_mob_=0
      imax_flex_=0
      nstr_=0
      idxcl_=0
    endif
    ! if more species than current allocated size, increase 
    ! allocation.
    ! ** CURRENTLY NOT POSSIBLE ** as nspcmx is assumed
    ! maximum specie number in many other locations
    if (nspec_.eq.allocsz) then
      stop "Maximum number of species is 16"
      allocate(tmp_species(allocsz), STAT=stat)
      if (stat.ne.0) stop "Memory allocation of TMPSPECIES failed"
      do idx=1,nspec_
        call move_specie(tmp_species(idx),species(idx))
      enddo
      deallocate(species, STAT=stat)
      if (stat.ne.0) stop "Memory deallocation of SPECIES failed"
      allocsz=allocsz+nspcmx
      allocate(species(allocsz), STAT=stat)
      if (stat.ne.0) stop "Memory allocation of SPECIES failed"
      do idx=1,nspec_
        call move_specie(species(idx),tmp_species(idx))
      enddo
      deallocate(tmp_species, STAT=stat)
      if (stat.ne.0) stop "Memory deallocation of TMPSPECIES failed"
    endif
    ! increment counter
    nspec_=nspec_+1

    ! initialise specie object
    call ctor__(species(nspec_))
    species(nspec_)%type_=notype ! set to an invalid value
    ! optional argument to specie can be specie type
    select case (svalue(1:3))
    case (fsmobl(1:3))
      species(nspec_)%type_=nmob_l
      mask_(6)=.true.
    case (fsflxd(1:3))
      species(nspec_)%type_=nflx_l
      mask_(6)=.true.
    case (fschon(1:3))
      species(nspec_)%type_=ncho_l
      mask_(6)=.true.
    case (fsfree(1:3))
      species(nspec_)%type_=nfre_l
      mask_(6)=.true.
    case default
    end select
    ! process section
    j1f0t7: do
      val_ = " "
      nme_ = " "
      call readnv(fid,nme_,val_,istat)
      ! exit on bad read or section 'end'
      if (istat.ne.0) return
      if (nme_.eq.fsend) exit j1f0t7
      ! looking for xri,xz,chemp,chex,ctarg,ntarg,type_,(nispc,xyzr)
      d9i3s0: select case (nme_)
      case (fsz) d9i3s0
        read(val_,'(D20.13)')species(nspec_)%valency_
        mask_(1)=.true.
      case (fsd) d9i3s0
        read(val_,'(D20.13)')species(nspec_)%radius_
        species(nspec_)%radius_=species(nspec_)%radius_/2
        mask_(2)=.true.
      case (fschex) d9i3s0
        read(val_,'(D20.13)')species(nspec_)%chem_excess_
        mask_(3)=.true.
      case (fsn) d9i3s0
        if (associated(species(nspec_)%xyz_)) &
             & call c8c5p4(""""//fsn//""" tag appears more than once for a specie")
        read(val_,*)species(nspec_)%input_count_
        y5s6w0: select case (species(nspec_)%type_)
        case (nmob_l, nflx_l) y5s6w0
          ! required information is x,y,z,update radii
          allocate(species(nspec_)%xyz_(3,species(nspec_)%input_count_))
          species(nspec_)%xyz_ = 0
          allocate(species(nspec_)%sxyzr_(4,species(nspec_)%input_count_))
          species(nspec_)%sxyzr_ = 0
        case (ncho_l, nfre_l) y5s6w0
          ! required information is x,y,z
          allocate(species(nspec_)%xyz_(3,species(nspec_)%input_count_))
          species(nspec_)%xyz_ = 0
          nullify(species(nspec_)%sxyzr_)
        case default y5s6w0
          ! assume is mobile type, but do not set mask(6) here.
          species(nspec_)%type_=nmob_l
          allocate(species(nspec_)%xyz_(3,species(nspec_)%input_count_))
          species(nspec_)%xyz_ = 0
          allocate(species(nspec_)%sxyzr_(4,species(nspec_)%input_count_))
          species(nspec_)%sxyzr_ = 0
        end select y5s6w0
        call n5v2h7
      case (fsname) d9i3s0
        read(val_,*)species(nspec_)%code_name_
        mask_(4)=.true.
      case (fstype) d9i3s0
        read(val_,*)type_
        j4b6e8: select case (type_)
        case (fsmobl(1:3)) j4b6e8
          species(nspec_)%type_=nmob_l
        case (fsflxd(1:3)) j4b6e8
          species(nspec_)%type_=nflx_l
        case (fschon(1:3)) j4b6e8
          species(nspec_)%type_=ncho_l
        case (fsfree(1:3)) j4b6e8
          if (associated(species(nspec_)%xyz_)) &
               & call c8c5p4(""""//fsn//""" tag in """//fsfree//""" type specie is an error")
          species(nspec_)%type_=nfre_l
        case default j4b6e8
          call c8c5p4("Unknown specie type "//val_)
        end select j4b6e8
        mask_(5)=.true.
      case (fsrtsp) d9i3s0
        read(val_,'(D20.13)')species(nspec_)%rate_specie_
        mask_(6)=.true.
      case (fsrtmv) d9i3s0
        read(val_,'(D20.13)')species(nspec_)%rate_move_
        mask_(7)=.true.
      case (fsrtex) d9i3s0
        read(val_,'(D20.13)')species(nspec_)%rate_exchange_
        mask_(8)=.true.
      case (fsrtgr) d9i3s0
        read(val_,'(D20.13)')species(nspec_)%rate_change_
        mask_(9)=.true.
      case (fsrtrg) d9i3s0
        call redflt(val_,species(nspec_)%rate_region_,nrgnmx,idx,istat)
        if (idx.ne.4.or.istat.ne.0) then
          call c8c5p4("Error reading four numbers from line """//fsrtrg//" "//trim(val_)//"""")
        endif
        mask_(10)=.true.
      case default d9i3s0
        call c8c5p4("Name "//nme_//" is not valid in specie section")
      end select d9i3s0
    enddo j1f0t7
    if (species(nspec_)%type_.eq.notype) then
      if (associated(species(nspec_)%xyz_)) then
        call c8c5p4("Using """//fsn//""" without specifying """//fstype//""" is an error.")
      else
        ! Assume is 'free' type
        species(nspec_)%type_ = nfre_l
      endif
    endif

    select case (species(nspec_)%type_)
    case (nmob_l, nflx_l, ncho_l) 
      if (dfeq(species(nspec_)%rate_specie_,0.D0)) species(nspec_)%rate_specie_=0.2D0
      species(nspec_)%rate_exchange_=0.D0
      if (dfeq(species(nspec_)%rate_move_,0.D0)) species(nspec_)%rate_move_=1.D0
      species(nspec_)%rate_change_=0.D0
      species(nspec_)%rate_region_=0.D0
      mask_(6)=.true.
      mask_(7)=.true.
      mask_(8)=.true.
      mask_(9)=.true.
      mask_(10)=.true.
    case default 
    end select

    d0q8g6: if (.not.all(mask_)) then
      if (debug) write (*,*)"! ",mask_
      ! allow ratreg to be missing iff ratgri is missing
      if (.not.all(mask_(1:9)).or.mask_(9)) call c8c5p4("Not all required tags were present.")
    endif d0q8g6

    call sortit

  contains

    ! Subroutine for mobile/flexible ions that can have either 4 or 7 elements per line.
    ! * 4 elements : 'x y z R' where centre point and start point are identical
    ! * 7 elements : ' x y z R X Y Z' where xyz is start point and XYZ is centre
    subroutine n5v2h7
      implicit none
      double precision, dimension(7) :: xyz_data_
      integer :: idx, jdx, istat
      !integer, parameter :: linemax=256
      !character(linemax) :: line_
      do idx=1,species(nspec_)%input_count_
        ! read line and strip comments.
        do
          istat=0
          read(unit=fid,fmt='(A)',iostat=istat) val_
          if (istat.ne.0) then
            call c8c5p4("Ion position specification error reading input file")
          endif
          val_=decmnt(val_)
          if (len_trim(val_).ne.0) exit
        end do
        call redflt(val_,xyz_data_,7,jdx,istat)
        select case (jdx)
        case (3)
          if (species(nspec_)%type_.eq.nmob_l.or.species(nspec_)%type_.eq.nflx_l) then
            call c8c5p4("Mobile/flexible ion specification error reading """//trim(val_)//"""")
          endif
          species(nspec_)%xyz_(1:3,idx) = xyz_data_(1:3)
        case (4)
          if (species(nspec_)%type_.eq.nmob_l.or.species(nspec_)%type_.eq.nflx_l) then
            species(nspec_)%xyz_(1:3,idx) = xyz_data_(1:3)
            species(nspec_)%sxyzr_(2:4,idx) = xyz_data_(1:3)
            species(nspec_)%sxyzr_(1,idx) = sqr(xyz_data_(4))
          else
            call c8c5p4("Solute/channel-only ion specification error reading """//trim(val_)//"""")
          endif
        case (7)
          if (species(nspec_)%type_.eq.nmob_l.or.species(nspec_)%type_.eq.nflx_l) then
            species(nspec_)%xyz_(1:3,idx) = xyz_data_(1:3)
            species(nspec_)%sxyzr_(2:4,idx) = xyz_data_(5:7)
            species(nspec_)%sxyzr_(1,idx) = sqr(xyz_data_(4))
          else
            call c8c5p4("Solute/channel-only ion specification error reading """//trim(val_)//"""")
          endif
        case default
          call c8c5p4("Ion specification error reading """//trim(val_)//"""")
        end select
      enddo
    end subroutine n5v2h7

    subroutine c8c5p4(msg)
      use strngs
      implicit none
      character(len=*), intent(in) :: msg
      write(unit=fidlog,fmt=*)
      write(unit=fidlog,fmt=*)"## Bad specie section in input ##"
      write(unit=fidlog,fmt=*)msg
      write(unit=fidlog,fmt=*)"## Bad specie section in input ##"
      call help_specie(fidlog)
      stop 1
    end subroutine c8c5p4

    subroutine sortit
      implicit none
      integer :: insert ! insert position
      integer :: idx ! loop index
      !  We sort species in order imax_mob_, imax_flex_, ncho, idxcl_, nfre
      !  nchonly ions appear before free ions
      insert=0
      select case (species(nspec_)%type_)
      case (nmob_l)
        imax_mob_=imax_mob_+1
        imax_flex_=imax_flex_+1
        nstr_=nstr_+1
        idxcl_=idxcl_+1
        insert=imax_mob_
      case (nflx_l)
        imax_flex_=imax_flex_+1
        nstr_=nstr_+1
        idxcl_=idxcl_+1
        insert=imax_flex_
      case (ncho_l)
        nstr_=nstr_+1
        idxcl_=idxcl_+1
        insert=nstr_
      case default
        if (species(nspec_)%code_name_.eq.'Cl') then
          idxcl_=idxcl_+1
          insert=idxcl_
        else
          ! free ion is already in its final position
          return
        endif
      end select

      if (nspec_.eq.1) return
      ! here we know where to insert specie into species
      do idx=nspec_,insert + 1,-1
        call swapspec(species(idx-1),species(idx))
      enddo
    end subroutine sortit

  end subroutine rdspec

  ! -----------------------------------------------------
  ! Read subspecie information section for a subspecie
  !
  ! @param fid : input unit number
  ! @param sname : the name value that caused this function to be called
  ! @param svalue : the value associated with the name (may be mepty string)
  !
  ! @pre sname=fssubs
  !
  subroutine rdsubs(fid,sname,svalue,istat)
    use strngs
    implicit none
    integer, intent(in) :: fid
    character(len=*), intent(in) :: sname,svalue
    integer, intent(out) :: istat
    logical, dimension(5) :: mask_
    character(32) :: nme_
    character(1024) :: val_
    integer :: subcount_
    if (dbc) then
      if (dbc_level.ge.dbc_require) then
        call require_readable(fid)
        if (sname.ne.fssubs) stop "Error: incorrect section name"
        if (len_trim(svalue).ne.0) stop "Error: section does not any parameters"
      endif
    endif
    mask_=.false.

    if (.not.allocated(subspecies)) then
      allocate(subspecies(nsltmx))
      nsubs_=0
    endif
    subcount_ = 0
    nsubs_=nsubs_+1
    call ctor_subs(subspecies(nsubs_))
    ! process section
    e2a4b7: do
      val_ = " "
      nme_ = " "
      call readnv(fid,nme_,val_,istat)
      ! exit on bad read or section 'end'
      if (istat.ne.0) return
      if (nme_.eq.fsend) exit e2a4b7
      ! looking for subss(name or cation)%code_name_,chex,ctarg,ratgrs
      w5k6s5: select case (nme_)
      case (fsname) w5k6s5
        read(val_,*)subspecies(nsubs_)%code_name_
        mask_(1)=.true.
      case (fsrtsw) w5k6s5
        read(val_,'(D20.13)')subspecies(nsubs_)%rate_swap_
        mask_(2)=.true.
      case (fsenth) w5k6s5
        read(val_,'(D20.13)')subspecies(nsubs_)%enthalpy_
        mask_(3)=.true.
      case (fsentr) w5k6s5
        read(val_,'(D20.13)')subspecies(nsubs_)%entropy_
        mask_(4)=.true.
      case (fsgrnd) w5k6s5
        read(val_,*)subspecies(nsubs_)%sub_names_(1)
        subcount_ = subcount_ + 1
        if (subcount_.eq.2) mask_(5)=.true.
      case (fsexct) w5k6s5
        read(val_,*)subspecies(nsubs_)%sub_names_(2)
        subcount_ = subcount_ + 1
        if (subcount_.eq.2) mask_(5)=.true.
      case default w5k6s5
        call d0s21r("Name "//nme_//" is not valid in sub-specie section")
      end select w5k6s5
    enddo e2a4b7

    d0q8g6: if (.not.all(mask_)) then
      call d0s21r("Not all required tags were present.")
    endif d0q8g6

  contains

    subroutine d0s21r(msg)
      use strngs
      implicit none
      character(len=*), intent(in) :: msg

      write(unit=fidlog,fmt=*)"## Bad subspecie section in input: ##"
      write(unit=fidlog,fmt=*)msg
      write(unit=fidlog,fmt=*)"## Bad subspecie section in input: ##"
      call help_subspecie(fidlog)
      stop 1
    end subroutine d0s21r

  end subroutine rdsubs

  ! ------------------------------------------------------------
  ! Convert the species data into standard specie form.
  !
  ! @codenotes
  ! actions
  !  * write data from rdtmp into std so that imax_mob_, imax_flex_ and 
  !  nchonly ions appear before free ions
  !  * set indxcl
  !  * map salt ions to their cation indices
  !  * calculate derived data (q,chempi)
  ! @endcodenotes
  subroutine rfspec(rl1, zlim)
    use strngs
    implicit none
    double precision, intent(in) :: rl1, zlim
    ! LOCALS
    integer :: iv    ! salt valency
    integer :: igc ! salt loop index
    integer :: ispec ! specie loop index
    character(*), parameter, dimension(5) :: ftyp=(/ "unkn", " mob", "flex", "chnl", "free" /)
    integer :: styp
    character(20) :: fltout

    ! map subspecie ions to their specie indices (species are
    ! sorted into their final positions during rdspec)
    if (nsubs_.ge.1) then
      f8k92y: do igc=1,nsubs_
        ! Find specie from first two letters of fsubs!
        subspecies(igc)%sub_indices_(1) = private_spcidx(subspecies(igc)%sub_names_(1))
        subspecies(igc)%sub_indices_(2) = private_spcidx(subspecies(igc)%sub_names_(2))
        if (subspecies(igc)%sub_indices_(1).eq.subspecies(igc)%sub_indices_(2)) then
          write(unit=fidlog,fmt='(1X,"Both subnames ",A2," and ",A2," map to the same index")') &
               subspecies(igc)%sub_names_(1),subspecies(igc)%sub_names_(2)
          stop "Error: Bad subspecie specification"
        endif
        a6r04k: if (subspecies(igc)%sub_indices_(1).eq.0) then
          write(unit=fidlog,fmt='(1X,"Could not find sub-specie",A2," for ",A2)')subspecies(igc)%sub_names_(1) &
               & ,subspecies(igc)%code_name_
          stop "Error: Bad subspecie specification"
        endif a6r04k
        if (species(subspecies(igc)%sub_indices_(1))%type_.ne.nfre_l) then
          write(unit=fidlog,fmt='(1X,"Sub specie",A2," for ",A2,"is not a solute ion")')subspecies(igc)%sub_names_(1) &
               & ,subspecies(igc)%code_name_
          stop "Error: Invalid non-solute ion in subspecie specification"
        endif
        a6r40k: if (subspecies(igc)%sub_indices_(2).eq.0) then
          write(unit=fidlog,fmt='(1X,"Could not find cation ",A2," for ",A2)')subspecies(igc)%sub_names_(2) &
               & ,subspecies(igc)%code_name_
          stop "Error: Bad subs specification"
        endif a6r40k
        if (species(subspecies(igc)%sub_indices_(2))%type_.ne.nfre_l) then
          write(unit=fidlog,fmt='(1X,"Sub specie ",A2," for ",A2," is not a solute ion")')subspecies(igc)%sub_names_(1) &
               & ,subspecies(igc)%code_name_
          stop "Error: Invalid non-solute ion in subspecie specification"
        endif

        if (subspecies(igc)%sub_indices_(1).eq.subspecies(igc)%sub_indices_(2)) then
          write(unit=fidlog,fmt='(1X,"The same specie ",A2," appears twice in ",A2)')subspecies(igc)%sub_names_(1) &
               & ,subspecies(igc)%code_name_
          stop "Error: Invalid use of a specie twice in one subspecie specification"
        endif
      enddo f8k92y
      ! Check that each specie only appears once
      if (nsubs_.ge.2) then
        do igc=1,nsubs_-1
          do ispec=igc+1,nsubs_
            if (subspecies(igc)%sub_indices_(1).eq.subspecies(ispec)%sub_indices_(1)) then
              write(unit=fidlog,fmt='(1X,"The same specie ",A2," appears in both ",A2," and ",A2," subspecies")') &
                   & subspecies(igc)%sub_names_(1), subspecies(igc)%code_name_, subspecies(ispec)%code_name_
              stop "Error: Invalid use of a specie in two subspecie specifications"
            endif
            if (subspecies(igc)%sub_indices_(1).eq.subspecies(ispec)%sub_indices_(2)) then
              write(unit=fidlog,fmt='(1X,"The same specie ",A2," appears in both ",A2," and ",A2," subspecies")') &
                   & subspecies(igc)%sub_names_(1), subspecies(igc)%code_name_, subspecies(ispec)%code_name_
              stop "Error: Invalid use of a specie in two subspecie specifications"
            endif
            if (subspecies(igc)%sub_indices_(2).eq.subspecies(ispec)%sub_indices_(1)) then
              write(unit=fidlog,fmt='(1X,"The same specie ",A2," appears in both ",A2," and ",A2," subspecies")') &
                   & subspecies(igc)%sub_names_(2), subspecies(igc)%code_name_, subspecies(ispec)%code_name_
              stop "Error: Invalid use of a specie in two subspecie specifications"
            endif
            if (subspecies(igc)%sub_indices_(2).eq.subspecies(ispec)%sub_indices_(2)) then
              write(unit=fidlog,fmt='(1X,"The same specie ",A2," appears in both ",A2," and ",A2," subspecies")') &
                   & subspecies(igc)%sub_names_(2), subspecies(igc)%code_name_, subspecies(ispec)%code_name_
              stop "Error: Invalid use of a specie in two subspecie specifications"
            endif
          enddo
        enddo
      endif
    endif

    ! --- salt specific parameters -------------------------
    ! map salt ions to their cation indices (species are
    ! sorted into their final positions during rdspec)
    f8k9y2: do igc=1,nsalt_
      ! Find specie from first two letters of fsalt!
      salts(igc)%cation_index_ = spcidx(salts(igc)%code_name_(1:2))
      a6r0k4: if (salts(igc)%cation_index_.eq.0) then
        write(unit=fidlog,fmt='(1X,"Could not find cation ",A2," for ",A4)')salts(igc)%code_name_ &
             & ,salts(igc)%code_name_
        stop "Error: No cation matches salt specification"
      endif a6r0k4
      ! Check this cation index is only used in one slat definition
      if (igc.gt.1) then
        do ispec=1,igc-1
          if (salts(igc)%cation_index_.eq.salts(ispec)%cation_index_) then
            write(unit=fidlog,fmt='(1X,"Error: cation ",A2," appears in two salt: ",A4," ",A4)')&
                 salts(igc)%code_name_,salts(igc)%code_name_,salts(ispec)%code_name_
            stop "Error: Cation in two salt specifications"
          endif
        enddo
      endif
    enddo f8k9y2

    ! -- generate q values --
    m0u9b7: do ispec=1,nspec_
      if (dfeq(species(ispec)%valency_,0.D0)) then
        write(unit=fidlog,fmt=*) "Warning: specie without charge is experimental."
        species(ispec)%red_charge_ = 0.D0
      else
        species(ispec)%red_charge_ = qstar() * species(ispec)%valency_
      endif
    enddo m0u9b7

    write(unit=fidlog,fmt='(72("-"))')
    write(unit=fidlog,fmt=*)"Specie type and salt data summary"
    write(unit=fidlog,fmt='(72("-"))')


    write(unit=fidlog,fmt=*)" Number of structural ion species = ", nstr_
    write(unit=fidlog,fmt='(52("-"))')
    write(unit=fidlog,fmt=*)"       mobile within sphere (mob) = ", imax_mob_
    write(unit=fidlog,fmt=*)"     mobile within channel (chnl) = ", nstr_-imax_flex_
    write(unit=fidlog,fmt=*)"                  flexible (flex) = ", imax_flex_-imax_mob_
    write(unit=fidlog,fmt='(52("-"))')
    write(unit=fidlog,fmt='(12x,"  xq [e]  ","  d [A]  "," type")')
    write(unit=fidlog,fmt='(52("-"))')
    o7a4o6: do ispec=1,nspec_
      styp=1 ! default to unknown
      if (mobile(ispec)) then
        styp=2
      else if (flexible(ispec)) then
        styp=3
      else if (chonly(ispec)) then
        styp=4
      else if (isfree(ispec)) then
        styp=5
      endif
      write(unit=fidlog,fmt='("  ion   ",A,5X,F5.2,1X,F8.4,3X,A)')species(ispec)%code_name_, &
           & species(ispec)%valency_, species(ispec)%radius_*2,ftyp(styp)
    enddo o7a4o6
    write(unit=fidlog,fmt='(52("-"))')

    ! ---- DEFINITION OF VARIABLES FOR INDIVIDUAL GC -----------

    ! - calculate individual ion concentrations from salts
    species(1:nspec_)%target_concentration_ = 0.D0
    d5j5c4: do igc=1,nsalt_
      ispec = salts(igc)%cation_index_ ! cation index
      iv=nint(species(ispec)%valency_)
      species(ispec)%target_concentration_ = species(ispec)%target_concentration_ + &
           & salts(igc)%target_concentration_
      species(idxcl_)%target_concentration_ = species(idxcl_)%target_concentration_ + &
           & iv * salts(igc)%target_concentration_
    enddo d5j5c4

    if (nsubs_.ne.0) then
      ! Massage concentrations for subspecies.
      call subs_init
    endif

    ! Output the finalised results
    write(unit=fidlog,fmt='(72("-"))')
    call output_specie_potentials(fidlog)
    write(unit=fidlog,fmt='(72("-"))')
    call output_salt_potentials(fidlog)
    write(unit=fidlog,fmt='(72("-"))')

    ! normalise rates.
    call normrate

    ! rescale center points.
    ! call rescale

    write(unit=fidlog,fmt='(72("-"))')
    call str(ionstr(), fltout)
    write(unit=fidlog,fmt=*)"TOTAL IONIC STRENGTH: ",trim(adjustl(fltout))
    write(unit=fidlog,fmt='(72("-"))')
    write(unit=fidlog,fmt=*)"Interpreted specie, salt and subspecie sections of input"
    write(unit=fidlog,fmt='(72("-"))')
    call echoin(fidlog)
    call ecsalt(fidlog)
    call ecsubs(fidlog)
    write(unit=fidlog,fmt='(72("-"))')

  contains

    ! ---- NORMALIZE TRIAL RATES -------------------------------
    subroutine normrate
      implicit none
      double precision :: ratsum    ! sum of rates (for normalisaton)
      integer :: ireg ! region loop index
      integer :: igc ! salt loop index
      integer :: ispec ! specie loop index

      write(unit=fidlog,fmt='(70("-"))')
      write(unit=fidlog,fmt=*)"Grand canonical update rates"
      write(unit=fidlog,fmt='(70("-"))')
      write(unit=fidlog,fmt='(1X,A3,2X,7("|",1X,A6,1X))')"ION","MOVE","EXCG","ADD","REG 1","REG 2","REG 3","REG 4"
      ! ---- normalize ratreg ------------------------------------
      m1z6t1: do ispec=idxcl_,nspec_
        ratsum = sum(species(ispec)%rate_region_)
        if (dfeq(ratsum,0.D0)) then
          write(unit=fidlog,fmt=*) "Error: sum of region add/delete rates for ",species(ispec)%code_name_," is zero"
          stop 1
        endif
        species(ispec)%rate_region_ = species(ispec)%rate_region_ / ratsum
      enddo m1z6t1

      ! ---- normalize ratspc ------------------------------------
      ratsum=sum(species(1:nspec_)%rate_specie_)
      if (dfeq(ratsum,0.D0)) then
        write(unit=fidlog,fmt=*) "Warning: sum of specie rates for moves (",fsrtsp,") is zero"
      else
        species(1:nspec_)%rate_specie_ =species(1:nspec_)%rate_specie_ / ratsum
      endif

      ! ---- normalize ratgri ------------------------------------
      ratsum=sum(species(idxcl_:nspec_)%rate_change_)
      if (dfeq(ratsum,0.D0)) then
        write(unit=fidlog,fmt=*) "Warning: sum of rates for individual add/delete (",fsrtgr,") is zero"
      else
        species(idxcl_:nspec_)%rate_change_ = species(idxcl_:nspec_)%rate_change_ / ratsum
      endif
      k6g5g1: do ispec=1,nspec_
        if (ispec.le.nstr_) then
          write(unit=fidlog,fmt='(1X,A2,3X,("|",1X,F7.4))')species(ispec)%code_name_, &
               & species(ispec)%rate_specie_
        else
          write(unit=fidlog,fmt='(1X,A2,3X,7("|",1X,F7.4))')species(ispec)%code_name_, &
               & species(ispec)%rate_specie_, species(ispec)%rate_exchange_, species(ispec)%rate_change_, &
               & (species(ispec)%rate_region_(ireg),ireg=1,nrgnmx)
        endif
      enddo k6g5g1

      ! ---- normalize salt ratgrs ------------------------------------
      ratsum=sum(salts(1:nsalt_)%rate_change_)
      if (dfeq(ratsum,0.D0)) then
        write(unit=fidlog,fmt=*) "Warning: sum of rates for salt add/delete (",fsrtgr,") is zero"
      else
        salts(1:nsalt_)%rate_change_ = salts(1:nsalt_)%rate_change_ / ratsum
      endif
      write(unit=fidlog,fmt='(1X,A4,1X,"|",1X,A7)')"SALT","CHANGE"
      k8p4f9: do igc=1,nsalt_
        write(unit=fidlog,fmt='(1X,A4,1X,("|",1X,F7.4))')salts(igc)%code_name_, salts(igc)%rate_change_
      enddo k8p4f9

      ! ---- normalize subspecie ratgrs ------------------------------------
      if (nsubs_.ge.1) then
        ratsum=sum(subspecies(1:nsubs_)%rate_swap_)
        if (dfeq(ratsum,0.D0)) then
          write(unit=fidlog,fmt=*) "Warning: sum of rates for subspecie swaps (",fsrtgr,") is zero"
        else
          subspecies(1:nsubs_)%rate_swap_ = subspecies(1:nsubs_)%rate_swap_ / ratsum
        endif
        write(unit=fidlog,fmt='(1X,A4,1X,"|",1X,A7)')"SUBS","CHANGE"
        k8p49f: do igc=1,nsubs_
          write(unit=fidlog,fmt='(1X,A2,"-",A2,("|",1X,F7.4))')subspecies(igc)%sub_names_(1), &
               & subspecies(igc)%sub_names_(2), subspecies(igc)%rate_swap_
        enddo k8p49f
      endif
    end subroutine normrate

    ! ---------------------------------------------------------------
    ! Scale localised ions centroid positions rsx, rsy and rsz to fit
    ! into the geometry.
    !
    subroutine rescale
      implicit none
      double precision :: rzmin, rzmax, rrmxsq, rtmp, rscale, zoffset
      integer :: ispec, jspec, idx
      !  -----------------------------------------------
      !  Handle MOBILE type structural ion centre-points
      !  -----------------------------------------------
      ! rescale mobile centre points
      rzmin =  zlim
      rzmax = -zlim
      ! Use jspec and rrmxsq to cache minimum radial
      jspec = 0
      rrmxsq = 0
      !! RESCALE RSX and RSY
      do ispec=1,imax_mob_
        ! We do not rescale 'flexible' ions in rsx, rsy as they are allowed anywhere 
        ! a free ion is permitted. We do scale them in rsz, but ignore them when
        ! calculating the z scale factor.
        !! We use r1**2 - radius**2 as the maximum squared radius purposefully instead of
        !! (r1 - radius)**2.
        !! rrmax is position inside cylinder where an arc of the length
        !! of the sphere diameter intersects with a radial at right angles
        !! ie rrmax ^ 2 + r_sphere ^ 2 == r_cyl ^ 2 { rrmxsq = rrmax ^ 2}
        rrmxsq = sqr(rl1)-sqr(species(ispec)%radius_)
        do idx=1,species(ispec)%input_count_
          rzmin = min(rzmin,species(ispec)%sxyzr_(4, idx) - species(ispec)%radius_)
          rzmax = max(rzmax,species(ispec)%sxyzr_(4, idx) + species(ispec)%radius_)
          rtmp = sqr(species(ispec)%sxyzr_(2, idx)) + sqr(species(ispec)%sxyzr_(3, idx))
          if (rrmxsq.lt.rtmp) then
            rscale=sqrt(rtmp/rrmxsq)
            write(fidlog,*)"WARNING: Mobile particle ",idx,"'s centre-point is outside maximum radius"
            write(fidlog,*)" within the channel, rescaling rsx,rsy by ",(1/rscale)
            species(ispec)%sxyzr_(2:3, idx) = species(ispec)%sxyzr_(2:3, idx) / rscale
          endif
        enddo
      enddo

      !! RESCALE and/or TRANSLATE RSZ
      if (rzmin.lt.-zlim.or.rzmax.gt.zlim) then
        zoffset = 0.D0
        rscale = 1.D0
        ! rescale only if interval is greater than
        if (rzmax-rzmin.gt.2*zlim) then
          rscale = (2 * zlim)/(rzmax - rzmin)
          ! offset along z axis
          zoffset = -(rzmin + zlim)
        else
          if (rzmax.gt.zlim) then
            zoffset = zlim - rzmax
          else if (rzmin.lt.-zlim) then
            zoffset = -(rzmin + zlim)
          end if
        endif
        zoffset = zoffset / rscale
        write(fidlog,*)"WARNING: Mobile particle centre-points beyond cylinder end-points, adjusting all"
        write(fidlog,*)" mobile _and_ flexible ion rsz by:"
        write(fidlog,*)" scale : ",rscale
        write(fidlog,*)" translate : ",zoffset
        do ispec=1,imax_flex_
          do idx=1,species(ispec)%input_count_
            species(ispec)%sxyzr_(4, idx) = (zoffset + species(ispec)%sxyzr_(4, idx)) * rscale
          enddo
        enddo
      endif
    end subroutine rescale

    ! -----------------------------------------------------
    ! Echo the input file to file-id fid
    !
    ! Write out the interpreted specie and salt data from the input
    ! file in the same format as an input file.  This includes program
    ! default values for optional input data and normalised rate
    ! values.
    subroutine echoin(fid)
      implicit none
      integer, intent(in) :: fid

      ! locals
      integer :: ispec
      double precision, dimension(0,0) :: dummy
      integer, parameter :: dumsz=0

      do ispec=1,nspec_
        ! for species with particles defined output x y z (r) coordinates
        if (species(ispec)%input_count_.ne.0) then
          if (associated(species(ispec)%xyz_)) then
            ! need to get this information from species
            call ecspec(fid,ispec,species(ispec)%xyz_,species(ispec)%input_count_)
          else
            call ecspec(fid,ispec,dummy,dumsz)
          endif
        else
          call ecspec(fid,ispec,dummy,dumsz)
        endif
      enddo

    end subroutine echoin
  end subroutine rfspec

  double precision function rsx1(ii)
    integer, intent(in) :: ii
    integer :: ispec, idx
    idx = ii
    do ispec=1,imax_mob_+imax_flex_
      if (idx.gt.0.and.idx.le.species(ispec)%input_count_) then
        rsx1 = species(ispec)%sxyzr_(2,idx)
        return
      else
        idx = idx - species(ispec)%input_count_
      endif
    enddo
    stop "Particle index given to rsx is not a localized particle"
  end function rsx1

  double precision function rsx2(ispec, idx)
    integer, intent(in) :: ispec, idx
    if (idx.gt.0.and.idx.le.species(ispec)%input_count_) then
      rsx2 = species(ispec)%sxyzr_(2,idx)
    else
      stop "Particle specie/index given to rsx is not a localized particle"
    endif
  end function rsx2

  double precision function rsy1(ii)
    integer, intent(in) :: ii
    integer :: ispec, idx
    idx = ii
    do ispec=1,imax_mob_+imax_flex_
      if (idx.gt.0.and.idx.le.species(ispec)%input_count_) then
        rsy1 = species(ispec)%sxyzr_(3,idx)
        return
      else
        idx = idx - species(ispec)%input_count_
      endif
    enddo
    stop "Particle index given to rsy is not a localized particle"
  end function rsy1

  double precision function rsy2(ispec, idx)
    integer, intent(in) :: ispec, idx
    if (idx.gt.0.and.idx.le.species(ispec)%input_count_) then
      rsy2 = species(ispec)%sxyzr_(3,idx)
    else
      stop "Particle specie/index given to rsx is not a localized particle"
    endif
  end function rsy2

  double precision function rsz1(ii)
    integer, intent(in) :: ii
    integer :: ispec, idx
    idx = ii
    do ispec=1,imax_mob_+imax_flex_
      if (idx.gt.0.and.idx.le.species(ispec)%input_count_) then
        rsz1 = species(ispec)%sxyzr_(4,idx)
        return
      else
        idx = idx - species(ispec)%input_count_
      endif
    enddo
    stop "Particle index given to rsz is not a localized particle"
  end function rsz1

  double precision function rsz2(ispec, idx)
    integer, intent(in) :: ispec, idx
    if (idx.gt.0.and.idx.le.species(ispec)%input_count_) then
      rsz2 = species(ispec)%sxyzr_(4,idx)
    else
      stop "Particle specie/index given to rsx is not a localized particle"
    endif
  end function rsz2

  double precision function rsr1(ii)
    integer, intent(in) :: ii
    integer :: ispec, idx
    idx = ii
    do ispec=1,imax_mob_+imax_flex_
      if (idx.gt.0.and.idx.le.species(ispec)%input_count_) then
        rsr1 = species(ispec)%sxyzr_(1,idx)
        return
      else
        idx = idx - species(ispec)%input_count_
      endif
    enddo
    stop "Particle index given to rsr is not a localized particle"
  end function rsr1

  double precision function rsr2(ispec, idx)
    integer, intent(in) :: ispec, idx
    if (idx.gt.0.and.idx.le.species(ispec)%input_count_) then
      rsr2 = species(ispec)%sxyzr_(1,idx)
    else
      stop "Particle specie/index given to rsr is not a localized particle"
    endif
  end function rsr2

  ! ------------------------------------------------------------
  ! Select a salt using the given random number
  !
  ! REQUIRE(0 <= rndnum <= 1)
  integer function select_salt(rndnum)
    implicit none
    double precision, intent(in) :: rndnum
    double precision :: rnd_local
    integer :: idx
    if (dbc) then
      !! ASSERT(0 <= rndnum <= 1)
      if (0.gt.rndnum.or.1.lt.rndnum) then
        stop "Invalid random number not in range (0,1) in select_salt"
      end if
    end if
    rnd_local = rndnum
    select_salt = nsalt_
    if (nsalt_.eq.1) return
    m2w0s3: do idx=1,nsalt_
      t7z4l0: if (rnd_local.lt.salts(idx)%rate_change_) then
        select_salt = idx
        return
      endif t7z4l0
      rnd_local = rnd_local - salts(idx)%rate_change_
    enddo m2w0s3
  end function select_salt

  ! ------------------------------------------------------------
  ! Select a specie for a GC step using the given random number
  !
  ! REQUIRE(0 <= rndnum <= 1)
  integer function select_gc_specie(rndnum)
    implicit none
    double precision, intent(in) :: rndnum
    integer :: idx
    double precision :: rnd_local
    rnd_local = rndnum
    if (dbc) then
      !! ASSERT(0 <= rndnum <= 1)
      if (0.gt.rndnum.or.1.lt.rndnum) then
        stop "Invalid random number not in range (0,1) in select_gc_specie"
      end if
    end if
    select_gc_specie=nspec_
    m2w0s3: do idx=idxcl_,nspec_
      t7z4l0: if (rnd_local.lt.species(idx)%rate_change_) then
        select_gc_specie = idx
        return
      endif t7z4l0
      rnd_local = rnd_local - species(idx)%rate_change_
    enddo m2w0s3
  end function select_gc_specie

  ! ------------------------------------------------------------
  ! Select a specie for a move step using the given random number
  !
  ! if isbulk is true only free species are selected for a move
  ! REQUIRE(0 <= rndnum <= 1)
  integer function select_mv_specie(rndnum, isbulk)
    implicit none
    double precision, intent(in) :: rndnum
    logical, intent(in) :: isbulk
    integer :: idx, istart
    double precision :: rnd_local
    rnd_local = rndnum
    if (dbc) then
      !! ASSERT(0 <= rndnum <= 1)
      if (0.gt.rndnum.or.1.lt.rndnum) then
        stop "Invalid random number not in range (0,1) in select_mv_specie"
      end if
    end if
    if (isbulk) then
      istart = idxcl_
      rnd_local = rnd_local * sum(species(idxcl_:nspec_)%rate_specie_)
    else
      istart = 1
    endif
    select_mv_specie=nspec_
    m2w0s3: do idx = istart, nspec_
      rnd_local = rnd_local - species(idx)%rate_specie_
      t7z4l0: if (rnd_local.le.0.D0) then
        select_mv_specie = idx
        return
      endif t7z4l0
    enddo m2w0s3
  end function select_mv_specie

  ! ------------------------------------------------------------
  ! Select a subspecie for a swap move using the given random number
  !
  ! REQUIRE(0 <= rndnum <= 1)
  integer function select_subs(rndnum)
    implicit none
    double precision, intent(in) :: rndnum
    integer :: idx
    double precision :: rnd_local
    rnd_local = rndnum
    if (dbc) then
      !! ASSERT(0 <= rndnum <= 1)
      if (0.gt.rndnum.or.1.lt.rndnum) then
        stop "Invalid random number not in range (0,1) in select_subs"
      end if
    end if
    select_subs=nsubs_
    if (nsubs_.eq.1) return
    m2w0s3: do idx=1,nsubs_
      t7z4l0: if (rnd_local.lt.subspecies(idx)%rate_swap_) then
        select_subs = idx
        return
      endif t7z4l0
      rnd_local = rnd_local - subspecies(idx)%rate_swap_
    enddo m2w0s3
  end function select_subs

  ! ------------------------------------------------------------
  ! Select a region using the given random number
  !
  ! REQUIRE(0 <= rndnum <= 1)
  integer function select_region(ispec, rndnum)
    implicit none
    integer, intent(in) :: ispec
    double precision, intent(in) :: rndnum
    integer :: idx
    double precision :: rnd_local
    rnd_local = rndnum
    if (dbc) then
      !! ASSERT(0 <= rndnum <= 1)
      if (0.gt.rndnum.or.1.lt.rndnum) then
        stop "Invalid random number not in range (0,1) in select_region"
      end if
      if (1.gt.ispec.or.nspec_.lt.ispec) then
        stop "Invalid specie in select_region"
      end if
    end if
    select_region=ibulk
    m2w0s3: do idx=izlim,ibulk
      t7z4l0: if (rnd_local.lt.species(ispec)%rate_region_(idx)) then
        select_region = idx
        return
      endif t7z4l0
      rnd_local = rnd_local - species(ispec)%rate_region_(idx)
    enddo m2w0s3
  end function select_region

  ! ------------------------------------------------------------
  ! Select a move subtype using the given random number
  !
  ! REQUIRE(0 <= rndnum <= 1)
  integer function select_move(ispec, rndnum)
    implicit none
    integer, intent(in) :: ispec
    double precision, intent(in) :: rndnum
    double precision :: rnd_local
    rnd_local = rndnum
    if (dbc) then
      !! ASSERT(0 <= rndnum <= 1)
      if (0.gt.rndnum.or.1.lt.rndnum) then
        stop "Invalid random number not in range (0,1) in select_move"
      end if
      if (1.gt.ispec.or.nspec_.lt.ispec) then
        stop "Invalid specie in select_move"
      end if
    end if
    select_move = jmpin_
    if (localized(ispec)) then
      select_move = sphr_l
    else if (.not.isfree(ispec).or.rnd_local.gt.species(ispec)%rate_exchange_) then
      if (ranff().lt.species(ispec)%rate_move_) then
        select_move = sphr_l
      else
        select_move = jump_l
      endif
    endif
  end function select_move

  ! --------------------------------------------------
  ! Convert specie name into an index
  !
  ! ENSURE(spcidx = 0 on error | 1 <= spcidx <= nspec)
  !
  ! In DBC mode this calls stop if no match
  ! Otherwise it returns 0 for no match
  integer function spcidx(a_name)
    implicit none
    character(2), intent(in) :: a_name
    integer :: ispec, igc
    spcidx=0
    if (nspec_.le.0) return
    do ispec=1,nspec_
      if (a_name.eq.species(ispec)%code_name_) then
        spcidx=ispec
        return
      endif
    enddo
    ! If not found, look for label as a subspecie
    ! and return ground state ion index
    do igc=1,nsubs_
      if (a_name.eq.subspecies(igc)%code_name_) then
        ! perform search for ground state ion specie.
        do ispec=1,nspec_
          if (subspecies(igc)%sub_names_(1).eq.species(ispec)%code_name_) then
            spcidx=ispec
            return
          endif
        enddo
      endif
    enddo
    if (dbc) then
      write(unit=fidlog,fmt=*)"Specie name ",a_name," is unknown"
      stop 1
    endif
  end function spcidx

  ! --------------------------------------------------
  ! Convert specie name into an index excluding subspecies
  !
  ! ENSURE(spcidx = 0 on error | 1 <= spcidx <= nspec)
  !
  ! In DBC mode this calls stop if no match
  ! Otherwise it returns 0 for no match
  integer function private_spcidx(a_name)
    implicit none
    character(2), intent(in) :: a_name
    integer :: ispec
    private_spcidx=0
    if (nspec_.le.0) return
    do ispec=1,nspec_
      if (a_name.eq.species(ispec)%code_name_) then
        private_spcidx=ispec
        return
      endif
    enddo
    if (dbc) then
      write(unit=fidlog,fmt=*)"Specie name ",a_name," is unknown"
      stop 1
    endif
  end function private_spcidx

   ! --------------------------------------------------
  ! iterate the structural ions
  !
  ! Get information about the position of particles that where
  ! given in the input file.  To enhance data-hiding this method
  ! gives the 'conf' module a well-defined interface for retrieving
  ! the position data.  Once the 'conf' module has retrieved all
  ! the position information it can call 'deltmp' to let the 'spec'
  ! module know that it can delete the input xyz date for each
  ! specie.
  !
  ! REQUIRE(1 <= ispec <= nspec)
  ! REQUIRE(1 <= idc <= input_count(ispec))
  subroutine struks(ispec,idx,ax,ay,az)
    implicit none
    ! current specie/particle
    integer, intent(in) :: ispec,idx
    ! output as particle data
    double precision, intent(out) :: ax,ay,az
    if (dbc) then
      if (ispec.lt.1.or.ispec.gt.nspec_) stop "Invalid specie index"
      if (idx.lt.1.or.idx.gt.species(ispec)%input_count_) stop "Invalid initial particle index"
    endif
    ax = species(ispec)%xyz_(1, idx)
    ay = species(ispec)%xyz_(2, idx)
    az = species(ispec)%xyz_(3, idx)
  end subroutine struks

  ! ----------------------------------------------------------------
  ! Massage concentrations for subspecies.
  ! NOTE: ions can only appear in one subspecie so the
  ! following code assumes each pair can be treated 
  ! independently.
  !
  ! back calculates salt concentration and possibly adds any missing
  ! salt definitions.
  subroutine subs_init
    implicit none
#include "require.h"
    integer :: idx1
    integer :: igc, igc1, igc2
    integer :: ground_spec, excit_spec
    double precision :: rtargi, exponent
    logical :: is_igc1, is_igc2

    if (dbc) then
      if (nsubs_.le.0) stop "subs_init called with nsubs_ < 1"
    endif
    do idx1=1,nsubs_
      ! calculate total.
      ground_spec = subspecies(idx1)%sub_indices_(1)
      rtargi = species(ground_spec)%target_concentration_
      excit_spec = subspecies(idx1)%sub_indices_(2)
      rtargi = rtargi + species(excit_spec)%target_concentration_

      ! If specie 1 is ground state then dH - dS.T_ > 0
      !
      !  P(excited) ~= 1 / (1 + exp(dH - dS.T))
      !
      ! assign dG based on dH and dS
      exponent = subspecies(idx1)%enthalpy_ * (1000.D0) - subspecies(idx1)%entropy_ * tmptur()
      if (isnan(exponent).or.exponent.lt.0.D0) then
        write(*,*)"ERROR: Ground state subspecie ",species(ground_spec)%code_name_," (of subspecie ", &
             & subspecies(idx1)%code_name_,")"
        write(*,*)" is less favoured than excited state ",species(excit_spec)%code_name_
        stop "Ground state less favoured than excited state."
      endif
      ! divide by kT
      exponent = exponent * beta() / avog

      ! calculate the relative concentrations
      species(excit_spec)%target_concentration_ =  rtargi / (1.D0 + exp(exponent))
      species(ground_spec)%target_concentration_ =  rtargi - species(excit_spec)%target_concentration_

      if (isnan(species(ground_spec)%target_concentration_)) then
        write(*,*)"Subspecie ",species(ground_spec)%code_name_," of ", &
             & subspecies(idx1)%code_name_," has undefined target concentration; aborting run."
        stop "Ion has zero concentration"
      endif
      if (dfeq(species(ground_spec)%target_concentration_,0.D0)) then
        write(*,*)"Subspecie ",species(ground_spec)%code_name_," of ", &
             & subspecies(idx1)%code_name_," has zero target concentration; aborting run."
        stop "Ion has zero concentration"
      endif
      if (isnan(species(excit_spec)%target_concentration_)) then
        write(*,*)"Subspecie ",species(excit_spec)%code_name_," of ", &
             & subspecies(idx1)%code_name_," has undefined target concentration; aborting run."
        stop "Ion has zero concentration"
      endif
      if (dfeq(species(excit_spec)%target_concentration_,0.D0)) then
        write(*,*)"Subspecie ",species(excit_spec)%code_name_," of ", &
             & subspecies(idx1)%code_name_," has zero target concentration; aborting run."
        write(*,*)"enthalpy = ",subspecies(idx1)%enthalpy_ * (1000.D0)
        write(*,*)"entropy = ",subspecies(idx1)%entropy_
        write(*,*)"temperature = ",tmptur()
        write(*,*)"exponent = ",exponent
        write(*,*)"exp(exponent) = ",exp(exponent)
        write(*,*)"[",species(ground_spec)%code_name_,"] = ",species(ground_spec)%target_concentration_
        write(*,*)"[",species(excit_spec)%code_name_,"] = ",species(excit_spec)%target_concentration_
        stop "Ion has zero concentration"
      endif
    enddo
    igc1=-1
    igc2=-1
    ! back calculate salt concentration and possibly add any missing salt definitions.
    do idx1=1, nsubs_
      ground_spec = subspecies(idx1)%sub_indices_(1)
      excit_spec = subspecies(idx1)%sub_indices_(2)
      is_igc1=.false.
      is_igc2=.false.
      do igc=1,nsalt_
        if (salts(igc)%cation_index_.eq.ground_spec) then
          is_igc1=.true.
          igc1=igc
        endif
        if (salts(igc)%cation_index_.eq.excit_spec) then
          is_igc2=.true.
          igc2=igc
        endif
      enddo
      if (.not.(is_igc1.and.is_igc2)) then
        ! both false is an error
        if (.not.is_igc1.and..not.is_igc2) then
          stop "Can not define a sub specie pair without any corresponding salts"
        endif
        ! have a missing salt
        if (.not.is_igc1) then
          nsalt_ = nsalt_ + 1
          igc1 = nsalt_
          if (igc1.gt.nsltmx) stop "Too many salt definitions"
          call assign_salt(salts(igc1), salts(igc2))
          salts(igc1)%code_name_(1:2) = species(ground_spec)%code_name_
          salts(igc1)%cation_index_ = ground_spec
          salts(igc1)%target_concentration_ = 0.D0
          is_igc1 = .true.
        else if (.not.is_igc2) then
          nsalt_ = nsalt_ + 1
          igc2 = nsalt_
          if (igc2.gt.nsltmx) stop "Too many salt definitions"
          call assign_salt(salts(igc2), salts(igc1))
          salts(igc2)%code_name_(1:2) = species(excit_spec)%code_name_
          salts(igc2)%cation_index_ = excit_spec
          salts(igc2)%target_concentration_ = 0.D0
          is_igc2 = .true.
        endif
      endif
      salts(igc1)%target_concentration_ = species(ground_spec)%target_concentration_
      salts(igc2)%target_concentration_ = species(excit_spec)%target_concentration_

      !! ! reset excess chemical potentials for excited state
      !! if (dfeq(species(ground_spec)%chem_excess_,species(excit_spec)%chem_excess_)) then
      !!   ! assume chem_excess is unknown/guess for excited state if
      !!   ! equal to ground state value. Calculate based on the energy
      !!   ! of hydration.
      !!   species(excit_spec)%chem_excess_ = chempi(ground_spec) &
      !!        & - log(rtargi) - exponent
      !! endif
    end do
  end subroutine subs_init

  ! ----------------------------------------------------------------
  ! Get the specie index of a sub-specie component
  !
  ! REQUIRE(isubs < nsubs_ & (idx = 1 | idx = 2))
  integer function subspecie_index(isubs,idx)
    implicit none
    integer, intent(in) :: isubs, idx
    subspecie_index = subspecies(isubs)%sub_indices_(idx)
  end function subspecie_index

  ! ----------------------------------------------------------------
  ! Exchange the contents of two specie objects
  subroutine swapspec(spc1, spc2)
    implicit none
    type (specie), intent(inout) :: spc1, spc2
    character(2) :: fspctmp
    double precision, dimension(:,:), pointer :: xyzrtmp

    call swap(spc1%radius_,  spc2%radius_)
    call swap(spc1%valency_,   spc2%valency_)
    call swap(spc1%red_charge_,   spc2%red_charge_)
    call swap(spc1%chem_excess_, spc2%chem_excess_)
    ! Swap xyz
    if (associated(spc1%xyz_)) then
      xyzrtmp => spc1%xyz_
    else
      nullify(xyzrtmp)
    endif
    if (associated(spc2%xyz_)) then
      spc1%xyz_ => spc2%xyz_
    else
      nullify(spc1%xyz_)
    endif
    if (associated(xyzrtmp)) then
      spc2%xyz_ => xyzrtmp
    else
      nullify(spc2%xyz_)
    endif
    ! swap sxyzr
    if (associated(spc1%sxyzr_)) then
      xyzrtmp => spc1%sxyzr_
    else
      nullify(xyzrtmp)
    endif
    if (associated(spc2%sxyzr_)) then
      spc1%sxyzr_ => spc2%sxyzr_
    else
      nullify(spc1%sxyzr_)
    endif
    if (associated(xyzrtmp)) then
      spc2%sxyzr_ => xyzrtmp
    else
      nullify(spc2%sxyzr_)
    endif

    call swap(spc1%type_, spc2%type_)
    call swap(spc1%input_count_, spc2%input_count_)
    ! swap code name
    fspctmp = spc1%code_name_
    spc1%code_name_ = spc2%code_name_
    spc2%code_name_ = fspctmp
    call swap(spc1%rate_specie_, spc2%rate_specie_)
    call swap(spc1%rate_exchange_, spc2%rate_exchange_)
    call swap(spc1%rate_move_, spc2%rate_move_)
    call swap(spc1%rate_change_, spc2%rate_change_)
    call swap(spc1%rate_region_, spc2%rate_region_)
  end subroutine swapspec

  pure elemental double precision function xri(ispec)
    implicit none
    integer, intent(in) :: ispec
    xri=species(ispec)%radius_
  end function xri

  pure elemental double precision function xz(ispec)
    implicit none
    integer, intent(in) :: ispec
    xz=species(ispec)%valency_
  end function xz

  pure elemental double precision function xq(ispec)
    implicit none
    integer, intent(in) :: ispec
    xq=species(ispec)%red_charge_
  end function xq

end module spec

