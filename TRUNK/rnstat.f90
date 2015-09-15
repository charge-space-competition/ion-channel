
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
! Running Statistics Module
! 
! Simple module for objects that hold the running average and
! variance of a variable
! -------------------------------------------------------------
! Types:
! type (runstat)
!    object for the running statistics.  This needs to be
!    initialised with rs_clear before use.
! -------------------------------------------------------------
! Methods:
! rs_init, rs_reset
!    intialise or reset the variable.
! rs_push
!    add a new estimate to the variable.
! rs_mean
!    get the current mean of the variable.
! rs_variance
!    get the current variance of the variable.
! rs_stddev
!    get the current standard deviation of the variable.
! rs_count
!    get the number of samples of the variable.
! -------------------------------------------------------------

module rnstat
  use const
  implicit none
  public
  type hist1array
    integer (kind=Cntr_K), allocatable, dimension(:) :: bincounters
    integer (kind=Cntr_K) :: samples
    integer :: size
    logical :: inpush
  end type hist1array

  type hist2array
    integer (kind=Cntr_K), allocatable, dimension(:,:) :: bincounters
    integer (kind=Cntr_K) :: samples
    integer, dimension(2) :: size
    logical :: inpush
  end type hist2array

  type hist3array
    integer (kind=Cntr_K), allocatable, dimension(:,:,:) :: bincounters
    integer (kind=Cntr_K) :: samples
    integer, dimension(3) :: size
    logical :: inpush
  end type hist3array

  type hist4array
    integer (kind=Cntr_K), allocatable, dimension(:,:,:,:) :: bincounters
    integer (kind=Cntr_K) :: samples
    integer, dimension(4) :: size
    logical :: inpush
  end type hist4array

  type runstat
    double precision :: mean, variance
    integer (kind=Cntr_K) :: count
  end type runstat

  ! initialise a histogram to a particular size
  interface hist_init
    module procedure hist1_init
    module procedure hist2_init
    module procedure hist3_init
    module procedure hist4_init
  end interface hist_init

  ! remove all the storage associated with a histogram
  interface hist_delete
    module procedure hist1_delete
    module procedure hist2_delete
    module procedure hist3_delete
    module procedure hist4_delete
  end interface hist_delete

  ! reset all samples in the the histogram to 0
  interface hist_reset
    module procedure hist1_reset
    module procedure hist2_reset
    module procedure hist3_reset
    module procedure hist4_reset
  end interface hist_reset

  ! add data to the given histogram position.  This method only increments the samples
  ! counter once for each time hist_end_sample is called.
  interface hist_push
    module procedure hist1_sample
    module procedure hist2_sample
    module procedure hist3_sample
    module procedure hist4_sample
  end interface hist_push

  ! Indicate that a sampling sequence is completed
  interface hist_end_push
    module procedure hist1_end_sample
    module procedure hist2_end_sample
    module procedure hist3_end_sample
    module procedure hist4_end_sample
  end interface hist_end_push

  ! Get the number of sampling sequences completed
  interface hist_count
    module procedure hist1_samples
    module procedure hist2_samples
    module procedure hist3_samples
    module procedure hist4_samples
  end interface hist_count

  ! Get the mean number of entries in a histogram bin
  interface hist_mean
    module procedure hist1_mean
    module procedure hist2_mean
    module procedure hist3_mean
    module procedure hist4_mean
  end interface hist_mean

  ! get the sum of all entries in a bin
  interface hist_sum
    module procedure hist1_sum
    module procedure hist2_sum
    module procedure hist3_sum
    module procedure hist4_sum
  end interface hist_sum
  
  interface rs_reset
    module procedure rs_init
  end interface rs_reset
  !public rs_init, rs_push, rs_count, rs_mean, rs_variance, rs_stddev

contains
  ! -------------------------------------------------------------
  ! Clear or initialise or reset average/variance pairs
  subroutine hist1_init(var, asize)
  implicit none
  type (hist1array), intent(inout) :: var
  integer, intent(in) :: asize
  var%samples = 0
  var%inpush=.false.
  if (allocated(var%bincounters)) then
    if (asize.ne.var%size) then
      var%size = 0
      deallocate(var%bincounters)
      if (0.ne.var%size) then
        allocate(var%bincounters(asize))
        var%size = asize
      endif
    endif
  else
    var%size = 0
    allocate(var%bincounters(asize))
    var%size = asize
  endif
  if (0.ne.var%size) var%bincounters = 0
  end subroutine hist1_init

  ! -------------------------------------------------------------
  ! Release any stored memory and reset size to 0
  subroutine hist1_delete(var)
  implicit none
  type (hist1array), intent(inout) :: var
  var%samples = 0
  var%inpush=.false.
  if (allocated(var%bincounters)) then
    var%size = 0
    deallocate(var%bincounters)
  endif
  end subroutine hist1_delete

  ! -------------------------------------------------------------
  ! Add to an average/variance pair at the given index
  subroutine hist1_sample(var, index, count)
  implicit none
  type (hist1array), intent(inout) :: var
  integer, intent(in) :: index
  integer, intent(in) :: count
  if (.not.allocated(var%bincounters)) then
    stop "Attempt to use runarray before initialising"
  endif
  if (.not.var%inpush) then
    var%samples=var%samples + 1
    var%inpush = .true.
  endif
  var%bincounters(index)=var%bincounters(index) + count
  end subroutine hist1_sample

  ! -------------------------------------------------------------
  ! Tell variable that you have finished pushing bin counts for
  ! this sample
  subroutine hist1_end_sample(var)
  implicit none
  type (hist1array), intent(inout) :: var
  var%inpush = .false.
  end subroutine hist1_end_sample

  ! -------------------------------------------------------------
  ! Get the current average
  function hist1_mean(var, index)
  implicit none
  type (hist1array), intent(inout) :: var
  integer, intent(in) :: index
  double precision :: hist1_mean
  if (var%samples.eq.0) then
    hist1_mean = 0.0D0
  else
    hist1_mean = dble(var%bincounters(index))/dble(var%samples)
  endif
  end function hist1_mean

  ! -------------------------------------------------------------
  ! Get the current average
  function hist1_sum(var, index)
  implicit none
  type (hist1array), intent(inout) :: var
  integer, intent(in) :: index
  integer (kind=Cntr_K) :: hist1_sum
  if (var%samples.eq.0) then
    hist1_sum = 0
  else
    hist1_sum = var%bincounters(index)
  endif
  end function hist1_sum

  ! -------------------------------------------------------------
  ! Get the current average
  subroutine hist1_reset(var)
  implicit none
  type (hist1array), intent(inout) :: var
  var%samples = 0
  var%inpush=.false.
  if (0.ne.var%size) var%bincounters = 0
  end subroutine hist1_reset

  ! -------------------------------------------------------------
  ! Get the current count
  function hist1_samples(var)
  implicit none
  type (hist1array), intent(inout) :: var
  integer (kind=Cntr_K) :: hist1_samples
  hist1_samples = var%samples
  end function hist1_samples

  ! -------------------------------------------------------------
  ! Clear or initialise or reset average/variance pairs
  subroutine hist2_init(var, asize1, asize2)
  implicit none
  type (hist2array), intent(inout) :: var
  integer, intent(in) :: asize1, asize2
  var%samples = 0
  var%inpush=.false.
  if (allocated(var%bincounters)) then
    if (asize1.ne.var%size(1).or.asize2.ne.var%size(2)) then
      var%size = 0
      deallocate(var%bincounters)
      allocate(var%bincounters(asize1,asize2))
      var%size(1) = asize1
      var%size(2) = asize2
    endif
  else
    var%size = 0
    allocate(var%bincounters(asize1,asize2))
    var%size(1) = asize1
    var%size(2) = asize2
  endif
  var%bincounters = 0
  end subroutine hist2_init

  ! -------------------------------------------------------------
  ! Release any stored memory and reset size to 0
  subroutine hist2_delete(var)
  implicit none
  type (hist2array), intent(inout) :: var
  var%samples = 0
  var%inpush=.false.
  if (allocated(var%bincounters)) then
    var%size = 0
    deallocate(var%bincounters)
  endif
  end subroutine hist2_delete

  ! -------------------------------------------------------------
  ! Reset average/variance pairs
  subroutine hist2_reset(var)
  implicit none
  type (hist2array), intent(inout) :: var
  var%samples = 0
  var%inpush=.false.
  var%bincounters = 0
  end subroutine hist2_reset
  ! -------------------------------------------------------------
  ! Add to an average/variance pair at the given index
  subroutine hist2_sample(var, i, j, count)
  implicit none
  type (hist2array), intent(inout) :: var
  integer, intent(in) :: i, j
  integer, intent(in) :: count
  if (.not.allocated(var%bincounters)) then
    stop "Attempt to use runarray before initialising"
  endif
  if (.not.var%inpush) then
    var%samples=var%samples + 1
    var%inpush = .true.
  endif
  var%bincounters(i,j)=var%bincounters(i,j) + count
  end subroutine hist2_sample

  ! -------------------------------------------------------------
  ! Tell variable that you have finished pushing bin counts for
  ! this sample
  subroutine hist2_end_sample(var)
  implicit none
  type (hist2array), intent(inout) :: var
  var%inpush = .false.
  end subroutine hist2_end_sample

  ! -------------------------------------------------------------
  ! Get the current average
  function hist2_mean(var, i, j)
  implicit none
  type (hist2array), intent(inout) :: var
  integer, intent(in) :: i, j
  double precision :: hist2_mean
  if (var%samples.eq.0) then
    hist2_mean = 0.0D0
  else
    hist2_mean = dble(var%bincounters(i, j))/dble(var%samples)
  endif
  end function hist2_mean

  ! -------------------------------------------------------------
  ! Get the current average
  function hist2_sum(var, i, j)
  implicit none
  type (hist2array), intent(inout) :: var
  integer, intent(in) :: i, j
  integer (kind=Cntr_K) :: hist2_sum
  if (var%samples.eq.0) then
    hist2_sum = 0
  else
    hist2_sum = var%bincounters(i, j)
  endif
  end function hist2_sum

  ! -------------------------------------------------------------
  ! Get the current count
  function hist2_samples(var)
  implicit none
  type (hist2array), intent(inout) :: var
  integer (kind=Cntr_K) :: hist2_samples
  hist2_samples = var%samples
  end function hist2_samples

  ! -------------------------------------------------------------
  ! Clear or initialise or reset bin counters
  subroutine hist3_init(var, asize1, asize2, asize3)
  implicit none
  type (hist3array), intent(inout) :: var
  integer, intent(in) :: asize1, asize2, asize3
  var%samples = 0
  var%inpush=.false.
  if (allocated(var%bincounters)) then
    if (asize1.ne.var%size(1).or.asize2.ne.var%size(2).or.asize3.ne.var%size(3)) then
      var%size = 0
      deallocate(var%bincounters)
      allocate(var%bincounters(asize1,asize2,asize3))
      var%size(1) = asize1
      var%size(2) = asize2
      var%size(3) = asize3
    endif
  else
    var%size = 0
    allocate(var%bincounters(asize1, asize2, asize3))
    var%size(1) = asize1
    var%size(2) = asize2
    var%size(3) = asize3
  endif
  var%bincounters = 0
  end subroutine hist3_init

  ! -------------------------------------------------------------
  ! Clear or initialise or reset bin counters
  subroutine hist3_reset(var)
  implicit none
  type (hist3array), intent(inout) :: var
  var%samples = 0
  var%inpush=.false.
  var%bincounters = 0
  end subroutine hist3_reset

  ! -------------------------------------------------------------
  ! Release any stored memory and reset size to 0
  subroutine hist3_delete(var)
  implicit none
  type (hist3array), intent(inout) :: var
  var%samples = 0
  var%inpush=.false.
  if (allocated(var%bincounters)) then
    var%size = 0
    deallocate(var%bincounters)
  endif
  end subroutine hist3_delete

  ! -------------------------------------------------------------
  ! Add to an average/variance pair at the given index
  subroutine hist3_sample(var, i, j, k, count)
  implicit none
  type (hist3array), intent(inout) :: var
  integer, intent(in) :: i,j,k
  integer, intent(in) :: count
  if (.not.allocated(var%bincounters)) then
    stop "Attempt to use runarray before initialising"
  endif
  if (.not.var%inpush) then
    var%samples=var%samples + 1
    var%inpush = .true.
  endif
  var%bincounters(i,j,k)=var%bincounters(i,j,k) + count
  end subroutine hist3_sample

  ! -------------------------------------------------------------
  ! Tell variable that you have finished pushing bin counts for
  ! this sample
  subroutine hist3_end_sample(var)
  implicit none
  type (hist3array), intent(inout) :: var
  var%inpush = .false.
  end subroutine hist3_end_sample

  ! -------------------------------------------------------------
  ! Get the current average
  function hist3_mean(var, i, j, k)
  implicit none
  type (hist3array), intent(inout) :: var
  integer, intent(in) :: i,j,k
  double precision :: hist3_mean
  if (var%samples.eq.0) then
    hist3_mean = 0.0D0
  else
    hist3_mean = dble(var%bincounters(i, j, k))/dble(var%samples)
  endif
  end function hist3_mean

  ! -------------------------------------------------------------
  ! Get the current average
  function hist3_sum(var, i, j, k)
  implicit none
  type (hist3array), intent(inout) :: var
  integer, intent(in) :: i,j,k
  integer (kind=Cntr_K) :: hist3_sum
  if (var%samples.eq.0) then
    hist3_sum = 0
  else
    hist3_sum = var%bincounters(i, j, k)
  endif
  end function hist3_sum

  ! -------------------------------------------------------------
  ! Get the current count
  function hist3_samples(var)
  implicit none
  type (hist3array), intent(inout) :: var
  integer (kind=Cntr_K) :: hist3_samples
  hist3_samples = var%samples
  end function hist3_samples

  ! -------------------------------------------------------------
  ! Clear or initialise or reset bin counters
  subroutine hist4_init(var, sz1, sz2, sz3, sz4)
  implicit none
  type (hist4array), intent(inout) :: var
  integer, intent(in) :: sz1, sz2, sz3, sz4
  var%samples = 0
  var%inpush=.false.
  if (allocated(var%bincounters)) then
    if (sz1.ne.var%size(1).or.sz2.ne.var%size(2).or.sz3.ne.var%size(3).or.sz4.ne.var%size(4)) then
      var%size = 0
      deallocate(var%bincounters)
      allocate(var%bincounters(sz1,sz2,sz3,sz4))
      var%size(1) = sz1
      var%size(2) = sz2
      var%size(3) = sz3
      var%size(4) = sz4
    endif
  else
    var%size = 0
    allocate(var%bincounters(sz1,sz2,sz3,sz4))
    var%size(1) = sz1
    var%size(2) = sz2
    var%size(3) = sz3
    var%size(4) = sz4
  endif
  var%bincounters = 0
  end subroutine hist4_init

  ! -------------------------------------------------------------
  ! Clear or initialise or reset bin counters
  subroutine hist4_reset(var)
  implicit none
  type (hist4array), intent(inout) :: var
  var%samples = 0
  var%inpush=.false.
  var%bincounters = 0
  end subroutine hist4_reset

  ! -------------------------------------------------------------
  ! Release any stored memory and reset size to 0
  subroutine hist4_delete(var)
  implicit none
  type (hist4array), intent(inout) :: var
  var%samples = 0
  var%inpush=.false.
  if (allocated(var%bincounters)) then
    var%size = 0
    deallocate(var%bincounters)
  endif
  end subroutine hist4_delete

  ! -------------------------------------------------------------
  ! Add to an average/variance pair at the given index
  subroutine hist4_sample(var, i, j, k, l, count)
  implicit none
  type (hist4array), intent(inout) :: var
  integer, intent(in) :: i,j,k,l
  integer, intent(in) :: count
  if (.not.allocated(var%bincounters)) then
    stop "Attempt to use runarray before initialising"
  endif
  if (.not.var%inpush) then
    var%samples=var%samples + 1
    var%inpush = .true.
  endif
  var%bincounters(i,j,k,l)=var%bincounters(i,j,k,l) + count
  end subroutine hist4_sample

  ! -------------------------------------------------------------
  ! Tell variable that you have finished pushing bin counts for
  ! this sample
  subroutine hist4_end_sample(var)
  implicit none
  type (hist4array), intent(inout) :: var
  var%inpush = .false.
  end subroutine hist4_end_sample

  ! -------------------------------------------------------------
  ! Get the current average
  function hist4_mean(var, i, j, k, l)
  implicit none
  type (hist4array), intent(inout) :: var
  integer, intent(in) :: i,j,k,l
  double precision :: hist4_mean
  if (var%samples.eq.0) then
    hist4_mean = 0.0D0
  else
    hist4_mean = dble(var%bincounters(i, j, k, l))/dble(var%samples)
  endif
  end function hist4_mean

  ! -------------------------------------------------------------
  ! Get the current average
  function hist4_sum(var, i, j, k, l)
  implicit none
  type (hist4array), intent(inout) :: var
  integer, intent(in) :: i,j,k,l
  integer (kind=Cntr_K) :: hist4_sum
  if (var%samples.eq.0) then
    hist4_sum = 0
  else
    hist4_sum = var%bincounters(i, j, k, l)
  endif
  end function hist4_sum

  ! -------------------------------------------------------------
  ! Get the current count
  function hist4_samples(var)
  implicit none
  type (hist4array), intent(inout) :: var
  integer (kind=Cntr_K) :: hist4_samples
  hist4_samples = var%samples
  end function hist4_samples

  ! -------------------------------------------------------------
  ! Clear or initialise or reset average/variance pair
  subroutine rs_init(var)
  implicit none
  type (runstat), intent(inout) :: var
  var%count = 0
  end subroutine rs_init

  ! -------------------------------------------------------------
  ! Add to an average/variance pair
  subroutine rs_push(var, val)
  implicit none
  type (runstat), intent(inout) :: var
  double precision, intent(in) :: val
  double precision :: av_old
  var%count=var%count + 1
  if (var%count.eq.1) then
    var%mean = val
    var%variance = 0.0D0
  else
    av_old = var%mean
    var%mean = av_old + (val - av_old)/dble(var%count)
    var%variance = var%variance + (val - av_old)*(val - var%mean)
  endif
  end subroutine rs_push

  ! -------------------------------------------------------------
  ! Get the current average
  function rs_mean(var)
  implicit none
  type (runstat), intent(inout) :: var
  double precision :: rs_mean
  if (var%count.eq.0) then
    rs_mean = 0.0D0
  else
    rs_mean = var%mean
  endif
  end function rs_mean

  ! -------------------------------------------------------------
  ! Get the current count
  function rs_count(var)
  implicit none
  type (runstat), intent(inout) :: var
  integer (kind=Cntr_K) :: rs_count
  rs_count = var%count
  end function rs_count

  ! -------------------------------------------------------------
  ! Variance of the variable
  function rs_variance(var)
  implicit none
  type (runstat), intent(inout) :: var
  double precision :: rs_variance
  if (var%count.le.1) then
    rs_variance = 0.0D0
  else
    rs_variance = var%variance/dble(var%count - 1)
  endif
  end function rs_variance

  ! -------------------------------------------------------------
  ! Variance of the variable
  function rs_stddev(var)
  implicit none
  type (runstat), intent(inout) :: var
  double precision :: rs_stddev
  if (var%count.le.1) then
    rs_stddev = 0.0D0
  else
    rs_stddev = sqrt(var%variance/dble(var%count - 1))
  endif
  end function rs_stddev


end module rnstat
