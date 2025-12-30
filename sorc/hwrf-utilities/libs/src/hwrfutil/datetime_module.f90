module datetime_module
  use, intrinsic :: iso_c_binding, only: C_INT64_T, C_DOUBLE
  implicit none
  private

  integer, parameter, public :: grads_date_len=15

  public :: epochtime, ymdhtime, ymdhmtime, grads_date, set_timezone_utc
  public :: epochnow, epoch_split_ymdhm

  ! ------------------------------------------------------------

contains

  ! ------------------------------------------------------------

  subroutine grads_date(out,epochtime)
    ! Generates a string date in GrADS format, given the specified
    ! epoch time.  First converts the epoch time to a YYYYMMDDHHMM
    ! integer(kind=C_INT64_T) via ymdhmtime.  Any rounding is handled by that
    ! routine.  Return value is a character(len=grads_date_len).
    !
    ! Grads date format: April 20, 2014 at 04:30 UTC would return
    ! the string "04:30Z20apr2014".
    real(kind=C_DOUBLE), intent(in) :: epochtime
    character(len=grads_date_len) :: out
    ! 04:30Z20apr2014
    ! 123456789012345
    integer(kind=C_INT64_T) :: ymdhm
    integer :: imonth,i
    character(len=3), parameter :: cmonth(12) = &
         (/ 'jan', 'feb', 'mar', 'apr', 'may', 'jun', &
            'jul', 'aug', 'sep', 'oct', 'nov', 'dec' /)
    out='X'
    ymdhm=ymdhmtime(epochtime)
    imonth=mod(ymdhm/100/100/100,100)
    if(imonth>12 .or. imonth<1) then
       write(0,*) 'ERROR: ended up with invalid month ',imonth+1
       stop 30
    endif
30  format(I2,':',I2,'Z',I2,A3,I4)
35  format(I0," => ",I0,',',I0,',',I0,',',A,',',I0)
    ! write(6,35) ymdhm,mod(ymdhm/100,100), &
    !      mod(ymdhm,100), mod(ymdhm/10000,100), &
    !      cmonth(imonth), ymdhm/100000000
    write(out,30) mod(ymdhm/100,100), &
         mod(ymdhm,100), mod(ymdhm/10000,100), &
         cmonth(imonth), ymdhm/100000000
    do i=1,grads_date_len
       if(out(i:i)==' ') then
          out(i:i)='0'
       endif
    enddo
  end subroutine grads_date

  ! ------------------------------------------------------------

  subroutine set_timezone_utc(ierr)
    ! Sets the timezone to UTC.  If this is done, it should be done as
    ! early in the program as possible.
    integer(kind=4) :: ierr
    call c_etimeinit(ierr)
  end subroutine set_timezone_utc

  ! ------------------------------------------------------------

  function epochnow()
    ! Returns the current unix epoch time (seconds since 1970 began)
    ! with the precision of one microsecond, on platforms that
    ! correctly support the C gettimeofday routine.
    real(kind=C_DOUBLE) :: epochnow
    call c_now(epochnow)
  end function epochnow

  ! ------------------------------------------------------------

  function epochtime(ymdh)
    ! epochtime -- returns the unix epoch time for a YYYYMMDDHH time.
    !          This is a simple wrapper around the C c_etime subroutine.
    !   ymdh -- an integer of the form YYYYMMDDHH specifying a date/time
    !   returns -- number of seconds since midnight, January 1, 1970 UTC
    implicit none
    logical, save :: inited=.false.
    integer, intent(in) :: ymdh
    real(kind=C_DOUBLE) :: epochtime
    integer(kind=C_INT64_T) :: temp
    integer :: err

    if(.not.inited) then
       call c_etimeinit(err)
       inited=.true.
    endif

    call c_etime(ymdh,temp)
    epochtime=real(temp,kind=C_DOUBLE)
  end function epochtime

  ! ------------------------------------------------------------

  function ymdhtime(epoch)
    ! ymdhtime -- returns an integer YYYYMMDDHH time for a unix epoch time
    !       This is a simple wrapper around the C c_ytime subroutine.
    !     epoch -- integer number of seconds since midnight, January 1, 1970 UTC
    !     returns -- the equivalent integer YYYYMMDDHH time
    implicit none
    real(kind=C_DOUBLE), intent(in) :: epoch
    integer :: ymdhtime
    integer(kind=C_INT64_T) :: temp
    temp=nint(epoch,kind=C_INT64_T)
    call c_ytime(temp,ymdhtime)
  end function ymdhtime

  ! ------------------------------------------------------------

  function ymdhmtime(epoch)
    ! ymdhtime -- returns an eight byte integer YYYYMMDDHHMM time for
    !       a unix epoch time This is a simple wrapper around the C
    !       c_ymdhm subroutine.
    !     epoch -- integer number of seconds since midnight, January
    !       1, 1970 UTC
    !     returns -- the equivalent integer YYYYMMDDHHMM time as an
    !       eight byte integer.
    implicit none
    real(kind=C_DOUBLE), intent(in) :: epoch
    integer(kind=C_INT64_T) :: ymdhmtime
    integer(kind=C_INT64_T) :: temp
    temp=nint(epoch,kind=C_INT64_T)
    call c_ymdhm(temp,ymdhmtime)
  end function ymdhmtime

  subroutine epoch_split_ymdhm(epoch,year,mon,day,hr,min)
    real(kind=C_DOUBLE), intent(in) :: epoch
    integer(kind=C_INT64_T) :: year,mon,day,hr,min,sec, ymdhm

    ymdhm=ymdhmtime(epoch)
    year=ymdhm/int(100000000,kind=C_INT64_T)
    mon=mod(ymdhm/int(1000000,kind=C_INT64_T),100)
    day=mod(ymdhm/int(10000,kind=C_INT64_T),100) 
    hr=mod(ymdhm/int(100,kind=C_INT64_T),100)
    min=mod(ymdhm,100)
  end subroutine epoch_split_ymdhm
  
end module datetime_module
