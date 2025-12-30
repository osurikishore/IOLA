PROGRAM edit_gfsfrctime
  !
  ! Edit hycom time of template gfs.domain.forcing.param.[ab]. 
  ! The purpose of this set of codes is to build dummy forcing files
  ! which are required for only reading for hurricane-ocean coupled run.
  !
  USE mod_za,ONLY : xcspmd,zaiost,zaiopf,zaiowr,zaiocl,zaiord,idm,jdm
  !
  IMPLICIT NONE
  REAL :: fldmin,fldmax,otime,oldtime,time,time1,dtime,time0
  REAL, DIMENSION(:,:), ALLOCATABLE :: var1
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: imsk,ip

  INTEGER, PARAMETER :: lua=77,lub=lua+1001,luaa=15,luas=16
  INTEGER, PARAMETER :: level1=0

  INTEGER :: m,n,num_records,nheads,num
  CHARACTER(LEN=79), DIMENSION(:), ALLOCATABLE :: preambl
  CHARACTER(LEN=20) :: labela,fvar
  character(len=80) :: oldfrcprefix, oldAfrc, oldBfrc, newAfrc, newBfrc
  character(len=10) :: pdy
  !
  ! Read input parameters: edit_gfsfrc.in
  READ(5,'(a)') labela
  READ(5,'(a)') oldfrcprefix
  READ(5,*) nheads
  READ(5,*) num_records
  READ(5,*) pdy

  ! set input and output filenames 
  fvar=trim(labela(4:9))
  oldAfrc=trim(oldfrcprefix)//'.'//trim(fvar)//'.a'
  oldBfrc=trim(oldfrcprefix)//'.'//trim(fvar)//'.b'
  newAfrc='forcing.'//trim(fvar)//'.a'
  newBfrc='forcing.'//trim(fvar)//'.b'

  !
  call hytime(pdy,time0)
  time1=time0 - (6.-1.44)/24.
 
  WRITE(*,*) 'nheads, num_records, time0 =',nheads, num_records,time0
  write(*,*) ' '
  write(*,*) 'old a file=',oldAfrc
  write(*,*) 'old b file=',oldBfrc
  write(*,*) 'new a file=',newAfrc
  write(*,*) 'new b file=',newBfrc

  !
  ! Get grid parameters and allocate arrays.
  CALL xcspmd  !input idm,jdm by use association

  m=idm 
  n=jdm

  ALLOCATE ( VAR1(idm,jdm),preambl(nheads),imsk(idm,jdm),ip(idm,jdm))

  !
  ! Initialize I/O uints.
  CALL zaiost ! initialize units

  !
  ! Read HYCOM mask from regional.depth.[ab]
  CALL mask_hycom_2(imsk)
  WRITE(*,*) 'imsk min, max = ',MINVAL(imsk),MAXVAL(imsk)

  !
  ! Open output files and write header in 'forcing.airtmp.b file.
  CALL zaiopf(newAfrc,'new' , lua)
  OPEN (unit=lub,file=newBfrc,status='new' ,action='write')

  OPEN (unit=10,file=oldBfrc,status='old' ,action='read')
  READ (10,'(a79)') preambl
  CLOSE(10)

  ! NOTE: header will have the date of cration of the original file
  ! no mention that the field was corrected
  WRITE(lub,'(a79)') preambl
  !
  ! Loop over records
  records: DO num=1,num_records
     !
     ! Read values of forcing variable 
     CALL get_abfld1(luaa,ip,VAR1,m,n,labela, &
          level1,oldAfrc,nheads,oldtime) 
     PRINT *,'num=',num
          if (num .eq. 1) otime=oldtime; dtime=time1-otime
          time=oldtime+dtime
     IF (num==num_records) THEN 
        CALL zaiocl(luaa)
     ENDIF
     !
     !
     ! Write update forcinng.wndspd.a values
     CALL zaiowr(VAR1,imsk,.FALSE., fldmin,fldmax, lua,.FALSE.)
     WRITE(lub,'(A,'' date,span,range = '',F10.2,'' 0 '',1P2E16.7)') &
          &       trim(labela),time,fldmin,fldmax

     !
  ENDDO records
  !
  CLOSE(lub)
  CALL zaiocl(lua)
  !
  DEALLOCATE (VAR1)
  DEALLOCATE(imsk)
  DEALLOCATE(ip)
  DEALLOCATE(preambl)
  !
  !
CONTAINS
  !
  !========================================================================
  !
  SUBROUTINE mask_hycom_2(imsk)
    !
    !   Read HYCOM mask and calculate HYCOM land/sea mask (land=0,sea=1)
    !
    INTEGER, DIMENSION(:,:), INTENT(out) :: imsk
    REAL, DIMENSION(:,:), ALLOCATABLE :: dc
    REAL :: dmin,dmax
    INTEGER :: i,j,nxhycom,nyhycom,jdm
    nxhycom=SIZE(imsk,1)
    nyhycom=SIZE(imsk,2)
    ALLOCATE (dc(1:nxhycom,1:nyhycom))
    imsk=0
!!!!    call zaiost 
    CALL zaiopf('regional.depth.a','old', 61)
    CALL zaiord(dc,imsk,.FALSE.,dmin,dmax, 61)
    CALL zaiocl(61)
    DO j=1,nyhycom
       DO i=1,nxhycom
          IF(dc(i,j)>1.) THEN
             imsk(i,j)=1  
          ELSE 
             imsk(i,j)=0 
          ENDIF
       END DO
    END DO
  END SUBROUTINE mask_hycom_2
END PROGRAM edit_gfsfrctime

!==========================================================================
!
   subroutine hytime(date,hdate)
    implicit none
    real,intent(out) :: hdate
    character, intent(in):: date*(*)
    integer :: iyear,imonth,iday,ihour,jdn,day1
    integer,save  :: day0=-1
    external jdn
! ---   model day is calendar days since 12/31/1900
    if(day0<0) call julday(day0,1900,12,31)
    read(date(1:10),'(i4,i2,i2,i2)') iyear,imonth,iday,ihour
    call julday(day1,iyear,imonth,iday)
    hdate=day1-day0+ihour/24.d0
    return
  end subroutine hytime
!
!==========================================================================
!
   subroutine julday(jdn,iyear,month,iday)
      implicit none
      integer,intent(in) :: iyear,month,iday
      integer,intent(out) :: jdn
      jdn  =    iday - 32075 &
 &          + 1461 * (iyear + 4800 + (month - 14) / 12) / 4 &
 &          + 367 * (month - 2 - (month -14) / 12 * 12) / 12 &
 &          - 3 * ((iyear + 4900 + (month - 14) / 12) / 100) / 4
      return
  end subroutine julday
!
!===========================================================================================

