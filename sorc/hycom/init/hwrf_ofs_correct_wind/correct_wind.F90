PROGRAM correct_wind
  !
  ! Correct GFS surface wind, and wind stress from forcing.tau[en]wd1.[ab] 
  ! and forcing.wndspd1.[ab] files.
  !
  ! Sea surface wind and wind stress are computed from an independent source
  ! and these will replace GFS wind and wind stresses. 
  ! 
  ! The result is written to forcing.tau[en]wd.[ab] and forcing.wndspd.[ab].
  !
  !--------
  !   Note:
  !--------
  !      Designed to accomodate wind data from a variety of sources.
  !
  !---------
  ! history:
  !---------
  !  5/27/2008 Modified so to include time-varying weight. 
  !	This is designed based on the input value (nframes=5). 
  !                       by hyun-sook kim 
  !
  !  10/24/2007  Modified 'correct_forcing.f90' by Hyun-Sook Kim
  !
  !------------------------------------------------------------------------------------
  
  USE mod_za,ONLY : xcspmd,zaiost,zaiopf,zaiowr,zaiocl,zaiord,idm,jdm

  IMPLICIT NONE

  REAL :: fldmin,fldmax,time
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: VAR3d
  REAL, DIMENSION(:,:), ALLOCATABLE :: VAR1,VAR2,VAR0,VAR00,var01
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: imsk,ip

  real :: arg1, arg2, w1, w2,ref_val 
  real, parameter :: pi=3.141592
  real, parameter :: vel_mask=12.0, cd=0.0025, small_vel=0.01
  
  integer :: i,j, nt4t0
  real, dimension(:,:), allocatable :: tmpvar
  integer, dimension(:,:), allocatable :: matbin
  
  INTEGER, PARAMETER :: lua=77,lub=75,lp=6
  INTEGER, PARAMETER :: level1=0
  
  INTEGER :: m,n,nv,nt
  INTEGER :: nframes, nrcrds
  INTEGER :: ihead, iheada
  
  CHARACTER(LEN=79), DIMENSION(:), ALLOCATABLE :: preambl
  CHARACTER(LEN=20) :: labela3v(3), labela
  CHARACTER(LEN=16) :: BFNAME3(3),AFNAME3(3),OBFNAME3(3)
  CHARACTER(LEN=16) :: BFNAME,AFNAME,OBFNAME
  
  CHARACTER(LEN=80) :: DFNAME   
  
!
! Read input parameters from 'correct_wind.in'.
  READ(5,'(a)') dfname
  READ(5,'(3(a16,1x))') bfname3
  READ(5,'(3(a16,1x))') afname3
  READ(5,'(3(a16,1x))') obfname3
  READ(5,*) ihead
  READ(5,'(3(a10,1x))') labela3v
  READ(5,*) nframes
  READ(5,*) nrcrds
  
!
! following 2 variables: fixed
!dhi   nframes=5		! fixed. caveat: time-interval=3-hr
  nt4t0=4		! t=0 for next cycle
  
!
! Get grid parameters and allocate arrays.
    CALL XCSPMD   			!input idm,jdm by use association

m=idm
n=jdm

! Initialize I/O uints.
    CALL ZAIOST  				! initialize units


ALLOCATE ( preambl(ihead) )
ALLOCATE ( ip(m,n) )
ALLOCATE ( imsk(m,n) )
ALLOCATE ( var0(m,n) )
ALLOCATE ( var00(m,n) )
ALLOCATE ( var01(m,n) )
ALLOCATE ( var1(m,n) )
ALLOCATE ( var2(m,n) )
ALLOCATE ( VAR3d(m,n,4) )
ALLOCATE ( tmpvar(m,n) )
ALLOCATE ( matbin(m,n) )
  
!
! Read HYCOM mask from regional.mask.[ab]
CALL MASK_HYCOM_2(imsk) 
WRITE(*,*) 'imsk min, max = ',MINVAL(imsk),MAXVAL(imsk)
 
!----------------- work on each variable ---------------------
!
! Open the individual output files and write header 
!     in e.g., 'forcing.wndspd.b' file.
!      1. forcing.wndspd.b
!      2. forcing.tauewd.b
!      3. forcing.taunwd.b
!  
DO NV=1,3

!
! read in input data from 'DFNAME'
    OPEN (unit=44,file=trim(DFNAME),status='old',form='unformatted')

    afname=afname3(nv)
    bfname=bfname3(nv)
    obfname=obfname3(nv)
      

    CALL ZAIOPF(AFNAME,'new',lua)

    OPEN (unit=lua,file=trim(BFNAME),status='new' ,action='write')

!get the preamble from old file
    OPEN (unit=10,file=trim(OBFNAME),status='old' ,action='read')
    READ (10,'(a79)') preambl
    CLOSE(10)

! NOTE: header will have the date of creation of the original file
! not mention that the field was corrected
     WRITE(lua,'(a79)') preambl
     iheada=ihead
     labela=trim(labela3v(nv))
     write(*,*)labela

     if(nv.eq.1) then
       ref_val=vel_mask
       write(*,*)'Cutoff velocity for background ',ref_val
     else
       ref_val=vel_mask*vel_mask*cd*1.4
       write(*,*)'Cutoff stress for background ',ref_val
     endif
!
! Loop over records
records: DO nt=1,nframes
     !
     ! Loop over each variable - (wind_speed, taux, and tauy)
     READ(44) VAR3d

     VAR0=var3d(:,:,1)    !fm
     VAR01=var3d(:,:,2)   !speed
     VAR1=VAR3d(:,:,nv+1)
     CALL get_abfld1(lub,ip,VAR2,m,n,labela, &
          level1,obfname,iheada,time) 

!------------------------------------------------
!    cutoff high values  of  GFS fields if close to  FM hurricane
          
     DO i=1,idm
       DO j=1,jdm
         if ( var01(i,j).gt.small_vel.or.var0(1,j).gt. 0.99) then
           if (var2(i,j).gt.ref_val) then
             var2(i,j)=ref_val
           elseif(var2(i,j).lt.-ref_val) then
             var2(i,j)=-ref_val
           endif
         endif
        ENDDO
     ENDDO

!------------------------------------------------
!    if a wind stress component, then rescale weight
     if (nv.ge.2) then      
       DO i=1,idm
       DO j=1,jdm
         var0(i,j)=min(1.0,var0(i,j)**2)
       ENDDO
       ENDDO
     endif
!------------------------------------------------
! melding VAR1 (H-wind) and VAR2 (GFS)
     
     DO i=1,idm
       DO j=1,jdm
	      Var1(i,j)=Var2(i,j)*(1-var0(i,j))+var0(i,j)*var1(i,j)
       ENDDO
     ENDDO
	          
     !
     ! 5: Write corrected wind entity.
     !--------------------------------
     CALL zaiowr(VAR1,imsk,.TRUE., fldmin,fldmax, lua,.FALSE.)

     WRITE(lua,'(a10," date,span,range = ",F10.2,'' 0 '',1P2E16.7)') &
          & labela(1:10),time,fldmin,fldmax
     WRITE(lp,'(a10," date,span,range = ",F10.2,'' 0 '',1P2E16.7)') &
          & labela(1:10),time,fldmin,fldmax

     !
!     write(1000*nv+nt,'(e18.6)')var1
!     close(1000*nv+nt)

!     write(5000*nv+nt,'(e18.6)')var3d
!     close(5000*nv+nt)
     
     ENDDO records
!


     call zaiocl(lub)
     CLOSE(lub)
!
! close files:

     CALL zaiocl(lua)     
     CLOSE(lua)
!
     CLOSE(44) 
ENDDO 

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
    INTEGER :: i,j,nxhycom,nyhycom
    nxhycom=SIZE(imsk,1)
    nyhycom=SIZE(imsk,2)
    ALLOCATE (dc(1:nxhycom,1:nyhycom))
    imsk=0
!!!!    call zaiost 
    CALL zaiopf('regional.mask.a','old', 61)
    CALL zaiord(dc,imsk,.FALSE.,dmin,dmax, 61)
    CALL zaiocl(61)
    DO j=1,nyhycom
       DO i=1,nxhycom
! old_HYCOM:        IF(dc(i,j)>0.) THEN
          IF(dc(i,j)<1E+20) THEN
             imsk(i,j)=1  
          ELSE 
             imsk(i,j)=0 
          ENDIF
       END DO
    END DO
  END SUBROUTINE mask_hycom_2
END PROGRAM correct_wind


!========================================================================
#if defined(AIX)
  subroutine flush(unit)
    integer unit
    call flush_(unit)
  end subroutine flush
#endif /* AIX */
#if defined(IFC)
  subroutine flush(unit)
    integer unit
  end subroutine flush
#endif /* IFC */
