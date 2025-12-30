PROGRAM concate_forcing
  !
  ! Add forcing1 and forcing2, to generate forcing3.
  ! It includes sorting.
  !    ./temp1/forcing.OFIELD.[ab] for times
  ! extracted from forcing file in 
  !    ./temp2/forcing.IFIELD.b
  ! write expanded file to
  !    ./forcing.OFIELD.[ab] 
  ! 
  !
  USE mod_za,ONLY : xcspmd,zaiost,zaiopf,zaiowr,zaiocl,zaiord,idm,jdm
  !
  IMPLICIT NONE

  REAL :: fldmin,fldmax,time,time_old,time_star
  REAL, DIMENSION(:,:), ALLOCATABLE :: T_a
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: imsk,ip
  INTEGER, PARAMETER :: lua=88,lub=lua+1001,lua1=15,lua2=16
  INTEGER :: m,n,ihead,knt1,knt,level1,num,iheada,iheads
  CHARACTER(LEN=79), DIMENSION(:), ALLOCATABLE :: preambl
  CHARACTER(LEN=11) :: labela,labels
  CHARACTER(LEN=10) :: ifile1='forcing1.b',ifile2='forcing2.b'
  CHARACTER(LEN=16) :: afile
  character(len=14) :: oname='forcing.123456'
  character(len=6)  :: ofield
  integer :: mark
  real , dimension(:), allocatable ::g_time
  real :: wei
  real ,parameter :: time_now=-1.
  integer, parameter :: nmax=200
  !
  ! Read input parameters
  read(5,'(a)') ofield
  read(5,'(a)') ifile1
  read(5,'(a)') ifile2
  READ(5,*) ihead

  ! set up label in a
  labela = '   '//ofield//':'

!
  ALLOCATE( preambl(ihead) )
  ALLOCATE(g_time(nmax))
!
!
! write(*,*),'file1=','temp1/'//TRIM(ifile1)
  OPEN (unit=10,file='./temp1/'//TRIM(ifile1),status='old' ,action='read')
  READ (10,'(a79)') preambl
!
  knt=0
  do 
    knt=knt+1
    read(10,'(28x,f11.2)',end=111) g_time(knt)
  enddo
111 continue
  knt1=knt-1
  write(*,*)'First target time for file 1 ',g_time(1)
  write(*,*)'Last  target time for file 1 ',g_time(knt1)
  write(*,*)'Number of records for file 1 = ', knt1
  close(10)
!
!
  OPEN (unit=20,file='./temp2/'//TRIM(ifile2),status='old' ,action='read')
  READ (20,'(a79)') preambl
!
  knt=knt1
  do 
    knt=knt+1
    read(20,'(28x,f11.2)',end=222) g_time(knt)
  enddo
222 continue
  write(*,*)'First target time for file 2 ',g_time(knt1+1)
  write(*,*)'Last  target time for file 2 ',g_time(knt-1)
  write(*,*)'Number of records for file 2 = ', knt-knt1-1
  close(20)

!
! --- output
!
  oname(9:9+6-1)=ofield
  ! 
  !
  ! Get grid parameters and allocate arrays.
  CALL xcspmd  !input idm,jdm by use association
  m=idm 
  n=jdm
  ALLOCATE ( imsk(idm,jdm),ip(idm,jdm),T_a(idm,jdm) )

  !
  ! Initialize I/O uints.
   CALL zaiost ! initialize units
  !
  ! Open output files and write header in 'forcing.airtmp.b file.
   afile=oname//'.a'
   CALL zaiopf(afile,'new' , lua)
   OPEN (unit=lub,file=oname//'.b',status='new' ,action='write')

  ! NOTE: header will have the date of creation of the original file
  ! no mention that the field was corrected
   WRITE(lub,'(a79)') preambl

  level1=0

  iheada=ihead  
  iheads=ihead
  
  !
  ! Loop over records
  num=0
  mark=0
  records1: DO while(mark<knt1)
     mark=mark+1
     !
     ! Read ofile
     CALL get_abfld1(lua1,ip,T_a,m,n,labela, &
          level1,'./temp1/'//ifile1,iheada,time) 

     if (time==g_time(mark)) then
       if (g_time(mark)>=maxval(g_time(:mark-1))) then
          CALL zaiowr(T_a,imsk,.FALSE., fldmin,fldmax, lua,.FALSE.)
          WRITE(lub,'(A,'' date,span,range = '',F10.2,'' 0 '',1P2E16.7)') &
              &       labela(1:10),time,fldmin,fldmax
          num=num+1
        endif
      endif
     !
  ENDDO records1
     !
  print*,'num=',num
  CALL zaiocl(lua1)

!
!- second records
!
  level1=0

  iheada=ihead
  iheads=ihead

  DEALLOCATE(T_a)
  ALLOCATE(T_a(idm,jdm))

  ! Loop over the 2nd records (later in time)
  mark=knt1
  knt=knt-1
  records2: DO while(mark<knt)
     mark=mark+1
     !
     ! Read ifile2
     call get_abfld1(lua2,ip,T_a,m,n,labela, &
          level1,'./temp2/'//ifile2,iheada,time)

     if (time==g_time(mark)) then
        if (g_time(mark)>=maxval(g_time(:knt1))) then
          CALL zaiowr(T_a,imsk,.FALSE., fldmin,fldmax, lua,.FALSE.)
          WRITE(lub,'(A,'' date,span,range = '',F10.2,'' 0 '',1P2E16.7)') &
          &       labela(1:10),time,fldmin,fldmax
          num=num+1
     print*,'2nd loop, num=',num
        endif
     endif
  ENDDO records2
  !
  CALL zaiocl(lua2)

  CLOSE(lub)
  CALL zaiocl(lua)
  !
  DEALLOCATE(T_a)
  DEALLOCATE(ip)
  DEALLOCATE(imsk)
  DEALLOCATE(preambl)
  !
  !
END PROGRAM concate_forcing
