      PROGRAM read_grid
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C ******************IMPORTANT ******************
C Made the following change on Sep 21 2006
C IGDS(03)=0 to IGDS(03)=203
c This is to reflect the new addition to the GRIB table ON388 which is:
c 203   Curivilinear Orthogonal Grid
C Made the following change on Mar 21 2007
C IGDS(03)=203 to IGDS(03)=204
c This is to reflect the new addition to the GRIB table ON388 which is:
c 204   Curivilinear Orthogonal Grid
c TO add subregions gridno got rid of :if gridno=120: statement from gds section 


C     The regional.grid file is used for HYCOM 
C
C
      INTEGER    IWI,JWI
      REAL*4     XFIN,YFIN,DX0,DY0
      REAL*4     XeIN,YeIN,dxin,dyin


      PARAMETER (DX0=1.0/30.0, DY0=1.0/30.0)
C
C     BATHYMETRY ARRAYS.
C
      INTEGER, ALLOCATABLE :: IP(:,:)
      REAL*4,  ALLOCATABLE :: PLON(:,:),PLAT(:,:)
      REAL*4,  ALLOCATABLE :: ULON(:,:),ULAT(:,:)
      REAL*4,  ALLOCATABLE :: VLON(:,:),VLAT(:,:)


      REAL*4  PLATIJ,YFMIN,YFMAX,XAMAX,XAMIN,YAMAX,YAMIN
      CHARACTER PREAMBL(5)*79
      CHARACTER*80 header_old
      character*80 cregion
      character*80 title
C
C     INTERPOLATION ARRAYS.
C
      INTEGER IBD(4)
C
C     SAVE  /NPROCS/
C
      REAL*4          COAST
      INTEGER         INTERP

      integer         periodic
C
      EXTERNAL LINEAR, AVERMS,MINMAX
C
      REAL*4     ZERO,RADIAN
      PARAMETER (ZERO=0.0, RADIAN=57.2957795)
C
      CHARACTER*80 CLINE
      INTEGER I,II,IWIX,J,L,LENGTH
      REAL*4  XFD,YFD,DXD,DYD,BLAND
      REAL*4  HMINA,HMINB,HMAXA,HMAXB
C
C --- MODEL ARRAYS.
C
      CALL XCSPMD  !define idm,jdm
      ALLOCATE( PLON(IDM,JDM) )
      ALLOCATE( PLAT(IDM,JDM) )
      ALLOCATE( ULON(IDM,JDM) )
      ALLOCATE( ULAT(IDM,JDM) )
      ALLOCATE( VLON(IDM,JDM) )
      ALLOCATE( VLAT(IDM,JDM) )

C
C     NAMELIST INPUT.
C
      CALL ZHOPEN(6, 'FORMATTED', 'UNKNOWN', 0)
C
C     Read from the NAMELIST file here

C
C     GRID INPUT.
C
      CALL ZAIOST
C
      CALL ZHOPNC(21, 'regional.grid.b', 'FORMATTED', 'OLD', 0)
      CALL ZAIOPF('regional.grid.a', 'OLD', 21)
C
      READ(21,*) ! skip idm
      READ(21,*) ! skip jdm
      READ(21,*) ! skip mapflg
      READ(21,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
      CALL ZAIORD(PLON,IP,.FALSE., HMINA,HMAXA, 21)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (plon):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF
C skip qlat, qlon
      READ(21,*) ! skip qlat
      READ(21,*) ! skip qlon
      CALL ZAIOSK (21)
      CALL ZAIOSK (21)
C
      READ(21,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
      CALL ZAIORD(PLAT,IP,.FALSE., HMINA,HMAXA, 21)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (plat):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      READ(21,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
      CALL ZAIORD(ULON,IP,.FALSE., HMINA,HMAXA, 21)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (ulon):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF

C
      READ(21,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
      CALL ZAIORD(ULAT,IP,.FALSE., HMINA,HMAXA, 21)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (ulat):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF

c
      READ(21,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
      CALL ZAIORD(VLON,IP,.FALSE., HMINA,HMAXA, 21)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (vlon):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF

C
      READ(21,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
      CALL ZAIORD(VLAT,IP,.FALSE., HMINA,HMAXA, 21)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (vlat):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF

      CLOSE(UNIT=21)
      CALL ZAIOCL(21)
C

      call grib_nomads(plon,plat,ulon,ulat,vlon,vlat,idm,jdm)
      stop
      end



      subroutine  grib_nomads(plon,plat,ulon,ulat,vlon,vlat,iim,jjm)
      
      real plon(iim,jjm),plat(iim,jjm),ulon(iim,jjm),ulat(iim,jjm)
     *     ,vlon(iim,jjm),vlat(iim,jjm)
      real bnv,gds4
      logical run_u,run_v,run_p
      data   run_u/.true. /, run_v/.true. /, run_p
     &  /.true. /
      integer yr,rmm,dd,fhr,hh,igridno
      real idlat,idlon

      write(6,*) "Reading Input "
      read(5,*) yr,rmm,dd,fhr,hh,idlat,idlon,igridno
      write(6,*) yr,rmm,dd,fhr,hh

C Changed the table # from 1 to 129: March 2007
C Changed plat plon parameter # to 207 208
C Changed ulat ulon parameter # to 192 193
C Changed vlat vlon parameter # to 188 189

      if (run_p) then
      call baopen(51,'pgrib',IRET)
      write(6,'(a25,i5)') 'after baopen iret = ',iret
      call parm_def(plon,208,51,129,1,1,3,0.,iim,jjm,
     *              yr,rmm,dd,fhr,hh,idlat,idlon,igridno)
      call parm_def(plat,207,51,129,1,1,3,0.,iim,jjm,
     *              yr,rmm,dd,fhr,hh,idlat,idlon,igridno)
      call baclose(iu,iret)
      endif

      if (run_u) then
      call baopen(52,'ugrib',IRET)
      write(6,'(a25,i5)') 'after baopen iret = ',iret
      call parm_def(ulon,193,52,129,1,1,3,0.,iim,jjm,
     *              yr,rmm,dd,fhr,hh,idlat,idlon,igridno)
      call parm_def(ulat,192,52,129,1,1,3,0.,iim,jjm,
     *              yr,rmm,dd,fhr,hh,idlat,idlon,igridno)
      call baclose(iu,iret)
      endif

      if (run_v) then
      call baopen(53,'vgrib',IRET)
      write(6,'(a25,i5)') 'after baopen iret = ',iret
      call parm_def(vlon,189,53,129,1,1,3,0.,iim,jjm,
     *              yr,rmm,dd,fhr,hh,idlat,idlon,igridno)
      call parm_def(vlat,188,53,129,1,1,3,0.,iim,jjm,
     *              yr,rmm,dd,fhr,hh,idlat,idlon,igridno)
      call baclose(iu,iret)
      endif

        return
        end
c--------------------------------------------------------------

      SUBROUTINE PARM_DEF(work,parmet2,iu,
     * pdstbver,level,value,dec,gds4,iim,jjm,
     * yr,rmm,dd,fhr,hh,idlat,idlon,igridno)

      integer im,jm
      real anv,bnv,gds4
      integer kpds(40),kgds(22),ibds(9),nf,nbm,parmet2,imjm
      integer pdstbver,dec,iim0,jjm0,iim1,jjm1,iiijjj
      integer imin,imax,jmin,jmax,ij
      real ilat0,ilat1,ilon0,ilon1
      real work(iim,jjm),maxx,minn 
      real, allocatable, dimension (:,:) ::arryn
      real, allocatable, dimension (:) ::arr1d
      logical*1, allocatable, dimension (:) ::LB1D

      read(22,*) iim0,jjm0,iii,jjj
      write(6,*) "From 22 ", iim0,jjm0,iii,jjj
      rewind(22)
      iim1=iim0+iii-1
      jjm1=jjm0+jjj-1
      write(*,*) 'iim1 and jjm1 are ',iim1,jjm1
      allocate (arryn(iii,jjj))

      iiijjj=iii*jjj
      allocate ( arr1d(iiijjj) )
      allocate ( LB1D(iiijjj) )
      idlat=idlat*1000
      idlon=idlon*1000
      write(6,*)"idlat,idlon",idlat,idlon

      ilat0=0.
      ilat1=0.
      ilon0=0.
      ilon1=0.
      ii=0
      jold=0
      jj=0
         do j=1,jjm
            do i=1,iim
              if ((i.ge.iim0).and.(i.le.iim1).and.(j.ge.jjm0)
     *            .and.(j.le.jjm1)) then
                 if(jold.lt.j) then
                   jold=j
                   jj=jj+1
                   ii=1
                 else
                   ii=ii+1
                 endif
                 arryn(ii,jj)=work(i,j)
              endif
            enddo
           enddo

      bnv=2.**100
         ij=0
         do j=1,jjj
            do i=1,iii
             ij=ij+1
                 arr1d(ij)=arryn(i,j)
                  if(arr1d(ij).eq.bnv) then
                     LB1D(ij)=.false.
                  else
                     LB1D(ij)=.true.
                  endif
            enddo
         enddo


       call pdsgds(yr,rmm,dd,fhr,
     1           hh,parmet2,pdstbver,level,value,dec,
     2           ilat0,ilat1,ilon0,ilon1,idlon,idlat,iii,jjj,
     3           igridno,gds4,KPDS,KGDS)

        call putgb(iu,iiijjj,KPDS,KGDS,LB1D,arr1d,iret)
        write(6,'(a25,i5)') 'after putgb iret = ',iret

      deallocate ( arryn )
      deallocate ( arr1d )
      deallocate ( LB1D )
         return
         end


      SUBROUTINE PDSGDS (iyr,imo,ida,fhr,tau,parm,pdstbver,
     2 level,value,dec,ilat0,ilat1,ilon0,ilon1,idlon,idlat,
     1 im,jm,igridno,gds4,KPDS,KGDS)
      integer KPDS(40),KGDS(22)
      integer fhr,tau,parm,igridno
      integer im,jm,dec,pdstbver
      real gds4,idlat,idlon 
      real ilat0,ilat1,ilon0,ilon1
      integer yr,rmm,dd,hh
c
C Changed the process ID number from 45 to 47: March 2007
      write(6,*) "pdstbver is ",pdstbver

      do i=1,40
       kpds(i)=0
      enddo
c
C     KPDS     - ARRAY CONTAINING PDS ELEMENTS.  (EDITION 1)
      kpds(1) =   7    ! ID OF CENTER
      kpds(2) =  47    ! GENERATING PROCESS ID NUMBER
      kpds(3) = igridno ! GRID DEFINITION
      kpds(4) = 192    ! GDS/BMS FLAG (RIGHT ADJ COPY OF OCTET 8)
      kpds(5) = parm   ! INDICATOR OF PARAMETER
      kpds(6) = level  ! TYPE OF LEVEL
      kpds(7) = value  ! HEIGHT/PRESSURE , ETC OF LEVEL
      kpds(8) = iyr    ! YEAR of century
      kpds(9) = imo    ! MONTH OF YEAR
      kpds(10)= ida    ! DAY OF MONTH
      kpds(11)= tau    ! HOUR OF DAY
      kpds(12)=   0    ! MINUTE OF HOUR
      kpds(13)=   1    ! INDICATOR OF FORECAST TIME UNIT
      kpds(14)= fhr    ! TIME RANGE 1
      kpds(15)=   0    ! TIME RANGE 2
      kpds(16)=   0    ! TIME RANGE FLAG
      kpds(17)=   0    ! NUMBER INCLUDED IN AVERAGE
      kpds(18)=   1    ! VERSION NR OF GRIB SPECIFICATION
      kpds(19)= pdstbver !VERSION NR OF PARAMETER TABLE
      kpds(20)=   0    ! NR MISSING FROM AVERAGE/ACCUMULATION
      kpds(21)=  21    ! CENTURY OF REFERENCE TIME OF DATA
      kpds(22)=  dec   ! UNITS DECIMAL SCALE FACTOR
      kpds(23)=   0    ! SUBCENTER NUMBER
      kpds(24)=   0    ! PDS BYTE 29, FOR NMC ENSEMBLE PRODUCTS
C                 128 IF FORECAST FIELD ERROR
C                  64 IF BIAS CORRECTED FCST FIELD
C                  32 IF SMOOTHED FIELD
C                  WARNING: CAN BE COMBINATION OF MORE THAN 1
      kpds(25)=   0    !  - PDS BYTE 30, NOT USED
C       (26-35)  - RESERVED
C       (36-N)   - CONSECUTIVE BYTES EXTRACTED FROM PROGRAM
C                  DEFINITION SECTION (PDS) OF GRIB MESSAGE
c

      print*,'first kpds:',kpds
C Changed kgds(1) from 203 204 : March 20 2007
C     KGDS     - ARRAY CONTAINING GDS ELEMENTS.
      kgds(1) =   204  ! DATA REPRESENTATION TYPE
      kgds(2) =    im  ! N(I) NR POINTS ON LATITUDE CIRCLE
      kgds(3) =    jm  ! N(J) NR POINTS ON LONGITUDE MERIDIAN
      kgds(4) = nint(gds4)*1000  ! LA(1) LATITUDE OF ORIGIN  (26.5 N)
      kgds(5) = ilon0  ! LO(1) LONGITUDE OF ORIGIN (83.0 W)
      kgds(6) =     8  ! RESOLUTION FLAG (RIGHT ADJ COPY OF OCTET 17)
      kgds(7) = ilat1  ! LA(2) LATITUDE OF EXTREME POINT (47.2 N)
      kgds(8) = ilon1  ! LO(2) LONGITUDE OF EXTREME POINT (50.0 N)
      kgds(9) =     0  ! DI LATITUDINAL DIRECTION OF INCREMENT
      kgds(10)=     0  ! DJ LONGITUDINAL DIRECTION INCREMENT
      kgds(11)=    64  ! SCANNING MODE FLAG (RIGHT ADJ COPY OF OCTET 28)
      kgds(19)=     0  ! NUMBER OF VERTICAL COORDINATE PARAMETERS
      kgds(20)=   255  ! OCTET NUMBER OF THE LIST OF VERTICAL COORDINATE
C                          PARAMETERS
C                        OR
C                         OCTET NUMBER OF THE LIST OF NUMBERS OF POINTS
C                          IN EACH ROW
C                        OR
C                         255 IF NEITHER ARE PRESENT
      kgds(21)=     0  ! FOR GRIDS WITH PL, NUMBER OF POINTS IN GRID
      kgds(22)=     0  ! NUMBER OF WORDS IN EACH ROW
      do i=12,18
        kgds(i) = 0    ! 12-18 not used
      enddo
      print*,'kpds:',kpds
c
       print*,'KGDS:',kgds
        return
      end

