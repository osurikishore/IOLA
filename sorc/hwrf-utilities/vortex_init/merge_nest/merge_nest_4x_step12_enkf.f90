      PROGRAM HWRF_MERGE_NEST_ENKF

!??????????????????????????????????????????????????????????
!     
! ABSTRACT: Create enemble domain used in covariance calucation
!     
! ORIGINAL AUTHOR: Mingjing Tong, NCEP/EMC, 2015
!
!     DECLARE VARIABLES
!
     IMPLICIT NONE

      INTEGER ITIM,IBGS,I,J,K,NX,NY,NZ,NZ1,IX,IY,IZ,IZ1
      INTEGER ICLAT,ICLON,I360,ILOC,JLOC,IUNIT
      REAL(4)  IVOBS,CLON0,CLAT0,ERR,DIF1
      REAL(4)  CLON_NHC,CLAT_NHC

! Variables from outer nest (d01)

      REAL(4) DLMD1,DPHD1,PT1,PDTOP1
      REAL(4) WBD1,SBD1,CLON1,CLAT1

      REAL(4), ALLOCATABLE :: T1(:,:,:),Q1(:,:,:)
      REAL(4), ALLOCATABLE :: U1(:,:,:),V1(:,:,:) 
      REAL(4), ALLOCATABLE :: Z1(:,:,:),P1(:,:,:)
      REAL(4), ALLOCATABLE :: HLON1(:,:),HLAT1(:,:),VLON1(:,:),VLAT1(:,:)
      REAL(4), ALLOCATABLE :: ETA1(:),ETA2(:)
      REAL(4), ALLOCATABLE :: PD1(:,:)
      REAL(4), ALLOCATABLE :: A101(:,:),B101(:,:)

! Variables from inner nest (d02)

      REAL(4) DLMD2,DPHD2,PT2,PDTOP2
      REAL(4) WBD2,SBD2,CLON2,CLAT2,wind_s
      REAL(4), ALLOCATABLE :: HLON2(:,:),HLAT2(:,:),VLON2(:,:),VLAT2(:,:)
      REAL(4), ALLOCATABLE :: T2(:,:,:),Q2(:,:,:)
      REAL(4), ALLOCATABLE :: U2(:,:,:),V2(:,:,:)
      REAL(4), ALLOCATABLE :: Z2(:,:,:),P2(:,:,:)
      REAL(4), ALLOCATABLE :: PD2(:,:),A102(:,:),B102(:,:),C102(:,:)
   
      CHARACTER*1 SN,EW

      READ(5,*)ITIM,IVOBS,IBGS,CLAT0,CLON0   ! CLAT0,CLON0 is the new domain center

! read in storm center

       read(11,11)ICLAT,SN,ICLON,EW
  11   format(33x,I3,A1,I5,A1)
       CLAT_NHC=ICLAT*0.1
       CLON_NHC=ICLON*0.1

       IF(SN.eq.'S')CLAT_NHC=-CLAT_NHC
       IF(EW.eq.'W')CLON_NHC=-CLON_NHC

      I360=180
      if(abs(CLON_NHC).gt.90.)then
         I360=360
      end if

      print*,'I360=',I360

      if(I360.eq.360) then
        IF(CLON0.GT.0.)CLON0=CLON0-360.
        if(CLON_NHC.gt.0.)CLON_NHC=CLON_NHC-360.
      endif

! READ PARENT DATA (d01)   ! 6 hour forecast data

      IUNIT=20+ITIM

      READ(IUNIT) NX,NY,NZ

      print*,'NX,NY,NZ=',NX,NY,NZ

      NZ1=NZ+1

      ALLOCATE ( T1(NX,NY,NZ),Q1(NX,NY,NZ) )
      ALLOCATE ( U1(NX,NY,NZ),V1(NX,NY,NZ) )
      ALLOCATE ( Z1(NX,NY,NZ1),P1(NX,NY,NZ1) )
      ALLOCATE ( HLON1(NX,NY),HLAT1(NX,NY),VLON1(NX,NY),VLAT1(NX,NY) )
      ALLOCATE ( PD1(NX,NY),A101(NX,NY),B101(NX,NY) )
      ALLOCATE ( ETA1(NZ1),ETA2(NZ1) )

      READ(IUNIT) DLMD1,DPHD1,CLON1,CLAT1
      READ(IUNIT) PT1,PDTOP1
      READ(IUNIT) T1
      READ(IUNIT) Q1
      READ(IUNIT) U1
      READ(IUNIT) V1
      READ(IUNIT) Z1
      READ(IUNIT) HLON1,HLAT1,VLON1,VLAT1
      READ(IUNIT) P1
      READ(IUNIT) PD1
      READ(IUNIT) ETA1
      READ(IUNIT) ETA2

      print*,'reading A101'
      READ(IUNIT) A101                        ! A101 = land sea mask, B101 = ZNT
      READ(IUNIT) B101
      print*,'finishing reading A101'

      if(I360.eq.360) then
        DO J=1,NY
        DO I=1,NX
          IF(HLON1(I,J).GT.0.)HLON1(I,J)=HLON1(I,J)-360.
          IF(VLON1(I,J).GT.0.)VLON1(I,J)=VLON1(I,J)-360.
        END DO
        END DO
      endif

      REWIND(IUNIT)
      WRITE(IUNIT) NX,NY,NZ,I360
      WRITE(IUNIT) DLMD1,DPHD1,CLON1,CLAT1
      WRITE(IUNIT) PT1,PDTOP1
      WRITE(IUNIT) T1
      WRITE(IUNIT) Q1
      WRITE(IUNIT) U1
      WRITE(IUNIT) V1
      WRITE(IUNIT) Z1
      WRITE(IUNIT) HLON1,HLAT1,VLON1,VLAT1
      WRITE(IUNIT) P1
      WRITE(IUNIT) PD1
      WRITE(IUNIT) ETA1
      WRITE(IUNIT) ETA2
      WRITE(IUNIT) A101                        ! A101 = land sea mask, B101 = ZNT
      WRITE(IUNIT) B101
 
      CLOSE(IUNIT)

      print*,'read in d01, HLON1,HLAT1,VLON1,VLAT1=',     &
              HLON1(1,1),HLAT1(1,1),VLON1(1,1),VLAT1(1,1)

      DEALLOCATE(T1,Q1,U1,V1)
      DEALLOCATE(Z1,P1)
      DEALLOCATE(PD1,A101,B101)
      DEALLOCATE(ETA1,ETA2)
!
! READ INNER NEST DATA (d02)   ! 6 hour forecast data

      IUNIT=30+ITIM

      READ(IUNIT) IX,IY,IZ              ! IZ==NZ

      print*,'IX,IY,IZ=',IX,IY,IZ

      IZ1=IZ+1

      ALLOCATE ( HLON2(IX,IY),HLAT2(IX,IY) )
      ALLOCATE ( VLON2(IX,IY),VLAT2(IX,IY) ) 
      ALLOCATE ( T2(IX,IY,IZ),Q2(IX,IY,IZ) )
      ALLOCATE ( U2(IX,IY,IZ),V2(IX,IY,IZ) )
      ALLOCATE ( Z2(IX,IY,IZ1),P2(IX,IY,IZ1) )
      ALLOCATE ( A102(IX,IY),B102(IX,IY),C102(IX,IY) )
      ALLOCATE ( PD2(IX,IY) )
      ALLOCATE ( ETA1(IZ1),ETA2(IZ1) )

      READ(IUNIT) DLMD2,DPHD2,CLON2,CLAT2
      READ(IUNIT) PT2,PDTOP2
      READ(IUNIT) T2
      READ(IUNIT) Q2
      READ(IUNIT) U2
      READ(IUNIT) V2
      READ(IUNIT) Z2
      READ(IUNIT) HLON2,HLAT2,VLON2,VLAT2
      READ(IUNIT) P2
      READ(IUNIT) PD2
      READ(IUNIT) ETA1
      READ(IUNIT) ETA2

      READ(IUNIT) A102
      READ(IUNIT) B102
      READ(IUNIT) C102

      if(I360.eq.360) then
        IF(CLON2.GT.0.)CLON2=CLON2-360.
        DO J=1,IY
        DO I=1,IX
          IF(HLON2(I,J).GT.0.)HLON2(I,J)=HLON2(I,J)-360.
          IF(VLON2(I,J).GT.0.)VLON2(I,J)=VLON2(I,J)-360.
        END DO
        END DO
      endif

      ERR=1.e20
      DO J=1,NY
      DO I=1,NX
        DIF1=abs(HLON2(1,1)-HLON1(I,J))+abs(HLAT2(1,1)-HLAT1(I,J))
        IF(DIF1.LT.ERR)THEN
          ILOC=I
          JLOC=J
          ERR=DIF1
        END IF
      END DO
      END DO

      WBD1=-(NX-1)*DLMD1                 ! PARENT wbd
      SBD1=-(NY/2)*DPHD1       ! PARENT SBD

      WBD2= WBD1 + (ILOC -1)*2.*DLMD1 + MOD(JLOC+1,2)*DLMD1
      SBD2= SBD1 + (JLOC -1)*DPHD1

      IUNIT=50+ITIM

      WRITE(IUNIT) IX,IY,IZ,I360
      WRITE(IUNIT) DLMD2,DPHD2,CLON2,CLAT2
      WRITE(IUNIT) PT2,PDTOP2,WBD2,SBD2
      WRITE(IUNIT) T2
      WRITE(IUNIT) Q2
      WRITE(IUNIT) U2
      WRITE(IUNIT) V2
      WRITE(IUNIT) Z2
      WRITE(IUNIT) HLON2,HLAT2,VLON2,VLAT2
      WRITE(IUNIT) P2
      WRITE(IUNIT) PD2
      WRITE(IUNIT) ETA1
      WRITE(IUNIT) ETA2

      CLOSE(IUNIT)

      DO J=1,IY
      DO I=1,IX
        wind_s=sqrt(U2(I,J,1)**2+V2(I,J,1)**2)+1.E-10
! Modification in case wind10 > lowest model level wind, especially if the lowest
! model level is lower than 10 m
        C102(I,J)=min(1.0,C102(I,J)/wind_s)
      END DO
      END DO

      IUNIT=60+ITIM

      WRITE(IUNIT) IX,IY,I360
      WRITE(IUNIT) HLON2,HLAT2,VLON2,VLAT2
      WRITE(IUNIT) A102
      WRITE(IUNIT) B102
      WRITE(IUNIT) C102

      CLOSE(IUNIT)

      DEALLOCATE(HLON1,HLAT1,VLON1,VLAT1)
      DEALLOCATE(HLON2,HLAT2,VLON2,VLAT2)
      DEALLOCATE(T2,Q2,U2,V2)
      DEALLOCATE(Z2,P2,PD2)
      DEALLOCATE(A102,B102,C102,ETA1,ETA2)

      END

