
      PROGRAM HWRF_FINAL_MERGE

!??????????????????????????????????????????????????????????
!     
! ABSTRACT: CREATE 6x DEGREE2 HIGH RESOLUTION DOMAIN
!     
!     DECLARE VARIABLES
!
!      IMPLICIT NONE

      INTEGER I,J,K,NX,NY,NZ,NST,IFLAG
      REAL(8) xxi,yyi
!
      PARAMETER (NST=5)
!      PARAMETER (IX=117,IY=225)
!      PARAMETER (NX=215,NY=431,NZ=42,NST=5)
!      PARAMETER (JX=393,JY=735)       ! fixed for 9 km resolution
      PARAMETER (GAMMA=6.5E-3,G=9.8,Rd=287.05,D608=0.608)
      PARAMETER (Cp=1004.)

! Variables on new outer nest hybrid coordinate  (GFS data)

      REAL(4) DLMD3,DPHD3,PT3,PDTOP3              ! use the new inner nest data
      REAL(4) WBD3,SBD3,CLON3,CLAT3

      REAL(4), ALLOCATABLE :: HLON3(:,:),HLAT3(:,:)
      REAL(4), ALLOCATABLE :: VLON3(:,:),VLAT3(:,:)
      REAL(4), ALLOCATABLE :: ETA1(:),ETA2(:)

      REAL(4), ALLOCATABLE :: T3(:,:,:),Q3(:,:,:)
      REAL(4), ALLOCATABLE :: U3(:,:,:),V3(:,:,:)
      REAL(4), ALLOCATABLE :: Z3(:,:,:),P3(:,:,:)
      REAL(4), ALLOCATABLE :: PD3(:,:)

! Variables old outer nest (6 hour WRF FORECAST)

      REAL(4) DLMD1,DPHD1,PT1,PDTOP1
      REAL(4) WBD1,SBD1,CLON1,CLAT1

      REAL(4), ALLOCATABLE :: HLON1(:,:),HLAT1(:,:)
      REAL(4), ALLOCATABLE :: VLON1(:,:),VLAT1(:,:)

      REAL(4), ALLOCATABLE :: T1(:,:,:),Q1(:,:,:)
      REAL(4), ALLOCATABLE :: U1(:,:,:),V1(:,:,:) 
      REAL(4), ALLOCATABLE :: Z1(:,:,:),P1(:,:,:)
      REAL(4), ALLOCATABLE :: PD1(:,:)

! Variables from 2x data

      REAL(4) DLMD2,DPHD2,PT2,PDTOP2
      REAL(4) WBD2,SBD2,CLON2,CLAT2

      REAL(4), ALLOCATABLE :: HLON2(:,:),HLAT2(:,:)
      REAL(4), ALLOCATABLE :: VLON2(:,:),VLAT2(:,:)

      REAL(4), ALLOCATABLE :: T2(:,:,:),Q2(:,:,:)
      REAL(4), ALLOCATABLE :: U2(:,:,:),V2(:,:,:)
      REAL(4), ALLOCATABLE :: Z2(:,:,:),P2(:,:,:)
      REAL(4), ALLOCATABLE :: PD2(:,:) 
   
!!!!!!!!!!!!!!!!11


      COEF1=Rd/Cp
      COEF3=Rd*GAMMA/G
      COEF2=1./COEF3

      GRD=G/Rd

      pi=4.*atan(1.)
      pi_deg=180./pi
      pi180=1./pi_deg

      DIST1=6.371E3*pi180

      READ(5,*)ITIM,KSTORM

      print*,'ITIM,KSTORM=',ITIM,KSTORM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! READ NEW OUTER NEST DATA    ! new data

      IUNIT=41
      READ(IUNIT) JX,JY,KZ
      REWIND(IUNIT)

      JX1=JX+1
      JY1=JY+1
      KZ1=KZ+1

      NZ=KZ

      print*,'JX,JY,KZ=',JX,JY,KZ

      ALLOCATE ( HLON3(JX,JY),HLAT3(JX,JY) )
      ALLOCATE ( VLON3(JX,JY),VLAT3(JX,JY) )
      ALLOCATE ( ETA1(KZ1),ETA2(KZ1) )

      ALLOCATE ( T3(JX,JY,KZ),Q3(JX,JY,KZ) )
      ALLOCATE ( U3(JX,JY,KZ),V3(JX,JY,KZ) )
      ALLOCATE ( Z3(JX,JY,KZ1),P3(JX,JY,KZ1) )
      ALLOCATE ( PD3(JX,JY) )

      ALLOCATE ( T1(JX,JY,KZ),Q1(JX,JY,KZ) )
      ALLOCATE ( U1(JX,JY,KZ),V1(JX,JY,KZ) )
      ALLOCATE ( Z1(JX,JY,KZ1),P1(JX,JY,KZ1) )
      ALLOCATE ( PD1(JX,JY) )
 
      ALLOCATE ( T2(JX,JY,KZ),Q2(JX,JY,KZ) )
      ALLOCATE ( U2(JX,JY,KZ),V2(JX,JY,KZ) )
      ALLOCATE ( Z2(JX,JY,KZ1),P2(JX,JY,KZ1) )
      ALLOCATE ( PD2(JX,JY) )

      T1=0.
      Q1=0.
      U1=0.
      V1=0.
      Z1=0.
      P1=0.
      PD1=0.


      DO KST=1,KSTORM

!      IUNIT=40+ITIM

      IUNIT=40+KST

      READ(IUNIT) JX,JY,KZ
      READ(IUNIT) DLMD3,DPHD3,CLON3,CLAT3
      READ(IUNIT) PT3,PDTOP3
      READ(IUNIT) T3
      READ(IUNIT) Q3
      READ(IUNIT) U3
      READ(IUNIT) V3
      READ(IUNIT) Z3
      READ(IUNIT) HLON3,HLAT3
      READ(IUNIT) P3
      READ(IUNIT) PD3
      READ(IUNIT) ETA1
      READ(IUNIT) ETA2

      CLOSE(IUNIT)

      T1=T1+T3
      Q1=Q1+Q3
      U1=U1+U3
      V1=V1+V3
      Z1=Z1+Z3
      P1=P1+P3
      PD1=PD1+PD3

      END DO    ! KST


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! READ ORG PARENT DATA    ! new data

      IUNIT=40

      READ(IUNIT) JX,JY,KZ
      READ(IUNIT) ! DLMD2,DPHD2,CLON2,CLAT2
      READ(IUNIT) ! PT2,PDTOP2
      READ(IUNIT) T2
      READ(IUNIT) Q2
      READ(IUNIT) U2
      READ(IUNIT) V2
      READ(IUNIT) Z2
      READ(IUNIT) ! HLON2,HLAT2
      READ(IUNIT) P2
      READ(IUNIT) PD2
      READ(IUNIT) ! ETA1
      READ(IUNIT) ! ETA2

      CLOSE(IUNIT)

      T3=T1-(KSTORM-1)*T2
      Q3=Q1-(KSTORM-1)*Q2
      U3=U1-(KSTORM-1)*U2
      V3=V1-(KSTORM-1)*V2
      Z3=Z1-(KSTORM-1)*Z2
      P3=P1-(KSTORM-1)*P2
      PD3=PD1-(KSTORM-1)*PD2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IUNIT=50+ITIM

      WRITE(IUNIT) JX,JY,KZ
      WRITE(IUNIT) DLMD3,DPHD3,CLON3,CLAT3
      WRITE(IUNIT) PT3,PDTOP3
      WRITE(IUNIT) T3
      WRITE(IUNIT) Q3
      WRITE(IUNIT) U3
      WRITE(IUNIT) V3
      WRITE(IUNIT) Z3
      WRITE(IUNIT) HLON3,HLAT3
      WRITE(IUNIT) P3
      WRITE(IUNIT) PD3
      WRITE(IUNIT) ETA1
      WRITE(IUNIT) ETA2

      CLOSE(IUNIT)

!      WRITE(61)((SLP3(I,J),I=1,JX),J=1,JY,2)
!      DO K=1,KZ+1
!        WRITE(61)((Z3(I,J,K),I=1,JX),J=1,JY,2)
!      END DO
!      DO K=1,KZ+1
!        WRITE(61)((P3(I,J,K),I=1,JX),J=1,JY,2)
!      END DO
!      DO K=1,KZ
!        WRITE(61)((T3(I,J,K),I=1,JX),J=1,JY,2)
!      END DO
!      DO K=1,KZ
!        WRITE(61)((Q3(I,J,K),I=1,JX),J=1,JY,2)
!      END DO
!      DO K=1,KZ
!        WRITE(61)((U3(I,J,K),I=1,JX),J=1,JY,2)
!      END DO
!      DO K=1,KZ
!        WRITE(61)((V3(I,J,K),I=1,JX),J=1,JY,2)
!      END DO
       END

