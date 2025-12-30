!-----------------------------------------------------------------------
      SUBROUTINE INTGRIB1(K1F,KGDS1F,K2F,KGDS2F,MRL,MRO,                &
     &                    IV,IP,IPOPT,KGDS1,K1,IB1,L1,F1,G1,KGDS2,K2,   &
     &                    IB2,L2,F2,G2,IRET)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:    INTGRIB1    INTERPOLATE FIELD
!   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
!
! ABSTRACT: INTERPOLATE FIELD.
!
! PROGRAM HISTORY LOG:
!   96-07-19  IREDELL
!
! USAGE:    CALL INTGRIB1(K1F,KGDS1F,K2F,KGDS2F,MRL,MRO,
!    &                    IV,IP,IPOPT,KGDS1,K1,IB1,L1,F1,G1,KGDS2,K2,
!    &                    IB2,L2,F2,G2,IRET)
!   INPUT ARGUMENTS:
!     K1F          INTEGER REGULARIZED INPUT DIMENSION
!     KGDS1F       INTEGER (200) REGULARIZED INPUT GRID PARAMETERS
!     K2F          INTEGER REGULARIZED OUTPUT DIMENSION
!     KGDS2F       INTEGER (200) REGULARIZED OUTPUT GRID PARAMETERS
!     MRL          INTEGER DIMENSION OF RLAT AND RLON
!     MRO          INTEGER DIMENSION OF CROT AND SROT
!     IV           INTEGER VECTOR FLAG
!     IP           INTEGER INTERPOLATION TYPE
!     IPOPT        INTEGER (20) INTERPOLATION OPTIONS
!     KGDS1        INTEGER (200) INPUT GRID PARAMETERS
!     K1           INTEGER INPUT DIMENSION
!     IB1          INTEGER INPUT BITMAP FLAG
!     L1           LOGICAL*1 (K1) INPUT BITMAP IF IB1=1
!     F1           REAL (K1) INPUT FIELD
!     G1           REAL (K1) INPUT Y-COMPONENT IF IV=1
!     KGDS2        INTEGER (200) OUTPUT GRID PARAMETERS
!     K2           INTEGER OUTPUT DIMENSION
!     IB2          INTEGER OUTPUT BITMAP FLAG
!     L2           LOGICAL*1 (K2) OUTPUT BITMAP
!     F2           REAL (K2) OUTPUT FIELD
!     G2           REAL (K2) OUTPUT Y-COMPONENT IF IV=1
!
! SUBPROGRAMS CALLED:
!   IPOLATES
!   IPOLATEV
!   IPXWAFS2
!   IPXWAFS3
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN
!
!$$$
      IMPLICIT NONE
      INTEGER :: K1,K2,K1F,K2F,MRL,MRO,IP,IB1,IV,KI,IB2,IRET,IB1F
      INTEGER :: IB2F
      INTEGER IPOPT(20)
      INTEGER KGDS1(200),KGDS2(200)
      LOGICAL*1 L1(K1),L2(K2)
      REAL F1(K1),F2(K2),G1(K1),G2(K2)
      INTEGER KGDS1F(200),KGDS2F(200)
      LOGICAL*1 L1F(K1F),L2F(K2F)
      REAL F1F(K1F),F2F(K2F),G1F(K1F),G2F(K2F)
      REAL RLAT(MRL),RLON(MRL),CROT(MRO),SROT(MRO)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  REGLR TO REGLR SCALAR
      BIGIF: IF(K1F.EQ.1.AND.K2F.EQ.1.AND.IV.EQ.0) THEN
        CALL IPOLATES(IP,IPOPT,KGDS1,KGDS2,K1,K2,1,IB1,L1,F1,           &
     &                KI,RLAT,RLON,IB2,L2,F2,IRET)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  IRREG TO REGLR SCALAR
      ELSEIF(K1F.NE.1.AND.K2F.EQ.1.AND.IV.EQ.0) THEN
         stop 10
        IF(IP.EQ.2) THEN
          CALL IPXWAFS3(1,K1,K1F,1,                                     &
     &                  KGDS1,IB1,L1,F1,KGDS1F,IB1F,L1F,F1F,IRET)
        ELSE
          CALL IPXWAFS2(1,K1,K1F,1,                                     &
     &                  KGDS1,IB1,L1,F1,KGDS1F,IB1F,L1F,F1F,IRET)
        ENDIF
        IF(IRET.EQ.0) THEN
          CALL IPOLATES(IP,IPOPT,KGDS1F,KGDS2,K1F,K2,1,IB1F,L1F,F1F,    &
     &                  KI,RLAT,RLON,IB2,L2,F2,IRET)
        ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  REGLR TO IRREG SCALAR
      ELSEIF(K1F.EQ.1.AND.K2F.NE.1.AND.IV.EQ.0) THEN
         stop 11
        CALL IPOLATES(IP,IPOPT,KGDS1,KGDS2F,K1,K2F,1,IB1,L1,F1,         &
     &                KI,RLAT,RLON,IB2F,L2F,F2F,IRET)
        IF(IRET.EQ.0) THEN
          IF(IP.EQ.2) THEN
            CALL IPXWAFS3(-1,K2,K2F,1,                                  &
     &                    KGDS2,IB2,L2,F2,KGDS2F,IB2F,L2F,F2F,IRET)
          ELSE
            CALL IPXWAFS2(-1,K2,K2F,1,                                  &
     &                    KGDS2,IB2,L2,F2,KGDS2F,IB2F,L2F,F2F,IRET)
          ENDIF
        ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  IRREG TO IRREG SCALAR
      ELSEIF(K1F.NE.1.AND.K2F.NE.1.AND.IV.EQ.0) THEN
         stop 12
        IF(IP.EQ.2) THEN
          CALL IPXWAFS3(1,K1,K1F,1,                                     &
     &                  KGDS1,IB1,L1,F1,KGDS1F,IB1F,L1F,F1F,IRET)
        ELSE
          CALL IPXWAFS2(1,K1,K1F,1,                                     &
     &                  KGDS1,IB1,L1,F1,KGDS1F,IB1F,L1F,F1F,IRET)
        ENDIF
        IF(IRET.EQ.0) THEN
          CALL IPOLATES(IP,IPOPT,KGDS1F,KGDS2F,K1F,K2F,1,IB1F,L1F,F1F,  &
     &                  KI,RLAT,RLON,IB2F,L2F,F2F,IRET)
          IF(IRET.EQ.0) THEN
            IF(IP.EQ.2) THEN
              CALL IPXWAFS3(-1,K2,K2F,1,                                &
     &                      KGDS2,IB2,L2,F2,KGDS2F,IB2F,L2F,F2F,IRET)
            ELSE
              CALL IPXWAFS2(-1,K2,K2F,1,                                &
     &                      KGDS2,IB2,L2,F2,KGDS2F,IB2F,L2F,F2F,IRET)
            ENDIF
          ENDIF
        ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  REGLR TO REGLR VECTOR
      ELSEIF(K1F.EQ.1.AND.K2F.EQ.1.AND.IV.NE.0) THEN
        CALL IPOLATEV(IP,IPOPT,KGDS1,KGDS2,K1,K2,1,IB1,L1,F1,G1,        &
     &                KI,RLAT,RLON,CROT,SROT,IB2,L2,F2,G2,IRET)
        IF(IRET.EQ.0.AND.KI.EQ.K2-1) THEN
          F2(K2)=0
          G2(K2)=0
        ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  IRREG TO REGLR VECTOR
      ELSEIF(K1F.NE.1.AND.K2F.EQ.1.AND.IV.NE.0) THEN
        IF(IP.EQ.2) THEN
          CALL IPXWAFS3(1,K1,K1F,1,                                     &
     &                  KGDS1,IB1,L1,F1,KGDS1F,IB1F,L1F,F1F,IRET)
          CALL IPXWAFS3(1,K1,K1F,1,                                     &
     &                  KGDS1,IB1,L1,G1,KGDS1F,IB1F,L1F,G1F,IRET)
        ELSE
          CALL IPXWAFS2(1,K1,K1F,1,                                     &
     &                  KGDS1,IB1,L1,F1,KGDS1F,IB1F,L1F,F1F,IRET)
          CALL IPXWAFS2(1,K1,K1F,1,                                     &
     &                  KGDS1,IB1,L1,G1,KGDS1F,IB1F,L1F,G1F,IRET)
        ENDIF
        IF(IRET.EQ.0) THEN
          CALL IPOLATEV(IP,IPOPT,KGDS1F,KGDS2,K1F,K2,1,                 &
     &                  IB1F,L1F,F1F,G1F,                               &
     &                  KI,RLAT,RLON,CROT,SROT,IB2,L2,F2,G2,IRET)
          IF(IRET.EQ.0.AND.KI.EQ.K2-1) THEN
            F2(K2)=0
            G2(K2)=0
          ENDIF
        ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  REGLR TO IRREG VECTOR
      ELSEIF(K1F.EQ.1.AND.K2F.NE.1.AND.IV.NE.0) THEN
         stop 13
        CALL IPOLATEV(IP,IPOPT,KGDS1,KGDS2F,K1,K2F,1,IB1,L1,F1,G1,      &
     &                KI,RLAT,RLON,CROT,SROT,IB2F,L2F,F2F,G2F,IRET)
        IF(IRET.EQ.0) THEN
          IF(IP.EQ.2) THEN
            CALL IPXWAFS3(-1,K2,K2F,1,                                  &
     &                    KGDS2,IB2,L2,F2,KGDS2F,IB2F,L2F,F2F,IRET)
            CALL IPXWAFS3(-1,K2,K2F,1,                                  &
     &                    KGDS2,IB2,L2,G2,KGDS2F,IB2F,L2F,G2F,IRET)
          ELSE
            CALL IPXWAFS2(-1,K2,K2F,1,                                  &
     &                    KGDS2,IB2,L2,F2,KGDS2F,IB2F,L2F,F2F,IRET)
            CALL IPXWAFS2(-1,K2,K2F,1,                                  &
     &                    KGDS2,IB2,L2,G2,KGDS2F,IB2F,L2F,G2F,IRET)
          ENDIF
        ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  IRREG TO IRREG VECTOR
      ELSEIF(K1F.NE.1.AND.K2F.NE.1.AND.IV.NE.0) THEN
         stop 14
        IF(IP.EQ.2) THEN
          CALL IPXWAFS3(1,K1,K1F,1,                                     &
     &                  KGDS1,IB1,L1,F1,KGDS1F,IB1F,L1F,F1F,IRET)
          CALL IPXWAFS3(1,K1,K1F,1,                                     &
     &                  KGDS1,IB1,L1,G1,KGDS1F,IB1F,L1F,G1F,IRET)
        ELSE
          CALL IPXWAFS2(1,K1,K1F,1,                                     &
     &                  KGDS1,IB1,L1,F1,KGDS1F,IB1F,L1F,F1F,IRET)
          CALL IPXWAFS2(1,K1,K1F,1,                                     &
     &                  KGDS1,IB1,L1,G1,KGDS1F,IB1F,L1F,G1F,IRET)
        ENDIF
        IF(IRET.EQ.0) THEN
          CALL IPOLATEV(IP,IPOPT,KGDS1F,KGDS2F,K1F,K2F,1,               &
     &                  IB1F,L1F,F1F,G1F,                               &
     &                  KI,RLAT,RLON,CROT,SROT,IB2F,L2F,F2F,G2F,IRET)
          IF(IRET.EQ.0) THEN
            IF(IP.EQ.2) THEN
              CALL IPXWAFS3(-1,K2,K2F,1,                                &
     &                      KGDS2,IB2,L2,F2,KGDS2F,IB2F,L2F,F2F,IRET)
              CALL IPXWAFS3(-1,K2,K2F,1,                                &
     &                      KGDS2,IB2,L2,G2,KGDS2F,IB2F,L2F,G2F,IRET)
            ELSE
              CALL IPXWAFS2(-1,K2,K2F,1,                                &
     &                      KGDS2,IB2,L2,F2,KGDS2F,IB2F,L2F,F2F,IRET)
              CALL IPXWAFS2(-1,K2,K2F,1,                                &
     &                      KGDS2,IB2,L2,G2,KGDS2F,IB2F,L2F,G2F,IRET)
            ENDIF
          ENDIF
        ENDIF
      ENDIF BIGIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END SUBROUTINE INTGRIB1
!-----------------------------------------------------------------------
      INTEGER FUNCTION LENGDSF(KGDS,KGDSF)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:    LENGDSF     RETURN THE LENGTH OF A FILLED GRID
!   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
!
! ABSTRACT: GIVEN A GRID DESCRIPTION SECTION (IN W3FI63 FORMAT),
!   RETURN THE GRID DESCRIPTION SECTION AND SIZE OF ITS REGULARIZED
!   COUNTERPART.  THAT IS, IF THE INPUT GRID IS REGULAR, THEN ITSELF
!   IS RETURNED ALONG WITH ITS GRID SIZE; HOWEVER IF THE INPUT GRID IS
!   ONLY QUASI-REGULAR (SUCH AS THE WAFS GRIDS), THEN ITS FILLED REGULAR
!   VERSION IS RETURNED ALONG WITH ITS FILLED GRID SIZE.
!
! PROGRAM HISTORY LOG:
!   96-07-19  IREDELL
!
! USAGE:    CALL LENGDSF(KGDS,KGDSF)
!   INPUT ARGUMENTS:
!     KGDS         INTEGER (200) GDS PARAMETERS IN W3FI63 FORMAT
!   OUTPUT ARGUMENTS:
!     KGDSF        INTEGER (200) REGULAR GDS PARAMETERS IN W3FI63 FORMAT
!     LENGDSF      INTEGER SIZE OF REGULARIZED GRID
!
! SUBPROGRAMS CALLED:
!   IPXWAFS
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN
!
!$$$
      IMPLICIT NONE
      INTEGER KGDS(200),KGDSF(200),DUM,IRET,DUMF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      BIGIF: IF(KGDS(1).EQ.201) THEN
        KGDSF=KGDS
        LENGDSF=KGDS(7)*KGDS(8)-KGDS(8)/2
      ELSEIF(KGDS(1).EQ.202) THEN
        KGDSF=KGDS
        LENGDSF=KGDS(7)*KGDS(8)
      ELSEIF(KGDS(19).EQ.0.AND.KGDS(20).NE.255) THEN
        CALL IPXWAFS(1,1,1,0,KGDS,DUM,KGDSF,DUMF,IRET)
        IF(IRET.EQ.0) THEN
          LENGDSF=KGDSF(2)*KGDSF(3)
        ELSE
          LENGDSF=0
        ENDIF
      ELSE
        KGDSF=KGDS
        LENGDSF=KGDS(2)*KGDS(3)
      ENDIF BIGIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END FUNCTION LENGDSF
