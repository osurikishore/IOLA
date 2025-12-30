      SUBROUTINE COF2GRD (LUN1,NC,JROMB,JCAP,KMAX,
     & IM,JM,LM,KB,KBETA,IMJM,JMP1,IMM1,IMJMM1,KHL00,KHH00,KNE,
     & KNW,KSW,KSE,KSL,KSLM1,LMM1,LMP1,IMT,JMT,IMAX,JMAX,
     & use_gfs_nemsio,use_gfs_netcdf,alat,alon,ib,jb,nxg,nyg)
C     
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C     .      .    .                                       .
C     SUBPROGRAM:    COF2GRD     CONVERT ONE RECORD OF SIGMA COEFF FILE
C     TO LAT/LON GRID
C     PRGMMR: ROGERS           ORG: W/NP22     DATE: 99-01-28
C     
C     Biju Thomas 2020/07/02: ADDED OPTION TO READ NetCDF GFS(V16) ANALYSIS FILE 
C     IN ADDITION TO NEMSIO/SPECTRAL FILES
C
C     ABSTRACT: CONVERT SIGMA COEFFICIENT RECORD TO GRID SPACE USING
C     SPLIB ROUTINES. THESE ROUTINES WILL RETURN A GLOBAL
C     LAT/LON GRID WHOSE RESOLUTION IS DETERMINED BY THE 
C     NUMBER OF GRID POINTS. THEN, THE RELEVENT SUBSET FOR
C     WHICH WE HAVE HIGH-RES OROGRAPHY IS EXTRACTED (DIMENSION
C     OF BOTH THE EXTRACTED GRID AND GLOBAL GRID SET IN 
C     parmanl FILE)
C     
C     PROGRAM HISTORY LOG:
C     99-01-28  ROGERS
C     
C     USAGE:    CALL COF2GRD(LUN1,NC,KMAX,JROMB, 
C     JCAP,XGRID,PGRID)
C     
C     INPUT ARGUMENT LIST:
C     LUN1     - FORTRAN UNIT FOR SIGMA FILE
C     NC       - LENGTH OF SIGMA RECORD = RES+1*RES+2
C     JROMB    - SPECTRAL DOMAIN SHAPE (0 FOR TRIANGULAR, 
C     1 FOR RHOMBOIDAL)
C     JCAP     - SPECTRAL TRUNCATION
C     
C     OUTPUT FILES:
C     KMAX     - NUMBER OF SIGMA LEVELS IN GLOBAL MODEL
C     XGRID    - ARRAY holding the IMAX x JMAX grids at KMAX levels
C     FOR FOUR VARIABLES:  1-Tv  2-U  3-V  4-q
C     PGRID    - ARRAY holding the IMAX x JMAX grids of Z* and p*
C     
C     ATTRIBUTES:
C     LANGUAGE: FORTRAN-90
C     MACHINE: CRAY C-90
C     
C$$$  

      USE SP_GRID_MODULE, ONLY: TGRID,UGRID,VGRID,QGRID,CWMGRID,
     &     PRESGRID,PINTGRID,PGRID, SPNLEVS, REAL_32
      USE SIGIO_MODULE
      use nemsio_module
      use read_gmv16_ncf
      IMPLICIT NONE
      TYPE(sigio_head):: HEAD
      TYPE(sigio_data):: DATA
      type(nemsio_gfile) :: gfile

      INTEGER :: NGUSED,NGBD,LTBGRI,LBDIM
      INTEGER :: I,J,IRET,L,IER,K,LLL,JCAP,KMAX,JROMB,LUN1,NC
      INTEGER :: LUN1HOLD,nthreads

      REAL    :: POLEI,POLEJ,XMESHL,ALONVT
      REAL    :: TPH0D,TLM0D,WBD,SBD

      LOGICAL :: POLA,NORTH
      logical :: use_gfs_nemsio, use_gfs_netcdf
      integer, parameter :: nnvcoord=2, nidvc=2, nidsl=1 ! For use_gfs_netcdf 
      integer :: ntrac,nframe,fieldsize,ib,jb,ic,jc,inemsio,jnemsio
      integer :: nxg,nyg,nzg,nzhg,atttype,attlen
      integer :: ii,jj
      real(4),allocatable :: fis(:),fit(:),vcoord(:,:,:)
      real,allocatable:: lat(:,:),lon(:,:)
      real,allocatable:: fld2d(:,:),fld3d(:,:,:)
      real,allocatable:: fld3d_cl(:,:,:), fld3d_gr(:,:,:), 
     &        fld3d_ic(:,:,:), fld3d_rw(:,:,:), fld3d_sn(:,:,:)
      real,allocatable:: pfull(:), phalf(:), ak(:), bk(:), abk(:,:)
      real :: dcenter_lat,dcenter_lon,r1,r2
      real :: alat(imax,jmax),alon(imax,jmax)

      INCLUDE "parmlbc"
C     
C     XMESHL = EAST-WEST GRID INCREMENT
C     POLEI  = SOUTH-NORTH GRID INCREMENT
C     POLEJ  = WESTERN BOUNDARY OF LAT/LON GRID
C     ALONVT = SOUTHERN BOUNDARY OF LAT/LON GRID
C     
      CHARACTER HOLDFIL*80,FILENAME*80
C     
      REAL(REAL_32),ALLOCATABLE,DIMENSION(:,:)::GRD
      REAL(REAL_32),ALLOCATABLE,DIMENSION(:,:)::MYGRD
      REAL(REAL_32),ALLOCATABLE,DIMENSION(:,:,:)::MYGRD2,PD
      REAL(REAL_32),ALLOCATABLE,DIMENSION(:)::DWORK,ZWORK
C     
      real :: PI

      print *,'entering cof2grd',jromb,jcap,imax,jmax,kmax
C     
      write(filename,633) LUN1
 633  format('fort.',I2.2)
      PI = 3.141592654

!     CALL sigio_srohdc(12,'gfsbc1',head,data,iret)
!     write(6,*) 'iret from sigio_srohdc: ', iret
!     write(6,*) 'head%levs: ', head%levs

      CLOSE(LUN1)

      write(6,*)'use_gfs_netcdf = ',use_gfs_netcdf,
     &                      ' use_gfs_nemsio = ', use_gfs_nemsio

      if ( use_gfs_nemsio ) then
!--- Initialize nemsio
        call nemsio_init(iret=iret)
!--- Opean a NEMSIO file
        call nemsio_open(gfile,trim(filename),'READ',iret=iret)
        call nemsio_getfilehead(gfile,iret=iret,jcap=head%jcap,
     &  dimx=head%lonb,dimy=head%latb,idsl=head%idsl,idvc=head%idvc,
     &  dimz=head%levs,ntrac=head%ntrac,nframe=nframe)
        head%nvcoord=2
        inemsio=head%lonb
        jnemsio=head%latb
        nxg=inemsio
        nyg=jnemsio
        nzg=head%levs
        call sigio_alhead(head,iret)

        fieldsize=(head%lonb+2*nframe)*(head%latb+2*nframe)
        allocate(fis(fieldsize),fit(fieldsize))
        allocate(lat(inemsio,jnemsio),lon(inemsio,jnemsio))
        allocate(vcoord(head%levs+1,3,2))

        call nemsio_getfilehead(gfile,iret=iret,vcoord=vcoord,
     &  lat=fis,lon=fit)
        do i=1,2
        do j=1,head%levs+1
          head%vcoord(j,i)=vcoord(j,i,1)
        enddo
        enddo
        do j=1,jnemsio
        do i=1,inemsio
          lat(i,j)=fis(i+(j-1)*head%lonb)
          lon(i,j)=fit(i+(j-1)*head%lonb)
        enddo
        enddo
        REWIND 45
        READ(45,*) dcenter_lat
        READ(45,*) dcenter_lon
        if (dcenter_lon .le. 0.) dcenter_lon=dcenter_lon+360.0
        do j=1,jnemsio-1
          r1=lat(1,j)
          r2=lat(1,j+1)
          if(dcenter_lat .ge. r2 .and. dcenter_lat .lt. r1) then
            jc=j
          endif
        enddo
!zhang latest test
        do i=1,inemsio-1
          r1=lon(i,1)
          r2=lon(i+1,1)
          if(dcenter_lon .ge. r1 .and. dcenter_lon .lt. r2) then
            ic=i
          endif
        enddo

        if(jc-int(jmax/2) .lt. 1 ) then
          jb=1
        endif

        if(jc-int(jmax/2) .ge. 1 .and. jc+int(jmax/2) .le. 
     &                                      jnemsio ) then
          jb=jc-int(jmax/2)
        endif

        if(jc+int(jmax/2) .gt. jnemsio) then
          jb=jnemsio-jmax
        endif

        ib=ic-int(imax/2)
        if(ib .lt. 1) ib=inemsio+ib
        print*, "ic=,jc=",ic,jc 
        print*, "imax=,jmax=",imax,jmax
        print*, "ib=,jb=",ib,jb 
        print*, "end=====",ib+imax,jb+jmax
!       print*, "lat==",lat(1,jb),lat(1,jc),lat(1,jb+jmax)
!       print*, "lon==",lon(ib,1),lon(ic,1),lon(ib+imax,1)
!       print*, "lat==",lat(imax,jb),lat(imax,jc),lat(imax,jb+jmax)
!       print*, "lon==",lon(ib,jmax),lon(ic,jmax),lon(ib+imax,jmax)
        do j=1,jmax
        do i=1,imax
          ii=i+ib-1
!         if(ii .gt. inemsio) ii=ii-inemsio
          alat(IMAX-I+1,JMAX-J+1)=lat(i,j-1+jb)
          if(ii .le. inemsio) then     
            alon(i,j)=360.-lon(ii,j)
          else
            alon(i,j)=-lon(ii-inemsio,j)
          endif
        enddo
        enddo

      else if(use_gfs_netcdf) then
        CALL get_dims(TRIM(filename), nxg, nyg, nzg, nzhg)
        allocate(lat(nxg, nyg),lon(nxg, nyg),vcoord(nzg+1,3,2))
        allocate(fld2d(nxg,nyg),fld3d(nxg,nyg,nzg))
        allocate(pfull(nzg), phalf(nzhg))
        CALL get_coord(TRIM(filename), nxg, nyg, nzg, nzhg, 'lon', 
     &       'lat', 'pfull', 'phalf', lon, lat, pfull, phalf)
        rewind 45
        read(45,*) dcenter_lat
        read(45,*) dcenter_lon
        if (dcenter_lon .le. 0.) dcenter_lon=dcenter_lon+360.0
        do j=1,nyg-1
          r1=lat(1,j)
          r2=lat(1,j+1)
          if(dcenter_lat .ge. r2 .and. dcenter_lat .lt. r1) then
            jc=j
          endif
        enddo
        do i=1,nxg-1
          r1=lon(i,1)
          r2=lon(i+1,1)
          if(dcenter_lon .ge. r1 .and. dcenter_lon .lt. r2) then
            ic=i
          endif
        enddo

        if(jc-int(jmax/2) .lt. 1 ) then
          jb=1
        endif

        if(jc-int(jmax/2) .ge. 1 .and. jc+int(jmax/2) .le.
     &                                      nyg ) then
          jb=jc-int(jmax/2)
        endif

        if(jc+int(jmax/2) .gt. nyg) then
          jb=nyg-jmax
        endif

        ib=ic-int(imax/2)
        if(ib .lt. 1) ib=nxg+ib
        do j=1,jmax
        do i=1,imax
          ii=i+ib-1
!         if(ii .gt. nxg) ii=ii-nxg
          alat(imax-i+1,jmax-j+1)=lat(i,j-1+jb)
          if(ii .le. nxg) then
            alon(i,j)=360.-lon(ii,j)
          else
            alon(i,j)=-lon(ii-nxg,j)
          endif
        enddo
        enddo

        call get_attrlen(TRIM(filename),'ak',atttype,attlen)
        allocate(ak(attlen),bk(attlen),abk(attlen,2))
        call get_attr(TRIM(filename),'ak',attlen,ak)
        call get_attr(TRIM(filename),'bk',attlen,bk)
!       reverse the order
        do l=1, attlen
          abk(l,1) = ak(attlen-l+1)
          abk(l,2) = bk(attlen-l+1)
        enddo
      else
        print*, "read in GFS data in sigio format"
        CALL sigio_srohdc(LUN1,trim(filename),head,data,iret)
        call summary()
      endif

      ALLOCATE(MYGRD(IMAX,JMAX))
      ALLOCATE(MYGRD2(IMAX,JMAX,nzg))
      ALLOCATE(PD(IMAX,JMAX,nzg))
      ALLOCATE(GRD(IMAX,JMAX))

      SPNLEVS=nzg
      ALLOCATE(TGRID(IMAX,JMAX,nzg))
      ALLOCATE(UGRID(IMAX,JMAX,nzg))
      ALLOCATE(VGRID(IMAX,JMAX,nzg))
      ALLOCATE(QGRID(IMAX,JMAX,nzg))
      ALLOCATE(CWMGRID(IMAX,JMAX,nzg))
      ALLOCATE(PRESGRID(IMAX,JMAX,nzg))
      ALLOCATE(PINTGRID(IMAX,JMAX,nzg+1))
      ALLOCATE(PGRID(IMAX,JMAX,3))

      write(6,*) 'iret from sigio_srohdc: ', iret
      write(6,*) 'head%levs: ', nzg
      write(6,*) 'head%jcap: ', head%jcap
      write(6,*) 'head%ntrac: ', head%ntrac
      write(6,*) 'From file: ',filename
      write(6,*) 'start COF2GRD',nxg,nyg

!     hs (surface topo)

      write(6,*) 'Process data%hs...'
      if ( use_gfs_nemsio ) then
        call nemsio_readrecv(gfile,'hgt','sfc',1,fis(:),iret=iret)
C$OMP PARALLEL DO PRIVATE(I,J,II,JJ)
! Subseting
        DO J=1,JMAX
           jj=j-1+jb
           DO I=1,IMAX
             ii=i-1+ib
             if(ii .lt. 1) ii=inemsio+ii
             if(ii .gt. inemsio) ii=ii-inemsio
             PGRID(I,JMAX-J+1,2)=fis(ii+(jj-1)*head%lonb)
           ENDDO
        ENDDO
C$OMP END PARALLEL DO

      else if(use_gfs_netcdf) then
        CALL get_ncf2d(TRIM(filename), nxg, nyg, 'hgtsfc', fld2d)
C$OMP PARALLEL DO PRIVATE(I,J,II,JJ)
! Subseting
        DO J=1,JMAX
           jj=j-1+jb
           DO I=1,IMAX
             ii=i-1+ib
             if(ii .lt. 1) ii=nxg+ii
             if(ii .gt. nxg) ii=ii-nxg
             PGRID(I,JMAX-J+1,2)=fld2d(ii,jj)
           ENDDO
        ENDDO
C$OMP END PARALLEL DO
      else
        CALL SPTRAN(0,head%JCAP,0,IMAX,JMAX,1,0,0,
     2     -IMAX,IMAX,
     1     0,0,0,0,1,data%hs,MYGRD(1,JMAX),MYGRD(1,1),1)

        DO J=1,JMAX 
        DO I=1,IMAX
            PGRID(I,J,2)=MYGRD(I,J)
        ENDDO
        ENDDO
      endif
      write(6,*) '  ... done.'

      write(6,*) 'data%hs onto PGRID'
      do J=1,JMAX,JMAX/40
         write(6,617) (PGRID(I,J,2),I=1,IMAX,IMAX/20)
      enddo

 617  format(30(f5.0,1x))

!     ps and midlayer P values


      write(6,*) 'Process data%ps...'
      if ( use_gfs_nemsio ) then
        call nemsio_readrecv(gfile,'pres','sfc',1,fis(:),iret=iret)
C$OMP PARALLEL DO PRIVATE(I,J,II,JJ)
! Subseting
        DO J=1,JMAX
          jj=j-1+jb
          DO I=1,IMAX
            ii=i-1+ib
            if(ii .lt. 1) ii=inemsio+ii
            if(ii .gt. inemsio) ii=ii-inemsio
            PGRID(I,JMAX-J+1,1)=fis(ii+(jj-1)*head%lonb)
          ENDDO
        ENDDO
C$OMP END PARALLEL DO
C$OMP PARALLEL DO PRIVATE(I,J)
        DO J=1,JMAX
        DO I=1,IMAX
           MYGRD(I,J)=PGRID(I,J,1)
        ENDDO
        ENDDO
C$OMP END PARALLEL DO

      else if(use_gfs_netcdf) then
        CALL get_ncf2d(TRIM(filename), nxg, nyg, 'pressfc', fld2d)
C$OMP PARALLEL DO PRIVATE(I,J,II,JJ)
! Subseting
        DO J=1,JMAX
           jj=j-1+jb
           DO I=1,IMAX
             ii=i-1+ib
             if(ii .lt. 1) ii=nxg+ii
             if(ii .gt. nxg) ii=ii-nxg
             PGRID(I,JMAX-J+1,1)=fld2d(ii,jj)
           ENDDO
        ENDDO
C$OMP END PARALLEL DO
C$OMP PARALLEL DO PRIVATE(I,J)
        DO J=1,JMAX
        DO I=1,IMAX
           MYGRD(I,J)=PGRID(I,J,1)
        ENDDO
        ENDDO
C$OMP END PARALLEL DO

      else
        CALL SPTRAN(0,head%JCAP,0,IMAX,JMAX,1,0,0,
     2     -IMAX,IMAX,
     1     0,0,0,0,1,data%ps,MYGRD(1,JMAX),MYGRD(1,1),1)
        write(6,*) '  ... done.'

C$OMP PARALLEL DO PRIVATE(I,J)
        DO J=1,JMAX 
        DO I=1,IMAX
           PGRID(I,J,1)=1000.*EXP(MYGRD(I,J))
           MYGRD(I,J)=PGRID(I,J,1)
        ENDDO
        ENDDO
C$OMP END PARALLEL DO
      endif

      write(6,*) 'data%ps onto PGRID'
      do J=1,JMAX,JMAX/40
         write(6,617) (PGRID(I,J,1)/100.,I=1,IMAX,IMAX/20)
      enddo

      write(6,*) 'Run sigio_modpr...'
      if (use_gfs_netcdf) then
        CALL sigio_modpr(IMAX*JMAX,IMAX*JMAX,nzg,
     +     nnvcoord, nidvc, nidsl, abk,
     +     iret,ps=MYGRD,pm=MYGRD2, pd=PD)
      else  
        CALL sigio_modpr(IMAX*JMAX,IMAX*JMAX,head%levs,
     +     head%nvcoord,head%idvc,head%idsl,          
     +     head%vcoord,iret,ps=MYGRD,pm=MYGRD2,
     +     pd=PD)
        write(6,*) '  ... done.'
      endif

!     temperatures

      write(6,*) 'before 3D read of T'
C$OMP PARALLEL DO PRIVATE(I,J)
        DO J=1,JMAX
          DO I=1,IMAX
             PINTGRID(I,J,1)=MYGRD(I,J)
          ENDDO
        ENDDO
C$OMP END PARALLEL DO

      write(6,*) 'Process data%t...'
      if ( use_gfs_nemsio ) then
        do l=1,head%levs
          call nemsio_readrecv(gfile,'tmp','mid layer',l,fis(:)
     &    ,iret=iret)
C$OMP PARALLEL DO PRIVATE(I,J,II,JJ)
! Subseting
          DO J=1,JMAX
            jj=j-1+jb
            DO I=1,IMAX
              ii=i-1+ib
              if(ii .lt. 1) ii=inemsio+ii
              if(ii .gt. inemsio) ii=ii-inemsio
              TGRID(I,JMAX-J+1,l)=fis(ii+(jj-1)*head%lonb)
            ENDDO
          ENDDO
C$OMP END PARALLEL DO
        enddo
      else if(use_gfs_netcdf) then
        CALL get_ncf3d(TRIM(filename), nxg, nyg, nzg, 'tmp', fld3d)
        do l=1,nzg
C$OMP PARALLEL DO PRIVATE(I,J,II,JJ)
! Subseting
           DO J=1,JMAX
             jj=j-1+jb
             DO I=1,IMAX
               ii=i-1+ib
               if(ii .lt. 1) ii=nxg+ii
               if(ii .gt. nxg) ii=ii-nxg
               TGRID(I,JMAX-J+1,l)=fld3d(ii,jj,nzg-l+1)
             ENDDO
           ENDDO
C$OMP END PARALLEL DO
         enddo
      else
        CALL SPTRAN(0,head%JCAP,0,IMAX,JMAX,head%levs,0,0,
     2        -IMAX,IMAX,
     1        0,0,0,0,1,data%t,TGRID(1,JMAX,1),TGRID(1,1,1),1)
        write(6,*) '  ... done.'
      endif

      DO L=1,nzg
C$OMP PARALLEL DO PRIVATE(I,J)
        DO J=1,JMAX
        DO I=1,IMAX
           PINTGRID(I,J,L+1)=PINTGRID(I,J,L)-PD(I,J,L)
           IF (L .eq. nzg .and. PINTGRID(I,J,L+1) .lt. 5.) then
              PINTGRID(I,J,L+1)=5.
           ENDIF
        ENDDO
        ENDDO
C$OMP END PARALLEL DO
      ENDDO
      I=1
      J=1
      DO L=1,nzg
         write(6,*) '1,1,PINT,PD:: ', 
     &       I,J,L,PINTGRID(I,J,L+1),PD(I,J,L)
      ENDDO

      write(6,*) 'after 3D read of T'
      call summary()

      print *,'ok after terrain coeffs'
C     
C     READ SFC PRESSURE COEFFICIENTS
C     

C$OMP PARALLEL DO PRIVATE(I,J,L)
      DO L = 1, nzg
      DO J = 1, JMAX
      DO I = 1, IMAX
         PRESGRID(I,J,L) =  MYGRD2(I,J,L)
      ENDDO
      ENDDO
      ENDDO
C$END PARALLEL DO
      I=1
      J=1
      DO L = 1, nzg
         write(6,*) 'L, PRESGRID(1,1,L): ', L, PRESGRID(1,1,L)
      ENDDO

      DEALLOCATE(GRD)
C     
C     READ DIVERGENCE AND VORTICITY COEFFICIENTS
C     
      write(6,*) 'Process data%d and data%z (velocity)...'
      if ( use_gfs_nemsio ) then
        DO L = 1, head%levs
          call nemsio_readrecv(gfile,'ugrd','mid layer',l,fis(:)
     &,   iret=iret)
          call nemsio_readrecv(gfile,'vgrd','mid layer',l,fit(:)
     &,   iret=iret)
C$OMP PARALLEL DO PRIVATE(I,J,II,JJ)
          DO J = 1, JMAX
            jj=j-1+jb
            DO I = 1, IMAX
               ii=i-1+ib
               if(ii .lt. 1) ii=inemsio+ii
               if(ii .gt. inemsio) ii=ii-inemsio
               UGRID(I,JMAX-J+1,L)=fis(ii+(jj-1)*head%lonb)
               VGRID(I,JMAX-J+1,L)=fit(ii+(jj-1)*head%lonb)
            ENDDO
          ENDDO
C$END PARALLEL DO
        ENDDO
      else if(use_gfs_netcdf) then
        CALL get_ncf3d(TRIM(filename), nxg, nyg, nzg, 'ugrd', fld3d)
        DO L = 1, nzg
C$OMP PARALLEL DO PRIVATE(I,J,II,JJ)
           DO J = 1, JMAX
             jj=j-1+jb
             DO I = 1, IMAX
               ii=i-1+ib
               if(ii .lt. 1) ii=nxg+ii
               if(ii .gt. nxg) ii=ii-nxg
               UGRID(I,JMAX-J+1,L)=fld3d(ii,jj,nzg-l+1)
             ENDDO
           ENDDO
C$END PARALLEL DO
        ENDDO
        CALL get_ncf3d(TRIM(filename), nxg, nyg, nzg, 'vgrd', fld3d)
        DO L = 1, nzg
C$OMP PARALLEL DO PRIVATE(I,J,II,JJ)
           DO J = 1, JMAX
             jj=j-1+jb
             DO I = 1, IMAX
               ii=i-1+ib
               if(ii .lt. 1) ii=nxg+ii
               if(ii .gt. nxg) ii=ii-nxg
               VGRID(I,JMAX-J+1,L)=fld3d(ii,jj,nzg-l+1)
             ENDDO
           ENDDO
C$END PARALLEL DO
        ENDDO
      else
        ALLOCATE(DWORK(NC))
        ALLOCATE(ZWORK(NC))
        CALL SPTRANV(0,head%jcap,0,IMAX,JMAX,head%levs,0,0,-IMAX,IMAX,
     1  0,0,0,0,1,data%d,data%z,UGRID(1,JMAX,1),UGRID(1,1,1),
     2  VGRID(1,JMAX,1),VGRID(1,1,1),1)
        DEALLOCATE(DWORK)
        DEALLOCATE(ZWORK)
      endif

      write(6,*) '  ... done.'
      DO L = 1, nzg
        write(6,*) 'U(50,50),V(50,50): ', 
     &      L, UGRID(50,50,L),VGRID(50,50,L)
      ENDDO
C       DO L = 1, nzg
C         CALL SPTRANV(0,head%jcap,0,IMAX,JMAX,1,0,0,-IMAX,IMAX,
C      1   0,0,0,0,1,data%d(:,L),data%z(:,L),UGRID(1,JMAX,L),UGRID(1,1,L),
C      2   VGRID(1,JMAX,L),VGRID(1,1,L),1)
C         write(6,*) 'U(50,50),V(50,50): ', 
C      &      L, UGRID(50,50,L),VGRID(50,50,L)
C       ENDDO
      print *,'ok after div/vort coeffs'

C     
C     READ SPECIFIC HUMIDITY COEFFICIENTS
C     
      write(6,*) 'Process data%q(:,:,1) (water vapor)...'
      if ( use_gfs_nemsio ) then
        do l=1,head%levs
          call nemsio_readrecv(gfile,'spfh','mid layer',l,fis(:)
     &,   iret=iret)
C$OMP PARALLEL DO PRIVATE(I,J,II,JJ)
          DO J = 1, JMAX
            jj=j-1+jb
            DO I = 1, IMAX
               ii=i-1+ib
               if(ii .lt. 1) ii=inemsio+ii
               if(ii .gt. inemsio) ii=ii-inemsio
               QGRID(I,JMAX-J+1,L)=fis(ii+(jj-1)*head%lonb)
            ENDDO
          ENDDO
C$END PARALLEL DO
        enddo

      else if(use_gfs_netcdf) then
        CALL get_ncf3d(TRIM(filename), nxg, nyg, nzg, 'spfh', fld3d)
        do l=1,nzg
C$OMP PARALLEL DO PRIVATE(I,J,II,JJ)
          DO J = 1, JMAX
            jj=j-1+jb
            DO I = 1, IMAX
               ii=i-1+ib
               if(ii .lt. 1) ii=nxg+ii
               if(ii .gt. nxg) ii=ii-nxg
               QGRID(I,JMAX-J+1,L)=fld3d(ii,jj,nzg-l+1)
            ENDDO
          ENDDO
C$END PARALLEL DO
        enddo
      else
        CALL SPTRAN(0,head%JCAP,0,IMAX,JMAX,head%levs,0,0,-IMAX,IMAX,
     1     0,0,0,0,1,data%q(:,:,1),QGRID(1,JMAX,1),QGRID(1,1,1),1)
      endif
      write(6,*) '  ... done.'
C       DO L = 1, nzg
C          CALL SPTRAN(0,head%JCAP,0,IMAX,JMAX,1,0,0,-IMAX,IMAX,
C      1        0,0,0,0,1,data%q(:,L,1),QGRID(1,JMAX,L),QGRID(1,1,L),1)
C       ENDDO
      print *,'ok after q coeffs'

C$OMP PARALLEL DO PRIVATE(I,J,K)
      DO K = 1, nzg
         DO J = 1, JMAX
            DO I = 1, IMAX
               QGRID(I,J,K) = AMAX1(QGRID(I,J,K),1.0E-12)
            ENDDO
         ENDDO
      ENDDO
C$OMP END PARALLEL DO

! Convert temperature to virtual temperature (nemsio output T rather Tv
C$OMP PARALLEL DO
C$OMP& PRIVATE(I,J,K)
       DO K=1,nzg
        DO J=1,JMAX
         DO I=1,IMAX
          TGRID(I,J,K)=TGRID(I,J,K)*(1.0+0.608*QGRID(I,J,K))
         ENDDO
        ENDDO
       ENDDO
C$OMP END PARALLEL DO
      do L=1,nzg
         write(6,*) 'q element 1...L,Q(10,180,L): ', L, QGRID(10,180,L)
      enddo
C     

      CWMGRID=-9999.
      write(6,*) 'Process data%q(:,:,3) (total condensate)...'
      if ( use_gfs_nemsio ) then
        do l=1,head%levs
!          call nemsio_readrecv(gfile,'clwmr','mid layer',l,fis(:)
!     &,   iret=iret)
          call nemsio_readrecv(gfile,'cld_amt','mid layer',l,fis(:)
     &,   iret=iret)
C$OMP PARALLEL DO PRIVATE(I,J,II,JJ)
          DO J = 1, JMAX
            jj=j-1+jb
            DO I = 1, IMAX
               ii=i-1+ib
               if(ii .lt. 1) ii=inemsio+ii
               if(ii .gt. inemsio) ii=ii-inemsio
               CWMGRID(I,JMAX-J+1,L)=fis(ii+(jj-1)*head%lonb)
            ENDDO
          ENDDO
C$END PARALLEL DO
        enddo

      else if(use_gfs_netcdf) then
        deallocate(fld2d, fld3d)
        allocate(fld3d_cl(nxg,nyg,nzg), fld3d_gr(nxg,nyg,nzg), 
     &  fld3d_ic(nxg,nyg,nzg),fld3d_rw(nxg,nyg,nzg),
     &  fld3d_sn(nxg,nyg,nzg))
        CALL get_ncf3d(TRIM(filename), nxg, nyg, nzg,'clwmr', fld3d_cl)
        CALL get_ncf3d(TRIM(filename), nxg, nyg, nzg, 'grle', fld3d_gr)
        CALL get_ncf3d(TRIM(filename), nxg, nyg, nzg, 'icmr', fld3d_ic)
        CALL get_ncf3d(TRIM(filename), nxg, nyg, nzg, 'rwmr', fld3d_rw)
        CALL get_ncf3d(TRIM(filename), nxg, nyg, nzg, 'snmr', fld3d_sn)
        do l=1,nzg
C$OMP PARALLEL DO PRIVATE(I,J,II,JJ)
          DO J = 1, JMAX
            jj=j-1+jb
            DO I = 1, IMAX
               ii=i-1+ib
               if(ii .lt. 1) ii=nxg+ii
               if(ii .gt. nxg) ii=ii-nxg
               CWMGRID(I,JMAX-J+1,L)=fld3d_cl(ii,jj,nzg-l+1)+
     &          fld3d_gr(ii,jj,nzg-l+1)+fld3d_ic(ii,jj,nzg-l+1)+
     &          fld3d_rw(ii,jj,nzg-l+1)+fld3d_sn(ii,jj,nzg-l+1)
c               CWMGRID(I,JMAX-J+1,L)=fld3d_cl(ii,jj,nzg-l+1)
            ENDDO
          ENDDO
C$END PARALLEL DO
        enddo
        deallocate(fld3d_cl,fld3d_gr,fld3d_ic,fld3d_rw,fld3d_sn)
      else
        CALL SPTRAN(0,head%JCAP,0,IMAX,JMAX,head%levs,0,0,-IMAX,IMAX,
     1     0,0,0,0,1,data%q(:,:,3),CWMGRID(1,JMAX,1),CWMGRID(1,1,1),1)
      endif
      write(6,*) '  ... done.'
C       DO L = 1, nzg
C          CALL SPTRAN(0,head%JCAP,0,IMAX,JMAX,1,0,0,-IMAX,IMAX,
C      1     0,0,0,0,1,data%q(:,L,3),CWMGRID(1,JMAX,L),CWMGRID(1,1,L),1)
C       ENDDO
      write(6,*) 'maxval CWMGRID: ', maxval(CWMGRID)
      do L=1,nzg
         write(6,*) 'q(3) ..L,CWMGRID(10,180,L): ',L,CWMGRID(10,180,L)
      enddo

      LUN1HOLD=LUN1+100
      WRITE(HOLDFIL,1000)LUN1
 1000 FORMAT('holdsig',I3.3)
c$$$      OPEN(UNIT=LUN1HOLD,FILE=HOLDFIL,FORM='UNFORMATTED',IOSTAT=IER)
c$$$C     
c$$$      write(6,*) 'TGRID(1,1,1): ', TGRID(1,1,1)
c$$$      WRITE(LUN1HOLD)TGRID
c$$$      write(6,*) 'UGRID(1,1,1): ', UGRID(1,1,1)
c$$$      WRITE(LUN1HOLD)UGRID
c$$$      write(6,*) 'VGRID(1,1,1): ', VGRID(1,1,1)
c$$$      WRITE(LUN1HOLD)VGRID
c$$$      write(6,*) 'QGRID(1,1,1): ', QGRID(1,1,1)
c$$$      WRITE(LUN1HOLD)QGRID
c$$$      write(6,*) 'CWMGRID(1,1,1): ', CWMGRID(1,1,1)
c$$$      WRITE(LUN1HOLD)CWMGRID
c$$$      write(6,*) 'PRESGRID(1,1,1): ', PRESGRID(1,1,1)
c$$$      WRITE(LUN1HOLD)PRESGRID
c$$$      write(6,*) 'PINTGRID(1,1,1): ', PINTGRID(1,1,1)
c$$$      WRITE(LUN1HOLD)PINTGRID
c$$$      write(6,*) 'PGRID(1,1,1): ', PGRID(1,1,1)
c$$$      WRITE(LUN1HOLD)PGRID
c$$$      CLOSE(LUN1HOLD)
c$$$      DEALLOCATE(TGRID)
c$$$      DEALLOCATE(UGRID)
c$$$      DEALLOCATE(VGRID)
c$$$      DEALLOCATE(QGRID)
c$$$      DEALLOCATE(CWMGRID)
c$$$      DEALLOCATE(PRESGRID)
c$$$      DEALLOCATE(PINTGRID)
c$$$      DEALLOCATE(PGRID)
c$$$      DEALLOCATE(PD)

      KMAX=nzg
      if ( use_gfs_nemsio ) then
        deallocate(fis)
        deallocate(fit)
        deallocate(vcoord)
        deallocate(lat,lon)
      endif
      if ( use_gfs_netcdf ) then
         deallocate(ak, bk, abk, pfull, phalf)
      endif

      write(6,*) 'end COF2GRD'
      if ( use_gfs_nemsio ) call nemsio_close(gfile,iret=iret)
      call summary()
C     


C Deallocate data:
      if ( .not. use_gfs_netcdf ) then
        CALL sigio_axdata(data,iret)
      endif


      RETURN
      END
