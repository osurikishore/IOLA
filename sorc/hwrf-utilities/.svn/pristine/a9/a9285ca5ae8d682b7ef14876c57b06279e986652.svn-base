      subroutine read_nemsio(filename,head,data,iret)

      use nemsio_module
      USE SIGIO_MODULE
      IMPLICIT NONE
      TYPE(sigio_head):: HEAD
      TYPE(sigio_data):: DATA
!
      type(nemsio_gfile) :: gfile
      integer im,jm,km,jcap,jrec,ntrac,nframe,nrec,fieldsize,vlev,iret
      integer iz,jz,kz,ll
      character(8) vname
      character(16) vlevtyp
      character filename*80
      real(4),allocatable :: fis(:),fit(:),vcoord(:,:,:)
      real(4),allocatable :: clwm(:,:),ugrd(:,:), vgrd(:,:)
      real(4),allocatable :: tmp(:,:),spfh(:,:)
      real,allocatable:: slat(:),clat(:),wlat(:)
!
!--- Initialize
      iret=0
      print*, "zhang xxxx iret===", iret, trim(filename)
      call nemsio_init(iret=iret)
      print*, "zhang yyyy iret====" ,iret 
!
!--- Opean a NEMSIO file
      call nemsio_open(gfile,trim(filename),'READ',iret=iret)
!
      print*, "zhang aaaa", trim(filename)
      print*, "zhang bbbb"
!--- Get dimension
      call nemsio_getfilehead(gfile,iret=iret,jcap=head%jcap,
     & dimx=head%lonb,dimy=head%latb,idsl=head%idsl,idvc=head%idvc,
     & dimz=head%levs,ntrac=head%ntrac,nframe=nframe)
      print*, "zhang cccc"
      im=head%lonb
      jm=head%latb
      km=head%levs
      jcap=head%jcap
      head%latf=jm
      head%lonf=im
      head%nvcoord=2
      call sigio_alhead(head,iret)
      call sigio_aldata(head,data,iret)
      print*, "zhang dddd", im,jm,head%ntrac,nframe,head%idsl,head%idvc
!
!---Allocate array
      fieldsize=(im+2*nframe)*(jm+2*nframe)
      allocate(fis(fieldsize))
      allocate(fit(fieldsize))
      allocate(vcoord(km+1,3,2))
      allocate(slat(jm),clat(jm),wlat(jm))
      call splat(4,jm,slat,wlat)
      clat=acos(sqrt(1-slat**2))*180./3.1415926
      slat=asin(slat)*180./3.1415926
      write(29,*) "--------slat-----"
      write(29,*) slat
      write(29,*) "--------clat-----"
      write(29,*) clat
      call nemsio_getfilehead(gfile,iret=iret,vcoord=vcoord)
      print*, "zhang fieldsize===", fieldsize
      do jz=1,2
      do kz=1,km+1
      print*, "zhang vcoord===",jz,kz 
      head%vcoord(kz,jz)=vcoord(kz,jz,1)
      print*, kz,head%vcoord(kz,jz)
      enddo
      enddo
!
!--- You could get the name, levtyp, and level of j-th record data field by calling:
!zz      jrec=1
!zz      call nemsio_getrechead(gfile,jrec,vname,vlevtyp,vlev,iret=iret)
!
!--- Get one data field out by data field record number
!zz      call nemsio_readrec(gfile,jrec,fis(:),iret=iret)
!
!--- Get one data field out by giving field name, levtyp, and level
!zz      call nemsio_readrecv(gfile,'hgt','sfc',1,fis(:),iret=iret)
! get sfchgt and ps on Gaussian grids, and convert them to spectrums
      call nemsio_readrecv(gfile,'hgt','sfc',1,fis(:),iret=iret)
      call sptez(0,jcap,4,im,jm,data%hs,fis,-1)
      call nemsio_readrecv(gfile,'pres','sfc',1,fis(:),iret=iret)
      fis=log(fis*1.e-3)
      call sptez(0,jcap,4,im,jm,data%ps,fis,-1)
      deallocate(fis)
      deallocate(fit)

      allocate(clwm(fieldsize,km))
      allocate(ugrd(fieldsize,km))
      allocate(vgrd(fieldsize,km))
      allocate(tmp(fieldsize,km))
      allocate(spfh(fieldsize,km))
      do ll=1,km
      print*, "ll=====",ll
      call nemsio_readrecv(gfile,'clwm','mid layer',ll,clwm(:,ll)
     &,iret=iret)
      call nemsio_readrecv(gfile,'ugrd','mid layer',ll,ugrd(:,ll)
     &,iret=iret)
      call nemsio_readrecv(gfile,'vgrd','mid layer',ll,vgrd(:,ll)
     &,iret=iret)
      call nemsio_readrecv(gfile,'tmp','mid layer',ll,tmp(:,ll)
     &,iret=iret)
      call nemsio_readrecv(gfile,'spfh','mid layer',ll,spfh(:,ll)
     &,iret=iret)
      enddo
!
!--- Close the NEMSIO file
      call nemsio_close(gfile,iret=iret)
!
!--- Finalize
      call nemsio_finalize()

! get variables on model vertical levels on Gaussian grids, and convert them to spectrums
C$OMP PARALLEL DO PRIVATE(ll,iret)
!
      do ll=1,km
      print*, "ll=====",ll
      call sptez(0,jcap,4,im,jm,data%q(1,ll,3),clwm(:,ll),-1)
      call sptezv(0,jcap,4,im,jm,data%d(1,ll),data%z(1,ll)
     &,ugrd(:,ll),vgrd(:,ll),-1)
      tmp(:,ll)=tmp(:,ll)*(1+(461.50/287.05-1)*spfh(:,ll))
      call sptez(0,jcap,4,im,jm,data%t(1,ll),tmp(:,ll),-1)
      call sptez(0,jcap,4,im,jm,data%q(1,ll,1),spfh(:,ll),-1)
      enddo
C$OMP END PARALLEL DO
      print*, "end of read nemsio"
      deallocate(clwm)
      deallocate(ugrd)
      deallocate(vgrd)
      deallocate(tmp)
      deallocate(spfh)
      deallocate(vcoord)
      return
!
      end
