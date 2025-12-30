PROGRAM wind2hycom
!---------------------------------------------------------------------------
! 
! WIND2HYCOM 
!    read in a set of wind properties on a given domain, and 
!    do horizontal interpolation and rotation of them on an OFS domain.
!
!-------
! input: 
!-------
!      from 'wind2hycom.in' file
!      		xyfname = grid (lon,lat,etc) and time info file
!      		binfname = (time,im,jm,3) in binary
!
!      lat,lon, & workflxs, 3-D (wind-speed,taux,tauy)
!
!--------
! output:
!--------
!      output:   hycomflxs, 3-D (wind-speed,taux,tauy)
!
!
!---------
! history:
!--------- &
!  10/16/2008 Extend the number of input fields from 3 to 4  
!  10/30/2007 Edited by Carlos
!  10/25/2007 by Hyun-Sook Kim
!    
!---------------------------------------------------------------------------

use mod_za
use mod_geom
implicit none
real, dimension(:), allocatable :: exhycom,eyhycom,exmrf,eymrf
real, dimension(:,:), allocatable :: qxhycom2d,qyhycom2d,anhycom2d
real, dimension(:,:), allocatable :: exhycom2d,eyhycom2d

integer :: i,nxhycom,nyhycom,dbgn,iii,jjj,idum,iu
real    :: cs,ss

character (len=128) :: binfname,file4hycom,grid_infofile

integer :: k, ntimes, inthour, startdate0,nxmrf, nymrf
logical :: flip
real, dimension(:,:,:), allocatable :: var3d
real, dimension(:,:,:), allocatable :: hycomflxs
real, dimension(:,:), allocatable :: hycomflxu,hycomflxv


!------------------------------------------------------------------------------
! Main part
!-----------
       
        call xcspmd
        nxhycom=idm
        nyhycom=jdm
        call zaiost
      	allocate ( &
 &     exhycom2d(nxhycom+1,nyhycom+1),eyhycom2d(nxhycom+1,nyhycom+1) &
 &    ,qxhycom2d(nxhycom+1,nyhycom+1),qyhycom2d(nxhycom+1,nyhycom+1) &
 &    ,anhycom2d(nxhycom,nyhycom)                                    &
 &         )
        allocate(hycomflxs(nxhycom,nyhycom,4))
        allocate(hycomflxu(nxhycom,nyhycom),hycomflxv(nxhycom,nyhycom))

      	call hycom_na(anhycom2d,qxhycom2d,qyhycom2d,dbgn)

  
! 
! READ generic grid (lon,lat)


  read(5,'(a)')binfname  !output binamy file name
  read(5,'(a)')file4hycom
  read(5,*)grid_infofile
  
  idum=101
  iu=12
  
  open(iu,file=grid_infofile,form='formatted',status='old')
  call read_grid_parm(ntimes, inthour, startdate0,nxmrf, nymrf,iu)
  
  allocate ( exmrf(nxmrf),eymrf(nymrf))
  call read_grid_lonlat(exmrf,eymrf,nxmrf,nymrf,flip,iu)
  
  close(iu)
  allocate ( var3d(nxmrf,nymrf,4) )

!sanity check ex monotonic increasing and ey monotonic decreasing  
  
! else stop

!
! Pre-compute interpolation indices and weights for multiple interpolations

      	call horiz_intp(exmrf,eymrf,qxhycom2d,qyhycom2d,exhycom2d,eyhycom2d,idum)

!
! READ var3D=(fm, |Wind|, Taux, Tauy)

  open(15,file=trim(binfname),status='new',form='unformatted')
  open(10,file=trim(file4hycom),status='old',form='unformatted')

! 
! iteration for number of total records  

  do k=1,ntimes
        read(10) var3d

!flip ud along 2nd index

        if(flip) call ud(var3d,nxmrf,nymrf,4)

!
        write(6,*)k,' time frame', ntimes,'/total'
      	do i=1,4
           call horiz_intp_(var3D(:,:,i),hycomflxs(:,:,i),    &
&                           exhycom2d(:,:),eyhycom2d(:,:))
      	enddo

!
! Vector variable
!- zonal component
      	i=3
      	hycomflxu=hycomflxs(:,:,i)

!- meridional component
      	i=4
!
! Rotating     
      	do jjj=1,nyhycom
           do iii=1,nxhycom
            	cs=cos(anhycom2d(iii,jjj))
            	ss=sin(anhycom2d(iii,jjj))
            	hycomflxv(iii,jjj)= hycomflxs(iii,jjj,i)*cs-    &
&                               hycomflxu(iii,jjj  )*ss  
           	hycomflxu(iii,jjj)= hycomflxs(iii,jjj,i)*ss+    &
&                               hycomflxu(iii,jjj  )*cs      

!cj               speed=max(minspd,sqrt(hycomflxu(iii,jjj)**2+hycomflxv(iii,jjj)**2))
!cj               cdval=min(cd(speed),cdmax)*speed
!here we change sign to compensate for a change in sign in the averaging loop
!cj               hycomflxu(iii,jjj)=-  cdval*hycomflxu(iii,jjj)
!cj               hycomflxv(iii,jjj)=-  cdval*hycomflxv(iii,jjj)
	   enddo
        enddo

!
! replace old values to new values after rotation
     	hycomflxs(:,:,3)=hycomflxu
     	hycomflxs(:,:,4)=hycomflxv
  
!
! final stage:
!- write out in binary

        write(15) hycomflxs
  end do

  close(10)
  close(15)
  
 end program wind2hycom

!---------------------------------------------------------------------------
 subroutine read_grid_parm(ntimes, inthour, startdate0,nxmrf, nymrf,iu)
   implicit none
   integer :: ntimes, inthour, startdate0,nxmrf, nymrf,iu
   read(iu,*)ntimes,inthour,startdate0
   read(iu,*)nxmrf, nymrf
   return
 end subroutine read_grid_parm
 
!---------------------------------------------------------------------------
 subroutine read_grid_lonlat(lon,lat,im,jm,flip,iu)
   implicit none
   integer          :: iu,im,jm,i
   real,dimension(im)::lon
   real,dimension(jm)::lat,tmp
   real,parameter:: radian=57.2957795
   real              :: dum
   integer           ::hj
   logical           ::flip

   do i=1,im
    read(iu,*)lon(i)
    lon(i)=lon(i)/radian
   enddo

   do i=1,jm
    read(iu,*)lat(i)
    lat(i)=lat(i)/radian
   enddo

   if(lat(2)>lat(1)) then
! flip
   flip=.true.
   hj=jm/2
   do i=1,hj
    	dum=lat(i)
    	lat(i)=lat(jm-i+1)
    	lat(jm-i+1)=dum
   enddo
   else
   	flip=.false.
   endif
 
  return
 end subroutine read_grid_lonlat
 
!---------------------------------------------------------------------------
 subroutine ud(a,n1,n,n2)
    integer idim
    real, dimension(n1,n,n2) :: a
    integer n,nh,i,n1,n2,j,k
    real ::slab

    nh=n/2

    do i=1,nh
       do k=1,n2
          do j=1,n1
             slab=a(j,i,k)
             a(j,i,k)=a(j,n-i+1,k)
             a(j,n-i+1,k)=slab
          enddo
       enddo
    enddo
    
    return
 end subroutine ud
       
