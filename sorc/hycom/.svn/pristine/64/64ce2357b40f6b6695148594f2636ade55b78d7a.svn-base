c####################################################################
c   WIND.F                                                          #
c           Originally written by a hurricane group of GSO/URI      #
c           for a feature model for a hurrican wind field.          #
c                                                                   #
c--------                                                           #
c purpose                                                           #
c--------                                                           #
c        1. computes wind (x,y) and stress (x,y) at a given time    #
c        2. modified so that produes a set of wind-speed &          #
c                    wind-stress components for HYCOM               #
c                                                                   #
c--------                                                           #
c history                                                           #
c--------                                                           #
c      10/25/2007 Modified by Hyun-Sook Kim, for easy control input #
c                 variables - filename, startdate, hour.            #
c      10/15/2008 signify in/out hurricane in the wind speed        #
c                 cd correction tbd                                 #
c####################################################################

      integer*4 startdate
      character*15 filename
      REAL LATMIN,LATMAX,LONGMIN,LONGMAX
      COMMON/sphere/LATMIN,LATMAX,LONGMIN,LONGMAX

      integer IM,JM
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: VAR3D
      REAL, DIMENSION(:,:), ALLOCATABLE :: wmag0,taux0,tauy0,fm
      REAL, DIMENSION(:), ALLOCATABLE :: xs,ys
      INTEGER I, J, startdate0, nrcds, k, l, idm,jdm
      REAL inthour
      character*80 binfname
      character*80 xyfname
      
      
chsk following are from control.in:
c
c      filename='ISABEL.message'
c      startdate=03091512
c      hour=0.
      
      read(5,*) filename
      read(5,*) idm
      read(5,*) jdm
      read(5,*) startdate0
      read(5,*) inthour
      read(5,*) nrcds
      read(5,*) xyfname
      read(5,*) binfname
      read(5,*) LATMIN
      read(5,*) LATMAX
      read(5,*) LONGMIN
      read(5,*) LONGMAX     	
      
c      inowind=73
      im=idm
      jm=jdm
      allocate ( wmag0(IM,JM) )
      allocate ( taux0(IM,JM) )
      allocate ( tauy0(IM,JM) )
      allocate ( fm(IM,JM) )
      allocate ( xs(IM) )
      allocate ( ys(JM) )
     
chsk following are the GSO HWRF-POM domain info:
c         extended domain in order to cover our HYCOM western hurrican domain
c      LATMIN=5.
c      LATMAX=51.
c      LONGMIN=-102.5
c      LONGMAX=-43.5

      do I=1,IM
         xs(I)=LONGMIN+(I-1)*(LONGMAX-LONGMIN)/(IM-1)
      end do
   
      do J=1,JM
         ys(J)=LATMIN+(J-1)*(LATMAX-LATMIN)/(JM-1)
      end do
	
chsk computes feature model wind properties:
c change this to be passed in
c
      write(*,*) 'nrecords=',nrcds

chsk  write out 4-D variable
      open(48,file=binfname,status='new',form='unformatted')

      allocate ( var3d(IM,JM,4) )
      do i=1,nrcds
         startdate=startdate0         
         hour=inthour*(i-1)
	 time=hour/24.  
c     	 if((time*24.0).lt.float(inowindh)) then
         call wind(filename,startdate,time,im,jm,taux0,tauy0,wmag0,fm)
c        end if
      
chsk  deposit estimated value onto a 4-D new workspace
c
         var3d(:,:,1)=fm
         var3d(:,:,2)=wmag0
         var3d(:,:,3)=taux0
         var3d(:,:,4)=tauy0
         write(48) var3d
	 
      end do
      close(48)    
      
chsk  write output 1
      open(38,file=xyfname,status='replace',form='formatted')
      write(38,'(1x,i4,1x,f4.2,1x,i012)') nrcds,inthour,startdate0
      write(38,'(1x,i4,1x,i4)') im,jm
      write(38,'(1x,f12.4)') xs
      write(38,'(1x,f12.4)') ys
      close(38)     

      
      stop
      end

c==============================================================================
      SUBROUTINE WIND(filename,startdate,time,im,jm,taux,tauy,
     $     wndspd,fmeld)
      INTEGER PLN
c      PARAMETER (IM=254,JM=225)
      INTEGER IM, JM
      PARAMETER(PLN=100)
      integer,parameter::reset_flag=0
c
      integer hour, lat, long, mx, rmw
      integer*4 startdate,date
      integer day,month,year
      integer garb(5),Rd1(4),Rd2(4)
      character*19 name
      character*15 filename
      character*1 ns,ew
c
      DIMENSION X(PLN),Y(PLN),TM(PLN),PRES(PLN),PRES0(PLN),
     *      RMAXa(PLN),WSPMAX(PLN),
     *      R18v(PLN,4),R26v(PLN,4),Rref18v(5),Rref26v(5),alphv(5)
      DIMENSION RAD(14),WS(14),RADM(14),WSM(14),ANGL(14)
      REAL CMP,T1,T2,F0,F1,L0,L1,REARTH,R,A7,B,E,DELP,x0,y0
      REAL DELTAX,DELTAX1,DELTAY,DELTAY1,DXDY,julday
      REAL WUSURF(IM,JM),WVSURF(IM,JM),TAUX(IM,JM),TAUY(IM,JM)
      REAL COR(IM,JM)
      REAL UWIND(IM,JM),VWIND(IM,JM)
      REAL LATMIN,LATMAX,LONGMIN,LONGMAX

      COMMON/sphere/LATMIN,LATMAX,LONGMIN,LONGMAX

c<--- HSK 10/25/2007
c<---     10/15/2008 add lower case variables (see arguments in the subroutine wind)
      real fmask,bc
      REAL WNDSPD(IM,JM),fmeld(im,jm)
c---> HSK       

      DATA RM,R0/60.E3,480.E3/,RHO_0/1024.E0/
      DATA RAD/0.,.4,.7,.8,.95,1.,1.35,2.7,4.05,5.4,6.75
     * ,8.1,10.8,13.5/
      DATA WS/0.,.1,.5,.8,.95,1.,.97,.72,.54,.44,.4,.36
     * ,.27,.23/
      DATA ANGL/0.,2.,4.,6.,7.,7.,14.,23.,24.,22.,
     * 21.,21.,21.,21./
      
      print*,'In subroutine WIND ...'
c
      REARTH=6371.E3
      WIND_SCALE=0.8
      ROA=1.28
      RMAX=50.E3
      PI=3.1415927
      E=exp(1.)
c
c----------------------- Reading message file --------------------
c 2/18/16 - new format:
c
c          1         2         3         4         5         6         7         8
c 12345678901234567890123456789012345678901234567890123456789012345678901234567890
c NHC  01L ALEX      20160114 0000 279N 0304W 045 067 0988 1012 0927 26 056 0111 0
c
c          9         0         1         2         3         4         5         6
c 12345678901234567890123456789012345678901234567890123456789012345678901234567890
c 222 0222 -999 M -999 0056 -999 -999 72 560N 0340W -999 -999 -999 -999
c
c 1-43 name,date,hour,lat,ns,long,ew
c 44-66 garb
c 67-73 mx,rmw
c 74-93 Rd1
c 94-115 Rd2

  17  format(A19,I8,1x,I4,1x,I3,a1,1x,I4,a1,
     * 1x,I3,1x,I3,1x,I4,1x,I4,1x,I4,
     * 1x,I2,1x,I3,
     * 1x,I4,1x,I4,1x,I4,1x,I4,
     * 3x,I4,1x,I4,1x,I4,1x,I4)
      print*,'reading file ',filename
      open(15,file=filename,status='old')
      end=0.
      I=0
      do while(end.eq.0)
        read(15,17) name,date,hour,lat,ns,long,ew,garb,mx,rmw,Rd1,Rd2
        write(6,17) name,date,hour,lat,ns,long,ew,garb,mx,rmw,Rd1,Rd2
        if(date.eq.0) goto 20
        I=I+1
c
        date=date*100+hour/100
        call date2day(year,julday,date)
c
        TM(i)=julday
        if(ns.eq.'S') then
           lat=-lat
        endif
        if(ew.eq.'W') then
           long=-long
        endif
        X(i)=-long/10.
        Y(i)=lat/10.
        PRES(i)=float(garb(3))
        PRES0(i)=float(garb(4))
        WSPMAX(i)=float(mx)
        RMAXa(i)=float(rmw)
c----------------------------------------------------------------
c to use info in any element in the set
        if(reset_flag.eq.1) then
c- for R18v	
	   if (Rd1(n).eq.-999) then
	      do k=n+1,4
	         if (Rd1(k).ne.-999) then
		     Rd1(n)=Rd1(k)
	         endif
   	      enddo
	   endif
	   if (Rd1(4).eq.-999.and.Rd1(1).ne.-999) Rd1(4)=Rd1(1)
c- for R26v
	   if (Rd2(n).eq.-999) then
	      do k=n+1,4
	         if (Rd2(k).ne.-999) then
		     Rd2(n)=Rd2(k)
	         endif
   	      enddo
	   endif
	   if (Rd2(4).eq.-999.and.Rd2(1).ne.-999) Rd2(4)=Rd2(1)
       else
         print*,'Do not reset flags'
         endif
c----------------------------------------------------------------
	do n=1,4
          n1=n+(1-mod(n,2))*sign(2,3-n)
          R18v(i,n)=Rd1(n1)
          R26v(i,n)=Rd2(n1)
          if(wspmax(i).le.26.or.R26v(i,n).le.RMAXa(i)) R26v(i,n)=-999
          if(wspmax(i).le.18.or.R18v(i,n).le.RMAXa(i)) R18v(i,n)=-999
          if(R26v(i,n).gt.R18v(i,n)) R26v(i,n)=-999
        end do
      end do
  20  end=1.
      ipmax=I
      print*,'Number of hurricane path datapoints read: ',ipmax
      print*,'tm=',(tm(i),i=1,ipmax)
      print*,'year=',year
      print*,'startdate=',startdate
      close(15)
c
c--------------------- Calculating starting day -----------------
c
      call date2day(year,julday,startdate)
      do i=1,ipmax
        TM(i)=TM(i)-julday
      end do
C++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  INTERPOLATION OF HURRICANE PATH TO DETERMINE THE CURRENT POSITION

      CMP=TM(ipmax)
      if(cmp.lt.time.or.time.lt.tm(1)) then
        print*,'NO HURRICANE PATH DATA FOR THIS TIME'
        print*,'   time=',time
        print*,'   tm(1)=',tm(1)
        print*,'   tm(ipmax)=',tm(ipmax)
        print*,'   day=',julday
        RETURN
      end if

      call verinterp(ipmax,1,TM,X,TIME,F0)
      call verinterp(ipmax,1,TM,Y,TIME,L0)
      call verinterp(ipmax,1,TM,PRES,TIME,PRES1)
      call verinterp(ipmax,1,TM,PRES0,TIME,PRES2)
      DELP=(PRES2-PRES1)*100.
      call verinterp(ipmax,1,TM,WSPMAX,TIME,WSMAX)
      WSMAX=WSMAX*WIND_SCALE
      WS18=18*WIND_SCALE
      WS26=26*WIND_SCALE
      call verinterp(ipmax,1,TM,RMAXa,TIME,RMAX)
      RMAX=RMAX*1.e3
      do n=1,4
        call interp1d(-999.0,ipmax,1,TM,R18v(1,n),TIME,Rref18v(n))
        if(Rref18v(n).ne.-999) Rref18v(n) = Rref18v(n)*1.e3
        call interp1d(-999.0,ipmax,1,TM,R26v(1,n),TIME,Rref26v(n))
        if(Rref26v(n).ne.-999) Rref26v(n) = Rref26v(n)*1.e3
        alphv(n) = (n-1)*pi/2
      end do
      do n=2,6
        n1=mod(n-1,4)+1
        nm1=mod(n-2,4)+1
        np1=mod(n,4)+1
        if(Rref18v(n1).eq.-999) then
          if(Rref18v(nm1).ne.-999) then
            if(Rref18v(np1).ne.-999) then
              Rref18v(n1)=0.5*(Rref18v(nm1)+Rref18v(np1))
            else
              Rref18v(n1)=Rref18v(nm1)
            end if
          else
            if(Rref18v(np1).ne.-999) then
              Rref18v(n1)=Rref18v(np1)
            else
              Rref18v(n1)=-999
            end if
          end if
        end if
        if(Rref26v(n1).eq.-999) then
          if(Rref26v(nm1).ne.-999) then
            if(Rref26v(np1).ne.-999) then
              Rref26v(n1)=0.5*(Rref26v(nm1)+Rref26v(np1))
            else
              Rref26v(n1)=Rref26v(nm1)
            end if
          else
            if(Rref26v(np1).ne.-999) then
              Rref26v(n1)=Rref26v(np1)
            else
              Rref26v(n1)=-999
            end if
          end if
        end if
      end do
c
      Rref18v(5) = Rref18v(1)
      Rref26v(5) = Rref26v(1)
      alphv(5) = alphv(4)+pi/2
c
      print*,'Time=',Time
      print*,'Current hurricane position (x,y): ',f0,l0
      print*,'WSMAX=',WSMAX,'; DELP=',DELP,'; RMAX=',RMAX
      print*,'Rref18v=',Rref18v
      print*,'Rref26v=',Rref26v
c
      x0=f0
      y0=l0
      F0=F0*2.*PI/360
      L0=L0*2.*PI/360
C--- F0,L0 ARE THE CARRENT (LONGITUDE,LATTITUDE) COORDINATES OF THE HURRICANE
C
C--- CALCULATING UTX AND UTY (HURRICANE SPEED)
c
      cmp=tm(ipmax)
      do i=1,ipmax
        if(abs(tm(i)-time).le.cmp) then
          cmp=abs(tm(i)-time)
          ii=i
        end if
      end do
      if((tm(ii)-time) .le. 0. .and. ii .ne. ipmax) then
        t1=tm(ii)
        t2=tm(ii+1)
        x1=x(ii)
        x2=x(ii+1)
        y1=y(ii)
        y2=y(ii+1)
      else
        t2=tm(ii)
        t1=tm(ii-1)
        x2=x(ii)
        x1=x(ii-1)
        y2=y(ii)
        y1=y(ii-1)
      end if
      deltax1=rearth*cos(l0)*(x2-x1)*2.*pi/360
      deltay1=rearth*(y2-y1)*2.*pi/360
      utx=deltax1/((t2-t1)*24.*3600.)
      uty=deltay1/((t2-t1)*24.*3600.)
c
      print*,'utx,uty: ',utx,uty
C
C--- CALCULATING PARAMETERS FOR WIND PROFILE FORMULA
c
      B=WSMAX**2*E*ROA/DELP
      A7=RMAX**B
      print*,'B= ',B
c
      do i=1,14
        RADM(I)=RMAX*RAD(I)
      end do
C
      DO 350 J=1,JM
      DO 351 I=1,IM

cj<--- initialize arrays
         taux(i,j)=0.0
         tauy(i,j)=0.0
         wusurf(i,j)=0.0
         wvsurf(i,j)=0.0
         cor(i,j)=4*pi*sin(L0)/(3600*24)
cj---> initialize arrays
C  CALCULATING U-WIND STRESS FOR I,J POINT OF U-VELOCITY GRID
      F1=LONGMIN+(I-1)*(LONGMAX-LONGMIN)/((IM-0.5)-1)
      F1=F1*2.*PI/360
      L1=LATMIN+(J-1)*(LATMAX-LATMIN)/(JM-1)
      L1=L1*2.*PI/360
      DELTAX=REARTH*COS(L0)*(F1-F0)
      deltax=abs(deltax)
      if (deltax .eq. 0.) deltax=1.e-8
      DELTAY=REARTH*(L1-L0)
      deltay=abs(deltay)
      DXDY=DELTAX*DELTAY
      R=SQRT(DELTAX**2+DELTAY**2)
      alpha=atan(abs(DELTAY/DELTAX))*sign(1.,DXDY) +
     1      (1 - sign(1.,DXDY))*pi/2 + (1 - sign(1.,DELTAY))*pi/2
      if(alpha.ge.pi/4) then
        alpha = alpha - pi/4
      else
        alpha = alpha - pi/4 + 2*pi
      end if
      call verinterp(5,1,alphv,Rref18v,alpha,Rref18)
      call verinterp(5,1,alphv,Rref26v,alpha,Rref26)

      call verinterp(14,1,radm,angl,R,RANGL)

c      if(R.GT.RADM(14).OR.R.EQ.0) then
c        k=14
c      else
c        k=1
c        do while(R.GE.RADM(k+1))
c          k=k+1
c        end do
c      end if

C    CALCULATING WIND SPEED 

      if(Rref18.le.0.and.Rref26.le.0) then
        WND=SQRT(A7*B*DELP*EXP(-A7/R**B)/(ROA*R**B)+R**2*
     1  COR(I,J)**2/4.)-R*COR(I,J)/2.
        if(R.le.RMAX) then
           fmask=1.0
           UTXa=UTX/2.
           UTYa=UTY/2.

        else
           fmask=exp( (RMAX-R)/(2*RMAX) )
           UTXa=0.0
           UTYa=0.0
        endif

      else
        call EXPWND1(R,RMAX,Rref18,Rref26,WSMAX,WS18,WS26,wnd,fmask)
        UTXa=0.
        UTYa=0.
      end if
c      RANGL=ANGL(K)*PI/180.
      RANGL=RANGL*PI/180.
c      WF=WND
c      WR=-WND*TAN(RANGL)
      WF=WND*COS(RANGL)
      WR=-WND*SIN(RANGL)
      WX=WR*(DELTAX)/R-WF*(DELTAY)/R+UTXa
      WY=WF*(DELTAX)/R+WR*(DELTAY)/R+UTYa
      WM=SQRT(WX**2+WY**2)
cRMY      IF(WM.LT.10.) CD=1.14*1.E-3
cRMY      IF(WM.GE.10.) CD=(0.49+.065*WM)*1.E-3
cRMY      IF(WM.LE.35.) then
cRMY      WUSURF(I,J)=WUSURF(I,J)-CD*ROA*WM*WX/RHO_0
cRMY      else
cRMY      WUSURF(I,J)=WUSURF(I,J)-(3.3368+
cRMY     1      (WM-34.0449)**0.3)*WX/(WM*RHO_0)
cRMY      end if
cRMY
cRMY ***** NEW WIND STRESS PARAMETERIZATION 2****
      IF(WM.LE.12.5) Z0=(0.0185/9.80665)*(0.001*WM*WM+0.028*WM)**2.
      IF(WM.GT.12.5) Z0=(0.085*WM-0.58)/1000.
cRMY_92pct      CD=(0.4/log(10./Z0))**2.0
      CD=0.92*(0.4/log(10./Z0))**2.0
      WUSURF(I,J)=WUSURF(I,J)-CD*ROA*WM*WX/RHO_0
cRMY ********************************************

C  CALCULATING V-WIND STRESS FOR I,J POINT OF V-VELOCITY GRID
      F1=LONGMIN+(I-1)*(LONGMAX-LONGMIN)/(IM-1)
      F1=F1*2.*PI/360
      L1=LATMIN+(J-1)*(LATMAX-LATMIN)/((JM-0.5)-1)
      L1=L1*2.*PI/360
      DELTAX=REARTH*COS(L0)*(F1-F0)
      if (deltax .eq. 0.) deltax=1.e-8
      DELTAY=REARTH*(L1-L0)
      DXDY=DELTAX*DELTAY
      R=SQRT(DELTAX**2+DELTAY**2)
      alpha=atan(abs(DELTAY/DELTAX))*sign(1.,DXDY)+
     1      (1 - sign(1.,DXDY))*pi/2 + (1 - sign(1.,DELTAY))*pi/2
      if(alpha.ge.pi/4) then
        alpha = alpha - pi/4
      else
        alpha = alpha - pi/4 + 2*pi
      end if
      call verinterp(5,1,alphv,Rref18v,alpha,Rref18)
      call verinterp(5,1,alphv,Rref26v,alpha,Rref26)

      call verinterp(14,1,radm,angl,R,RANGL)

C    CALCULATING WIND SPEED 

      if(Rref18.le.0.and.Rref26.le.0) then
        WND=SQRT(A7*B*DELP*EXP(-A7/R**B)/(ROA*R**B)+R**2*
     1  COR(I,J)**2/4.)-R*COR(I,J)/2.
        UTXa=UTX/2.
        UTYa=UTY/2.

        if(R.le.RMAX) then
           fmask=1.0
           UTXa=UTX/2.
           UTYa=UTY/2.

        else
           fmask=exp( (RMAX-R)/(2*RMAX) )
           UTXa=0.0
           UTYa=0.0
        endif

      else
        call EXPWND1(R,RMAX,Rref18,Rref26,WSMAX,WS18,WS26,wnd,fmask)
        UTXa=0.
        UTYa=0.
      end if
c      RANGL=ANGL(K)*PI/180.
      RANGL=RANGL*PI/180.
      WF=WND*COS(RANGL)
      WR=-WND*SIN(RANGL)
      WX=WR*(DELTAX)/R-WF*(DELTAY)/R+UTXa
      WY=WF*(DELTAX)/R+WR*(DELTAY)/R+UTYa
      WM=SQRT(WX**2+WY**2)
cRMY      IF(WM.LT.10.) CD=1.14*1.E-3
cRMY      IF(WM.GE.10.) CD=(0.49+.065*WM)*1.E-3
cRMY      IF(WM.le.35.) then
cRMY      WVSURF(I,J)=WVSURF(I,J)-CD*ROA*WM*WY/RHO_0
cRMY      else
cRMY      WVSURF(I,J)=WVSURF(I,J)-(3.3368+
cRMY     1      (WM-34.0449)**0.3)*WY/(WM*RHO_0)
cRMY      end if
cRMY
cRMY ***** NEW WIND STRESS PARAMETERIZATION 2****
      IF(WM.LE.12.5) Z0=(0.0185/9.80665)*(0.001*WM*WM+0.028*WM)**2.
      IF(WM.GT.12.5) Z0=(0.085*WM-0.58)/1000.
cRMY_92pct      CD=(0.4/log(10./Z0))**2.0
      CD=0.92*(0.4/log(10./Z0))**2.0
      WVSURF(I,J)=WVSURF(I,J)-CD*ROA*WM*WY/RHO_0
cRMY ********************************************

C  CALCULATING WIND STRESS FOR I,J POINT OF DEPTH GRID
      F1=LONGMIN+(I-1)*(LONGMAX-LONGMIN)/(IM-1)
      F1=F1*2.*PI/360.
      L1=LATMIN+(J-1)*(LATMAX-LATMIN)/(JM-1)
      L1=L1*2.*PI/360.
      DELTAX=REARTH*COS(L0)*(F1-F0)
      if (deltax .eq. 0.) deltax=1.e-8
      DELTAY=REARTH*(L1-L0)
      DXDY=DELTAX*DELTAY
      R=SQRT(DELTAX**2+DELTAY**2)
      alpha=atan(abs(DELTAY/DELTAX))*sign(1.,DXDY)+
     1      (1 - sign(1.,DXDY))*pi/2 + (1 - sign(1.,DELTAY))*pi/2
      if(alpha.ge.pi/4) then
        alpha = alpha - pi/4
      else
        alpha = alpha - pi/4 + 2*pi
      end if
      call verinterp(5,1,alphv,Rref18v,alpha,Rref18)
      call verinterp(5,1,alphv,Rref26v,alpha,Rref26)

      call verinterp(14,1,radm,angl,R,RANGL)

C    CALCULATING WIND SPEED 

      if(Rref18.le.0.and.Rref26.le.0) then
        WND=SQRT(A7*B*DELP*EXP(-A7/R**B)/(ROA*R**B)+R**2*
     1  COR(I,J)**2/4.)-R*COR(I,J)/2.
        UTXa=UTX/2.
        UTYa=UTY/2.
        if(R.le.RMAX) then
           fmask=1.0
           UTXa=UTX/2.
           UTYa=UTY/2.

        else
           fmask=exp( (RMAX-R)/(2*RMAX) )
           UTXa=0.0
           UTYa=0.0
        endif
      else
        call EXPWND1(R,RMAX,Rref18,Rref26,WSMAX,WS18,WS26,wnd,fmask)
        UTXa=0.
        UTYa=0.
      end if
c      RANGL=ANGL(K)*PI/180.
      RANGL=RANGL*PI/180.
      WF=WND*COS(RANGL)
      WR=-WND*SIN(RANGL)
      WX=WR*(DELTAX)/R-WF*(DELTAY)/R+UTXa
      WY=WF*(DELTAX)/R+WR*(DELTAY)/R+UTYa
c      WINDX(i,j)=WX
c      WINDY(i,j)=WY
      WM=SQRT(WX**2+WY**2)
cRMY      IF(WM.LT.10.) CD=1.14*1.E-3
cRMY      IF(WM.GE.10.) CD=(0.49+.065*WM)*1.E-3
cRMY      IF(WM.GT.35.) THEN
cRMY      TMAX=3.3368+(WM-34.0449)**0.3
cRMY      TAUX(I,J)=TAUX(I,J)+TMAX*WX/WM
cRMY      TAUY(I,J)=TAUY(I,J)+TMAX*WY/WM
cRMY      ELSE
cRMY      TAUY(I,J)=TAUY(I,J)+CD*ROA*WM*WY
cRMY      TAUX(I,J)=TAUX(I,J)+CD*ROA*WM*WX
cRMY
cRMY ***** NEW WIND STRESS PARAMETERIZATION 2****
      IF(WM.LE.12.5) Z0=(0.0185/9.80665)*(0.001*WM*WM+0.028*WM)**2.
      IF(WM.GT.12.5) Z0=(0.085*WM-0.58)/1000.
cRMY_92pct      CD=(0.4/log(10./Z0))**2.0
      CD=0.92*(0.4/log(10./Z0))**2.0
      TAUY(I,J)=TAUY(I,J)+CD*ROA*WM*WY
      TAUX(I,J)=TAUX(I,J)+CD*ROA*WM*WX
cRMY ********************************************

C      WTSURF(I,J)=30.*WM/(4000*1024.)

cRMY      END IF

      UWIND(I,J)=WX
      VWIND(I,J)=WY

c<--- HSK 10/25/2007      
      WNDSPD(I,J)=SQRT(WX**2 + WY**2)
      fmeld(i,j)=fmask
c---> HSK
      
 351  CONTINUE
 350  CONTINUE

c      print*,'Exiting WIND ...'
c      
c      OPEN(39,FILE='TXY.dat',STATUS='UNKNOWN',form='unformatted')
c      WRITE(39) WUSURF
c      WRITE(39) WVSURF
c      CLOSE(39)
      

c      OPEN(49,FILE='TXY0.dat',STATUS='UNKNOWN',form='unformatted')
c      WRITE(49) TAUX
c      WRITE(49) TAUY
c      CLOSE(49)
      
      
c      OPEN(59,FILE='WIND.dat',STATUS='UNKNOWN',form='unformatted')
c      WRITE(59) UWIND
c      WRITE(59) VWIND
c      CLOSE(59)
           
      RETURN
      END
c==============================================================================
        SUBROUTINE DATE2DAY(year,julday,date)
      integer*4 date
      integer dat2day(12),dat2dayl(12),day,month,year,hour
      real julday
      real*8 tmp
      data dat2day/31,28,31,30,31,30,31,31,30,31,30,31/
      data dat2dayl/31,29,31,30,31,30,31,31,30,31,30,31/

      year=int(date/1000000.)
      month=nint(100*(date/1000000.-int(date/1000000.)))
      julday=0
        if(mod(year,4).eq.0) then
      do n=1,month-1
      julday=julday+dat2dayl(n)
      end do
        else
      do n=1,month-1
      julday=julday+dat2day(n)
      end do
        end if
      julday=julday+nint(100*(date/10000.-int(date/10000.)))
      hour=date-nint(date/100.)*100
      julday=julday+float(hour)/24.

      return
      end
c==============================================================================
        SUBROUTINE DAY2DATE(year,julday,date)
      integer*4 date
      integer dat2day(12),dat2dayl(12),day,month,year,year1,hour
      real julday,julday1
      real*8 tmp
      data dat2day/31,28,31,30,31,30,31,31,30,31,30,31/
      data dat2dayl/31,29,31,30,31,30,31,31,30,31,30,31/

      if(int(julday).gt.365+int(1./(mod(year,4)*100.+1.))) then
      julday1=julday-365-int(1./(mod(year,4)*100.+1.))
      year1=year+1
      else
      julday1=julday
      year1=year
      end if
      day=0
      n=1
        if(mod(year1,4).eq.0) then
      do while(day+dat2dayl(n).lt.int(julday1))
      day=day+dat2dayl(n)
      n=n+1
      end do
        else
      do while(day+dat2day(n).lt.int(julday1))
      day=day+dat2day(n)
      n=n+1
      end do
        end if
      month=n
      day=int(julday1-day)
      hour=nint((julday1-int(julday1))*24.)
      date=year1*1000000+month*10000+day*100+hour
c
      return
      end
c==============================================================================
      SUBROUTINE VERINTERP(n,ni,x,y,xi,yi)
c--------------------------------------------------------------
c   This subroutine determines ni values yi at the points xi
c   interpolating between the n values y at the points x
c--------------------------------------------------------------
      real x(n),y(n),xi(ni),yi(ni),tmp(500)
      real cmp
      integer n,ni,ii

            do i=1,ni
      if((xi(i).gt.x(n).and.xi(i).gt.x(1)).or.
     *   (xi(i).lt.x(1).and.xi(i).lt.x(n))) then
        if(xi(i).gt.x(n).and.xi(i).gt.x(1)) then
            if(x(n).gt.x(1)) then
              yi(i)=y(n)
            else
              yi(i)=y(1)
            end if
        else
            if(x(n).gt.x(1)) then
              yi(i)=y(1)
            else
              yi(i)=y(n)
            end if
        end if
      else
        do j=1,n-1
        tmp(j)=(xi(i)-x(j))*(xi(i)-x(j+1))
        end do
c
         do j=1,n-1
           if(tmp(j).le.0) then
             ii=j
           end if
         end do
c
        yi(i)=(y(ii)*abs(x(ii+1)-xi(i))+y(ii+1)*abs(xi(i)-x(ii)))
     1        /abs(x(ii+1)-x(ii))
      end if
            end do

      return
      end
c==============================================================================
      SUBROUTINE INTERP1d(mask,n,ni,x,y,xi,yi)
c--------------------------------------------------------------
c   This subroutine determines ni values yi at the points xi
c   interpolating between the n values y at the points x
c   values equal to mask are ignored
c--------------------------------------------------------------
      real x(n),y(n),xi(ni),yi(ni),tmp(500)
      real cmp,mask
      integer n,ni,ii
c
      do i=1,ni
        if((xi(i).gt.x(n).and.xi(i).gt.x(1)).or.
     *  (xi(i).lt.x(1).and.xi(i).lt.x(n))) then
          if(xi(i).gt.x(n).and.xi(i).gt.x(1)) then
            if(x(n).gt.x(1)) then
              yi(i)=y(n)
            else
              yi(i)=y(1)
            end if
          else
            if(x(n).gt.x(1)) then
              yi(i)=y(1)
            else
              yi(i)=y(n)
            end if
          end if
        else
          do j=1,n-1
            tmp(j)=(xi(i)-x(j))*(xi(i)-x(j+1))
          end do
          do j=1,n-1
            if(tmp(j).le.0) ii=j
          end do
          if(y(ii).eq.mask.or.y(ii+1).eq.mask) then
            yi(i)=mask
          else
            yi(i)=(y(ii)*abs(x(ii+1)-xi(i))+y(ii+1)*abs(xi(i)-x(ii)))/
     1      abs(x(ii+1)-x(ii))
          end if
        end if
      end do
c
      return
      end

c==============================================================================
      subroutine  EXPWND1(R,RMAX,Rref18,Rref26,WSMAX,WS18,WS26,
     $     expwnd,fmask)
      real Rref18,Rref26,R,RMAX,WSMAX,WS18,WS26
      real b,expwnd
      real fmask
c
       r1=0.5*(Rref18+Rref26)
       WS=0.5*(WS18+WS26)
       if(Rref18.le.0.) then
       r1=Rref26
       WS=WS26
       end if
       if(Rref26.le.0.) then
       r1=Rref18
       WS=WS18
       end if
c
         if(R.GE.RMAX) then
         b=(RMAX-r1)/log(WS/WSMAX)
         expwnd=WSMAX*exp((RMAX-R)/b)
         fmask=exp( (RMAX-R)/b )
         else
         expwnd=R*WSMAX/RMAX
         fmask=1.0
         end if
c
      return
      end

