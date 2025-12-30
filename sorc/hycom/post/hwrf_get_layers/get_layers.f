      program get_layers
      implicit none
      integer k1,k2
      real nav,bdiff,cdiff,biso
      real             ::hnav
      real, allocatable, dimension(:) ::a_in
      real, allocatable, dimension(:) ::b_in
      real, allocatable, dimension(:) ::c_in

      real, allocatable, dimension(:) ::d_in        ! mix_depth
      real, allocatable, dimension(:) ::sst         ! SST
      integer, allocatable, dimension(:) ::iwet	    ! wet point

      real, allocatable, dimension(:) ::ohc

      real, allocatable::b(:),bt(:),bd(:),bs(:),bz(:)
      real, allocatable::c(:),ct(:),cd(:),cs(:),c_at_bz(:)
      real, allocatable::a(:),at(:),ab(:),a1(:),a2(:)
      real  one,zero,p
      integer :: i,j,ip,k,kp,imjm,ios,idum,jdum,im,jm,idm,jdm,i1,j1
      integer,parameter::ic=1,jc=1
      integer,parameter::i0=1,j0=1
      character*100 fil
      integer indx
      real pref,dz,tw
 
      real, parameter:: s_heat=3850, rho_w=1025

      one=1.0
      zero=0.0
      k=0
      read(5,*) idm,jdm,k1,k2,nav,bdiff,cdiff,biso,pref
      hnav=0.5*nav
      im=idm
      jm=jdm
      i1=im
      j1=jm
      imjm=im*jm
      allocate (a_in(imjm), b_in(imjm), c_in(imjm), d_in(imjm))
      allocate(a(imjm),at(imjm),ab(imjm),a1(imjm),a2(imjm))
      allocate(b(imjm),bt(imjm),bd(imjm),bs(imjm),bz(imjm))
      allocate(c(imjm),ct(imjm),cd(imjm),cs(imjm),c_at_bz(imjm))
      
      allocate (sst(imjm))
      allocate (ohc(imjm))
      allocate (iwet(imjm))

c -- initialize
      do kp=1,imjm
c interface
         a(kp)=zero
         at(kp)=zero
         ab(kp)=zero
         a1(kp)=zero
         a2(kp)=zero
c scalar 1: WTMPC
         b(kp)=zero
         bt(kp)=zero
         bd(kp)=nav
         bs(kp)=nav
         bz(kp)=nav
c scalar 2: DENSITY
         c(kp)=zero
         ct(kp)=zero
         cd(kp)=nav
         cs(kp)=zero
         c_at_bz(kp)=nav
c wet point: 
         iwet(kp)=1
c ohc 
         ohc(kp)=zero
         sst(kp)=nav
      enddo
      
      do
         read(5,'(a)',iostat=ios)fil
     
         if(ios.ne.0) exit
         k=k+1

c -- read record for VAR1 (LAYTH)
         open(10,file=fil,form='formatted')
         read(10,*)idum,jdum
         do j=1,jm
         do i=1,im
           ip=i+(j-1)*im
           read(10,*)a_in(ip)

c <-- identify wet point (just one time run)
           if(a_in(ip).eq.nav)then
              iwet(ip)=0
           endif
c -->
         enddo
         enddo
         close(10)

         read(5,'(a)',iostat=ios)fil
         if(ios.ne.0) exit
c -- read record for VAR2 (WTMPC)
         open(11,file=fil,form='formatted')
         read(11,*)idum,jdum
         do j=1,jm
         do i=1,im
         ip=i+(j-1)*im
         read(11,*)b_in(ip)
           if(k.eq.1)then
              if(b_in(ip).ne.nav)then
                 sst(ip)=b_in(ip)
              endif
           endif
         enddo
         enddo
         close(11)

         read(5,'(a)',iostat=ios)fil
         if(ios.ne.0) exit
c -- read record for VAR3 (DENS)
         open(12,file=fil,form='formatted')
         read(12,*)idum,jdum
         do j=1,jm
         do i=1,im
         ip=i+(j-1)*im
         read(12,*)c_in(ip)
         enddo
         enddo
         close(12)

         read(5,'(a)',iostat=ios)fil
         if(ios.ne.0) exit
c -- read record for VAR4 (MIX_DPTH)
         open(13,file=fil,form='formatted')
         read(13,*)idum,jdum
         do j=1,jm
         do i=1,im
         ip=i+(j-1)*im
         read(13,*)d_in(ip)
         enddo
         enddo
         close(13)

         do kp=1,imjm
            if(iwet(kp).eq.1) then

c -- interface
               at(kp)=a(kp)
               a(kp)=a(kp)+a_in(kp) !layer depth
               if(k.eq.k1) a1(kp)=a(kp)
               if(k.eq.k2) a2(kp)=a(kp)
c -- scalars
c -- scalar 2: DENSITY
               ct(kp)=c(kp)
               c(kp)=c_in(kp)
               if(k.eq.1) then
                  cs(kp)=c(kp)
               else
                 if(cd(kp).eq.nav) then
                   if(c(kp)-cs(kp).ge.cdiff) then
                     p=(cs(kp)-ct(kp)-cdiff)/(c(kp)-ct(kp))
                     p=min(max(p,zero),one)
                     cd(kp)=at(kp)+p*(a(kp)-at(kp))
                   endif
                 endif
               endif

c -- scalar 1: WTMPC
               bt(kp)=b(kp)
               b(kp)=b_in(kp)
               if(k.eq.1) then
c                bs(kp)=b(kp)
               else
                 if(bd(kp).eq.nav) then
                   if(bs(kp)-b(kp).ge.bdiff) then
                     p=(bs(kp)-bt(kp)-bdiff)/(b(kp)-bt(kp))
                     p=min(max(p,zero),one)
                     bd(kp)=at(kp)+p*(a(kp)-at(kp))
                   endif
                 endif
                 if(bz(kp).eq.nav) then
                   if(b(kp).le.biso.and.bt(kp).gt.biso) then
                     p=(biso-b(kp))/(bt(kp)-b(kp))
                     p=min(max(p,zero),one)
                     bz(kp)=at(kp)+p*(a(kp)-at(kp))
                     c_at_bz(kp)=ct(kp)+p*(c(kp)-ct(kp))
                   elseif(b(kp).ge.biso.and.bt(kp).lt.biso) then
                     p=(biso-b(kp))/(bt(kp)-b(kp))
                     p=min(max(p,zero),one)
                     bz(kp)=at(kp)+p*(a(kp)-at(kp))
                     c_at_bz(kp)=ct(kp)+p*(c(kp)-ct(kp))
                   endif
                 endif
               endif


c -- scalar x: OHC
              tw=(b(kp)+bt(kp))/2.
              if (tw.ge.26.) then
                 dz=a(kp)-at(kp)
                 ohc(kp)=ohc(kp)+s_heat*rho_w*abs(dz)*(tw-26.)
              endif      


             endif		! if for wet
           enddo		! do for kp
         enddo                  ! DO

         write(6,*)k,' no of frames'
         write(6,*)'last file '
         write(6,'(a)')fil
c
c  -- write results
c

c         write(20,'(8(e18.6))')
c     $     sst(kp),cs(kp),d_in(kp),bd(kp),cd(kp),bz(kp),c_at_bz(kp),
c     $     ohc(kp)   
c-----------------------------------------------------------------------------------------------
c 1      2        3           4               5                  6    7          8
c sst, ssdensity, model-MLD,  MLD w/ dT=0.25, MLD w/ dSigma=0.1, Z26, Sigma @26. OHC*heat*rho
c-----------------------------------------------------------------------------------------------

      open(21,status='unknown',form='unformatted')
      write (21) sst,cs,d_in,bd,cd,bz,c_at_bz,ohc
      write(6,*)k,' no of frames'
      write(6,*)'last file '
      write(6,'(a)')fil
      stop
      end program get_layers
