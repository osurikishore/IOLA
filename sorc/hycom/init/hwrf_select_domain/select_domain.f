      program select_domain
      implicit none
      integer,parameter::idm=4500,jdm=3298
      real lon(idm,jdm),lat(idm,jdm)
      integer num_dom,nd(99),is(99),js(99),ie(99),je(99)
      integer in(99),fin(99)
      real lonmin(99),lonmax(99),latmin(99),latmax(99)
      character*11 app(99)
      integer num_pt, hsf
      real lont(99),latt(99),t(99),x(99),y(99),d0
      integer ix(99),jy(99),iout(99)
      integer k,l,mark,iu,debug,prev,same,basin
      debug=1

c
c -- get previous domain identifier
c
      iu=8
      read(iu,*)prev, basin
      close(iu)
      write(6,*)'Previous domain number identifier', prev
      write(6,*)'Storm basin number identifier', basin
c
c -- get subdomain boundaries in the big domain
c
      iu=10
      call get_domains(num_dom,nd,is,js,ie,je,app,
     &             lonmin,lonmax,latmin,latmax,iu)
      if(debug.gt.0) then
        write(6,*)num_dom,' sub-domains:'
        do l=1,num_dom
         write(6,'(5i6,4(1x,f10.5))')nd(l),is(l),js(l),ie(l),je(l),
     &         lonmin(l),lonmax(l),latmin(l),latmax(l)
        enddo
      endif
c
c -- get track
c      
      iu=12
      call get_track(num_pt,lont,latt,t,iu)
      
c 
c -- convert lon if WPac
c
      if(basin.eq.3) then
        do k=1,num_pt
          lont(k)=lont(k)*(-1)-360
        enddo
      endif

      if(debug.gt.0) then
        do k=1,num_pt
          write(6,*)lont(k),latt(k)
        enddo
      endif

c
c -- determine scores for points inside the domain
c
      if(num_pt.gt.0) then
         
        write(6,'(a,a)')
     &        ' dom no., no. track points in,'
        do l=1,num_dom
         if(basin.eq.0.or.basin.eq.nd(l)/10) then
          in(l)=0
          fin(l)=0
          do k=1,num_pt
            if(lont(k)+360.le.lonmax(l).and.
     &           lont(k)+360.ge.lonmin(l).and.
     &           latt(k).le.latmax(l).and.
     &           latt(k).ge.latmin(l))then
                in(l)=in(l)+1
                if (k.eq.1) then
                    fin(l)=1
                endif
            endif
          enddo
         endif 
         write(6,'(2i6)')nd(l),in(l)
        enddo
c
c -- select domain containing the entire track, and if it was
c    the last domain; otherwise select the domain with the
c    most points in the domain. Preference is first found.
c

        mark=0
        hsf=0
        do l=1,num_dom

c check if previous domain has all points
          if(nd(l).eq.prev.and.
     &       in(l).eq.num_pt) then 
              mark=l
              exit
          endif
c check if this domain looks promising
          if(in(l).gt.hsf) then
              hsf=in(l)
              mark=l
          endif

c check if this domain looks even better
          if(in(l).gt.hsf.and.hsf.ne.0) then
              hsf=in(l)
              mark=l
            endif
c         endif

        enddo

        same=0
        if (mark.eq.prev) then
          same=1
        endif
c
c --> hsk's doing to keep previous domain when points are tied.
        if (hsf.eq.in(prev)) then
           same=1
           mark=prev
        endif
c <-- hsk's

        if(in(mark).ge.1) then
          write(15,'(a,i2)')'export domain_no=',nd(mark)
          write(15,*)'export Application=',trim(app(mark))
          if(prev.eq.0) then
            write(15,*)'export same_domain=-1'
          endif
          if(same.eq.0.and.prev.ne.0) then
            write(15,*)'export same_domain=0'
          endif
          if(same.eq.1.and.prev.ne.0) then
            write(15,*)'export same_domain=1'
          endif
          write(6,*)'Select domain ',nd(mark)
          
        else
          write(15,'(a,i1)')'export domain_no=',0
          write(15,*)'export Application=NONE'
          write(6,*)'Select NO domain '

        endif
         
      else
        write(6,*)'Track data is missing'
        write(15,'(a,i2)')'export domain_no=',-1
        write(15,*)'export Application=NONE'
      endif

      close(15)
      stop
      end

      subroutine get_domains(num_dom,nd,is,js,ie,je,app,
     &                    lonmin,lonmax,latmin,latmax,iu)
      implicit none
      integer num_dom,is(99),js(99),ie(99),je(99),nd(99)
      real lonmin(99),lonmax(99),latmin(99),latmax(99)
      character*11 app(99)
      integer k,iu,ios,nd_,is_,js_,ie_,je_
      real lonmin_,lonmax_,latmin_,latmax_
      character*11 app_
      num_dom=0
      do 
      read(iu,'(i2,4i5,1x,a,4f5.0)',iostat=ios)nd_,is_,js_,ie_,je_,app_
     &                                ,lonmin_,lonmax_,latmin_,latmax_
         if(ios.ne.0)exit
         num_dom=num_dom+1
         nd(num_dom)=nd_
         is(num_dom)=is_
         ie(num_dom)=ie_
         js(num_dom)=js_
         je(num_dom)=je_
         app(num_dom)=trim(app_)
         lonmin(num_dom)=lonmin_
         lonmax(num_dom)=lonmax_
         latmin(num_dom)=latmin_
         latmax(num_dom)=latmax_
      enddo
      close(iu)
      return
      end
      subroutine get_track(num_pt,lont,latt,t,iu)
      implicit none
      integer num_pt,iu
      real lont(99),latt(99),t(99)
      integer ios
      real dumlon,dumlat,dumt
      num_pt=0
      do 
         read(iu,*,iostat=ios)dumt,dumlon,dumlat
         if(ios.ne.0) exit
         num_pt=num_pt+1
         lont(num_pt)=dumlon
         latt(num_pt)=dumlat
         t(num_pt)=dumt
      enddo
      close(iu)
      return
      end

      subroutine get_grid(lon,lat,idm,jdm,iu)
      integer idm,jdm,iu
      real lon(idm,jdm),lat(idm,jdm)
      integer i,j
      do j=1,jdm
        do i=1,idm
          read(iu,*)lon(i,j),lat(i,j)
        enddo
      enddo
      write(6,*)lon(1,1),lon(idm,jdm),' lon at 1,1 idm,jdm'
      write(6,*)lat(1,1),lat(idm,jdm),' lat at 1,1 idm,jdm'
      return
      end

