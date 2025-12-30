      program rivers
      use mod_za  ! HYCOM I/O interface
      implicit none
c
c     extract a subregion flds from a full region HYCOM 2.0 flds
c
      character*79         :: preambl(5)
      character*80         :: cline,fld_name
      character*79         :: cline_out
      character*128        :: flnm_in,flnm_fld,stat_in
      integer              :: i,j,l,ni,nf,k,ns
      real                 :: hmina,hminb,hmaxa,hmaxb,rtmp,val
      integer              :: num_stations
      
      integer, allocatable :: m_in(:,:)
      real,    allocatable :: a_in(:,:)
      real,    allocatable :: st(:,:)
      real,    allocatable :: line(:)
      integer, allocatable :: iw(:),jw(:)
c
      real,    parameter   :: hspval=0.5*2.0**100
      real,    parameter   :: spval=     2.0**100

c
      call xcspmd

      allocate( m_in(idm,jdm) )
      allocate( a_in(idm,jdm) )

      call zaiost
      call blkdat(num_stations,stat_in,fld_name,
     &            flnm_in,flnm_fld,cline_out)

      allocate( st(num_stations,12) )
      allocate (line(14))
      allocate( iw(num_stations) )
      allocate( jw(num_stations) )
      



c
c
c     open input and output files.
c
c
      ns = 12
      l  = len_trim(stat_in)
      open (unit=ns,file= stat_in(1:l),form='formatted',
     .      status='old',action='read')
      write(6,'(a,a,a,i4)') 
     &   'opened file: ',stat_in(1:l),'  on unit: ',ns

      ni = 14
      l  = len_trim(flnm_in)
      open (unit=ni,file= flnm_in(1:l)//'.b',form='formatted',
     .      status='old',action='read')
      call zaiopf( flnm_in(1:l)//'.a','old', ni)
      write(6,'(a,a,a,i4)') 
     &   'opened file: ',flnm_in(1:l)//'.a','  on unit: ',ni
c
      nf = 15
      l  = len_trim(flnm_fld)
      open (unit=nf,file= flnm_fld(1:l)//'.b',form='formatted',
     .      status='new',action='write')
      call zaiopf( flnm_fld(1:l)//'.a','new', nf)
      write(6,'(a,a,a,i4)') 
     &   'opened file: ',flnm_fld(1:l)//'.a','  on unit: ',nf

c
c    load river outflow monthly data
c
      do k=1,num_stations
        read(ns,*)line
        iw(k)=line(1)
        jw(k)=line(2)
        do i=1,12
          st(k,i)=line(2+i)
        enddo
      enddo
      close(ns)
           

c
c     process the header
c
      read( ni,'(a79)') preambl
      write(lp,*)
      write(lp,'(a79)') preambl
      write(lp,*)
      call flush(lp)
c

c
c     input the domain bathymetry.
c

      read( ni,'(a)') cline
      l = index(cline,'=')
      read (cline(l+1:),*)  hminb,hmaxb
      call zaiord(a_in,m_in,.false., hmina,hmaxa, ni)
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - full domain .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        stop
      endif
c
c     construct mask
c
      do j= 1,jdm
         do i= 1,idm
            if     (a_in(i,j).lt.hspval) then
               m_in(i,j) = 1
            else
               m_in(i,j) = 0
            endif
         enddo
      enddo
c

c
c--- write header
c
      preambl(1)=cline_out
      write(nf,'(a)')    preambl

c 
c -- write monthly averages
c

      do k=1,12
        do j= 1,jdm
          do i= 1,idm
            if     (m_in(i,j).eq.1) then
              a_in(i,j) = 0.0
            else
              a_in(i,j) = spval
            endif
          enddo
        enddo
        do l=1,num_stations
          i=iw(l)
          j=jw(l)
          val=st(l,k)
          a_in(i,j)=val+a_in(i,j)
        enddo
c
        call zaiowr(a_in,m_in,.true., hmina,hmaxa, nf, .false.)

        l = index(cline,':')
        fld_name=cline(1:l)

        write(nf,'(a,i5,2e18.7)')
     $       fld_name(1:l)//' month, range = ',
     $       k,hmina,hmaxa
        call flush(nf)
        write(lp,'(a,i5,2e18.7)')'min,max water flux (m/s) = ',
     $       k,hmina,hmaxa

        call flush(lp)
      enddo

      stop

c
      end program rivers

      subroutine blkdat(num_stations,stat_in,fld_name,
     &                  flnm_in,flnm_fld,cline_out)
      use mod_xc  ! HYCOM communication interface
      implicit none
      character*128 :: flnm_in,flnm_out,flnm_ref,flnm_fld,stat_in
      character*79  :: cline_out
      character*80  :: fld_name

      integer       :: num_stations
c
c --- read blkdat.input for topog subregion.
c
      integer       :: irefi,irefo,jrefi,jrefo
c
c --- 'fld_name'  = name of the field
c --- 'stat_in'   = input file with river ouflows
c --- 'flnm_in'   = input   region bathymetry filename
c --- 'flnm_fld'  = output   region fld filename
c --- 'cline_out' = output title line (replaces preambl(5))
c
      read( *,'(a)')      fld_name
      write(6,'(a)') trim(fld_name)
      read( *,'(a)')      stat_in
      write(6,'(a)') trim(stat_in)
      read( *,'(a)')      flnm_in
      write(6,'(a)') trim(flnm_in)
      read( *,'(a)')      flnm_fld
      write(6,'(a)') trim(flnm_fld)
      read( *,'(a)') cline_out
      write(6,'(a)') cline_out
      write(6,*)
c --- 'num-stations' = number of stations in stat_in file

      call blkini(num_stations,'num-stations')

      call flush(6)
      end subroutine blkdat

      subroutine blkinr(rvar,cvar,cfmt)
      implicit none
c
      real      rvar
      character cvar*6,cfmt*(*)
c
c     read in one real value
c
      character*6 cvarin
c
      read( *,*) rvar,cvarin
      write(6,cfmt) cvarin,rvar
      call flush(6)
c
      if     (cvar.ne.cvarin) then
        write(6,*) 
        write(6,*) 'error in blkinr - input ',cvarin,
     +                      ' but should be ',cvar
        write(6,*) 
        call flush(6)
        stop
      endif
      return
      end subroutine blkinr

      subroutine blkini(ivar,cvar)
      implicit none
c
      integer     ivar
      character*6 cvar
c
c     read in one integer value
c
      character*6 cvarin
c
      read( *,*) ivar,cvarin
      write(6,6000) cvarin,ivar
      call flush(6)
c
      if     (cvar.ne.cvarin) then
        write(6,*) 
        write(6,*) 'error in blkini - input ',cvarin,
     +                      ' but should be ',cvar
        write(6,*) 
        call flush(6)
        stop
      endif
      return
 6000 format(a6,' =',i6)
      end subroutine blkini

      subroutine blkinl(lvar,cvar)
      implicit none
c
      logical     lvar
      character*6 cvar
c
c     read in one logical value
c     due to a SGI bug for logical I/O: read in an integer 0=F,1=T
c
      character*6 cvarin
      integer     ivar
c
      read( *,*) ivar,cvarin
      lvar = ivar .ne. 0
      write(6,6000) cvarin,lvar
      call flush(6)
c
      if     (cvar.ne.cvarin) then
        write(6,*) 
        write(6,*) 'error in blkinr - input ',cvarin,
     +                      ' but should be ',cvar
        write(6,*) 
        call flush(6)
        stop
      endif
      return
 6000 format(a6,' =',l6)
      end subroutine blkinl
