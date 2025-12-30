#define byteswap8(v) \
  ior(iand(ishft(v,-56),Z'FF'),\
  ior(iand(ishft(v,-40),Z'FF00'),\
  ior(iand(ishft(v,-24),Z'FF0000'),\
  ior(iand(ishft(v, -8),Z'FF000000'),\
  ior(iand(ishft(v,  8),Z'FF00000000'),\
  ior(iand(ishft(v, 24),Z'FF0000000000'),\
  ior(iand(ishft(v, 40),Z'FF000000000000'),\
      iand(ishft(v, 56),Z'FF00000000000000') )))))))

#define char3toint8(c,i) \
  ior(ishft(ichar(c((i):(i)),kind=8),16), \
  ior(ishft(ichar(c((i)+1:(i)+1),kind=8),8), \
  ichar(buf((i)+2:(i)+2),kind=8)))

#define char2toint8(c,i) \
  ior(ishft(ichar(c((i):(i)),kind=8),8),ichar(c((i)+1:(i)+1),kind=8))

module module_regrid_merge
  implicit none

  logical, parameter :: big_endian = ichar(transfer(1,'1')) == 0

  logical(kind=1), pointer, private ::  bitmask_storage(:) => NULL()

  integer, parameter :: maxbit=32

  type gribrecord
     logical(kind=1), pointer :: mask(:) => NULL()
     real, pointer :: data(:)    => NULL()
     character(len=:), pointer :: cgds => NULL()
     character(len=:), pointer :: cpds => NULL()
     integer, pointer :: gds(:) => NULL()
     integer, pointer :: gdsf(:) => NULL()
     integer, pointer :: pds(:) => NULL()
     integer(kind=8) :: secloc(0:4)=(/0,0,0,0,0/)
     integer(kind=8) :: reclen=0
     integer :: npoints=0,ilengds=0,ilengdsf=0
  end type gribrecord

  type gribfile
     type(gribrecord), pointer :: records(:) => NULL()
     character(len=:), pointer :: filename => NULL()
     logical, pointer :: done(:) => NULL()
     integer :: unit=0, nrecords=0
     integer(kind=8) :: writepos=0
  end type gribfile

  type regrid_merge
     type(gribfile), pointer :: files(:) => NULL()
     integer, pointer :: gds(:)
     character(len=:), pointer :: cgds
     integer :: gdslen=0, nopened=0, igrid=0, gridwind=0
  end type regrid_merge

  ! --------------------------------------------------------------------

contains

  ! --------------------------------------------------------------------

  logical function same_pds(cpds1,cpds2)
    character(len=*), intent(in) :: cpds1,cpds2

    same_pds = cpds1(4:4)==cpds2(4:4) .and. cpds1(9:21)==cpds2(9:21) &
         .and. cpds1(6:6)==cpds2(6:6) .and. len(cpds1)==len(cpds2)

    if(same_pds) then
       if(len(cpds1)>40) then
          same_pds=cpds1(41:)==cpds2(41:)
       endif
    endif
  end function same_pds

  ! --------------------------------------------------------------------

  logical function other_vector_component(cpds1,cpds2)
    character(len=*), intent(in) :: cpds1,cpds2
    logical :: same_pds

    other_vector_component=(ichar(cpds1(9:9))==33 .and. ichar(cpds2(9:9))==34) &
         .or. (ichar(cpds1(9:9))==34 .and. ichar(cpds2(9:9))==33)

    if(.not.other_vector_component) return

    same_pds = cpds1(4:4)==cpds2(4:4) .and. cpds1(10:21)==cpds2(10:21) &
         .and. cpds1(6:6)==cpds2(6:6) .and. len(cpds1)==len(cpds2)

    if(same_pds) then
       if(len(cpds1)>40) then
          same_pds=cpds1(41:)==cpds2(41:)
       endif
    endif

    other_vector_component=same_pds

  end function other_vector_component

  ! --------------------------------------------------------------------

  subroutine copy_pds(cf,crec,of,orec,igrid)
    implicit none
    type(gribfile), pointer :: cf,of
    integer, intent(in) :: crec,orec,igrid
    type(gribrecord), pointer :: cr,or
    integer :: i

    cr=>cf%records(crec)
    or=>of%records(orec)
    
    allocate(or%pds(200))
    or%pds=0

    do i=1,size(cr%pds)
       or%pds(i)=cr%pds(i)
    enddo

    or%pds(3)=igrid

!    print 10,of%filename,orec,or%pds
10  format(A,' record #',I0,' PDS=',200(I0,', '))
  end subroutine copy_pds

  ! --------------------------------------------------------------------

  subroutine make_space_for_record(gf,n)
    implicit none
    type(gribfile), pointer :: gf
    type(gribrecord), pointer :: newrecords(:)
    logical, pointer :: newdone(:)
    integer, intent(in) :: n
    integer :: i
    if(.not.associated(gf%records)) then
       allocate(gf%records(max(n+9,10)))
       allocate(gf%done(max(n+9,10)))
       gf%done=.false.
    endif
    if(n>size(gf%records)) then
       allocate(newrecords(max(n+9,2*size(gf%records)+1)))
       allocate(newdone(max(n+9,2*size(gf%records)+1)))
       newdone=.false.
       do i=1,size(gf%records)
          newdone(i)=gf%done(i)
          newrecords(i)%mask=>gf%records(i)%mask
          newrecords(i)%data=>gf%records(i)%data
          newrecords(i)%cgds=>gf%records(i)%cgds
          newrecords(i)%cpds=>gf%records(i)%cpds
          newrecords(i)%gds=>gf%records(i)%gds
          newrecords(i)%gdsf=>gf%records(i)%gdsf
          newrecords(i)%pds=>gf%records(i)%pds
          newrecords(i)%secloc=gf%records(i)%secloc
          newrecords(i)%reclen=gf%records(i)%reclen
          newrecords(i)%npoints=gf%records(i)%npoints
          newrecords(i)%ilengds=gf%records(i)%ilengds
          newrecords(i)%ilengdsf=gf%records(i)%ilengdsf
       end do
       deallocate(gf%records)
       deallocate(gf%done)
       gf%records=>newrecords
       gf%done=>newdone
    end if
  end subroutine make_space_for_record

  ! --------------------------------------------------------------------

  subroutine write_record(fp,irec,igrid,icomp,istat)
    implicit none

    ! interface
    !   SUBROUTINE W3FI72(ITYPE,FLD,IFLD,IBITL,&
    !  &                  IPFLAG,ID,PDS,&
    !  &                  IGFLAG,IGRID,IGDS,ICOMP,&
    !  &                  IBFLAG,IBMAP,IBLEN,IBDSFL,&
    !  &                  NPTS,KBUF,ITOT,JERR)

    !     integer :: itype,ibitl,ipflag,igflag,igrid,igds,,icomp,ibflag,
    !     REAL            FLD(*)
    !   INTEGER         IFLD(*)
    !   INTEGER         IBMAP(*)      INTEGER         ID(*)
    !   CHARACTER * 1   PDS(*)
    ! end interface

    type(gribfile), pointer :: fp
    integer, intent(in) :: irec
    integer, intent(inout) :: istat
    integer, intent(in) :: igrid, icomp

    character(len=:), pointer :: gribdat
    integer, pointer :: ibitmask(:)
    real, pointer :: adata(:)
    type(gribrecord),pointer :: rc
    character(len=4) :: recsize
    integer :: kfo

    real :: fmin,fmax
    integer :: ipds(200),igds(200),nbit,itot,kbm,i,kf
    integer :: kpds(200),kgds(200),ibds(200)
    character :: pds(400)
    real :: sum, avg

    nullify(gribdat,ibitmask,adata,rc)

    rc=>fp%records(irec)

    kf=rc%npoints
    kpds=0
    kgds=0
    kpds(1:size(rc%pds))=rc%pds
    kgds(1:size(rc%gds))=rc%gds

    kpds(4)=192

    ! Convert KPDS to IPDS format for input to w3fi72
    call R63W72(kpds,kgds,ipds,igds)
    ibds=0

    ! Count valid datapoints
    kbm=kf
    if(ipds(7)/=0) then
       allocate(ibitmask(kf))
       kbm=0
       !$OMP PARALLEL DO PRIVATE(I) REDUCTION(+:KBM)
       do i=1,kf
          if(rc%mask(i)) then
             ibitmask(i)=1
             kbm=kbm+1
             sum=sum+rc%data(i)
          else
             ibitmask(i)=0
          endif
       enddo
       if(kbm==rc%npoints) then
          !write(0,*) 'NO BITMASK - ALL POINTS ON'
          ipds(7)=0
       else
          avg=sum/kbm
          !$OMP PARALLEL DO PRIVATE(I) REDUCTION(+:KBM)
          do i=1,kf
             if(.not.rc%mask(i)) then
                rc%data(i)=avg
             endif
          enddo
       endif
    else
       !write(0,*) 'NO BITMASK - IPDS(7)==0'
    endif

    ! Get number of bits and round data.
    allocate(adata(kf))
    if(kbm==0) then
       ! No valid data.
       !$OMP PARALLEL DO PRIVATE(I)
       do i=1,kf
          adata(i)=0
       enddo
    else
       call GETBIT(ipds(7),0,ipds(25),kf,ibitmask,rc%data,adata,fmin,fmax,nbit)
       nbit=min(nbit,maxbit)
    endif

    allocate(character(len=size(rc%mask)+size(rc%data)*5+1000) :: gribdat)

    call w3fi72(0,adata,0,nbit,0,ipds,pds, &
                1,255,igds,0,0,ibitmask,kf,ibds,&
                kfo,gribdat,itot,istat)

    if(istat/=0) then
       write(0,10) fp%filename,irec,'w3fi72 failed',istat
10     format(A,' (',I0,'): ',A,' status=',I0)
       return
    endif

20  format(A,' (',I0,'): write ',I0,' bytes plus head/tail with ',I0,' in mask')
    !print 20, fp%filename,irec,itot,kbm

    recsize(1:1)=char(iand(ishft(itot,-24),255))
    recsize(2:2)=char(iand(ishft(itot,-16),255))
    recsize(3:3)=char(iand(ishft(itot, -8),255))
    recsize(4:4)=char(iand(      itot,     255))

    write(fp%unit) recsize,gribdat(1:itot),recsize
    fp%writepos=fp%writepos+int(itot,kind=8)+8_8
    fp%nrecords=irec
    istat=0

    if(associated(gribdat)) deallocate(gribdat)
    if(associated(adata)) deallocate(adata)
    if(associated(ibitmask)) deallocate(ibitmask)
  end subroutine write_record

  ! --------------------------------------------------------------------
  
  subroutine read_record(fp,irec,istat)
    type(gribfile), pointer :: fp
    type(gribrecord),pointer :: rc
    integer, intent(in) :: irec
    integer, intent(inout) :: istat
    character(len=:), pointer :: gribdat
    integer :: kptr(200),kret
    integer :: kgds(200), kpds(200), flags, npoints
    real, pointer :: data(:)
    logical :: gdsflag, bmsflag
    integer(kind=8) :: bms_start, unused, bitmask_size,bms_len,ones, avail,i
    real :: dmin,dmax
    nullify(gribdat,data)

    kpds=0
    kgds=0

    rc=>fp%records(irec)
    istat=1

    allocate(character(len=rc%reclen) :: gribdat)
    read(fp%unit,pos=rc%secloc(0)) gribdat

    flags=ichar(gribdat(16:16))
    gdsflag=iand(flags,128)/=0
    bmsflag=iand(flags,64)/=0

    if(bmsflag) then
       bms_start=rc%secloc(3)-rc%secloc(0)+1
       bms_len=rc%secloc(4)-rc%secloc(3)
       unused=ichar(gribdat(bms_start+3:bms_start+3))
       avail=char2toint8(gribdat,bms_start+4)
       if(avail/=0) then
          write(0,20) fp%filename,'predefined bitmasks are unsupported'
          goto 1000
       endif

       bitmask_size=(bms_len-6)*8-unused
    endif

    if(.not.associated(bitmask_storage)) then
       allocate(bitmask_storage(1048576*24*8)) ! theoretical maximum gridpoint count
    endif

    kgds=0
    kptr=0
    CALL FI631(gribdat,KPTR,KPDS,KRET)
    if(kret/=0) then
       write(0,20) fp%filename,'FI631 failed (cannot find GRIB1 record start)'
       goto 1000
    endif

    CALL FI632(gribdat,KPTR,KPDS,KRET)
    if(kret/=0) then
       write(0,20) fp%filename,'FI632 failed (cannot split record into sections)'
       goto 1000
    endif

    IF (IAND(KPDS(4),128).NE.0) THEN
       CALL FI633(gribdat,KPTR,KGDS,KRET)
       if(kret/=0) then
          write(0,20) fp%filename,'FI633 failed (cannot get grid information)'
          goto 1000
       endif
    ENDIF

    CALL FI634(gribdat,KPTR,KPDS,KGDS,bitmask_storage,KRET)
    if(kret/=0 .and. kret/=9) then
       write(0,20) fp%filename,'FI634 failed (mask processing)'
       goto 1000
    endif
    npoints=kptr(10)
    rc%npoints=npoints

    allocate(rc%mask(npoints))
    rc%mask(1:npoints)=bitmask_storage(1:npoints)

    allocate(rc%data(npoints))

    kptr=0
    call W3FI63(gribdat,kpds,kgds,bitmask_storage,rc%data,kptr,kret)
    if(kret/=0) then
       write(0,20) fp%filename,'w3fi63 failed (cannot decode grib message)'
       goto 1000
    endif
    
    if(kptr(10)/=rc%npoints) then
       write(0,20) fp%filename,'w3fi63 produced different point count than FI634'
       goto 1000
    endif


    dmin=9e19
    dmax=-9e19
    if(kptr(12)/=0) then
       if(.not. bmsflag) then
          write(0,20) fp%filename,'ERROR: FI634 wanted a bitmask but section 1 said none was needed'
          goto 1000
       endif
       ones=0
       do i=1,rc%npoints
          if(rc%mask(i)) then
             ones=ones+1
             dmin=min(dmin,rc%data(i))
             dmax=max(dmax,rc%data(i))
          endif
       enddo
       !print 12, fp%filename,kptr(10),ones,minval(rc%data),maxval(rc%data)
    else
       do i=1,rc%npoints
          dmin=min(dmin,rc%data(i))
          dmax=max(dmax,rc%data(i))
       enddo
       !print 10, fp%filename,kptr(10),minval(rc%data),maxval(rc%data)
    end if

    if(associated(gribdat)) deallocate(gribdat)

    allocate(rc%gds(200))
    rc%gds=kgds
    allocate(rc%pds(200))
    rc%pds=kpds

    istat=0
    return

1000 continue ! error handling
    if(associated(gribdat)) deallocate(gribdat)
    if(associated(data)) deallocate(data)

20  format(A,': ERROR: ',A)
10  format(A,': npoints=',I0,' (no bitmask) min=',F0.3,' max=',F0.3)
12  format(A,': npoints=',I0,' bitmask on=',I0,' min=',F0.3,' max=',F0.3)
  end subroutine read_record

  ! --------------------------------------------------------------------

  subroutine make_record(fp,irec,ogds,istat)
    !! Allocate space and make PDS for target record.
    implicit none

    interface
       INTEGER FUNCTION LENGDSF(KGDS,KGDSF)
         INTEGER KGDS(200),KGDSF(200)
       END FUNCTION LENGDSF
    end interface

    interface
       integer function lengds(gds)
         integer, intent(in) :: gds(200)
       end function lengds
    end interface

    type(gribfile), pointer :: fp !! Pointer to output file
    integer, intent(in) :: ogds(200) !! Output GDS
    integer, intent(in) :: irec !! record index within fp
    integer, intent(inout) :: istat !! error status (0 on success)

    integer :: npoints

    type(gribrecord),pointer :: rc
    istat=1
    rc=>fp%records(irec)

    allocate(rc%gds(max(size(ogds),200)))
    allocate(rc%gdsf(200))
    rc%gds=0
    rc%gds(1:size(ogds))=ogds
    rc%ilengds=lengds(rc%gds)
    
    rc%ilengdsf=lengdsf(rc%gds(1:200),rc%gdsf)
    npoints=max(rc%ilengdsf,rc%ilengds)
    rc%npoints=npoints

    allocate(rc%mask(npoints+1000))
    allocate(rc%data(npoints+1000))

    rc%mask=.false.
    rc%data=0

    !print 10, fp%filename,npoints

    istat=0
10  format(A,': allocated record with ',I0,' points for output')
  end subroutine make_record

  ! --------------------------------------------------------------------

  subroutine free_record(fp,irec,istat)
    type(gribfile), pointer :: fp
    type(gribrecord),pointer :: rc
    integer, intent(in) :: irec
    integer, intent(inout) :: istat
    istat=0
    rc=>fp%records(irec)
    if(associated(rc%mask)) then
       deallocate(rc%mask)
       nullify(rc%mask)
    endif
    if(associated(rc%data)) then
       deallocate(rc%data)
       nullify(rc%data)
    endif
    if(associated(rc%gds)) then
       deallocate(rc%gds)
       nullify(rc%gds)
    endif
    if(associated(rc%cgds)) then
       deallocate(rc%cgds)
       nullify(rc%cgds)
    endif
    if(associated(rc%gdsf)) then
       deallocate(rc%gdsf)
       nullify(rc%gdsf)
    endif
    if(associated(rc%pds)) then
       deallocate(rc%pds)
       nullify(rc%pds)
    endif
    if(associated(rc%cpds)) then
       deallocate(rc%cpds)
       nullify(rc%cpds)
    endif
  end subroutine free_record

  ! --------------------------------------------------------------------

  subroutine interp_scalar_record(inf,inrec,ouf,ourec,mergeflag,istat)
    type(gribfile), pointer :: inf,ouf
    integer,intent(in) :: inrec,ourec
    logical,intent(in) :: mergeflag
    integer,intent(inout) :: istat

    call interp_record_impl(inf,inrec,ouf,ourec,mergeflag,istat)
  end subroutine interp_scalar_record

  ! --------------------------------------------------------------------

  subroutine interp_velocity_record(inf,inrec,ouf,ourec,mergeflag,istat,invec,ouvec)
    type(gribfile), pointer :: inf,ouf
    integer,intent(in) :: inrec,ourec
    integer,intent(in) :: invec,ouvec
    logical,intent(in) :: mergeflag
    integer,intent(inout) :: istat

    integer :: inrec1,inrec2
    integer :: ourec1,ourec2

    ! Make sure velocity components are in the right order for iplib.

    if(ichar(inf%records(invec)%cpds(9:9))==33) then
       inrec2=inrec
       inrec1=invec
       ourec2=ourec
       ourec1=ouvec
    else
       inrec1=inrec
       inrec2=invec
       ourec1=ourec
       ourec2=ouvec
    endif

    call interp_record_impl(inf,inrec1,ouf,ourec1,mergeflag,istat,inrec2,ourec2)

  end subroutine interp_velocity_record

  ! --------------------------------------------------------------------

  subroutine interp_record_impl(inf,inrec,ouf,ourec,mergeflag,istat,invec,ouvec)
    implicit none

    interface
       SUBROUTINE INTGRIB1(K1F,KGDS1F,K2F,KGDS2F,MRL,MRO,&
            IV,IP,IPOPT,KGDS1,K1,IB1,L1,F1,G1,KGDS2,K2, &
            IB2,L2,F2,G2,IRET) 
         IMPLICIT NONE
         INTEGER IPOPT(20),K1F,K2F,MRL,MRO,IV,IP,K1,IB1,K2,IB2,IRET
         INTEGER KGDS1F(200),KGDS2F(200)
         LOGICAL*1 L1(K1),L2(K2)
         REAL F1(K1),F2(K2),G1(K1),G2(K2)
         INTEGER KGDS1(200),KGDS2(200)
       END SUBROUTINE INTGRIB1
    end interface

    type(gribfile), pointer :: inf,ouf
    integer,intent(in) :: inrec,ourec
    integer,optional,intent(in) :: invec,ouvec
    logical,intent(in) :: mergeflag
    integer,intent(inout) :: istat

    real :: mx,mn

    type(gribrecord), pointer :: ir,or,iv,ov
    real, pointer :: data2(:), data3(:)
    logical(kind=1), pointer :: mask2(:)
    integer :: i,no,ni,im
    integer :: ivflag, ipopt(20),mro,ip,ibmp,obmp,k1f,k2f,mrl

    nullify(ir,or,iv,ov,data2,data3,mask2)

    istat=1

    ip=0 ! interpolation type (0=bilinear)
    ivflag=0 ! vector flag (0=scalar)
    mro=1 ! size of crot and srot (1 if iv=0)
    ipopt=0 ! interpolation options

    ibmp=1 ! we always make a bitmask
    obmp=1 ! we always make a bitmask

    ir=>inf%records(inrec)
    or=>ouf%records(ourec)

    if(present(invec)) then
       iv=>inf%records(invec)
       ivflag=1
    else
       iv=>ir
       ivflag=0
    endif
    if(present(ouvec)) then
       ov=>ouf%records(ouvec)
       ivflag=1
    else
       ov=>or
       ivflag=0
    endif

    if(present(invec) .neqv. present(ouvec)) then
       write(0,*) "INTERNAL ERROR: PRESENT(INVEC)=",present(invec),&
                  ' PRESENT(OUVEC)=',present(ouvec)
       stop 23
    endif

    ni=ir%npoints
    no=or%npoints

    if(mergeflag) then
       allocate(data2(or%npoints+1000))
       if(present(ouvec)) then
          allocate(data3(ov%npoints+1000))
       else
          data3=>data2
       endif
       allocate(mask2(or%npoints+1000))
    else
       data2=>or%data
       data3=>ov%data
       mask2=>or%mask
    endif

    ! Two loops to ensure initialization is in the write thread for
    ! each memory location.

    !$OMP PARALLEL DO PRIVATE(I)
    do i=1,no
       mask2(i)=.false.
    enddo
    !$OMP PARALLEL DO PRIVATE(I)
    do i=1,1000
       mask2(no+i)=.false.
    enddo

    k1f=ir%ilengdsf
    if(ir%ilengdsf==ir%ilengds) k1f=1
    k2f=or%ilengdsf
    if(or%ilengdsf==or%ilengds) k2f=1
    mrl=max(or%ilengds,or%ilengdsf)

    if(ivflag/=0) then
       mro=mrl
    endif

    if(.not.associated(iv%data)) then
       write(0,*) 'no iv%data'
       stop 33
    endif
    if(.not.associated(data3)) then
       write(0,*) 'no data3'
       stop 36
    endif

    istat=999    
    call INTGRIB1(k1f,ir%gdsf,k2f,or%gdsf,&
         mrl,mro,ivflag,ip,ipopt,&
         ir%gds,ni,ibmp,ir%mask,ir%data,iv%data,&
         or%gds,no,obmp,mask2,  data2,  data3 , &
         istat)
    if(istat/=0) return

    if(mergeflag) then
       !$OMP PARALLEL DO PRIVATE(I)
       do i=1,no
          if(mask2(i)) then
             or%data(i)=data2(i)
             or%mask(i)=.true.
          endif
       enddo
       if(present(ouvec)) then
          !$OMP PARALLEL DO PRIVATE(I)
          do i=1,no
             if(mask2(i)) then
                ov%data(i)=data3(i)
                ov%mask(i)=.true.
             endif
          enddo
          deallocate(data3)
       endif
       deallocate(data2)
       deallocate(mask2)
    else
       !$OMP PARALLEL DO PRIVATE(I)
       do i=1,no
          if(mask2(i)) then
             ov%mask(i)=.true.
          endif
       enddo

    endif

    mx=-9e19
    mn=9e19
    im=0
    if(present(ouvec)) then
       !$OMP PARALLEL DO PRIVATE(I) REDUCTION(max:mx) REDUCTION(min:mn) REDUCTION(+:im)
       do i=1,no
          if(or%mask(i) .and. ov%mask(i)) then
             im=im+1
             mx=max(mx,sqrt(or%data(i)*or%data(i) + ov%data(i)*ov%data(i)))
             mn=min(mn,sqrt(or%data(i)*or%data(i) + ov%data(i)*ov%data(i)))
          endif
       enddo
    else
       !$OMP PARALLEL DO PRIVATE(I) REDUCTION(max:mx) REDUCTION(min:mn) REDUCTION(+:im)
       do i=1,no
          if(or%mask(i)) then
             im=im+1
             mx=max(mx,or%data(i))
             mn=min(mn,or%data(i))
          endif
       enddo
    endif

    if(present(ouvec)) then
       if(mergeflag) then
          !print 21,ouf%filename,ourec,ouvec,inf%filename,inrec,invec,mn,mx,no,im
       else
          !print 11,ouf%filename,ourec,ouvec,inf%filename,inrec,invec,mn,mx,no,im
       endif
    else
       if(mergeflag) then
          !print 20,ouf%filename,ourec,inf%filename,inrec,mn,mx,no,im
       else
          !print 10,ouf%filename,ourec,inf%filename,inrec,mn,mx,no,im
       endif
    endif
10 format(A,':',I0,': filled with ',A,':',I0,' min=',F0.3,' max=',F0.3,' npoints=',I0,' mask=',I0)
20 format(A,':',I0,': merged with ',A,':',I0,' min=',F0.3,' max=',F0.3,' npoints=',I0,' mask=',I0)

11 format(A,':',I0,'&',I0,': filled with ',A,':',I0,'&',I0,' min=',F0.3,' max=',F0.3,' npoints=',I0,' mask=',I0)
21 format(A,':',I0,'&',I0,': merged with ',A,':',I0,'&',I0,' min=',F0.3,' max=',F0.3,' npoints=',I0,' mask=',I0)
  end subroutine interp_record_impl

  ! --------------------------------------------------------------------

  subroutine interp_files(rm,istat)
    implicit none
    type(regrid_merge), intent(inout) :: rm
    integer,intent(inout) :: istat

    integer :: cfile,sfile,ofile, crec,srec,orec,nfiles
    integer :: cvec,svec,ovec, next_orec
    type(gribfile), pointer :: cf,sf,of
    type(gribrecord), pointer :: cr,sr,or, cv,sv,ov

    nullify(cf,sf,of,cr,sr,or,cv,sv,ov)

    next_orec=1
    nfiles=size(rm%files)-1
    of=>rm%files(nfiles+1)
    controlfile: do cfile=1,nfiles
       cf=>rm%files(cfile)
       controlrec: do crec=1,cf%nrecords
          nullify(cv,ov)
          if(cf%done(crec)) cycle controlrec
          cr=>cf%records(crec)
          if(cr%reclen<16) then
             cf%done(crec)=.true.
             cycle controlrec
          endif
          !print 10, cf%filename,crec,cf%nrecords
          cf%done(crec)=.true.

          orec=next_orec
          next_orec=next_orec+1

          !print 20, of%filename,orec

          call read_record(cf,crec,istat)
          if(istat/=0) cycle controlrec

          if(ichar(cr%cpds(9:9))==33 .or. ichar(cr%cpds(9:9))==34) then
             ! Find other velocity component in control file
             nullify(cv)
             find_cvec: do cvec=crec+1,cf%nrecords
                if(cvec==crec) cycle find_cvec
                if(cf%done(cvec)) cycle find_cvec
                cv=>cf%records(cvec)
                if(cv%reclen<16) then
                   cf%done(cvec)=.true.
                   cycle find_cvec
                endif

                if(.not.other_vector_component(cr%cpds,cv%cpds)) cycle find_cvec
                exit find_cvec
             enddo find_cvec
             if(.not.associated(cv)) then
                !print 70, cf%filename,crec
                cycle controlrec
             endif
             cf%done(cvec)=.true.
             !print 11, cf%filename,cvec,cf%nrecords
             call read_record(cf,cvec,istat)
             if(istat/=0) cycle controlrec
          else
             cvec=-1
             nullify(cv)
          endif

          call make_space_for_record(of,orec)
          if(istat/=0) cycle controlrec
          call make_record(of,orec,rm%gds,istat)
          if(istat/=0) cycle controlrec
          call copy_pds(cf,crec,of,orec,rm%igrid)

          if(associated(cv)) then
             ! VECTOR interpolation of first file
             ovec=next_orec
             next_orec=next_orec+1
             !print 21,of%filename,ovec
             call make_space_for_record(of,ovec)
             ov=>of%records(ovec)
             if(istat/=0) cycle controlrec
             call make_record(of,ovec,rm%gds,istat)
             if(istat/=0) cycle controlrec
             call copy_pds(cf,cvec,of,ovec,rm%igrid)
             call interp_velocity_record(cf,crec,of,orec,.false.,istat,cvec,ovec)
             if(istat/=0) cycle controlrec
          else
             ! SCALAR interpolation of first file
             nullify(ov)
             call interp_scalar_record(cf,crec,of,orec,.false.,istat)
             if(istat/=0) cycle controlrec
          endif

          secondfile: do sfile=cfile+1,nfiles
             nullify(sv)
             sf=>rm%files(sfile)
             interprec: do srec=1,sf%nrecords
                if(sf%done(srec)) then
                   !print 50, sf%filename,srec
                   cycle interprec
                endif
                sr=>sf%records(srec)
                if(sr%reclen<16) then
                   write(0,40) sf%filename,srec,sr%reclen
                   sf%done(srec)=.true.
                   cycle controlrec
                endif

                if(srec/=crec) then
                   write(0,60) cf%filename,crec,sf%filename,srec
                endif

                if(.not.associated(sr%cpds)) then
                   write(0,80) sf%filename,srec
                   stop 9
                endif
                if(.not.associated(cr%cpds)) then
                   write(0,80) cf%filename,crec
                   stop 10
                endif
                if(.not.same_pds(sr%cpds,cr%cpds)) cycle interprec

                sf%done(srec)=.true.
                !print 15, sf%filename,srec,sf%nrecords


                if(associated(cv)) then
                   ! Find other velocity component in control file
                   nullify(sv)
                   find_svec: do svec=srec+1,sf%nrecords
                      if(svec==srec) cycle find_svec
                      if(sf%done(svec)) cycle find_svec
                      sv=>sf%records(svec)
                      if(sv%reclen<16) then
                         sf%done(svec)=.true.
                         cycle find_svec
                      endif

                      if(.not.other_vector_component(sr%cpds,sv%cpds)) then
                         cycle find_svec
                      endif
                      exit find_svec
                   enddo find_svec
                   if(.not.associated(sv)) then
                      !print 70, sf%filename,srec
                      cycle interprec
                   endif
                   sf%done(svec)=.true.
                   !print 16, sf%filename,svec,sf%nrecords
                   call read_record(sf,svec,istat)
                   if(istat/=0) cycle interprec
                endif


                call read_record(sf,srec,istat)
                if(istat/=0) cycle interprec

                if(associated(sv)) then
                   ! VECTOR interpolation of second file
                   call interp_velocity_record(sf,srec,of,orec,.true.,istat,svec,ovec)
                   if(istat/=0) cycle controlrec
                else
                   call interp_scalar_record(sf,srec,of,orec,.true.,istat)
                endif
                if(istat/=0) then
                   write(0,*) 'interpolation failed with status ',istat
                   stop 9
                   cycle interprec
                endif

                call free_record(sf,srec,istat)
                if(associated(sv)) call free_record(sf,svec,istat)
                cycle secondfile
             end do interprec
          enddo secondfile

          call write_record(of,orec,rm%igrid,rm%gridwind,istat)
          if(associated(ov)) then
             call write_record(of,ovec,rm%igrid,rm%gridwind,istat)
          endif

          call free_record(cf,crec,istat)
          call free_record(of,orec,istat)
          if(associated(cv)) call free_record(cf,cvec,istat)
          if(associated(ov)) call free_record(of,ovec,istat)
       end do controlrec
    enddo controlfile
10  format(A,': input #1 record #',I0,' of ',I0)
11  format(A,': input #1 record #',I0,' of ',I0,' (velocity second component)')
15  format(A,': input #2 record #',I0,' of ',I0)
16  format(A,': input #2 record #',I0,' of ',I0,' (velocity second component)')
20  format(A,': target record #',I0)
21  format(A,': target record #',I0,' (velocity second component)')
40  format(A,': record #',I0,' has invalid size ',I0)
50  format(A,': record #',I0,' is already done')
60  format('WARNING: ',A,' using record #',I0,' but using ',A,' record ',I0)
70  format(A,': record #',I0,' is a vector component with no partner')
80  format(A,': record ',I0,' has no CPDS')
  end subroutine interp_files

  ! --------------------------------------------------------------------

  subroutine index_files(rm,istat)
    implicit none
    type(regrid_merge), intent(inout) :: rm
    integer,intent(inout) :: istat

    integer :: ifile,nfiles,unit
    integer(kind=8), parameter :: bufsize=128
    character(len=:), pointer :: buf,filename
    integer(kind=8) :: idx, gpos, buflen, pdslen, bmslen
    integer(kind=8) :: recsize, gribver, filesize, left, gdslen
    integer :: irecord

    logical :: gdsflag, bmsflag
    integer :: flags

    type(gribfile), pointer :: gf
    type(gribrecord), pointer :: gr

    istat=0

    allocate(character(len=bufsize) :: buf)
    
    nfiles=size(rm%files)-1
    fileloop: do ifile=1,nfiles
       gf=>rm%files(ifile)
       unit=gf%unit
       filename=>gf%filename
       inquire(unit=unit,size=filesize)

       !print '(A,": unit ",I0)',filename,unit

       idx=1
       irecord=0
       readloop: do while(idx<filesize)
          left=filesize-idx+1
          buflen=min(bufsize,left)
          read(unit,pos=idx,iostat=istat) buf(1:buflen)
          if(istat/=0) exit readloop
          gpos=index(buf(1:buflen),'GRIB',kind=8)
          if(gpos<0) then
             idx=idx+buflen
             cycle readloop
          else
             irecord=irecord+1
          endif

          gpos=idx+gpos-1
          read(unit,pos=gpos,iostat=istat) buf(1:16)
          if(istat/=0) goto 2000
          gribver=ichar(buf(8:8))
          if(gribver==1) then
             recsize=char3toint8(buf,5)
          elseif(gribver==2) then
             recsize=transfer(buf(9:16),recsize)
             if(big_endian) recsize=byteswap8(recsize)
             goto 1000
          else
             WRITE(0,40) filename,gribver
             istat=3
             return
          endif
          if(recsize<22) goto 1000

          flags=ichar(buf(16:16))
          gdsflag=iand(flags,128)/=0
          bmsflag=iand(flags,64)/=0

          call make_space_for_record(gf,irecord)
          gr=>gf%records(irecord)
          gr%secloc(0)=gpos
          gr%reclen=recsize

          ! Read PDS
          read(unit,pos=gpos+8) buf(1:3)
          pdslen=char3toint8(buf,1)
          allocate(character(len=pdslen) :: gr%cpds)
          read(unit,pos=gpos+8) gr%cpds
          gr%secloc(1)=gpos+8

          if(gdsflag) then
             ! Read GDS
             read(unit,pos=gpos+8+pdslen) buf(1:3)
             gdslen=char3toint8(buf,1)
             allocate(character(len=pdslen) :: gr%cgds)
             read(unit,pos=gpos+8+pdslen) gr%cgds
             gr%secloc(2)=gpos+8+pdslen
          else
             gdslen=0
             gr%secloc(2)=gpos+8+pdslen
          endif

          if(bmsflag) then
             gr%secloc(3)=gpos+8+pdslen+gdslen
             read(unit,pos=gr%secloc(3)) buf(1:3)
             bmslen=char3toint8(buf,1)
             gr%secloc(4)=gr%secloc(3)+bmslen
          else
             bmslen=0
             gr%secloc(3)=gpos+8+pdslen+gdslen
             gr%secloc(4)=gpos+8+pdslen+gdslen+bmslen
          endif

          ! Read end of record
          read(unit,pos=gpos+recsize-4) buf(1:4)
          if(buf(1:4)/='7777') then
             write(0,60) filename,buf(1:4)
          endif

          idx=gpos+recsize-1
          cycle readloop
!          print 30,filename,recsize,pdslen,gdslen

1000      continue ! jump here when we discard the record
          irecord=irecord-1
          idx=gpos+recsize-1 ! seek to end of GRIB record
          cycle readloop

2000      continue ! jump here to skip the rest of the file
          irecord=irecord-1
          exit readloop
       enddo readloop
       gf%nrecords=irecord       
       do irecord=1,gf%nrecords
          if(.not.associated(gf%records(irecord)%cpds)) then
             write(0,*) gf%filename,irecord
             stop 2
          endif
       enddo
       write(0,50) gf%filename,gf%nrecords
    enddo fileloop

80  format(A,': record ',I0,' has no CPDS')
70  format(A,': error in FI633 unpacking GDS')
60  format(A,': ERROR: record end string is wrong: "',A,'" instead of "7777"')
50  format(A,': found ',I0,' records')    
40  format(A,': unrecognized GRIB version ',I0,'.  Abort.')
30  format(A,': record points=',I0,' bytes=',I0,' pds size=',I0,' gds size=',I0)
  end subroutine index_files

  ! --------------------------------------------------------------------

  subroutine add_infile(rm,filename,istat)
    implicit none
    character*(*), intent(in) :: filename
    type(regrid_merge), intent(inout) :: rm
    integer,intent(inout) :: istat

    integer :: ifile

    ifile=rm%nopened+1
    if(ifile==size(rm%files)) then
       write(0,19) filename,size(rm%files)
       istat=3
       return
    endif
    
    open(newunit=rm%files(ifile)%unit,access='stream',status='old',&
         file=filename,iostat=istat)
    allocate(character(len=len(filename)) :: rm%files(ifile)%filename)
    rm%files(ifile)%filename=filename
    if(istat==0) then
       !print 10,filename,ifile
       rm%nopened=ifile
    else
       write(0,11) filename
       stop 18
    endif

10  format(A,': opened file #',I0,' for read')
11  format(A,': cannot open for read')
19  format(A,': too many input files (nfiles=',I0,' includes output)')
  end subroutine add_infile

  ! --------------------------------------------------------------------

  subroutine add_outfile(rm,filename,istat)
    character*(*), intent(in) :: filename
    type(regrid_merge), intent(inout) :: rm
    integer,intent(inout) :: istat

    integer :: ifile

    ifile=rm%nopened+1
    if(ifile/=size(rm%files)) then
       write(0, 19) filename,ifile,size(rm%files)
       istat=1
       return
    endif

    open(newunit=rm%files(ifile)%unit,access='stream',status='replace',&
         file=filename,iostat=istat)
    allocate(character(len=len(filename)) :: rm%files(ifile)%filename)
    rm%files(ifile)%filename=filename
    
    if(istat==0) then
       rm%nopened=ifile
       write(0,10) filename,ifile
    else
       write(0,11) filename
       stop 32
    endif

10  format(A,': opened file #',I0,' for write')
11  format(A,': cannot open for write')
19  format(A,': too few input files (received ',I0,' of ',I0,')')
  end subroutine add_outfile

  ! --------------------------------------------------------------------

  subroutine init_regrid_merge(rm,carg,nfiles,istat)
    implicit none
    character*(*), intent(in) :: carg
    type(regrid_merge), intent(inout) :: rm
    integer,intent(out) :: istat
    integer,intent(in) :: nfiles
    character(len=400) :: gds
    integer kgdsi(200),karg(100),igi,iret,lgds,mi
    integer, external :: lengds
    igi=-1
    kgdsi=0
    kgdsi(20)=255
    karg(1)=-1
    karg(2:100)=kgdsi(1:99)
    call FPARSEI(carg(1:len_trim(carg)),100,karg)
    igi=karg(1)
    if(igi>0 .and. igi<255) then
       call makgds(igi,kgdsi,gds,lgds,iret)
       if(iret/=0) igi=-1
       rm%igrid=igi
    else
       rm%igrid=255
       kgdsi(1:99)=karg(2:100)
    endif
    IF(IGI.LT.0.OR.IGI.EQ.0.OR.IGI.GT.255) THEN
       write(0,33) 'invalid output grid',carg
       istat=1
       return
    ENDIF
    mi=lengds(kgdsi)
    if(mi<=0) then
       write(0,33) 'unsupported output grid',carg
       istat=2
       return
    end if
    if(iand(kgdsi(6),8)==0) then
       rm%gridwind=0
    else
       rm%gridwind=1
    endif

    allocate(rm%files(nfiles))
    allocate(rm%gds(200))
    rm%gds=kgdsi

    istat=0

    !print 44,rm%gds

33  format('regrid_merge: ',A,' ',A)
44  format('KGDS = ',200(I0,', '))
  end subroutine init_regrid_merge

end module module_regrid_merge
