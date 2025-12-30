program atcf_prob
  use atcf_prob_module

  type(multimodel) :: mm
  integer :: ierr

  call mm%read_namelist('atcf_prob.nml','atcf_prob.nml.out',ierr)
  if(ierr/=0) then
100  format('atcf_prob.nml: error reading namelist (',I0,')')
     call error_handler()
  endif

  call mm%read_atcf(5000,ierr)
  if(ierr/=0) then
200  format('ERROR: no valid atcf files found for any models.')
     call error_handler()
  endif

  irec=0
  do while(mm%next_time(ftime))
     irec=irec+1
     call write_grads_dat(mm,'multimodel_prob.dat',irec)
  enddo

  call write_grads_ctl(mm,'multimodel_prob.ctl','multimodel_prob.dat',irec)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine error_handler()
1100 format(A)
    print 1100, 'EPIC FAIL.  Giving up.  Will produce no output.'
    print 1100, 'Check preceding messages for more details.'
    print 1100, 'Will now exit with status 1.'
    stop 1
  end subroutine error_handler

  subroutine write_grads_dat(mm,datfile,irec)
    class(multimodel) :: mm
    character*(*) :: datfile
    real(kind=4), pointer :: outdat(:,:)
    integer :: i,j

    allocate(outdat(mm%nx,mm%ny))

    !$OMP PARALLEL DO PRIVATE(i,j)
    do j=1,mm%ny
       do i=1,mm%nx
          if(.not.mm%prob%mask(i,j,1)) then
             outdat(i,j)=GRADS_BADVAL
          else
             outdat(i,j)=mm%prob%rdata(i,j,1)
          endif
       enddo
    enddo
    !$OMP END PARALLEL DO

    open(201,file=trim(datfile),access='DIRECT',recl=mm%nx*mm%ny)
    write(201,rec=mm%nx*mm%ny) outdat
    close(201)

    deallocate(outdat)
  end subroutine write_grads_dat

  subroutine write_grads_ctl(mm,ctlfile,dset,irec)
    use datetime_module, only: grads_date_len, grads_date
    character(len=grads_date_len) :: stratime
    integer, intent(in) :: irec
    class(multimodel) :: mm
    character*(*) :: dset, ctlfile

    stratime=' '
    call grads_date(stratime,mm%atime)

    open(200,file=trim(ctlfile),form='FORMATTED')
    write(200,10) dset
    write(200,20) 0,-999.
    write(200,30) mm%nx,mm%lons(1),mm%lons(2)-mm%lons(1)
    write(200,40) mm%ny,mm%lats(1),mm%lats(2)-mm%lats(1)
    write(200,50) ' '
    write(200,60) irec,trim(stratime),nint(mm%time_step)
    write(200,70) 1
    write(200,80) 'track_prob','track probability'
    write(200,90) ' '
    close(200)

10  format("DSET  ^",A)
20  format("TITLE Multi-Model Track Probability: Fcst Hour ",I0,/,"UNDEF ",F0.1)
30  format("XDEF ",I0," LINEAR ",F0.5,"  ",F0.5)
40  format("YDEF ",I0," LINEAR ",F0.5,"  ",F0.5)
50  format("ZDEF 1 LEVELS 1000",A)
60  format("TDEF ",I0," LINEAR ",A," ",I0,"hr")
70  format("VARS ",I0)
80  format(A10,      "       0 99 ",A)
90  format('ENDVARS',A)
  end subroutine write_grads_ctl
end program atcf_prob
