program hwrf_regrid_merge
  use module_regrid_merge
  implicit none
  type(regrid_merge) :: rm
  integer :: nargs,iarg,arglen,istat
  character(len=:), pointer :: carg
  nargs=command_argument_count()
  argloop: do iarg=1,nargs
     call get_command_argument(iarg,length=arglen)
     allocate(character(len=arglen) :: carg)
     call get_command_argument(iarg,carg)
     if(iarg==1) then
        call init_regrid_merge(rm,carg,nargs-1,istat)
     elseif(iarg<nargs) then
        call add_infile(rm,carg,istat)
     else
        call add_outfile(rm,carg,istat)
     endif
     deallocate(carg)
  enddo argloop
  call index_files(rm,istat)
  call interp_files(rm,istat)

  if(istat/=0) stop 1
end program hwrf_regrid_merge
