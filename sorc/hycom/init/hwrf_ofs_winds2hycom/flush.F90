#if defined(AIX)
subroutine flush(iunit)
implicit none
integer iunit
!
! --- wrapper for flush system call under AIX.
!
integer*4 iunit4
!
iunit4=iunit
call flush_(iunit4)
return
end subroutine flush
#endif /* AIX */
#if defined(IFC)
subroutine flush(iunit)
implicit none
integer iunit
integer*4 iunit4
!
! --- disable the flush system call under Intel's IFC compiler.
!
iunit4=iunit
return
end subroutine flush
#endif /* IFC */
