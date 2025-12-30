program test_ud
integer :: idim
real,dimension(2,3,4)::a
integer :: i,j,k

idim=1

do i=1,2
do j=1,3
do k=1,4
   a(i,j,k)=j

enddo
enddo
enddo

call ud(a,2,3,4)
do i=1,2
do j=1,3
do k=1,4
write(*,*)i,j,k, a(i,j,k)
enddo
enddo
enddo
stop
end program test_ud

subroutine ud(a,n1,n,n2)
    integer idim
    real, dimension(n1,n,n2) :: a
    integer n,hn,i,n1,n2,j,k
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
