      program qc_river
!
!     collects from decoded river data the information needed for
!     selected rivers, and also assigns a qc value.
!
      implicit none
      CHARACTER*8, allocatable,  dimension(:) :: name,ident
      integer count,ic,nrivers,iriv,icc
      real,    allocatable, dimension(:) :: clat,clon,dchgd,rstgd
      real dchgd_iriv,rstgd_iriv
      integer, allocatable, dimension(:) :: ii, jj, qc
     
 50   format(i6)
 51   format(a8,1x,i4,1x,i4)
 60   format(a8,1x,f8.2,1x,f8.2,1x,f8.3,1x,f8.3)
 70   format(a8,1x,i4,1x,i4,1x,f8.3,1x,f8.3,1x,i1)

! 
! read data
!
      read(11,50) count 
      allocate(name(count),clat(count),clon(count),
     &     dchgd(count),rstgd(count))
      do ic = 1,count
         read(10,60)name(ic),clat(ic),clon(ic),dchgd(ic),rstgd(ic)
      enddo 

!
! read table of rivers
!           
      read(13,50) nrivers


      allocate(ident(nrivers),ii(nrivers),jj(nrivers),qc(nrivers))
            
      do iriv = 1,nrivers
         read(13,51) ident(iriv),ii(iriv),jj(iriv)
      enddo
      
     
            
!
! find the rivers in the data using the identifier
!
           
      do iriv = 1,nrivers
         qc(iriv) = 0
         dchgd_iriv = 0.0
         rstgd_iriv = 0.0
         icc = 0  
         do ic = 1,count
              
            if (ident(iriv) == name(ic)) then
               icc = icc +1
                  dchgd_iriv = dchgd(ic)+ dchgd_iriv
                  rstgd_iriv = rstgd(ic)+ rstgd_iriv
            endif
         enddo
         
         
         if (icc == 0) then
           dchgd_iriv = 0.0
            rstgd_iriv = 0.0
            qc(iriv)=0
          else
             dchgd_iriv = dchgd_iriv/float(icc)
            rstgd_iriv = rstgd_iriv/float(icc)
            qc(iriv)=1
          endif
!
! write only admissible stations
!
         write(60,70) ident(iriv),ii(iriv),jj(iriv),
     &        dchgd_iriv,rstgd_iriv,qc(iriv)
               
      enddo
                
  
      deallocate(ident,ii,jj,qc,name,clat,clon,dchgd,rstgd)

      stop
      end program qc_river
          
             
