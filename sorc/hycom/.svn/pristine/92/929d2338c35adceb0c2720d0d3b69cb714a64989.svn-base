      CHARACTER*8  SUBSET,name                                   
      Real*8 DATES(5),RTIM(6)
      Real*8 data1(2,500),data2(4,500),data3(5)
      real*8 xlocat(5)
      integer two,three,count,nocount
      DATA LUNIN /11/                                                   
      DATA BMISS /10E10/                                                
      data two/2/,three/3/
      equivalence(name,rpid)
                                                                        
c
c
C  OPEN AND READ THRU THE INPUT BUFR FILE                               
C  --------------------------------------                               
                                                                        
      CALL OPENBF(LUNIN,'IN',LUNIN)                                     
      IPRT = 0                                                          
                                                                        
C  READ THROUGH THE MESSAGES/SUBSETS IN THE FILE                        
C  ---------------------------------------------                        
      count = 0                                                                  
      nocount=0
      
c     write(6,*) 'Before the While loop',LUNIN,SUBSET,IDATE
      DO WHILE(IREADMG(LUNIN,SUBSET,IDATE).EQ.0)                       
c     write(6,*) 'Subset ',subset,'  idate ',idate
c     write(6,'('' finished ireaderm.'')')
       DO WHILE(IREADSB(LUNIN).EQ.0)                                    
c        write(6,*) ' Lunin  ',lunin
c        write(6,'('' finished ireaders.'')')
                                                                        
C  READ THE INTERNAL DATE AND CHECK FOR REALISM                         
C  --------------------------------------------                         
                                                                        
      CALL UFBINT(LUNIN,DATES,5,1,IRET,'YEAR MNTH DAYS HOUR MINU')
      IYR = NINT(DATES(1))
      MON = NINT(DATES(2))
      IDAY = NINT(DATES(3))                                           
      IHR = NINT(DATES(4))                                            
      MINS = NINT(DATES(5))                                            
c                                               
       call UFBINT(LUNIN,rpid,1,1,iret,'RPID')
c
      CALL UFBINT(LUNIN,XLOCAT,5,1,IRET,'CLAT CLON SELV')
      clat=xlocat(1)
      clon=xlocat(2)
      selv=xlocat(3)
c
c... Get Temp Data
c
      CALL UFBINT(LUNIN,data3,5,1,IRET,'TMDB SST1')
c
      if(data3(1).eq.bmiss)then
          tmdb = -99.99
       else
       tmdb = data3(1) -273.15
      endif
c
      if(data3(2).eq.bmiss)then
         sst = -99.99
        else
        sst = data3(2) -273.15
      endif
c
c
c  Get discharge data
c
      CALL UFBINT(LUNIN,data1,2,500,nlev,'DCHG DDDC')
c
      if(data1(1,1).eq.bmiss.or.data1(1,1).eq.00.0)then
          dchgd = -99.99
       else
        dchgd = data1(1,1)
      endif
c
      if(data1(2,1).eq.bmiss)then
         dchgi = -99.99
        else
        dchgi = data1(2,1)
      endif
c
c  Get River stage data
c
      CALL UFBint(LUNIN,data1,2,500,nlev1,'RSHM DDRS')
c
      if(data1(1,1).eq.BMISS.or.data1(1,1).eq.00.0)then
          rstgd = -99.99
        else
        rstgd = data1(1,1)
      endif
c
      if(data1(2,1).eq.bmiss)then
         rstgdi = -99.99
        else
        rstgdi = data1(2,1)
      endif
c
c Get Stage hgt elev above ngvd 1929, strm velocity, salinity
c
      CALL UFBINT(LUNIN,data3,5,1,IRET,'RSH29 STRV SALN')
c
      if(data3(1).eq.bmiss)then
          RSH29 = -99.99
       else
       rsh29 = data3(1)
      endif
c
      if(data3(2).eq.bmiss)then
         strv= -99.99
        else
         strv = data3(2)
      endif
c
      if(data3(3).eq.BMISS)then
          saln = -99.99
        else
        saln = data3(3)
      endif
c
c Get precipitation data
c
      CALL UFBINT(LUNIN,data1,2,500,nlev2,'TOPC DDPC')
c
       if(data1(1,1).ge.bmiss)then
           precp = -99.99
       else
         precp = data1(1,1)
       endif
c
      if(data1(2,1).eq.bmiss)then
        precpi = -99.99
      else
        precpi = data1(2,1)
      endif
c
c,..  Get Wind data
c
      CALL UFBint(LUNIN,data2,4,500,nlev3,
     *            'WDIR WSPD DDWD DDWS')
c
       if(data2(1,1).eq.bmiss)then
           wdir = -99.99
        else
          wdir = data2(1,1)
       endif
c
      if (data2(2,1).eq.bmiss)then
           wspd = -99.99
        else
         wspd = data2(2,1)
      endif
c                                                                       
c.. Print out level information
c

c
c           write(10,60)name,clat,clon,iyr,mon,iday,ihr,
c    *          mins,dchgd,rstgd,precp,tmdb,sst,wspd,wdir
c  60 format(a8,1x,2f8.2,1x,i4,4i2,1x,f8.3,f7.3,1x,f10.6,4f7.2)
c
       if(dchgd.ne.-99.99) then
           write(60,60)name,clat,clon,dchgd,rstgd
     $      ,iyr,mon,iday,ihr
   60      format(a8,1x,f8.2,1x,f8.2,1x,f8.3,1x,f8.3,4i5)

           count = count +1

        else
           nocount=nocount+1
        endif
c
       ENDDO
      ENDDO    

           write(61,50) count
           write(61,50) nocount
 50        format(i6)
C  WHEN FINISHED MAKE SURE ALL BUFFERS ARE FLUSHED THEN EXIT            
C  ---------------------------------------------------------            
                                                                        
c
      CALL CLOSBF(LUNIN)                                                
c
      PRINT*,'*** PROCESSING ENDED NORMALLY ***'                        
      STOP                                                              
      END                                                               
