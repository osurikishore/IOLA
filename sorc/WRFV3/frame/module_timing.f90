


MODULE module_timing

   INTEGER, PARAMETER, PRIVATE :: cnmax = 30
   INTEGER, PRIVATE :: cn = 0 
   REAL, PRIVATE    :: elapsed_seconds , elapsed_seconds_total = 0

   REAL(kind=8) :: epoch_seconds_hires(cnmax)

CONTAINS

   SUBROUTINE init_module_timing
      
      
      call init_hires_timer()
      cn = 0
   END SUBROUTINE init_module_timing


   SUBROUTINE start_timing
     use module_wrf_error, only: silence

      IMPLICIT NONE

      if(silence/=0) return
      cn = cn + 1
      IF ( cn .gt. cnmax ) THEN
        CALL wrf_error_fatal3("<stdin>",30,&
'module_timing: clock nesting error (too many nests)' )
        RETURN
      ENDIF

      call hires_timer(epoch_seconds_hires(cn))

   END SUBROUTINE start_timing

   SUBROUTINE end_timing(string)
      CHARACTER *(*) :: string
      CALL end_timing_fmt(string)
   END SUBROUTINE end_timing

   SUBROUTINE end_timing_fmt ( string,fmt)
     use module_wrf_error, only: silence, stderrlog, buffered
   
      IMPLICIT NONE
      character*(*), intent(in), optional :: fmt
      REAL(kind=8) :: now_hires
      CHARACTER *(*) :: string
      character*512 :: buf

      if(silence/=0) return

      IF ( cn .lt. 1 ) THEN
        CALL wrf_error_fatal3("<stdin>",56,&
'module_timing: clock nesting error, cn<1' ) 
      ELSE IF ( cn .gt. cnmax ) THEN
        CALL wrf_error_fatal3("<stdin>",59,&
'module_timing: clock nesting error, cn>cnmax' ) 
      ENDIF

      call hires_timer(now_hires)
      
      elapsed_seconds = REAL(now_hires-epoch_seconds_hires(cn))
      elapsed_seconds_total = elapsed_seconds_total + elapsed_seconds

      
3031 format("Timing for ",A,": ",F10.5," elapsed seconds")

      if(buffered/=0) then
         if(present(fmt)) then
            write(buf,fmt=fmt) TRIM(string),elapsed_seconds
         else
            write(buf,3031) TRIM(string),elapsed_seconds
         endif
         call wrf_message(buf)
      else
         if(present(fmt)) then
            if(stderrlog/=0) &
                 write(0,fmt=fmt) TRIM(string),elapsed_seconds
            write(6,fmt=fmt) TRIM(string),elapsed_seconds
         else
            if(stderrlog/=0) &
                 write(0,3031) TRIM(string),elapsed_seconds
            write(6,3031) TRIM(string),elapsed_seconds
         endif
      endif





      cn = cn - 1

   END SUBROUTINE end_timing_fmt

   FUNCTION now_time() result(timef)
     
     
     
     
     implicit none
     real*8 :: timef
     call hires_timer(timef)
   END FUNCTION now_time

END MODULE module_timing

