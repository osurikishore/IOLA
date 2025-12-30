
!******************************************************************************
! 
! Diagnostic post-processor (diapost) for WRF-NMM & its varieties
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This program is a post-processor for WRF-NMM & its varieties (HWRF etc.)
! with a diagnostic package designed for hurricane research.  Starting from 
! the model's NETCDF output on a rotated E-grid at sigma-p levels, variables 
! are extracted and converted to an A-grid for analysis & display with GrADS, 
! they are also interpolated and output to pressure and height levels. 
! 
! The program produces text data files and descriptions for hurricane tracks: 
! easy_track.msg describes the data easy_track.txt of a simple storm tracker, 
! atcf_track.msg describes atcf_track.txt with same info in different format, 
! trker_dim.msg describes trker_dim.txt for dimensions of GFDL tracker fields. 
! 
! The program also produces binary data files with GrADS control files, e.g., 
! sigma_lev.dat contains sigma-level data with sigma_rot.ctl for display on 
! the native rotated grid, and sigma_std.ctl for interpolated standard grid;
! while press_lev.dat, press_rot.ctl & press_std.ctl are for pressure levels,
! and height_lev.dat, height_rot.ctl & height_std.ctl are for height levels. 
! 
! Special data sets are also produced for evaluation and demonstration, 
! including basic_rot.dat, basic_rot.ctl & basic_std.ctl for the HFIP test, 
! and wind_swath.dat, swath_rot.ctl & swath_std.ctl for hires wind-swath 
! when the processed data are on moving nests. Some CFADs are also provided 
! with cfads_demo.dat & cfads_demo.ctl. 
! 
! The program (diapost) is driven by the shell script Diapost, which has 
! a list of 0/1 switches for choosing the types of output and diagnostics. 
! For debugging purpose, the sigma_level switch has an extra option of -1 
! to output the prognostic variables only; it works with all other switches 
! turned off except the wind_swath switch. 
! 
! -----------------------------------------------------------------------------
! 
! 2008/12/01, initial version (1.0) by Kevin Yeh <Kevin.Yeh@noaa.gov> 
! -- original code: John Michalakes' diffwrf.F90 for io_netcdf 
! -- include a simple storm tracker modified from Gopal's code 
! 
!______________________________________________________________________________


PROGRAM diapost !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

   USE  wrf_data
   IMPLICIT NONE
#include "wrf_status_codes.h"
#include "netcdf.inc"
   Integer, external :: iargc 
   Logical, external :: iveceq 
   Character(len=2)  :: storm_num 
   Character(len=4)  :: staggering 
   Character(len=3)  :: ordering, ord 
   Character(len=80) :: SysInfo, dimnames(3) 
   Character(len=160):: filename, argument 
   Character(len=19) :: DateStr, VarName, variable(99) 
   Integer:: env=0               ! 1 for removing the mean motion
   Integer:: iorder  =  4        !* 2/4: interpolation order 
   Integer:: itermax =  10       !* max number of iterations 
   Real(4):: refine  = .10       !* refining threshhold (deg) 
   Real(4):: aspect  = .80       !* output domain aspect ratio 
   Real(4):: r_slope =  20       !* max storm axis slope (dr/dz) 
   Real(4):: r_storm = 2.0       !* reference storm radius (deg) 
   Real(4):: azimuth = 3.0       !* azimuthal/radial dim ratio 
   Integer:: kmb  = 4            !* number of HFIP basic levels 
   Real(4):: pbasic(4)           !* HFIP basic pressure levels (mb) 
      DATA   pbasic / 850, 700, 500, 200 / 
   Integer:: kman = 9            !* number of mandatory levels 
   Real(4):: pman(9)             !* mandatory pressure  levels (mb) 
      DATA   pman / 1000, 850, 700, 500, 400, 300, 200, 150, 100 / 
   Integer:: kmp  = 37           !* # of 25mb pressure  levels 
   Real(4):: plevel(37)          !* 25mb pressure levels (mb) 
      DATA   plevel / 1000, 975, 950, 925, 900, 875, 850, 825, &
                       800, 775, 750, 725, 700, 675, 650, 625, &
                       600, 575, 550, 525, 500, 475, 450, 425, &
                       400, 375, 350, 325, 300, 275, 250, 225, &
                       200, 175, 150, 125, 100  / 
   Integer:: kmz  = 60           !*  # of  height levels 
   Real(4):: zlevel(60)          !* output height levels (m) 
      DATA   zlevel / 35,   70,   & 
               100,  200,  300,  400,  500,  600,  700,  800,  900,  &
              1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800,  &
              3000, 3300, 3600, 3900, 4200, 4500, 4800, 5100, 5400, 5700,  &
              6000, 6300, 6600, 6900, 7200, 7500, 7800, 8100, 8400, 8700,  &
              9000, 9500,10000,10500,11000,11500,12000,12500,13000,13500,  &
             14000,14500,15000,15500,16000,16500,17000,17500,18000   / 
   Real(4)::  undef  = 1.E+25    !* undefined value for GrADS
   Real(4):: rd_air  = 287.053   !* dry-air gas  constant (J/kg*K) 
   Real(4):: cp_air  = 1005.65   !* dry-air heat capacity (J/kg*K)
   Real(4):: kd_air  = 0.28544   !* dry-air constant (Rd/Cp: N.D.)
   Real(4):: kelv0c  = 273.15    !* Kelvin temperature (K) at 0'C 
   Real(4):: vapr0c  = 6.1121    !* water vapor pressure (mb) at 0'C 
   Real(4):: heat0c  = 2.5E+6    !* w.vapor latent heat (J/kg) at 0'C 
   Real(4):: viscos  = 1.5E-5    !* kinetic viscosity of air (/s) 
   Real(4):: karman  = 0.4000    !* PBL von Karman constant 
   Real(4):: r_earth = 6371220.  !* average earth radius (m) 
   Real(4):: g_earth = 9.80616   !* average earth gravity  (m/s^2) 
   Real(4):: o_earth = 7.292E-5  !* average earth rotation (rad/s) 
   Real(8):: pi_cnst = 3.141592653589793238462643383279502884 
   Real(8):: deg2rad = 3.141592653589793238462643383279502884/180. 
   Logical:: sigma_level, press_level, height_level, debug_data, defined 
   Logical:: wind_swath, basic_level, diagnostics, height_cylin, tilt_cylin 
   Logical:: derive_vars, press_vars, height_vars, bad_cylin, define_u10 
   Integer:: msigma, mpress, mheight, mswath, mbasic, mdiagno, mzcylin, mtilt 
   Integer:: dh0, Status, Status_next_time,Status_next_var,mdebug,ntime,process_freq_sec
   Integer:: WrfType, start_index(4), end_index(4), ndim, icross, ierror 
   Integer:: domain_num, start_date, ifcst_hour, ifcst_minute, lat10, lon10 
   Integer:: mvars, imax, jmax, kmax, ide, jde, kde, im, jm, km, ibd, jbd 
   Integer:: nvars, imx, jmx, ims, jms, imr, jmr, im3, jm3, k2d, imc, jmc 
   Integer:: i, j, k, m, n,  ii, jj, i0, i1, i2, j0, j1, j2, k0, itmp, jtmp, kstart 
   Integer:: imid, jmid, m0, m1, m2, m3, m4, n0, n1, n2, n3, n4, iold, jold
   Integer:: kkont,mkont,akont,bkont,kkk
   Integer:: iloc,jloc,ioff,joff,iter_cen,iter_eqt,mrange,mbins,kout,iter,iloc1,jloc1 
   Real(4):: ETIME, atime(2), btime(2), cpu_time_a, cpu_time_b, user_time 
   Real(4):: ctr_lat, ctr_lon, dlat, dlon, yrlen, xrs11, yrs11, xrs31, vcut 
   Real(4):: fcst_time, storm_lat, storm_lon, vmax, pmax, pmin, pcut, d01cut 
   Real(4):: olamda, otheta, alamda, atheta, dlamda, dtheta, r_cylin, d02cut 
   Real(4):: xrr1, xrr2, xrr8, xrr9, yrr1, yrr2, yrr8, yrr9, percent, binsize 
   Real(4):: xxr1, xxr_, yyr1, yyr_, xxs1, xxs_, yys1, yys_, dxx, dyy, vrange 
   Real(4):: eps, ykappam, qqq, ttt, tt0, ttx, tts, tts_error, scale, rho 
   Real(4):: ff(4,4), wtx(4), wty(4), aa1, aa2, bb1, bb2, ees, ee0, z02, z10 
   Real(4):: aaa, bbb, ccc, ddd, eee, fff, xxx, yyy, zzz, www, rrr, uur,vvr, zz0  
   Real(4):: fcor,vg,atemp,btemp,ctemp,dtemp,etemp,ftemp,gtemp,htemp,utemp,vtemp,huur,hvvr          ! gopal's doing  
   Real(4):: uenv,venv
   Real(4):: lat1(20000),lon1(20000) 

!
!=============================================================================================================
!  HDOBS SECTION: if HDOBS are available open the file, proceed as usual
!  otherwise
!
!  README: HDOBS comparison requires 
!  (1) file 1: out_hdob_YYYYMODYHR copied as hdobs.dat (one directory per file - does not do multiple windows) 
!  (2) file 2: obshead.dat:: maxobs,flight_start,flight_end,# eye penetration, HWRF simulation day (1/2/3/4/5)
!  (3) output files: flight_track*.txt, verify_vars1*.txt and verify_vars2*.txt (verified at forecast hour) 
!==============================================================================================================
   Character(len=160)     :: filename1,filename2,filename3
   Integer                :: maxobs,flight_start,flight_end,nday  ! (e.g. format: 980 075030 153000 2)
   Integer, parameter     :: nrad=4        !radius of obs influence in grid units 
   Integer                :: proceed,local_clock,iset1,iset2,flegs
   Integer, Allocatable   :: hdtime(:)
   Real(4)                :: obs1,obs2,obs3,obs4,obs5,obs6,obs7,obs8,obs9,obs10,obs11,obs12,obs13
   Real(4)                :: latmin,latmax,lonmin,lonmax
   Real(4), Allocatable   :: hdlat(:),hdlon(:),hdppp(:),hdgpm(:),hdslp(:),hdrain(:)
   Real(4), Allocatable   :: hdtt(:),hdtd(:),hdwd(:),hdws(:)
   
!
!* Generalized routine includes E2A and B2A. This is Thiago's doings
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!* Variables used to determine the generating model and data grid type. Thiago 11/12/2014.
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   integer,parameter :: GRID_TYPE_ROTATEDLLE=1      !Rotated Lat-Lon projection with Arakawa E-grid staggering.
   integer,parameter :: GRID_TYPE_ROTATEDLLB=2      !Rotated Lat-Lon projection with Arakawa B-grid staggering.
   integer :: igridtype                             !The ID of the projection and grid type of the underlying data.
   Character(len=:),allocatable :: smodelname       !The name of the model which generated this data set.
 
!* 2D variables on rotated A-grid 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   Real(4), Allocatable::  zs(:,:),  sm(:,:),  p0(:,:), prc(:,:), p0B(:,:)  ! for BEST_MSLP 
   Real(4), Allocatable::  ps(:,:), ths(:,:),  ts(:,:),  qs(:,:) 
   Real(4), Allocatable::  u10(:,:), v10(:,:), t10(:,:), q10(:,:) 
   Real(4), Allocatable::  sst(:,:), ssq(:,:), t02(:,:), q02(:,:) 
   Real(4), Allocatable::  tw(:,:),  qw(:,:),  uf(:,:),  z0(:,:), z0t(:,:) 
   Real(4), Allocatable::  tx(:,:), ty(:,:)

   
!* 3D variables on sigma-p levels 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   Real(4), Allocatable:: uu(:,:,:), vv(:,:,:), ww(:,:,:) 
   Real(4), Allocatable:: tt(:,:,:), td(:,:,:), qt(:,:,:) 
   Real(4), Allocatable:: pp(:,:,:), qq(:,:,:), cw(:,:,:)
   Real(4), Allocatable:: gh(:,:,:), rh(:,:,:), dp(:,:,:) 
   Real(4), Allocatable:: f_ice(:,:,:),f_rain(:,:,:),f_rimef(:,:,:)
   Real(4), Allocatable::                            mc_heat(:,:,:)
   Real(4), Allocatable:: tn(:,:,:), qr(:,:,:), qc(:,:,:), qi(:,:,:), qn(:,:,:) ! RAB 2/27/12

!* 3D diffusion variables on sigma-p levels
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   Real(4), Allocatable::  difu(:,:,:),  difv(:,:,:),km_diffus(:,:,:),kh_diffus(:,:,:)  
   Real(4), Allocatable:: hdifu(:,:,:), hdifv(:,:,:) 

!* Coordinates and general workspace 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   Real(4), Allocatable::  xrr(:,:),   yrr(:,:),  xrs(:,:),  yrs(:,:) 
   Real(4), Allocatable::  xsr(:,:),   ysr(:,:),  xss(:,:),  yss(:,:) 
   Real(4), Allocatable::  xxr(:,:),   yyr(:,:),  xxs(:,:),  yys(:,:) 
   Real(4), Allocatable:: wcos(:,:),  wsin(:,:),  xrlen(:),  td2(:,:)  
   Real(4), Allocatable:: u10x(:,:),  v10x(:,:),  p0x(:,:),  p02(:,:) 
   Real(4), Allocatable:: uus(:,:,:), vvs(:,:,:), urs(:,:),  vrs(:,:) 
   Real(4), Allocatable:: dat(:,:,:), wks(:,:,:),  wk(:,:)  
   Real(4), Allocatable:: txx(:,:),  tyx(:,:)
   
!* 3D workspace for pressure and height levels 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   Real(4), Allocatable:: uu_(:,:,:), vv_(:,:,:), ww_(:,:,:) 
   Real(4), Allocatable:: tt_(:,:,:), td_(:,:,:), qt_(:,:,:) 
   Real(4), Allocatable:: pp_(:,:,:), qq_(:,:,:),cw_(:,:,:)
   Real(4), Allocatable:: gh_(:,:,:), rh_(:,:,:) 
   Real(4), Allocatable:: ux_(:,:,:), vx_(:,:,:) 
   Real(4), Allocatable:: f_ice_(:,:,:),f_rain_(:,:,:),f_rimef_(:,:,:)
   Real(4), Allocatable::                              mc_heat_(:,:,:)

!* 3D diffusion variables on height levels (gopal's doing) 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   Real(4), Allocatable::  difu_(:,:,:),  difv_(:,:,:),km_diffus_(:,:,:),kh_diffus_(:,:,:)
   Real(4), Allocatable:: hdifu_(:,:,:), hdifv_(:,:,:)

!* 3D workspace for triple grid refinement 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   Real(4), Allocatable:: uu3(:,:,:), vv3(:,:,:), ww3(:,:,:) 
   Real(4), Allocatable:: tt3(:,:,:), td3(:,:,:), qt3(:,:,:) 
   Real(4), Allocatable:: pp3(:,:,:), qq3(:,:,:), cw3(:,:,:) 
   Real(4), Allocatable:: gh3(:,:,:), rh3(:,:,:) 
   Real(4), Allocatable:: ux3(:,:,:), vx3(:,:,:) 
   Real(4), Allocatable:: f_ice3(:,:,:),f_rain3(:,:,:),f_rimef3(:,:,:)
   Real(4), Allocatable::                              mc_heat3(:,:,:)  

   Real(4), Allocatable:: uu3_(:,:,:), vv3_(:,:,:), ww3_(:,:,:) 
   Real(4), Allocatable:: tt3_(:,:,:), td3_(:,:,:), qt3_(:,:,:) 
   Real(4), Allocatable:: pp3_(:,:,:), qq3_(:,:,:), cw3_(:,:,:) 
   Real(4), Allocatable:: gh3_(:,:,:), rh3_(:,:,:) 
   Real(4), Allocatable:: ux3_(:,:,:), vx3_(:,:,:) 
   Real(4), Allocatable:: f_ice3_(:,:,:),f_rain3_(:,:,:),f_rimef3_(:,:,:)
   Real(4), Allocatable::                                mc_heat3_(:,:,:)

!* 3D diffusion variables for triple grid refinement (gopal's doing) 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   Real(4), Allocatable::  difu3(:,:,:),  difv3(:,:,:),km_diffus3(:,:,:),kh_diffus3(:,:,:)
   Real(4), Allocatable:: hdifu3(:,:,:), hdifv3(:,:,:)

!* Workspace for variables on cylindrical grid 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   Integer, Allocatable::  icen(:),    jcen(:) 
   Real(4), Allocatable::  rcos(:),    rsin(:),   rhc(:,:,:) 
   Real(4), Allocatable:: ppc(:,:,:), qqc(:,:,:), cwc(:,:,:) 
   Real(4), Allocatable:: ttc(:,:,:), tdc(:,:,:), qtc(:,:,:) 
   Real(4), Allocatable:: uuc(:,:,:), vvc(:,:,:), wwc(:,:,:) 
   Real(4), Allocatable:: uxc(:,:,:), vxc(:,:,:), mwc(:,:,:),qwc(:,:,:) 
   Real(4), Allocatable:: f_icec(:,:,:),f_rainc(:,:,:),f_rimefc(:,:,:)
   Real(4), Allocatable::                              mc_heatc(:,:,:)

!* 3D workspace for height levels (gopal's doing)
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   Real(4), Allocatable:: u1d(:),v1d(:),t1d(:),rad1d(:)
   Real(4), Allocatable:: uu2d(:,:,:),vv2d(:,:,:),ww2d(:,:,:),pp2d(:,:,:),tt2d(:,:,:),vort2d(:,:,:),th2d(:,:,:)
   Real(4), Allocatable:: upert(:,:,:),vpert(:,:,:),wpert(:,:,:),vopert(:,:,:),thpert(:,:,:)
   Real(4), Allocatable:: dvdz(:,:,:),edrvf(:,:,:),eddvdz(:,:,:)
   Real(4), Allocatable:: mean_var1(:,:,:),rms_var1(:,:,:),wstat2d(:,:),wstat(:,:,:)
   Real(4), Allocatable:: mean_var2(:,:,:),rms_var2(:,:,:),tstat2d(:,:),tstat(:,:,:)
   Real(4), Allocatable:: lfluxz(:),hfluxz(:),xstress(:),ystress(:),stress(:)
   Real(4), Allocatable:: warmxz(:,:),uvortxz(:,:)
   Real(4), Allocatable:: udiffc(:,:,:), vdiffc(:,:,:),uadc(:,:,:),tadc(:,:,:)
   Real(4), Allocatable:: hudiffc(:,:,:),hvdiffc(:,:,:)
   Real(4), Allocatable:: vadc(:,:,:), tan_grd_wind(:,:,:)
   Real(4), Allocatable:: tpert(:,:,:), vv_verify(:,:,:),pgr(:,:,:),vgr(:,:,:),ror(:,:,:)
   Real(4), Allocatable:: kmc(:,:,:),khc(:,:,:)

!* For SEE computation 
!* ~~~~~~~~~~~~~~~~~~~
   Real(4), Allocatable:: zdiff2d(:,:,:),hdiff2d(:,:,:)
   Real(4), Allocatable:: thdot(:,:,:),edth1(:,:,:),edth2(:,:,:),edth3(:,:,:),fthdiff(:,:,:),fmomdiff(:,:,:)
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

   Real(4), Allocatable:: wwfr(:,:) 
   
!           *     *     *     *     *     *     *     *     *     * 
   
   cpu_time_a = ETIME(atime) 
   CALL set_wrf_debug_level( 2 )
   CALL ext_ncd_ioinit(SysInfo,Status)
   iter_cen = 0 ; iter_eqt = 0 

!**************************************! 
!  I/O preparation for data processing ! 
!______________________________________! 
   
!* Get arguments from script
!* ~~~~~~~~~~~~~~~~~~~~~~~~~ 
   CALL getarg( 1,filename)
   CALL getarg( 2,argument) ; READ(argument,*) msigma
   CALL getarg( 3,argument) ; READ(argument,*) mpress
   CALL getarg( 4,argument) ; READ(argument,*) mheight
   CALL getarg( 5,argument) ; READ(argument,*) mswath
   CALL getarg( 6,argument) ; READ(argument,*) mbasic
   CALL getarg( 7,argument) ; READ(argument,*) mzcylin
   CALL getarg( 8,argument) ; READ(argument,*) mtilt
   CALL getarg( 9,argument) ; READ(argument,*) mdiagno
   CALL getarg(10,argument) ; READ(argument,*) mdebug
   CALL getarg(11,argument) ; READ(argument,*) domain_num
   CALL getarg(12,argument) ; READ(argument,*) storm_num
   CALL getarg(13,argument) ; READ(argument,*) start_date
   CALL getarg(14,argument) ; READ(argument,*) r_cylin
   CALL getarg(15,argument) ; READ(argument,*) d02cut
   CALL getarg(16,argument) ; READ(argument,*) d01cut
   CALL getarg(17,argument) ; READ(argument,*) ntime
   CALL getarg(18,argument) ; READ(argument,*) process_freq_sec 
   
   derive_vars  = ( mswath  >= 0 ) 
   sigma_level  = ( msigma  == 1 ) 
   press_level  = ( mpress  == 1 ) 
   basic_level  = ( mbasic  == 1 ) 
   height_level = ( mheight == 1 ) 
   debug_data   = ( mdebug  == 1 ) 
   diagnostics  = ( mdiagno == 1 ) 
   tilt_cylin   = ( mtilt   == 1 ) 
   wind_swath   = ( mswath  /= 0 ) ! .AND. ( domain_num > 1 ) 
   height_cylin = ( mzcylin == 1 ) ! .AND. ( domain_num > 1 )
   height_vars  = ( height_level ) .OR.  ( height_cylin ) 
   press_vars   = ( press_level  ) 
   
!* Variables to be extracted 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~ 
   variable( 1) = 'U10'
   variable( 2) = 'V10'
   variable( 3) = 'MSLP'
   variable( 4) = 'HLAT'
   variable( 5) = 'HLON'
   variable( 6) = 'XTIME'
   variable( 7) = 'FIS'
   variable( 8) = 'SM'
   variable( 9) = 'SST'
   variable(10) = 'QS'
   variable(11) = 'THS'
   variable(12) = 'Z0'
   variable(13) = 'USTAR'
   variable(14) = 'TWBS'
   variable(15) = 'QWBS'
   variable(16) = 'PREC'
   variable(17) = 'U'
   variable(18) = 'V'
   variable(19) = 'W'
   variable(20) = 'T'
   variable(21) = 'Q'
   variable(22) = 'PINT'
   variable(23) = 'CWM'

! budget related (gopal's doing)

   variable(24) = 'DUDT'           !'DUBLDT' per second
   variable(25) = 'DVDT'           !'DVBLDT' per second
   variable(26) = 'TAUX'
   variable(27) = 'TAUY'
   variable(28) = 'F_ICE'
   variable(29) = 'F_RAIN'
   variable(30) = 'F_RIMEF'
   variable(31) = 'KM'
   variable(32) = 'KH'
   variable(33) = 'DUHDIFFDT'
   variable(34) = 'DVHDIFFDT'
   variable(35) = 'T_CLOUD'  !heating per time step
   variable(36) = 'BEST_MSLP'

   mvars = 6 ; IF (derive_vars) mvars = 36

!* Open the data files to write to 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   OPEN(51,file='easy_track.msg',status='unknown')
   OPEN(52,file='easy_track.csv',position='append')
   OPEN(53,file='atcf_track.msg',status='unknown')
   OPEN(54,file='atcf_track.csv',position='append')
   OPEN(55,file='maxmin.csv',position='append')
   OPEN(56,file='fluxes.csv',position='append')
   OPEN(57,file='uvenv.csv',position='append')
   OPEN(20,file='shear.csv',position='append')
   OPEN(21,file='wstats.csv',position='append')
   OPEN(22,file='tstats.csv',position='append')

!  OPEN(30,file='synthetic.dat',position='append')
!  OPEN(23,file='extreme.txt',position='append')
!  OPEN(24,file='updraft.txt',position='append')

   if ( wind_swath ) then 
      if ( debug_data ) then 
      OPEN(61,file='swath_rot.ctl',form='formatted')
      OPEN(62,file='swath_rot.dat',form='unformatted',position='append')
      endif 
   OPEN(63,file='swath_std.ctl',form='formatted')
   OPEN(64,file='swath_std.dat',form='unformatted',position='append')
   endif 

   if ( sigma_level ) then 
   OPEN(100,file='sigma_lev.dat',form='unformatted',position='append')
   OPEN(101,file='sigma_rot.ctl',form='formatted')
   OPEN(102,file='sigma_std.ctl',form='formatted')
   endif 

   if ( press_level ) then 
   OPEN(200,file='press_lev.dat',form='unformatted',position='append')
   OPEN(201,file='press_rot.ctl',form='formatted')
   OPEN(202,file='press_std.ctl',form='formatted')
   endif 

   if ( basic_level ) then 
      if ( debug_data ) then 
      OPEN(251,file='basic_rot.ctl',form='formatted')
      OPEN(252,file='basic_rot.dat',form='unformatted',position='append')
      endif 
   OPEN(253,file='basic_std.ctl',form='formatted')
   OPEN(254,file='basic_std.dat',form='unformatted',position='append')
   OPEN(255,file='basic_std.msg',status='unknown')
   OPEN(256,file='basic_std.txt',position='append')
   endif 
   
   if ( height_level ) then 
   OPEN(300,file='height_lev.dat',form='unformatted',position='append')
   OPEN(301,file='height_rot.ctl',form='formatted')
   OPEN(302,file='height_std.ctl',form='formatted')
   endif 

   if ( height_cylin ) then 
   OPEN(310,file='height_cyl.dat',form='unformatted',position='append')
   OPEN(311,file='height_cyl.ctl',form='formatted')
   OPEN(312,file='height_cpr.ctl',form='formatted')
   endif 

   if ( diagnostics ) then 
   OPEN(331,file='cfads_demo.ctl',form='formatted')
   OPEN(332,file='cfads_demo.dat',form='unformatted',position='append')
   endif 

!* Open the data file to read from 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   CALL ext_ncd_open_for_read( trim(filename), 0, 0, "", dh0, Status )

   PRINT* 
   PRINT*, 'Reading ', trim(filename) 
   IF ( Status /= 0 ) THEN
      PRINT*, 'Error opening', filename, ' Status = ', Status
      STOP
   ENDIF
   PRINT* 
!
!* Generalized routine includes E2A and B2A. This is Thiago's doings
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!* Get the model name and determine the underlying data projection and grid type. Thiago 11/12/2014.
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   WRITE(*,*)"Searching for global netCDF attribute 'MODELNAME'."
   CALL ext_ncd_get_dom_ti_char( dH0,"MODELNAME",argument,ierror) !Search for the model name attribute.
   IF(ierror==0) THEN   !Model name attribute was found.
       IF(trim(argument)=="NMMB") THEN
           WRITE(*,*)"Processing data file produced by the 'NMMB' model on rotated B-grid."
           smodelname="NMMB"                !This data set was produced by the "NMMB" model.
           igridtype=GRID_TYPE_ROTATEDLLB   !NMMB uses rotated LAT-LON projection on B-grid staggering.
       ELSE
           WRITE(*,*)"Error: NETCDF data sets from model '",trim(argument),"'are not currently supported by Diapost."
           STOP
       ENDIF
   ELSE !<-- Default
       WRITE(*,*)"Cannot find global netCDF attribute 'MODELNAME'. Assume HWRF on rotated E-grid."
       smodelname="HWRF"                    !Default, assume HWRF model.
       igridtype=GRID_TYPE_ROTATEDLLE       !Default, assume rotated LAT-LON projection on E-grid staggering.
   ENDIF

!* Dimensions & resolutions of data grids 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   CALL ext_ncd_get_dom_ti_integer( dh0,'WEST-EAST_GRID_DIMENSION',   &
                                         imax, 1, icross, ierror      )
   CALL ext_ncd_get_dom_ti_integer( dh0,'SOUTH-NORTH_GRID_DIMENSION', &
                                         jmax, 1, icross, ierror      )
   CALL ext_ncd_get_dom_ti_integer( dh0,'BOTTOM-TOP_GRID_DIMENSION',  &
                                         kmax, 1, icross, ierror      )
   CALL ext_ncd_get_dom_ti_real( dh0, 'CEN_LAT', ctr_lat, 1, icross, ierror )
   CALL ext_ncd_get_dom_ti_real( dh0, 'CEN_LON', ctr_lon, 1, icross, ierror )
   CALL ext_ncd_get_dom_ti_real( dh0, 'DX', dlon, 1, icross, ierror )
   CALL ext_ncd_get_dom_ti_real( dh0, 'DY', dlat, 1, icross, ierror ) 
   
   ide = imax - 1 ; im = 2*ide - 1    !* X dimensions of E-grid & A-grid 
   jde = jmax - 1 ; jm =   jde        !* Y dimensions of E-grid & A-grid 
   kde = kmax - 1 ; km =   kde        !* Number of sigma levels (layers) 

!
!* Generalized routine includes E2A and B2A. This is Thiago's doings
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  For rotated LAT-LON projections on B grid, set "im=ide" because the wind grid will
!  be interpolated to the mass grid. The mass grid is considered to be unstaggered.
!
   IF(igridtype==GRID_TYPE_ROTATEDLLB) im=ide   !Thiago 11/12/2014.

   
!* Allocate memory for 2D & sigma-level variables 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   Allocate(  zs(im,jm),  sm(im,jm),  p0(im,jm), p0B(im,jm), prc(im,jm) )  ! p0B is BEST_MSLP
   Allocate(  ps(im,jm), ths(im,jm),  ts(im,jm),  qs(im,jm) ) 
   Allocate( u10(im,jm), v10(im,jm), t10(im,jm), q10(im,jm) ) 
   Allocate( sst(im,jm), ssq(im,jm), t02(im,jm), q02(im,jm) ) 
   Allocate(  tw(im,jm),  qw(im,jm),  uf(im,jm),  z0(im,jm), z0t(im,jm) ) 
   Allocate( uu(im,jm,km), vv(im,jm,km), ww(im,jm,km) ) 
   Allocate( tt(im,jm,km), td(im,jm,km), qt(im,jm,km) ) 
   Allocate( pp(im,jm,km), qq(im,jm,km), cw(im,jm,km) ) 
   Allocate( gh(im,jm,km), rh(im,jm,km) ) 
   Allocate( tx(im,jm), ty(im,jm))
   Allocate( f_ice(im,jm,km),f_rain(im,jm,km),f_rimef(im,jm,km))
   Allocate(                                  mc_heat(im,jm,km))
   Allocate( tn(im,jm,km), qc(im,jm,km), qr(im,jm,km), qi(im,jm,km), qn(im,jm,km) )   ! RAB 2/21/12
   
!* Budget computation (gopal's doing)
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   Allocate(  difu(im,jm,km),  difv(im,jm,km), km_diffus(im,jm,km), kh_diffus(im,jm,km))
   Allocate( hdifu(im,jm,km), hdifv(im,jm,km))
   
!           *     *     *     *     *     *     *     *     *     * 

   CALL ext_ncd_get_next_time( dh0, DateStr, Status_next_time )
   
   DO WHILE ( Status_next_time .EQ. 0 ) !==================================> 

!*************************************! 
!  Read the original data from E-grid !  
!_____________________________________! 
   
   Allocate( u10x(ide,jde), v10x(ide,jde), dp(im,jm,km) ) 
   Allocate( txx(ide,jde), tyx(ide,jde))

   CALL ext_ncd_get_next_var( dh0, VarName, Status_next_var )

   DO WHILE ( Status_next_var .EQ. 0 ) !-----------------------------------> 

!*    Get the variable information 
!*    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      start_index= 1 ; end_index= 1
      call ext_ncd_get_var_info( dh0, VarName, ndim, ordering, staggering, &
                                 start_index, end_index, WrfType, ierror   )

      IF ( WrfType /= WRF_REAL .AND. WrfType /= WRF_DOUBLE) THEN
         CALL ext_ncd_get_next_var( dh0, VarName, Status_next_var )
         CYCLE
      ENDIF

      do m = 1, mvars 
      IF ( VarName == variable(m) ) THEN  !* - - - - - - - - - - - - - - 

      IF ( ndim == 3 )  ord = 'XYZ'
      IF ( ndim == 2 )  ord = 'XY'
      IF ( ndim == 1 )  ord = 'Z'
      IF ( ndim == 0 )  ord = '0'

      Allocate( dat( end_index(1),end_index(2),end_index(3)   ) )
      Allocate( wks( end_index(1),end_index(2),end_index(3)-1 ) )

!*    Read original data on rotated E-grid 
!*    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL ext_ncd_read_field( dh0, DateStr, trim(VarName), dat , WRF_REAL,  &
                        0,0,0,ord,staggering,dimnames,start_index,end_index, & 
                        start_index,end_index,start_index,end_index, ierror  ) 

      if ( ierror /= 0 ) then
         PRINT*, 'Error in reading data:'
         PRINT*, 'ndim, end_index(1:3) =', ndim, end_index(1:3) 
         STOP 
      endif

      WRITE(*,'(1x,A9,1x,I1,2(1x,I3),1x,I2,1x,A3,1x,A)')             &
            VarName, ndim, end_index(1), end_index(2), end_index(3), &
            trim(ordering), trim(DateStr)

!*    Convert E-grid data to enhanced A-grid 
!*    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      IF ( trim(VarName)== 'XTIME') fcst_time = dat(1,1,1)  !* (min) 
      IF ( trim(VarName)== 'HLAT' )yrs11 = dat(1,1,1)
      IF ( trim(VarName)== 'HLON' ) THEN
              xrs11 = dat(1,1,1) ; xrs31 = dat(2,1,1) 
      ENDIF 

!*
!* Generalized routine includes E2A and B2A. This is Thiago's doings
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
      IF ( trim(VarName)== 'SM' ) & 
      CALL Destagger_Grid_Scalar( sm(1,1), im, jm, 1, dat(1,1,1), ide, jde, 1 , igridtype)  !Thiago 11/12/2014.
           !CALL E2A_Grid_Scalar( sm(1,1), im, jm, 1, dat(1,1,1), ide, jde, 1 ) 
   
      IF ( trim(VarName)== 'FIS' ) THEN
      CALL Destagger_Grid_Scalar( zs(1,1), im, jm, 1, dat(1,1,1), ide, jde, 1 , igridtype)  !Thiago 11/12/2014. 
           !CALL E2A_Grid_Scalar( zs(1,1), im, jm, 1, dat(1,1,1), ide, jde, 1 ) 
      zs(:,:) = zs(:,:)/g_earth  !* (m) 
      ENDIF 
   
      IF ( trim(VarName)== 'MSLP' ) THEN
      CALL Destagger_Grid_Scalar( p0(1,1), im, jm, 1, dat(1,1,1), ide, jde, 1 ,igridtype)  !Thiago 11/12/2014. 
           !CALL E2A_Grid_Scalar( p0(1,1), im, jm, 1, dat(1,1,1), ide, jde, 1 ) 
      p0(:,:) = p0(:,:)*.01      !* (mb) 
      ENDIF 

      IF ( trim(VarName)== 'BEST_MSLP' ) THEN  ! HWRF outputs BEST_MSLP for MSLP after X hours of forecast 
      CALL Destagger_Grid_Scalar( p0B(1,1), im, jm, 1, dat(1,1,1), ide, jde, 1 ,igridtype)  !Thiago 11/12/2014. 
           !CALL E2A_Grid_Scalar( p0B(1,1), im, jm, 1, dat(1,1,1), ide, jde, 1 ) 
      p0B(:,:) = p0B(:,:)*.01      !* (mb) 
      ENDIF

      do i=1,im; do j=1,jm
        p0(i,j) = max(p0(i,j),p0B(i,j))   ! When BEST_MSLP is available use that
      enddo; enddo

      IF ( trim(VarName)== 'PREC' ) & 
      CALL Destagger_Grid_Scalar( prc(1,1), im, jm, 1, dat(1,1,1), ide, jde, 1 , igridtype) !Thiago 11/12/2014.      
           !CALL E2A_Grid_Scalar( prc(1,1), im, jm, 1, dat(1,1,1), ide, jde, 1 ) 
   
      IF ( trim(VarName)== 'QS' ) &
      CALL Destagger_Grid_Scalar( qs(1,1), im, jm, 1, dat(1,1,1), ide, jde, 1 , igridtype)  !Thiago 11/12/2014. 
           !CALL E2A_Grid_Scalar( qs(1,1), im, jm, 1, dat(1,1,1), ide, jde, 1 ) 
   
      IF ( trim(VarName)== 'THS' ) & 
      CALL Destagger_Grid_Scalar( ths(1,1), im, jm, 1, dat(1,1,1), ide, jde, 1 , igridtype) !Thiago 11/12/2014.
           !CALL E2A_Grid_Scalar( ths(1,1), im, jm, 1, dat(1,1,1), ide, jde, 1 ) 
   
      IF ( trim(VarName)== 'SST' ) &
      CALL Destagger_Grid_Scalar( sst(1,1), im, jm, 1, dat(1,1,1), ide, jde, 1 , igridtype) !Thiago 11/12/2014. 
           !CALL E2A_Grid_Scalar( sst(1,1), im, jm, 1, dat(1,1,1), ide, jde, 1 ) 
   
      IF ( trim(VarName)== 'Z0' ) & 
      CALL Destagger_Grid_Scalar( z0t(1,1), im, jm, 1, dat(1,1,1), ide, jde, 1 , igridtype) !Thiago 11/12/2014.
           !CALL E2A_Grid_Scalar( z0(1,1), im, jm, 1, dat(1,1,1), ide, jde, 1 ) 
   
      IF ( trim(VarName)== 'USTAR' ) &
      CALL Destagger_Grid_Scalar( uf(1,1), im, jm, 1, dat(1,1,1), ide, jde, 1 , igridtype)  !Thiago 11/12/2014. 
           !CALL E2A_Grid_Scalar( uf(1,1), im, jm, 1, dat(1,1,1), ide, jde, 1 ) 
   
      IF ( trim(VarName)== 'TWBS' ) &
      CALL Destagger_Grid_Scalar( tw(1,1), im, jm, 1, dat(1,1,1), ide, jde, 1 , igridtype)  !Thiago 11/12/2014. 
           !CALL E2A_Grid_Scalar( tw(1,1), im, jm, 1, dat(1,1,1), ide, jde, 1 ) 
   
      IF ( trim(VarName)== 'QWBS' ) & 
      CALL Destagger_Grid_Scalar( qw(1,1), im, jm, 1, dat(1,1,1), ide, jde, 1 , igridtype)  !Thiago 11/12/2014.
           !CALL E2A_Grid_Scalar( qw(1,1), im, jm, 1, dat(1,1,1), ide, jde, 1 ) 
   
      IF ( trim(VarName)== 'U10' ) u10x(:,:) = dat(:,:,1) 
      IF ( trim(VarName)== 'V10' ) v10x(:,:) = dat(:,:,1) 

      IF ( trim(VarName)== 'TAUX' ) txx(:,:) = dat(:,:,1)
      IF ( trim(VarName)== 'TAUY' ) tyx(:,:) = dat(:,:,1)

      IF ( trim(VarName)== 'U' ) & 
      CALL Destagger_Grid_Vector( uu(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde , igridtype) !Thiago 11/12/2014.
           !CALL E2A_Grid_Vector( uu(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde ) 
   
      IF ( trim(VarName)== 'V' ) & 
      CALL Destagger_Grid_Vector( vv(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde , igridtype) !Thiago 11/12/2014.
           !CALL E2A_Grid_Vector( vv(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde ) 
   
      IF ( trim(VarName)== 'W' ) THEN
          IF(smodelname=="NMMB") THEN
              !wks=dat   !Simply copy "dat" into temporary variable "wks".
              CALL Destagger_Grid_Scalar( ww(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde , igridtype) !Thiago 11/12/2014.
          ELSE  !<-- For HWRF.
              dat(:,:,kde+1) = dat(:,:,kde)
              DO k = 1, kde ; wks(:,:,k) = .5*(dat(:,:,k)+dat(:,:,k+1)) ; ENDDO
              CALL Destagger_Grid_Scalar( ww(1,1,1), im, jm, km, wks(1,1,1), ide, jde, kde , igridtype) !Thiago 11/12/2014.
          ENDIF
          !CALL Destagger_Grid_Scalar( ww(1,1,1), im, jm, km, wks(1,1,1), ide, jde, kde , 9999) !igridtype) !Thiago 11/12/2014.
          !CALL E2A_Grid_Scalar( ww(1,1,1), im, jm, km, wks(1,1,1), ide, jde, kde ) 
      ENDIF
 
      IF ( trim(VarName)== 'T' ) & 
      CALL Destagger_Grid_Scalar( tt(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde , igridtype) !Thiago 11/12/2014.
           !CALL E2A_Grid_Scalar( tt(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde ) 
   
      IF ( trim(VarName)== 'Q' ) & 
      CALL Destagger_Grid_Scalar( qq(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde , igridtype) !Thiago 11/12/2014.
           !CALL E2A_Grid_Scalar( qq(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde ) 
   
      IF ( trim(VarName)== 'PINT' ) THEN !* ( Pa -> mb ) 
      DO k = 1, kde ; wks(:,:,k) = .005*(dat(:,:,k)+dat(:,:,k+1)) ; ENDDO
      CALL Destagger_Grid_Scalar( pp(1,1,1), im, jm, km, wks(1,1,1), ide, jde, kde , igridtype) !Thiago 11/12/2014.
           !CALL E2A_Grid_Scalar( pp(1,1,1), im, jm, km, wks(1,1,1), ide, jde, kde ) 
      DO k = 1, kde ; wks(:,:,k) =  .005*(dat(:,:,k)-dat(:,:,k+1)) ; ENDDO
      CALL Destagger_Grid_Scalar( dp(1,1,1), im, jm, km, wks(1,1,1), ide, jde, kde , igridtype) !Thiago 11/12/2014.
           !CALL E2A_Grid_Scalar( dp(1,1,1), im, jm, km, wks(1,1,1), ide, jde, kde )
      CALL Destagger_Grid_Scalar( ps(1,1), im, jm, 1, dat(1,1,1), ide, jde, 1 , igridtype)      !Thiago 11/12/2014.
           !CALL E2A_Grid_Scalar( ps(1,1), im, jm, 1, dat(1,1,1), ide, jde, 1 ) 
      ps(:,:) = ps(:,:)*.01
      ENDIF
       
      IF ( trim(VarName)== 'CWM' ) & 
      CALL Destagger_Grid_Scalar( cw(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde , igridtype) !Thiago 11/12/2014.
           !CALL E2A_Grid_Scalar( cw(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde ) 
 
      IF ( trim(VarName)== 'F_ICE' ) &
      CALL Destagger_Grid_Scalar( f_ice(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde , igridtype) !Thiago 11/12/2014.
           !CALL E2A_Grid_Scalar( f_ice(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde )

      IF ( trim(VarName)== 'F_RAIN' ) &
      CALL Destagger_Grid_Scalar( f_rain(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde , igridtype) !Thiago 11/12/2014.
           !CALL E2A_Grid_Scalar( f_rain(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde)

      IF ( trim(VarName)== 'F_RIMEF' ) &
      CALL Destagger_Grid_Scalar( f_rimef(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde , igridtype) !Thiago 11/12/2014.
           !CALL E2A_Grid_Scalar( f_rimef(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde)

      IF ( trim(VarName)== 'T_CLOUD' ) &
      CALL Destagger_Grid_Scalar( mc_heat(1,1,1), im, jm, km, dat(1,1,1), ide,jde, kde , igridtype) ! Gopal
           !CALL E2A_Grid_Scalar( mc_heat(1,1,1), im, jm, km, dat(1,1,1), ide,
           !jde, kde)

!      IF ( trim(VarName)== 'TRAIN' ) &
!      CALL Destagger_Grid_Scalar( tn(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde , igridtype) !Thiago 11/12/2014.
!           !CALL E2A_Grid_Scalar( tn(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde ) 

      IF ( trim(VarName)== 'QCLOUD' ) &
      CALL Destagger_Grid_Scalar( qc(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde , igridtype) !Thiago 11/12/2014.
           !CALL E2A_Grid_Scalar( qc(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde ) 

      IF ( trim(VarName)== 'QRAIN' ) &
      CALL Destagger_Grid_Scalar( qr(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde , igridtype) !Thiago 11/12/2014.
           !CALL E2A_Grid_Scalar( qr(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde ) 

      IF ( trim(VarName)== 'QICE' ) &
      CALL Destagger_Grid_Scalar( qi(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde , igridtype) !Thiago 11/12/2014.
           !CALL E2A_Grid_Scalar( qi(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde ) 

      IF ( trim(VarName)== 'QSNOW' ) &
      CALL Destagger_Grid_Scalar( qn(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde , igridtype) !Thiago 11/12/2014.
           !CALL E2A_Grid_Scalar( qn(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde )

!     For diffusion budget (gopal's doing)

      IF ( trim(VarName)== 'DUDT' ) &
      CALL Destagger_Grid_Scalar( difu(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde , igridtype) !Thiago 11/12/2014.
           !CALL E2A_Grid_Scalar( difu(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde )

      IF ( trim(VarName)== 'DUBLDT' ) &
      CALL Destagger_Grid_Scalar( difu(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde , igridtype) !Thiago 11/12/2014.
           !CALL E2A_Grid_Scalar( difu(1,1,1), im, jm, km, dat(1,1,1), ide, jde,
           !kde )

      IF ( trim(VarName)== 'DVDT' ) &
      CALL Destagger_Grid_Scalar( difv(1,1,1), im, jm, km, dat(1,1,1), ide, jde,kde , igridtype) !Thiago 11/12/2014.
           !CALL E2A_Grid_Scalar( difv(1,1,1), im, jm, km, dat(1,1,1), ide, jde,
           !kde )

      IF ( trim(VarName)== 'DVBLDT' ) &
      CALL Destagger_Grid_Scalar( difv(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde , igridtype) !Thiago 11/12/2014.
           !CALL E2A_Grid_Scalar( difv(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde )

      IF ( trim(VarName)== 'KM' ) &
      CALL Destagger_Grid_Scalar( km_diffus(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde , igridtype) !Thiago 11/12/2014.
           !CALL E2A_Grid_Scalar( km_diffus(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde )

      IF ( trim(VarName)== 'KH' ) &
      CALL Destagger_Grid_Scalar( kh_diffus(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde , igridtype) !Thiago 11/12/2014.
           !CALL E2A_Grid_Scalar( kh_diffus(1,1,1), im, jm, km, dat(1,1,1), ide, jde, kde )

!     Horizontal diffusion (gopal's doing)

      IF ( trim(VarName)== 'DUHDIFFDT' ) &
      CALL Destagger_Grid_Scalar( hdifu(1,1,1), im, jm, km, dat(1,1,1), ide, jde,kde , igridtype) !Thiago 11/12/2014.
           !CALL E2A_Grid_Scalar( hdifu(1,1,1), im, jm, km, dat(1,1,1), ide, jde,
           !kde )

      IF ( trim(VarName)== 'DVHDIFFDT' ) &
      CALL Destagger_Grid_Scalar( hdifv(1,1,1), im, jm, km, dat(1,1,1), ide,jde,kde , igridtype) !Thiago 11/12/2014.
           !CALL E2A_Grid_Scalar( hdifv(1,1,1), im, jm, km, dat(1,1,1), ide, jde,
           !kde )

 
      DeAllocate( dat, wks ) 

      ENDIF !* for VarName - - - - - - - - - - - - - - - - - - - - - - - 
      enddo

      CALL ext_ncd_get_next_var( dh0, VarName, Status_next_var ) 

   ENDDO !<----------------------------------------------------------------- 
!
!* Special E-to-A-grid conversion for 10m winds 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!* Generalized routine includes E2A and B2A. This is Thiago's doings
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IF    ( fcst_time == 0 .or. igridtype==GRID_TYPE_ROTATEDLLB) THEN                        !Thiago 11/12/2014.
      CALL Destagger_Grid_Vector( u10(1,1), im, jm, 1, u10x(1,1), ide, jde, 1 , igridtype)  !Thiago 11/12/2014.
      CALL Destagger_Grid_Vector( v10(1,1), im, jm, 1, v10x(1,1), ide, jde, 1 , igridtype)  !Thiago 11/12/2014.
      !CALL E2A_Grid_Vector( u10(1,1), im, jm, 1, u10x(1,1), ide, jde, 1 ) 
      !CALL E2A_Grid_Vector( v10(1,1), im, jm, 1, v10x(1,1), ide, jde, 1 ) 
   ELSE !( fcst_time /= 0 )! 
      CALL E2A_Debug_Vector(u10(1,1), im, jm, 1, u10x(1,1), ide, jde, 1 )
      CALL E2A_Debug_Vector(v10(1,1), im, jm, 1, v10x(1,1), ide, jde, 1 )
   ENDIF

   DeAllocate( u10x, v10x ) 
   DeAllocate( txx, tyx )
   
!**************************************************! 
!  Coordinates & transformation of rotated A-grids ! 
!__________________________________________________! 

!* Rotated grid's rotated & standard coordinates (degree) 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   Allocate(  xrr(im,jm),  yrr(im,jm), xrs(im,jm), yrs(im,jm) ) 
   Allocate( wcos(im,jm), wsin(im,jm), urs(im,jm), vrs(im,jm), xrlen(jm) ) 
   
   IF ( xrs11 /= xrs31 ) THEN !* Earth latlon defined ------------------ 
      Call Convert_LatLon( 1, 1, 1, xrs11, yrs11, xxx, yyy, ctr_lon, ctr_lat ) 
      xrr(1,1) = nint(xxx/dlon)*dlon ; yrr(1,1) = nint(yyy/dlat)*dlat 
   ELSE  !* Earth latlon undefined, assume data on d01 ----------------- 
      xrr(1,1) = 0.5*(1-im)*dlon     ; yrr(1,1) = 0.5*(1-jm)*dlat 
   ENDIF !* ------------------------------------------------------------ 
   DO i = 1, im ; xrr(i,:) = xrr(1,1) + (i-1)*dlon ; ENDDO 
   DO j = 1, jm ; yrr(:,j) = yrr(1,1) + (j-1)*dlat ; ENDDO 

   Call Convert_LatLon( 2, im, jm, xrs, yrs, xrr, yrr, ctr_lon, ctr_lat ) 

!* Rotated grid's interval lengths (m) 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   yrlen = dlat*deg2rad*r_earth ; xxx = dlon*deg2rad*r_earth 
   DO j = 1, jm ; xrlen(j) = xxx*cos( yrr(1,j)*deg2rad ) ; ENDDO 
   
!* Matrix elements for converting vector components 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!  uur =   wcos*uus + wsin*vvs  !  uus = wcos*uur - wsin*vvr 
!  vvr = - wsin*uus + wcos*vvs  !  vvs = wsin*uur + wcos*vvr 
   
   olamda = ctr_lon*deg2rad 
   otheta = ctr_lat*deg2rad 
   aaa = cos(otheta) ; bbb = sin(otheta) 
   DO j = 1, jm ; ccc = 1./cos(yrr(1,j)*deg2rad)
   DO i = 1, im         !* poles not allowed ... 
      xxx = xrs(i,j)*deg2rad - olamda ; yyy = yrs(i,j)*deg2rad 
      wcos(i,j)= ( aaa*cos(yyy) + bbb*sin(yyy)*cos(xxx) )*ccc 
      wsin(i,j)= - bbb*sin(xxx)*ccc 
   ENDDO 
   ENDDO 

!************************************************! 
!  Derive diagnostic variables on sigma-p levels ! 
!________________________________________________! 

   IF ( derive_vars ) THEN  !* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   
!* Surface air temperature (K) 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   DO j = 1, jm ; DO i = 1, im 
      ts(i,j) = ths(i,j)*( .001*ps(i,j) )**kd_air 
   ENDDO ; ENDDO 
   
!* Cut off SST & SSQ over land 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ssq(:,:) = qs(:,:) 
   DO j = 1, jm ; DO i = 1, im 
      IF ( sm(i,j) < 1. ) THEN 
         sst(i,j) = undef ; ssq(i,j) = undef 
      ENDIF 
   ENDDO ; ENDDO 

!* Fix signs of surface heat fluxes 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   tw(:,:) = - tw(:,:) ; qw(:,:) = - qw(:,:) 

!* Geopotential height (m) 
!* ~~~~~~~~~~~~~~~~~~~~~~~ 

!  WRITE(0,*)'------------ GPM HEIGHTS AT THE DOMAIN CENTER----------'
   ccc = rd_air/g_earth 
   DO j = 1, jm 
   DO i = 1, im 
         aaa = ccc*(1+.608*qq(i,j,1))*tt(i,j,1)*dp(i,j,1)/pp(i,j,1) 
         gh(i,j,1) = zs(i,j) + aaa 
!         IF(i==im/2 .and. j==jm/2)WRITE(0,*)gh(i,j,1),pp(i,j,1),qq(i,j,1),tt(i,j,1),dp(i,j,1)
      DO k = 2, km 
         bbb = ccc*(1+.608*qq(i,j,k))*tt(i,j,k)*dp(i,j,k)/pp(i,j,k) 
         gh(i,j,k) = gh(i,j,k-1) + aaa + bbb 
	 aaa = bbb 
!         IF(i==im/2 .and.j==jm/2)WRITE(0,*)gh(i,j,k),pp(i,j,k),qq(i,j,k),tt(i,j,k),dp(i,j,k)
      ENDDO 
   ENDDO 
   ENDDO 
   
!...............................................................................
 
!* PBL wind, temperature & humidity 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   aaa = .0185/g_earth !* momentum roughness coeff over ocean 
   bbb = .520*.001*aaa !* adjusted momentum roughness for 30m wind 
   vcut = 0.1 ; define_u10 = .false. 

   IF ( define_u10 ) THEN !* re-define uf, z0, z0t, u10, v10 - - - - - - 
   
   eee = 4.7*karman*g_earth 
!  if ( fcst_time == 0. ) then 
   DO j = 1, jm ; DO i = 1, im 
      www = sqrt( uu(i,j,1)**2 + vv(i,j,1)**2 ) ; www = max(www,vcut) 
      zz0 = z0(i,j) ; IF ( sm(i,j) == 1.) zz0 = bbb*www**2 !* ocean roughness 
      ttt = tt(i,j,1)*( 1000./pp(i,j,1) )**kd_air ; zzz = gh(i,j,1) - zs(i,j) 
!     uf(i,j) = ( karman*www - eee*(ttt-ths(i,j))*zzz/(ttt*www) )/log(zzz/zz0) 
      uf(i,j) = ( karman*www )/log(zzz/zz0) 
   ENDDO ; ENDDO 
!  endif 
   
   z0t(:,:) = z0(:,:) 
   DO j = 1, jm ; DO i = 1, im 
   if ( sm(i,j) == 1. ) then 
      z0(i,j) = aaa*uf(i,j)**2  !* momentum (z0) & thermal (z0t) roughness *!
     z0t(i,j) = z0(i,j)/exp( 2.67*sqrt(sqrt(z0(i,j)*uf(i,j)/viscos)) - 2.57 ) 
   endif 
   ENDDO ; ENDDO 
   
!           *     *     *     *     *     *     *     *     *     * 
   
   ddd = 1.00 !* 1.27 ? 
   eee = - 3.*log(2.) + .5*pi_cnst 
   DO j = 1, jm ; DO i = 1, im 
      www = sqrt( uu(i,j,1)**2 + vv(i,j,1)**2 )
      www = max(www,vcut) ; uur = uu(i,j,1)/www ; vvr = vv(i,j,1)/www 
!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      rho =   ps(i,j)*100./( rd_air*(1+.608*qs(i,j))*ts(i,j) ) 
      zzz = - 10*g_earth*karman*tw(i,j)/(rho*cp_air)/(ths(i,j)*uf(i,j)**3) 
      IF ( zzz < 0. ) THEN 
         xxx =   sqrt( sqrt( 1. - 16.*zzz ) ) 
         z10 = 2.*log(1.+xxx) + log(1.+xxx**2) - 2.*atan(xxx) + eee 
      ELSEIF ( zzz <= .5 ) THEN 
         z10 = - 5.*zzz 
      ELSEIF ( zzz <= 10 ) THEN 
         z10 = - 7.*log(zzz) - 4.25/zzz + .5/zzz**2 - .852  
      ELSE ! ( zzz  > 10 ) 
         z10 = log(zzz) - .7585*zzz - 12.11 
      ENDIF 
      www = ddd*uf(i,j)/karman*( log(10./z0(i,j)) - z10 ) 
      u10(i,j) = www*uur  ;  v10(i,j) = www*vvr 
   ENDDO ; ENDDO 
   
   ENDIF  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
   
!* Pielke's PBL temperature & humidity  
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   aaa = 0.40 !* von Karman constant (k) 
   bbb = 1.05 !* mixing height (beta) 
   ccc = 1.00/(aaa*bbb) ; ddd = .5*pi_cnst 
   
   IF ( fcst_time > 0.) THEN 
   DO j = 1, jm ; DO i = 1, im 
      rho = 100.*ps(i,j)/( rd_air*(1+.608*qs(i,j))*ts(i,j) ) 
      xxx = tw(i,j)/(rho*cp_air) ; ttt = - ccc*xxx/uf(i,j) 
      yyy = qw(i,j)/(rho*heat0c) ; qqq = - ccc*yyy/uf(i,j) 
      rrr = - g_earth*aaa*xxx/(ths(i,j)*uf(i,j)**3)  !* 1/L
      IF ( rrr > 0. ) THEN 
         z02 = - 8.00*2.*rrr ; z10 = - 8.00*10*rrr 
      ELSE !* Pielke 2002 *! 
         z02 = 2.*log( .5*(1.+sqrt(1-11.6*2.*rrr)) ) 
         z10 = 2.*log( .5*(1.+sqrt(1-11.6*10*rrr)) ) 
      ENDIF 
      xxx = log(2./z0(i,j)) ; yyy = log(10/z0(i,j)) 
      q02(i,j) =  qs(i,j) + qqq*( xxx - z02 ) 
      q10(i,j) =  qs(i,j) + qqq*( yyy - z10 ) 
      t02(i,j) = ths(i,j) + ttt*( xxx - z02 ) 
      t10(i,j) = ths(i,j) + ttt*( yyy - z10 ) 
!     --------------------------------------- 
      zzz = gh(i,j,1) - zs(i,j) ; z02 = (zzz-2.)/zzz ; z10 = (zzz-10)/zzz 
      t02(i,j) = t02(i,j)*(.001*(z02*ps(i,j)+(1.-z02)*pp(i,j,1)))**kd_air 
      t10(i,j) = t10(i,j)*(.001*(z10*ps(i,j)+(1.-z10)*pp(i,j,1)))**kd_air 
   ENDDO ; ENDDO 
   ENDIF 
   
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   
!* Dewpoint temperature (K) on sigma-p levels 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   CALL Dewpoint_RH( pp, qq, tt, td, rh, im, jm, km ) 
   
!* Equivalent potential temperature (K) on sigma-p levels 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   tt0 = 273.15   !* ice-point temperature (K) 
   ee0 = 6.1121   !* ice-point vapor pressure (mb) 
   aa1 = 17.502   !*  factor  for water vapor pressure 
   bb1 = 240.97   !* constant for water vapor pressure 
   aa2 = 22.452   !*  factor  for  ice  vapor pressure 
   bb2 = 272.55   !* constant for  ice  vapor pressure 
   eps = .62198   !* ratio of vapor & air molecular weights
   tts_error = .1 !* error constraint for "tts" computation 
   
   iter_eqt = 0 
   DO k = 1, km ; DO j = 1, jm ; DO i = 1, im 
      ttt =     tt(i,j,k) 
      rrr =     qq(i,j,k)/( 1. - qq(i,j,k) ) 
      eee = rrr*pp(i,j,k)/( eps + rrr  ) 
      ykappam =  1./( kd_air*(1.-.245*rrr) ) 
      ees = eee 
      ttx = ttt 
      DO n = 0, itermax 
         IF ( ees == ee0 ) THEN 
            tts = tt0 
         ELSEIF ( ees > ee0 ) THEN 
            xxx = log(ees/ee0) 
            tts = bb1*xxx/(aa1-xxx) + tt0 
         ELSE  !( ees < ee0 ) 
            xxx = log(ees/ee0) 
            tts = bb2*xxx/(aa2-xxx) + tt0 
         ENDIF 
         IF ( tts > ttx - tts_error ) THEN 
            EXIT 
         ELSE 
            ees = eee*exp( log(tts/ttt)*ykappam ) 
            ttx = tts 
         ENDIF 
      ENDDO 
      iter_eqt = max( iter_eqt, n ) 
      aaa = ttt*( (1000.+1607.8*rrr)/pp(i,j,k) )**kd_air 
      bbb = aaa*exp( (3132.9-2.37*tts)*rrr/tts ) 
      qt(i,j,k) = bbb*exp( exp(1.62*log(bbb)+14.3*log(tts)-96.) ) 
   ENDDO ; ENDDO ; ENDDO 
   
!  PRINT*, '------------------------------------------------------------' 
!  PRINT*, 'Max Iterations for Eqv.Pot.Temp =', iter_eqt 
!  PRINT*, '------------------------------------------------------------' 
   
   ENDIF !* ( derive_vars ) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   
   DeAllocate( dp ) 
   
!* Forecast time and domain 
!* ~~~~~~~~~~~~~~~~~~~~~~~~ 
   n = nint(fcst_time) ; ifcst_hour = n/60 ; ifcst_minute = Mod(n,60) 
   
   WRITE(*,'(/, " Processing Domain: ",  I3)') domain_num 
   WRITE(*,'(   " Forecast Hour Min: ", 2I3)') ifcst_hour, ifcst_minute 
   
!******************************************************! 
!  a. Storm tracker by centroid of MSLP on moving nest ! 
!______________________________________________________! 
   
   IF ( domain_num == 1 ) THEN 
      pmin = MinVal( p0(:,:) )  !* (mb) 
      vmax = 1.94384449*MaxVal( sqrt(u10(:,:)**2+v10(:,:)**2) ) !* (knot) 
      DO j = 1, jm ; DO i = 1, im 
       if(pmin .eq. p0(i,j))then
         iloc=i;jloc=j;storm_lon = xrs(iloc,jloc); storm_lat = yrs(iloc,jloc)
       endif
      ENDDO;ENDDO
!     iloc = -999 ; jloc = -999 ; storm_lon = -999 ; storm_lat = -999 
   ENDIF 

   IF ( domain_num > 1 ) THEN !* Centroid Method ------------------- 
   
   ioff = nint(r_storm/dlon) ; ibd = 3 
   joff = nint(r_storm/dlat) ; jbd = 3 
   iloc = (im+1)/2  ;  jloc = (jm+1)/2 ; ierror = 0 
   storm_lon = xrs(iloc,jloc); storm_lat = yrs(iloc,jloc) 

   DO n = 0, itermax  !* - - - - - - - - - - - - - - - - - - - - - - 
   
   i1 = iloc - ioff ; i1 = max( 1+ibd, i1  ) 
   i2 = iloc + ioff ; i2 = min( i2, im-ibd ) 
   j1 = jloc - joff ; j1 = max( 1+jbd, j1  ) 
   j2 = jloc + joff ; j2 = min( j2, jm-jbd ) 
   
   pmax = MaxVal( p0(i1:i2,j1:j2) ) !* (mb) 
   pmin = MinVal( p0(i1:i2,j1:j2) ) !* (mb) 
   pcut = .5*( pmax + pmin ) 
 
   aaa = 0 ; bbb = 0 ; ccc = 0 ; ddd = 0 ; eee = 0 
   DO j = j1, j2 ; DO i = i1, i2 
      www = max( pcut - p0(i,j), 0. )
      aaa = aaa + www*i ; bbb = bbb + www*j ; eee = eee + www 
      ccc = ccc + www*xrs(i,j) ; ddd = ddd + www*yrs(i,j) 
   ENDDO ; ENDDO
   
   IF ( eee == 0. ) EXIT 
   iold = iloc ; jold = jloc 
   iloc = nint(aaa/eee); jloc = nint(bbb/eee) 
   storm_lon = ccc/eee ; storm_lat = ddd/eee 
   IF ( (iloc-iold)**2 + (jloc-jold)**2 <= ierror**2 ) EXIT 
   
   ENDDO  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   iter_cen = n 
   
   WRITE(*,'(" MSLP Iterations for Storm Tracking =", I9)'), iter_cen 
   WRITE(*,'(" Storm Center by the MSLP: lon, lat =", 2F12.5)') &
               storm_lon, storm_lat 

!   vmax = 1.94384449*MaxVal( sqrt(u10(:,:)**2+v10(:,:)**2) ) !* (knot) 

   vmax = MaxVal( sqrt(u10(:,:)**2+v10(:,:)**2) )

   ENDIF !* -------------------------------------------------------- 

   WRITE(0,*)'iloc=',iloc,'jloc=',jloc,'pmin=',pmin     !p0(iloc,jloc)

   WRITE(52,1009) &         ! gopal's change
   ifcst_hour,storm_lon,storm_lat,pmin,vmax 

1009 FORMAT(I5,',',F11.4,',',F10.4,',',F11.4,',',F12.4)    !,',',F12.4)

   lat10 = nint(storm_lat*10) ; lon10 = nint(-storm_lon*10) 

!
!       recompute model storm location based on pmin and store as iloc1
!       and iloc2 (earlier computation is lost in iteration)
!
        do jj=1,jm; do ii=1,im
         if(p0(ii,jj) .eq. pmin)then; iloc1=ii; jloc1=jj; endif
        enddo; enddo
!
        write(0,*)'NEW iloc and jloc',iloc1,jloc1
!

   WRITE(54,'(A2,A1,A3,A1,I11,A1,I3.2,A1,A5,A1,I4,A1,I4,A2,I5,A2,I4,A1,I5)') & 
        'AL',',', storm_num,',', start_date,',', 99,',', ' HWRX',',',  & 
        ifcst_hour,',', lat10,'N,', lon10,'W,', nint(vmax),',', nint(pmin) 
   

   WRITE(51,*) &
  '      Date_Time      i_a  j_a  Storm_Lon  Storm_Lat  10m_V(kt)   MSLP(mb)',&
  '-------------------  ---  ---  ---------  ---------  ---------  ---------' 
   
   write(53,*) '                                                             '
   write(53,*) 'Column 1. Storm basin: AL = Atlantic, EP = Eastern Pacific.  '
   write(53,*) 'Column 2. Annual cyclone (storm) number: 1-99.               '
   write(53,*) 'Column 3. Model start date: YYYYMMDDHH.                      '
   write(53,*) 'Column 4. Technique sorting number: 00-99.                   '
   write(53,*) 'Column 5. Technique: HWRX = Experimental HWRF (HWRFx).       '
   write(53,*) 'Column 6. Forecast   period    (hour).                       '
   write(53,*) 'Column 7. Latitude  (tenths of degree).                      '
   write(53,*) 'Column 8. Longitude (tenths of degree).                      '
   write(53,*) 'Column 9. Maximum sustained 10m wind speed (knot).           '
   write(53,*) 'Column10. Minimum sea-level pressure (mb).                   '
   write(53,*) '                                                             '
   
   CLOSE(51) ; CLOSE(52) ; CLOSE(53) ; CLOSE(54) 

!   
!* Output data grid information 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   i0  = (im+1)/2   ; j0  = (jm+1)/2   ; xxx = sqrt( dlon**2 + dlat**2 ) 
   aaa = xrs(i0,j0) ; bbb = yrs(i0,j0) ; ccc = ctr_lon ; ddd = ctr_lat 
   imr = im ; jmr = jm ; dxx = dlon ; dyy = dlat 
   IF ( dlon > refine )  THEN 
      imr = 3*im-2 ; jmr = 3*jm-2 ; dxx = dlon/3. ; dyy = dlat/3. 
   ENDIF 
   WRITE(*,'(" Processed  Domain Center: lon, lat =", 2F12.5)') aaa, bbb 
   WRITE(*,'(" Origin of  Rotated  Grid: lon, lat =", 2F12.5)') ccc, ddd 
   WRITE(*,'(" Native E-Grid Resolution: lon, lat =", 2F12.5)') xxx, xxx 
   WRITE(*,'(" Output A-Grid Resolution: lon, lat =", 2F12.5)') dxx, dyy 
   WRITE(*,'(" Output A-Grid Dimensions:  X, Y, Z =", 3I10)')   imr, jmr, km 
   
!******************************************************! 
!  b. Embed moving nest in outer domain for wind swath ! 
!______________________________________________________! 

   IF ( wind_swath ) THEN !* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!* Dimensions of canvas for moving nest composite 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   imx = (nint(d01cut/dlon)/2)*2+1 ; jmx = (nint(aspect*imx)/2)*2+1 
   Allocate( u10x(imx,jmx), v10x(imx,jmx), p0x(imx,jmx) ) 

   IF ( debug_data ) THEN !* --------------------------------------- 

!* Embed moving nest in extended rotated grid 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   u10x(:,:) = undef ; v10x(:,:) = undef ; p0x(:,:) = undef 
   i0 = (imx-1)/2 + nint(xrr(1,1)/dlon) 
   j0 = (jmx-1)/2 + nint(yrr(1,1)/dlat) 

   IF ( i0+im > imx .OR. j0+jm > jmx ) THEN !* - - - - - - - - - - - 
   
   PRINT*, '*ATTN: Outer domain is too small for wind swath on rotated  grid !' 
   
   ELSE !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   
   DO j = 1, jm ; j1 = j0 + j 
   DO i = 1, im ; i1 = i0 + i 
       p0x(i1,j1) =  p0(i,j) 
      u10x(i1,j1) = u10(i,j) 
      v10x(i1,j1) = v10(i,j) 
   ENDDO 
   ENDDO 

   ENDIF  !* ( ERROR ) - - - - - - - - - - - - - - - - - - - - - - - 
   
   WRITE(62) p0x ; WRITE(62) u10x ; WRITE(62) v10x 
   nvars =  3 

!* GrADS control file for moving nest on rotated grid 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   xxr1 = .5*(1-imx)*dlon ; xxr_ = dlon  !* Rotated longitudes 
   yyr1 = .5*(1-jmx)*dlat ; yyr_ = dlat  !* Rotated latitudes 
   
   WRITE(61,'("                                                      ")') 
   WRITE(61,'(" DSET    swath_rot.dat                                ")') 
   WRITE(61,'(" TITLE   Moving nest data on rotated latlon grid      ")') 
   WRITE(61,'(" OPTIONS sequential big_endian                        ")') 
   WRITE(61,'(" UNDEF",E11.2)') undef 
   WRITE(61,'("                                                      ")') 
   WRITE(61,'(" XDEF", I7,"   LINEAR", 2F11.5 )') imx, xxr1, xxr_ 
   WRITE(61,'(" YDEF", I7,"   LINEAR", 2F11.5 )') jmx, yyr1, yyr_ 
   WRITE(61,'(" ZDEF", I7,"   LINEAR          1          1 ")')   1 
   WRITE(61,'(" TDEF", I7,"   LINEAR   00:00Z01JAN0001  1yr")') ntime 
   WRITE(61,'("                                                      ")') 
   WRITE(61,'(" VARS", I3)') nvars
   WRITE(61,'(" p0        0   0    mean-sea-level pressure (mb)      ")') 
   WRITE(61,'(" u10       0   0    rotated 10m x-wind (m/s)          ")') 
   WRITE(61,'(" v10       0   0    rotated 10m y-wind (m/s)          ")') 
   WRITE(61,'(" ENDVARS                                              ")') 
   WRITE(61,'("                                                      ")') 
   
   ENDIF !* ( debug_data ) ----------------------------------------- 
   
!           *     *     *     *     *     *     *     *     *     * 
   
!* Standard grid that circumscribes the moving rotated grid 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   aaa = MinVal( xrs( 1,1:jm) )        ; ccc = MinVal( yrs(1:im,1 ) ) 
   bbb = MaxVal( xrs(im,1:jm) )        ; ddd = MaxVal( yrs(1:im,jm) ) 
   ims = (nint((bbb-aaa)/dlon)/2)*2+1  ; jms = (nint((ddd-ccc)/dlat)/2)*2+1
   Allocate( xsr(ims,jms), ysr(ims,jms), xss(ims,jms), yss(ims,jms) ) 
   
   xss(1,1) = ctr_lon + ( nint((.5*(aaa+bbb)-ctr_lon)/dlon) - (ims-1)/2 )*dlon 
   yss(1,1) = ctr_lat + ( nint((.5*(ccc+ddd)-ctr_lat)/dlat) - (jms-1)/2 )*dlat 
   DO i = 1, ims ; xss(i,:) = xss(1,1) + (i-1)*dlon ; ENDDO 
   DO j = 1, jms ; yss(:,j) = yss(1,1) + (j-1)*dlat ; ENDDO 

   Call Convert_LatLon( 1, ims, jms, xss, yss, xsr, ysr, ctr_lon, ctr_lat ) 
   
!* Convert rotated wind components to standard ones 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   urs(:,:) = wcos(:,:)*u10(:,:) - wsin(:,:)*v10(:,:) 
   vrs(:,:) = wsin(:,:)*u10(:,:) + wcos(:,:)*v10(:,:) 
   
!* Embed moving nest in extended standard grid 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   u10x(:,:) = undef ; v10x(:,:) = undef ; p0x(:,:) = undef 
   xrr1 = xrr(1,1) ; xrr2 = xrr(2,1) ; xrr8 = xrr(im-1,1) 
   yrr1 = yrr(1,1) ; yrr2 = yrr(1,2) ; yrr8 = yrr(1,jm-1) 
    i0 = (imx-1)/2 + nint( (xss(1,1)-ctr_lon)/dlon ) 
    j0 = (jmx-1)/2 + nint( (yss(1,1)-ctr_lat)/dlat ) 
   
   IF ( i0+ims > imx .OR. j0+jms > jmx ) THEN !* - - - - - - - - - - 
   
   PRINT*, '*ATTN: Outer domain is too small for wind swath on standard grid !' 
   PRINT*, '...... Wind swath data is *NOT* defined for this time level ......' 
   
   ELSEIF ( iorder == 2 ) THEN  !* - - - - - - - - - - - - - - - - - 

   DO j = 1, jms ; j1 = j0 + j 
   DO i = 1, ims ; i1 = i0 + i 
      xxx = xsr(i,j) ;  yyy = ysr(i,j) 
      IF ( yrr2 <= yyy .AND. yyy < yrr8 ) THEN 
      if ( xrr2 <= xxx .AND. xxx < xrr8 ) then 
             m  = 1 + int((xxx-xrr1)/dlon) ;  m1 = m + 1 
             n  = 1 + int((yyy-yrr1)/dlat) ;  n1 = n + 1 
            aaa = ( xxx - xrr(m,n) )/dlon  ; bbb = 1 - aaa 
            ccc = ( yyy - yrr(m,n) )/dlat  ; ddd = 1 - ccc 
          p0x(i1,j1) = bbb*( ddd* p0(m ,n) + ccc* p0(m ,n1) ) &
                     + aaa*( ddd* p0(m1,n) + ccc* p0(m1,n1) ) 
         u10x(i1,j1) = bbb*( ddd*urs(m ,n) + ccc*urs(m ,n1) ) &
                     + aaa*( ddd*urs(m1,n) + ccc*urs(m1,n1) ) 
         v10x(i1,j1) = bbb*( ddd*vrs(m ,n) + ccc*vrs(m ,n1) ) &
                     + aaa*( ddd*vrs(m1,n) + ccc*vrs(m1,n1) ) 
      endif 
      ENDIF 
   ENDDO 
   ENDDO 

   ELSEIF ( iorder == 4 ) THEN  !* - - - - - - - - - - - - - - - - - 

   DO j = 1, jms ; j1 = j0 + j 
   DO i = 1, ims ; i1 = i0 + i 
      xxx = xsr(i,j) ;  yyy = ysr(i,j) 
      IF ( yrr2 <= yyy .AND. yyy < yrr8 ) THEN 
      if ( xrr2 <= xxx .AND. xxx < xrr8 ) then 
            m2 =  1 + int((xxx-xrr1)/dlon);    m1 = m2-1 ; m4 = m2+2 
            n2 =  1 + int((yyy-yrr1)/dlat);    n1 = n2-1 ; n4 = n2+2 
           aaa =  (xxx - xrr(m2,n2))/dlon ;   bbb = (yyy - yrr(m2,n2))/dlat 
         wtx(1)= -(aaa)*(aaa-1)*(aaa-2)/6 ; wtx(2)= (aaa+1)*(aaa-1)*(aaa-2)/2
         wtx(3)= -(aaa+1)*(aaa)*(aaa-2)/2 ; wtx(4)= (aaa+1)*(aaa)*(aaa-1)/6
         wty(1)= -(bbb)*(bbb-1)*(bbb-2)/6 ; wty(2)= (bbb+1)*(bbb-1)*(bbb-2)/2
         wty(3)= -(bbb+1)*(bbb)*(bbb-2)/2 ; wty(4)= (bbb+1)*(bbb)*(bbb-1)/6
!        - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
         ff= p0(m1:m4,n1:n4); CALL Cubic_Interpo( p0x(i1,j1),ff,wtx,wty,undef) 
         ff=urs(m1:m4,n1:n4); CALL Cubic_Interpo(u10x(i1,j1),ff,wtx,wty,undef) 
         ff=vrs(m1:m4,n1:n4); CALL Cubic_Interpo(v10x(i1,j1),ff,wtx,wty,undef) 
      endif 
      ENDIF 
   ENDDO 
   ENDDO 

   ENDIF !* ( iorder ) - - - - - - - - - - - - - - - - - - - - - - - 
   
   WRITE(64) p0x ; WRITE(64) u10x ; WRITE(64) v10x 
   nvars = 3 

!* GrADS control file for moving nest on standard grid 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   xxs1 = .5*(1-imx)*dlon + ctr_lon ; xxs_ = dlon !* d01 west  bdry 
   yys1 = .5*(1-jmx)*dlat + ctr_lat ; yys_ = dlat !* d01 south bdry 
   
   WRITE(63,'("                                                      ")') 
   WRITE(63,'(" DSET    swath_std.dat                                ")') 
   WRITE(63,'(" TITLE   Moving nest data on standard latlon grid     ")') 
   WRITE(63,'(" OPTIONS sequential big_endian                        ")') 
   WRITE(63,'(" UNDEF",E11.2)') undef 
   WRITE(63,'("                                                      ")') 
   WRITE(63,'(" XDEF", I7,"   LINEAR", 2F11.5 )') imx, xxs1, xxs_ 
   WRITE(63,'(" YDEF", I7,"   LINEAR", 2F11.5 )') jmx, yys1, yys_ 
   WRITE(63,'(" ZDEF", I7,"   LINEAR          1          1 ")')   1 
   WRITE(63,'(" TDEF", I7,"   LINEAR   00:00Z01JAN0001  1yr")') ntime 
   WRITE(63,'("                                                      ")') 
   WRITE(63,'(" VARS", I3)') nvars
   WRITE(63,'(" p0        0   0    mean-sea-level pressure (mb)      ")') 
   WRITE(63,'(" u10       0   0    10m zonal wind comp (m/s)         ")') 
   WRITE(63,'(" v10       0   0    10m meridional wind (m/s)         ")') 
   WRITE(63,'(" ENDVARS                                              ")') 
   WRITE(63,'("                                                      ")') 
   
   CLOSE(61) ; CLOSE(62) ; CLOSE(63) ; CLOSE(64) 
   DeAllocate( u10x, v10x, p0x, xsr, ysr, xss, yss ) 
   
   ENDIF !* ( wind_swath ) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!******************************************************! 
!  1. Output sigma-p-level variables on rotated A-grid ! 
!______________________________________________________! 

   IF ( sigma_level ) THEN !* ------------------------------------------ 
   
!* Output variables on sigma-p levels 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   IF ( dlon <= refine ) THEN !* - - - - - - - - - - - - - - - - - - 

   WRITE(100) xrs ; WRITE(100) yrs 
   WRITE(100) zs  ; WRITE(100) ps 
   WRITE(100) ts  ; WRITE(100) qs 
   WRITE(100) p0  ; WRITE(100) prc
   WRITE(100) u10 ; WRITE(100) v10 
   WRITE(100) t10 ; WRITE(100) q10 
   WRITE(100) sst ; WRITE(100) ssq
   WRITE(100) tw  ; WRITE(100) qw 
   WRITE(100) uf  ; WRITE(100) z0  ;  WRITE(100) z0t 
   do k =  1, km  ; WRITE(100) uu(:,:,k) ; enddo 
   do k =  1, km  ; WRITE(100) vv(:,:,k) ; enddo 
   do k =  1, km  ; WRITE(100) ww(:,:,k) ; enddo 
   do k =  1, km  ; WRITE(100) tt(:,:,k) ; enddo 
   do k =  1, km  ; WRITE(100) td(:,:,k) ; enddo 
   do k =  1, km  ; WRITE(100) qt(:,:,k) ; enddo 
   do k =  1, km  ; WRITE(100) gh(:,:,k) ; enddo 
   do k =  1, km  ; WRITE(100) pp(:,:,k) ; enddo 
   do k =  1, km  ; WRITE(100) qq(:,:,k) ; enddo 
   do k =  1, km  ; WRITE(100) cw(:,:,k) ; enddo 
   do k =  1, km  ; WRITE(100) rh(:,:,k) ; enddo 
   
   nvars = 30 ; dxx = dlon ; dyy = dlat 
    imr  = im ; jmr =  jm  ; ims =  im ; jms = (nint(aspect*ims)/2)*2 + 1 

! gopal's doing for synthetic vortex
!
!   WRITE(0,*)'------------ SYNTHETIC OBS ---------------'
!   Do j=jloc-1,jloc+1
!    Do i=iloc-1,iloc+1
!      WRITE(0,*)'ILOC= ',ILOC,'JLOC= ',JLOC
!      WRITE(30,*)'Start Date= ',start_date,' forecast hour= ',ifcst_hour
!      WRITE(30,*)'LAT= ',xrs(i,j),' LON= ',yrs(i,j)
!      Do k=1,km
!       WRITE(30,1001)gh(i,j,k),pp(i,j,k),uu(i,j,k),vv(i,j,k),tt(i,j,k),td(i,j,k)
!      ENDDO
!     ENDDO
!   ENDDO
!
!1001 FORMAT(6F10.4)
!

   
!* Refine the grid for low resolutions 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   ELSE !* ( dlon > refine ) - - - - - - - - - - - - - - - - - - - - 
   
   im3 = 3*im - 2 ; jm3 = 3*jm - 2 ; k2d = 19 
   Allocate( uu3(im3,jm3,km), vv3(im3,jm3,km), ww3(im3,jm3,km) ) 
   Allocate( tt3(im3,jm3,km), td3(im3,jm3,km), qt3(im3,jm3,km) ) 
   Allocate( pp3(im3,jm3,km), qq3(im3,jm3,km), cw3(im3,jm3,km)) 
   Allocate( gh3(im3,jm3,km), rh3(im3,jm3,km), wks(im3,jm3,k2d))
   Allocate( f_ice3(im3,jm3,km),f_rain3(im3,jm3,km),f_rimef3(im3,jm3,km))
   Allocate(                                        mc_heat3(im3,jm3,km)) 
!
   Allocate(  difu3(im3,jm3,km),  difv3(im3,jm3,km), km_diffus3(im3,jm3,km), kh_diffus3(im3,jm3,km))   ! gopal's doing for the vertical diffusion term
   Allocate( hdifu3(im3,jm3,km), hdifv3(im3,jm3,km))                                                 ! gopal's doing for the horizontal diffusion term
   
   CALL Refine_Grid3x( xrs(1,1), wks(1,1,1),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( yrs(1,1), wks(1,1,2),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  zs(1,1), wks(1,1,3),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  ps(1,1), wks(1,1,4),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  ts(1,1), wks(1,1,5),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  qs(1,1), wks(1,1,6),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  p0(1,1), wks(1,1,7),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( prc(1,1), wks(1,1,8),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( u10(1,1), wks(1,1,9),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( v10(1,1), wks(1,1,10), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( t10(1,1), wks(1,1,11), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( q10(1,1), wks(1,1,12), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( sst(1,1), wks(1,1,13), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( ssq(1,1), wks(1,1,14), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  tw(1,1), wks(1,1,15), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  qw(1,1), wks(1,1,16), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  uf(1,1), wks(1,1,17), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  z0(1,1), wks(1,1,18), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( z0t(1,1), wks(1,1,19), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( uu, uu3, undef, im3, jm3, im, jm, km )
   CALL Refine_Grid3x( vv, vv3, undef, im3, jm3, im, jm, km )
   CALL Refine_Grid3x( ww, ww3, undef, im3, jm3, im, jm, km )
   CALL Refine_Grid3x( tt, tt3, undef, im3, jm3, im, jm, km )
   CALL Refine_Grid3x( td, td3, undef, im3, jm3, im, jm, km )
   CALL Refine_Grid3x( qt, qt3, undef, im3, jm3, im, jm, km )
   CALL Refine_Grid3x( gh, gh3, undef, im3, jm3, im, jm, km )
   CALL Refine_Grid3x( pp, pp3, undef, im3, jm3, im, jm, km )
   CALL Refine_Grid3x( qq, qq3, undef, im3, jm3, im, jm, km )
   CALL Refine_Grid3x( cw, cw3, undef, im3, jm3, im, jm, km )
   CALL Refine_Grid3x( rh, rh3, undef, im3, jm3, im, jm, km )
   CALL Refine_Grid3x( f_ice, f_ice3, undef, im3, jm3, im, jm, km )
   CALL Refine_Grid3x( f_rain, f_rain3, undef, im3, jm3, im, jm, km )
   CALL Refine_Grid3x( f_rimef, f_rimef3, undef, im3, jm3, im, jm, km )
   CALL Refine_Grid3x( mc_heat, mc_heat3, undef, im3, jm3, im, jm, km )

!* For diffusion terms (gopal's doing)
!*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   CALL Refine_Grid3x( difu, difu3, undef, im3, jm3, im, jm, km )
   CALL Refine_Grid3x( difv, difv3, undef, im3, jm3, im, jm, km )
   CALL Refine_Grid3x( km_diffus, km_diffus3, undef, im3, jm3, im, jm, km )  
   CALL Refine_Grid3x( kh_diffus, kh_diffus3, undef, im3, jm3, im, jm, km )
   CALL Refine_Grid3x( hdifu, hdifu3, undef, im3, jm3, im, jm, km )
   CALL Refine_Grid3x( hdifv, hdifv3, undef, im3, jm3, im, jm, km )
 
!* Positive-definiteness constraint 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   wks(:,:,6)  = max( wks(:,:,6),  0. )
   wks(:,:,8)  = max( wks(:,:,8),  0. )
   wks(:,:,12) = max( wks(:,:,12), 0. )
   wks(:,:,14) = max( wks(:,:,14), 0. )
   wks(:,:,17) = max( wks(:,:,17), 0. )
   wks(:,:,18) = max( wks(:,:,18), 0. )
   wks(:,:,19) = max( wks(:,:,19), 0. )
   cw3(:,:,:)  = max( cw3(:,:,:),  0. )
   rh3(:,:,:)  = min( rh3(:,:,:), 100.)
   td3(:,:,:)  = min( td3(:,:,:), tt3(:,:,:) )
   f_ice3(:,:,:)  = max( f_ice3(:,:,:),  0. )
   f_rain3(:,:,:)  = max( f_rain3(:,:,:),  0. )
   f_rimef3(:,:,:)  = max( f_rimef3(:,:,:),  0. )
   mc_heat3(:,:,:)  = max( mc_heat3(:,:,:),  0. )

!* Diffusion terms (gopal's doing)
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   difu3(:,:,:)  = max( difu3(:,:,:),  0. )
   difv3(:,:,:)  = max( difv3(:,:,:),  0. )
   km_diffus3(:,:,:)  = max( km_diffus3(:,:,:),  0. )
   kh_diffus3(:,:,:)  = max( kh_diffus3(:,:,:),  0. )
   hdifu3(:,:,:)  = max( hdifu3(:,:,:),  0. )
   hdifv3(:,:,:)  = max( hdifv3(:,:,:),  0. )
   
   do k = 1, k2d ; WRITE(100) wks(:,:,k) ; enddo 
   do k = 1, km  ; WRITE(100) uu3(:,:,k) ; enddo 
   do k = 1, km  ; WRITE(100) vv3(:,:,k) ; enddo 
   do k = 1, km  ; WRITE(100) ww3(:,:,k) ; enddo 
   do k = 1, km  ; WRITE(100) tt3(:,:,k) ; enddo 
   do k = 1, km  ; WRITE(100) td3(:,:,k) ; enddo 
   do k = 1, km  ; WRITE(100) qt3(:,:,k) ; enddo 
   do k = 1, km  ; WRITE(100) gh3(:,:,k) ; enddo 
   do k = 1, km  ; WRITE(100) pp3(:,:,k) ; enddo 
   do k = 1, km  ; WRITE(100) qq3(:,:,k) ; enddo 
   do k = 1, km  ; WRITE(100) cw3(:,:,k) ; enddo 
   do k = 1, km  ; WRITE(100) rh3(:,:,k) ; enddo 
   
   nvars = k2d + 11 ;  dxx = dlon/3. ; dyy = dlat/3. 
   imr = im3 ; jmr = jm3 ; ims = im3 ; jms = (nint(aspect*ims)/2)*2 + 1 
   
   DeAllocate( uu3, vv3, ww3, tt3, td3, qt3, gh3, pp3, qq3, cw3, rh3, wks) 
   DeAllocate( f_ice3,f_rain3,f_rimef3)
   DeAllocate(                mc_heat3)

   DeAllocate(  difu3,  difv3, km_diffus3,kh_diffus3)   ! gopal's doing for the diffusion term
   DeAllocate( hdifu3, hdifv3)
   
   ENDIF !* ( refine ) - - - - - - - - - - - - - - - - - - - - - - - 

!           *     *     *     *     *     *     *     *     *     * 
   
!* GrADS control file for sigma-p-level data on rotated grid 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   xxr1 = .5*(1-imr)*dxx ; xxr_ = dxx !* Lagrangian coord for moving nest 
   yyr1 = .5*(1-jmr)*dyy ; yyr_ = dyy !* and for fixed outer domain (d01) 
   
   WRITE(101,'("                                                      ")') 
   WRITE(101,'(" DSET    sigma_lev.dat                                ")') 
   WRITE(101,'(" TITLE   Sigma-p-level data on rotated latlon grid    ")') 
   WRITE(101,'(" OPTIONS sequential big_endian                        ")') 
   WRITE(101,'(" UNDEF",E11.2)') undef 
   WRITE(101,'("                                                      ")') 
   WRITE(101,'(" XDEF", I7,"   LINEAR", 2F11.5 )') imr, xxr1, xxr_ 
   WRITE(101,'(" YDEF", I7,"   LINEAR", 2F11.5 )') jmr, yyr1, yyr_ 
   WRITE(101,'(" ZDEF", I7,"   LINEAR          1          1 ")')  km 
   WRITE(101,'(" TDEF", I7,"   LINEAR   00:00Z01JAN0001  1yr")') ntime 
   WRITE(101,'("                                                      ")') 
   WRITE(101,'(" VARS", I3)') nvars
   WRITE(101,'(" xx        0   0    Earth longitude (deg)             ")') 
   WRITE(101,'(" yy        0   0    Earth latitude  (deg)             ")') 
   WRITE(101,'(" zs        0   0    Earth topography (m)              ")') 
   WRITE(101,'(" ps        0   0    surface air   pressure (mb)       ")') 
   WRITE(101,'(" ts        0   0    surface air temperature (K)       ")') 
   WRITE(101,'(" qs        0   0    surface specific humidity (kg/kg) ")') 
   WRITE(101,'(" p0        0   0    mean-sea-level pressure (mb)      ")') 
   WRITE(101,'(" prc       0   0    precipitation (m)                 ")') 
   WRITE(101,'(" u10       0   0    rotated 10m x-wind (m/s)          ")') 
   WRITE(101,'(" v10       0   0    rotated 10m y-wind (m/s)          ")') 
   WRITE(101,'(" t10       0   0    10m air temperature (K)           ")') 
   WRITE(101,'(" q10       0   0    10m specific humidity (kg/kg)     ")') 
   WRITE(101,'(" sst       0   0    sea-surface temperature (K)       ")') 
   WRITE(101,'(" ssq       0   0    sea-sfc specific humidity (kg/kg) ")') 
   WRITE(101,'(" tw        0   0    sensible heat flux (W/m^2)        ")') 
   WRITE(101,'(" qw        0   0    latent   heat flux (W/m^2)        ")') 
   WRITE(101,'(" uf        0   0    friction velocity  (m/s)          ")') 
   WRITE(101,'(" z0        0   0    momentum roughness height (m)     ")') 
   WRITE(101,'(" z0t       0   0    thermal  roughness height (m)     ")') 
   WRITE(101,'(" uu  ", I7,"   0    rotated x-wind (m/s)              ")') km 
   WRITE(101,'(" vv  ", I7,"   0    rotated y-wind (m/s)              ")') km 
   WRITE(101,'(" ww  ", I7,"   0    vertical  wind (m/s)              ")') km 
   WRITE(101,'(" tt  ", I7,"   0    air temperature (K)               ")') km 
   WRITE(101,'(" td  ", I7,"   0    dewpoint temperature (K)          ")') km 
   WRITE(101,'(" qt  ", I7,"   0    equiv. pot. temperature (K)       ")') km 
   WRITE(101,'(" gh  ", I7,"   0    geopotential height (m)           ")') km 
   WRITE(101,'(" pp  ", I7,"   0    air  pressure  (mb)               ")') km 
   WRITE(101,'(" qq  ", I7,"   0    specific humidity (kg/kg)         ")') km 
   WRITE(101,'(" cw  ", I7,"   0    total condensate  (kg/kg)         ")') km 
   WRITE(101,'(" rh  ", I7,"   0    relative humidity (%)             ")') km 
   WRITE(101,'(" ENDVARS                                              ")') 
   WRITE(101,'("                                                      ")') 
   
!* GrADS control file for sigma-level data on standard grid 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   xxs1 = .5*(1-ims)*dxx + ctr_lon ; xxs_ = dxx !* d01 west  bdry 
   yys1 = .5*(1-jms)*dyy + ctr_lat ; yys_ = dyy !* d01 south bdry 
   
   WRITE(102,'("                                                      ")') 
   WRITE(102,'(" DSET    sigma_lev.dat                                ")') 
   WRITE(102,'(" TITLE   Sigma-p-level data on standard latlon grid   ")') 
   WRITE(102,'(" OPTIONS sequential big_endian                        ")') 
   WRITE(102,'(" UNDEF",E11.2)') undef 
   WRITE(102,'("                                                      ")') 
   WRITE(102,'(" PDEF", I7,  I6,   A10,   4F11.5)') & 
                       imr, jmr, 'ETA.U', ctr_lon, ctr_lat, dxx,  dyy 
   WRITE(102,'(" XDEF", I7,"   LINEAR", 2F11.5 )') ims, xxs1, xxs_ 
   WRITE(102,'(" YDEF", I7,"   LINEAR", 2F11.5 )') jms, yys1, yys_ 
   WRITE(102,'(" ZDEF", I7,"   LINEAR          1          1 ")')  km 
   WRITE(102,'(" TDEF", I7,"   LINEAR   00:00Z01JAN0001  1yr")') ntime 
   WRITE(102,'("                                                      ")') 
   WRITE(102,'(" VECTORPAIRS   u10,v10  uu,vv                         ")') 
   WRITE(102,'(" VARS", I3)') nvars
   WRITE(102,'(" xx        0   0    Earth longitude (deg)             ")') 
   WRITE(102,'(" yy        0   0    Earth latitude  (deg)             ")') 
   WRITE(102,'(" zs        0   0    Earth topography (m)              ")') 
   WRITE(102,'(" ps        0   0    surface air   pressure (mb)       ")') 
   WRITE(102,'(" ts        0   0    surface air temperature (K)       ")') 
   WRITE(102,'(" qs        0   0    surface specific humidity (kg/kg) ")') 
   WRITE(102,'(" p0        0   0    mean-sea-level pressure (mb)      ")') 
   WRITE(102,'(" prc       0   0    precipitation (m)                 ")') 
   WRITE(102,'(" u10       0   0    10m zonal wind comp (m/s)         ")') 
   WRITE(102,'(" v10       0   0    10m meridional wind (m/s)         ")') 
   WRITE(102,'(" t10       0   0    10m air temperature (K)           ")') 
   WRITE(102,'(" q10       0   0    10m specific humidity (kg/kg)     ")') 
   WRITE(102,'(" sst       0   0    sea-surface temperature (K)       ")') 
   WRITE(102,'(" ssq       0   0    sea-sfc specific humidity (kg/kg) ")') 
   WRITE(102,'(" tw        0   0    sensible heat flux (W/m^2)        ")') 
   WRITE(102,'(" qw        0   0    latent   heat flux (W/m^2)        ")') 
   WRITE(102,'(" uf        0   0    friction velocity  (m/s)          ")') 
   WRITE(102,'(" z0        0   0    momentum roughness height (m)     ")') 
   WRITE(102,'(" z0t       0   0    thermal  roughness height (m)     ")') 
   WRITE(102,'(" uu  ", I7,"   0    zonal wind comp (m/s)             ")') km 
   WRITE(102,'(" vv  ", I7,"   0    meridional wind (m/s)             ")') km 
   WRITE(102,'(" ww  ", I7,"   0    vertical  wind (m/s)              ")') km 
   WRITE(102,'(" tt  ", I7,"   0    air temperature (K)               ")') km 
   WRITE(102,'(" td  ", I7,"   0    dewpoint temperature (K)          ")') km 
   WRITE(102,'(" qt  ", I7,"   0    equiv. pot. temperature (K)       ")') km 
   WRITE(102,'(" gh  ", I7,"   0    geopotential height (m)           ")') km 
   WRITE(102,'(" pp  ", I7,"   0    air  pressure  (mb)               ")') km 
   WRITE(102,'(" qq  ", I7,"   0    specific humidity (kg/kg)         ")') km 
   WRITE(102,'(" cw  ", I7,"   0    total condensate  (kg/kg)         ")') km 
   WRITE(102,'(" rh  ", I7,"   0    relative humidity (%)             ")') km 
   WRITE(102,'(" ENDVARS                                              ")') 
   WRITE(102,'("                                                      ")') 
   
   CLOSE(100) ; CLOSE(101) ; CLOSE(102) 
   
   ENDIF !* ( sigma_level ) -------------------------------------------- 
   
!*********************************************************! 
!  2. Interpolate variables from sigma to pressure levels ! 
!_________________________________________________________! 

   IF ( press_vars ) THEN !* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   
!* Allocate work arrays for pressure levels 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   Allocate( uu_(im,jm,kmp), vv_(im,jm,kmp), ww_(im,jm,kmp) ) 
   Allocate( tt_(im,jm,kmp), td_(im,jm,kmp), qt_(im,jm,kmp) ) 
   Allocate( gh_(im,jm,kmp), qq_(im,jm,kmp), cw_(im,jm,kmp)) 
   Allocate( rh_(im,jm,kmp) )
   Allocate(f_ice_(im,jm,kmp),f_rain_(im,jm,kmp),f_rimef_(im,jm,kmp))
   Allocate(                                     mc_heat_(im,jm,kmp))
   
!* Vertical interpolation to pressure levels (mb) 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   uu_(:,:,:) = undef ; vv_(:,:,:) = undef ; ww_(:,:,:) = undef 
   tt_(:,:,:) = undef ; td_(:,:,:) = undef ; qt_(:,:,:) = undef 
   gh_(:,:,:) = undef ; qq_(:,:,:) = undef ; cw_(:,:,:) = undef 
   rh_(:,:,:) = undef 
   f_ice_(:,:,:) = undef; f_rain_(:,:,:) = undef;f_rimef_(:,:,:) = undef
                                                 mc_heat_(:,:,:) = undef
   
   DO j = 1, jm 
   DO i = 1, im ; k0 = 1 
      DO m = 1, kmp ; xxx = plevel(m) 
      if ( pp(i,j,1) >= xxx .AND. xxx > pp(i,j,km) ) then 
         DO k = k0, km-1 ; IF ( xxx > pp(i,j,k+1) ) EXIT ; ENDDO 
         aaa = (xxx-pp(i,j,k+1))/(pp(i,j,k)-pp(i,j,k+1)) ; bbb = 1 - aaa 
         uu_(i,j,m) = aaa*uu(i,j,k) + bbb*uu(i,j,k+1) 
         vv_(i,j,m) = aaa*vv(i,j,k) + bbb*vv(i,j,k+1) 
         ww_(i,j,m) = aaa*ww(i,j,k) + bbb*ww(i,j,k+1) 
         tt_(i,j,m) = aaa*tt(i,j,k) + bbb*tt(i,j,k+1) 
         td_(i,j,m) = aaa*td(i,j,k) + bbb*td(i,j,k+1) 
         qt_(i,j,m) = aaa*qt(i,j,k) + bbb*qt(i,j,k+1) 
         gh_(i,j,m) = aaa*gh(i,j,k) + bbb*gh(i,j,k+1) 
         qq_(i,j,m) = aaa*qq(i,j,k) + bbb*qq(i,j,k+1) 
         cw_(i,j,m) = aaa*cw(i,j,k) + bbb*cw(i,j,k+1) 
         rh_(i,j,m) = aaa*rh(i,j,k) + bbb*rh(i,j,k+1) 
         f_ice_(i,j,m) = aaa*f_ice(i,j,k) + bbb*f_ice(i,j,k+1)
         f_rain_(i,j,m) = aaa*f_rain(i,j,k) + bbb*f_rain(i,j,k+1)
         f_rimef_(i,j,m) = aaa*f_rimef(i,j,k) + bbb*f_rimef(i,j,k+1)
         mc_heat_(i,j,m) = aaa*mc_heat(i,j,k) + bbb*mc_heat(i,j,k+1)
         k0 = k 
      endif 
      ENDDO 
   ENDDO 
   ENDDO 

!
!  shear Computation. Hua's doing 
!
!  Deep Layer (850-200)

   ioff=18
   joff=18

   aaa=sum(uu_(iloc-ioff:iloc+ioff,jloc-joff:jloc+joff,33))/((2*ioff+1)*(2*ioff+1))- &
       sum(uu_(iloc-ioff:iloc+ioff,jloc-joff:jloc+joff,7))/((2*ioff+1)*(2*ioff+1))
   bbb=sum(vv_(iloc-ioff:iloc+ioff,jloc-joff:jloc+joff,33))/((2*ioff+1)*(2*ioff+1))- &
       sum(vv_(iloc-ioff:iloc+ioff,jloc-joff:jloc+joff,7))/((2*ioff+1)*(2*ioff+1))
   ccc=sqrt(aaa**2+bbb**2) 

   WRITE(20,1012)ifcst_hour,ifcst_minute,aaa,bbb,ccc,pmin !p0(iloc,jloc)

1012 FORMAT(I5,',',2(F11.4,','),F11.4,',',2F11.4)
!
!===========================================================================================================================
!  HDOBS SECTION: if HDOBS are available open the file, proceed as usual otherwise
!  -------------------------------------------------------------------------------
!  README: HDOBS comparison requires 
!  (1) file 1: out_hdob_YYYYMODYHR copied as hdobs.dat (one directory per file - does not do multiple windows) 
!  (2) file 2: obshead.dat:: maxobs,flight_start,flight_end,# eye penetration, HWRF simulation day starting at 1
!  (3) output files: flight_track*.txt, verify_vars1*.txt and verify_vars2*.txt (verified at forecast hour) 
!============================================================================================================================
!
      open(401,file='obshead.dat',form='formatted',status='old')
      read(401,*,IOSTAT=proceed)maxobs,flight_start,flight_end,flegs,nday  ! (e.g.: 980 075030 153000 3 2)
!
      if(maxobs .le. 1)then
        write(0,*)'WARNING: obshead format: maxobs,flight_start,flight_end,flegs,nday'
        stop
      endif
      if(flegs .le. 2)then
         write(0,*)'WARNING: CHECK EYE PENETRATIONS; MUST BE GREATER THAN 1'
         write(0,*)'WARNING: obshead format: maxobs,flight_start,flight_end,flegs,nday' 
      endif
      if(nday .lt. 1)then
         write(0,*)'WARNING: CHECK obshead.dat; nday MUST BE GREATER THAN 1'
         write(0,*)'WARNING: obshead format: maxobs,flight_start,flight_end,flegs,nday'
      endif
!      
      allocate(hdtime(maxobs),hdlat(maxobs),hdlon(maxobs))
      allocate(hdppp(maxobs),hdgpm(maxobs),hdslp(maxobs),hdrain(maxobs))
      allocate(hdtt(maxobs),hdtd(maxobs),hdwd(maxobs),hdws(maxobs))
!
      if(proceed .lt. 0)goto 113
       if(flight_end .lt. flight_start)then
         local_clock=flight_end+240000      ! clock correction
       else
         local_clock=flight_end
       endif
       iset1=(nday-1)*24 + int(flight_start/10000) ! convert flight time to model forecast hour               
       iset2=(nday-1)*24 + local_clock/10000       ! convert flight time to model forecast hour        
!
       write(0,*)'HDOBS TIME WINDOW IN MODEL FORECAST HOURS',iset1,ifcst_hour,iset2
       if(ifcst_hour .gt. iset1 .and. ifcst_hour .lt. iset2)then ! hour window
        filename1='flight_track_'//DateStr//'.txt'
        filename2='verify_varsA_'//DateStr//'.txt'
        filename3='verify_varsB_'//DateStr//'.txt'
        OPEN(402,file='hdobs.dat',form='formatted',status='old')
        OPEN(403,file=filename1,form='formatted',status='unknown',position='append')
        OPEN(404,file=filename2,form='formatted',status='unknown',position='append')
        OPEN(405,file=filename3,form='formatted',status='unknown',position='append')
        write(0,*)'------------------------------------------'
        write(0,*)'READING HDOBS'
        do i=1,maxobs
           read(402,*)obs1,obs2,obs3,obs4,obs5,obs6,obs7,obs8,obs9,obs10,obs11,obs12,obs13
           hdtime(i)   = obs1
           if(flight_end .lt. flight_start .and. obs1 .lt. 120000)hdtime(i)=obs1+240000
           hdlat(i)  = obs2/100; hdlon(i)  = obs3/100
           hdppp(i)  = obs4; hdgpm(i)  = obs5; hdslp(i)  = obs6
           hdtt(i)   = obs7/10; hdtd(i)   = obs8/10
           hdwd(i)   = int(obs9/1000)
           hdws(i)   = (obs9-int(obs9/1000)*1000)*0.514444  ! wind speed in m/s
           hdrain(i) = obs12
           write(0,1013)hdtime(i),hdlat(i),hdlon(i),hdslp(i),hdtt(i),hdtd(i),hdws(i)
        enddo
        write(0,*)'------------------------------------------'

!         do jj=1,jm; do ii=1,im
!          if(p0(ii,jj) .lt. 970.00)write(408,1017)ii,jj,yrs(ii,jj),xrs(ii,jj),p0(ii,jj),sqrt((uu_(ii,jj,14)**2)+(vv_(ii,jj,14)**2))
!         enddo; enddo

!
!       Storm-Centric Transformation in model space (adjust flight path to
!       forecasted MSLP) 
!
         do kkont=1,flegs-1      ! adjustment is done for each flight leg 
!
           kkk=int(maxobs/(flegs-1))
           i1=1 + kkk*(kkont-1)
           i2=i1+kkk-1
           aaa=MinVal(hdslp(i1:i2)) ! determine the minimum MSLP for that flight leg
           write(0,*)
           write(0,*)'ADJUST FLIGHT TO MODEL SPACE BASED EYE PENETRATION #',kkont
           do i=i1,i2                                            
            if(hdslp(i) .eq. aaa)then
             bbb=hdlon(i)-xrs(iloc1,jloc1)
             ccc=hdlat(i)-yrs(iloc1,jloc1)
             write(0,*)'OBSERVED and FORECASTED MINUMUM MSLP',aaa,pmin
             write(0,*)'OBSERVED   LAT/LON',hdlon(i),hdlat(i)
             write(0,*)'FORECASTED LAT/LON',xrs(iloc1,jloc1),yrs(iloc1,jloc1)
             write(0,*)'CORRECTION  FACTOR',bbb,ccc
             if(bbb .ge. 1.00 .or. ccc .ge. 1.00)write(0,*)'WARNING: TRACK DEVIATION BETWEEN FORECAST AND HDOBS > 1.0 DEGREES'
             if(aaa .ge. 1000.0)then
              write(0,*)'WARNING: SOMETHING IS WRONG WITH OBSERVED MINIMUM MSLP IN FLIGHT LEG',kkont
              write(0,*)'WARNING: CHECK obshead.dat AND REDUCE THE FLIGHT LEGS'
              stop
             endif
            endif
           enddo
!
           do i=i1,i2
            hdlon(i)=hdlon(i)-bbb
            hdlat(i)=hdlat(i)-ccc
           enddo
!
         enddo
!
          write(403,*)'model forecast hour =',ifcst_hour
          write(404,*)'model forecast hour =',ifcst_hour   
          write(405,*)'model forecast hour =',ifcst_hour 
          do i=1,maxobs   ! find where the obs lat lon lies in the model grid box
!          write(0,*)'--------------------','obs number',i
!          write(0,*)'flight lat lon',hdlon(i),hdlat(i)
           do jj=nrad,jm-nrad; do ii=nrad,im-nrad
              if(hdlat(i) .ge. yrs(ii,jj-nrad) .and. hdlat(i) .le. yrs(ii,jj+nrad) .and. &
                 hdlon(i) .ge. xrs(ii-nrad,jj) .and. hdlon(i) .le. xrs(ii+nrad,jj))then
                 do k=1,kmp  ! get the height that matches observations 
                   if(hdgpm(i) .ge. gh_(ii,jj,k) .and. hdgpm(i) .le. gh_(ii,jj,k+1))kkk=k
                 enddo
                 write(403,1014)i,hdtime(i),hdlon(i),hdlat(i),xrs(ii,jj),yrs(ii,jj),hdgpm(i),0.5*(gh_(ii,jj,kkk)+gh_(ii,jj,kkk+1))
                 write(404,1014)i,hdtime(i),xrs(ii,jj),yrs(ii,jj),hdslp(i),p0(ii,jj),hdtt(i),(tt_(ii,jj,kkk)-273.14),hdws(i),  &
                                sqrt((uu_(ii,jj,kkk)**2)+(vv_(ii,jj,kkk)**2))
                 if(hdtd(i) .ne. -99.9)then    ! skip flagged data
                  if(hdrain(i) .eq. -999)hdrain(i)=0
                  write(405,1014)i,hdtime(i),hdtd(i),(td_(ii,jj,kkk)-273.14),hdrain(i),prc(ii,jj)*1000.
                 endif
              endif
           enddo; enddo
          enddo            ! go to next observation
!
       endif   ! for observation hour window
113   write(0,*)'SKIPPING OBS BECAUSE hdobs.dat IS NOT AVAILABLE' ! if hdobs and obshead.dat not available simply proceed with diapost  
!
1013  format(I7,2x,6F12.2)
1014  format(I4,2x,I6,2x,12(F10.4,2x))
1017  format(2(I4,2x),12(F10.4,2x))
      close(401)
      close(402) 
      deallocate(hdtime,hdlat,hdlon)
      deallocate(hdppp,hdgpm,hdslp,hdrain)
      deallocate(hdtt,hdtd,hdwd,hdws)
!
   ENDIF !* ( press_vars ) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!           *     *     *     *     *     *     *     *     *     * 
   
   IF ( press_level ) THEN !* ------------------------------------------ 
   
!* Output variables on pressure levels 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   IF ( dlon <= refine ) THEN !* - - - - - - - - - - - - - - - - - - 

   WRITE(200) xrs ; WRITE(200) yrs 
   WRITE(200) zs  ; WRITE(200) ps 
   WRITE(200) ts  ; WRITE(200) qs 
   WRITE(200) p0  ; WRITE(200) prc
   WRITE(200) u10 ; WRITE(200) v10 
   WRITE(200) t10 ; WRITE(200) q10 
   WRITE(200) sst ; WRITE(200) ssq
   WRITE(200) tw  ; WRITE(200) qw 
   WRITE(200) uf  ; WRITE(200) z0 
   do k =  1, kmp ; WRITE(200) uu_(:,:,k) ; enddo 
   do k =  1, kmp ; WRITE(200) vv_(:,:,k) ; enddo 
   do k =  1, kmp ; WRITE(200) ww_(:,:,k) ; enddo 
   do k =  1, kmp ; WRITE(200) tt_(:,:,k) ; enddo 
   do k =  1, kmp ; WRITE(200) td_(:,:,k) ; enddo 
   do k =  1, kmp ; WRITE(200) qt_(:,:,k) ; enddo 
   do k =  1, kmp ; WRITE(200) gh_(:,:,k) ; enddo 
   do k =  1, kmp ; WRITE(200) qq_(:,:,k) ; enddo 
   do k =  1, kmp ; WRITE(200) cw_(:,:,k) ; enddo 
   do k =  1, kmp ; WRITE(200) rh_(:,:,k) ; enddo 
   
   nvars = 28 ; dxx = dlon ; dyy = dlat 
    imr  = im ; jmr =  jm  ; ims =  im ; jms = (nint(aspect*ims)/2)*2 + 1 
   
!* Refine the grid for low resolutions 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   ELSE !* ( dlon > refine ) - - - - - - - - - - - - - - - - - - - - 
   
   im3 = 3*im - 2 ; jm3 = 3*jm - 2 ; k2d = 18 
   Allocate( uu3(im3,jm3,kmp), vv3(im3,jm3,kmp), ww3(im3,jm3,kmp) ) 
   Allocate( tt3(im3,jm3,kmp), td3(im3,jm3,kmp), qt3(im3,jm3,kmp) ) 
   Allocate( gh3(im3,jm3,kmp), qq3(im3,jm3,kmp), cw3(im3,jm3,kmp)) 
   Allocate( rh3(im3,jm3,kmp), wks(im3,jm3,k2d) ) 
   Allocate( f_ice3(im3,jm3,kmp),f_rain3(im3,jm3,kmp),f_rimef3(im3,jm3,kmp))
   Allocate(                                          mc_heat3(im3,jm3,kmp))
   
   CALL Refine_Grid3x( xrs(1,1), wks(1,1,1),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( yrs(1,1), wks(1,1,2),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  zs(1,1), wks(1,1,3),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  ps(1,1), wks(1,1,4),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  ts(1,1), wks(1,1,5),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  qs(1,1), wks(1,1,6),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  p0(1,1), wks(1,1,7),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( prc(1,1), wks(1,1,8),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( u10(1,1), wks(1,1,9),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( v10(1,1), wks(1,1,10), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( t10(1,1), wks(1,1,11), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( q10(1,1), wks(1,1,12), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( sst(1,1), wks(1,1,13), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( ssq(1,1), wks(1,1,14), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  tw(1,1), wks(1,1,15), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  qw(1,1), wks(1,1,16), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  uf(1,1), wks(1,1,17), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  z0(1,1), wks(1,1,18), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( uu_, uu3, undef, im3, jm3, im, jm, kmp )
   CALL Refine_Grid3x( vv_, vv3, undef, im3, jm3, im, jm, kmp )
   CALL Refine_Grid3x( ww_, ww3, undef, im3, jm3, im, jm, kmp )
   CALL Refine_Grid3x( tt_, tt3, undef, im3, jm3, im, jm, kmp )
   CALL Refine_Grid3x( td_, td3, undef, im3, jm3, im, jm, kmp )
   CALL Refine_Grid3x( qt_, qt3, undef, im3, jm3, im, jm, kmp )
   CALL Refine_Grid3x( gh_, gh3, undef, im3, jm3, im, jm, kmp )
   CALL Refine_Grid3x( qq_, qq3, undef, im3, jm3, im, jm, kmp )
   CALL Refine_Grid3x( cw_, cw3, undef, im3, jm3, im, jm, kmp )
   CALL Refine_Grid3x( rh_, rh3, undef, im3, jm3, im, jm, kmp )
   CALL Refine_Grid3x( f_ice_, f_ice3, undef, im3, jm3, im, jm, kmp )
   CALL Refine_Grid3x( f_rain_, f_rain3, undef, im3, jm3, im, jm, kmp )
   CALL Refine_Grid3x( f_rimef_, f_rimef3, undef, im3, jm3, im, jm, kmp )
   CALL Refine_Grid3x( mc_heat_, mc_heat3, undef, im3, jm3, im, jm, kmp )
   
!* Positive-definiteness constraint 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   wks(:,:,6)  = max( wks(:,:,6),  0. )
   wks(:,:,8)  = max( wks(:,:,8),  0. )
   wks(:,:,12) = max( wks(:,:,12), 0. )
   wks(:,:,14) = max( wks(:,:,14), 0. )
   wks(:,:,17) = max( wks(:,:,17), 0. )
   wks(:,:,18) = max( wks(:,:,18), 0. )
   cw3(:,:,:)  = max( cw3(:,:,:),  0. )
   rh3(:,:,:)  = min( rh3(:,:,:), 100.)
   td3(:,:,:)  = min( td3(:,:,:), tt3(:,:,:) )
   f_ice3(:,:,:)  = max( f_ice3(:,:,:),  0. )
   f_rain3(:,:,:)  = max( f_rain3(:,:,:),  0. )
   f_rimef3(:,:,:)  = max( f_rimef3(:,:,:),  0. )
   mc_heat3(:,:,:)  = max( mc_heat3(:,:,:),  0. )
   
   do k = 1, 18  ; WRITE(200) wks(:,:,k) ; enddo 
   do k = 1, kmp ; WRITE(200) uu3(:,:,k) ; enddo 
   do k = 1, kmp ; WRITE(200) vv3(:,:,k) ; enddo 
   do k = 1, kmp ; WRITE(200) ww3(:,:,k) ; enddo 
   do k = 1, kmp ; WRITE(200) tt3(:,:,k) ; enddo 
   do k = 1, kmp ; WRITE(200) td3(:,:,k) ; enddo 
   do k = 1, kmp ; WRITE(200) qt3(:,:,k) ; enddo 
   do k = 1, kmp ; WRITE(200) gh3(:,:,k) ; enddo 
   do k = 1, kmp ; WRITE(200) qq3(:,:,k) ; enddo 
   do k = 1, kmp ; WRITE(200) cw3(:,:,k) ; enddo 
   do k = 1, kmp ; WRITE(200) rh3(:,:,k) ; enddo 
   
   nvars = 28 ; dxx = dlon/3. ; dyy = dlat/3. 
   imr = im3 ; jmr = jm3 ; ims = im3 ; jms = (nint(aspect*ims)/2)*2 + 1 
   
   DeAllocate( uu3, vv3, ww3, tt3, td3, qt3, gh3, qq3, cw3, rh3, wks ) 
   DeAllocate( f_ice3,f_rain3,f_rimef3)
   DeAllocate(                mc_heat3)
   
   ENDIF !* ( refine ) - - - - - - - - - - - - - - - - - - - - - - - 

!           *     *     *     *     *     *     *     *     *     * 
   
!* GrADS control file for pressure-level data on rotated grid 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   xxr1 = .5*(1-imr)*dxx ; xxr_ = dxx !* Lagrangian coord for moving nest 
   yyr1 = .5*(1-jmr)*dyy ; yyr_ = dyy !* and for fixed outer domain (d01) 
   
   WRITE(201,'("                                                      ")') 
   WRITE(201,'(" DSET    press_lev.dat                                ")') 
   WRITE(201,'(" TITLE   Pressure-level data on rotated latlon grid   ")') 
   WRITE(201,'(" OPTIONS sequential big_endian                        ")') 
   WRITE(201,'(" UNDEF",E11.2)') undef 
   WRITE(201,'("                                                      ")') 
   WRITE(201,'(" XDEF", I7,"   LINEAR", 2F11.5 )') imr, xxr1, xxr_ 
   WRITE(201,'(" YDEF", I7,"   LINEAR", 2F11.5 )') jmr, yyr1, yyr_ 
   WRITE(201,'(" ZDEF", I7,"   LEVELS")') kmp 
   write(201,*)  plevel 
   WRITE(201,'(" TDEF", I7,"   LINEAR   00:00Z01JAN0001  1yr")') ntime 
   WRITE(201,'("                                                      ")') 
   WRITE(201,'(" VARS", I3)') nvars
   WRITE(201,'(" xx        0   0    Earth longitude (deg)             ")') 
   WRITE(201,'(" yy        0   0    Earth latitude  (deg)             ")') 
   WRITE(201,'(" zs        0   0    Earth topography (m)              ")') 
   WRITE(201,'(" ps        0   0    surface air   pressure (mb)       ")') 
   WRITE(201,'(" ts        0   0    surface air temperature (K)       ")') 
   WRITE(201,'(" qs        0   0    surface specific humidity (kg/kg) ")') 
   WRITE(201,'(" p0        0   0    mean-sea-level pressure (mb)      ")') 
   WRITE(201,'(" prc       0   0    precipitation (m)                 ")') 
   WRITE(201,'(" u10       0   0    10m rotated x-wind (m/s)          ")') 
   WRITE(201,'(" v10       0   0    10m rotated y-wind (m/s)          ")') 
   WRITE(201,'(" t10       0   0    10m air temperature (K)           ")') 
   WRITE(201,'(" q10       0   0    10m specific humidity (kg/kg)     ")') 
   WRITE(201,'(" sst       0   0    sea-surface temperature (K)       ")') 
   WRITE(201,'(" ssq       0   0    sea-sfc specific humidity (kg/kg) ")') 
   WRITE(201,'(" tw        0   0    sensible heat flux (W/m^2)        ")') 
   WRITE(201,'(" qw        0   0    latent   heat flux (W/m^2)        ")') 
   WRITE(201,'(" uf        0   0    friction velocity  (m/s)          ")') 
   WRITE(201,'(" z0        0   0    roughness height   (m)            ")') 
   WRITE(201,'(" uu  ", I7,"   0    rotated x-wind (m/s)              ")') kmp 
   WRITE(201,'(" vv  ", I7,"   0    rotated y-wind (m/s)              ")') kmp 
   WRITE(201,'(" ww  ", I7,"   0    vertical  wind (m/s)              ")') kmp 
   WRITE(201,'(" tt  ", I7,"   0    air temperature (K)               ")') kmp 
   WRITE(201,'(" td  ", I7,"   0    dewpoint temperature (K)          ")') kmp 
   WRITE(201,'(" qt  ", I7,"   0    equiv. pot. temperature (K)       ")') kmp 
   WRITE(201,'(" gh  ", I7,"   0    geopotential height (m)           ")') kmp 
   WRITE(201,'(" qq  ", I7,"   0    specific humidity (kg/kg)         ")') kmp 
   WRITE(201,'(" cw  ", I7,"   0    total condensate  (kg/kg)         ")') kmp 
   WRITE(201,'(" rh  ", I7,"   0    relative humidity (%)             ")') kmp 
   WRITE(201,'(" ENDVARS                                              ")') 
   WRITE(201,'("                                                      ")') 
   
!* GrADS control file for pressure-level data on standard grid 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   xxs1 = .5*(1-ims)*dxx + ctr_lon ; xxs_ = dxx !* d01 west  bdry 
   yys1 = .5*(1-jms)*dyy + ctr_lat ; yys_ = dyy !* d01 south bdry 
   
   WRITE(202,'("                                                      ")') 
   WRITE(202,'(" DSET    press_lev.dat                                ")') 
   WRITE(202,'(" TITLE   Pressure-level data on standard latlon grid  ")') 
   WRITE(202,'(" OPTIONS sequential big_endian                        ")') 
   WRITE(202,'(" UNDEF",E11.2)') undef 
   WRITE(202,'("                                                      ")') 
   WRITE(202,'(" PDEF", I7,  I6,   A10,   4F11.5)') & 
                       imr, jmr, 'ETA.U', ctr_lon, ctr_lat, dxx,  dyy 
   WRITE(202,'(" XDEF", I7,"   LINEAR", 2F11.5 )') ims, xxs1, xxs_ 
   WRITE(202,'(" YDEF", I7,"   LINEAR", 2F11.5 )') jms, yys1, yys_ 
   WRITE(202,'(" ZDEF", I7,"   LEVELS")') kmp 
   write(202,*)  plevel 
   WRITE(202,'(" TDEF", I7,"   LINEAR   00:00Z01JAN0001  1yr")') ntime 
   WRITE(202,'("                                                      ")') 
   WRITE(202,'(" VECTORPAIRS   u10,v10  uu,vv                         ")') 
   WRITE(202,'(" VARS", I3)') nvars
   WRITE(202,'(" xx        0   0    Earth longitude (deg)             ")') 
   WRITE(202,'(" yy        0   0    Earth latitude  (deg)             ")') 
   WRITE(202,'(" zs        0   0    Earth topography (m)              ")') 
   WRITE(202,'(" ps        0   0    surface air   pressure (mb)       ")') 
   WRITE(202,'(" ts        0   0    surface air temperature (K)       ")') 
   WRITE(202,'(" qs        0   0    surface specific humidity (kg/kg) ")') 
   WRITE(202,'(" p0        0   0    mean-sea-level pressure (mb)      ")') 
   WRITE(202,'(" prc       0   0    precipitation (m)                 ")') 
   WRITE(202,'(" u10       0   0    10m zonal wind comp (m/s)         ")') 
   WRITE(202,'(" v10       0   0    10m meridional wind (m/s)         ")') 
   WRITE(202,'(" t10       0   0    10m air temperature (K)           ")') 
   WRITE(202,'(" q10       0   0    10m specific humidity (kg/kg)     ")') 
   WRITE(202,'(" sst       0   0    sea-surface temperature (K)       ")') 
   WRITE(202,'(" ssq       0   0    sea-sfc specific humidity (kg/kg) ")') 
   WRITE(202,'(" tw        0   0    sensible heat flux (W/m^2)        ")') 
   WRITE(202,'(" qw        0   0    latent   heat flux (W/m^2)        ")') 
   WRITE(202,'(" uf        0   0    friction velocity  (m/s)          ")') 
   WRITE(202,'(" z0        0   0    roughness height   (m)            ")') 
   WRITE(202,'(" uu  ", I7,"   0    zonal wind comp (m/s)             ")') kmp 
   WRITE(202,'(" vv  ", I7,"   0    meridional wind (m/s)             ")') kmp 
   WRITE(202,'(" ww  ", I7,"   0    vertical  wind (m/s)              ")') kmp 
   WRITE(202,'(" tt  ", I7,"   0    air temperature (K)               ")') kmp 
   WRITE(202,'(" td  ", I7,"   0    dewpoint temperature (K)          ")') kmp 
   WRITE(202,'(" qt  ", I7,"   0    equiv. pot. temperature (K)       ")') kmp 
   WRITE(202,'(" gh  ", I7,"   0    geopotential height (m)           ")') kmp 
   WRITE(202,'(" qq  ", I7,"   0    specific humidity (kg/kg)         ")') kmp 
   WRITE(202,'(" cw  ", I7,"   0    total condensate  (kg/kg)         ")') kmp 
   WRITE(202,'(" rh  ", I7,"   0    relative humidity (%)             ")') kmp 
   WRITE(202,'(" ENDVARS                                              ")') 
   WRITE(202,'("                                                      ")') 
   
   CLOSE(200) ; CLOSE(201) ; CLOSE(202) 
   
   ENDIF !* ( press_level ) -------------------------------------------- 

   IF ( press_vars ) THEN 
   DeAllocate( uu_, vv_, ww_, tt_, td_, qt_, gh_, qq_, cw_, rh_) 
   DeAllocate( f_ice_,f_rain_,f_rimef_)
   DeAllocate(                mc_heat_)
   ENDIF 
   
!********************************************! 
!  2x. Output HFIP basic pressure level data ! 
!____________________________________________! 

   IF ( basic_level ) THEN !* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   
   xxx = d02cut ; IF ( domain_num == 1 ) xxx = d01cut 
   imx = (nint(xxx/dlon)/2)*2+1 ; jmx = (nint(aspect*imx)/2)*2+1 
   
!* Allocate work arrays for HFIP test 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   Allocate( uu_(im,jm,4), vv_(im,jm,4) ) 
   Allocate( tt_(im,jm,4), td_(im,jm,4) ) 
   Allocate( gh_(im,jm,4), vx_(im,jm,4) ) 
   Allocate( uus(im,jm,4), vvs(im,jm,4) ) 
   Allocate( p02(im,jm), td2(im,jm), wk(im,jm) ) 
   
!* Vertical interpolation to HFIP levels (mb) 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   uu_(:,:,:) = undef ; vv_(:,:,:) = undef 
   tt_(:,:,:) = undef ; td_(:,:,:) = undef 
   gh_(:,:,:) = undef ; vx_(:,:,:) = undef 
   
   DO j = 1, jm 
   DO i = 1, im ; k0 = 1 
      DO m = 1, 4 ; xxx = pbasic(m) 
      if ( pp(i,j,1) >= xxx .AND. xxx > pp(i,j,km) ) then 
         DO k = k0, km-1 ; IF ( xxx > pp(i,j,k+1) ) EXIT ; ENDDO 
         aaa = (xxx-pp(i,j,k+1))/(pp(i,j,k)-pp(i,j,k+1)) ; bbb = 1 - aaa 
         uu_(i,j,m) = aaa*uu(i,j,k) + bbb*uu(i,j,k+1) 
         vv_(i,j,m) = aaa*vv(i,j,k) + bbb*vv(i,j,k+1) 
         tt_(i,j,m) = aaa*tt(i,j,k) + bbb*tt(i,j,k+1) 
         td_(i,j,m) = aaa*td(i,j,k) + bbb*td(i,j,k+1) 
         gh_(i,j,m) = aaa*gh(i,j,k) + bbb*gh(i,j,k+1) 
         k0 = k 
      endif 
      ENDDO 
   ENDDO 
   ENDDO 

!* Absolute voritcity (/s)
!* ~~~~~~~~~~~~~~~~~~~~~~~ 
   DO j = 2, jm-1 
   DO i = 2, im-1 
      ccc = 2*o_earth*Sin(yrs(i,j)*deg2rad) !* Coriolis parameter 
      DO k = 1, 4 
      aaa = max( vv_(i+1,j,k), vv_(i-1,j,k), uu_(i,j-1,k), uu_(i,j+1,k) ) 
      IF ( aaa < undef ) THEN 
         vx_(i,j,k) = .5*( vv_(i+1,j,k) - vv_(i-1,j,k) )/xrlen(j) & 
                    + .5*( uu_(i,j-1,k) - uu_(i,j+1,k) )/yrlen + ccc 
      ENDIF 
      ENDDO 
   ENDDO 
   ENDDO 

!* 2m Dewpoint temperature (td2) 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   if ( fcst_time /= 0.) then 
   DO j = 1, jm ; DO i = 1, im 
      zzz = gh(i,j,1)-zs(i,j) ; aaa = (zzz-2)/zzz 
      p02(i,j) =  aaa*ps(i,j) + (1-aaa)*pp(i,j,1)
   ENDDO ; ENDDO  !* 2m air pressure *! 
   CALL Dewpoint_RH( p02(1,1), q02(1,1), t02(1,1), td2(1,1), wk(1,1), im,jm,1 ) 
   endif 
   
!* 10m winds in Earth coordinates (urs, vrs) 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   urs(:,:) = wcos(:,:)*u10(:,:) - wsin(:,:)*v10(:,:) 
   vrs(:,:) = wsin(:,:)*u10(:,:) + wcos(:,:)*v10(:,:) 
   
!* 3D winds in Earth coordinates (uus, vvs) 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   DO k = 1, 4 ; DO j = 1, jm ; DO i = 1, im  
      IF ( max( uu_(i,j,k), vv_(i,j,k) ) == undef ) THEN 
         uus(i,j,k) = undef ; vvs(i,j,k) = undef 
      ELSE 
         uus(i,j,k) = wcos(i,j)*uu_(i,j,k) - wsin(i,j)*vv_(i,j,k) 
         vvs(i,j,k) = wsin(i,j)*uu_(i,j,k) + wcos(i,j)*vv_(i,j,k) 
      ENDIF 
   ENDDO ; ENDDO ; ENDDO 
   
   IF ( dlon <= refine ) THEN !* no refine - - - - - - - - - - - - - 

   dxx = dlon ; dyy = dlat 
   imr = im ; jmr = jm ; ims = imx ; jms = jmx 
   Allocate( uu3(imr,jmr,4), vv3(imr,jmr,4) ) 
   Allocate( tt3(imr,jmr,4), td3(imr,jmr,4) ) 
   Allocate( gh3(imr,jmr,4), vx3(imr,jmr,4) ) 
   Allocate( pp3(imr,jmr,6), qq3(imr,jmr,2) ) 
   Allocate( xxr(imr,jmr),   yyr(imr,jmr)   ) 
   Allocate( xxs(imr,jmr),   yys(imr,jmr)   ) 
   
!* Copy rotated grid's coordinates (degree) 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   xxr(:,:) = xrr(:,:) ; yyr(:,:) = yrr(:,:) 
   xxs(:,:) = xrs(:,:) ; yys(:,:) = yrs(:,:) 
   
!* Define HFIP variables without refinement 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   pp3(:,:,1) =  p0(:,:)   ; pp3(:,:,2) = prc(:,:) 
   pp3(:,:,3) = t02(:,:)   ; pp3(:,:,4) = td2(:,:) 
   pp3(:,:,5) = urs(:,:)   ; pp3(:,:,6) = vrs(:,:) 
   uu3(:,:,:) = uus(:,:,:) ; vv3(:,:,:) = vvs(:,:,:) 
   tt3(:,:,:) = tt_(:,:,:) ; td3(:,:,:) = td_(:,:,:) 
   gh3(:,:,:) = gh_(:,:,:) ; vx3(:,:,:) = vx_(:,:,:) 
   qq3(:,:,1) =  zs(:,:)   ; qq3(:,:,2) =  ps(:,:) 
   
   ELSE !* Refine low resolution - - - - - - - - - - - - - - - - - - 
   
   dxx = dlon/3.   ; dyy = dlat/3. 
   imr = 3*im - 2  ; jmr = 3*jm - 2 
   ims = 3*imx - 2 ; jms = 3*jmx - 2 
   Allocate( uu3(imr,jmr,4), vv3(imr,jmr,4) ) 
   Allocate( tt3(imr,jmr,4), td3(imr,jmr,4) ) 
   Allocate( gh3(imr,jmr,4), vx3(imr,jmr,4) ) 
   Allocate( pp3(imr,jmr,6), qq3(imr,jmr,2) ) 
   Allocate( xxr(imr,jmr),   yyr(imr,jmr)   ) 
   Allocate( xxs(imr,jmr),   yys(imr,jmr)   ) 
   
!* Refined rotated grid's coordinates (degree) 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   DO i = 1, imr ; xxr(i,:) = xrr(1,1) + (i-1)*dxx ; ENDDO 
   DO j = 1, jmr ; yyr(:,j) = yrr(1,1) + (j-1)*dyy ; ENDDO 
   Call Convert_LatLon( 2, imr, jmr, xxs, yys, xxr, yyr, ctr_lon, ctr_lat ) 

!* Define HFIP variables with triple refinement 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   CALL Refine_Grid3x(  p0(1,1), pp3(1,1,1), undef, imr, jmr, im, jm, 1 )
   CALL Refine_Grid3x( prc(1,1), pp3(1,1,2), undef, imr, jmr, im, jm, 1 )
   CALL Refine_Grid3x( t02(1,1), pp3(1,1,3), undef, imr, jmr, im, jm, 1 )
   CALL Refine_Grid3x( td2(1,1), pp3(1,1,4), undef, imr, jmr, im, jm, 1 )
   CALL Refine_Grid3x( urs(1,1), pp3(1,1,5), undef, imr, jmr, im, jm, 1 )
   CALL Refine_Grid3x( vrs(1,1), pp3(1,1,6), undef, imr, jmr, im, jm, 1 )
   
   CALL Refine_Grid3x( uus, uu3, undef, imr, jmr, im, jm, 4 )
   CALL Refine_Grid3x( vvs, vv3, undef, imr, jmr, im, jm, 4 )
   CALL Refine_Grid3x( tt_, tt3, undef, imr, jmr, im, jm, 4 )
   CALL Refine_Grid3x( td_, td3, undef, imr, jmr, im, jm, 4 )
   CALL Refine_Grid3x( gh_, gh3, undef, imr, jmr, im, jm, 4 )
   CALL Refine_Grid3x( vx_, vx3, undef, imr, jmr, im, jm, 4 )
   CALL Refine_Grid3x(  zs(1,1), qq3(1,1,1), undef, imr, jmr, im, jm, 1 )
   CALL Refine_Grid3x(  ps(1,1), qq3(1,1,2), undef, imr, jmr, im, jm, 1 )
   
!* Positive-definiteness constraint 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   pp3(:,:,2) = max( pp3(:,:,2), 0. )
   pp3(:,:,4) = min( pp3(:,:,4), pp3(:,:,3) )
   td3(:,:,:) = min( td3(:,:,:), tt3(:,:,:) )

   ENDIF !* ( refine ) - - - - - - - - - - - - - - - - - - - - - - - 

   IF ( debug_data ) THEN !* --------------------------------------- 

!* Output HFIP test data on rotated grid 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   do n = 1, 6 ; WRITE(252) pp3(:,:,n) ; enddo 
   do k = 1, 4 ; WRITE(252) uu3(:,:,k) ; enddo 
   do k = 1, 4 ; WRITE(252) vv3(:,:,k) ; enddo 
   do k = 1, 4 ; WRITE(252) tt3(:,:,k) ; enddo 
   do k = 1, 4 ; WRITE(252) td3(:,:,k) ; enddo 
   do k = 1, 4 ; WRITE(252) gh3(:,:,k) ; enddo 
   do k = 1, 4 ; WRITE(252) vx3(:,:,k) ; enddo 
   do n = 1, 2 ; WRITE(252) qq3(:,:,n) ; enddo 
   WRITE(252) xxs ; WRITE(252) yys 
   nvars = 16 
   
!* GrADS control file for HFIP test data on rotated grid 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   xxr1 = .5*(1-imr)*dxx ; xxr_ = dxx !* Lagrangian coord for moving nest 
   yyr1 = .5*(1-jmr)*dyy ; yyr_ = dyy !* and  for  the fixed domain (d01) 
   
   WRITE(251,'("                                                      ")') 
   WRITE(251,'(" DSET    basic_rot.dat                                ")') 
   WRITE(251,'(" TITLE   HFIP p-level data on rotated latlon grid     ")') 
   WRITE(251,'(" OPTIONS sequential big_endian                        ")') 
   WRITE(251,'(" UNDEF",E11.2)') undef 
   WRITE(251,'("                                                      ")') 
   WRITE(251,'(" XDEF", I7,"   LINEAR", 2F11.5 )') imr, xxr1, xxr_ 
   WRITE(251,'(" YDEF", I7,"   LINEAR", 2F11.5 )') jmr, yyr1, yyr_ 
   WRITE(251,'(" ZDEF", I7,"   LEVELS")') kmb 
   WRITE(251,*)  pbasic 
   WRITE(251,'(" TDEF", I7,"   LINEAR   00:00Z01JAN0001  1yr")') ntime 
   WRITE(251,'("                                                      ")') 
   WRITE(251,'(" VARS", I3)') nvars
   WRITE(251,'(" p0        0   0    mean-sea-level pressure (mb)      ")') 
   WRITE(251,'(" prc       0   0    precipitation (m)                 ")') 
   WRITE(251,'(" t02       0   0    2-m air temperature (K)           ")') 
   WRITE(251,'(" td2       0   0    2-m dewpoint temperature (K)      ")') 
   WRITE(251,'(" u10       0   0    10m zonal wind comp (m/s)         ")') 
   WRITE(251,'(" v10       0   0    10m meridional wind (m/s)         ")') 
   WRITE(251,'(" uu  ", I7,"   0    zonal wind comp (m/s)             ")') kmb 
   WRITE(251,'(" vv  ", I7,"   0    meridional wind (m/s)             ")') kmb 
   WRITE(251,'(" tt  ", I7,"   0    air temperature  (K)              ")') kmb 
   WRITE(251,'(" td  ", I7,"   0    dewpoint temperature (K)          ")') kmb 
   WRITE(251,'(" gh  ", I7,"   0    geopotential height  (m)          ")') kmb 
   WRITE(251,'(" vx  ", I7,"   0    absolute vorticity  (/s)          ")') kmb 
   WRITE(251,'(" zs        0   0    Earth topography (m)              ")') 
   WRITE(251,'(" ps        0   0    surface pressure (mb)             ")') 
   WRITE(251,'(" xx        0   0    Earth longitude (deg)             ")') 
   WRITE(251,'(" yy        0   0    Earth latitude  (deg)             ")') 
   WRITE(251,'(" ENDVARS                                              ")') 
   WRITE(251,'("                                                      ")') 
   WRITE(251,'("*NOTE: Winds are in Earth coordinates -- NOT rotated! ")') 
   WRITE(251,'("                                                      ")') 
   
   ENDIF !* ( debug_data ) ----------------------------------------- 
   
!           *     *     *     *     *     *     *     *     *     * 
   
!* Prescribed standard grid's coordinates (degree) 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   Allocate( xsr(ims,jms), ysr(ims,jms), xss(ims,jms), yss(ims,jms) ) 
   
   i0 = (im+1)/2 ; j0 = (jm+1)/2 
   xss(1,1) = ctr_lon + ( nint((xrs(i0,j0)-ctr_lon)/dxx) - (ims-1)/2 )*dxx
   yss(1,1) = ctr_lat + ( nint((yrs(i0,j0)-ctr_lat)/dyy) - (jms-1)/2 )*dyy
   DO i = 1, ims ; xss(i,:) = xss(1,1) + (i-1)*dxx ; ENDDO 
   DO j = 1, jms ; yss(:,j) = yss(1,1) + (j-1)*dyy ; ENDDO 
   
   Call Convert_LatLon( 1, ims, jms, xss, yss, xsr, ysr, ctr_lon, ctr_lat ) 
   
!* Interpolate refined data to standard grid 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   Allocate( uu3_(ims,jms,4), vv3_(ims,jms,4) ) 
   Allocate( tt3_(ims,jms,4), td3_(ims,jms,4) ) 
   Allocate( gh3_(ims,jms,4), vx3_(ims,jms,4) ) 
   Allocate( pp3_(ims,jms,6), qq3_(ims,jms,2) ) 
   
   uu3_(:,:,:) = undef ; vv3_(:,:,:) = undef 
   tt3_(:,:,:) = undef ; td3_(:,:,:) = undef 
   gh3_(:,:,:) = undef ; vx3_(:,:,:) = undef 
   pp3_(:,:,:) = undef ; qq3_(:,:,:) = undef 
   
   ibd = 4 ; IF (dlon > refine) ibd = 12 !* bad boundary points 
   xrr1 = xxr(1,1) ; xrr2 = xxr(1+ibd,1) ; xrr8 = xxr(imr-ibd,1) 
   yrr1 = yyr(1,1) ; yrr2 = yyr(1,1+ibd) ; yrr8 = yyr(1,jmr-ibd) 
   
   DO j = 1, jms ; DO i = 1, ims !* -------------------------------> 
   
      xxx = xsr(i,j) ;  yyy = ysr(i,j) 
   IF ( yrr2 <= yyy .AND. yyy < yrr8 ) THEN 
   if ( xrr2 <= xxx .AND. xxx < xrr8 ) then !* - - - - - - - - - - - 
   
      m2  =  1 + int((xxx-xrr1)/dxx) ;    m1  = m2-1 ; m4 = m2+2 
      n2  =  1 + int((yyy-yrr1)/dyy) ;    n1  = n2-1 ; n4 = n2+2 
     aaa  = ( xxx - xxr(m2,n2) )/dxx ;   bbb  = ( yyy - yyr(m2,n2) )/dyy 
   wtx(1) = -(aaa)*(aaa-1)*(aaa-2)/6 ; wtx(2) = (aaa+1)*(aaa-1)*(aaa-2)/2 
   wtx(3) = -(aaa+1)*(aaa)*(aaa-2)/2 ; wtx(4) = (aaa+1)*(aaa)*(aaa-1)/6 
   wty(1) = -(bbb)*(bbb-1)*(bbb-2)/6 ; wty(2) = (bbb+1)*(bbb-1)*(bbb-2)/2 
   wty(3) = -(bbb+1)*(bbb)*(bbb-2)/2 ; wty(4) = (bbb+1)*(bbb)*(bbb-1)/6 
   
   DO k = 1, 6 
   ff = pp3(m1:m4,n1:n4,k); CALL Cubic_Interpo( pp3_(i,j,k),ff,wtx,wty,undef ) 
   ENDDO 
   
   DO k = 1, 4 
   ff = uu3(m1:m4,n1:n4,k); CALL Cubic_Interpo( uu3_(i,j,k),ff,wtx,wty,undef ) 
   ff = vv3(m1:m4,n1:n4,k); CALL Cubic_Interpo( vv3_(i,j,k),ff,wtx,wty,undef ) 
   ff = tt3(m1:m4,n1:n4,k); CALL Cubic_Interpo( tt3_(i,j,k),ff,wtx,wty,undef ) 
   ff = td3(m1:m4,n1:n4,k); CALL Cubic_Interpo( td3_(i,j,k),ff,wtx,wty,undef ) 
   ff = gh3(m1:m4,n1:n4,k); CALL Cubic_Interpo( gh3_(i,j,k),ff,wtx,wty,undef ) 
   ff = vx3(m1:m4,n1:n4,k); CALL Cubic_Interpo( vx3_(i,j,k),ff,wtx,wty,undef ) 
   ENDDO 
   
   DO k = 1, 2 
   ff = qq3(m1:m4,n1:n4,k); CALL Cubic_Interpo( qq3_(i,j,k),ff,wtx,wty,undef ) 
   ENDDO 

   endif  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   ENDIF 
   
   ENDDO ; ENDDO !* <----------------------------------------------- 

!* Positive-definiteness constraint 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   pp3_(:,:,2) = max( pp3_(:,:,2), 0. )
   pp3_(:,:,4) = min( pp3_(:,:,4), pp3_(:,:,3) )
   td3_(:,:,:) = min( td3_(:,:,:), tt3_(:,:,:) )
   
!* Output HFIP test data on standard grid 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   do n = 1, 6 ; WRITE(254) pp3_(:,:,n) ; enddo 
   do k = 1, 4 ; WRITE(254) uu3_(:,:,k) ; enddo 
   do k = 1, 4 ; WRITE(254) vv3_(:,:,k) ; enddo 
   do k = 1, 4 ; WRITE(254) tt3_(:,:,k) ; enddo 
   do k = 1, 4 ; WRITE(254) td3_(:,:,k) ; enddo 
   do k = 1, 4 ; WRITE(254) gh3_(:,:,k) ; enddo 
   do k = 1, 4 ; WRITE(254) vx3_(:,:,k) ; enddo 
   do n = 1, 2 ; WRITE(254) qq3_(:,:,n) ; enddo 
   WRITE(254) xss ; WRITE(254) yss 
   nvars = 16 
   
!* GrADS control file for HFIP data on standard grid 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   xxs_ = dxx ; yys_ = dyy 
   IF ( domain_num == 1 ) THEN 
      xxs1 = xss(1,1) ; yys1 = yss(1,1) 
   ELSE !* Lagrangian coordinates for moving nest 
      xxs1 = .5*(1-ims)*dxx ; yys1 = .5*(1-jms)*dyy 
   ENDIF 
   
   WRITE(253,'("                                                      ")') 
   WRITE(253,'(" DSET    basic_std.dat                                ")') 
   WRITE(253,'(" TITLE   HFIP p-level data on standard latlon grid    ")') 
   WRITE(253,'(" OPTIONS sequential big_endian                        ")') 
   WRITE(253,'(" UNDEF",E11.2)') undef 
   WRITE(253,'("                                                      ")') 
   WRITE(253,'(" XDEF", I7,"   LINEAR", 2F11.5 )') ims, xxs1, xxs_ 
   WRITE(253,'(" YDEF", I7,"   LINEAR", 2F11.5 )') jms, yys1, yys_ 
   WRITE(253,'(" ZDEF", I7,"   LEVELS")') kmb 
   WRITE(253,*)  pbasic 
   WRITE(253,'(" TDEF", I7,"   LINEAR   00:00Z01JAN0001  1yr")') ntime 
   WRITE(253,'("                                                      ")') 
   WRITE(253,'(" VARS", I3)') nvars
   WRITE(253,'(" p0        0   0    mean-sea-level pressure (mb)      ")') 
   WRITE(253,'(" prc       0   0    precipitation (m)                 ")') 
   WRITE(253,'(" t02       0   0    2-m air temperature (K)           ")') 
   WRITE(253,'(" td2       0   0    2-m dewpoint temperature (K)      ")') 
   WRITE(253,'(" u10       0   0    10m zonal wind comp (m/s)         ")') 
   WRITE(253,'(" v10       0   0    10m meridional wind (m/s)         ")') 
   WRITE(253,'(" uu  ", I7,"   0    zonal wind comp (m/s)             ")') kmb 
   WRITE(253,'(" vv  ", I7,"   0    meridional wind (m/s)             ")') kmb 
   WRITE(253,'(" tt  ", I7,"   0    air temperature  (K)              ")') kmb 
   WRITE(253,'(" td  ", I7,"   0    dewpoint temperature (K)          ")') kmb 
   WRITE(253,'(" gh  ", I7,"   0    geopotential height  (m)          ")') kmb 
   WRITE(253,'(" vx  ", I7,"   0    absolute vorticity  (/s)          ")') kmb 
   WRITE(253,'(" zs        0   0    Earth topography (m)              ")') 
   WRITE(253,'(" ps        0   0    surface pressure (mb)             ")') 
   WRITE(253,'(" xx        0   0    Earth longitude (deg)             ")') 
   WRITE(253,'(" yy        0   0    Earth latitude  (deg)             ")') 
   WRITE(253,'(" ENDVARS                                              ")') 
   WRITE(253,'("                                                      ")') 
   
!           *     *     *     *     *     *     *     *     *     * 
   
!* HFIP data dimensions on standard grid 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   WRITE(256,'(I10,I5,I3,I3,1X,2I5,2F9.4,2F10.4,4F8.4)')  & 
   start_date, ifcst_hour, ifcst_minute, domain_num, ims, jms,  & 
   yss(1,jms), yss(1,1), xss(1,1), xss(ims,1), dxx, dyy, dlon, dlat 

!* Message file for HFIP data dimensions 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   write(255,*) '                                                         '
   write(255,*) 'Column  Description of the Data                   Example'
   write(255,*) '------  -------------------------------------  ----------'
   write(255,*) '   1    Starting date:                         2008090800'
   write(255,*) '   2    Forecast hour:                                120'
   write(255,*) '   3    Forecast minute:                               30'
   write(255,*) '   4    Domain number:                                  1'
   write(255,*) '   5    Number of grid points in X-direction:         997'
   write(255,*) '   6    Number of grid points in Y-direction:         721'
   write(255,*) '   7    Northern  most  latitude   (degree):      40.2000'
   write(255,*) '   8    Southern  most  latitude   (degree):      -3.0000'
   write(255,*) '   9    Western   most  longitude  (degree):    -110.0000'
   write(255,*) '  10    Eastern   most  longitude  (degree):     -50.0000'
   write(255,*) '  11    Output resolution in X-dir (degree):       0.0600'
   write(255,*) '  12    Output resolution in Y-dir (degree):       0.0600'
   write(255,*) '  13    Model  resolution in X-dir (degree):       0.1800'
   write(255,*) '  14    Model  resolution in Y-dir (degree):       0.1800'
   write(255,*) '                                                         '
   
   CLOSE(251) ; CLOSE(252) ; CLOSE(253) ; CLOSE(254) 
   CLOSE(255) ; CLOSE(256) 
   
   DeAllocate( uu_,  vv_,  tt_,  td_,  gh_,  vx_,  uus,  vvs, p02, td2, wk ) 
   DeAllocate( uu3,  vv3,  tt3,  td3,  gh3,  vx3,  pp3,  qq3  ) 
   DeAllocate( uu3_, vv3_, tt3_, td3_, gh3_, vx3_, pp3_, qq3_ ) 
   DeAllocate( xxr, yyr, xxs, yys, xsr, ysr, xss, yss ) 
   
   ENDIF !* ( basic_level ) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   
!*******************************************************! 
!  3. Interpolate variables from sigma to height levels ! 
!_______________________________________________________! 

   IF ( height_vars ) THEN !* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!* Allocate work arrays for height levels 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   Allocate( uu_(im,jm,kmz), vv_(im,jm,kmz), ww_(im,jm,kmz) ) 
   Allocate( tt_(im,jm,kmz), td_(im,jm,kmz), qt_(im,jm,kmz) ) 
   Allocate( pp_(im,jm,kmz), qq_(im,jm,kmz), cw_(im,jm,kmz) ) 
   Allocate( ux_(im,jm,kmz), vx_(im,jm,kmz), rh_(im,jm,kmz) )
   Allocate( f_ice_(im,jm,kmz),f_rain_(im,jm,kmz),f_rimef_(im,jm,kmz))
   Allocate(                                      mc_heat_(im,jm,kmz)) 

!* Diffusion terms (gopal's doing)
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   Allocate(  difu_(im,jm,kmz),  difv_(im,jm,kmz), km_diffus_(im,jm,kmz), kh_diffus_(im,jm,kmz) )
   Allocate( hdifu_(im,jm,kmz), hdifv_(im,jm,kmz))
 
   uu_(:,:,:) = undef ; vv_(:,:,:) = undef ; ww_(:,:,:) = undef 
   tt_(:,:,:) = undef ; td_(:,:,:) = undef ; qt_(:,:,:) = undef 
   pp_(:,:,:) = undef ; qq_(:,:,:) = undef ; cw_(:,:,:) = undef 
   ux_(:,:,:) = undef ; vx_(:,:,:) = undef ; rh_(:,:,:) = undef 
   f_ice_(:,:,:) = undef; f_rain_(:,:,:) = undef; f_rimef_(:,:,:) = undef
                                                  mc_heat_(:,:,:) = undef

   difu_(:,:,:)=undef; difv_(:,:,:)=undef   ! gopal's doing for diffusion terms
   km_diffus_(:,:,:)=undef; kh_diffus_(:,:,:)=undef
   hdifu_(:,:,:)=undef; hdifv_(:,:,:)=undef
!

!* Take care of height zlevel(1) that could be less than lowest geopotential height 


!* Vertical interpolation to height levels (m) 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   DO j = 1, jm 
   DO i = 1, im ; k0 = 1 
      DO m = 1, kmz ; xxx = zlevel(m)  
       if ( gh(i,j,1) >= xxx)xxx=gh(i,j,1) ! fix for filling zlevel if those are below gh(1) 
       if ( gh(i,j,1) <= xxx .AND. xxx < gh(i,j,km) ) then 
         DO k = k0, km-1 ; IF ( xxx < gh(i,j,k+1) ) EXIT ; ENDDO
         aaa = (gh(i,j,k+1)-xxx)/(gh(i,j,k+1)-gh(i,j,k)) ; bbb = 1 - aaa 
         uu_(i,j,m) = aaa*uu(i,j,k) + bbb*uu(i,j,k+1) 
         vv_(i,j,m) = aaa*vv(i,j,k) + bbb*vv(i,j,k+1) 
         ww_(i,j,m) = aaa*ww(i,j,k) + bbb*ww(i,j,k+1) 
         tt_(i,j,m) = aaa*tt(i,j,k) + bbb*tt(i,j,k+1) 
         td_(i,j,m) = aaa*td(i,j,k) + bbb*td(i,j,k+1) 
         qt_(i,j,m) = aaa*qt(i,j,k) + bbb*qt(i,j,k+1) 
         pp_(i,j,m) = aaa*pp(i,j,k) + bbb*pp(i,j,k+1) 
         qq_(i,j,m) = aaa*qq(i,j,k) + bbb*qq(i,j,k+1) 
         cw_(i,j,m) = aaa*cw(i,j,k) + bbb*cw(i,j,k+1) 
         rh_(i,j,m) = aaa*rh(i,j,k) + bbb*rh(i,j,k+1) 
         f_ice_(i,j,m) = aaa*f_ice(i,j,k) + bbb*f_ice(i,j,k+1)
         f_rain_(i,j,m) = aaa*f_rain(i,j,k) + bbb*f_rain(i,j,k+1)
         f_rimef_(i,j,m) = aaa*f_rimef(i,j,k) + bbb*f_rimef(i,j,k+1)
         mc_heat_(i,j,m) = aaa*mc_heat(i,j,k) + bbb*mc_heat(i,j,k+1)

! gopal's doing for the diffusion term
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         km_diffus_(i,j,m) = aaa*km_diffus(i,j,k) + bbb*km_diffus(i,j,k+1)
         kh_diffus_(i,j,m) = aaa*kh_diffus(i,j,k) + bbb*kh_diffus(i,j,k+1)

         difu_(i,j,m) = aaa*difu(i,j,k) + bbb*difu(i,j,k+1)
         difv_(i,j,m) = aaa*difv(i,j,k) + bbb*difv(i,j,k+1)

         hdifu_(i,j,m) = aaa*hdifu(i,j,k) + bbb*hdifu(i,j,k+1)
         hdifv_(i,j,m) = aaa*hdifv(i,j,k) + bbb*hdifv(i,j,k+1)

         k0 = k 
         kstart=k0
       else 
         WRITE(0,*)'PROBLEMS WITH VERTICAL INTERPOLATION',I,J,K
         WRITE(0,*)gh(i,j,1),xxx,gh(i,j,km)
       endif
      ENDDO 
   ENDDO 
   ENDDO 

!* Divergence & vorticity  
!* ~~~~~~~~~~~~~~~~~~~~~~~
   DO k = 1, kmz ; DO j = 2, jm-1 ; DO i = 2, im-1 
      aaa = max( uu_(i+1,j,k), uu_(i-1,j,k), vv_(i,j+1,k), vv_(i,j-1,k) )  
      IF ( aaa < undef ) THEN 
         ux_(i,j,k) = .5*( uu_(i+1,j,k) - uu_(i-1,j,k) )/xrlen(j) & 
                    + .5*( vv_(i,j+1,k) - vv_(i,j-1,k) )/yrlen
      ENDIF 
      bbb = max( vv_(i+1,j,k), vv_(i-1,j,k), uu_(i,j-1,k), uu_(i,j+1,k) ) 
      IF ( bbb < undef ) THEN 
         vx_(i,j,k) = .5*( vv_(i+1,j,k) - vv_(i-1,j,k) )/xrlen(j) & 
                    + .5*( uu_(i,j-1,k) - uu_(i,j+1,k) )/yrlen
      ENDIF 
   ENDDO ; ENDDO ; ENDDO 

   ENDIF !* ( height_vars ) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!           *     *     *     *     *     *     *     *     *     * 
   
   IF ( height_level ) THEN !* ----------------------------------------- 

!* Output variables on height levels 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   IF ( dlon <= refine ) THEN !* - - - - - - - - - - - - - - - - - - 

   WRITE(300) xrs ; WRITE(300) yrs 
   WRITE(300) zs  ; WRITE(300) ps 
   WRITE(300) ts  ; WRITE(300) qs 
   WRITE(300) p0  ; WRITE(300) prc
   WRITE(300) u10 ; WRITE(300) v10 
   WRITE(300) t10 ; WRITE(300) q10 
   WRITE(300) sst ; WRITE(300) ssq
   WRITE(300) tw  ; WRITE(300) qw 
   WRITE(300) uf  ; WRITE(300) z0 
   do k =  1, kmz ; WRITE(300) uu_(:,:,k)  ; enddo 
   do k =  1, kmz ; WRITE(300) vv_(:,:,k)  ; enddo 
   do k =  1, kmz ; WRITE(300) ww_(:,:,k)  ; enddo 
   do k =  1, kmz ; WRITE(300) tt_(:,:,k)  ; enddo 
   do k =  1, kmz ; WRITE(300) vx_(:,:,k)  ; enddo 
   do k =  1, kmz ; WRITE(300) qt_(:,:,k)  ; enddo 
   do k =  1, kmz ; WRITE(300) pp_(:,:,k)  ; enddo 
   do k =  1, kmz ; WRITE(300) qq_(:,:,k)  ; enddo 
   do k =  1, kmz ; WRITE(300) cw_(:,:,k)  ; enddo 
   do k =  1, kmz ; WRITE(300) rh_(:,:,k)  ; enddo
!  do k =  1, kmz ; WRITE(300) qt_(:,:,k) - qt_(im,jm,k) ; enddo !thetae pert change
!  do k =  1, kmz ; WRITE(300) tt_(:,:,k) - tt_(im,jm,k) ; enddo 
   do k =  1, kmz ; WRITE(300) f_ice_(:,:,k); enddo
   do k =  1, kmz ; WRITE(300) f_rain_(:,:,k); enddo
   do k =  1, kmz ; WRITE(300) f_rimef_(:,:,k); enddo
   
   nvars = 31 ; dxx = dlon ; dyy = dlat 
    imr  = im ; jmr =  jm  ; ims =  im ; jms = (nint(aspect*ims)/2)*2 + 1 
   
!* Refine the grid for low resolutions 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   ELSE !* ( dlon > refine ) - - - - - - - - - - - - - - - - - - - - 
   
   im3 = 3*im - 2 ; jm3 = 3*jm - 2 ; k2d = 18 
   Allocate( uu3(im3,jm3,kmz), vv3(im3,jm3,kmz), ww3(im3,jm3,kmz) ) 
   Allocate( tt3(im3,jm3,kmz), td3(im3,jm3,kmz), qt3(im3,jm3,kmz) ) 
   Allocate( pp3(im3,jm3,kmz), qq3(im3,jm3,kmz), cw3(im3,jm3,kmz)) 
   Allocate( rh3(im3,jm3,kmz), wks(im3,jm3,k2d) )
   Allocate( f_ice3(im3,jm3,kmz),f_rain3(im3,jm3,kmz),f_rimef3(im3,jm3,kmz))
   Allocate(                                          mc_heat3(im3,jm3,kmz))
!
   Allocate(  difu3(im3,jm3,kmz),  difv3(im3,jm3,kmz), km_diffus3(im3,jm3,kmz),kh_diffus3(im3,jm3,kmz))  ! gopal's doing for diffusion terms
   Allocate( hdifu3(im3,jm3,kmz), hdifv3(im3,jm3,kmz))
   
   CALL Refine_Grid3x( xrs(1,1), wks(1,1,1),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( yrs(1,1), wks(1,1,2),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  zs(1,1), wks(1,1,3),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  ps(1,1), wks(1,1,4),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  ts(1,1), wks(1,1,5),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  qs(1,1), wks(1,1,6),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  p0(1,1), wks(1,1,7),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( prc(1,1), wks(1,1,8),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( u10(1,1), wks(1,1,9),  undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( v10(1,1), wks(1,1,10), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( t10(1,1), wks(1,1,11), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( q10(1,1), wks(1,1,12), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( sst(1,1), wks(1,1,13), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( ssq(1,1), wks(1,1,14), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  tw(1,1), wks(1,1,15), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  qw(1,1), wks(1,1,16), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  uf(1,1), wks(1,1,17), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x(  z0(1,1), wks(1,1,18), undef, im3, jm3, im, jm, 1 )
   CALL Refine_Grid3x( uu_, uu3, undef, im3, jm3, im, jm, kmz )
   CALL Refine_Grid3x( vv_, vv3, undef, im3, jm3, im, jm, kmz )
   CALL Refine_Grid3x( ww_, ww3, undef, im3, jm3, im, jm, kmz )
   CALL Refine_Grid3x( tt_, tt3, undef, im3, jm3, im, jm, kmz )
   CALL Refine_Grid3x( td_, td3, undef, im3, jm3, im, jm, kmz )
   CALL Refine_Grid3x( qt_, qt3, undef, im3, jm3, im, jm, kmz )
   CALL Refine_Grid3x( pp_, pp3, undef, im3, jm3, im, jm, kmz )
   CALL Refine_Grid3x( qq_, qq3, undef, im3, jm3, im, jm, kmz )
   CALL Refine_Grid3x( cw_, cw3, undef, im3, jm3, im, jm, kmz )
   CALL Refine_Grid3x( rh_, rh3, undef, im3, jm3, im, jm, kmz )
   CALL Refine_Grid3x( f_ice_, f_ice3, undef, im3, jm3, im, jm, kmz )
   CALL Refine_Grid3x( f_rain_, f_rain3, undef, im3, jm3, im, jm, kmz )
   CALL Refine_Grid3x( f_rimef_, f_rimef3, undef, im3, jm3, im, jm, kmz )
   CALL Refine_Grid3x( mc_heat_, mc_heat3, undef, im3, jm3, im, jm, kmz )
!
   CALL Refine_Grid3x( difu_, difu3, undef, im3, jm3, im, jm, kmz )  ! gopal's doing for vertical diffusion terms
   CALL Refine_Grid3x( difv_, difv3, undef, im3, jm3, im, jm, kmz )
   CALL Refine_Grid3x( km_diffus_, km_diffus3, undef, im3, jm3, im, jm, kmz )  ! gopal's doing for vertical diffusion terms
   CALL Refine_Grid3x( kh_diffus_, kh_diffus3, undef, im3, jm3, im, jm, kmz )

   CALL Refine_Grid3x( hdifu_, hdifu3, undef, im3, jm3, im, jm, kmz )  ! gopal's doing for horizontal diffusion terms
   CALL Refine_Grid3x( hdifv_, hdifv3, undef, im3, jm3, im, jm, kmz )
   
!* Positive-definiteness constraint 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   wks(:,:,6)  = max( wks(:,:,6),  0. )
   wks(:,:,8)  = max( wks(:,:,8),  0. )
   wks(:,:,12) = max( wks(:,:,12), 0. )
   wks(:,:,14) = max( wks(:,:,14), 0. )
   wks(:,:,17) = max( wks(:,:,17), 0. )
   wks(:,:,18) = max( wks(:,:,18), 0. )
   cw3(:,:,:)  = max( cw3(:,:,:),  0. )
   rh3(:,:,:)  = min( rh3(:,:,:), 100.)
   td3(:,:,:)  = min( td3(:,:,:), tt3(:,:,:) )
   f_ice3(:,:,:)  = max( f_ice3(:,:,:),  0. )
   f_rain3(:,:,:)  = max( f_rain3(:,:,:),  0. )
   f_rimef3(:,:,:)  = max( f_rimef3(:,:,:),  0. )
   mc_heat3(:,:,:)  = max( mc_heat3(:,:,:),  0. )
!
   difu3(:,:,:)  = max( difu3(:,:,:),  0. )      ! gopal's doing for vertical diffusion term
   difv3(:,:,:)  = max( difv3(:,:,:),  0. )
   km_diffus3(:,:,:)  = max( km_diffus3(:,:,:),  0. )
   kh_diffus3(:,:,:)  = max( kh_diffus3(:,:,:),  0. )

   hdifu3(:,:,:)  = max( hdifu3(:,:,:),  0. )      ! gopal's doing for horizontal diffusion term
   hdifv3(:,:,:)  = max( hdifv3(:,:,:),  0. )

   
   do k = 1, 18  ; WRITE(300) wks(:,:,k) ; enddo 
   do k = 1, kmz ; WRITE(300) uu3(:,:,k) ; enddo 
   do k = 1, kmz ; WRITE(300) vv3(:,:,k) ; enddo 
   do k = 1, kmz ; WRITE(300) ww3(:,:,k) ; enddo 
   do k = 1, kmz ; WRITE(300) tt3(:,:,k) ; enddo 
   do k = 1, kmz ; WRITE(300) td3(:,:,k) ; enddo 
   do k = 1, kmz ; WRITE(300) qt3(:,:,k) ; enddo 
   do k = 1, kmz ; WRITE(300) pp3(:,:,k) ; enddo 
   do k = 1, kmz ; WRITE(300) qq3(:,:,k) ; enddo 
   do k = 1, kmz ; WRITE(300) cw3(:,:,k) ; enddo 
   do k = 1, kmz ; WRITE(300) rh3(:,:,k) ; enddo 
!   do k =  1, kmz ; WRITE(300) qt3(:,:,k) - qt3(im,jm,k) ; enddo   !thetae pert change
!   do k =  1, kmz ; WRITE(300) tt3(:,:,k) - tt3(im,jm,k) ; enddo
   do k =  1, kmz ; WRITE(300) f_ice3(:,:,k); enddo 
   do k =  1, kmz ; WRITE(300) f_rain3(:,:,k); enddo
   do k =  1, kmz ; WRITE(300) f_rimef3(:,:,k); enddo
  
   nvars = 31 ; dxx = dlon/3. ; dyy = dlat/3. 
   imr = im3 ; jmr = jm3 ; ims = im3 ; jms = (nint(aspect*ims)/2)*2 + 1 
   
   DeAllocate( uu3, vv3, ww3, tt3, td3, qt3, pp3, qq3, cw3, rh3, wks) 
   DeAllocate( difu3, difv3, km_diffus3,kh_diffus3 )
   DeAllocate( hdifu3, hdifv3)
   DeAllocate( f_ice3,f_rain3,f_rimef3) 
   DeAllocate(                mc_heat3)
   
   ENDIF !* ( refine ) - - - - - - - - - - - - - - - - - - - - - - - 

!           *     *     *     *     *     *     *     *     *     * 
   
!* GrADS control file for height-level data on rotated grid 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   xxr1 = .5*(1-imr)*dxx ; xxr_ = dxx !* Lagrangian coord for moving nest 
   yyr1 = .5*(1-jmr)*dyy ; yyr_ = dyy !* and for fixed outer domain (d01) 
   
   WRITE(301,'("                                                      ")') 
   WRITE(301,'(" DSET    height_lev.dat                               ")') 
   WRITE(301,'(" TITLE   Height-level data on rotated latlon grid     ")') 
   WRITE(301,'(" OPTIONS sequential big_endian                        ")') 
   WRITE(301,'(" UNDEF",E11.2)') undef 
   WRITE(301,'("                                                      ")') 
   WRITE(301,'(" XDEF", I7,"   LINEAR", 2F11.5 )') imr, xxr1, xxr_ 
   WRITE(301,'(" YDEF", I7,"   LINEAR", 2F11.5 )') jmr, yyr1, yyr_ 
   WRITE(301,'(" ZDEF", I7,"   LEVELS")') kmz 
   WRITE(301,*)  zlevel 
   WRITE(301,'(" TDEF", I7,"   LINEAR   00:00Z01JAN0001  1yr")') ntime 
   WRITE(301,'("                                                      ")') 
   WRITE(301,'(" VARS", I3)') nvars
   WRITE(301,'(" xx        0   0    Earth longitude (deg)             ")') 
   WRITE(301,'(" yy        0   0    Earth latitude  (deg)             ")') 
   WRITE(301,'(" zs        0   0    Earth topography (m)              ")') 
   WRITE(301,'(" ps        0   0    surface air   pressure (mb)       ")') 
   WRITE(301,'(" ts        0   0    surface air temperature (K)       ")') 
   WRITE(301,'(" qs        0   0    surface specific humidity (kg/kg) ")') 
   WRITE(301,'(" p0        0   0    mean-sea-level pressure (mb)      ")') 
   WRITE(301,'(" prc       0   0    precipitation (m)                 ")') 
   WRITE(301,'(" u10       0   0    rotated 10m x-wind (m/s)          ")') 
   WRITE(301,'(" v10       0   0    rotated 10m y-wind (m/s)          ")') 
   WRITE(301,'(" t10       0   0    10m air temperature (K)           ")') 
   WRITE(301,'(" q10       0   0    10m specific humidity (kg/kg)     ")') 
   WRITE(301,'(" sst       0   0    sea-surface temperature (K)       ")') 
   WRITE(301,'(" ssq       0   0    sea-sfc specific humidity (kg/kg) ")') 
   WRITE(301,'(" tw        0   0    sensible heat flux (W/m^2)        ")') 
   WRITE(301,'(" qw        0   0    latent   heat flux (W/m^2)        ")') 
   WRITE(301,'(" uf        0   0    friction velocity  (m/s)          ")') 
   WRITE(301,'(" z0        0   0    roughness height   (m)            ")') 
   WRITE(301,'(" uu  ", I7,"   0    rotated x-wind (m/s)              ")') kmz 
   WRITE(301,'(" vv  ", I7,"   0    rotated y-wind (m/s)              ")') kmz 
   WRITE(301,'(" ww  ", I7,"   0    vertical  wind (m/s)              ")') kmz 
   WRITE(301,'(" tt  ", I7,"   0    air temperature (K)               ")') kmz 
   WRITE(301,'(" vx  ", I7,"   0    vorticity (/s)          ")') kmz 
   WRITE(301,'(" qt  ", I7,"   0    equiv. pot. temperature (K)       ")') kmz 
   WRITE(301,'(" pp  ", I7,"   0    air  pressure  (mb)               ")') kmz 
   WRITE(301,'(" qq  ", I7,"   0    specific humidity (kg/kg)         ")') kmz 
   WRITE(301,'(" cw  ", I7,"   0    total condensate  (kg/kg)         ")') kmz 
   WRITE(301,'(" rh  ", I7,"   0    relative humidity (%)             ")') kmz 
!   WRITE(301,'("tpert", I7,"   0    temperature perturbation (%)      ")') kmz
   WRITE(301,'("f_ice", I7,"   0     ice fraction      ")') kmz
   WRITE(301,'("f_rain", I7,"   0    rain fraction      ")') kmz
   WRITE(301,'("f_rimef", I7,"   0   rime factor      ")') kmz
   WRITE(301,'(" ENDVARS                                              ")') 
   WRITE(301,'("                                                      ")') 
   
!* GrADS control file for height-level data on standard grid 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   xxs1 = .5*(1-ims)*dxx + ctr_lon ; xxs_ = dxx !* d01 west  bdry 
   yys1 = .5*(1-jms)*dyy + ctr_lat ; yys_ = dyy !* d01 south bdry 
   
   WRITE(302,'("                                                      ")') 
   WRITE(302,'(" DSET    height_lev.dat                               ")') 
   WRITE(302,'(" TITLE   Height-level data on standard latlon grid    ")') 
   WRITE(302,'(" OPTIONS sequential big_endian                        ")') 
   WRITE(302,'(" UNDEF",E11.2)') undef 
   WRITE(302,'("                                                      ")') 
   WRITE(302,'(" PDEF", I7,  I6,   A10,   4F11.5)') & 
                       imr, jmr, 'ETA.U', ctr_lon, ctr_lat, dxx,  dyy 
   WRITE(302,'(" XDEF", I7,"   LINEAR", 2F11.5 )') ims, xxs1, xxs_ 
   WRITE(302,'(" YDEF", I7,"   LINEAR", 2F11.5 )') jms, yys1, yys_ 
   WRITE(302,'(" ZDEF", I7,"   LEVELS")') kmz 
   WRITE(302,*)  zlevel 
   WRITE(302,'(" TDEF", I7,"   LINEAR   00:00Z01JAN0001  1yr")') ntime 
   WRITE(302,'("                                                      ")') 
   WRITE(302,'(" VECTORPAIRS   u10,v10  uu,vv                         ")') 
   WRITE(302,'(" VARS", I3)') nvars
   WRITE(302,'(" xx        0   0    Earth longitude (deg)             ")') 
   WRITE(302,'(" yy        0   0    Earth latitude  (deg)             ")') 
   WRITE(302,'(" zs        0   0    Earth topography (m)              ")') 
   WRITE(302,'(" ps        0   0    surface air   pressure (mb)       ")') 
   WRITE(302,'(" ts        0   0    surface air temperature (K)       ")') 
   WRITE(302,'(" qs        0   0    surface specific humidity (kg/kg) ")') 
   WRITE(302,'(" p0        0   0    mean-sea-level pressure (mb)      ")') 
   WRITE(302,'(" prc       0   0    precipitation (m)                 ")') 
   WRITE(302,'(" u10       0   0    10m zonal wind comp (m/s)         ")') 
   WRITE(302,'(" v10       0   0    10m meridional wind (m/s)         ")') 
   WRITE(302,'(" t10       0   0    10m air temperature (K)           ")') 
   WRITE(302,'(" q10       0   0    10m specific humidity (kg/kg)     ")') 
   WRITE(302,'(" sst       0   0    sea-surface temperature (K)       ")') 
   WRITE(302,'(" ssq       0   0    sea-sfc specific humidity (kg/kg) ")') 
   WRITE(302,'(" tw        0   0    sensible heat flux (W/m^2)        ")') 
   WRITE(302,'(" qw        0   0    latent   heat flux (W/m^2)        ")') 
   WRITE(302,'(" uf        0   0    friction velocity  (m/s)          ")') 
   WRITE(302,'(" z0        0   0    roughness height   (m)            ")') 
   WRITE(302,'(" uu  ", I7,"   0    rotated x-wind (m/s)              ")') kmz 
   WRITE(302,'(" vv  ", I7,"   0    rotated y-wind (m/s)              ")') kmz 
   WRITE(302,'(" ww  ", I7,"   0    vertical  wind (m/s)              ")') kmz 
   WRITE(302,'(" tt  ", I7,"   0    air temperature (K)               ")') kmz 
   WRITE(302,'(" vx  ", I7,"   0    vorticity (/s)          ")') kmz 
   WRITE(302,'(" qt  ", I7,"   0    equiv. pot. temperature (K)       ")') kmz 
   WRITE(302,'(" pp  ", I7,"   0    air  pressure  (mb)               ")') kmz 
   WRITE(302,'(" qq  ", I7,"   0    specific humidity (kg/kg)         ")') kmz 
   WRITE(302,'(" cw  ", I7,"   0    total condensate  (kg/kg)         ")') kmz 
   WRITE(302,'(" rh  ", I7,"   0    relative humidity (%)             ")') kmz 
!   WRITE(302,'("tpert", I7,"   0    temperature perturbation (%)      ")') kmz
   WRITE(302,'("f_ice", I7,"   0     ice fraction      ")') kmz
   WRITE(302,'("f_rain", I7,"   0    rain fraction      ")') kmz
   WRITE(302,'("f_rimef", I7,"   0   rime factor      ")') kmz
   WRITE(302,'(" ENDVARS                                              ")') 
   WRITE(302,'("                                                      ")') 
   
   CLOSE(300) ; CLOSE(301) ; CLOSE(302) 
   
   ENDIF !* ( height_level ) ------------------------------------------- 

!********************************************! 
!  3a. Height-level data on cylindrical grid ! 
!____________________________________________! 

   IF ( height_cylin ) THEN !* ----------------------------------------- 

!* Variables on cylindrical grid with south-polar projection 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   jmc = nint(r_cylin/dlat) !*   radial  dim: j=0 -> center of cylinder 
   imc = nint(azimuth*jmc)  !* azimuthal dim: i=0 -> south; + clockwise 
   dlamda = 2*pi_cnst/imc   !* azimuthal resolution  ( dr = yrlen ) 
   k2d = 20                 !* total # of 2D fields 
   
   Allocate( icen(0:kmz), jcen(0:kmz), rcos(0:imc), rsin(0:imc) ) 
   Allocate( ppc(0:imc,0:jmc,kmz), qqc(0:imc,0:jmc,kmz), cwc(0:imc,0:jmc,kmz) ) 
   Allocate( ttc(0:imc,0:jmc,kmz), tdc(0:imc,0:jmc,kmz), qtc(0:imc,0:jmc,kmz) ) 
   Allocate( uuc(0:imc,0:jmc,kmz), vvc(0:imc,0:jmc,kmz), wwc(0:imc,0:jmc,kmz) ) 
   Allocate( uxc(0:imc,0:jmc,kmz), vxc(0:imc,0:jmc,kmz), mwc(0:imc,0:jmc,kmz), qwc(0:imc,0:jmc,kmz)) 
   Allocate( rhc(0:imc,0:jmc,kmz), wks(0:imc,0:jmc,k2d) ) 
   Allocate( f_icec(0:imc,0:jmc,kmz),f_rainc(0:imc,0:jmc,kmz),f_rimefc(0:imc,0:jmc,kmz))
   Allocate(                                                  mc_heatc(0:imc,0:jmc,kmz))
!
!* Advection quantities (gopal's doing)
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   Allocate( lfluxz(0:jmc),hfluxz(0:jmc),warmxz(0:jmc,kmz),uvortxz(0:jmc,kmz)) 
   Allocate( xstress(0:jmc),ystress(0:jmc),stress(0:jmc) )
   Allocate( u1d(kmz), v1d(kmz), t1d(kmz), rad1d(kmz))
   Allocate( uu2d(0:imc,0:jmc,kmz),vv2d(0:imc,0:jmc,kmz),ww2d(0:imc,0:jmc,kmz),pp2d(0:imc,0:jmc,kmz))
   Allocate( tt2d(0:imc,0:jmc,0:jmc),vort2d(0:imc,0:jmc,kmz),th2d(0:imc,0:jmc,kmz))
   Allocate( upert(0:imc,0:jmc,kmz),vpert(0:imc,0:jmc,kmz),wpert(0:imc,0:jmc,kmz),vopert(0:imc,0:jmc,kmz),thpert(0:imc,0:jmc,kmz))
   Allocate( dvdz(0:imc,0:jmc,kmz),edrvf(0:imc,0:jmc,kmz),eddvdz(0:imc,0:jmc,kmz))
   Allocate( wstat2d(0:imc,0:jmc),tstat2d(0:imc,0:jmc))
   Allocate( mean_var1(0:imc,0:jmc,kmz),rms_var1(0:imc,0:jmc,kmz),wstat(0:imc,0:jmc,kmz))
   Allocate( mean_var2(0:imc,0:jmc,kmz),rms_var2(0:imc,0:jmc,kmz),tstat(0:imc,0:jmc,kmz))
   Allocate( uadc(0:imc,0:jmc,kmz), vadc(0:imc,0:jmc,kmz), tadc(0:imc,0:jmc,kmz))
   Allocate(  udiffc(0:imc,0:jmc,kmz),  vdiffc(0:imc,0:jmc,kmz))
   Allocate( hudiffc(0:imc,0:jmc,kmz), hvdiffc(0:imc,0:jmc,kmz))
   Allocate( tan_grd_wind(0:imc,0:jmc,kmz))
   Allocate( tpert(0:imc,0:jmc,kmz), vv_verify(0:imc,0:jmc,kmz))
   Allocate( pgr(0:imc,0:jmc,kmz), vgr(0:imc,0:jmc,kmz), ror(0:imc,0:jmc,kmz))
   Allocate( kmc(0:imc,0:jmc,kmz), khc(0:imc,0:jmc,kmz))
!
!* For SEE computation
!* ~~~~~~~~~~~~~~~~~~~
   Allocate(zdiff2d(0:imc,0:jmc,kmz),hdiff2d(0:imc,0:jmc,kmz)) 
   Allocate(thdot(0:imc,0:jmc,kmz),edth1(0:imc,0:jmc,kmz),edth2(0:imc,0:jmc,kmz),edth3(0:imc,0:jmc,kmz))
   Allocate(fthdiff(0:imc,0:jmc,kmz),fmomdiff(0:imc,0:jmc,kmz))

   ppc(:,:,:) = undef ; qqc(:,:,:) = undef ; cwc(:,:,:) = undef 
   ttc(:,:,:) = undef ; tdc(:,:,:) = undef ; qtc(:,:,:) = undef 
   uuc(:,:,:) = undef ; vvc(:,:,:) = undef ; wwc(:,:,:) = undef 
   uxc(:,:,:) = undef ; vxc(:,:,:) = undef ; mwc(:,:,:) = undef ; qwc(:,:,:) = undef 
   rhc(:,:,:) = undef ; wks(:,:,:) = undef 
   f_icec(:,:,:) = undef; f_rainc(:,:,:) = undef; f_rimefc(:,:,:) = undef
                                                  mc_heatc(:,:,:) = undef
!
!* Advection quantities (gopal's doing)
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
   lfluxz(:) = undef; hfluxz(:) = undef; warmxz(:,:) = undef; uvortxz(:,:) = undef 
   xstress(:)=undef; ystress(:)=undef; stress(:)=undef
   u1d(:) = undef; v1d(:) = undef; t1d(:) = undef; rad1d(:) = undef 
   uu2d(:,:,:) =undef; vv2d(:,:,:) =undef; ww2d(:,:,:) =undef; pp2d(:,:,:) =undef
   tt2d(:,:,:) =undef; vort2d(:,:,:) =undef; th2d(:,:,:) =undef
   dvdz(:,:,:) =undef; edrvf(:,:,:) =undef; eddvdz(:,:,:) =undef
   upert(:,:,:)=undef;  vpert(:,:,:)=undef; wpert(:,:,:)=undef; vopert(:,:,:)=undef;thpert(:,:,:)=undef
   mean_var1(:,:,:) =undef; rms_var1(:,:,:) =undef; wstat2d(:,:) =undef; wstat(:,:,:) =undef
   mean_var2(:,:,:) =undef; rms_var2(:,:,:) =undef; tstat2d(:,:) =undef; tstat(:,:,:) =undef
   uadc(:,:,:) = undef ; vadc(:,:,:) = undef; tadc(:,:,:) = undef
   udiffc(:,:,:) = undef ;  vdiffc(:,:,:) = undef; kmc(:,:,:) = undef ; khc(:,:,:) = undef
   hudiffc(:,:,:) = undef ; hvdiffc(:,:,:) = undef
   tan_grd_wind(:,:,:) = undef
   tpert(:,:,:) = undef ; vv_verify(:,:,:) = undef
   pgr(:,:,:) = undef ; vgr(:,:,:) = undef ; ror(:,:,:) = undef

!
!* For SEE computation
!* ~~~~~~~~~~~~~~~~~~~
   zdiff2d(:,:,:)=undef;hdiff2d(:,:,:)=undef
   thdot(:,:,:)=undef;edth1(:,:,:)=undef;edth2(:,:,:)=undef;edth3(:,:,:)=undef 
   fthdiff(:,:,:)=undef;fmomdiff(:,:,:)=undef
   
!
!           *     *     *     *     *     *     *     *     *     * 
   
!* 3D Storm center for tilted cylinder axis 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   icen(:) = iloc ; jcen(:) = jloc !* Storm center by MSLP 
   
   IF ( tilt_cylin ) THEN  !* 3D Storm center by vorticity - - - - - 

   ioff = nint(.5*r_storm/dlon) ; ibd = 3 
   joff = nint(.5*r_storm/dlat) ; jbd = 3 
   iold = iloc  ;  jold = jloc  ; zzz = 0. 
 
   DO k = 1, kmz 
      aaa = 0 ; bbb = 0 ; eee = 0 
      i1 = icen(k-1) - ioff ; i1 = max( 1+ibd, i1  ) 
      i2 = icen(k-1) + ioff ; i2 = min( i2, im-ibd ) 
      j1 = jcen(k-1) - joff ; j1 = max( 1+jbd, j1  ) 
      j2 = jcen(k-1) + joff ; j2 = min( j2, jm-jbd ) 
      DO j = j1, j2 ; DO i = i1, i2 
         www = max( vx_(i,j,k), 0. ) 
         IF ( www == undef ) www = 0. 
         aaa = aaa + www*i ; bbb = bbb + www*j ; eee = eee + www 
      ENDDO ; ENDDO
      IF ( eee == 0. ) THEN 
         icen(k) = icen(k-1) ; jcen(k) = jcen(k-1) 
      ELSE 
         itmp = nint(aaa/eee) ; jtmp = nint(bbb/eee) 
          m0  = int( (zlevel(k)-zzz)*r_slope/yrlen ) 
         icen(k) = max( iold-m0, min(itmp, iold+m0) ) 
         jcen(k) = max( jold-m0, min(jtmp, jold+m0) ) 
         IF ( icen(k)/=icen(k-1) .OR. jcen(k)/=jcen(k-1) ) THEN 
            iold = icen(k) ; jold = jcen(k) ; zzz = zlevel(k) 
         ENDIF
      ENDIF
   ENDDO 
 
!  PRINT*, '------------------------------------------------------------' 
!  DO k = 0, kmz ; PRINT*, 'k, icen, jcen =', k, icen(k), jcen(k) ; ENDDO 
!  PRINT*, '------------------------------------------------------------' 
 
   ENDIF !* ( tilt_cylin ) - - - - - - - - - - - - - - - - - - - - - 

!* Rotated latlon to cylindrical coord transform 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!  outward radial wind: uuc =   rcos*uur + rsin*vvr 
!  CCW tangential wind: vvc = - rsin*uur + rcos*vvr 
   
   DO i = 0, imc !* South-Polar Projection 
       xxx = 1.5*pi_cnst - i*dlamda  
      rcos(i) = cos(xxx) ; rsin(i) = sin(xxx) 
   ENDDO 

!* Remove Mean flow for further analysis (works for 1 hourly or less outputs)  
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!
   IF ( domain_num > 1 .and. env==1) THEN
     OPEN(52,file='easy_track.csv',status='old')
      DO mkont=0,ifcst_hour
       READ(52,*)kkont,atemp,btemp
       lon1(mkont)=atemp
       lat1(mkont)=btemp
      ENDDO
      IF(ifcst_hour .EQ. 0)WRITE(0,*)'forecast hour,    uenv (m/s),   venv (m/s)'
      IF(ifcst_hour .GE. 1)THEN
       WRITE(0,*)lon1(ifcst_hour),lon1(ifcst_hour-1),lat1(ifcst_hour),lat1(ifcst_hour-1)
       uenv=(lon1(ifcst_hour)-lon1(ifcst_hour-1))*111130/process_freq_sec
       venv=(lat1(ifcst_hour)-lat1(ifcst_hour-1))*111130/process_freq_sec 
       WRITE(57,*)ifcst_hour,ifcst_minute,uenv,venv
      ENDIF 
     CLOSE(52)

     DO j = 1, jm
      DO i = 1, im
        u10(i,j)=u10(i,j)-uenv     ! u10
        v10(i,j)=v10(i,j)-venv     ! v10
      ENDDO
     ENDDO

     DO k = 1, kmz
      DO j = 1, jm
       DO i = 1, im
        uu_(i,j,k)=uu_(i,j,k)-uenv   ! uu in height coordinate
        vv_(i,j,k)=vv_(i,j,k)-venv   ! vv in height coordinate
       ENDDO
      ENDDO
     ENDDO

   ENDIF ! (this code is commented for basin scale)


! Store Coriolis at the domain center for relative vorticity computation

   fcor=2*o_earth*Sin(yrs(iloc,jloc)*deg2rad) !* Coriolis parameter 

   WRITE(0,*)
   WRITE(0,*)'-------- Coriolis at the center of the domain --------'
   WRITE(0,*)fcor
   WRITE(0,*)'------------------------------------------------------'
   WRITE(0,*)

!* Interpolate variables to cylindrical grid 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   ibd = 3 ; i1 = 1 + ibd ; i2 = im - ibd 
   jbd = 3 ; j1 = 1 + jbd ; j2 = jm - jbd 
   bad_cylin = .false.    ; kout = km + 1 

   DO j = 1, jmc 
   DO i = 1, imc ; k = 0 
      xxx = j*rcos(i) ; iloc = int(xxx) ; IF(xxx<0)iloc=iloc-1; aaa = xxx - iloc ; bbb = 1 - aaa
      yyy = j*rsin(i) ; jloc = int(yyy) ; IF(yyy<0)jloc=jloc-1; ccc = yyy - jloc ; ddd = 1 - ccc

!        - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
         m = icen(k) + iloc ; m1 = m + 1 ; n = jcen(k) + jloc ; n1 = n + 1 
         IF ( (i1 <= m).AND.(m1 <= i2).AND.(j1 <= n).AND.(n1 <= j2) ) THEN 
            wks(i,j, 1) = bbb*( ddd*xrs(m ,n) + ccc*xrs(m ,n1) ) &
                        + aaa*( ddd*xrs(m1,n) + ccc*xrs(m1,n1) ) 
            wks(i,j, 2) = bbb*( ddd*yrs(m ,n) + ccc*yrs(m ,n1) ) &
                        + aaa*( ddd*yrs(m1,n) + ccc*yrs(m1,n1) ) 
            wks(i,j, 3) = bbb*( ddd* zs(m ,n) + ccc* zs(m ,n1) ) &
                        + aaa*( ddd* zs(m1,n) + ccc* zs(m1,n1) ) 
            wks(i,j, 4) = bbb*( ddd* ps(m ,n) + ccc* ps(m ,n1) ) &
                        + aaa*( ddd* ps(m1,n) + ccc* ps(m1,n1) ) 
            wks(i,j, 5) = bbb*( ddd* ts(m ,n) + ccc* ts(m ,n1) ) &
                        + aaa*( ddd* ts(m1,n) + ccc* ts(m1,n1) ) 
            wks(i,j, 6) = bbb*( ddd* qs(m ,n) + ccc* qs(m ,n1) ) &
                        + aaa*( ddd* qs(m1,n) + ccc* qs(m1,n1) ) 
            wks(i,j, 7) = bbb*( ddd* p0(m ,n) + ccc* p0(m ,n1) ) &
                        + aaa*( ddd* p0(m1,n) + ccc* p0(m1,n1) ) 
            wks(i,j, 8) = bbb*( ddd*prc(m ,n) + ccc*prc(m ,n1) ) &
                        + aaa*( ddd*prc(m1,n) + ccc*prc(m1,n1) ) 
                   uur  = bbb*( ddd*u10(m ,n) + ccc*u10(m ,n1) ) &
                        + aaa*( ddd*u10(m1,n) + ccc*u10(m1,n1) ) 
                   vvr  = bbb*( ddd*v10(m ,n) + ccc*v10(m ,n1) ) &
                        + aaa*( ddd*v10(m1,n) + ccc*v10(m1,n1) ) 
            wks(i,j, 9) =   rcos(i)*uur + rsin(i)*vvr  !*   radial   wind 
            wks(i,j,10) = - rsin(i)*uur + rcos(i)*vvr  !* tangential wind 
            wks(i,j,11) = bbb*( ddd*t10(m ,n) + ccc*t10(m ,n1) ) &
                        + aaa*( ddd*t10(m1,n) + ccc*t10(m1,n1) ) 
            wks(i,j,12) = bbb*( ddd*q10(m ,n) + ccc*q10(m ,n1) ) &
                        + aaa*( ddd*q10(m1,n) + ccc*q10(m1,n1) ) 
            IF ( max(sst(m,n),sst(m,n1),sst(m1,n),sst(m1,n1)) < undef ) THEN 
            wks(i,j,13) = bbb*( ddd*sst(m ,n) + ccc*sst(m ,n1) ) &
                        + aaa*( ddd*sst(m1,n) + ccc*sst(m1,n1) ) 
            wks(i,j,14) = bbb*( ddd*ssq(m ,n) + ccc*ssq(m ,n1) ) &
                        + aaa*( ddd*ssq(m1,n) + ccc*ssq(m1,n1) ) 
            ENDIF 
            wks(i,j,15) = bbb*( ddd* tw(m ,n) + ccc* tw(m ,n1) ) &
                        + aaa*( ddd* tw(m1,n) + ccc* tw(m1,n1) ) 
            wks(i,j,16) = bbb*( ddd* qw(m ,n) + ccc* qw(m ,n1) ) &
                        + aaa*( ddd* qw(m1,n) + ccc* qw(m1,n1) ) 
            wks(i,j,17) = bbb*( ddd* uf(m ,n) + ccc* uf(m ,n1) ) &
                        + aaa*( ddd* uf(m1,n) + ccc* uf(m1,n1) ) 
            wks(i,j,18) = bbb*( ddd* z0(m ,n) + ccc* z0(m ,n1) ) &
                        + aaa*( ddd* z0(m1,n) + ccc* z0(m1,n1) ) 

!           taux and tauy (frictional stresses)  ! gopal's doing

                   uur  = bbb*( ddd*tx(m ,n) + ccc*tx(m ,n1) ) &
                        + aaa*( ddd*tx(m1,n) + ccc*tx(m1,n1) )
                   vvr  = bbb*( ddd*ty(m ,n) + ccc*ty(m ,n1) ) &
                        + aaa*( ddd*ty(m1,n) + ccc*ty(m1,n1) )
            wks(i,j,19) =   rcos(i)*uur + rsin(i)*vvr  !*   radial   wind 
            wks(i,j,20) = - rsin(i)*uur + rcos(i)*vvr  !* tangential wind

         ENDIF 
!        - - - - - - - Interpolate 2D fields above / 3D fields below - - - 
      DO k = 1, kmz 
         m = icen(k) + iloc ; m1 = m + 1 ; n = jcen(k) + jloc ; n1 = n + 1 
         defined = (i1 <= m).AND.(m1 <= i2) .AND. (j1 <= n).AND.(n1 <= j2) 
         IF ( .not. defined ) THEN 
            bad_cylin = .true. ; kout = min(k,kout) ; WRITE(0,*)'BAD CYLINDER',i,j,k;EXIT 
         ENDIF 
         eee = max( pp_(m,n,k), pp_(m,n1,k), pp_(m1,n,k), pp_(m1,n1,k) ) 
         IF ( eee < undef ) THEN !* ------------------------------------ 
            ppc(i,j,k) = bbb*( ddd*pp_(m ,n,k) + ccc*pp_(m ,n1,k) ) & 
                       + aaa*( ddd*pp_(m1,n,k) + ccc*pp_(m1,n1,k) ) 
            qqc(i,j,k) = bbb*( ddd*qq_(m ,n,k) + ccc*qq_(m ,n1,k) ) & 
                       + aaa*( ddd*qq_(m1,n,k) + ccc*qq_(m1,n1,k) ) 
            cwc(i,j,k) = bbb*( ddd*cw_(m ,n,k) + ccc*cw_(m ,n1,k) ) & 
                       + aaa*( ddd*cw_(m1,n,k) + ccc*cw_(m1,n1,k) ) 
            f_icec(i,j,k) = bbb*( ddd*f_ice_(m ,n,k) + ccc*f_ice_(m ,n1,k) ) &
                       + aaa*( ddd*f_ice_(m1,n,k) + ccc*f_ice_(m1,n1,k) )
            f_rainc(i,j,k) = bbb*( ddd*f_rain_(m ,n,k) + ccc*f_rain_(m ,n1,k) ) &
                       + aaa*( ddd*f_rain_(m1,n,k) + ccc*f_rain_(m1,n1,k) )
            f_rimefc(i,j,k) = bbb*( ddd*f_rimef_(m ,n,k) + ccc*f_rimef_(m ,n1,k) ) &
                       + aaa*( ddd*f_rimef_(m1,n,k) + ccc*f_rimef_(m1,n1,k) )
            mc_heatc(i,j,k) = bbb*( ddd*mc_heat_(m ,n,k) + ccc*mc_heat_(m ,n1,k)) &
                       + aaa*( ddd*mc_heat_(m1,n,k) + ccc*mc_heat_(m1,n1,k) )
            ttc(i,j,k) = bbb*( ddd*tt_(m ,n,k) + ccc*tt_(m ,n1,k) ) & 
                       + aaa*( ddd*tt_(m1,n,k) + ccc*tt_(m1,n1,k) ) 
            tdc(i,j,k) = bbb*( ddd*td_(m ,n,k) + ccc*td_(m ,n1,k) ) & 
                       + aaa*( ddd*td_(m1,n,k) + ccc*td_(m1,n1,k) ) 
            rhc(i,j,k) = bbb*( ddd*rh_(m ,n,k) + ccc*rh_(m ,n1,k) ) & 
                       + aaa*( ddd*rh_(m1,n,k) + ccc*rh_(m1,n1,k) ) 
            qtc(i,j,k) = bbb*( ddd*qt_(m ,n,k) + ccc*qt_(m ,n1,k) ) & 
                       + aaa*( ddd*qt_(m1,n,k) + ccc*qt_(m1,n1,k) ) 
            uxc(i,j,k) = bbb*( ddd*ux_(m ,n,k) + ccc*ux_(m ,n1,k) ) & 
                       + aaa*( ddd*ux_(m1,n,k) + ccc*ux_(m1,n1,k) ) 
            vxc(i,j,k) = bbb*( ddd*vx_(m ,n,k) + ccc*vx_(m ,n1,k) ) & 
                       + aaa*( ddd*vx_(m1,n,k) + ccc*vx_(m1,n1,k) ) 
            wwc(i,j,k) = bbb*( ddd*ww_(m ,n,k) + ccc*ww_(m ,n1,k) ) &
                       + aaa*( ddd*ww_(m1,n,k) + ccc*ww_(m1,n1,k) )

!           IF(wwc(i,j,k) .GE. 30.)THEN
!               WRITE(0,*)'LARGE W'
!               WRITE(0,*)xxx,iloc,yyy,jloc
!               WRITE(0,*)aaa,bbb,ccc,ddd
!               WRITE(0,*)bbb,( ddd*ww_(m,n,k) + ccc*ww_(m,n1,k) )
!               WRITE(0,*)aaa,( ddd*ww_(m1,n,k)+ ccc*ww_(m1,n1,k) )
!               WRITE(0,*)m,n,n1,k,ww_(m ,n,k),ww_(m ,n1,k)
!               WRITE(0,*)i,j,k,'                   ',wwc(i,j,k)
!               WRITE(0,*)m1,n,n1,k,ww_(m1,n,k),ww_(m1,n1,k)
!            ENDIF

                  uur  = bbb*( ddd*uu_(m ,n,k) + ccc*uu_(m ,n1,k) ) & 
                       + aaa*( ddd*uu_(m1,n,k) + ccc*uu_(m1,n1,k) ) 
                  vvr  = bbb*( ddd*vv_(m ,n,k) + ccc*vv_(m ,n1,k) ) & 
                       + aaa*( ddd*vv_(m1,n,k) + ccc*vv_(m1,n1,k) ) 
            uuc(i,j,k) =   rcos(i)*uur + rsin(i)*vvr !*   radial   wind 
            vvc(i,j,k) = - rsin(i)*uur + rcos(i)*vvr !* tangential wind 
!           - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            mwc(i,j,k) =   wwc(i,j,k)*ppc(i,j,k)*100. & 
                       / ( rd_air*(1.+.608*qqc(i,j,k))*ttc(i,j,k) )
!
! gopal's doing for the budget computation
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
            qwc(i,j,k) =  1000.*qqc(i,j,k)*mwc(i,j,k)    
        
            kmc(i,j,k) = bbb*( ddd*km_diffus_(m ,n,k) + ccc*km_diffus_(m ,n1,k) ) &
                       + aaa*( ddd*km_diffus_(m1,n,k) + ccc*km_diffus_(m1,n1,k) ) 

            khc(i,j,k) = bbb*( ddd*kh_diffus_(m ,n,k) + ccc*kh_diffus_(m ,n1,k) ) &
                       + aaa*( ddd*kh_diffus_(m1,n,k) + ccc*kh_diffus_(m1,n1,k) )

            uur  = bbb*( ddd*difu_(m ,n,k) + ccc*difu_(m ,n1,k) ) &
                 + aaa*( ddd*difu_(m1,n,k) + ccc*difu_(m1,n1,k) )
            vvr  = bbb*( ddd*difv_(m ,n,k) + ccc*difv_(m ,n1,k) ) &
                 + aaa*( ddd*difv_(m1,n,k) + ccc*difv_(m1,n1,k) )  

            udiffc(i,j,k) =   rcos(i)*uur + rsin(i)*vvr !*   radial vertical diffusion 
            vdiffc(i,j,k) = - rsin(i)*uur + rcos(i)*vvr !*   tangential vertical diffusion

            huur  = bbb*( ddd*hdifu_(m ,n,k) + ccc*hdifu_(m ,n1,k) ) &
                  + aaa*( ddd*hdifu_(m1,n,k) + ccc*hdifu_(m1,n1,k) )
            hvvr  = bbb*( ddd*hdifv_(m ,n,k) + ccc*hdifv_(m ,n1,k) ) &
                  + aaa*( ddd*hdifv_(m1,n,k) + ccc*hdifv_(m1,n1,k) )

            hudiffc(i,j,k) =   rcos(i)*huur + rsin(i)*hvvr !*   radial horizontal diffusion 
            hvdiffc(i,j,k) = - rsin(i)*huur + rcos(i)*hvvr !*   tangential horizontal diffusion
!
! any species based on the above computations can be add as new variables here Bob!
!

         ELSE
          WRITE(0,*)'PROBLEMS WITH HORIZONTAL INTERPOLATION',I,J,K
          WRITE(0,*)pp_(m,n,k),pp_(m,n1,k),pp_(m1,n,k),pp_(m1,n1,k)
         ENDIF !* ------------------------------------------------------ 
      ENDDO 
   ENDDO 
   ENDDO 

!
!  gopal's changes to include radial  and tangential advection
!  dlamda = 2*pi_cnst/imc   !* azimuthal resolution  ( dr = yrlen )
!
   t1d=0.0
   do k=1,kmz
!   t1d(k) = qtc(1,jmc,k)    !thetae pert change
    t1d(k) = ttc(1,jmc,k)
   enddo
!
!  Compute domain average profiles of u and v 
!
!    write(0,*)'---------------------------------------------------------'
!    write(0,*)'Domain Averaged Profiles of Radial and Tangential Winds'
    do k = 1, kmz-2
     utemp=0.0
     vtemp=0.0
     do j = 1, jmc-1
      do i = 1, imc-1
        utemp = utemp + uuc(i,j,k)
        vtemp = vtemp + vvc(i,j,k)
      enddo
     enddo
        u1d(k)=utemp/((imc-1)*(jmc-1))
        v1d(k)=vtemp/((imc-1)*(jmc-1))
!       write(0,*)k,u1d(k),v1d(k)
    enddo
!   write(0,*)'---------------------------------------------------------'

!
!  Compute temperature perturbation (local and wrt mesoscale environment) 
!

   tpert=0.0
   do k = 1, kmz-2
    do j = 1, jmc-1
     do i = 1, imc-1
       tpert(i,j,k) = ttc(i,j,k) - t1d(k)
     enddo
    enddo
   enddo

!
!  Compute azimuthally averaged wind, thermodynamic and vorticity fields 
!
   uu2d=0.0
   vv2d=0.0
   ww2d=0.0
   tt2d=0.0
   pp2d=0.0
   vort2d=0.0
   th2d=0.0
   do k = 1, kmz-2
    do j = 1, jmc
     atemp=0.0
     btemp=0.0
     ctemp=0.0
     dtemp=0.0
     etemp=0.0
     ftemp=0.0
     gtemp=0.0
     do i = 1, imc-1
      atemp = atemp + uuc(i,j,k)
      btemp = btemp + vvc(i,j,k)
      ctemp = ctemp + wwc(i,j,k)
      dtemp = dtemp + ppc(i,j,k)
      etemp = etemp + ttc(i,j,k)
      ftemp = ftemp + vxc(i,j,k)
      gtemp = gtemp + qtc(i,j,k)
     enddo
     do i = 1, imc-1 ! kept as 3d variable for plotting purpose
      uu2d(i,j,k)   = atemp/(imc-1)
      vv2d(i,j,k)   = btemp/(imc-1)
      ww2d(i,j,k)   = ctemp/(imc-1)
      pp2d(i,j,k)   = dtemp/(imc-1)
      tt2d(i,j,k)   = etemp/(imc-1)
      vort2d(i,j,k) = ftemp/(imc-1) + fcor   ! Absolute vorticity
      th2d(i,j,k)   = gtemp/(imc-1)          ! Potential temperature for SEE
     enddo
    enddo
   enddo

!  Terms for mean budgets. Mean and Eddy flux terms are computed here

   dvdz=0.0
   do k = 2, kmz-2
    do j = 1, jmc
     do i = 1, imc-1
      dvdz(i,j,k)=0.5*(vv2d(i,j,k+1)-vv2d(i,j,k-1))/(zlevel(k+1)-zlevel(k))
     enddo
    enddo
   enddo

!  Compute perturbation of tangential wind, and w

   upert=0.0;vpert=0.0;wpert=0.0;vopert=0.0;thpert=0.0
   do k = 1, kmz-2
    do j = 1, jmc
     do i = 1, imc-1
      upert(i,j,k) = uuc(i,j,k) - uu2d(i,j,k)
      vpert(i,j,k) = vvc(i,j,k) - vv2d(i,j,k)
      wpert(i,j,k) = wwc(i,j,k) - ww2d(i,j,k) 
      vopert(i,j,k)= vxc(i,j,k) - (vort2d(i,j,k) -fcor) 
      thpert(i,j,k)= qtc(i,j,k) - th2d(i,j,k)
     enddo
    enddo
   enddo
   
   edrvf=0.0   ! Eddy radial flux of absolute vorticity
   eddvdz=0.0  ! vertical advection of mean tangential momentum
   do k = 2, kmz-2
    do j = 1, jmc
     atemp=0.0
     btemp=0.0
     do i = 1, imc-1
      atemp = atemp + upert(i,j,k)*vopert(i,j,k)
      btemp = btemp + 0.5*wpert(i,j,k)*(vpert(i,j,k+1)-vpert(i,j,k-1))/(zlevel(k+1)-zlevel(k))
     enddo
     do i = 1, imc-1 ! kept as 3d variable for plotting purpose
       edrvf(i,j,k) = atemp/(imc-1)  ! Eddy radial flux of absolute vorticity
       eddvdz(i,j,k)= btemp/(imc-1)  ! vertical advection of mean tangential momentum
     enddo
    enddo
   enddo


!* For SEE computation
!* ~~~~~~~~~~~~~~~~~~~
   thdot=0.0
   edth1=0.0
   edth2=0.0
   edth3=0.0
   fthdiff=0.0
   fmomdiff=0.0
   do k = 2, kmz-2
    do j = 1, jmc        ! remember j is in radial
     atemp=0.0
     btemp=0.0
     ctemp=0.0
     dtemp=0.0
     etemp=0.0
     ftemp=0.0
     do i = 1, imc-1    ! remember i is along azimuth
      atemp = atemp + mc_heatc(i,j,k)
      btemp = btemp + 0.5*upert(i,j,k)*(thpert(i,j+1,k)-thpert(i,j-1,k))/yrlen               ! remember j is along radial
      ctemp = ctemp + 0.5*(vpert(i,j,k)/(j*yrlen))*(thpert(i+1,j,k)-thpert(i-1,j,k))/dlamda  ! remember i is along azimuth
      dtemp = dtemp + 0.5*wpert(i,j,k)*(thpert(i,j,k+1)-thpert(i,j,k-1))/(zlevel(k+1)-zlevel(k)) 
     enddo
     do i = 1, imc-1 ! kept as 3d variable for plotting purpose
       thdot(i,j,k)  = atemp/(imc-1)
       edth1(i,j,k)  = btemp/(imc-1)  ! Mean Eddy radial advection of potential temperature perturbation   
       edth2(i,j,k)  = ctemp/(imc-1)  ! Mean Tangential advection of potential temperature perturbation 
       edth3(i,j,k)  = dtemp/(imc-1)  ! Mean Vertical advection of potential temperature perturbation
     enddo
    enddo
   enddo

!
! Diffsion terms in SEE Momentum Equation

   do k = 2, kmz-2
    do j = 1, jmc       ! remember j is in radial
     atemp=0.0
     btemp=0.0
     do i = 1, imc-1    ! remember i is along azimuth
        atemp = atemp + vdiffc(i,j,k)
        btemp = btemp + hvdiffc(i,j,k) 
     enddo
     do i = 1, imc-1 ! kept as 3d variable for plotting purpose
        zdiff2d(i,j,k) = atemp/(imc-1)    
        hdiff2d(i,j,k) = btemp/(imc-1)  
     enddo
    enddo
   enddo

!----------------- End of mean budget computation -------------

!  determine the radius of maximum wind at each height

   do k = 1, kmz-2
    atemp = Maxval(vv2d(1,2:jmc-1,K))
    do j = 1, jmc
     if(vv2d(1,j,k) == atemp)rad1d(k)=j*yrlen/1000.
    enddo
   enddo

!  btemp=66.785 - 0.09102*MaxVal(vv2d(1,2:jmc-1,1)) + 1.0619*(storm_lat -25.0)

   btemp=u1d(15)-u1d(45)   ! shear across 1.5 km and 10 km
   ctemp=v1d(15)-v1d(45)   ! shear across 1.5 km and 10 km
   atemp=sqrt(btemp*btemp+ctemp*ctemp)
   dtemp=(180.0/pi_cnst)*atan2(btemp,ctemp)
   WRITE(55,1010)ifcst_hour,ifcst_minute,MaxVal(vv2d(1,2:jmc-1,1:20)),  &
                 MinVal(pp2d(1,2:jmc-1,1)),rad1d(1),btemp,ctemp,atemp,dtemp

1010 format(2i5,',',3(f12.2,','),3(f12.2,','),f12.2)

!
!  Compute azimuthal average for W and temperature perturbation
!
    mean_var1=0.0
    mean_var2=0.0
    do k = 1, kmz-2
     do j = 1, jmc-1
      atemp = 0.0
      btemp = 0.0
      do i = 1, imc-1
       atemp = atemp + wwc(i,j,k) 
       btemp = btemp + tpert(i,j,k) 
      enddo
      do i = 1, imc-1
       mean_var1(i,j,k) = atemp/(imc-1)     ! kept as 3d variable for plotting purpose
       mean_var2(i,j,k) = btemp/(imc-1)
      enddo
     enddo
    enddo
!
!  Compute rmse along the azimuthal for W and temperature perturbation
!
   rms_var1=0.0
   rms_var2=0.0
   do k = 1, kmz-1
    do j = 1, jmc-1 
     atemp = 0.0
     btemp = 0.0
     do i = 1, imc-1
       atemp = atemp + (wwc(i,j,k)-mean_var1(i,j,k))**2   
       btemp = btemp + (tpert(i,j,k)-mean_var2(i,j,k))**2  
     enddo 
     do i = 1, imc-1
       rms_var1(i,j,k)=sqrt(atemp/(imc-1))  ! kept as 3d variable for plotting purpose            
       rms_var2(i,j,k)=sqrt(btemp/(imc-1))
     enddo
    enddo
   enddo

! Track the number of extreme events discounting the same system (tower), i.e., column maximum
! The assumption here is that each updraft is associated with thermal perturbations


   wstat2d=0.0
   tstat2d=0.0
   do j = 1, jmc-1
    do i=1, imc-1
       wstat2d(i,j)=MaxVal( wwc(i,j,1:kmz-10))   ! 10 levels from the top is neglected
       tstat2d(i,j)=MaxVal( tpert(i,j,1:kmz-10)) ! 10 levels from the top is neglected
    enddo
   enddo

!  Find ratio of updrafts between > 5 m/s to those < 5 m/s 

   mkont=0
   kkont=0
   do j = 1, jmc-1
    do i=1, imc-1
      if(wstat2d(i,j) .ge. 5.0)THEN
         mkont=mkont+1
!       WRITE(23,*)ifcst_hour,i,j,wstat2d(i,j),tstat2d(i,j)
       else if(wstat2d(i,j) .gt. 0.5 .and. wstat2d(i,j) .lt. 5.0)then
         kkont=kkont+1
!        WRITE(24,*)ifcst_hour,i,j,wstat2d(i,j),tstat2d(i,j)
       endif
    enddo
   enddo

! Find the ratio of extreme events to the total number of events

   atemp=0.
   if(kkont .gt. 0)atemp=float(mkont)/(kkont+mkont)
   WRITE(21,104)ifcst_hour,ifcst_minute,MaxVal(mean_var1(:,:,1:50)),MaxVal(rms_var1(:,:,1:50)),MaxVal(wstat2d),mkont,kkont,atemp*100.
   WRITE(22,104)ifcst_hour,ifcst_minute,MaxVal(mean_var2(:,:,1:50)),MaxVal(rms_var2(:,:,1:50)),MaxVal(tstat2d),mkont,kkont,atemp*100.

104 FORMAT(2I5,',',3(F12.4,','),2(I10,','),F8.1)

!  maximum of extreme events in a given azimuth (variation along r-z)

   atemp=0.
   btemp=0.
   do k = 1, kmz-1
    do j = 1, jmc-1
       atemp=MaxVal( wwc(1:imc-1,j,k))
       btemp=MaxVal( tpert(1:imc-1,j,k))
       do i=1, imc-1
        wstat(i,j,k)=atemp
        tstat(i,j,k)=btemp
       enddo
    enddo
   enddo

!  Momentum Budgets Computation

   do k = 1, kmz-1
    do j = 2, jmc-1
     do i = 2, imc-1
!
      ror(i,j,k) = ppc(i,j,k)*100./( rd_air*(1+.608*qqc(i,j,k))*ttc(i,j,k))              ! rho
!
      pgr(i,j,k) = -50.*(ppc(i,j+1,k) - ppc(i,j-1,k))/yrlen                              ! radial pgf
!
      vgr(i,j,k) = (vvc(i,j,k)**2/(j*yrlen)) + fcor*vvc(i,j,k)                           ! centrepetal force + fv
!     vgr(i,j,k) = j*yrlen                                                               ! for Paul
!
      uadc(i,j,k) = 0.5*uuc(i,j,k)*(uuc(i,j+1,k)-uuc(i,j-1,k))/yrlen &                   ! radial advect
                  - (vvc(i,j,k)/(j*yrlen))*0.5*(uuc(i+1,j,k)-uuc(i-1,j,k))/dlamda
!
      vadc(i,j,k) = 0.5*uuc(i,j,k)*(vvc(i,j+1,k)-vvc(i,j-1,k))/yrlen      &
                  - (vvc(i,j,k)/(j*yrlen))*0.5*(vvc(i+1,j,k)-vvc(i-1,j,k))/dlamda

     tan_grd_wind(i,j,k) = 0.5*(1.0/(ror(i,j,k)*j*yrlen))*(ppc(i+1,j,k)-ppc(i-1,j,k))/dlamda &
                          - (uuc(i,j,k)*vvc(i,j,k)/(j*yrlen)) - fcor*uuc(i,j,k)

!     vv_verify(i,j,k)   = (j*yrlen)*(vvc(i,j,k) + fcor*j*yrlen/2.0)        ! angular momentum

      vv_verify(i,j,k)   = - (uuc(i,j,k)*vvc(i,j,k)/(j*yrlen)) - fcor*uuc(i,j,k)   ! gencor

      tadc(i,j,k) = 0.5*uuc(i,j,k)*(ttc(i,j+1,k)-ttc(i,j-1,k))/yrlen    !&           ! radial temperature advect
!                 - (vvc(i,j,k)/(j*yrlen))*0.5*(ttc(i+1,j,k)-ttc(i-1,j,k))/dlamda

     enddo
    enddo
   enddo

! Azimuthal Averages of latent, sensible, core warming and radial vorticity advection

   lfluxz=0.0
   hfluxz=0.0
   xstress=0.0
   ystress=0.0
   do j = 1, jmc-1
     atemp=0.0
     btemp=0.0
     ctemp=0.0
     dtemp=0.0
     do i = 1, imc-1
       atemp = atemp + wks(i,j,16)
       btemp = btemp + wks(i,j,15) 
       ctemp = ctemp + wks(i,j,19)
       dtemp = dtemp + wks(i,j,20)
     enddo
     lfluxz(j) = atemp/(imc-1)
     hfluxz(j) = btemp/(imc-1)
     xstress(j) = ctemp/(imc-1)
     ystress(j) = dtemp/(imc-1)
     stress(j)=sqrt(xstress(j)**2 + ystress(j)**2)
   enddo

   WRITE(56,1011)ifcst_hour,ifcst_minute,MaxVal( lfluxz(2:jmc-1) ),MaxVal( hfluxz(2:jmc-1)),MaxVal(stress(2:jmc-1))

1011 format(2I3,',',F10.2,',',F10.2)

   warmxz=0.0
   uvortxz=0.0
   do k = 1, kmz-2
    do j = 1, jmc-1
     atemp=0.0
     btemp=0.0
     do i = 1, imc-1
       atemp = atemp + tpert(i,j,k)
       btemp = btemp + tan_grd_wind(i,j,k) + vadc(i,j,k) + vdiffc(i,j,k)     !uuc(i,j,k)*vxc(i,j,k)    ! ud(vort)/dr 
     enddo
     warmxz(j,k) = atemp/(imc-1)
     uvortxz(j,k)= btemp/(imc-1)
    enddo
   enddo


!   WRITE(55,1010)ifcst_hour,MaxVal( vv2d(1,2:jmc-1,1) ),MaxVal( lfluxz(2:jmc-1) ),MaxVal( hfluxz(2:jmc-1) ),MaxVal( warmxz(2:jmc-1,1:50) ),3600.*Maxval( uvortxz(2:jmc-1,1:20) ),3600.*Maxval( uvortxz(2:jmc-1,21:50) ),rad1d(1)

   
!* Cylindrical Boundary conditions 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   wks(0,:, :) = wks(imc,:,:) ; ii= icen(0); jj= jcen(0)
   wks(:,0, 1) = xrs(ii,jj)   ; wks(:,0, 2) = yrs(ii,jj)
   wks(:,0, 3) =  zs(ii,jj)   ; wks(:,0, 4) =  ps(ii,jj)
   wks(:,0, 5) =  ts(ii,jj)   ; wks(:,0, 6) =  qs(ii,jj)
   wks(:,0, 7) =  p0(ii,jj)   ; wks(:,0, 8) = prc(ii,jj)
   wks(:,0,11) = t10(ii,jj)   ; wks(:,0,12) = q10(ii,jj)
   wks(:,0,13) = sst(ii,jj)   ; wks(:,0,14) = ssq(ii,jj)
   wks(:,0,15) =  tw(ii,jj)   ; wks(:,0,16) =  qw(ii,jj)
   wks(:,0,17) =  uf(ii,jj)   ; wks(:,0,18) =  z0(ii,jj)
   wks(:,0,19) =  tx(ii,jj)   ; wks(:,0,20) =  ty(ii,jj)  
 
   DO k = 1, kmz ; ii = icen(k) ; jj = jcen(k) 
      ppc(0,:,k) = ppc(imc,:,k) ; ppc(:,0,k) = pp_(ii,jj,k) 
      qqc(0,:,k) = qqc(imc,:,k) ; qqc(:,0,k) = qq_(ii,jj,k) 
      cwc(0,:,k) = cwc(imc,:,k) ; cwc(:,0,k) = cw_(ii,jj,k) 
      f_icec(0,:,k) = f_icec(imc,:,k) ; f_icec(:,0,k) = f_ice_(ii,jj,k)
      f_rainc(0,:,k) = f_rainc(imc,:,k) ; f_rainc(:,0,k) = f_rain_(ii,jj,k)
      f_rimefc(0,:,k) = f_rimefc(imc,:,k) ; f_rimefc(:,0,k) = f_rimef_(ii,jj,k)
      mc_heatc(0,:,k) = mc_heatc(imc,:,k) ; mc_heatc(:,0,k) = mc_heat_(ii,jj,k)
      ttc(0,:,k) = ttc(imc,:,k) ; ttc(:,0,k) = tt_(ii,jj,k) 
      tdc(0,:,k) = tdc(imc,:,k) ; tdc(:,0,k) = td_(ii,jj,k) 
      rhc(0,:,k) = rhc(imc,:,k) ; rhc(:,0,k) = rh_(ii,jj,k) 
      qtc(0,:,k) = qtc(imc,:,k) ; qtc(:,0,k) = qt_(ii,jj,k) 
      uxc(0,:,k) = uxc(imc,:,k) ; uxc(:,0,k) = ux_(ii,jj,k) 
      vxc(0,:,k) = vxc(imc,:,k) ; vxc(:,0,k) = vx_(ii,jj,k) 
      wwc(0,:,k) = wwc(imc,:,k) ; wwc(:,0,k) = ww_(ii,jj,k) 
      uuc(0,:,k) = uuc(imc,:,k) 
      vvc(0,:,k) = vvc(imc,:,k) 
      mwc(0,:,k) = mwc(imc,:,k) 
      qwc(0,:,k) = qwc(imc,:,k) 

!* Mean Budget Quantities
!* ~~~~~~~~~~~~~~~~~~~~~~

   uu2d(0,:,k)  = 0.0; uu2d(:,0,k)  = 0.0 
   vv2d(0,:,k)  = 0.0; vv2d(:,0,k)  = 0.0
   ww2d(0,:,k)  = 0.0; ww2d(:,0,k)  = 0.0
   vort2d(0,:,k)= 0.0; vort2d(:,0,k)= 0.0
   th2d(0,:,k)= 0.0; th2d(:,0,k)= 0.0
   pp2d(0,:,k)  = pp2d(imc-1,:,k); pp2d(:,0,k)  = pp_(ii,jj,k)
   tt2d(0,:,k)  = tt2d(imc-1,:,k); tt2d(:,0,k)  = tt_(ii,jj,k) 
   dvdz(0,:,k)= 0.0; dvdz(:,0,k)= 0.0
   edrvf(0,:,k)= 0.0; edrvf(:,0,k)= 0.0 
   eddvdz(0,:,k)= 0.0;eddvdz(:,0,k)= 0.0
   upert(0,:,k)= 0.0; upert(:,0,k)= 0.0
   vpert(0,:,k)= 0.0; vpert(:,0,k)= 0.0
   wpert(0,:,k)= 0.0; wpert(:,0,k)= 0.0
   vopert(0,:,k)= 0.0; vopert(:,0,k)= 0.0
   thpert(0,:,k)= 0.0; thpert(:,0,k)= 0.0

!* For SEE computation
!* ~~~~~~~~~~~~~~~~~~~ 
   zdiff2d(0,:,k)=0.0;zdiff2d(:,0,k)=0.0
   hdiff2d(:,0,k)=0.0;hdiff2d(0,:,k)=0.0
   thdot(0,:,k)=0.0;thdot(:,0,k)=0.0
   edth1(0,:,k)=0.0;edth1(:,0,k)=0.0
   edth2(0,:,k)=0.0;edth2(:,0,k)=0.0
   edth3(0,:,k)=0.0;edth3(:,0,k)= 0.0
   fthdiff(0,:,k)=0.0;fthdiff(:,0,k)=0.0
   fmomdiff(0,:,k)=0.0;fmomdiff(:,0,k)=0.0   

!* Advection quantities (gopal's doing)
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      uadc(0,:,k) = uadc(imc,:,k) ; uadc(:,0,k) = 0.0
      tadc(0,:,k) = tadc(imc,:,k) ; tadc(:,0,k) = 0.0 
      vadc(0,:,k) = vadc(imc,:,k) ; vadc(:,0,k) = 0.0 
      tan_grd_wind(0,:,k) = tan_grd_wind(imc,:,k) ; tan_grd_wind(:,0,k) = 0.0 
      tpert(0,:,k) = tpert(imc,:,k) ; tpert(:,0,k) = 0.0
      vv_verify(0,:,k) = vv_verify(imc,:,k) ; vv_verify(:,0,k) = 0.0
      pgr(0,:,k) = pgr(imc,:,k) ; pgr(:,0,k) = 0.0
      vgr(0,:,k) = vgr(imc,:,k) ; vgr(:,0,k) = 0.0
      ror(0,:,k) = ror(imc,:,k) ; ror(:,0,k) = 0.0

!      udiffc(0,:,k) = udiffc(imc,:,k) ; udiffc(:,0,k) = udiffc(ii,jj,k)
!      vdiffc(0,:,k) = vdiff(imc,:,k) ; vdiff(:,0,k) = vdiff(ii,jj,k)

      kmc(0,:,k) = kmc(imc,:,k) ; kmc(:,0,k) = 0.0
      khc(0,:,k) = khc(imc,:,k) ; khc(:,0,k) = 0.0

   ENDDO 
   
!* Output height-level data on cylindrical grid 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   do k = 1, k2d ; WRITE(310) wks(:,:,k) ; enddo 
   do k = 1, kmz ; WRITE(310) ppc(:,:,k) ; enddo 
   do k = 1, kmz ; WRITE(310) qqc(:,:,k) ; enddo 
   do k = 1, kmz ; WRITE(310) cwc(:,:,k) ; enddo 
   do k = 1, kmz ; WRITE(310) ttc(:,:,k) ; enddo 
   do k = 1, kmz ; WRITE(310) tdc(:,:,k) ; enddo 
   do k = 1, kmz ; WRITE(310) qtc(:,:,k) ; enddo 
   do k = 1, kmz ; WRITE(310) uuc(:,:,k) ; enddo 
   do k = 1, kmz ; WRITE(310) vvc(:,:,k) ; enddo 
   do k = 1, kmz ; WRITE(310) wwc(:,:,k) ; enddo 
   do k = 1, kmz ; WRITE(310) uxc(:,:,k) ; enddo 
   do k = 1, kmz ; WRITE(310) vxc(:,:,k) ; enddo 
   do k = 1, kmz ; WRITE(310) mwc(:,:,k) ; enddo 
   do k = 1, kmz ; WRITE(310) qwc(:,:,k) ; enddo
   do k = 1, kmz ; WRITE(310) rhc(:,:,k) ; enddo 

!
!* gopal's addition for momentum budget studies 
!*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   do k = 1, kmz ; WRITE(310) uu2d(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) vv2d(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) ww2d(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) pp2d(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) tt2d(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) vort2d(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) th2d(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) upert(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) vpert(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) wpert(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) vopert(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) thpert(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) dvdz(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) edrvf(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) eddvdz(:,:,k); enddo
!
!* gopal's addition for radial and tangential advection and diffusion studies
!*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   do k = 1, kmz ; WRITE(310) uadc(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) ror(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) pgr(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) vgr(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) udiffc(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) vadc(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) tan_grd_wind(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) vdiffc(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) tpert(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) vv_verify(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) tadc(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) mean_var1(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) rms_var1(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) wstat(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) mean_var2(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) rms_var2(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) tstat(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) f_icec(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) f_rainc(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) f_rimefc(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) kmc(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) khc(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) hudiffc(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) hvdiffc(:,:,k); enddo

!* For SEE computation 
!*
   do k = 1, kmz ; WRITE(310) thdot(:,:,k); enddo   
   do k = 1, kmz ; WRITE(310) edth1(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) edth2(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) edth3(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) zdiff2d(:,:,k); enddo
   do k = 1, kmz ; WRITE(310) hdiff2d(:,:,k); enddo

   nvars = k2d + 59           ! k2d=20
   
   DeAllocate( icen, jcen, rcos, rsin, rhc, wks )
   DeAllocate( ppc, qqc, cwc, ttc, tdc, qtc, uuc, vvc, wwc, uxc, vxc, mwc, qwc )
   DeAllocate( f_icec, f_rainc, f_rimefc )
   DeAllocate(                  mc_heatc )

   DeAllocate( uadc, vadc, tadc)  ! gopal's doing
   DeAllocate( udiffc, vdiffc)
   DeAllocate(hudiffc,hvdiffc)
   DeAllocate( tan_grd_wind)
   DeAllocate( tpert, vv_verify)
   DeAllocate( u1d,v1d,t1d,rad1d,pgr,vgr,ror)
   DeAllocate( upert,vpert,wpert,vopert,thpert)
   DeAllocate( uu2d,vv2d,ww2d,pp2d,tt2d,vort2d,th2d)
   DeAllocate( dvdz,edrvf,eddvdz)
   DeAllocate( mean_var1,rms_var1,wstat2d,wstat)
   DeAllocate( mean_var2,rms_var2,tstat2d,tstat)
   DeAllocate( lfluxz,hfluxz,warmxz,uvortxz)
   DeAllocate(kmc,khc)

!* For SEE computation
!* ~~~~~~~~~~~~~~~~~~~
   DeAllocate(zdiff2d,hdiff2d)
   DeAllocate(thdot,edth1,edth2,edth3,fthdiff,fmomdiff)
   
!           *     *     *     *     *     *     *     *     *     * 
   
!* GrADS control file for height-cylinder data on rectangular grid 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   dxx = 360./float(imc) !* resolution in azimuthal dir (deg) 
   dyy =  yrlen/1000.    !* resolution in  radial   dir (km) 
   
   WRITE(311,'("                                                      ")') 
   WRITE(311,'(" DSET    height_cyl.dat                               ")') 
   WRITE(311,'(" TITLE   Height-cylinder data on rectangular grid     ")') 
   WRITE(311,'(" OPTIONS sequential big_endian                        ")') 
   WRITE(311,'(" UNDEF",E11.2)') undef 
   WRITE(311,'("                                                      ")') 
   WRITE(311,'(" XDEF", I7,"   LINEAR", 2F10.3)') imc+1, -180., dxx 
   WRITE(311,'(" YDEF", I7,"   LINEAR", 2F10.3)') jmc+1,   0. , dyy 
   WRITE(311,'(" ZDEF", I7,"   LEVELS")') kmz 
   WRITE(311,*)  zlevel 
   WRITE(311,'(" TDEF", I7,"   LINEAR   00:00Z01JAN0001  1yr")') ntime 
   WRITE(311,'("                                                      ")') 
   WRITE(311,'(" VARS", I3)') nvars
   WRITE(311,'(" xx        0   0    Earth longitude (deg)             ")') 
   WRITE(311,'(" yy        0   0    Earth latitude  (deg)             ")') 
   WRITE(311,'(" zs        0   0    Earth topography (m)              ")') 
   WRITE(311,'(" ps        0   0    surface air   pressure (mb)       ")') 
   WRITE(311,'(" ts        0   0    surface air temperature (K)       ")') 
   WRITE(311,'(" qs        0   0    surface specific humidity (kg/kg) ")') 
   WRITE(311,'(" p0        0   0    mean-sea-level pressure (mb)      ")') 
   WRITE(311,'(" prc       0   0    precipitation (m)                 ")') 
   WRITE(311,'(" u10       0   0    10m   radial   wind (m/s)         ")') 
   WRITE(311,'(" v10       0   0    10m tangential wind (m/s)         ")') 
   WRITE(311,'(" t10       0   0    10m air temperature (K)           ")') 
   WRITE(311,'(" q10       0   0    10m specific humidity (kg/kg)     ")') 
   WRITE(311,'(" sst       0   0    sea-surface temperature (K)       ")') 
   WRITE(311,'(" ssq       0   0    sea-sfc specific humidity (kg/kg) ")') 
   WRITE(311,'(" tw        0   0    sensible heat flux (W/m^2)        ")') 
   WRITE(311,'(" qw        0   0    latent   heat flux (W/m^2)        ")') 
   WRITE(311,'(" uf        0   0    friction velocity  (m/s)          ")') 
   WRITE(311,'(" z0        0   0    roughness height   (m)            ")') 
   WRITE(311,'(" tx        0   0    radial stress (kg/m/s^2)          ")')
   WRITE(311,'(" ty        0   0    tangential stress (kg/m/s^2)      ")')
   WRITE(311,'(" pp  ", I7,"   0    air  pressure  (mb)               ")') kmz 
   WRITE(311,'(" qq  ", I7,"   0    specific humidity (kg/kg)         ")') kmz 
   WRITE(311,'(" cw  ", I7,"   0    total condensate  (kg/kg)         ")') kmz 
   WRITE(311,'(" tt  ", I7,"   0    air  temperature  (K)             ")') kmz 
   WRITE(311,'(" td  ", I7,"   0    dewpoint temperature (K)          ")') kmz 
   WRITE(311,'(" qt  ", I7,"   0    equiv. pot. temperature (K)       ")') kmz 
   WRITE(311,'(" uu  ", I7,"   0    outward radial wind (m/s)         ")') kmz 
   WRITE(311,'(" vv  ", I7,"   0    CCW tangential wind (m/s)         ")') kmz 
   WRITE(311,'(" ww  ", I7,"   0    vertical  wind (m/s)              ")') kmz 
   WRITE(311,'(" ux  ", I7,"   0    horizontal divergence (/s)        ")') kmz 
   WRITE(311,'(" vx  ", I7,"   0    relative vorticity (/s)           ")') kmz 
   WRITE(311,'(" mwc ", I7,"   0    vertical mass flux (kg/m^2/s)     ")') kmz 
   WRITE(311,'(" qwc ", I7,"   0    vertical moisture flux (g/m^2/s)  ")') kmz
   WRITE(311,'(" rh  ", I7,"   0    relative humidity (%)             ")') kmz 

!
!  gopal's addition for mean momentum budget
!
   WRITE(311,'("ubar ", I7,"   0    mean radial wind  (m/s^2)         ")') kmz
   WRITE(311,'("vbar ", I7,"   0    mean tangential wind  (m/s^2)     ")') kmz 
   WRITE(311,'("wbar ", I7,"   0    mean vertical wind  (m/s^2)         ")') kmz
   WRITE(311,'("pbar ", I7,"   0    mean pressure wind  (mb)         ")') kmz
   WRITE(311,'("tbar ", I7,"   0    mean temperature  (k)         ")') kmz
   WRITE(311,'("vxbar ", I7,"   0   mean vorticity  (m/s^2)         ")') kmz
   WRITE(311,'("thbar ", I7,"   0   mean potential temperature  (K)         ")') kmz
   WRITE(311,'("upert ", I7,"   0   upert  (m/s^2)         ")') kmz
   WRITE(311,'("vpert ", I7,"   0   upert  (m/s^2)         ")') kmz
   WRITE(311,'("wpert ", I7,"   0   upert  (m/s^2)         ")') kmz
   WRITE(311,'("vopert ", I7,"   0  vopert (/s^1)         ")') kmz
   WRITE(311,'("thpert ", I7,"   0  thpert (K)         ")') kmz
   WRITE(311,'("dvdz  ", I7,"   0   d<v>/dz  (1/s)                  ")') kmz
   WRITE(311,'("edrvf  ", I7,"   0  eddy radial vorticity flux   (m/s^2) ")') kmz
   WRITE(311,'("eddvdz  ", I7,"  0  eddy vertical advection of v (m/s^2) ")') kmz
!
!  gopal's addition for computing imbalances
!
   WRITE(311,'("radv ", I7,"   0    radial advection  (m/s^2)         ")') kmz
   WRITE(311,'("ror  ", I7,"   0    moist density (kg/m^3)            ")') kmz
   WRITE(311,'("pgr  ", I7,"   0    radial PGF (pa/m)                 ")') kmz
   WRITE(311,'("vgr  ", I7,"   0    centrepetal force + fv (m/s^2)    ")') kmz
   WRITE(311,'("rdif ", I7,"   0    radial diffusion wind  (m/s^2)    ")') kmz
   WRITE(311,'("tadv ", I7,"   0    tangential advection (m/s^2)      ")') kmz
   WRITE(311,'("tgrd ", I7,"   0    tangential gradient wind (m/s^2)      ")') kmz
   WRITE(311,'("tdif ", I7,"   0    tangential diffusion (m/s^2)      ")') kmz
   WRITE(311,'("tprt ", I7,"   0    temperature perturbation (k)      ")') kmz
   WRITE(311,'("gencor", I7,"   0    generalized Coriolis (m/s^2)          ")') kmz
   WRITE(311,'("tadv ", I7,"   0    radial temp advect (k/s)          ")') kmz
   WRITE(311,'("mean_var1",I7,"  0    azimuthan mean of var1 (m/s)        ")') kmz
   WRITE(311,'("rms_var1", I7,"  0    rms of var1 (m/s)                   ")') kmz
   WRITE(311,'("wstat", I7,"  0    stats of w (m/s)                   ")') kmz
   WRITE(311,'("mean_var2",I7,"  0    azimuthan mean of tpert (k)        ")') kmz
   WRITE(311,'("rms_var2", I7,"  0    rms of tpert (k)                   ")') kmz
   WRITE(311,'("tstat", I7,"  0    stats of tpert (k)                   ")') kmz
   WRITE(311,'("f_ice", I7,"  0    ice fraction           ")') kmz
   WRITE(311,'("f_rain", I7,"  0   rain fraction           ")') kmz
   WRITE(311,'("f_rimef", I7,"  0  rime factor           ")') kmz 
   WRITE(311,'("kmc", I7,"  0    Eddy Diffusivity for Momentum (m^2/s)   ")') kmz
   WRITE(311,'("khc", I7,"  0    Eddy Diffusivity for Heat (m^2/s)   ")') kmz
   WRITE(311,'("hrdif ", I7,"   0    lateral radial diffusion wind  (m/s^2)    ")') kmz
   WRITE(311,'("htdif ", I7,"   0    lateral tangential diffusion (m/s^2)      ")') kmz
!
!* For SEE computation
!*
   WRITE(311,'("thdot", I7,"   0 Mean heating (K per time step) ")')kmz
   WRITE(311,'("edth1", I7,"   0 Term 1- eddy term in SEE (K/s) ")')kmz 
   WRITE(311,'("edth2", I7,"   0 Term 2- eddy term in SEE (K/s) ")')kmz
   WRITE(311,'("edth3", I7,"   0 Term 3- eddy term in SEE (K/s) ")')kmz
   WRITE(311,'("zdiff", I7,"   0 Mean Diffusion along Z   (m^2/s)")')kmz 
   WRITE(311,'("vdiff", I7,"   0 Mean Diffusion along XY  (m/s^2)")') kmz

   WRITE(311,'(" ENDVARS                                              ")') 
   WRITE(311,'("                                                      ")') 
   
!* GrADS control file for height-cylinder data on polar projection 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   dxx = 360./float(imc) !* resolution in azimuthal dir (deg) 
   dyy =  90./float(jmc) !* resolution in  radial   dir (deg) 
   
   WRITE(312,'("                                                      ")') 
   WRITE(312,'(" DSET    height_cyl.dat                               ")') 
   WRITE(312,'(" TITLE   Height-cylinder data on polar projection     ")') 
   WRITE(312,'(" OPTIONS sequential big_endian                        ")') 
   WRITE(312,'(" UNDEF",E11.2)') undef 
   WRITE(312,'("                                                      ")') 
   WRITE(312,'(" XDEF", I7,"   LINEAR", 2F10.3)') imc+1, -180., dxx 
   WRITE(312,'(" YDEF", I7,"   LINEAR", 2F10.3)') jmc+1,  -90., dyy 
   WRITE(312,'(" ZDEF", I7,"   LEVELS")') kmz 
   WRITE(312,*)  zlevel 
   WRITE(312,'(" TDEF", I7,"   LINEAR   00:00Z01JAN0001  1yr")') ntime 
   WRITE(312,'("                                                      ")') 
   WRITE(312,'(" VARS", I3)') nvars
   WRITE(312,'(" xx        0   0    Earth longitude (deg)             ")') 
   WRITE(312,'(" yy        0   0    Earth latitude  (deg)             ")') 
   WRITE(312,'(" zs        0   0    Earth topography (m)              ")') 
   WRITE(312,'(" ps        0   0    surface air   pressure (mb)       ")') 
   WRITE(312,'(" ts        0   0    surface air temperature (K)       ")') 
   WRITE(312,'(" qs        0   0    surface specific humidity (kg/kg) ")') 
   WRITE(312,'(" p0        0   0    mean-sea-level pressure (mb)      ")') 
   WRITE(312,'(" prc       0   0    precipitation (m)                 ")') 
   WRITE(312,'(" u10       0   0    10m   radial   wind (m/s)         ")') 
   WRITE(312,'(" v10       0   0    10m tangential wind (m/s)         ")') 
   WRITE(312,'(" t10       0   0    10m air temperature (K)           ")') 
   WRITE(312,'(" q10       0   0    10m specific humidity (kg/kg)     ")') 
   WRITE(312,'(" sst       0   0    sea-surface temperature (K)       ")') 
   WRITE(312,'(" ssq       0   0    sea-sfc specific humidity (kg/kg) ")') 
   WRITE(312,'(" tw        0   0    sensible heat flux (W/m^2)        ")') 
   WRITE(312,'(" qw        0   0    latent   heat flux (W/m^2)        ")') 
   WRITE(312,'(" uf        0   0    friction velocity  (m/s)          ")') 
   WRITE(312,'(" z0        0   0    roughness height   (m)            ")') 
   WRITE(312,'(" tx        0   0    radial stress (kg/m/s^2)          ")')
   WRITE(312,'(" ty        0   0    tangential stress (kg/m/s^2)      ")')
   WRITE(312,'(" pp  ", I7,"   0    air  pressure  (mb)               ")') kmz 
   WRITE(312,'(" qq  ", I7,"   0    specific humidity (kg/kg)         ")') kmz 
   WRITE(312,'(" cw  ", I7,"   0    total condensate  (kg/kg)         ")') kmz 
   WRITE(312,'(" tt  ", I7,"   0    air  temperature  (K)             ")') kmz 
   WRITE(312,'(" td  ", I7,"   0    dewpoint temperature (K)          ")') kmz 
   WRITE(312,'(" qt  ", I7,"   0    equiv. pot. temperature (K)       ")') kmz 
   WRITE(312,'(" uu  ", I7,"   0    outward radial wind (m/s)         ")') kmz 
   WRITE(312,'(" vv  ", I7,"   0    CCW tangential wind (m/s)         ")') kmz 
   WRITE(312,'(" ww  ", I7,"   0    vertical  wind (m/s)              ")') kmz 
   WRITE(312,'(" ux  ", I7,"   0    horizontal divergence (/s)        ")') kmz 
   WRITE(312,'(" vx  ", I7,"   0    relative vorticity (/s)           ")') kmz 
   WRITE(312,'(" mwc ", I7,"   0    vertical mass flux (kg/m^2/s)     ")') kmz 
   WRITE(312,'(" qwc ", I7,"   0    vertical moisture flux (g/m^2/s)  ")') kmz
   WRITE(312,'(" rh  ", I7,"   0    relative humidity (%)             ")') kmz 

!
!  gopal's addition for mean momentum budget
!
   WRITE(312,'("ubar ", I7,"   0    mean radial wind  (m/s^2)         ")') kmz
   WRITE(312,'("vbar ", I7,"   0    mean tangential wind  (m/s^2)         ")') kmz
   WRITE(312,'("wbar ", I7,"   0    mean vertical wind  (m/s^2)         ")') kmz
   WRITE(312,'("pbar ", I7,"   0    mean pressure wind  (mb)         ")') kmz
   WRITE(312,'("tbar ", I7,"   0    mean temperature  (k)         ")') kmz
   WRITE(312,'("vxbar ", I7,"   0   mean vorticity  (m/s^2)         ")') kmz
   WRITE(312,'("thbar ", I7,"   0   mean potential temperature  (K)         ")')kmz
   WRITE(312,'("upert ", I7,"   0   upert  (m/s^2)         ")') kmz
   WRITE(312,'("vpert ", I7,"   0   upert  (m/s^2)         ")') kmz
   WRITE(312,'("wpert ", I7,"   0   upert  (m/s^2)         ")') kmz
   WRITE(312,'("vopert ", I7,"   0   upert  (m/s^2)         ")') kmz
   WRITE(312,'("thpert ", I7,"   0  thpert (K)         ")') kmz
   WRITE(312,'("dvdz  ", I7,"   0   d<v>/dz  (1/s)                  ")') kmz
   WRITE(312,'("edrvf  ", I7,"   0  eddy radial vorticity flux   (m/s^2) ")') kmz
   WRITE(312,'("eddvdz  ", I7,"  0  eddy vertical advection of v (m/s^2) ")') kmz
!
!
!  gopal's addition for computing imbalances
!
   WRITE(312,'("radv ", I7,"   0    radial advection  (m/s^2)         ")') kmz
   WRITE(312,'("ror  ", I7,"   0    moist density (kg/m^3)            ")') kmz
   WRITE(312,'("pgr  ", I7,"   0    radial PGF (pa/m)                 ")') kmz
   WRITE(312,'("vgr  ", I7,"   0    centrepetal force + fv (m/s^2)    ")') kmz
   WRITE(312,'("rdif ", I7,"   0    radial diffusion wind  (m/s^2)    ")') kmz
   WRITE(312,'("tadv ", I7,"   0    tangential advection (m/s^2)      ")') kmz
   WRITE(312,'("tgrd ", I7,"   0    tangential gradient wind (m/s^2)      ")') kmz
   WRITE(312,'("tdif ", I7,"   0    tangential diffusion (m/s^2)      ")') kmz
   WRITE(312,'("tprt ", I7,"   0    temperature perturbation (k)      ")') kmz
   WRITE(312,'("gencor", I7,"   0    generalized Coriolis (m/s^2)          ")') kmz
   WRITE(312,'("tadv ", I7,"   0    radial temp advect (k/s)          ")') kmz
   WRITE(312,'("mean_var1",I7,"  0    azimuthan mean of var1 (m/s)        ")') kmz
   WRITE(312,'("rms_var1", I7,"  0    rms of var1 (m/s)                   ")') kmz
   WRITE(312,'("wstat", I7,"  0    stats of w (m/s)                   ")') kmz
   WRITE(312,'("mean_var2",I7,"  0    azimuthan mean of tpert (k)        ")') kmz
   WRITE(312,'("rms_var2", I7,"  0    rms of tpert (k)                   ")') kmz
   WRITE(312,'("tstat", I7,"  0    stats of tpert (k)                   ")') kmz
   WRITE(312,'("f_ice", I7,"  0    ice fraction           ")') kmz
   WRITE(312,'("f_rain", I7,"  0   rain fraction           ")') kmz
   WRITE(312,'("f_rimef", I7,"  0  rime factor           ")') kmz
   WRITE(312,'("kmc", I7,"  0    Eddy Diffusivity for Momentum (m^2/s)   ")') kmz
   WRITE(312,'("khc", I7,"  0    Eddy Diffusivity for Heat (m^2/s)   ")') kmz
   WRITE(312,'("hrdif ", I7,"   0    lateral radial diffusion wind  (m/s^2)")') kmz
   WRITE(312,'("htdif ", I7,"   0    lateral tangential diffusion (m/s^2)")') kmz
!
!* For SEE computation
!*
   WRITE(312,'("thdot", I7,"   0 Mean heating (K per time step) ")')kmz
   WRITE(312,'("edth1", I7,"   0 Term 1- eddy term in SEE (K/s) ")')kmz
   WRITE(312,'("edth2", I7,"   0 Term 2- eddy term in SEE (K/s) ")')kmz
   WRITE(312,'("edth3", I7,"   0 Term 3- eddy term in SEE (K/s) ")')kmz
   WRITE(312,'("zdiff", I7,"   0 Mean Diffusion along Z   (m^2/s)")')kmz
   WRITE(312,'("vdiff", I7,"   0 Mean Diffusion along XY  (m/s^2)")')kmz

   WRITE(312,'(" ENDVARS                                              ")') 
   WRITE(312,'("                                                      ")') 
   
   CLOSE(310) ; CLOSE(311) ; CLOSE(312) 
   
   ENDIF !* ( height_cylin ) ------------------------------------------- 


!******************************************************! 
!  3c. Constant Frequency by Altitude Diagrams (CFADs) ! 
!______________________________________________________! 

   IF ( diagnostics ) THEN  !* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!* Define W frequency (distribution) 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   vrange = 20. ; binsize = .5 
    mbins = nint(vrange/binsize) 
   Allocate( wwfr(-mbins:mbins,1:km) ) 

   wwfr(:,:) = 0. 
   percent = 100./(im*jm) 
   DO k = 1, km 
      DO j = 1, jm ; DO i = 1, im 
         n = nint( ww(i,j,k)/binsize ) 
         m = max( -mbins, min(n,mbins) ) 
         wwfr(m,k) = wwfr(m,k) + 1. 
      ENDDO ; ENDDO 
      DO m = -mbins, mbins 
         wwfr(m,k) = wwfr(m,k)*percent 
      ENDDO 
   ENDDO 
   
   nvars = 1 
   WRITE(332) wwfr 
   DeAllocate(wwfr) 
   
!* GrADS control file for CFADs 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   imax = 2*mbins + 1 ; jmax = km 
   xxr1 = -vrange ; dxx = binsize 
   
   WRITE(331,'("                                                      ")') 
   WRITE(331,'(" DSET    cfads_demo.dat                               ")') 
   WRITE(331,'(" TITLE   Constant Frequency by Altitude Diagrams      ")') 
   WRITE(331,'(" OPTIONS sequential big_endian                        ")') 
   WRITE(331,'(" UNDEF",E11.2)') undef 
   WRITE(331,'("                                                      ")') 
   WRITE(331,'(" XDEF", I7,"   LINEAR", 2F11.5)') imax,  xxr1, dxx 
   WRITE(331,'(" YDEF", I7,"   LINEAR          1          1 ")') jmax 
   WRITE(331,'(" ZDEF", I7,"   LINEAR          1          1 ")')   1 
   WRITE(331,'(" TDEF", I7,"   LINEAR   00:00Z01JAN0001  1yr")') ntime 
   WRITE(331,'("                                                      ")') 
   WRITE(331,'(" VARS", I3)') nvars
   WRITE(331,'(" wwfr      0   0    vertical velocity distribution    ")') 
   WRITE(331,'(" ENDVARS                                              ")') 
   WRITE(331,'("                                                      ")') 
   
   CLOSE(331) ; CLOSE(332) 
   
   ENDIF  !* ( diagnostics ) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   
   IF ( height_vars ) THEN 
   DeAllocate( uu_, vv_, ww_, tt_, td_, qt_, pp_, qq_, cw_, ux_, vx_, rh_) 
   DeAllocate( difu_,difv_,km_diffus_,kh_diffus_)
   DeAllocate( hdifu_,hdifv_)
   DeAllocate( f_ice_,f_rain_,f_rimef_)
   DeAllocate(                mc_heat_)
   ENDIF 
   
   DeAllocate( xrr, yrr, xrs, yrs, wcos, wsin, urs, vrs, xrlen ) 
   CALL ext_ncd_get_next_time( dh0, DateStr, Status_next_time ) 

   ENDDO !<======== next_time ============================================== 

   DeAllocate( zs, sm, p0, prc, ps, ths, ts, qs, u10, v10, p0B)
   DeAllocate( t02, q02, t10, q10, sst, ssq, tw, qw, uf, z0, z0t ) 
   DeAllocate( uu, vv, ww, tt, td, qt, pp, qq, cw, gh, rh) 
   DeAllocate( f_ice,f_rain,f_rimef)
   DeAllocate(              mc_heat)
   CALL ext_ncd_ioclose(dh0, Status)
   
   cpu_time_b = ETIME(btime) ;  n  = itermax 
   aaa = btime(1) - atime(1) ; bbb = cpu_time_b - cpu_time_a 
   WRITE(*,'(" Elapsed  CPU  Time (sec): usr, tot =", 2F12.5)') aaa, bbb 
   if ( iter_cen > itermax ) &  
   WRITE(*,'(" Storm Center has NOT converged after ",I2," iterations!")') n
   if ( iter_eqt > itermax ) &  
   WRITE(*,'(" Eqv.Pot.Temp has NOT converged after ",I2," iterations!")') n
   if ( bad_cylin ) &  
   WRITE(*,'(" Cylinder falls OUT domain at level =", I9)') kout 
   WRITE(*,'(" Normal Termination  -- ",/)') 
   
END PROGRAM diapost !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!==============================================================================

!******************************************************************************
! 
!  Buck's dewpoint temperature & relative humidity 
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!  Calculate dewpoint temperature & relative humidity from pressure and
!  specific humidity with Buck's equation for saturation vapor pressure of 
!  water (Buck 1981). The approximated dewpoint is contrained by air temp. 
!  
!______________________________________________________________________________
   
   SUBROUTINE Dewpoint_RH( pp, qq, tt, td, rh, im, jm, km ) 
   
   Implicit None 
!  -------------------------------------------------------------------- [input] 
   Integer  im, jm, km       !* x, y, z dimension 
   Real(4)  pp(im,jm,km)     !* air pressure (mb) 
   Real(4)  qq(im,jm,km)     !* specific humidity (kg/kg) 
   Real(4)  tt(im,jm,km)     !* air temperature (K) 
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - [output] 
   Real(4)  td(im,jm,km)     !* dewpoint temperature (K) 
   Real(4)  rh(im,jm,km)     !* relative humidity (%) 
!  -------------------------------------------------------------------- [local] 
   Real(4)  t0c, e0c, eps, aaa, bbb, eee, ees, rrr, xxx, yyy, ttt, ttd 
   Integer  i, j, k 
   
!           *     *     *     *     *     *     *     *     *     * 
   
!* Dewpoint temperature (K) on sigma-p levels 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   t0c = 273.15  !* Kelvin temperature (K) at 0'C 
   e0c = 6.1121  !* water vapor pressure (mb) at 0'C 
   aaa = 17.502  !* factor for water vapor pressure 
   bbb = 240.97  !* constant for water vapor pressure 
   eps = .62198  !* ratio of vapor/air molecular weights
   
   DO k = 1, km !* ----------------------------------------------------> 
   
   DO j = 1, jm ; DO i = 1, im 
   
      rrr = qq(i,j,k)/( 1. - qq(i,j,k) ) !* mixing ratio 
   
      eee = rrr*pp(i,j,k)/( eps + rrr  ) ; xxx = log(eee/e0c) 
   
      ttd = bbb*xxx/( aaa - xxx ) + t0c  !* Buck's dewpoint (K) 
   
      td(i,j,k) = min( ttd, tt(i,j,k) )  !* constrained dewpoint 
   
      ttt = tt(i,j,k) - t0c ; ees = e0c*exp( aaa*ttt/(ttt+bbb) ) 
   
      rh(i,j,k) = min(100.*eee/ees,100.) !* relative humidity 
   
   ENDDO ; ENDDO 

   ENDDO !* ( k = 1, km ) <---------------------------------------------
   
   END SUBROUTINE Dewpoint_RH 
   
!==============================================================================

!******************************************************************************
! 
!  Convert between standard & rotated latlon coordinates 
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!  This subroutine converts latlon coordinates between two coordinate systems, 
!  the standard geographic system {S} & the rotated spherical coord system {R} 
!  obtained by rotating the origin of {S} to (olat,olon) as the new origin: 
!  (0,0) -> (olat,olon) = (0,0)'. The poles are not allowed for conversion. 
!  
!  mode = 1 : convert standard coord (alat,alon) to rotated  coord (blat,blon) 
!  mode = 2 : convert rotated  coord (blat,blon) to standard coord (alat,alon) 
!  
!______________________________________________________________________________
   
   SUBROUTINE Convert_LatLon(mode, im, jm, slon, slat, rlon, rlat, olon, olat)
   
   Implicit None 
!  -------------------------------------------------------------------- [input] 
   Integer  mode, im, jm      !* mode= 1/2: Std -> Rot / Rot -> Std coord 
   Real(4)  olon, olat        !* origin of rotated grid in standard coord 
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - [in/out] 
   Real(4)  slon(im,jm), slat(im,jm) !* standard latlon coordinates (deg) 
   Real(4)  rlon(im,jm), rlat(im,jm) !* rotated  latlon coordinates (deg) 
!  -------------------------------------------------------------------- [local] 
   Real(8)  alon(im,jm), alat(im,jm) !* standard latlon coordinates (deg) 
   Real(8)  blon(im,jm), blat(im,jm) !* rotated  latlon coordinates (deg) 
   Real(8)  pi_cnst, deg2rad, eps1, sqrt2_, xcos, xsin, ycos, ysin 
   Real(8)  olamda, otheta, alamda, atheta, blamda, btheta 
   Real(8)  aaa, bbb, ccc, ddd, eee, xxx, yyy, zzz 
   Integer   i, j 
   
!           *     *     *     *     *     *     *     *     *     * 
   
   pi_cnst = ACOS(-1.) ; deg2rad = ACOS(-1.)/180. 
     eps1  = 1. - 1.D-6 ;  sqrt2_ = 1./sqrt(2.) 
   
!* Enhance precision for trigonometric functions 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   alon(:,:) = slon(:,:) ; alat(:,:) = slat(:,:) 
   blon(:,:) = rlon(:,:) ; blat(:,:) = rlat(:,:) 
   olamda = olon*deg2rad ; otheta = olat*deg2rad 
      aaa = Dcos(otheta) ; bbb = Dsin(otheta) 
   
!* Convert standard Earth latlon to rotated Model latlon 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   IF ( mode == 1 ) THEN  !* - - - - - - - - - - - - - - - - - - - - 
   
   DO j = 1, jm 
   DO i = 1, im 
      yyy = alat(i,j)*deg2rad 
      xxx = alon(i,j)*deg2rad - olamda 
     xcos = Dcos(xxx) ; xsin = Dsin(xxx) 
     ycos = Dcos(yyy) ; ysin = Dsin(yyy) 
!    - - - - - - - - - - - - - - - - - - 
      zzz = aaa*ysin - bbb*ycos*xcos 
      zzz = max(-eps1, min(zzz, eps1)) !* Poles undefined 
     btheta = DAsin(zzz) ; eee = 1./Dcos(btheta) 
      ccc = ( bbb*ysin + aaa*ycos*xcos )*eee 
      ddd = (                ycos*xsin )*eee 
      IF (  ccc > sqrt2_ ) THEN 
         blamda = DAsin(ddd) 
      ELSEIF (  ccc  < -sqrt2_ ) THEN 
         blamda = pi_cnst*Sign(1.D+0,ddd) - DAsin(ddd) 
      ELSE ! ( |ccc| <= sqrt2_ ) ! 
         blamda = ACOS(ccc)*Sign(1.D+0,ddd) 
      ENDIF 
      blat(i,j) = btheta/deg2rad 
      blon(i,j) = blamda/deg2rad 
   ENDDO 
   ENDDO 
   
!* Convert rotated Model latlon to standard Earth latlon
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   ELSE !* ( mode == 2 ) - - - - - - - - - - - - - - - - - - - - - - 
   
   DO j = 1, jm 
   DO i = 1, im 
      yyy = blat(i,j)*deg2rad 
      xxx = blon(i,j)*deg2rad 
     xcos = Dcos(xxx) ; xsin = Dsin(xxx) 
     ycos = Dcos(yyy) ; ysin = Dsin(yyy) 
!    - - - - - - - - - - - - - - - - - - 
      zzz = aaa*ysin + bbb*ycos*xcos 
      zzz = max(-eps1, min(zzz, eps1)) !* Poles undefined 
     atheta = DAsin(zzz) ; eee = 1./Dcos(atheta) 
      ccc = (-bbb*ysin + aaa*ycos*xcos )*eee 
      ddd = (                ycos*xsin )*eee 
      IF (  ccc > sqrt2_ ) THEN 
         alamda = DAsin(ddd) + olamda 
      ELSEIF (  ccc  < -sqrt2_ ) THEN 
         alamda = pi_cnst*Sign(1.D+0,ddd) - DAsin(ddd) + olamda 
      ELSE ! ( |ccc| <= sqrt2_ ) ! 
         alamda = ACOS(ccc)*Sign(1.D+0,ddd) + olamda 
      ENDIF 
      alat(i,j) = atheta/deg2rad 
      alon(i,j) = alamda/deg2rad 
   ENDDO 
   ENDDO 
   
   ENDIF  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   
!* Recover single precision for output 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   slon(:,:) = alon(:,:) ; slat(:,:) = alat(:,:) 
   rlon(:,:) = blon(:,:) ; rlat(:,:) = blat(:,:) 
   
   END SUBROUTINE Convert_LatLon 
   
!==============================================================================
   
!******************************************************************************
! 
!  Monotonic bi-cubic interpolation to the mid interval 
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  
!  -----(1)----------(2)---x------(3)----------(4)----- 
!  
!______________________________________________________________________________
   
   SUBROUTINE Cubic_Interpo( hhh, ff, wtx, wty, undef )
   
   Implicit None 
!  -------------------------------------------------------------------- [input] 
   Real(4):: ff(4,4), undef  !* given values on grid points 
   Real(4):: wtx(4), wty(4)  !* cubic interpolation weights 
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - [output] 
   Real(4):: hhh             !* interpolated value 
!  -------------------------------------------------------------------- [local] 
   Real(4):: gg(4), aaa, bbb, xxx, yyy 
   Integer:: i, j 
   
!           *     *     *     *     *     *     *     *     *     * 
   
   DO i = 1, 4 
      IF ( MaxVal(ff(i,:)) == undef ) THEN 
         gg(i) = undef 
      ELSE
         aaa = min( ff(i,2), ff(i,3) ) ;  bbb  = max( ff(i,2), ff(i,3) ) 
         yyy = SUM(wty(1:4)*ff(i,1:4)) ; gg(i) = max( aaa, min(yyy, bbb) ) 
      ENDIF 
   ENDDO 
   
   IF ( MaxVal(gg(:)) == undef ) THEN 
      hhh = undef 
   ELSE
      aaa = min( gg(2), gg(3) )   ; bbb = max( gg(2), gg(3) ) 
      xxx = SUM(wtx(1:4)*gg(1:4)) ; hhh = max( aaa, min(xxx, bbb) ) 
   ENDIF 
   
   END SUBROUTINE Cubic_Interpo 
   
!==============================================================================
   
!******************************************************************************
! 
!  Monotonic cubic grid refinement with triple resolution 
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  
!  ----(1)-----------(2)---|---|---(3)-----------(4)----- 
!  
!______________________________________________________________________________
   
   SUBROUTINE Refine_Grid3x( pp, pp3, undef, im3, jm3, im, jm, km )
   
   Implicit  None 
!  -------------------------------------------------------------------- [input] 
   Integer:: iorder = 4    !* 2/4 
   Integer:: im3, jm3, im, jm, km 
   Real(4):: pp(im,jm,km), pp3(im3,jm3,km), undef 
!  -------------------------------------------------------------------- [local] 
   Real(4):: aaa, bbb, xxx, yyy, a1_3, b2_3, awt(4), bwt(4), ff(4) 
   Integer:: i, j, k, m, n, j1, j2, j3, j4, n0, n1, n2 
   Integer:: m0, m1, m2, mm1, mm2, mm3, mm4 
   
!           *     *     *     *     *     *     *     *     *     * 
   
   SELECT CASE (iorder)  !================================================ 

!> ====================================
   CASE (2)  !* linear interpolation *! 
!> ====================================

!> ====================================
   CASE (4)  !* cubic  interpolation *! 
!> ====================================

    a1_3  =  1./3.  ;  b2_3  =  2./3. 
   awt(1) = -5./81. ; awt(2) = 20./27. ; awt(3) = 10./27. ; awt(4) = -4./81.  
   bwt(1) = -4./81. ; bwt(2) = 10./27. ; bwt(3) = 20./27. ; bwt(4) = -5./81.  
   
   DO k = 1, km !* ----------------------------------------------------> 
   
!* Refine the grid in Y-direction 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   DO i = 1, im ; m0 = 3*(i-1) + 1 
   
      j = jm ; n0 = 3*(j-1) + 1 
      pp3(m0,n0,k) = pp(i,j,k) 
   
      DO j = 1, jm-1, jm-2 ; j2 = j ; j3 = j+1 
         n0 = 3*(j-1) + 1  ; n1 = n0+1 ; n2 = n0+2 
         ff(2:3) = pp(i,j2:j3,k) ; pp3(m0,n0,k) = ff(2) 
         IF ( max(ff(2), ff(3)) == undef ) THEN 
            pp3(m0,n1,k) = undef ; pp3(m0,n2,k) = undef 
	 ELSE
            pp3(m0,n1,k) = b2_3*ff(2) + a1_3*ff(3) 
            pp3(m0,n2,k) = a1_3*ff(2) + b2_3*ff(3) 
         ENDIF 
      ENDDO 
   
      DO j = 2, jm-2  ; j1 = j-1 ; j4 = j+2 
         n0 = 3*(j-1) + 1 ; n1 = n0+1 ; n2 = n0+2 
         ff(1:4) = pp(i,j1:j4,k) ; pp3(m0,n0,k) = ff(2) 
	 aaa = min(ff(2), ff(3)) ; bbb = max(ff(2), ff(3)) 
         IF ( MaxVal(ff(:)) == undef ) THEN 
            pp3(m0,n1,k) = undef ; pp3(m0,n2,k) = undef 
	 ELSE
            xxx = SUM(awt(:)*ff(:)) ; pp3(m0,n1,k) = max( aaa, min(xxx, bbb) ) 
            yyy = SUM(bwt(:)*ff(:)) ; pp3(m0,n2,k) = max( aaa, min(yyy, bbb) ) 
         ENDIF 
      ENDDO 
   
   ENDDO 
   
!* Refine the grid in X-direction 
!* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   DO n = 1, jm3 
   
      DO i = 1, im-1, im-2 ; m0 = 3*(i-1) + 1 
	 m1 = m0+1 ; m2 = m0+2 ; mm2 = m0 ; mm3 = m0+3 
         ff(2) = pp3(mm2,n,k) ; ff(3) = pp3(mm3,n,k) 
         IF ( max(ff(2),ff(3)) == undef ) THEN 
            pp3(m1,n,k) = undef ; pp3(m2,n,k) = undef 
	 ELSE
            pp3(m1,n,k) = b2_3*ff(2) + a1_3*ff(3) 
            pp3(m2,n,k) = a1_3*ff(2) + b2_3*ff(3) 
         ENDIF 
      ENDDO 
   
      DO i = 2, im-2 
         m0 = 3*(i-1) + 1 ; m1 = m0+1 ; m2 = m0+2 
         mm1 = m0-3 ; mm2 = m0 ; mm3 = m0+3 ; mm4 = m0+6 
         ff(1) =  pp3(mm1,n,k) ; ff(2) =  pp3(mm2,n,k) 
         ff(3) =  pp3(mm3,n,k) ; ff(4) =  pp3(mm4,n,k) 
	 aaa = min(ff(2), ff(3)) ; bbb = max(ff(2), ff(3)) 
         IF ( MaxVal(ff(:)) == undef ) THEN 
            pp3(m1,n,k) = undef ; pp3(m2,n,k) = undef 
	 ELSE
	    xxx = SUM(awt(:)*ff(:)) ; pp3(m1,n,k) = max( aaa, min(xxx, bbb) ) 
	    yyy = SUM(bwt(:)*ff(:)) ; pp3(m2,n,k) = max( aaa, min(yyy, bbb) ) 
         ENDIF 
      ENDDO 
   
   ENDDO 
   
   ENDDO !* <-----------------------------------------------------------
   
   END SELECT  !========================================================== 

   END SUBROUTINE Refine_Grid3x 
   
!==============================================================================
   
!******************************************************************************
! 
!  Convert a scalar variable from E to A-grid 
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  This subroutine converts a scalar variable from E to A-grid with doubled 
!  grid points assuming that im = 2*ime-1, jm = jme (ime, jme = odd numbers). 
!  
!______________________________________________________________________________
   
   SUBROUTINE E2A_Grid_Scalar( hh, im, jm, km, ee, ime, jme, kme ) 
   
   Implicit None 
!  -------------------------------------------------------------------- [input] 
   Integer  im , jm , km     !* x, y, z-dimension of target A-grid 
   Integer  ime, jme, kme    !* x, y, z-dimension of source E-grid 
   Real(4)  ee(ime,jme,kme)  !* original data on the source E-grid 
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - [output] 
   Real(4)  hh(im,jm,km)     !* extended data on the target A-grid 
!  -------------------------------------------------------------------- [local] 
   Integer  ii, i, j, k 
   
!           *     *     *     *     *     *     *     *     *     * 
   
   DO k = 1, km !* ----------------------------------------------------> 
   
   DO j = 2, jm-1, 2 
      DO i = 2, im-1, 2 ; ii = i/2 
         hh(i,j,k) = ee(ii,j,k) 
      ENDDO 
      i = 1  ; ii = (i+1)/2 
      hh(i,j,k) = .50*( ee(ii,j-1,k) + ee(ii,j+1,k) ) 
      DO i = 3, im-2, 2 ; ii = (i+1)/2 
         hh(i,j,k) = .25*( ee(ii-1,j,k)+ee(ii,j,k)+ee(ii,j-1,k)+ee(ii,j+1,k) ) 
      ENDDO 
      i = im ; ii = (i+1)/2 
      hh(i,j,k) = .50*( ee(ii,j-1,k) + ee(ii,j+1,k) ) 
   ENDDO 
   
      j = 1 
      DO i = 1, im, 2 ; ii = (i+1)/2 
         hh(i,j,k) = ee(ii,j,k) 
      ENDDO 
      DO i = 2, im-1, 2 ; ii = i/2 
         hh(i,j,k) = .50*( ee(ii+1,j,k) + ee(ii,j,k) ) 
      ENDDO 
   
   DO j = 3, jm-2, 2 
      DO i = 1, im, 2 ; ii = (i+1)/2 
         hh(i,j,k) = ee(ii,j,k) 
      ENDDO 
      DO i = 2, im-1, 2 ; ii = i/2 
         hh(i,j,k) = .25*( ee(ii+1,j,k)+ee(ii,j,k)+ee(ii,j-1,k)+ee(ii,j+1,k) ) 
      ENDDO 
   ENDDO 
   
      j = jm 
      DO i = 1, im, 2 ; ii = (i+1)/2 
         hh(i,j,k) = ee(ii,j,k) 
      ENDDO 
      DO i = 2, im-1, 2 ; ii = i/2 
         hh(i,j,k) = .50*( ee(ii+1,j,k) + ee(ii,j,k) ) 
      ENDDO 

   ENDDO !* ( k = 1, km ) <---------------------------------------------
   
   END SUBROUTINE E2A_Grid_Scalar 
   
!==============================================================================

!******************************************************************************
! 
!  Convert a vector component from E to A-grid 
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  This subroutine converts a vector component from E to A-grid with doubled 
!  grid points assuming that im = 2*ime-1, jm = jme (ime, jme = odd numbers). 
!  
!______________________________________________________________________________
   
   SUBROUTINE E2A_Grid_Vector( uu, im, jm, km, ee, ime, jme, kme ) 
   
   Implicit None 
!  -------------------------------------------------------------------- [input] 
   Integer  im , jm , km     !* x, y, z-dimension of target A-grid 
   Integer  ime, jme, kme    !* x, y, z-dimension of source E-grid 
   Real(4)  ee(ime,jme,kme)  !* original data on the source E-grid 
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - [output] 
   Real(4)  uu(im,jm,km)     !* extended data on the target A-grid 
!  -------------------------------------------------------------------- [local] 
   Integer  ii, i, j, k 
   
!           *     *     *     *     *     *     *     *     *     * 
   
   DO k = 1, km !* ----------------------------------------------------> 
   
   DO j = 2, jm-1, 2 
      DO i = 1, im, 2 ; ii = (i+1)/2 
         uu(i,j,k) = ee(ii,j,k) 
      ENDDO 
      DO i = 2, im-1, 2 ; ii = i/2 
         uu(i,j,k) = .25*( ee(ii+1,j,k)+ee(ii,j,k)+ee(ii,j-1,k)+ee(ii,j+1,k) ) 
      ENDDO 
   ENDDO 
   
      j = 1 
      DO i = 2, im-1, 2 ; ii = i/2 
         uu(i,j,k) = ee(ii,j,k) 
      ENDDO 
      i = 1  ; ii = (i+1)/2 
      uu(i,j,k) = .25*( 3*ee(ii,j,k)-ee(ii+1,j,k)+3*ee(ii,j+1,k)-ee(ii,j+3,k) ) 
      DO i = 3, im-2, 2 ; ii = (i+1)/2 
         uu(i,j,k) = .50*( ee(ii-1,j,k) + ee(ii,j,k) ) 
      ENDDO 
      i = im ; ii = (i+1)/2 
      uu(i,j,k) = .25*(3*ee(ii-1,j,k)-ee(ii-2,j,k)+3*ee(ii,j+1,k)-ee(ii,j+3,k)) 
   
   DO j = 3, jm-2, 2 
      DO i = 2, im-1, 2 ; ii = i/2 
         uu(i,j,k) = ee(ii,j,k) 
      ENDDO 
      i = 1  ; ii = (i+1)/2 
      uu(i,j,k) = .50*( ee(ii,j-1,k) + ee(ii,j+1,k) ) 
      DO i = 3, im-2, 2 ; ii = (i+1)/2 
         uu(i,j,k) = .25*( ee(ii-1,j,k)+ee(ii,j,k)+ee(ii,j-1,k)+ee(ii,j+1,k) ) 
      ENDDO 
      i = im ; ii = (i+1)/2 
      uu(i,j,k) = .50*( ee(ii,j-1,k) + ee(ii,j+1,k) ) 
   ENDDO 
   
      j = jm
      DO i = 2, im-1, 2 ; ii = i/2 
         uu(i,j,k) = ee(ii,j,k) 
      ENDDO 
      i = 1  ; ii = (i+1)/2 
      uu(i,j,k) = .25*( 3*ee(ii,j,k)-ee(ii+1,j,k)+3*ee(ii,j-1,k)-ee(ii,j-3,k) ) 
      DO i = 3, im-2, 2 ; ii = (i+1)/2 
         uu(i,j,k) = .50*( ee(ii-1,j,k) + ee(ii,j,k) ) 
      ENDDO 
      i = im ; ii = (i+1)/2 
      uu(i,j,k) = .25*(3*ee(ii-1,j,k)-ee(ii-2,j,k)+3*ee(ii,j-1,k)-ee(ii,j-3,k)) 
   
   ENDDO !* ( k = 1, km ) <---------------------------------------------
   
   END SUBROUTINE E2A_Grid_Vector 
   
!==============================================================================
   
!******************************************************************************
! 
!  Convert a mis-located vector component from E to A-grid 
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  This subroutine converts a mis-located vector component from E to A-grid 
!  with doubled grid points assuming im = 2*ime-1, jm = jme (ime, jme = odd). 
!  
!______________________________________________________________________________
   
   SUBROUTINE E2A_Debug_Vector( hh, im, jm, km, ee, ime, jme, kme ) 
   
   Implicit None 
!  -------------------------------------------------------------------- [input] 
   Integer  im , jm , km     !* x, y, z-dimension of target A-grid 
   Integer  ime, jme, kme    !* x, y, z-dimension of source E-grid 
   Real(4)  ee(ime,jme,kme)  !* original data on the source E-grid 
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - [output] 
   Real(4)  hh(im,jm,km)     !* extended data on the target A-grid 
!  -------------------------------------------------------------------- [local] 
   Integer  ii, i, j, k 
   
!           *     *     *     *     *     *     *     *     *     * 
   
   DO k = 1, km !* ----------------------------------------------------> 
   
   DO j = 2, jm-1, 2 
      DO i = 2, im-1, 2 ; ii = i/2 
         hh(i,j,k) = ee(ii,j,k) 
      ENDDO 
      DO i = 3, im-2, 2 ; ii = (i+1)/2 
         hh(i,j,k) = .25*( ee(ii-1,j,k)+ee(ii,j,k)+ee(ii,j-1,k)+ee(ii,j+1,k) ) 
      ENDDO 
      i = im ; ii = (i+1)/2 
      hh(i,j,k) = .50*( ee(ii,j-1,k) + ee(ii,j+1,k) ) 
   ENDDO 
   
      j = 1 
      DO i = 3, im, 2 ; ii = (i+1)/2 
         hh(i,j,k) = ee(ii,j,k) 
      ENDDO 
      i = 2 ; ii = i/2 
      hh(i,j,k) = .25*(3*ee(ii+1,j,k)-ee(ii+2,j,k)+3*ee(ii,j+1,k)-ee(ii,j+3,k)) 
      DO i = 4, im-1, 2 ; ii = i/2 
         hh(i,j,k) = .50*( ee(ii+1,j,k) + ee(ii,j,k) ) 
      ENDDO 
   
   DO j = 3, jm-2, 2 
      DO i = 3, im, 2 ; ii = (i+1)/2 
         hh(i,j,k) = ee(ii,j,k) 
      ENDDO 
      i = 2 ; ii = i/2 
      hh(i,j,k) = .50*( ee(ii,j-1,k) + ee(ii,j+1,k) ) 
      DO i = 4, im-1, 2 ; ii = i/2 
         hh(i,j,k) = .25*( ee(ii+1,j,k)+ee(ii,j,k)+ee(ii,j-1,k)+ee(ii,j+1,k) ) 
      ENDDO 
   ENDDO 
   
      j = jm 
      DO i = 3, im, 2 ; ii = (i+1)/2 
         hh(i,j,k) = ee(ii,j,k) 
      ENDDO 
      i = 2 ; ii = i/2 
      hh(i,j,k) = .25*(3*ee(ii+1,j,k)-ee(ii+2,j,k)+3*ee(ii,j-1,k)-ee(ii,j-3,k)) 
      DO i = 4, im-1, 2 ; ii = i/2 
         hh(i,j,k) = .50*( ee(ii+1,j,k) + ee(ii,j,k) ) 
      ENDDO 
   
      DO j = 1, jm 
         hh(1,j,k) = 2*hh(2,j,k) - hh(3,j,k) 
      ENDDO 

   ENDDO !* ( k = 1, km ) <---------------------------------------------
   
   END SUBROUTINE E2A_Debug_Vector 

!
!* Generalized routine includes E2A and B2A. This is Thiago's doings
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!******************************************************************************
!
!  Destagger a scalar variable via interpolation from an arbitrary grid to A-grid
!
!  Author: Thiago.Quirino@noaa.gov, DOC/OAR/AOML/HRD
!______________________________________________________________________________

   SUBROUTINE Destagger_Grid_Scalar( hh, im, jm, km, ee, ime, jme, kme , igridtype )

   Implicit None
!  -------------------------------------------------------------------- [constants]
   integer,parameter :: GRID_TYPE_ROTATEDLLE=1  !Rotated Lat-Lon projection with Arakawa E-grid staggering.
   integer,parameter :: GRID_TYPE_ROTATEDLLB=2  !Rotated Lat-Lon projection with Arakawa B-grid staggering.
!  -------------------------------------------------------------------- [input]
   Integer,intent(in) :: igridtype              !* grid projection and staggering type.
   Integer,intent(in) :: im , jm , km           !* x, y, z-dimension of target A-grid
   Integer,intent(in) :: ime, jme, kme          !* x, y, z-dimension of source grid
   Real(4),intent(in) :: ee(ime,jme,kme)        !* original data on the source grid
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - [output]
   Real(4),intent(inout) :: hh(im,jm,km)        !* extended data on the target A-grid

   !Decide which interpolation subroutine to invoke based on the grid type.
   IF(igridtype==GRID_TYPE_ROTATEDLLB) THEN
       call B2A_Grid_Scalar( hh, im, jm, km, ee, ime, jme, kme )    !For rotated Lat-Lon projection on Arakawa B-grid staggering for NMMB.
   ELSE
       call E2A_Grid_Scalar( hh, im, jm, km, ee, ime, jme, kme )    !Default, rotated Lat-Lon projection on Arakawa E-grid staggering for HWRF.
   ENDIF

   END SUBROUTINE Destagger_Grid_Scalar

!==============================================================================

!******************************************************************************
!
!  Destagger a vector component via interpolation from an arbitrary grid to A-grid
!
!  Author: Thiago.Quirino@noaa.gov, DOC/OAR/AOML/HRD
!______________________________________________________________________________

   SUBROUTINE Destagger_Grid_Vector( hh, im, jm, km, ee, ime, jme, kme , igridtype )

   Implicit None
!  -------------------------------------------------------------------- [constants]
   integer,parameter :: GRID_TYPE_ROTATEDLLE=1  !Rotated Lat-Lon projection with Arakawa E-grid staggering.
   integer,parameter :: GRID_TYPE_ROTATEDLLB=2  !Rotated Lat-Lon projection with Arakawa B-grid staggering.
!  -------------------------------------------------------------------- [input]
   Integer,intent(in) :: igridtype              !* grid projection and staggering type.
   Integer,intent(in) :: im , jm , km           !* x, y, z-dimension of target A-grid
   Integer,intent(in) :: ime, jme, kme          !* x, y, z-dimension of source grid
   Real(4),intent(in) :: ee(ime,jme,kme)        !* original data on the source grid
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - [output]
   Real(4),intent(inout) :: hh(im,jm,km)        !* extended data on the target A-grid

   !Decide which interpolation subroutine to invoke based on the grid type.
   IF(igridtype==GRID_TYPE_ROTATEDLLB) THEN
       call B2A_Grid_Vector( hh, im, jm, km, ee, ime, jme, kme )    !For rotated Lat-Lon projection on Arakawa B-grid staggering for NMMB.
   ELSE
       call E2A_Grid_Vector( hh, im, jm, km, ee, ime, jme, kme )    !Default, rotated Lat-Lon projection on Arakawa E-grid staggering for HWRF.
   ENDIF

   END SUBROUTINE Destagger_Grid_Vector

!==============================================================================
   
!******************************************************************************
!
!  Convert a scalar variable from B to A-grid
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  This subroutine converts a scalar variable from B to A-grid. It assumes that
!  im = ime, jm = jme, km=kme.
!
!  Author: Thiago.Quirino@noaa.gov, DOC/OAR/AOML/HRD
!______________________________________________________________________________

   SUBROUTINE B2A_Grid_Scalar( hh, im, jm, km, ee, ime, jme, kme )

   Implicit None
!  -------------------------------------------------------------------- [input]
   Integer  im , jm , km     !* x, y, z-dimension of target A-grid
   Integer  ime, jme, kme    !* x, y, z-dimension of source E-grid
   Real(4)  ee(ime,jme,kme)  !* original data on the source E-grid
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - [output]
   Real(4)  hh(im,jm,km)     !* extended data on the target A-grid
!  -------------------------------------------------------------------- [local]
   Integer  ii, i, j, k

!           *     *     *     *     *     *     *     *     *     *
   hh=ee    !B-grid scalar variables (mass grid) are assumed to be unstaggered.

   END SUBROUTINE B2A_Grid_Scalar

!==============================================================================

!******************************************************************************
!
!  Convert a vector component from B to A-grid
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  This subroutine converts a vector component from B to A-grid. It assumes that
!  im = ime, jm = jme, km=kme. It interpolates the wind grid onto the mass grid.
!
!  Author: Thiago.Quirino@noaa.gov, DOC/OAR/AOML/HRD
!______________________________________________________________________________

   SUBROUTINE B2A_Grid_Vector( hh, im, jm, km, ee, ime, jme, kme )

   Implicit None
!  -------------------------------------------------------------------- [input]
   Integer  im , jm , km     !* x, y, z-dimension of target A-grid
   Integer  ime, jme, kme    !* x, y, z-dimension of source E-grid
   Real(4)  ee(ime,jme,kme)  !* original data on the source E-grid
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - [output]
   Real(4)  hh(im,jm,km)     !* extended data on the target A-grid
!  -------------------------------------------------------------------- [local]
   Integer  i, j, k

!           *     *     *     *     *     *     *     *     *     *
   DO k = 1, km !* ---------------------------------------------------->
       !Interpolate non-border grid points.
       DO j = 2, jm-1
           DO i = 2, im-1
               hh(i,j,k) = .25*( ee(i-1,j-1,k)+ee(i-1,j,k)+ee(i,j-1,k)+ee(i,j,k) )
           ENDDO
       ENDDO

       !Interpolate non-corner points in first row
       j=1
       DO i = 2, im-1
           hh(i,j,k) = .25*( 3*ee(i-1,j,k)-ee(i-1,j+1,k) + 3*ee(i,j,k)-ee(i,j+1,k) )
       ENDDO

       !Interpolate non-corner points in last row
       j=jm
       DO i = 2, im-1
           hh(i,j,k) = .25*( 3*ee(i-1,j-1,k)-ee(i-1,j-2,k) + 3*ee(i,j-1,k)-ee(i,j-2,k) )
       ENDDO

       !Interpolate non-corner points in first column
       i=1
       DO j = 2, jm-1
           hh(i,j,k) = .25*( 3*ee(i,j-1,k)-ee(i+1,j-1,k) + 3*ee(i,j,k)-ee(i+1,j,k) )
       ENDDO

       !Interpolate non-corner points in last column
       i=im
       DO j = 2, jm-1
           hh(i,j,k) = .25*( 3*ee(i-1,j-1,k)-ee(i-2,j-1,k) + 3*ee(i-1,j,k)-ee(i-2,j,k) )
       ENDDO

       !Interpolate the SW corner point
       i=1 ; j=1
       hh(i,j,k) = 0.25*( 3*hh(i+1,j,k)-hh(i+2,j,k) + 3*hh(i,j+1,k)-hh(i,j+2,k) )

       !Interpolate the SE corner point
       i=im ; j=1
       hh(i,j,k) = 0.25*( 3*hh(i-1,j,k)-hh(i-2,j,k) + 3*hh(i,j+1,k)-hh(i,j+2,k) )

       !Interpolate the NW corner point
       i=1 ; j=jm
       hh(i,j,k) = 0.25*( 3*hh(i+1,j,k)-hh(i+2,j,k) + 3*hh(i,j-1,k)-hh(i,j-2,k) )

       !Interpolate the NE corner point
       i=im ; j=jm
       hh(i,j,k) = 0.25*( 3*hh(i-1,j,k)-hh(i-2,j,k) + 3*hh(i,j-1,k)-hh(i,j-2,k) )

   ENDDO !* ( k = 1, km ) <---------------------------------------------

   END SUBROUTINE B2A_Grid_Vector

!==============================================================================

   
!==============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 

MODULE read_util_module

#ifdef crayx1
#define iargc ipxfargc
#endif

CONTAINS !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#ifdef crayx1
   SUBROUTINE getarg(i, harg)
      implicit none
      character(len=*) :: harg
      integer :: ierr, ilen, i
      CALL pxfgetarg(i, harg, ilen, ierr)
      return
   END SUBROUTINE getarg
#endif

!------------------------------------------------------------------------------

   SUBROUTINE arguments(v2file, lmore)

   implicit none
   character(len=*) :: v2file
   character(len=120) :: harg
   logical :: lmore
   integer :: ierr, i, numarg
   integer, external :: iargc

   numarg = iargc()

   i = 1
   lmore = .false.

   DO while ( i < numarg)
      CALL getarg(i, harg)
      print*, 'harg = ', trim(harg)
      IF (harg == "-v") THEN
         i = i + 1
         lmore = .true.
      ELSEIF (harg == "-h") THEN
         CALL help
      ENDIF
   ENDDO

   CALL getarg(i,harg)
   v2file = harg

   END SUBROUTINE arguments

!------------------------------------------------------------------------------

   SUBROUTINE help

   implicit none
   character(len=120) :: cmd
   call getarg(0, cmd)

   write(*,'(/,"Usage: ", A, " [-v] v2file ")') trim(cmd)
   write(*,'(8x, "-v     : Print extra info")')
   write(*,'(8x, "v3file : MM5v3 file name to read.")')
   write(*,'(8x, "-h     : print this help message and exit.",/)')
   STOP

   END SUBROUTINE help

END MODULE read_util_module  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!==============================================================================

logical function wrf_dm_on_monitor()
  wrf_dm_on_monitor=.true.
end function wrf_dm_on_monitor

LOGICAL FUNCTION iveceq( a, b, n )
   implicit none
   integer n
   integer a(n), b(n)
   integer i
   iveceq = .true.
   do i = 1,n
      if ( a(i) .ne. b(i) ) iveceq = .false.
   enddo
   return
END FUNCTION iveceq

!* Stop for routine called by module_wrf_error
!* (used by netcdf implementation of IO api)

SUBROUTINE wrf_abort
   STOP
END SUBROUTINE wrf_abort

SUBROUTINE get_current_time_string( time_str )
   CHARACTER(LEN=*), INTENT(OUT) :: time_str
   time_str = ''
END SUBROUTINE get_current_time_string

SUBROUTINE get_current_grid_name( grid_str )
   CHARACTER(LEN=*), INTENT(OUT) :: grid_str
   grid_str = ''
END SUBROUTINE get_current_grid_name

!============================================================================== 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
