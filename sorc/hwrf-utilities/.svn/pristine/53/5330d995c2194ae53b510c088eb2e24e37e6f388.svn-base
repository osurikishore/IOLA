module distance_module
  implicit none
  private
  public :: greatarc

contains

  real function greatarc(lat1,lon1,lat2,lon2)
    ! greatarc -- gets the great arc distance (along-Earth distance)
    !       between two points.  The earth radius used in the
    !       calculation is the average of the radius at the two
    !       points, using Requator and flattening from the
    !       constants_module.
    !     lat1,lon1 -- the first point's latitude
    !     lat2,lon2 -- the second point's latitude
    !     returns -- the distance in meters.
    use constants_module, only: Requator,pi,DEGRAD,flattening
    implicit none
    real, intent(in) :: lat1,lon1, lat2,lon2
    real :: rlat1,rlon1, rlat2,rlon2
    real :: Rearth1,Rearth2
    real, parameter :: deg2rad=DEGRAD
    real, parameter :: flattening_inv=1/flattening

    rlat1=lat1*deg2rad  ;  rlon1=lon1*deg2rad
    rlat2=lat2*deg2rad  ;  rlon2=lon2*deg2rad

    Rearth1=Requator*(1-sin(rlat1)**2/flattening_inv)
    Rearth2=Requator*(1-sin(rlat2)**2/flattening_inv)

    greatarc=(Rearth1+Rearth2)*asin(min(1.0,sqrt( &
         sin((rlat1-rlat2)/2)**2+ &
         cos(rlat1)*cos(rlat2)*sin((rlon1-rlon2)/2)**2)))
  end function greatarc
end module distance_module
