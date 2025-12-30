MODULE read_gmv16_ncf

!=======================================================================

!$$$ PROGRAM DOCUMENTATION BLOCK
!
! Reading the GFS/GDAS V16 NetCDF Fields
! Email: biju.thomas@noaa.gov
!
!=======================================================================

USE netcdf
IMPLICIT NONE

CONTAINS
  SUBROUTINE get_dims(filename, nx, ny, nzf, nzh)
    CHARACTER(LEN = *), INTENT(IN) :: filename
    INTEGER, INTENT(OUT) :: nx, ny, nzf, nzh
    CHARACTER(LEN = 50) :: dname
    INTEGER :: ncid,dimlen
    INTEGER :: ndims, nvars, ngatts, unlimdimid
    LOGICAL :: file_exists

    INQUIRE(FILE=TRIM(filename), EXIST=file_exists)
    IF (.NOT. file_exists) THEN
      PRINT*, TRIM(filename),' NOT EXISTS'
      STOP
    ENDIF
    CALL check(NF90_OPEN(TRIM(filename), NF90_NOWRITE, ncid),'getdim:open')
    CALL check(NF90_INQUIRE(ncid,ndims,nvars,ngatts,unlimdimid),'getdim:inquire')
    PRINT*,'ndims = ',ndims
    IF(ndims < 4) THEN
      PRINT*, 'ndims < 4', ' ndims = ',ndims
      STOP
    ENDIF
    CALL check(NF90_INQUIRE_DIMENSION(ncid,1,dname,nx),'getdim:inquire_dimension 1')
    CALL check(NF90_INQUIRE_DIMENSION(ncid,2,dname,ny),'getdim:inquire_dimension 2')
    CALL check(NF90_INQUIRE_DIMENSION(ncid,3,dname,nzf),'getdim:inquire_dimension 3')
    CALL check(NF90_INQUIRE_DIMENSION(ncid,4,dname,nzh),'getdim:inquire_dimension 4')
    CALL check(NF90_CLOSE(ncid),'getdim:close')
  END SUBROUTINE get_dims

  SUBROUTINE get_attrlen(filename,attname,atttype,attlen)
    CHARACTER(LEN = *), INTENT(IN) :: filename, attname
    INTEGER, INTENT(OUT) :: attlen, atttype
    INTEGER :: ncid

    CALL check(NF90_OPEN(TRIM(filename), NF90_NOWRITE, ncid),'get_attrlen:open') 
    CALL check(NF90_INQUIRE_ATTRIBUTE(ncid, NF90_GLOBAL, attname,  &
                    xtype = atttype, len = attlen),'get_attrlen:inquire_attribute')
    CALL check(NF90_CLOSE(ncid),'get_attrlen:close')
  END SUBROUTINE get_attrlen

  SUBROUTINE get_attr(filename,attname,attlen,attval)
    CHARACTER(LEN = *), INTENT(IN) :: filename, attname
     INTEGER, INTENT(IN) :: attlen
    REAL, INTENT(OUT) :: attval(attlen)
    INTEGER :: ncid

    CALL check(NF90_OPEN(TRIM(filename), NF90_NOWRITE, ncid),'get_attrlen:open')
    CALL check(NF90_GET_ATT(ncid, NF90_GLOBAL, attname,  &
                    attval),'get_attr:inquire_get_att')
    CALL check(NF90_CLOSE(ncid),'get_attrlen:close')
  END SUBROUTINE get_attr

  SUBROUTINE get_coord(filename, nx, ny, nzf, nzh, xname, yname, zname, zhname, &
                                                   xgrid, ygrid, zgrid, zhgrid)
    CHARACTER(LEN = *), INTENT(IN) :: filename, xname, yname, zname, zhname
    INTEGER, INTENT(IN) :: nx, ny, nzf, nzh
    REAL, INTENT(OUT) :: xgrid(nx,ny), ygrid(nx,ny), zgrid(nzf), zhgrid(nzh)
    INTEGER :: ncid, varid

    CALL check(NF90_OPEN(TRIM(filename), NF90_NOWRITE, ncid),'get_coord:open')
    CALL check(NF90_INQ_VARID(ncid,TRIM(xname),varid),'get_coord:inq_varid x')
    CALL check(NF90_GET_VAR(ncid,varid,xgrid),'get_coord:get_var x')
    CALL check(NF90_INQ_VARID(ncid,TRIM(yname),varid),'get_coord:inq_varid y')
    CALL check(NF90_GET_VAR(ncid,varid,ygrid),'get_coord:get_var y')
    CALL check(NF90_INQ_VARID(ncid,TRIM(zname),varid),'get_coord:inq_varid z')
    CALL check(NF90_GET_VAR(ncid,varid,zgrid),'get_coord:get_var z')
    CALL check(NF90_INQ_VARID(ncid,TRIM(zhname),varid),'get_coord:inq_varid zh')
    CALL check(NF90_GET_VAR(ncid,varid,zhgrid),'get_coord:get_var zh')
    CALL check(NF90_CLOSE(ncid),'get_coord:close')
  END SUBROUTINE get_coord

  SUBROUTINE get_ncf2d(filename, nx, ny, fldname, field2d)
    CHARACTER(LEN = *), INTENT(IN) :: filename,fldname
    INTEGER, INTENT(IN) :: nx, ny
    REAL, INTENT(OUT) :: field2d(nx,ny)
    INTEGER :: ncid, varid

    CALL check(NF90_OPEN(TRIM(filename), NF90_NOWRITE, ncid),'read_ncf2d:open')
    CALL check(NF90_INQ_VARID(ncid,TRIM(fldname),varid),'read_ncf2d:inq_varid')
    CALL check(NF90_GET_VAR(ncid,varid,field2d),'read_ncf2d:get_var')
    CALL check(NF90_CLOSE(ncid),'read_ncf2d:close')
  END SUBROUTINE get_ncf2d

  SUBROUTINE get_ncf3d(filename, nx, ny, nz, fldname, field3d)
    CHARACTER(LEN = *), INTENT(IN) :: filename,fldname
    INTEGER, INTENT(IN) :: nx, ny, nz
    REAL, INTENT(OUT) :: field3d(nx,ny,nz)
    INTEGER :: ncid, varid

    CALL check(NF90_OPEN(TRIM(filename), NF90_NOWRITE, ncid),'read_ncf3d:open')
    CALL check(NF90_INQ_VARID(ncid,TRIM(fldname),varid),'read_ncf3d:inq_varid')
    CALL check(NF90_GET_VAR(ncid,varid,field3d),'read_ncf3d:get_var')
    CALL check(NF90_CLOSE(ncid),'read_ncf3d:close')
  END SUBROUTINE get_ncf3d

  SUBROUTINE check(index, mssg)
    CHARACTER(LEN = *), INTENT(IN) :: mssg
    INTEGER, INTENT(IN) :: index

    IF( index /= NF90_NOERR) THEN
      PRINT*, mssg, ' NF90_NOERR = ', NF90_STRERROR(index)
      PRINT*, mssg, ' EXITING'
      STOP
    ENDIF
  END SUBROUTINE check

END MODULE read_gmv16_ncf

