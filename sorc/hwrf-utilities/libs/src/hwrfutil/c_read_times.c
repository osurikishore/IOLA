#include <stdlib.h>
#include <netcdf.h>
#include <stdint.h>
#include <stdio.h>

void c_read_times(uint64_t *ncid, uint64_t *varid, int32_t *times_i1,uint64_t *ntimes,uint64_t *timelen,uint64_t *err) {
  const size_t start[5]={ 1,1,1,1,0 };
  size_t counts[5]={ 18,1,1,1,1};
  fprintf(stderr,"%lld %lld %lld %lld\n",start[0], start[1],counts[0],counts[1]);
  *err=nc_get_vara_int(*ncid,*varid,start,counts,times_i1);
  printf("read  %d\n",times_i1[0]);
}
void c_read_times_(uint64_t *ncid, uint64_t *varid, int32_t *times_i1,uint64_t *ntimes,uint64_t *timelen,uint64_t *err) {
  return c_read_times(ncid,varid,times_i1,ntimes,timelen,err);
}
void c_read_times__(uint64_t *ncid, uint64_t *varid, int32_t *times_i1,uint64_t *ntimes,uint64_t *timelen,uint64_t *err) {
  return c_read_times(ncid,varid,times_i1,ntimes,timelen,err);
}
void C_READ_TIMES(uint64_t *ncid, uint64_t *varid, int32_t *times_i1,uint64_t *ntimes,uint64_t *timelen,uint64_t *err) {
  return c_read_times(ncid,varid,times_i1,ntimes,timelen,err);
}
void C_READ_TIMES_(uint64_t *ncid, uint64_t *varid, int32_t *times_i1,uint64_t *ntimes,uint64_t *timelen,uint64_t *err) {
  return c_read_times(ncid,varid,times_i1,ntimes,timelen,err);
}
void C_READ_TIMES__(uint64_t *ncid, uint64_t *varid, int32_t *times_i1,uint64_t *ntimes,uint64_t *timelen,uint64_t *err) {
  return c_read_times(ncid,varid,times_i1,ntimes,timelen,err);
}
