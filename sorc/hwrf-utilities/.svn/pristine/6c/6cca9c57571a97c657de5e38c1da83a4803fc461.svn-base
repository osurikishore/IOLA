#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <ctype.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

/* Convenience macro for defining synonyms of c_stat: */
#define DEF_C_STAT(CSNEW) \
void CSNEW(const char *filename,uint64_t *lstatflag,uint64_t *ierr, \
  uint64_t *quiet, uint64_t *dev, uint64_t *ino, uint64_t *mode, \
  uint64_t *justperms, uint64_t *nlink, uint64_t *uid, \
  uint64_t *gid, uint64_t *rdev, uint64_t *size, uint64_t *blksize,\
  uint64_t *blocks, double *atime, double *mtime, double *ctime,\
  uint64_t *flags) {\
 c_stat(filename,lstatflag,ierr,quiet,dev,ino,mode,justperms,nlink,uid,\
        gid,rdev,size,blksize,blocks,atime,mtime,ctime,flags); \
}

void c_stat(const char *filename,uint64_t *lstatflag,uint64_t *ierr,
            uint64_t *quiet, uint64_t *dev, uint64_t *ino, uint64_t *mode,
            uint64_t *justperms, uint64_t *nlink, uint64_t *uid, 
            uint64_t *gid, uint64_t *rdev, uint64_t *size, uint64_t *blksize,
            uint64_t *blocks, double *atime, double *mtime, double *ctime,
            uint64_t *flags) {
  struct stat s;
  int iflag;
  char errbuf[100];
  *quiet=0;

  if(*lstatflag) {
  if(!*quiet)
    if(lstat(filename,&s)) {
      if(!*quiet)  {
        errbuf[0]='\0';
        strerror_r(errno,errbuf,100);
        fprintf(stderr,"%s: cannot lstat: %s\n",filename,errbuf);
      }
      *ierr=errno;
      if(*ierr==0) *ierr=-1;
      return;
    }
  } else { 
  if(!*quiet)
    if(stat(filename,&s)) {
      if(!*quiet) {
        errbuf[0]='\0';
        strerror_r(errno,errbuf,100);
        fprintf(stderr,"%s: cannot stat: %s\n",filename,errbuf);
      }
      *ierr=errno;
      if(*ierr==0) *ierr=-1;
      return;
    }
  }
  
  *ierr=0;

  *dev=(uint64_t)s.st_dev;
  *ino=(uint64_t)s.st_ino;
  *mode=(uint64_t)s.st_mode;
  *justperms=(uint64_t)(s.st_mode & 07777);
  *nlink=(uint64_t)s.st_nlink;
  *uid=(uint64_t)s.st_uid;
  *gid=(uint64_t)s.st_gid;
  *rdev=(uint64_t)s.st_rdev;
  *size=(uint64_t)s.st_size;
  *blksize=(uint64_t)s.st_blksize;
  *blocks=(uint64_t)s.st_blocks;
  *atime=(double)s.st_atime;
  *mtime=(double)s.st_mtime;
  *ctime=(double)s.st_ctime;

  iflag=0;
  if(S_ISREG(s.st_mode))  iflag|=1<<0;
  if(S_ISDIR(s.st_mode))  iflag|=1<<1;
  if(S_ISCHR(s.st_mode))  iflag|=1<<2;
  if(S_ISBLK(s.st_mode))  iflag|=1<<3;
  if(S_ISFIFO(s.st_mode)) iflag|=1<<4;
  if(S_ISLNK(s.st_mode))  iflag|=1<<5;
  if(S_ISSOCK(s.st_mode)) iflag|=1<<6;

  *flags = iflag;
}

/* Define synonyms using all common Fortran name mangling schemes */
DEF_C_STAT(c_stat_);
DEF_C_STAT(c_stat__);
DEF_C_STAT(C_STAT);
DEF_C_STAT(C_STAT_);
DEF_C_STAT(C_STAT__);
