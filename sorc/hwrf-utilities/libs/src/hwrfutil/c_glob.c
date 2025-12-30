#include <glob.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

/* Return statuses */
static const int F_GLOB_NOSPACE=1, F_GLOB_ABORTED=2, F_GLOB_NOMATCH=4,
  F_GLOB_UNKNOWN_ERROR=8;

/* Flags */
static const int F_GLOB_ERR=1, F_GLOB_MARK=2, F_GLOB_NOSORT=4,
  F_GLOB_DOOFFS=8;

/* GNU Extension Flags */
static const int F_GLOB_PERIOD=16, F_GLOB_BRACE=32, F_GLOB_TILDE=64,
  F_GLOB_TILDE_CHECK=128, F_GLOB_ONLYDIR=256;

/* Local Extensions */
static const int F_QUIET=512;

/* #define C_GLOB_START_ALIAS(NAME)                \ */
/*   void NAME(const char *pattern,,uint64_t *flags,uint64_t *obj, \ */
/*             uint64_t *idx,uint64_t *offs,uint64_t *len,         \ */
/*             uint64_t *ierr) {                                   \ */
/*     c_glob_start(pattern,flags,obj,idx,offs,len,ierr);          \ */
/*   } */

void c_glob_free(void **obj) {
  glob_t *pglob=(glob_t*)*obj;
  if(pglob) {
    globfree(pglob);
    free(pglob);
    *obj=NULL;
  }
}

void c_glob_start(const char *pattern,uint64_t *flags,void **obj,
                   uint64_t *idx,uint64_t *offs,uint64_t *len,
                   uint64_t *ierr) {
  glob_t *pglob=(glob_t*)malloc(sizeof(glob_t));
  int err,cflag;

  *ierr=0;
  *idx=0;
  *len=0;
  *offs=0;

  /* Allocate space for globbing */
  if(!pglob) {
    if(!*flags&F_QUIET) 
      fprintf(stderr,"glob: cannot allocate %lld bytes of memory: %s\n",
              (long long)sizeof(glob_t),strerror(errno));
    *ierr=F_GLOB_NOSPACE;
    return;
  }

  /* Get the C copy of the fortran flags */
  cflag=0;
  if(*flags&F_GLOB_ERR)         cflag|=GLOB_ERR;
  if(*flags&F_GLOB_MARK)        cflag|=GLOB_MARK;
  if(*flags&F_GLOB_NOSORT)      cflag|=GLOB_NOSORT;
  if(*flags&F_GLOB_DOOFFS)      cflag|=GLOB_DOOFFS;

#ifdef GLOB_PERIOD
  if(*flags&F_GLOB_PERIOD)      cflag|=GLOB_PERIOD;
#endif

#ifdef GLOB_BRACE
  if(*flags&F_GLOB_BRACE)       cflag|=GLOB_BRACE;
#endif

#ifdef GLOB_TILDE
  if(*flags&F_GLOB_TILDE)       cflag|=GLOB_TILDE;
#endif

#ifdef GLOB_TILDE_CHECK
  if(*flags&F_GLOB_TILDE_CHECK) cflag|=GLOB_TILDE_CHECK;
#endif

#ifdef GLOB_ONLYDIR
  if(*flags&F_GLOB_ONLYDIR)     cflag|=GLOB_ONLYDIR;
#endif

  /* Now do some globbing */
  err=glob(pattern,cflag,NULL,pglob);
  if(err) {
    globfree(pglob);
    switch(err) {
    case(GLOB_NOSPACE):
      if(!*flags&F_QUIET)
        fprintf(stderr,"glob: cannot allocate memory: %s\n",strerror(errno));
      *ierr=F_GLOB_NOSPACE;
      return;
    case(GLOB_NOMATCH):
      if(!*flags&F_QUIET)
        fprintf(stderr,"glob: %s: no match\n",pattern);
      *ierr=F_GLOB_NOMATCH;
      return;
    case(GLOB_ABORTED):
      if(!*flags&F_QUIET)
        fprintf(stderr,"glob: %s: read error: %s\n",pattern,strerror(errno));
      *ierr=F_GLOB_ABORTED;
      return;
    default:
      if(!*flags&F_QUIET)
        fprintf(stderr,"glob: %s: unknown error code %d from glob\n",pattern,err);
      *ierr=F_GLOB_UNKNOWN_ERROR;
      return;
    }
  }
  
  /* Glob was successful.  Store results. */
  *obj=(void*)pglob;
  *idx=0;
  *len=(uint64_t)pglob->gl_pathc;
  *ierr=0;
}

void c_glob_len(glob_t *pglob,uint64_t *idx,int64_t *len) {
  size_t i=*idx - 1;
  if(i>=pglob->gl_pathc) {
    *len=0;
    return;
  }
  *len=strlen(pglob->gl_pathv[i]);
}

void c_glob_fill(glob_t *pglob,uint64_t *idx,uint64_t *len,char *buffer) {
  size_t ix=*idx-1, L=*len;
  size_t i;
  const char *str=pglob->gl_pathv[ix];
  for(i=0;i<L && str[i];i++) buffer[i]=str[i];
  for(   ;i<L          ;i++) 
    buffer[i]=' '; /* Just a safeguard; should never get here */
}

/* C_GLOB_START_ALIAS(c_glob_start_) */
/* C_GLOB_START_ALIAS(c_glob_start__) */
/* C_GLOB_START_ALIAS(C_GLOB_START) */
/* C_GLOB_START_ALIAS(C_GLOB_START_) */
/* C_GLOB_START_ALIAS(C_GLOB_START__) */
