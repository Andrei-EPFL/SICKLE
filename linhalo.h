/**********************************************************
**                                                       **
**      Generate a linear halo catalogue                 **
**      Author: Cheng Zhao <zhaocheng03@gmail.com>       **
**                                                       **
**********************************************************/

#ifndef _LINHALO_H_
#define _LINHALO_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <fftw3.h>
#include <math.h>
#include <time.h>
#include "fftw_define.h"
#include "define.h"
#include "load_conf.h"



typedef struct {
  double x[3];
  size_t dens;
} HALOS;

typedef struct {
  size_t val;
  size_t rank;
  size_t index;
} INDEX;

int read_pk(const char *, double **, double **, size_t *);

int init_pk(double *, double *, const size_t, const CONF *);

int init_halos(const long int, HALOS **);

int init_index_arr(const long int, INDEX **);

int save_halo(const char *, HALOS *, const size_t);

double return_x(double );

int init_field(const int, FFT_CMPLX **, FFT_PLAN *);

int gauss_ran_field(const CONF *, const double *, const double *,
    const size_t, FFT_PLAN *, FFT_CMPLX *);

void qsort_dens_asc(FFT_CMPLX *, INDEX *, const size_t, int);

void qsort_dens_desc(FFT_CMPLX *, INDEX *, const size_t, int);

void qsort_idx_desc(INDEX *, const size_t );

void select_dens(FFT_CMPLX *, HALOS *, INDEX *, const int, const size_t, const double);

int save_dens(const char *, FFT_CMPLX *, const int);

size_t binary_search(INDEX *, size_t, size_t, size_t, size_t);

size_t binary_search_pos(FFT_CMPLX *, INDEX *, long int, long int, FFT_REAL, size_t);

size_t cnt_strcpy(char *, const char *, const size_t);

#endif
