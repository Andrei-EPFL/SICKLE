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
  FFT_REAL dens;  // density value
  size_t idx;  // index of the density value in the density field
} MAX_DENS;

int read_pk(const char *, double **, double **, size_t *);

int init_pk(double *, double *, const size_t, const CONF *);

int init_halos(const long int, HALOS **);

int init_max_dens(const long int, MAX_DENS **);

int save_halo(const char *, char *, HALOS *, const size_t);

double return_x(double );

int init_field(const int, FFT_CMPLX **, FFT_PLAN *);

int gauss_ran_field(const CONF *, const double *, const double *,
    const size_t, FFT_PLAN *, FFT_CMPLX *, const double);

void qsort_dens_asc(MAX_DENS *, const size_t);

void qsort_dens_desc(MAX_DENS *, const size_t);

void select_dens(FFT_CMPLX *, HALOS *, MAX_DENS *, const int, const size_t, const double);

int save_dens(const char *, FFT_CMPLX *, const int);

size_t binary_search(FFT_CMPLX *, size_t *, size_t, size_t, FFT_REAL, size_t);

size_t binary_search_pos(FFT_CMPLX *, size_t *, size_t, size_t, FFT_REAL, size_t);

size_t cnt_strcpy(char *, const char *, const size_t);

#endif
