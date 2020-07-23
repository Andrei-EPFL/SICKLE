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
#include <fftw3.h>
#include <math.h>
#include <time.h>
#include "define.h"
#include "load_conf.h"


typedef struct {
  double x[3];
} HALOS;

int read_pk(const char *, double **, double **, size_t *);

int init_pk(double *, double *, const size_t, const CONF *);

int init_halos(const long int, HALOS **);

int save_halo(const char *, HALOS *, const size_t);


double return_x(double );
void function_test();

#ifdef DOUBLE_PREC

int init_field(const int, fftw_complex **, fftw_plan *);

int gauss_ran_field(const CONF *, const double *, const double *,
    const size_t, fftw_plan *, fftw_complex *);

void qsort_dens(fftw_complex *, const size_t);

void select_dens(fftw_complex *, HALOS *, const int, const size_t, const double);

int save_dens(const char *, fftw_complex *, const int);

#else

int init_field(const int, fftwf_complex **, fftwf_plan *);

int gauss_ran_field(const CONF *, const double *, const double *,
    const size_t, fftwf_plan *, fftwf_complex *);

void qsort_dens(fftwf_complex *, const size_t);

void select_dens(fftwf_complex *, HALOS *, const int, const size_t, const double);

int save_dens(const char *, fftwf_complex *, const int);

#endif

size_t cnt_strcpy(char *, const char *, const size_t);

#endif
