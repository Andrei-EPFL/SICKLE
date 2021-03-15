#include "linhalo.h"
#include "fftw_define.h"
//#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_spline.h>
#ifdef OMP
#include <omp.h>
#endif

int init_pk(double *k, double *P, const size_t Nk, const CONF *conf) {
  size_t i;

  if (k[0] > 2 * M_PI / conf->Lbox) {
    P_ERR("the minimum k of the input power spectrum is too large.\n");
    return ERR_RANGE;
  }
  if (k[Nk - 1] < M_PI * conf->Ngrid / conf->Lbox) {
    P_ERR("the maximum k of the input power spectrum is too small.\n");
    return ERR_RANGE;
  }

  for (i = 0; i < Nk; i++) {
    k[i] = log(k[i]);
    P[i] = log(P[i]);
  }

  printf("  Power spectrum interpolated.\n");
  return 0;
}

int init_field(const int Ngrid, FFT_CMPLX **mesh, FFT_PLAN *plan)
{
  size_t size = (size_t) Ngrid * Ngrid * Ngrid;

  printf("  Initialising the density field mesh... \n");
  fflush(stdout);

  *mesh = FFT_MALLOC(sizeof(FFT_CMPLX) * size);
  if (!(*mesh)) {
    P_ERR("failed to allocate memory for the density field.\n");
    return ERR_MEM;
  }

#ifdef OMP
  if (FFT_INIT_OMP() == 0) {
    P_ERR("failed to initialize FFTW with OpenMP.\n");
    return ERR_OTHER;
  }
  FFT_PLAN_OMP(omp_get_max_threads());
#endif
  *plan = FFT_PLAN_C2C(Ngrid, Ngrid, Ngrid, *mesh, *mesh,
      FFTW_BACKWARD, FFTW_ESTIMATE);

  printf("\r  The density field is initialised.\n");
  return 0;
}

int init_halos(const long int Nh, HALOS **halos) {
  printf("  Initialising halos... \n");
  fflush(stdout);
  
  //(ptr) = (type *) malloc(sizeof(type) * (n));
  *halos = malloc( sizeof(HALOS) * Nh );
  
  if (!(*halos)) {
    P_ERR("failed to allocate memory for the halos.\n");
    return ERR_MEM;
  }
  return 0;
}

int init_max_dens(const long int Nh,  MAX_DENS **max_dens) {
  printf("  Initialising array of indices for the %ld largest density values ... \n", 8 * Nh);
  fflush(stdout);
  
  //(ptr) = (type *) malloc(sizeof(type) * (n));
  *max_dens = malloc( sizeof(MAX_DENS) * 8 * Nh );
  
  if (!(*max_dens)) {
    P_ERR("failed to allocate memory for the array of indices.\n");
    return ERR_MEM;
  }
  return 0;
}

int gauss_ran_field(const CONF *conf, const double *lnk, const double *lnP,
    const size_t Nk, FFT_PLAN *fp, FFT_CMPLX *mesh, const double factor) {
  int i, j, k, i_2, j_2, k_2;
  size_t idx, idx_2;
  double fac, norm, fac2, ki, kj, kk, ksq, kv, P, theta;
  double sqrtfactor = sqrt(factor);
  int Ng;
  gsl_rng *r;
  gsl_interp_accel *acc;
  gsl_spline *spline;

  Ng = conf->Ngrid;
  fac = 2 * M_PI / conf->Lbox;
  norm = pow(conf->Lbox, -1.5);
  acc = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline, Nk);
  gsl_spline_init(spline, lnk, lnP, Nk);

  printf("\r  Generating Fourier modes ... \n");
  fflush(stdout);

  mesh[0][0] = mesh[0][1] = 0;

/*generate the first half of the mesh with random numbers*/
#ifdef OMP
#pragma omp parallel private(r)
{
#endif
  r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, conf->seed + omp_get_thread_num());
#ifdef OMP
#pragma omp for private(j,i,ki,kj,kk,idx,idx_2,i_2,k_2,j_2,ksq,kv,P,fac2, theta)
#endif
  /* Loop in k,j,i order to reduce cache miss. */
  for (k = 0; k <= Ng / 2; k++) {
    kk = (k <= Ng / 2) ? k : k - Ng;
    kk *= fac;
    for (j = 0; j < Ng; j++) {
      kj = (j <= Ng / 2) ? j : j - Ng;
      kj *= fac;
      for (i = 0; i < Ng; i++) {
        ki = (i <= Ng / 2) ? i : i - Ng;
        ki *= fac;
        
        idx = MESH_IDX(Ng, i, j, k);
          
        if (i == 0 && j == 0 && k == 0) {
            mesh[idx][0] = mesh[idx][1] = 0;
            continue;
        }

        ksq = ki * ki + kj * kj + kk * kk;
        kv = log(ksq) * 0.5;    /* log(sqrt(ksq)) */
        P = gsl_spline_eval(spline, kv, acc);
        P = exp(P * 0.5) * sqrtfactor;       /* sqrt(P) */
        
        fac2 = P * norm;
        theta = gsl_ran_flat(r, 0, 2*M_PI);
        if(k == 0 || k == Ng / 2){
          if(j == 0 || j == Ng / 2){
            if(i == 0 || i == Ng / 2){
              //points: j=0,N/2 ; j=0,N/2 ; i=0,N/2
              mesh[idx][0] = 1;
              mesh[idx][1] = 0;
              mesh[idx][0] *= fac2;
            }
            else if(i < Ng / 2){
              //lines: k=0,N/2 ; j=0,N/2 ; 0<i<N/2
              mesh[idx][0] = cos(theta);
              mesh[idx][1] = sin(theta);
              mesh[idx][0] *= fac2;
              mesh[idx][1] *= fac2;

              //lines: k=0,N/2 ; j=0,N/2 ; i>N/2
              i_2 = Ng - i;
              j_2 = j;
              k_2 = k;
              idx_2 = MESH_IDX(Ng, i_2, j_2, k_2);
              mesh[idx_2][0] = mesh[idx][0];
              mesh[idx_2][1] = - mesh[idx][1];
            }
            else{
              //lines: k=0,N/2 ; j=0,N/2 ; i>N/2
              continue;
            }
          }
          else if(j < Ng / 2){
            //planes: k=0,N/2 ; 0<j<N/2 ; i all
            mesh[idx][0] = cos(theta);
            mesh[idx][1] = sin(theta);
            mesh[idx][0] *= fac2;
            mesh[idx][1] *= fac2;

            //planes: k=0,N/2 ; j>N/2 ; i all
            i_2 = (i == 0) ? 0 : Ng - i;
            j_2 = Ng - j;
            k_2 = k;
            idx_2 = MESH_IDX(Ng, i_2, j_2, k_2);
            mesh[idx_2][0] = mesh[idx][0];
            mesh[idx_2][1] = - mesh[idx][1];
          }
          else{
            //planes: k=0,N/2 ; j>N/2 ; i all
            continue;
          } 
        }
        else{
          //volume: 0<k<N/2 ; j all ; i all
          mesh[idx][0] = cos(theta);
          mesh[idx][1] = sin(theta);
          mesh[idx][0] *= fac2;
          mesh[idx][1] *= fac2;

          //volume: k>N/2 ; j all ; i all
          i_2 = (i == 0) ? 0 : Ng - i;
          j_2 = (j == 0) ? 0 : Ng - j;
          k_2 = Ng - k;
          idx_2 = MESH_IDX(Ng, i_2, j_2, k_2);
          mesh[idx_2][0] = mesh[idx][0];
          mesh[idx_2][1] = - mesh[idx][1];
        }
      }
    }
  }
  gsl_rng_free(r);
#ifdef OMP
}
#endif

  printf("\r  Density field in Fourier space is generated.\n");

  /*for (k = 0; k < Ng; k++) {
    printf("k=%d\n", k);
    for (j = Ng - 1; j >= 0; j--) {
      printf("j=%d: ", j);
      for (i = 0; i < Ng; i++) {
        idx = MESH_IDX(Ng, i, j, k);  
        printf("(i=%d : %f + %fi) ", i, mesh[idx][0], mesh[idx][1]);
      }
      printf("\n");
    }
  }*/

  printf("\r  Executing FFT ... \n");
  fflush(stdout);

  FFT_EXEC_C2C(*fp);
  FFT_DESTROY(*fp);
#ifdef OMP
  FFT_CLEAN_OMP();
#endif

  printf("\r  FFT finished successfully.\n");
  fflush(stdout);
  /*for (k = 0; k < Ng; k++) {
    printf("k=%d\n", k);
    for (j = Ng - 1; j >= 0; j--) {
      printf("j=%d: ", j);
      for (i = 0; i < Ng; i++) {
        idx = MESH_IDX(Ng, i, j, k);  
        printf("(i=%d : %f + %fi) ", i, mesh[idx][0], mesh[idx][1]);
      }
      printf("\n");
    }
  }*/

  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  return 0;
}

int save_dens(const char *fname, FFT_CMPLX *mesh, const int Ng) {
  FILE *fp;
  size_t i, Ntot;

  printf("\n  Filename : %s.\n  Preparing for data ... ", fname);
  fflush(stdout);

  Ntot = (size_t) Ng * Ng * Ng;

  if (!(fp = fopen(fname, "w"))) {
    P_ERR("cannot write to file `%s'.\n", fname);
    return ERR_FILE;
  }

  for(i = 0; i<Ntot; i++) {
    fprintf(fp, " %ld %f", i, mesh[i][0]);
  }
  printf("\r  The density field is saved.\n");
  fclose(fp);
  return 0;
}


int save_dens_binary(const char *fname, FFT_CMPLX *mesh, const int Ng)
{
  FILE *fp;
  size_t i, Ntot;

  printf("\n  Filename : %s.\n  Preparing for data ... ", fname);
  fflush(stdout);

  Ntot = (size_t) Ng * Ng * Ng;
#ifdef OMP
#pragma omp parallel for
#endif
  for (i = 0; i < Ntot; i++) mesh[i >> 1][i % 2] = mesh[i][1];

  if (!(fp = fopen(fname, "w"))) {
    P_ERR("cannot write to file `%s'.\n", fname);
    return ERR_FILE;
  }

  if (fwrite(mesh, sizeof(FFT_REAL) * Ntot, 1, fp) != 1) {
    P_EXT("failed to write density field to file `%s'.\n", fname);
    return ERR_FILE;
  }

  printf("\r  The density field is saved.\n");
  fclose(fp);
  return 0;
}
