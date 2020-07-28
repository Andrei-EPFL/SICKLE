#include "linhalo.h"
#include <ctype.h>//TEST
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


/* Quicksort the fist N imaginary elements of mesh. */
#ifdef DOUBLE_PREC
void qsort_dens_asc(fftw_complex *mesh, const size_t N)
#else
void qsort_dens_asc(fftwf_complex *mesh, const size_t N, int index)
#endif
{
  //const int index = 1; //This index tells us if we sort the imaginary(1) part or real(0) one 
  const int M = 7;
  const int NSTACK = 64;
  long istack[NSTACK];
  long i, ir, j, k, jstack, l;
  real a, tmp;

  jstack = -1;
  l = 0;
  ir = N - 1;

  for (;;) {
    if (ir - l < M) {
      for (j = l + 1; j <= ir; j++) {
        a = mesh[j][index];
        for (i = j - 1; i >= l; i--) {
          if (mesh[i][index] <= a) break;
          mesh[i + 1][index] = mesh[i][index];
        }
        mesh[i + 1][index] = a;
      }
      if (jstack < 0) break;
      ir = istack[jstack--];
      l = istack[jstack--];
    }
    else {
      k = (l + ir) >> 1;
      SWAP(mesh[k][index], mesh[l + 1][index], tmp);
      if (mesh[l][index] > mesh[ir][index]) {
        SWAP(mesh[l][index], mesh[ir][index], tmp);
      }
      if (mesh[l + 1][index] > mesh[ir][index]) {
        SWAP(mesh[l + 1][index], mesh[ir][index], tmp);
      }
      if (mesh[l][index] > mesh[l + 1][index]) {
        SWAP(mesh[l][index], mesh[l + 1][index], tmp);
      }
      i = l + 1;
      j = ir;
      a = mesh[l + 1][index];
      for (;;) {
        do i++; while (mesh[i][index] < a);
        do j--; while (mesh[j][index] > a);
        if (j < i) break;
        SWAP(mesh[i][index], mesh[j][index], tmp);
      }
      mesh[l + 1][index] = mesh[j][index];
      mesh[j][index] = a;
      jstack += 2;
      if (jstack >= NSTACK) {
        P_EXT("NSTACK for qsort is too small.\n");
        exit(ERR_OTHER);
      }
      if (ir - i + 1 >= j - 1) {
        istack[jstack] = ir;
        istack[jstack - 1] = i;
        ir = j - 1;
      }
      else {
        istack[jstack] = j - 1;
        istack[jstack - 1] = l;
        l = i;
      }
    }
  }
}

#ifdef DOUBLE_PREC
void qsort_dens_desc(fftw_complex *mesh, const size_t N)
#else
void qsort_dens_desc(fftwf_complex *mesh, const size_t N, int index)
#endif
{
  //const int index = 1; //This index tells us if we sort the imaginary(1) part or real(0) one 
   const int M = 7;
  const int NSTACK = 64;
  long istack[NSTACK];
  long i, ir, j, k, jstack, l;
  real a, tmp;

  jstack = -1;
  l = 0;
  ir = N - 1;

  for (;;) {
    if (ir - l < M) {
      for (j = l + 1; j <= ir; j++) {
        a = mesh[j][index];
        for (i = j - 1; i >= l; i--) {
          if (mesh[i][index] >= a) break;
          mesh[i + 1][index] = mesh[i][index];
        }
        mesh[i + 1][index] = a;
      }
      if (jstack < 0) break;
      ir = istack[jstack--];
      l = istack[jstack--];
    }
    else {
      k = (l + ir) >> 1;
      SWAP(mesh[k][index], mesh[l + 1][index], tmp);
      if (mesh[l][index] < mesh[ir][index]) {
        SWAP(mesh[l][index], mesh[ir][index], tmp);
      }
      if (mesh[l + 1][index] < mesh[ir][index]) {
        SWAP(mesh[l + 1][index], mesh[ir][index], tmp);
      }
      if (mesh[l][index] < mesh[l + 1][index]) {
        SWAP(mesh[l][index], mesh[l + 1][index], tmp);
      }
      i = l + 1;
      j = ir;
      a = mesh[l + 1][index];
      for (;;) {
        do i++; while (mesh[i][index] > a);
        do j--; while (mesh[j][index] < a);
        if (j < i) break;
        SWAP(mesh[i][index], mesh[j][index], tmp);
      }
      mesh[l + 1][index] = mesh[j][index];
      mesh[j][index] = a;
      jstack += 2;
      if (jstack >= NSTACK) {
        P_EXT("NSTACK for qsort is too small.\n");
        exit(ERR_OTHER);
      }
      if (ir - i + 1 >= j - 1) {
        istack[jstack] = ir;
        istack[jstack - 1] = i;
        ir = j - 1;
      }
      else {
        istack[jstack] = j - 1;
        istack[jstack - 1] = l;
        l = i;
      }
    }
  }
}


double return_x(double y) {
  if(y < 0.5) return sqrt(2 * y) - 1;
  else return 1 - sqrt(2 * (1 - y));
}


// void function_test() {
//   gsl_rng *r;
//   r = gsl_rng_alloc(gsl_rng_mt19937);
//   gsl_rng_set(r, 49823);
//   for(int i=0; i<10; i++) {
//     double y = gsl_ran_flat(r, 0, 1);
//     printf("%lf\n", return_x(y));
//   }
  
//   gsl_rng_free(r);

// }

/******************************************************************************
Function `cic`:
  Assign particles to the mesh with the Could-In-Cell scheme.
Arguments:
  * `data`:     structure for the data sample;
  * `ndata`:    number of objects in the sample;
  * `Ng`:       number of grid cells per box side;
  * `bmin`:     lower boundaries of the box;
  * `Lbox`:     side length of the box;
  * `rho`:      the mesh.
******************************************************************************/

#ifdef DOUBLE_PREC
static void cic(fftw_complex *mesh, double x, double y, double z, const int Ng,
    const double Lbox)
#else
static void cic(fftwf_complex *mesh, double x, double y, double z, const int Ng,
    const double Lbox)
#endif
{
  size_t Ntot;
  Ntot = (size_t) Ng * Ng * Ng;

  double mesh_x = x * Ng / Lbox;
  double mesh_y = y * Ng / Lbox;
  double mesh_z = z * Ng / Lbox;

  int x0 = (int) mesh_x;
  int y0 = (int) mesh_y;
  int z0 = (int) mesh_z;

  printf("\nidx:%f %f %f\n", x,y,z);
  
  /* Weights for neighbours. */
  double wx1 = mesh_x - x0;
  double wx0 = 1 - wx1;
  double wy1 = mesh_y - y0;
  double wy0 = 1 - wy1;
  double wz1 = mesh_z - z0;
  double wz0 = 1 - wz1;

  int x1 = (x0 == Ng - 1) ? 0 : x0 + 1;
  int y1 = (y0 == Ng - 1) ? 0 : y0 + 1;
  int z1 = (z0 == Ng - 1) ? 0 : z0 + 1;

  //mesh[MESH_IDX(Ng,x0,y0,z0)][0] -= wx0 * wy0 * wz0;
  mesh[MESH_IDX(Ng,x0,y0,z1)][0] -= wx0 * wy0 * wz1;
  printf(" \nVAL:%d \n", MESH_IDX(Ng,x0,y0,z1));
  
  // mesh[MESH_IDX(Ng,x0,y1,z0)][0] -= wx0 * wy1 * wz0;
  // mesh[MESH_IDX(Ng,x0,y1,z1)][0] -= wx0 * wy1 * wz1;
  // mesh[MESH_IDX(Ng,x1,y0,z0)][0] -= wx1 * wy0 * wz0;
  // mesh[MESH_IDX(Ng,x1,y0,z1)][0] -= wx1 * wy0 * wz1;
  // mesh[MESH_IDX(Ng,x1,y1,z0)][0] -= wx1 * wy1 * wz0;
  // mesh[MESH_IDX(Ng,x1,y1,z1)][0] -= wx1 * wy1 * wz1;
  //insertion_sort(mesh, Ntot, (size_t)MESH_IDX(Ng,x0,y0,z0));
  
}


#ifdef DOUBLE_PREC
void select_dens(fftw_complex *mesh, HALOS *halos, const int Ng,
    const size_t Nh, const double Lbox)
#else
void select_dens(fftwf_complex *mesh, HALOS *halos, const int Ng,
    const size_t Nh, const double Lbox)
#endif
{
  size_t i, j, k, idx, Ntot, idx_val;
  double x, y, z;
  real tmp;

  Ntot = (size_t) Ng * Ng * Ng;

  gsl_rng *r;
  r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, 49823);  

#ifdef OMP
#pragma omp parallel for
#endif
  for (i = 0; i < 8 * Nh; i++) mesh[i][1] = mesh[i][0];
  
  qsort_dens_asc(mesh, 8 * Nh, 1);
  
  for (i = 8 * Nh; i < Ntot; i++) {
    if (mesh[i][0] > mesh[0][1]) {
      mesh[0][1] = mesh[i][0];
      for (j = 0;;) {
        k = (j << 1) + 1;
        if (k > 8 * Nh - 1) break;
        if (k != 8 * Nh - 1 && mesh[k][1] > mesh[k + 1][1]) ++k;
        if (mesh[j][1] <= mesh[k][1]) break;
        SWAP(mesh[k][1], mesh[j][1], tmp);
        j = k;
      }
    }
  }
  qsort_dens_desc(mesh, 8 * Nh, 1);
  qsort_dens_desc(mesh, Ntot, 0);
  
  
  // for (int u = 0; u < Nh; u++) {

  //   idx = 0;
  //   idx_val = mesh[idx][1];
  //   k = idx_val/(Ng * Ng);
  //   j = (idx_val - k * Ng * Ng) / Ng;
  //   i = idx_val - k * Ng * Ng - j * Ng;
  //   printf("\n%d %d %f\n", idx_val, MESH_IDX(Ng, i, j, k), mesh[idx_val][0]);
  //   x = gsl_ran_flat(r, 0, 1);
  //   y = gsl_ran_flat(r, 0, 1);
  //   z = gsl_ran_flat(r, 0, 1);
  //   halos[u].x[0] = Lbox * (i + return_x(x)) / Ng;
  //   halos[u].x[1] = Lbox * (j + return_x(y)) / Ng;
  //   halos[u].x[2] = Lbox * (k + return_x(z)) / Ng;
    
  //   //cic(mesh, halos[u].x[0], halos[u].x[1], halos[u].x[2], Ng, Lbox);
  // }
  
  gsl_rng_free(r);
  
}


int save_halo(const char *fname, HALOS *halos, const size_t Nh) {
  FILE *fp;
  int i, m, n;
  char *buf, *end, *cache;

  printf("\n  Filename : %s.\n", fname);
  if (!(fp = fopen(fname, "w"))) {
    P_ERR("cannot write to file `%s'.\n", fname);
    return ERR_FILE;
  }
  
#ifdef OMP
#pragma omp parallel private(buf,end,cache,n) shared(fp)
{
#endif
  buf = calloc(CHUNK, sizeof *buf);
  cache = calloc(MAX_LEN_LINE, sizeof *cache);
  if (!buf || !cache) {
    P_EXT("failed to allocate memory for writing outputs.\n");
    exit(ERR_MEM);
  }
  end = buf;

#ifdef OMP
#pragma omp for private(m)
#endif
  for (i = 0; i < Nh; i++) {
    n = snprintf(cache, MAX_LEN_LINE,
        OFMT_REAL " " OFMT_REAL " " OFMT_REAL "\n", halos[i].x[0], halos[i].x[0], halos[i].x[0]);
    if (n < 0 || n >= MAX_LEN_LINE) {
      P_EXT(FMT_KEY(MAX_LEN_LINE)
          " in `define.h' is not large enough.\n");
      exit(ERR_STRING);
    }

    if (end - buf + n < CHUNK) {        /* there is still space in buf */
      m = cnt_strcpy(end, cache, n + 1);
      if (m >= n + 1) {
        P_EXT("unexpected error for writing line:\n%s\n", cache);
        exit(ERR_STRING);
      }
      end += m;
    }
    else {                              /* write buf to file */
#ifdef OMP
#pragma omp critical
      {
#endif
      if (fwrite(buf, sizeof(char) * (end - buf), 1, fp) != 1) {
        P_EXT("failed to write to output:\n%s\n", cache);
        exit(ERR_FILE);
      }
      fflush(fp);
#ifdef OMP
      }
#endif
      m = cnt_strcpy(buf, cache, n + 1);
      if (m >= n + 1) {
        P_EXT("unexpected error for writing line:\n%s\n", cache);
        exit(ERR_STRING);
      }
      end = buf + m;
    }    
  }

  if ((n = end - buf) > 0) {
#ifdef OMP
#pragma omp critical
    {
#endif
    if (fwrite(buf, sizeof(char) * n, 1, fp) != 1) {
      P_EXT("failed to write to output:\n%s\n", cache);
      exit(ERR_FILE);
    }
    fflush(fp);
#ifdef OMP
    }
#endif
  }

  free(buf);
  free(cache);

#ifdef OMP
}
#endif

  fclose(fp);
  return 0;
}

size_t cnt_strcpy(char *dest, const char *src, const size_t num) {
  size_t i = 0;
  while (i < num - 1 && src[i] != '\0') {
    dest[i] = src[i];
    i++;
  }
  dest[i] = '\0';
  while (src[i] != '\0') i++;
  return i;
}
