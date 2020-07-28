#include "linhalo.h"
#include <ctype.h>//TEST
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


/* Quicksort the fist N imaginary elements of mesh. */
void qsort_dens_asc(FFT_CMPLX *mesh, size_t *index_arr, const size_t N, int index) {
  //const int index = 1; //This index tells us if we sort the imaginary(1) part or real(0) one 
  const int M = 7;
  const int NSTACK = 64;
  long istack[NSTACK];
  long i, ir, j, k, jstack, l;
  size_t a, tmp;

  jstack = -1;
  l = 0;
  ir = N - 1;

  for (;;) {
    if (ir - l < M) {
      for (j = l + 1; j <= ir; j++) {
        a = index_arr[j];
        for (i = j - 1; i >= l; i--) {
          if (mesh[index_arr[i]][index] <= mesh[a][index]) break;
          index_arr[i + 1] = index_arr[i];
        }
        index_arr[i + 1] = a;
      }
      if (jstack < 0) break;
      ir = istack[jstack--];
      l = istack[jstack--];
    }
    else {
      k = (l + ir) >> 1;
      SWAP(index_arr[k], index_arr[l + 1], tmp);
      if (mesh[index_arr[l]][index] > mesh[index_arr[ir]][index]) {
        SWAP(index_arr[l], index_arr[ir], tmp);
      }
      if (mesh[index_arr[l + 1]][index] > mesh[index_arr[ir]][index]) {
        SWAP(index_arr[l + 1], index_arr[ir], tmp);
      }
      if (mesh[index_arr[l]][index] > mesh[index_arr[l + 1]][index]) {
        SWAP(index_arr[l], index_arr[l + 1], tmp);
      }
      i = l + 1;
      j = ir;
      a = index_arr[l + 1];
      for (;;) {
        do i++; while (mesh[index_arr[i]][index] < mesh[a][index]);
        do j--; while (mesh[index_arr[j]][index] > mesh[a][index]);
        if (j < i) break;
        SWAP(index_arr[i], index_arr[j], tmp);
      }
      index_arr[l + 1] = index_arr[j];
      index_arr[j] = a;
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

void qsort_dens_desc(FFT_CMPLX *mesh, size_t *index_arr, const size_t N, int index) {
  //const int index = 1; //This index tells us if we sort the imaginary(1) part or real(0) one 
   const int M = 7;
  const int NSTACK = 64;
  long istack[NSTACK];
  long i, ir, j, k, jstack, l;
  size_t a, tmp;

  jstack = -1;
  l = 0;
  ir = N - 1;

  for (;;) {
    if (ir - l < M) {
      for (j = l + 1; j <= ir; j++) {
        a = index_arr[j];
        for (i = j - 1; i >= l; i--) {
          if (mesh[index_arr[i]][index] >= mesh[a][index]) break;
          index_arr[i + 1] = index_arr[i];
        }
        index_arr[i + 1] = a;
      }
      if (jstack < 0) break;
      ir = istack[jstack--];
      l = istack[jstack--];
    }
    else {
      k = (l + ir) >> 1;
      SWAP(index_arr[k], index_arr[l + 1], tmp);
      if (mesh[index_arr[l]][index] < mesh[index_arr[ir]][index]) {
        SWAP(index_arr[l], index_arr[ir], tmp);
      }
      if (mesh[index_arr[l + 1]][index] < mesh[index_arr[ir]][index]) {
        SWAP(index_arr[l + 1], index_arr[ir], tmp);
      }
      if (mesh[index_arr[l]][index] < mesh[index_arr[l + 1]][index]) {
        SWAP(index_arr[l], index_arr[l + 1], tmp);
      }
      i = l + 1;
      j = ir;
      a = index_arr[l + 1];
      for (;;) {
        do i++; while (mesh[index_arr[i]][index] > mesh[a][index]);
        do j--; while (mesh[index_arr[j]][index] < mesh[a][index]);
        if (j < i) break;
        SWAP(index_arr[i], index_arr[j], tmp);
      }
      index_arr[l + 1] = index_arr[j];
      index_arr[j] = a;
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

size_t binary_search(FFT_CMPLX *mesh, size_t *index_arr, size_t l, size_t r, FFT_REAL x, size_t notfound) {
  size_t m;
  
  while (l <= r) { 
    m = l + (r - l) / 2; 
    // Check if x is present at mid 
    if (mesh[index_arr[m]][0] == x) 
      return m; 

    // If x smaller, ignore left half 
    if (mesh[index_arr[m]][0] > x) 
      l = m + 1; 

    // If x is larger, ignore right half 
    else
      r = m - 1; 
  } 

  // if we reach here, then element was 
  // not present
  return notfound; 
} 

double return_x(double y) {
  if(y < 0.5) return sqrt(2 * y) - 1;
  else return 1 - sqrt(2 * (1 - y));
}

static void cic(FFT_CMPLX *mesh, size_t *index_arr, double x, double y, double z, const int Ng, const size_t Nh, const double Lbox) {
  //size_t Ntot;
  //Ntot = (size_t) Ng * Ng * Ng;
    
  FFT_REAL mesh_x = x * Ng / Lbox;
  FFT_REAL mesh_y = y * Ng / Lbox;
  FFT_REAL mesh_z = z * Ng / Lbox;

  int x0 = (int) mesh_x;
  int y0 = (int) mesh_y;
  int z0 = (int) mesh_z;
  printf("\n %d %d %d", (int) x0, (int) y0, (int) z0);
  printf("\n %f %f %f", mesh_x, mesh_y, mesh_z);
  //printf("\nidx:%f %f %f\n", x,y,z);
  
  /* Weights for neighbours. */
  FFT_REAL wx1 = mesh_x - x0;
  FFT_REAL wx0 = 1 - wx1;
  FFT_REAL wy1 = mesh_y - y0;
  FFT_REAL wy0 = 1 - wy1;
  FFT_REAL wz1 = mesh_z - z0;
  FFT_REAL wz0 = 1 - wz1;

  int x1 = (x0 == Ng - 1) ? 0 : x0 + 1;
  int y1 = (y0 == Ng - 1) ? 0 : y0 + 1;
  int z1 = (z0 == Ng - 1) ? 0 : z0 + 1;

  size_t idx1 = binary_search(mesh, index_arr, (size_t)0, 8 * Nh - 1, mesh[MESH_IDX(Ng,x0,y0,z0)][0], 8 * Nh);
  if(idx1 != 8 * Nh) {
    mesh[MESH_IDX(Ng,x0,y0,z0)][0] -= wx0 * wy0 * wz0;
    
    if(mesh[index_arr[idx1]][0] < )
  }
  

  // size_t idx2 = binary_search(mesh, index_arr, (size_t)0, 8 * Nh - 1, mesh[MESH_IDX(Ng,x0,y0,z1)][0]);
   
  // size_t idx3 = binary_search(mesh, index_arr, (size_t)0, 8 * Nh - 1, mesh[MESH_IDX(Ng,x0,y1,z0)][0]);
   
  // size_t idx4 = binary_search(mesh, index_arr, (size_t)0, 8 * Nh - 1, mesh[MESH_IDX(Ng,x0,y1,z1)][0]);
   

  // size_t idx5 = binary_search(mesh, index_arr, (size_t)0, 8 * Nh - 1, mesh[MESH_IDX(Ng,x1,y0,z0)][0]);
   
  // size_t idx6 = binary_search(mesh, index_arr, (size_t)0, 8 * Nh - 1, mesh[MESH_IDX(Ng,x1,y0,z1)][0]);
   
  // size_t idx7 = binary_search(mesh, index_arr, (size_t)0, 8 * Nh - 1, mesh[MESH_IDX(Ng,x1,y1,z0)][0]);
   
  // size_t idx8 = binary_search(mesh, index_arr, (size_t)0, 8 * Nh - 1, mesh[MESH_IDX(Ng,x1,y1,z1)][0]);
  // mesh[MESH_IDX(Ng,x0,y0,z1)][0] -= wx0 * wy0 * wz1;
  // printf("\n%f", mesh[index_arr[0]][0]);
  // printf("\nId: %d Value: %f and Position %d;", (int)MESH_IDX(Ng,x0,y0,z1), mesh[MESH_IDX(Ng,x0,y0,z1)][0], (int)idx2);
  
  // mesh[MESH_IDX(Ng,x0,y1,z0)][0] -= wx0 * wy1 * wz0; 
  // printf("\n%f", mesh[index_arr[0]][0]);
  // printf("\nId: %d Value: %f and Position %d;", (int)MESH_IDX(Ng,x0,y1,z0), mesh[MESH_IDX(Ng,x0,y1,z0)][0], (int)idx3);
  
  // mesh[MESH_IDX(Ng,x0,y1,z1)][0] -= wx0 * wy1 * wz1;
  // printf("\n%f", mesh[index_arr[0]][0]);
  // printf("\nId: %d Value: %f and Position %d;", (int)MESH_IDX(Ng,x0,y1,z1), mesh[MESH_IDX(Ng,x0,y1,z1)][0], (int)idx4);
  
  // mesh[MESH_IDX(Ng,x1,y0,z0)][0] -= wx1 * wy0 * wz0;
  // printf("\n%f", mesh[index_arr[0]][0]);
  // printf("\nId: %d Value: %f and Position %d;", (int)MESH_IDX(Ng,x1,y0,z0), mesh[MESH_IDX(Ng,x1,y0,z0)][0], (int)idx5);

  // mesh[MESH_IDX(Ng,x1,y0,z1)][0] -= wx1 * wy0 * wz1;
  // printf("\n%f", mesh[index_arr[0]][0]);
  // printf("\nId: %d Value: %f and Position %d;", (int)MESH_IDX(Ng,x1,y0,z1), mesh[MESH_IDX(Ng,x1,y0,z1)][0], (int)idx6);

  // mesh[MESH_IDX(Ng,x1,y1,z0)][0] -= wx1 * wy1 * wz0;
  // printf("\n%f", mesh[index_arr[0]][0]);
  // printf("\nId: %d Value: %f and Position %d;", (int)MESH_IDX(Ng,x1,y1,z0), mesh[MESH_IDX(Ng,x1,y1,z0)][0], (int)idx7);

  // mesh[MESH_IDX(Ng,x1,y1,z1)][0] -= wx1 * wy1 * wz1;
  // printf("\n%f", mesh[index_arr[0]][0]);
  // printf("\nId: %d Value: %f and Position %d;", (int)MESH_IDX(Ng,x1,y1,z1), mesh[MESH_IDX(Ng,x1,y1,z1)][0], (int)idx8);
}



void select_dens(FFT_CMPLX *mesh, HALOS *halos, size_t *index_arr, const int Ng,
    const size_t Nh, const double Lbox) {
  size_t i, j, k, idx_max, Ntot, idx_val;
  double x, y, z;
  FFT_REAL tmp;

  Ntot = (size_t) Ng * Ng * Ng;

  gsl_rng *r;
  r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, 49823);  

#ifdef OMP
#pragma omp parallel for
#endif
  for (i = 0; i < 8 * Nh; i++) index_arr[i] = i;
  
  qsort_dens_asc(mesh, index_arr, 8 * Nh, 0);
    
  for (i = 8 * Nh; i < Ntot; i++) {
    if (mesh[i][0] > mesh[index_arr[0]][0]) {
      index_arr[0] = i;
      for (j = 0;;) {
        k = (j << 1) + 1;
        if (k > 8 * Nh - 1) break;
        if (k != 8 * Nh - 1 && mesh[index_arr[k]][0] > mesh[index_arr[k + 1]][0]) ++k;
        if (mesh[index_arr[j]][0] <= mesh[index_arr[k]][0]) break;
        SWAP(index_arr[k], index_arr[j], tmp);
        j = k;
      }
    }
  }
  qsort_dens_desc(mesh, index_arr, 8 * Nh, 0);

  for(i = 0; i < 8*Nh; i++)
  {
    printf("\n%d %f %d %f", (int)i, mesh[i][0], (int)index_arr[i], mesh[index_arr[i]][0] );
  }

  
  for (size_t u = 0; u < 1; u++) {
    printf("\n");
    idx_max = 0;
    idx_val = index_arr[idx_max];
    k = idx_val/(size_t)(Ng * Ng);
    j = (idx_val - (size_t)k * Ng * Ng) / (size_t)Ng;
    i = idx_val - (size_t)k * Ng * Ng - (size_t)j * Ng;
    
    x = gsl_ran_flat(r, 0, 1);
    y = gsl_ran_flat(r, 0, 1);
    z = gsl_ran_flat(r, 0, 1);
    
    halos[u].x[0] = Lbox * (i + return_x(x)) / Ng;
    halos[u].x[1] = Lbox * (j + return_x(y)) / Ng;
    halos[u].x[2] = Lbox * (k + return_x(z)) / Ng;
    
    cic(mesh, index_arr, halos[u].x[0], halos[u].x[1], halos[u].x[2], Ng, Nh, Lbox);
  }
  
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
