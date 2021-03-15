#include "linhalo.h"
#include <ctype.h>//TEST
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* Quicksort the fist N imaginary elements of mesh. */
void qsort_dens_asc(MAX_DENS *max_dens, const size_t N) {
  //const int index = 1; //This index tells us if we sort the imaginary(1) part or real(0) one 
  const int M = 7;
  const int NSTACK = 128;
  long istack[NSTACK];
  long i, ir, j, k, jstack, l;
  MAX_DENS a, tmp;

  jstack = -1;
  l = 0;
  ir = N - 1;

  for (;;) {
    if (ir - l < M) {
      for (j = l + 1; j <= ir; j++) {
        a = max_dens[j];
        for (i = j - 1; i >= l; i--) {
          if (max_dens[i].dens <= a.dens) break;
          max_dens[i + 1] = max_dens[i];
        }
        max_dens[i + 1] = a;
      }
      if (jstack < 0) break;
      ir = istack[jstack--];
      l = istack[jstack--];
    }
    else {
      k = (l + ir) >> 1;
      SWAP(max_dens[k], max_dens[l + 1], tmp);
      if (max_dens[l].dens > max_dens[ir].dens) {
        SWAP(max_dens[l], max_dens[ir], tmp);
      }
      if (max_dens[l + 1].dens > max_dens[ir].dens) {
        SWAP(max_dens[l + 1], max_dens[ir], tmp);
      }
      if (max_dens[l].dens > max_dens[l + 1].dens) {
        SWAP(max_dens[l], max_dens[l + 1], tmp);
      }
      i = l + 1;
      j = ir;
      a = max_dens[l + 1];
      for (;;) {
        do i++; while (max_dens[i].dens < a.dens);
        do j--; while (max_dens[j].dens > a.dens);
        if (j < i) break;
        SWAP(max_dens[i], max_dens[j], tmp);
      }
      max_dens[l + 1] = max_dens[j];
      max_dens[j] = a;
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

static void max_heapify(MAX_DENS *max_dens, size_t i, const size_t n, size_t *index_arr) 
{
  size_t l, r, largest;
  MAX_DENS tmp;
  
  while (i < n) {
    l = 2*i + 1;
    r = 2*i + 2;
    largest = i;
    index_arr[max_dens[i].idx] = i;
    
    if (l < n) {
      index_arr[max_dens[l].idx] = l;
      if (max_dens[l].dens > max_dens[largest].dens) largest = l;
    }

    if (r < n) {
      index_arr[max_dens[r].idx] = r;
      if (max_dens[r].dens > max_dens[largest].dens) largest = r; 
    }
    
    if (largest != i) {
      SWAP(max_dens[i], max_dens[largest], tmp);
      index_arr[max_dens[i].idx] = i;
      index_arr[max_dens[largest].idx] = largest;
      i = largest;
    }
    else {
      break;
    }
  }
} 

double return_x(double y) {
  if(y < 0.5) return sqrt(2 * y) - 1;
  else return 1 - sqrt(2 * (1 - y));
}

static void part_cic(MAX_DENS *max_dens, size_t *index_arr, size_t idx_mesh, FFT_REAL delta_mesh, const size_t Nh, const size_t Ntot) {
  size_t idx_max_dens = index_arr[idx_mesh];
  if (idx_max_dens != Ntot + 1) {
    max_dens[idx_max_dens].dens -= delta_mesh;
    max_heapify(max_dens, idx_max_dens, 8*Nh, index_arr);
  }
}

static void cic(MAX_DENS *max_dens, size_t *index_arr, double x, double y, double z, const int Ng, const size_t Nh, const double Lbox) {
  
  size_t Ntot = (size_t) Ng * Ng * Ng;
  FFT_REAL mesh_x = x * Ng / Lbox;
  FFT_REAL mesh_y = y * Ng / Lbox;
  FFT_REAL mesh_z = z * Ng / Lbox;
  
  int x0 = (int) mesh_x;
  int y0 = (int) mesh_y;
  int z0 = (int) mesh_z;
  
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
  
  part_cic(max_dens, index_arr, MESH_IDX(Ng, x0, y0, z0), wx0 * wy0 * wz0, Nh, Ntot);
  part_cic(max_dens, index_arr, MESH_IDX(Ng, x0, y0, z1), wx0 * wy0 * wz1, Nh, Ntot);
  part_cic(max_dens, index_arr, MESH_IDX(Ng, x0, y1, z0), wx0 * wy1 * wz0, Nh, Ntot);
  part_cic(max_dens, index_arr, MESH_IDX(Ng, x0, y1, z1), wx0 * wy1 * wz1, Nh, Ntot);
  part_cic(max_dens, index_arr, MESH_IDX(Ng, x1, y0, z0), wx1 * wy0 * wz0, Nh, Ntot);
  part_cic(max_dens, index_arr, MESH_IDX(Ng, x1, y0, z1), wx1 * wy0 * wz1, Nh, Ntot);
  part_cic(max_dens, index_arr, MESH_IDX(Ng, x1, y1, z0), wx1 * wy1 * wz0, Nh, Ntot);
  part_cic(max_dens, index_arr, MESH_IDX(Ng, x1, y1, z1), wx1 * wy1 * wz1, Nh, Ntot);
}

void select_dens(FFT_CMPLX *mesh, HALOS *halos, MAX_DENS *max_dens, const int Ng, const size_t Nh, const double Lbox) {
  size_t i, j, k, u, Ntot;
  size_t idx_max, idx_mesh;
  size_t n = 8*Nh;
  
  double x, y, z;
  MAX_DENS tmp;

  Ntot = (size_t) Ng * Ng * Ng;

  gsl_rng *r;
  r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, 49823);  
  
  printf("\r Take the 8Nhalo values and their indices from the mesh ... \n");
  fflush(stdout);
#ifdef OMP
#pragma omp parallel for
#endif
  for (i = 0; i < 8 * Nh; i++) {
    max_dens[i].idx = i;
    max_dens[i].dens = mesh[i][0];
  }
  
  printf("\r Sort the 8Nhalo values from the max_dens in ascending order ... \n");
  fflush(stdout);

  qsort_dens_asc(max_dens, 8 * Nh);
  
  printf("\r Obtain the 8Nhalo largest values and their indices from the mesh ... \n");
  fflush(stdout);
  
  for (i = 8 * Nh; i < Ntot; i++) {
    if (mesh[i][0] > max_dens[0].dens) {
      max_dens[0].dens = mesh[i][0];
      max_dens[0].idx = i;
      for (j = 0;;) {
        k = (j << 1) + 1;
        if (k > 8 * Nh - 1) break;
        if (k != 8 * Nh - 1 && max_dens[k].dens > max_dens[k + 1].dens) ++k;
        if (max_dens[j].dens <= max_dens[k].dens) break;
        SWAP(max_dens[k], max_dens[j], tmp);
        j = k;
      }
    }
  }
  
  printf("\nThe minimum density in the array of 8Nhalo largest values is %f \n", max_dens[0].dens);
  
  printf("\r Cast the memory allocated for the mesh to an array of indices and initialize it with SIZE_MAX ... \n");
  fflush(stdout);
  
  size_t *index_arr = (size_t *) mesh;

#ifdef OMP
#pragma omp parallel for
#endif
  for (i = 0; i < Ntot; i++) {
    index_arr[i] = Ntot + 1;
  }
  
  printf("\r Create a Max-Heap with 8Nhalo LARGEST values from the max_dens... \n");
  fflush(stdout);
  
  // This function basically builds max heap from min heap 
  // Start from bottommost and rightmost 
  // internal mode and heapify all internal 
  // modes in bottom up way 

  
  for (u = (n-2)/2; u >= 0; --u) {
    max_heapify(max_dens, u, n, index_arr);
    if (u==0)
      break;
  }

  printf("The maximum density in the mesh is %f \n", max_dens[0].dens);
  fflush(stdout);
  
  printf("\r Inverse CIC and populate the mesh with halos ... \n Start the time counter");
  fflush(stdout);
  
  time_t start = time(NULL);
  time_t end;
  float seconds;
  for (u = 0; u < Nh; u++) {
    
    idx_max = 0;
    idx_mesh = max_dens[idx_max].idx;
    
    k = idx_mesh/(size_t)(Ng * Ng);
    j = (idx_mesh - (size_t)k * Ng * Ng) / (size_t)Ng;
    i = idx_mesh - (size_t)k * Ng * Ng - (size_t)j * Ng;
    
    x = gsl_ran_flat(r, 0, 1);
    y = gsl_ran_flat(r, 0, 1);
    z = gsl_ran_flat(r, 0, 1);
    
    halos[u].x[0] = Lbox * (i + return_x(x)) / Ng;
    halos[u].x[1] = Lbox * (j + return_x(y)) / Ng;
    halos[u].x[2] = Lbox * (k + return_x(z)) / Ng;
    halos[u].dens = idx_mesh;
    
    halos[u].x[0] = (halos[u].x[0] < 0) ? (Lbox + halos[u].x[0]) : halos[u].x[0];
    halos[u].x[1] = (halos[u].x[1] < 0) ? (Lbox + halos[u].x[1]) : halos[u].x[1];
    halos[u].x[2] = (halos[u].x[2] < 0) ? (Lbox + halos[u].x[2]) : halos[u].x[2];
    cic(max_dens, index_arr, halos[u].x[0], halos[u].x[1], halos[u].x[2], Ng, Nh, Lbox);
    
    if(u%500000==0 || u == Nh - 1) {
      end = time(NULL);
      seconds = (float)(end - start);
      printf("\n Progress... Halo: %ld/%ld after %f seconds", (long int)u, (long int)Nh, seconds);
      fflush(stdout);
    }
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
        OFMT_REAL " " OFMT_REAL " " OFMT_REAL " " OFMT_REAL "\n", halos[i].x[0], halos[i].x[1], halos[i].x[2], (double)halos[i].dens);
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
