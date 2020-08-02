#include "linhalo.h"
#include <ctype.h>//TEST
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


/* Quicksort the fist N imaginary elements of mesh. */
void qsort_dens_asc(FFT_CMPLX *mesh, size_t *index_arr, const size_t N, int index) {
  //const int index = 1; //This index tells us if we sort the imaginary(1) part or real(0) one 
  const int M = 7;
  const int NSTACK = 128;
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
  const int NSTACK = 128;
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
    //printf("\n(%d %d %d)", (int)l, (int)m ,(int)r);
    
    // Check if x is present at mid 
    
    if (mesh[index_arr[m]][0] == x) {
      return m; 
    }
    // If x smaller, ignore left half 
    if (mesh[index_arr[m]][0] > x) { 
      l = m + 1;
    }
    // If x is larger, ignore right half 
    else
    {
      if (m == 0) {
        //printf("\nCACACACA");
        return notfound;
      }
      r = m - 1;
      
    }     
  }   

  // if we reach here, then element was 
  // not present
  return notfound; 
} 

size_t binary_search_pos(FFT_CMPLX *mesh, size_t *index_arr, size_t l, size_t r, FFT_REAL x, size_t notfound) {
  size_t m;
  size_t max_index = r;
  size_t min_index = l;
  
  while (l <= r) { 
    m = l + (r - l) / 2;
    // Check if x is present at mid 
    if (mesh[index_arr[m]][0] == x) {
      return m;
    }

    // If x smaller, ignore left half 
    if (mesh[index_arr[m]][0] > x) {
      l = m + 1;
      if(l <= max_index && mesh[index_arr[l]][0] < x) {
        return m;
      }
      
      if(l > max_index) {
        return notfound;
      }
    }
    // If x is larger, ignore right half 
    else {
      r = m - 1;
      if(r >= min_index && mesh[index_arr[r]][0] > x) {
        return r;
      }
      
      if(r < min_index) {
        return notfound;
      }
    }
  } 

  // if we reach here, then element was 
  // not present
  return notfound; 
} 

double return_x(double y) {
  if(y < 0.5) return sqrt(2 * y) - 1;
  else return 1 - sqrt(2 * (1 - y));
}

static void part_cic(FFT_CMPLX *mesh, size_t *index_arr, size_t idx_mesh, FFT_REAL delta_mesh, const size_t Nh, long int *max_index) {
  //printf("\nIN Before sub ;part_cic5: %ld %0.10lf, %0.10lf ", (long int) idx_mesh, mesh[idx_mesh][0], delta_mesh);
  
  //mesh[idx_mesh][0] = mesh[idx_mesh][0] - delta_mesh;
  //printf("\nIN After sub; part_cic5: %ld %0.10lf, %0.10lf ", (long int) idx_mesh, mesh[idx_mesh][0], delta_mesh);
  
  size_t idx_new_pos;
  size_t idx_idx_idx;
  size_t tmp;
  // idx_mesh = 263;
  // Search in the array of length max_index+1 with the max_index+1 largest elements from the mesh, for the neighbour of the maximum element of the mesh
  idx_idx_idx = binary_search(mesh, index_arr, 0, *max_index, mesh[idx_mesh][0], 8 * Nh);
  // printf("\n IN part cic 5: %d %d %lf %lf", (long int) idx_idx_idx, (long int)  index_arr[idx_idx_idx], (double)  mesh[index_arr[idx_idx_idx]][0], mesh[idx_mesh][0]);
  
//   // If the neighbour is in this array then.
  if(idx_idx_idx != 8 * Nh) {
    mesh[idx_mesh][0] = mesh[idx_mesh][0] - delta_mesh;
    // printf("\n IN part cic 5: %d %d %lf %lf", (long int) idx_idx_idx, (long int)  index_arr[idx_idx_idx], (double)  mesh[index_arr[idx_idx_idx]][0], mesh[idx_mesh][0]);
      
    if(mesh[index_arr[idx_idx_idx]][0] < mesh[index_arr[idx_idx_idx + 1]][0]) {
      // Search the new position of the modified element of the mesh.
      idx_new_pos = binary_search_pos(mesh, index_arr, idx_idx_idx + 1, *max_index, mesh[index_arr[idx_idx_idx]][0], 8 * Nh);
      if(idx_new_pos != 8 * Nh) {
        // printf("\n 2nd binary part cic 5: %d %d %lf %lf", (long int) idx_new_pos, (long int)  index_arr[idx_new_pos], (double)  mesh[index_arr[idx_new_pos]][0], mesh[idx_mesh][0]);
        tmp = index_arr[idx_idx_idx];
        memmove(index_arr + idx_idx_idx, index_arr + idx_idx_idx + 1, (idx_new_pos - idx_idx_idx)*sizeof(size_t));
        index_arr[idx_new_pos] = tmp;
      }
      else {
        // printf("\n 3nd binary part cic 5: %d %d %lf %lf", (long int) idx_new_pos, (long int)  index_arr[idx_new_pos], (double)  mesh[index_arr[idx_new_pos]][0], mesh[idx_mesh][0]);
        memmove(index_arr + idx_idx_idx, index_arr + idx_idx_idx + 1, (*max_index - idx_idx_idx)*sizeof(size_t));
      }
    }
  }
  *max_index = *max_index - 1;
  // for(size_t i = 0; i <= *max_index; i++)
  // {
  //   printf("\n%d %lf %lf", (int)i, (double)index_arr[i], (double)mesh[index_arr[i]][0]);//, (double)index_arr[8*Nh-100 + i].index, (double)index_arr[index_arr[8*Nh-100 + i].index] );
  // }

}

static void cic(FFT_CMPLX *mesh, size_t *index_arr, double x, double y, double z, const int Ng, const size_t Nh, const double Lbox, long int *max_index) {
    
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

  // printf("\n%ld %0.10lf ", (long int) MESH_IDX(Ng, x0, y0, z0), mesh[MESH_IDX(Ng, x0, y0, z0)][0]);
  // printf("\n%ld %0.10lf ", (long int) MESH_IDX(Ng, x0, y0, z1), mesh[MESH_IDX(Ng, x0, y0, z1)][0]);
  // printf("\n%ld %0.10lf ", (long int) MESH_IDX(Ng, x0, y1, z0), mesh[MESH_IDX(Ng, x0, y1, z0)][0]);
  // printf("\n%ld %0.10lf ", (long int) MESH_IDX(Ng, x0, y1, z1), mesh[MESH_IDX(Ng, x0, y1, z1)][0]);
  // printf("\n%ld %0.10lf ", (long int) MESH_IDX(Ng, x1, y0, z0), mesh[MESH_IDX(Ng, x1, y0, z0)][0]);
  // printf("\n%ld %0.10lf ", (long int) MESH_IDX(Ng, x1, y0, z1), mesh[MESH_IDX(Ng, x1, y0, z1)][0]);
  // printf("\n%ld %0.10lf ", (long int) MESH_IDX(Ng, x1, y1, z0), mesh[MESH_IDX(Ng, x1, y1, z0)][0]);
  // printf("\n%ld %0.10lf ", (long int) MESH_IDX(Ng, x1, y1, z1), mesh[MESH_IDX(Ng, x1, y1, z1)][0]);
  
  
  //printf("\nBefore:part_cic1: %ld %0.10lf ", (long int) MESH_IDX(Ng, 0, 0, 0), mesh[0][0]);
  //mesh[0][0] = mesh[0][0] - 1;
  //printf("\nAfter:part_cic1: %ld %0.10lf ", (long int) MESH_IDX(Ng, 0, 0, 0), mesh[0][0]);
  
  //printf("\nBefore:part_cic1: %ld %0.10lf %0.10lf", (long int) MESH_IDX(Ng, x0, y0, z0), mesh[MESH_IDX(Ng, x0, y0, z0)][0], wx0 * wy0 * wz0 );
  part_cic(mesh, index_arr, MESH_IDX(Ng, x0, y0, z0), wx0 * wy0 * wz0, Nh, max_index);
  //printf("\naAfter:part_cic1: %ld %0.10lf %0.10lf", (long int) MESH_IDX(Ng, x0, y0, z0), mesh[MESH_IDX(Ng, x0, y0, z0)][0], wx0 * wy0 * wz0);
  
  part_cic(mesh, index_arr, MESH_IDX(Ng, x0, y0, z1), wx0 * wy0 * wz1, Nh, max_index);
  // //printf("\npart_cic2: %ld %0.10lf ", (long int) MESH_IDX(Ng, x0, y0, z1), mesh[MESH_IDX(Ng, x0, y0, z1)][0]);
  
  part_cic(mesh, index_arr, MESH_IDX(Ng, x0, y1, z0), wx0 * wy1 * wz0, Nh, max_index);
  // //printf("\npart_cic3: %ld %0.10lf ", (long int) MESH_IDX(Ng, x0, y1, z0), mesh[MESH_IDX(Ng, x0, y1, z0)][0]);
  //printf("\n\n");
  
  //for (size_t s = 0; s <= *max_index; s++)
  //  printf("%ld %ld %ld %0.10lf %ld \n", (long int) s, (long int) index_arr[s].val, (long int) index_arr[s].rank, mesh[index_arr[s].val][0], (long int) (*max_index));
  
  //printf("\n\n");
  part_cic(mesh, index_arr, MESH_IDX(Ng, x0, y1, z1), wx0 * wy1 * wz1, Nh, max_index);
  // //printf("\npart_cic4: %ld %0.10lf ", (long int) MESH_IDX(Ng, x0, y1, z1), mesh[MESH_IDX(Ng, x0, y1, z1)][0]);
  
  //for (size_t s = 0; s <= *max_index; s++)
  //  printf("%ld %ld %ld %0.10lf %ld \n", (long int) s, (long int) index_arr[s].val, (long int) index_arr[s].rank, mesh[index_arr[s].val][0], (long int) (*max_index));
  
  part_cic(mesh, index_arr, MESH_IDX(Ng, x1, y0, z0), wx1 * wy0 * wz0, Nh, max_index);
  
  //printf("\n\n");
  //for (size_t s = 0; s <= *max_index; s++)
  //  printf("%ld %ld %ld %0.10lf %ld \n", (long int) s, (long int) index_arr[s].val, (long int) index_arr[s].rank, mesh[index_arr[s].val][0], (long int) (*max_index));
  
  //printf("\npart_cic5: %ld %0.10lf, %0.10lf ", (long int) MESH_IDX(Ng, x1, y0, z0), mesh[MESH_IDX(Ng, x1, y0, z0)][0], wx1 * wy0 * wz0);
  
  part_cic(mesh, index_arr, MESH_IDX(Ng, x1, y0, z1), wx1 * wy0 * wz1, Nh, max_index);
  //printf("\npart_cic6: %ld %0.10lf, %0.10lf ", (long int) MESH_IDX(Ng, x1, y0, z1), mesh[MESH_IDX(Ng, x1, y0, z1)][0], wx1 * wy0 * wz1);
  
  part_cic(mesh, index_arr, MESH_IDX(Ng, x1, y1, z0), wx1 * wy1 * wz0, Nh, max_index);
  // //printf("\npart_cic7: %ld %0.10lf ", (long int) MESH_IDX(Ng, x1, y1, z0), mesh[MESH_IDX(Ng, x1, y1, z0)][0]);
  
  part_cic(mesh, index_arr, MESH_IDX(Ng, x1, y1, z1), wx1 * wy1 * wz1, Nh, max_index);
  // //printf("\npart_cic8: %ld %0.10lf ", (long int) MESH_IDX(Ng, x1, y1, z1), mesh[MESH_IDX(Ng, x1, y1, z1)][0]);
  
}

void select_dens(FFT_CMPLX *mesh, HALOS *halos, size_t *index_arr, const int Ng, const size_t Nh, const double Lbox) {
  size_t i, j, k, idx_max, Ntot, idx_val;
  long int max_index;
  double x, y, z;
  size_t tmp;

  Ntot = (size_t) Ng * Ng * Ng;

  gsl_rng *r;
  r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, 49823);  
  
  printf("\r Take the indices of 8Nhalo values from the mesh ... \n");
  fflush(stdout);
#ifdef OMP
#pragma omp parallel for
#endif
  for (i = 0; i < 8 * Nh; i++) index_arr[i] = i;

  
  printf("\r IndexSort the 8Nhalo values from the mesh in ascending order ... \n");
  fflush(stdout);

  qsort_dens_asc(mesh, index_arr, 8 * Nh, 0);
  
  printf("\r Obtain the indices of the 8Nhalo largest values from the mesh ... \n");
  fflush(stdout);
  
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

  printf("\r IndexSort the 8Nhalo LARGEST values from the mesh in descending order... \n");
  fflush(stdout);
  
  qsort_dens_desc(mesh, index_arr, 8 * Nh, 0);



  // for(i = 0; i < 8 * Nh; i++)
  // {
  //   printf("\n%d %lf %0.10lf", (int)i, (double)index_arr[i], mesh[index_arr[i]][0]);//, (double)index_arr[8*Nh-100 + i].index, (double)index_arr[index_arr[8*Nh-100 + i].index] );
  // }
  printf("\nThe maximum density in the mesh is %f ", mesh[index_arr[0]][0]);
  printf("\nThe minimum density in the array of 8Nhalo largest values is %f \n", mesh[index_arr[8*Nh - 1]][0]);
  fflush(stdout);
  //printf("\nPozitie: %lf %lf", (double)binary_search(index_arr, 0, 8 * Nh - 1, MESH_IDX(Ng, 2, 2, 1), 8 * Nh), (double) MESH_IDX(Ng, 2, 2, 1));
  //FILE *ofile;
  //ofile = fopen("./output/test.txt", "w+");
  
  printf("\r Inverse CIC and populate the mesh with halos ... \n Start the time counter");
  fflush(stdout);
  
  max_index = 8 * Nh - 1;
  time_t start = time(NULL);
  time_t end;
  float seconds;
  for (size_t u = 0; u < Nh; u++) {
    
    if(u%10000==0) {
      end = time(NULL);
      seconds = (float)(end - start);
      printf("\n Progress... Halo: %ld/%ld after %f seconds", (long int)u, (long int)Nh, seconds);
      fflush(stdout);
    }

    idx_max = 0;
    idx_val = index_arr[idx_max];
    
    //fprintf(ofile, "######\n");
    //for (size_t s = 0; s < max_index; i++)
    //  fprintf(ofile, "%ld %ld %0.10lf %ld \n", (long int) s, (long int) index_arr[s], mesh[index_arr[s]][0], (long int) (max_index));
    //fprintf(ofile, "\n");
    k = idx_val/(size_t)(Ng * Ng);
    j = (idx_val - (size_t)k * Ng * Ng) / (size_t)Ng;
    i = idx_val - (size_t)k * Ng * Ng - (size_t)j * Ng;
    
    // printf("\nidx_val=%d i=%d j=%d k=%d mesh_idx=%d", idx_val, i ,j ,k, MESH_IDX(Ng, i,j , k));
    x = gsl_ran_flat(r, 0, 1);
    y = gsl_ran_flat(r, 0, 1);
    z = gsl_ran_flat(r, 0, 1);
    
    halos[u].x[0] = Lbox * (i + return_x(x)) / Ng;
    halos[u].x[1] = Lbox * (j + return_x(y)) / Ng;
    halos[u].x[2] = Lbox * (k + return_x(z)) / Ng;
    halos[u].dens = idx_val;
    
    halos[u].x[0] = (halos[u].x[0] < 0) ? (Lbox + halos[u].x[0]) : halos[u].x[0];
    halos[u].x[1] = (halos[u].x[1] < 0) ? (Lbox + halos[u].x[1]) : halos[u].x[1];
    halos[u].x[2] = (halos[u].x[2] < 0) ? (Lbox + halos[u].x[2]) : halos[u].x[2];
    //printf("\n\n %d Before cic: %ld %0.10lf ", (int) u, (long int) index_arr[idx_max], mesh[idx_val][0]);
    cic(mesh, index_arr, halos[u].x[0], halos[u].x[1], halos[u].x[2], Ng, Nh, Lbox, &max_index);
    //printf("\nAfter cic: %ld %0.10lf ", (long int) index_arr[idx_max], mesh[idx_val][0]);
    //max_index = max_index - 8;
  }
  //fclose(ofile);
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
