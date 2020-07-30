#include "linhalo.h"
#include <ctype.h>//TEST
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


/* Quicksort the fist N imaginary elements of mesh. */
void qsort_dens_asc(FFT_CMPLX *mesh, INDEX *index_arr, const size_t N, int index) {
  //const int index = 1; //This index tells us if we sort the imaginary(1) part or real(0) one 
  const int M = 7;
  const int NSTACK = 128;
  long istack[NSTACK];
  long i, ir, j, k, jstack, l;
  INDEX a, tmp;

  jstack = -1;
  l = 0;
  ir = N - 1;

  for (;;) {
    if (ir - l < M) {
      for (j = l + 1; j <= ir; j++) {
        a = index_arr[j];
        for (i = j - 1; i >= l; i--) {
          if (mesh[index_arr[i].val][index] <= mesh[a.val][index]) break;
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
      if (mesh[index_arr[l].val][index] > mesh[index_arr[ir].val][index]) {
        SWAP(index_arr[l], index_arr[ir], tmp);
      }
      if (mesh[index_arr[l + 1].val][index] > mesh[index_arr[ir].val][index]) {
        SWAP(index_arr[l + 1], index_arr[ir], tmp);
      }
      if (mesh[index_arr[l].val][index] > mesh[index_arr[l + 1].val][index]) {
        SWAP(index_arr[l], index_arr[l + 1], tmp);
      }
      i = l + 1;
      j = ir;
      a = index_arr[l + 1];
      for (;;) {
        do i++; while (mesh[index_arr[i].val][index] < mesh[a.val][index]);
        do j--; while (mesh[index_arr[j].val][index] > mesh[a.val][index]);
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

void qsort_dens_desc(FFT_CMPLX *mesh, INDEX *index_arr, const size_t N, int index) {
  //const int index = 1; //This index tells us if we sort the imaginary(1) part or real(0) one 
   const int M = 7;
  const int NSTACK = 128;
  long istack[NSTACK];
  long i, ir, j, k, jstack, l;
  INDEX a, tmp;

  jstack = -1;
  l = 0;
  ir = N - 1;

  for (;;) {
    if (ir - l < M) {
      for (j = l + 1; j <= ir; j++) {
        a = index_arr[j];
        for (i = j - 1; i >= l; i--) {
          if (mesh[index_arr[i].val][index] >= mesh[a.val][index]) break;
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
      if (mesh[index_arr[l].val][index] < mesh[index_arr[ir].val][index]) {
        SWAP(index_arr[l], index_arr[ir], tmp);
      }
      if (mesh[index_arr[l + 1].val][index] < mesh[index_arr[ir].val][index]) {
        SWAP(index_arr[l + 1], index_arr[ir], tmp);
      }
      if (mesh[index_arr[l].val][index] < mesh[index_arr[l + 1].val][index]) {
        SWAP(index_arr[l], index_arr[l + 1], tmp);
      }
      i = l + 1;
      j = ir;
      a = index_arr[l + 1];
      for (;;) {
        do i++; while (mesh[index_arr[i].val][index] > mesh[a.val][index]);
        do j--; while (mesh[index_arr[j].val][index] < mesh[a.val][index]);
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

void qsort_idx_desc(INDEX *index_arr, const size_t N) {
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
        a = index_arr[j].index;
        for (i = j - 1; i >= l; i--) {
          if (index_arr[index_arr[i].index].val >= index_arr[a].val) break;
          index_arr[i + 1].index = index_arr[i].index;
        }
        index_arr[i + 1].index = a;
      }
      if (jstack < 0) break;
      ir = istack[jstack--];
      l = istack[jstack--];
    }
    else {
      k = (l + ir) >> 1;
      SWAP(index_arr[k].index, index_arr[l + 1].index, tmp);
      if (index_arr[index_arr[l].index].val < index_arr[index_arr[ir].index].val) {
        SWAP(index_arr[l].index, index_arr[ir].index, tmp);
      }
      if (index_arr[index_arr[l + 1].index].val < index_arr[index_arr[ir].index].val) {
        SWAP(index_arr[l + 1].index, index_arr[ir].index, tmp);
      }
      if (index_arr[index_arr[l].index].val < index_arr[index_arr[l + 1].index].val) {
        SWAP(index_arr[l].index, index_arr[l + 1].index, tmp);
      }
      i = l + 1;
      j = ir;
      a = index_arr[l + 1].index;
      for (;;) {
        do i++; while (index_arr[index_arr[i].index].val > index_arr[a].val);
        do j--; while (index_arr[index_arr[j].index].val < index_arr[a].val);
        if (j < i) break;
        SWAP(index_arr[i].index, index_arr[j].index, tmp);
      }
      index_arr[l + 1].index = index_arr[j].index;
      index_arr[j].index = a;
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

size_t binary_search(INDEX *index_arr, size_t l, size_t r, size_t x, size_t notfound) {
  size_t m;
  
  while (l <= r) { 
    m = l + (r - l) / 2; 
    //printf("\n(%d %d %d)", (int)l, (int)m ,(int)r);
    
    // Check if x is present at mid 
    
    if (index_arr[index_arr[m].index].val == x)
      return m; 
    // If x smaller, ignore left half 
    if (index_arr[index_arr[m].index].val > x) 
      l = m + 1; 

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

size_t binary_search_pos(FFT_CMPLX *mesh, INDEX *index_arr, long int l, long int r, FFT_REAL x, size_t notfound) {
  long int m;
  long int max_index = r;
  long int min_index = l;
  
  while (l <= r) { 
    m = l + (r - l) / 2; 
    // Check if x is present at mid 
    if (mesh[index_arr[m].val][0] == x) 
      return m; 

    // If x smaller, ignore left half 
    if (mesh[index_arr[m].val][0] > x) {
      l = m + 1;
      if(l <= max_index && mesh[index_arr[l].val][0] < x)
        return m;
      
      if(l > max_index)
        return notfound;
    }
    // If x is larger, ignore right half 
    else
      r = m - 1;
      if(r >= min_index && mesh[index_arr[r].val][0] > x)
        return r;
      
      if(r < min_index)
        return notfound;
  } 

  // if we reach here, then element was 
  // not present
  return notfound; 
} 

double return_x(double y) {
  if(y < 0.5) return sqrt(2 * y) - 1;
  else return 1 - sqrt(2 * (1 - y));
}

static void part_cic(FFT_CMPLX *mesh, INDEX *index_arr, size_t idx_mesh, FFT_REAL delta_mesh, const size_t Nh, long int *max_index) {
  //printf("\nIN Before sub ;part_cic5: %ld %0.10lf, %0.10lf ", (long int) idx_mesh, mesh[idx_mesh][0], delta_mesh);
  
  //mesh[idx_mesh][0] = mesh[idx_mesh][0] - delta_mesh;
  //printf("\nIN After sub; part_cic5: %ld %0.10lf, %0.10lf ", (long int) idx_mesh, mesh[idx_mesh][0], delta_mesh);
  
  size_t idx_idx;
  size_t idx_idx_idx;
  size_t idx_idx1;
  size_t min_rank;
  INDEX tmp;

  // Search in the array of length max_index+1 with the max_index+1 largest elements from the mesh, for the neighbour of the maximum element of the mesh
  idx_idx_idx = binary_search(index_arr, 0, *max_index, idx_mesh, 8 * Nh);
  //printf("\n IN part cic 5: %d %d %d", (long int) idx_idx_idx, (long int)  index_arr[idx_idx_idx].index, (long int)  index_arr[index_arr[idx_idx_idx].index].val);
  min_rank = index_arr[*max_index].rank;
  
  // If the neighbour is in this array then.
  if(idx_idx_idx != 8 * Nh) {
    idx_idx = index_arr[idx_idx_idx].index;
    mesh[idx_mesh][0] = mesh[idx_mesh][0] - delta_mesh;
    
    if(idx_idx + 1 < *max_index && mesh[index_arr[idx_idx].val][0] < mesh[index_arr[idx_idx+1].val][0]) {
      // Search the new position of the modified element of the mesh.
      idx_idx1 = binary_search_pos(mesh, index_arr, (long int) idx_idx + 1, *max_index, mesh[index_arr[idx_idx].val][0], 8 * Nh);
      if(idx_idx1 != 8 * Nh) {
        tmp = index_arr[idx_idx];
        memmove(index_arr + idx_idx, index_arr + idx_idx + 1, (idx_idx1 - idx_idx)*sizeof(INDEX));
        index_arr[idx_idx1] = tmp;
        min_rank = index_arr[*max_index].rank;
      }
      else {
        min_rank = index_arr[idx_idx].rank;
        memmove(index_arr + idx_idx, index_arr + idx_idx + 1, (*max_index - idx_idx)*sizeof(INDEX));
      }
    }
    else if (idx_idx + 1 == *max_index && mesh[index_arr[idx_idx].val][0] < mesh[index_arr[idx_idx+1].val][0]) {
      min_rank = index_arr[idx_idx].rank;
      index_arr[idx_idx] = index_arr[idx_idx + 1];
    }
  }
  *max_index = *max_index - 1;

#ifdef OMP
#pragma omp parallel for
#endif
  for (size_t i = 0; i <= *max_index; i++) {
    if(index_arr[i].rank > min_rank) {
      index_arr[i].rank = index_arr[i].rank - 1;
    }
    index_arr[index_arr[i].rank].index = i;
    
  }
}

static void cic(FFT_CMPLX *mesh, INDEX *index_arr, double x, double y, double z, const int Ng, const size_t Nh, const double Lbox, long int *max_index) {
    
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

void select_dens(FFT_CMPLX *mesh, HALOS *halos, INDEX *index_arr, const int Ng, const size_t Nh, const double Lbox) {
  size_t i, j, k, idx_max, Ntot, idx_val;
  long int max_index;
  double x, y, z;
  INDEX tmp;

  Ntot = (size_t) Ng * Ng * Ng;

  gsl_rng *r;
  r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, 49823);  

#ifdef OMP
#pragma omp parallel for
#endif
  for (i = 0; i < 8 * Nh; i++) {
    index_arr[i].val = i;
    index_arr[i].rank = 8 * Nh - 1 - i;
    index_arr[i].index = 0;
  }

  
  qsort_dens_asc(mesh, index_arr, 8 * Nh, 0);
    
  for (i = 8 * Nh; i < Ntot; i++) {
    if (mesh[i][0] > mesh[index_arr[0].val][0]) {
      index_arr[0].val = i;
      for (j = 0;;) {
        k = (j << 1) + 1;
        if (k > 8 * Nh - 1) break;
        if (k != 8 * Nh - 1 && mesh[index_arr[k].val][0] > mesh[index_arr[k + 1].val][0]) ++k;
        if (mesh[index_arr[j].val][0] <= mesh[index_arr[k].val][0]) break;
        SWAP(index_arr[k], index_arr[j], tmp);
        j = k;
      }
    }
  }
  qsort_dens_desc(mesh, index_arr, 8 * Nh, 0);

#ifdef OMP
#pragma omp parallel for
#endif
  for (i = 0; i < 8 * Nh; i++) index_arr[i].index = i;

  qsort_idx_desc(index_arr, 8 * Nh);

#ifdef OMP
#pragma omp parallel for
#endif
  for (i = 0; i < 8 * Nh; i++) index_arr[index_arr[i].index].rank = i;


  // for(i = 0; i < 8 * Nh; i++)
  // {
  //   printf("\n%d %lf %lf %lf %lf", (int)i, (double)index_arr[i].val, (double)index_arr[index_arr[i].index].val, (double) index_arr[i].rank, mesh[index_arr[i].val][0]);//, (double)index_arr[8*Nh-100 + i].index, (double)index_arr[index_arr[8*Nh-100 + i].index].val );
  // }
  printf("\nThe maximum density is %f ", mesh[index_arr[0].val][0]);
  printf("\nThe minimum density is %f ", mesh[index_arr[8*Nh - 1].val][0]);
  fflush(stdout);
  //printf("\nPozitie: %lf %lf", (double)binary_search(index_arr, 0, 8 * Nh - 1, MESH_IDX(Ng, 2, 2, 1), 8 * Nh), (double) MESH_IDX(Ng, 2, 2, 1));
  //FILE *ofile;
  //ofile = fopen("./output/test.txt", "w+");

  // max_index = 8 * Nh - 1;
  // for (size_t u = 0; u < Nh; u++) {
  //   idx_max = 0;
  //   idx_val = index_arr[idx_max].val;
    
  //   //fprintf(ofile, "######\n");
  //   //for (size_t s = 0; s < max_index; i++)
  //   //  fprintf(ofile, "%ld %ld %0.10lf %ld \n", (long int) s, (long int) index_arr[s].val, mesh[index_arr[s].val][0], (long int) (max_index));
  //   //fprintf(ofile, "\n");
  //   if(u%1000==0)
  //     printf("\n%ld %ld", (long int)u, (long int) max_index);
  //   k = idx_val/(size_t)(Ng * Ng);
  //   j = (idx_val - (size_t)k * Ng * Ng) / (size_t)Ng;
  //   i = idx_val - (size_t)k * Ng * Ng - (size_t)j * Ng;
    
  //   // printf("\nidx_val=%d i=%d j=%d k=%d mesh_idx=%d", idx_val, i ,j ,k, MESH_IDX(Ng, i,j , k));
  //   x = gsl_ran_flat(r, 0, 1);
  //   y = gsl_ran_flat(r, 0, 1);
  //   z = gsl_ran_flat(r, 0, 1);
    
  //   halos[u].x[0] = Lbox * (i + return_x(x)) / Ng;
  //   halos[u].x[1] = Lbox * (j + return_x(y)) / Ng;
  //   halos[u].x[2] = Lbox * (k + return_x(z)) / Ng;
  //   halos[u].dens = idx_val;
    
  //   halos[u].x[0] = (halos[u].x[0] < 0) ? (Lbox + halos[u].x[0]) : halos[u].x[0];
  //   halos[u].x[1] = (halos[u].x[1] < 0) ? (Lbox + halos[u].x[1]) : halos[u].x[1];
  //   halos[u].x[2] = (halos[u].x[2] < 0) ? (Lbox + halos[u].x[2]) : halos[u].x[2];
  //   //printf("\n\n %d Before cic: %ld %0.10lf ", (int) u, (long int) index_arr[idx_max].val, mesh[idx_val][0]);
  //   cic(mesh, index_arr, halos[u].x[0], halos[u].x[1], halos[u].x[2], Ng, Nh, Lbox, &max_index);
  //   //printf("\nAfter cic: %ld %0.10lf ", (long int) index_arr[idx_max].val, mesh[idx_val][0]);
  //   //max_index = max_index - 8;
  // }
  // //fclose(ofile);
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
