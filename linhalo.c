#include "load_conf.h"
#include "linhalo.h"
#include <unistd.h>

int main(int argc, char *argv[]) {
  int ecode;
  double *k, *P;
  size_t Nk;
  HALOS *halos;
  CONF conf;
#ifdef DOUBLE_PREC
  fftw_complex *mesh;
  fftw_plan fp;
#else
  fftwf_complex *mesh;
  fftwf_plan fp;
#endif
  time_t startall = time(NULL);
  time_t start = time(NULL);

  printf("Loading configurations ... ");
  fflush(stdout);
  if ((ecode = read_conf(argc, argv, &conf))) {
    P_EXT("failed to load the configurations.\n");
    return ecode;
  }
  if ((ecode = check_conf(&conf))) {
    P_EXT("please check your configurations.\n");
    return ecode;
  }
  print_conf(&conf);
  printf(FMT_DONE);

  printf("Processing the linear power spectrum ... ");
  fflush(stdout);
  if ((ecode = read_pk(conf.pkfile, &k, &P, &Nk))) {
    P_EXT("please check your input linear power spectrum.\n");
    return ecode;
  }
  if ((ecode = init_pk(k, P, Nk, &conf))) {
    P_EXT("please check your input linear power spectrum.\n");
    return ecode;
  }
  printf(FMT_DONE);
  
  if ((ecode = init_halos(conf.Nhalo, &halos))) {
    P_EXT("failed to generate the halos.\n");
    return ecode;
  }
  printf(FMT_DONE);
  
  printf("Generating the density field and the inverted density field... ");
  fflush(stdout);
  if ((ecode = init_field(conf.Ngrid, &mesh, &fp))) {
    P_EXT("failed to generate the density field.\n");
    return ecode;
  }
  time_t end = time(NULL);
  float seconds = (float)(end - start);
  printf("Up to there it took %f seconds\n", seconds);
  fflush(stdout);
    
  start = time(NULL);
  if ((ecode = gauss_ran_field(&conf, k, P, Nk, &fp, mesh))) {
    P_EXT("failed to generate the density field.\n");
    return ecode;
  }
  free(k);
  free(P);
  printf(FMT_DONE);
  end = time(NULL);
  seconds = (float)(end - start);
  printf("Up to there it took %f seconds\n", seconds);
  fflush(stdout);
  
  start = time(NULL);
  
  printf("Populating with haloes the normal mesh... ");
  fflush(stdout);
  select_dens(mesh, halos, conf.Ngrid, conf.Nhalo, conf.Lbox);
  fflush(stdout);
  printf(FMT_DONE);
  end = time(NULL);
  seconds = (float)(end - start);
  printf("Up to there it took %f seconds\n", seconds);
  fflush(stdout);
  
  // printf("Saving haloes ... ");
  // fflush(stdout);
  // if ((ecode = save_halo(conf.output, halos, conf.Nhalo))) {
  //   P_EXT("failed to generate the halo catalogue.\n");
  //   return ecode;
  // }
  free(halos);
  printf(FMT_DONE);

  // if (conf.savedm) {
  //   printf("Saving the density field ... ");
  //   fflush(stdout);
  //   if ((ecode = save_dens(conf.dmout, mesh, conf.lowNg))) {
  //     P_EXT("failed to save the density field.\n");
  //     return ecode;
  //   }
  //   printf(FMT_DONE);
  // }

#ifdef DOUBLE_PREC
  fftw_free(mesh);
#else
  fftwf_free(mesh);
#endif
  time_t endall = time(NULL);
  printf("Everything took %f seconds\n", (float)(endall - startall));
  return 0;
}
