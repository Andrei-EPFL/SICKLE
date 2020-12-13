#include "load_conf.h"
#include "linhalo.h"
#include "fftw_define.h"
#include <unistd.h>

int main(int argc, char *argv[]) {
  int ecode;
  double *k, *P;
  size_t Nk;
  MAX_DENS *max_dens;
  HALOS *halos;
  CONF conf;

  FFT_CMPLX *mesh;
  FFT_PLAN fp;

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
  
  /*Initizalizations */
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
    P_EXT("failed to initialize the halos.\n");
    return ecode;
  }
  printf(FMT_DONE);

  if ((ecode = init_max_dens(conf.Nhalo, &max_dens))) {
    P_EXT("failed to initialize the array of indices.\n");
    return ecode;
  }
  printf(FMT_DONE);
  
  printf("Generating the density field... ");
  fflush(stdout);
  if ((ecode = init_field(conf.Ngrid, &mesh, &fp))) {
    P_EXT("failed to initialize the density field.\n");
    return ecode;
  }
  time_t end = time(NULL);
  float seconds = (float)(end - start);
  printf("1. Up to this point %f seconds passed.\n The counter is set to zero again.\n", seconds);
  fflush(stdout);
  
  /* Generation of the gaussian random field */
  start = time(NULL);
  if ((ecode = gauss_ran_field(&conf, k, P, Nk, &fp, mesh, conf.factor))) {
    P_EXT("failed to generate the density field.\n");
    return ecode;
  }
  free(k);
  free(P);

  printf(FMT_DONE);
  if (conf.savedm) {
    printf("Saving the density field ... ");
    fflush(stdout);
    if ((ecode = save_dens(conf.dmout, mesh, conf.Ngrid))) {
      P_EXT("failed to save the density field.\n");
      return ecode;
    }
    printf(FMT_DONE);
  }
  end = time(NULL);
  seconds = (float)(end - start);
  printf("2. Up to this point %f seconds passed.\n The counter is set to zero again.\n", seconds);
  fflush(stdout);


  /*Populate the mesh with halos*/
  start = time(NULL);
  printf("Populating the mesh with haloes... ");
  fflush(stdout);
  select_dens(mesh, halos, max_dens, conf.Ngrid, conf.Nhalo, conf.Lbox);
  printf(FMT_DONE);
  end = time(NULL);
  seconds = (float)(end - start);
  printf("\n3. Up to this point %f seconds passed.\n", seconds);
  fflush(stdout);
  
  /* Save halos */
  printf("Saving haloes ... ");
  fflush(stdout);
  if ((ecode = save_halo(conf.output, halos, conf.Nhalo))) {
    P_EXT("failed to generate the halo catalogue.\n");
    return ecode;
  }
  free(halos);
  printf(FMT_DONE);

  FFT_FREE(mesh);
  time_t endall = time(NULL);
  printf("The code has finished in %f seconds \n", (float)(endall - startall));
  return 0;
}
