#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "load_conf.h"

int read_conf(const int argc, char *const *argv, CONF *conf) {
  int ecode, optidx;
  coca *cfg;
  const int nfunc = 2;
  const coca_funcs funcs[2] = {
    {   'h',    "help",         usage,          argv[0] },
    {   't',    "template",     temp_conf,      NULL    }
  };

  cfg = coca_init();
  if (!cfg) {
    P_ERR("failed to allocate memory for configurations.\n");
    return ERR_MEM;
  }

  ecode = coca_set_param(cfg, 'c', "conf",   "CONFIG", COCA_VERBOSE_WARN);
  CHECK_COCA(ecode, CONFIG);
  ecode = coca_set_param(cfg, 'i', "input",  "PK_INPUT", COCA_VERBOSE_WARN);
  CHECK_COCA(ecode, PK_INPUT);
  ecode = coca_set_param(cfg, 'g', "grid",   "NGRID", COCA_VERBOSE_WARN);
  CHECK_COCA(ecode, NGRID);
  ecode = coca_set_param(cfg, 'b', "box",    "BOXSIZE", COCA_VERBOSE_WARN);
  CHECK_COCA(ecode, BOXSIZE);
  ecode = coca_set_param(cfg, 's', "seed",   "SEED", COCA_VERBOSE_WARN);
  CHECK_COCA(ecode, SEED);
  ecode = coca_set_param(cfg, 'n', "num",    "NHALO", COCA_VERBOSE_WARN);
  CHECK_COCA(ecode, NHALO);
  ecode = coca_set_param(cfg, 'l', "mksurvey", "MKSURVEY", COCA_VERBOSE_WARN);
  CHECK_COCA(ecode, MKSURVEY)
  ecode = coca_set_param(cfg, 'o', "output", "OUTPUT", COCA_VERBOSE_WARN);
  CHECK_COCA(ecode, OUTPUT);
  ecode = coca_set_param(cfg, 'a', "factor", "FACTOR", COCA_VERBOSE_WARN);
  CHECK_COCA(ecode, OUTPUT);
  ecode = coca_set_param(cfg, 'd', "savedm", "SAVE_DENS", COCA_VERBOSE_WARN);
  CHECK_COCA(ecode, SAVE_DENS);
  ecode = coca_set_param(cfg, 'm', "dmout",  "DENS_OUT", COCA_VERBOSE_WARN);
  CHECK_COCA(ecode, DENS_OUT);
  ecode = coca_set_param(cfg, 'f', "force",  "FORCE", COCA_VERBOSE_WARN);
  CHECK_COCA(ecode, FORCE);

  ecode = coca_read_opts(cfg, argc, argv, funcs, nfunc, &optidx,
      100, COCA_VERBOSE_WARN);
  if (ecode) {
    P_ERR("failed to read command line options.\n");
    return ecode;
  }

  ecode = coca_get_var(cfg, "CONFIG", COCA_DTYPE_STR, &(conf->config));
  if (ecode == COCA_ERR_NOT_SET) {
    if (safe_strcpy(conf->config, DEFAULT_CONF_FILE, COCA_MAX_VALUE_LEN)) {
      P_ERR("failed to set the default configuration file.\n");
      return ERR_STRING;
    }
  }
  else CHECK_COCA(ecode, CONFIG);

  ecode = coca_read_file(cfg, conf->config, 1, COCA_VERBOSE_WARN);
  if (ecode) {
    P_ERR("failed to read the configuration file.\n");
    return ecode;
  }

  ecode = coca_get_var(cfg, "PK_INPUT", COCA_DTYPE_STR, &(conf->pkfile));
  CHECK_COCA(ecode, PK_INPUT);
  ecode = coca_get_var(cfg, "NGRID", COCA_DTYPE_INT, &(conf->Ngrid));
  CHECK_COCA(ecode, NGRID);
  ecode = coca_get_var(cfg, "BOXSIZE", COCA_DTYPE_DBL, &(conf->Lbox));
  CHECK_COCA(ecode, BOXSIZE);
  ecode = coca_get_var(cfg, "FACTOR", COCA_DTYPE_DBL, &(conf->factor));
  CHECK_COCA(ecode, BOXSIZE);
  ecode = coca_get_var(cfg, "SEED", COCA_DTYPE_LONG, &(conf->seed));
  CHECK_COCA(ecode, SEED);
  ecode = coca_get_var(cfg, "NHALO", COCA_DTYPE_LONG, &(conf->Nhalo));
  CHECK_COCA(ecode, NHALO);
  ecode = coca_get_var(cfg, "OUTPUT", COCA_DTYPE_STR, &(conf->output));
  CHECK_COCA(ecode, OUTPUT);
  ecode = coca_get_var(cfg, "SAVE_DENS", COCA_DTYPE_BOOL, &(conf->savedm));
  CHECK_COCA(ecode, SAVE_DENS);
  ecode = coca_get_var(cfg, "MKSURVEY", COCA_DTYPE_STR, &(conf->mksurvey));
  if (ecode == COCA_ERR_NOT_SET) conf->mksurvey[0] = '\0';
  else CHECK_COCA(ecode, MKSURVEY);
  ecode = coca_get_var(cfg, "DENS_OUT", COCA_DTYPE_STR, &(conf->dmout));
  if (ecode == COCA_ERR_NOT_SET) conf->dmout[0] = '\0';
  else CHECK_COCA(ecode, DENS_OUT);
  ecode = coca_get_var(cfg, "FORCE", COCA_DTYPE_BOOL, &(conf->force));
  if (ecode == COCA_ERR_NOT_SET) conf->force = false;
  else CHECK_COCA(ecode, FORCE);

  coca_destroy(cfg);
  return 0;
}

int check_conf(CONF *conf) {
  int ecode;

  if ((ecode = check_input(conf->pkfile, "PK_INPUT"))) return ecode;

  if (conf->Ngrid <= 1) {
    P_ERR(FMT_KEY(NGRID) " is not set correctly.\n");
    return ERR_RANGE;
  }
  
  if (conf->Lbox <= 0) {
    P_ERR(FMT_KEY(BOXSIZE) " is not set correctly.\n");
    return ERR_RANGE;
  }
  if (conf->factor <= 0) {
    P_ERR(FMT_KEY(FACTOR) " is not set correctly.\n");
    return ERR_RANGE;
  }
  if (conf->seed < 0) {
    P_ERR(FMT_KEY(SEED) " must be non-negative.\n");
    return ERR_RANGE;
  }
  if (conf->Nhalo <= 0) {
    P_ERR(FMT_KEY(NHALO) " is not set correctly.\n");
    return ERR_RANGE;
  }
  if (conf->Nhalo > (size_t) conf->Ngrid * conf->Ngrid * conf->Ngrid / 2) {
    P_ERR(FMT_KEY(NHALO) " is larger than " FMT_KEY(NGRID) "^3 / 2.\n");
    return ERR_RANGE;
  }
  if (conf->Nhalo > (size_t) conf->Ngrid * conf->Ngrid * conf->Ngrid) {
    P_ERR(FMT_KEY(NHALO) " is larger than " FMT_KEY(NGRID) "^3.\n");
    return ERR_RANGE;
  }

  if ((ecode = check_output(conf->output, "OUTPUT", conf->force))) return ecode;
  if (conf->savedm)
    if ((ecode = check_output(conf->dmout, "DENS_OUT", conf->force)))
      return ecode;
  if (*conf->mksurvey && (ecode = check_input(conf->mksurvey, "MKSURVEY")))
    return ecode;
  return 0;
}

int check_input(const char *fname, const char *dscp) {
  if (fname[0] == '\0' || fname[0] == ' ') {
    P_ERR("the input " FMT_KEY(%s) " is not set.\n", dscp);
    return ERR_FILE;
  }
  if (access(fname, R_OK)) {
    P_ERR("cannot open " FMT_KEY(%s) ": `%s'.\n", dscp, fname);
    return ERR_FILE;
  }
  return 0;
}

int check_output(const char *fname, const char *dscp, const int force) {
  char confirm, *end;
  mystr path;
  int cnt = 0;

  if (fname[0] == '\0' || fname[0] == ' ') {
    P_ERR("the output " FMT_KEY(%s) " is not set.\n", dscp);
    return ERR_FILE;
  }
  if (!access(fname, F_OK) && force == 0) {     // If the output file exists.
    P_WRN("the output " FMT_KEY(%s) " `%s' exists.\n", dscp, fname);
    do {
      if ((++cnt) == TIMEOUT) {
        P_ERR("too many failed inputs.\n");
        return ERR_INPUT;
      }
      fprintf(stderr, "Are you going to overwrite it? (y/n): ");
      if (scanf("%c", &confirm) != 1) continue;
      while(getchar() != '\n');         // Ignore invalid inputs.
    }
    while (confirm != 'y' && confirm != 'n');
    if(confirm == 'n') {
      P_ERR("cannot write to the file.\n");
      return ERR_FILE;
    }
  }

  if (!access(fname, F_OK) && access(fname, W_OK)) {
    P_ERR("cannot write to file `%s'.\n", fname);
    return ERR_FILE;
  }

  // Check the path of the output file.
  if (safe_strcpy(path, fname, COCA_MAX_VALUE_LEN)) {
    P_ERR("failed to get the path of the output file.\n");
    return ERR_STRING;
  }
  if ((end = strrchr(path, '/')) != NULL) {
    *(end + 1) = '\0';
    if(access(path, X_OK)) {
      P_ERR("cannot access the path `%s'.\n", path);
      return ERR_FILE;
    }
  }

  return 0;
}

void print_conf(const CONF *conf) {
  printf("\n  Configuration file: %s\n", conf->config);
  printf("  " FMT_KEY(PK_INPUT) " = %s\n", conf->pkfile);
  printf("  " FMT_KEY(NGRID) " = %d\n", conf->Ngrid);
  printf("  " FMT_KEY(BOXSIZE) " = %g\n", conf->Lbox);
  printf("  " FMT_KEY(FACTOR) " = %g\n", conf->factor);
  printf("  " FMT_KEY(SEED) " = %ld\n", conf->seed);
  printf("  " FMT_KEY(NHALO) " = %ld\n", conf->Nhalo);
  if (*conf->mksurvey) printf("  " FMT_KEY(MKSURVEY) " = %s\n", conf->mksurvey);
  printf("  " FMT_KEY(OUTPUT) " = %s\n", conf->output);
  printf("  " FMT_KEY(SAVE_DENS) " = %d\n", conf->savedm);
  if (conf->savedm)
    printf("  " FMT_KEY(DENS_OUT) " = %s\n", conf->dmout);
  printf("  " FMT_KEY(FORCE) " = %d\n", conf->force);
}

void usage(void *pname) {
  printf("Usage: %s [OPTION [VALUE]]\n\
Generate a halo catalogue given an input linear matter power spectrum.\n\n\
  -h, --help\n\
        Display this message and exit.\n\
  -t, --template\n\
        Display a template configuration file and exit.\n\
  -c, --conf\n\
        Set the configuration file (default: `%s').\n\
  -i, --input\n\
        Set the input file for the linear power spectrum.\n\
  -g, --grid\n\
        Set the dimension of grids for the linear density field.\n\
  -a, --factor\n\
        Set the factor which multiplies the linear power spectrum.\n\
  -b, --box\n\
        Set the side length of the simulation box.\n\
  -s, --seed\n\
        Set the random seed.\n\
  -n, --num\n\
        Set the number of output haloes.\n\
  -o, --output\n\
        Set the output file for the halo catalogue.\n\
  -d, --savedm\n\
        Save the low resolution dark matter density field.\n\
  -m, --dmout\n\
        Set the output file for the dark matter density field.\n\
  -f, --force\n\
        Force overwriting existing files without notifications.\n\
Consult the -t option for more information on the configuraion.\n\
Report bugs to <zhaocheng03@gmail.com>.\n",
      (char *) pname, DEFAULT_CONF_FILE);
  exit(0);
}

void temp_conf(void *arg) {
  printf("# Configuration file (default: `%s').\n\
# NOTE that command line options have priority over this file.\n\
# FORMAT: KEYWORD = VALUE # COMMENT\n\
# See https://github.com/cheng-zhao/COCA for the detailed format.\n\
\n\
PK_INPUT        = \n\
        # Input linear matter power spectrum.\n\
        # The first 2 columns must be {k, P(k)}.\n\
NGRID           = \n\
        # The number of grids on each side for the simulation box.\n\
BOXSIZE         = \n\
        # The side length of the simulation box.\n\
FACTOR          = \n\
        # The factor which multiplies the linear power spectrum.\n\
SEED            = \n\
        # The random seed, must be non-negative.\n\
NHALO           = \n\
        # The number of haloes to be populated.\n\
OUTPUT          = \n\
        # The output file for the haloes.\n\
SAVE_DENS       = \n\
        # True for saving the lower resolution dark matter density.\n\
DENS_OUT        = \n\
        # The output file for the dark matter density field.\n\
FORCE           = \n\
        # True for force overwriting the output (default: %s).\n",
      DEFAULT_CONF_FILE, DEFAULT_FORCE);
  exit(0);
}
