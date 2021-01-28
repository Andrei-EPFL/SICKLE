#ifndef _LOAD_CONF_H_
#define _LOAD_CONF_H_

#include "define.h"
#include "coca.h"

typedef struct {
  bool force;
  bool savedm;
  int Ngrid;
  double factor;
  double Lbox;
  long Nhalo;
  long seed;
  mystr config;
  mystr pkfile;
  mystr mksurvey;
  mystr output;
  mystr dmout;
} CONF;

int read_conf(const int, char * const *, CONF *);

int check_conf(CONF *);

int check_input(const char *, const char *);

int check_output(const char *, const char *, const int);

void print_conf(const CONF *);

void temp_conf(void *);

void usage(void *);

#endif
