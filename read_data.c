#include "linhalo.h"
#include <ctype.h>

int read_pk(const char *fname, double **k, double **P, size_t *num) {
  FILE *fp;
  char *buf, *p, *end, *endl;
  size_t n, cnt, nrest;

  printf("\n  Filename: %s\n  Counting lines ...", fname);
  fflush(stdout);

  if (!(fp = fopen(fname, "rb"))) {
    P_ERR("cannot open file: `%s'\n", fname);
    return ERR_FILE;
  }

  n = 0;
  buf = malloc(sizeof *buf * CHUNK);
  if (!buf) {
    P_ERR("failed to allocate memory for reading the file.\n");
    return ERR_MEM;
  }

  while ((cnt = fread(buf, sizeof(char), CHUNK, fp))) {
    p = buf;
    end = p + cnt;
    while ((p = memchr(p, '\n', end - p))) {
      ++p;
      ++n;
    }
  }

  printf("\r  Number of records: %ld.\n", n);
  *k = malloc(sizeof(double *) * n);
  *P = malloc(sizeof(double *) * n);
  
  if (!(*k) || !(*P)) {
    P_ERR("failed to allocate memory for the power spectrum.\n");
    return ERR_MEM;
  }

  printf("  Reading ... ");
  fflush(stdout);

  fseek(fp, 0, SEEK_SET);
  n = nrest = 0;
  while ((cnt = fread(buf + nrest, sizeof(char), CHUNK - nrest, fp))) {
    p = buf;
    end = p + nrest + cnt;
    if (cnt < CHUNK - nrest) *end = '\n';       // add '\n' to the last line
    while ((endl = memchr(p, '\n', end - p))) {
      *endl = '\0';     // replace '\n' by '\0' for processing lines
      while (isspace(p[0])) ++p;                // remove leading whitespaces
      if (p[0] == COMMENT || p[0] == '\0') {    // COMMENT or empty line
        p = endl + 1;
        continue;
      }

      if (sscanf(p, "%lf %lf", *k + n, *P + n) != 2) {
        P_ERR("failed to read line: %s\n", p);
        return ERR_FILE;
      }

      ++n;
      p = endl + 1;
    }

    nrest = end - p;
    memmove(buf, p, nrest);
  }

  *num = n;
  fclose(fp);
  free(buf);

  printf("\r  %ld valid lines recorded.\n", n);
  return 0;
}
