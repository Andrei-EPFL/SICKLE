/*******************************************************************************
*
* COCA: Command line Option and Configuration file Analyser
*
* Load configurations from command line options and plain text files.
*
* Github repository:
*       https://github.com/cheng-zhao/COCA
*
* Copyright (c) 2019 Cheng Zhao <zhaocheng03@gmail.com>
* 
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
* 
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "coca.h"

coca *coca_init(void) {
  coca *cfg;
  cfg = calloc(1, sizeof *cfg);
  if (!cfg) return NULL;
  cfg->next = NULL;
  cfg->src = COCA_SRC_NULL;
  return cfg;
}

void coca_destroy(coca *cfg) {
  if (!cfg) return;
  coca_destroy(cfg->next);
  free(cfg);
}

int coca_set_param(coca *cfg, const char opt, const char *long_opt,
    const char *name, const int verbose) {
  coca *conf = NULL;
  int i;
  bool short_valid, long_valid;

  /* Arguments validation. */
  if (!cfg) return COCA_ERR_INIT;
  if (!name || name[0] == '\0') return COCA_ERR_ARG;
  short_valid = long_valid = false;

  for (i = 0; i < COCA_MAX_KEY_LEN; i++) {
    if (name[i] == '\0') break;
    if (!isalnum(name[i]) && name[i] != '_') return COCA_ERR_ARG;
  }
  if (name[i] != '\0') return COCA_ERR_ARG;

  if (isalpha(opt)) short_valid = true;
  else if (opt != 0 && verbose >= COCA_VERBOSE_WARN)
    COCA_MSG_WARN("omitting invalid short option of parameter: %s\n", name);

  if (long_opt && long_opt[0] != '\0') {
    for (i = 0; i < COCA_MAX_KEY_LEN; i++) {
      if (long_opt[i] == '\0') {
        long_valid = true;
        break;
      }
      if (!isgraph(long_opt[i]) || long_opt[i] == COCA_CMD_ASSIGN) {
        if (verbose >= COCA_VERBOSE_WARN)
          COCA_MSG_WARN("omitting invalid long option of parameter: %s\n",
              name);
        break;
      }
    }
    if (i == COCA_MAX_KEY_LEN - 1 && !long_valid) return COCA_ERR_ARG;
  }

  /* Check duplicates and intialise a new node if necessary. */
  if (cfg->name[0] != '\0') {           /* not a fresh linked list */
    do {
      if (!conf) conf = cfg;
      else conf = conf->next;

      if (!strncmp(name, conf->name, COCA_MAX_KEY_LEN)) return COCA_ERR_EXIST;
      if (short_valid && opt == conf->opt) return COCA_ERR_EXIST;
      if (long_valid && !strncmp(long_opt, conf->long_opt, COCA_MAX_KEY_LEN))
        return COCA_ERR_EXIST;
    }
    while (conf->next);

    if (!(conf->next = coca_init())) return COCA_ERR_MEMORY;
    conf = conf->next;
  }
  else conf = cfg;

  if (short_valid) conf->opt = opt;
  if (long_valid) {
    if (safe_strcpy(conf->long_opt, long_opt, COCA_MAX_KEY_LEN))
      return COCA_ERR_STRING;
  }
  if (safe_strcpy(conf->name, name, COCA_MAX_KEY_LEN))
    return COCA_ERR_STRING;

  return 0;
}

int coca_check_funcs(coca *cfg, const coca_funcs *funcs, const int nfunc,
    coca_func_valid *fval, const int verbose) {
  int i, m, n;
  coca *conf;

  for (n = 0; n < nfunc; n++) {
    if (isalpha(funcs[n].opt)) {
      fval[n].opt = funcs[n].opt;
      /* Check duplicates. */
      for (m = 0; m < n; m++)
        if (fval[m].opt == fval[n].opt) return COCA_ERR_EXIST;
      for (conf = cfg; conf != NULL; conf = conf->next)
        if (fval[n].opt == conf->opt) return COCA_ERR_EXIST;
    }
    else if (funcs[n].opt != 0 && verbose >= COCA_VERBOSE_WARN)
      COCA_MSG_WARN("omitting invalid function short option with code %d\n",
          funcs[n].opt);

    if (funcs[n].long_opt && funcs[n].long_opt[0] != '\0') {
      for (i = 0; i < COCA_MAX_KEY_LEN; i++) {
        if (funcs[n].long_opt[i] == '\0') {
          fval[n].long_opt = funcs[n].long_opt;
          /* Check duplicates. */
          for (m = 0; m < n; m++) {
            if (fval[m].long_opt && !strncmp(fval[m].long_opt,
                  fval[n].long_opt, COCA_MAX_KEY_LEN))
              return COCA_ERR_EXIST;
          }
          for (conf = cfg; conf != NULL; conf = conf->next) {
            if (!strncmp(fval[n].long_opt, conf->long_opt, COCA_MAX_KEY_LEN))
              return COCA_ERR_EXIST;
          }
          break;
        }
        if (!isgraph(funcs[n].long_opt[i]) ||
            funcs[n].long_opt[i] == COCA_CMD_ASSIGN) {
          if (verbose >= COCA_VERBOSE_WARN)
            COCA_MSG_WARN("omitting invalid function long option: %s\n",
                funcs[n].long_opt);
          break;
        }
      }
      if (i == COCA_MAX_KEY_LEN - 1 && !fval[n].long_opt)
        return COCA_ERR_ARG;
    }

    if (!(fval[n].opt || fval[n].long_opt)) return COCA_ERR_ARG;
  }
  return 0;
}

int coca_read_opts(coca *cfg, const int argc, char *const *argv,
    const coca_funcs *funcs, const int nfunc, int *optidx,
    const int prior, const int verbose) {
  int i, j, m, n, idx;
  bool with_arg = false;
  coca *conf;
  char *optarg;
  coca_func_valid *fval = NULL;

  if (prior <= COCA_SRC_NULL) return COCA_ERR_ARG;
  if (funcs == NULL && nfunc) return COCA_ERR_ARG;
  /* Do not check cfg, for command line options calling functions only. */

  /* Validation of function options. */
  if (nfunc) {
    if (!(fval = calloc(nfunc, sizeof *fval))) return COCA_ERR_MEMORY;
    if ((n = coca_check_funcs(cfg, funcs, nfunc, fval, verbose))) {
      free(fval);
      return n;
    }
  }

  /* Parsing command line options. */
  for (idx = 0, i = 1; i < argc; i++) {
    if (argv[i][0] != COCA_CMD_FLAG) {
      idx = i;
      break;
    }

    if (!(COCA_IS_OPT(argv[i]))) {
      if (verbose >= COCA_VERBOSE_WARN)
        COCA_MSG_WARN("invalid option: %s\n", argv[i]);
      continue;
    }

    optarg = NULL;
    with_arg = false;
    for (j = 0; j < COCA_MAX_VALUE_LEN; j++) {
      if (argv[i][j] == '\0' || argv[i][j] == COCA_CMD_ASSIGN) break;
    }
    if (argv[i][j] == COCA_CMD_ASSIGN) {
      with_arg = true;
      optarg = argv[i] + j + 1;
      argv[i][j] = '\0';
    }
    else if (argv[i][j] == '\0') {
      j = i + 1;
      if (j < argc && !(COCA_IS_OPT(argv[j]))) optarg = argv[j];
    }
    else {
      if (nfunc) free(fval);
      return COCA_ERR_STRING;
    }

    if (argv[i][1] != COCA_CMD_FLAG) {          /* short option */
      for (n = m = 0; n < nfunc; n++) {
        if (fval[n].opt == argv[i][1]) {
          if (optarg && verbose >= COCA_VERBOSE_WARN)
            COCA_MSG_WARN("omitting argument of option: %s\n", argv[i]);
          if (!fval[n].called) {
            funcs[n].func_ptr(funcs[n].func_arg);
            fval[n].called = true;
          }
          else if (verbose >= COCA_VERBOSE_WARN)
            COCA_MSG_WARN("duplicate option: %s\n", argv[i]);
          m = 1;
          break;
        }
      }
      if (m) continue;

      for (conf = cfg; conf != NULL; conf = conf->next) {
        if (conf->opt == argv[i][1]) {
          if (conf->src < prior) {
            if (optarg &&
                safe_strcpy(conf->value, optarg, COCA_MAX_VALUE_LEN)) {
              if (nfunc) free(fval);
              return COCA_ERR_STRING;
            }
          }
          else if (conf->src == prior && verbose >= COCA_VERBOSE_WARN)
            COCA_MSG_WARN("duplicate option: %s\n", argv[i]);
          conf->src = prior;
          break;
        }
      }
      if (!conf && verbose >= COCA_VERBOSE_WARN)
        COCA_MSG_WARN("invalid option: %s\n", argv[i]);
    }
    else if (argv[i][2]) {                         /* long option */
      for (n = m = 0; n < nfunc; n++) {
        if (fval[n].long_opt[n] &&
            !strncmp(fval[n].long_opt, argv[i] + 2, COCA_MAX_KEY_LEN)) {
          if (optarg && verbose >= COCA_VERBOSE_WARN)
            COCA_MSG_WARN("omitting argument of option: %s\n", argv[i]);
          if (!fval[n].called) {
            funcs[n].func_ptr(funcs[n].func_arg);
            fval[n].called = true;
          }
          else if (verbose >= COCA_VERBOSE_WARN)
            COCA_MSG_WARN("duplicate option: %s\n", argv[i]);
          m = 1;
          break;
        }
      }
      if (m) continue;

      for (conf = cfg; conf != NULL; conf = conf->next) {
        if (!strncmp(conf->long_opt, argv[i] + 2, COCA_MAX_KEY_LEN)) {
          if (conf->src < prior) {
            if (optarg &&
                safe_strcpy(conf->value, optarg, COCA_MAX_VALUE_LEN)) {
              if (nfunc) free(fval);
              return COCA_ERR_STRING;
            }
          }
          else if (conf->src == prior && verbose >= COCA_VERBOSE_WARN)
            COCA_MSG_WARN("duplicate option: %s\n", argv[i]);
          conf->src = prior;
          break;
        }
      }
      if (!conf && verbose >= COCA_VERBOSE_WARN)
        COCA_MSG_WARN("invalid option: %s\n", argv[i]);
    }
    else {                                      /* parser termination */
      idx = i + 1;
      break;
    }

    if (optarg && !with_arg) ++i;
  }

  if (idx == 0) idx = i + 1;
  *optidx = idx;
  if (idx < argc && verbose >= COCA_VERBOSE_DEBUG) {
    COCA_MSG_DEBUG("unused command line options:\n ");
    while (idx < argc) fprintf(stderr, " %s", argv[idx++]);
    fprintf(stderr, "\n");
  }

  if (nfunc) free(fval);
  return 0;
}

int coca_read_file(coca *cfg, const char *fname, const int prior,
    const int verbose) {
  FILE *fp;
  char line[COCA_MAX_LINE_LEN];
  char keyword[COCA_MAX_KEY_LEN];
  mystr value;
  int code;
  coca *conf;

  if (!cfg) return COCA_ERR_INIT;
  if (prior <= COCA_SRC_NULL) return COCA_ERR_ARG;
  if (!(fp = fopen(fname, "r"))) return COCA_ERR_FILE;
  memset(line, 0, sizeof(char) * COCA_MAX_LINE_LEN);
  while (fgets(line, COCA_MAX_LINE_LEN, fp) != NULL) {
    if (line[COCA_MAX_LINE_LEN - 1] != '\0' &&
        line[COCA_MAX_LINE_LEN - 1] != '\n' &&
        line[COCA_MAX_LINE_LEN - 1] != EOF) return COCA_ERR_STRING;
    code = coca_parse_line(line, keyword, value);

    if (code == COCA_PARSE_ERROR) {
      COCA_MSG_ERR("failed to parse line:\n  %s", line);
      return COCA_ERR_PARSE;
    }
    if (code == COCA_PARSE_WARN && verbose >= COCA_VERBOSE_WARN)
      COCA_MSG_WARN("failed to parse line:\n  %s", line);
    if (code == COCA_PARSE_DEBUG && verbose >= COCA_VERBOSE_DEBUG)
      COCA_MSG_DEBUG("failed to parse line:\n  %s", line);

    if (code == COCA_PARSE_DONE) {
      for (conf = cfg; conf != NULL; conf = conf->next) {
        if (!strncmp(conf->name, keyword, COCA_MAX_KEY_LEN)) {
          if (conf->src < prior) {
            if (safe_strcpy(conf->value, value, COCA_MAX_VALUE_LEN))
              return COCA_ERR_STRING;
          }
          else if (conf->src == prior && verbose >= COCA_VERBOSE_WARN)
            COCA_MSG_WARN("duplicate entry: %s\n", keyword);
          conf->src = prior;
          break;
        }
      }
      if (conf == NULL && verbose >= COCA_VERBOSE_WARN)
        COCA_MSG_WARN("unused parameter: %s\n", keyword);
    }
  }
  fclose(fp);
  return 0;
}

int coca_parse_line(const char *line, char *key, mystr value) {
  int i, ik, iv, state;
  char c, quote;

  if (!(line && key && value)) return COCA_PARSE_ERROR;
  ik = iv = 0;
  quote = '\0';         /* for distinguishing single/double quotation mark */
  state = COCA_PARSE_START;

  for (i = 0; i < COCA_MAX_LINE_LEN; i++) {
    c = line[i];
    if (c == '\0' || c == '\n' || c == EOF) {
      key[ik] = value[iv] = '\0';
      if (state == COCA_PARSE_VALUE || state == COCA_PARSE_VALUE_END ||
          state == COCA_PARSE_ARRAY_END)
        return COCA_PARSE_DONE;
      else if (state == COCA_PARSE_QUOTE || state == COCA_PARSE_ARRAY ||
          state == COCA_PARSE_ARRAY_QUOTE)
        return COCA_PARSE_ERROR;
      else if (state == COCA_PARSE_START) return COCA_PARSE_EMPTY;
      else return COCA_PARSE_DEBUG;
    }

    switch (state) {
      case COCA_PARSE_START:
        if (isalpha(c) || c == '_') {
          key[ik++] = c;
          state = COCA_PARSE_KEYWORD;
        }
        else if (c == COCA_SYM_COMMENT) return COCA_PARSE_EMPTY;
        else if (!isblank(c)) return COCA_PARSE_WARN;
        break;
      case COCA_PARSE_KEYWORD:
        if (c == COCA_SYM_EQUAL) {
          key[ik] = '\0';
          state = COCA_PARSE_VALUE_START;
        }
        else if (isblank(c)) {
          key[ik] = '\0';
          state = COCA_PARSE_EQUAL;
        }
        else if (isalnum(c) || c == '_') {
          if (ik >= COCA_MAX_KEY_LEN - 1) return COCA_PARSE_ERROR;
          key[ik++] = c;
        }
        else return COCA_PARSE_WARN;
        break;
      case COCA_PARSE_EQUAL:
        if (c == COCA_SYM_EQUAL) state = COCA_PARSE_VALUE_START;
        else if (!isblank(c)) return COCA_PARSE_WARN;
        break;
      case COCA_PARSE_VALUE_START:
        if (!isblank(c)) {
          if (c == '"' || c == '\'') {
            quote = c;
            state = COCA_PARSE_QUOTE;
          }
          else if (c == COCA_SYM_ARRAY_START) {
            value[iv++] = c;
            state = COCA_PARSE_ARRAY;
          }
          else if (c == COCA_SYM_COMMENT) return COCA_PARSE_DEBUG;
          else if (isprint(c)) {
            value[iv++] = c;
            state = COCA_PARSE_VALUE;
          }
          else return COCA_PARSE_WARN;
        }
        break;
      case COCA_PARSE_VALUE:
        if (isblank(c)) {
          value[iv] = '\0';
          state = COCA_PARSE_VALUE_END;
        }
        else if (isgraph(c)) {
          if (c == COCA_SYM_COMMENT) {
            value[iv] = '\0';
            return COCA_PARSE_DONE;
          }
          else {
            if (iv >= COCA_MAX_VALUE_LEN - 1) return COCA_PARSE_ERROR;
            value[iv++] = c;
          }
        }
        else return COCA_PARSE_WARN;
        break;
      case COCA_PARSE_QUOTE:
        if (c == quote) {
          value[iv] = '\0';
          quote = '\0';
          state = COCA_PARSE_VALUE_END;
        }
        else {
          if (iv >= COCA_MAX_VALUE_LEN - 1) return COCA_PARSE_ERROR;
          value[iv++] = c;
        }
        break;
      case COCA_PARSE_ARRAY:
        if (c == COCA_SYM_COMMENT) return COCA_PARSE_WARN;
        else if (c == '"' || c == '\'') {
          if (iv >= COCA_MAX_VALUE_LEN - 1) return COCA_PARSE_ERROR;
          quote = c;
          value[iv++] = c;
          state = COCA_PARSE_ARRAY_QUOTE;
        }
        else if (c == COCA_SYM_ARRAY_END) {
          if (iv >= COCA_MAX_VALUE_LEN - 1) return COCA_PARSE_ERROR;
          value[iv++] = c;
          value[iv] = '\0';
          state = COCA_PARSE_ARRAY_END;
        }
        else {
          if (iv >= COCA_MAX_VALUE_LEN - 1) return COCA_PARSE_ERROR;
          value[iv++] = c;
        }
        break;
      case COCA_PARSE_VALUE_END:
      case COCA_PARSE_ARRAY_END:
        if (c == COCA_SYM_COMMENT) return COCA_PARSE_DONE;
        else if (isgraph(c)) return COCA_PARSE_ERROR;
        break;
      case COCA_PARSE_ARRAY_QUOTE:
        if (c == quote) {
          quote = '\0';
          state = COCA_PARSE_ARRAY;
        }
        if (iv >= COCA_MAX_VALUE_LEN - 1) return COCA_PARSE_ERROR;
        value[iv++] = c;
        break;
      default:
        return COCA_PARSE_WARN;
    }
  }

  if (i == COCA_MAX_LINE_LEN) return COCA_PARSE_ERROR;
  else return COCA_PARSE_WARN;
}

int coca_parse_array(const mystr value, mystr **array, int *size) {
  int i, iv, state, n;
  char c, quote;

  if (!value) return COCA_PARSE_ERROR;
  quote = '\0';
  state = COCA_PARSE_START;

  /* Count number of elements in the array. */
  for (i = n = 0; i < COCA_MAX_VALUE_LEN; i++) {
    if (state == COCA_PARSE_ARRAY_DONE) break;
    c = value[i];
    if (c == '\0') {
      if (state != COCA_PARSE_ARRAY_END) return COCA_PARSE_ERROR;
      break;
    }

    switch (state) {
      case COCA_PARSE_START:
        if (c == COCA_SYM_ARRAY_START) state = COCA_PARSE_VALUE_START;
        else if (!isblank(c)) {                 /* not an array */
          return COCA_PARSE_ERROR;
        }
        break;
      case COCA_PARSE_VALUE_START:
        if (c == COCA_SYM_ARRAY_SEP) return COCA_PARSE_ERROR;
        else if (c == COCA_SYM_ARRAY_END) return COCA_PARSE_ERROR;
        else if (c == COCA_SYM_COMMENT) return COCA_PARSE_ERROR;
        else if (c == '"' || c == '\'') {
          quote = c;
          state = COCA_PARSE_QUOTE;
        }
        else if (isgraph(c)) state = COCA_PARSE_VALUE;
        else if (!isblank(c)) return COCA_PARSE_ERROR;
        break;
      case COCA_PARSE_VALUE:
        if (isblank(c)) state = COCA_PARSE_VALUE_END;
        else if (c == COCA_SYM_ARRAY_SEP) {   /* new array element */
          n++;
          state = COCA_PARSE_VALUE_START;
        }
        else if (c == COCA_SYM_ARRAY_END) state = COCA_PARSE_ARRAY_END;
        else if (c == COCA_SYM_COMMENT) return COCA_PARSE_ERROR;
        else if (!isgraph(c)) return COCA_PARSE_ERROR;
        break;
      case COCA_PARSE_QUOTE:
        if (c == quote) {
          quote = '\0';
          state = COCA_PARSE_VALUE_END;
        }
        break;
      case COCA_PARSE_VALUE_END:
        if (c == COCA_SYM_ARRAY_SEP) {        /* New array element */
          n++;
          state = COCA_PARSE_VALUE_START;
        }
        else if (c == COCA_SYM_ARRAY_END) state = COCA_PARSE_ARRAY_END;
        else if (!isblank(c)) return COCA_PARSE_ERROR;
        break;
      case COCA_PARSE_ARRAY_END:
        if (c == COCA_SYM_COMMENT) state = COCA_PARSE_ARRAY_DONE;
        else if (isgraph(c)) return COCA_PARSE_ERROR;
        break;
      default:
        return COCA_PARSE_ERROR;
    }
  }

  *size = ++n;
  *array = malloc(sizeof(mystr) * n);
  if (!(*array)) return COCA_ERR_MEMORY;

  n = iv = 0;
  quote = '\0';
  state = COCA_PARSE_START;

  /* Record elements in the array. */
  for (i = 0; i < COCA_MAX_VALUE_LEN; i++) {
    c = value[i];
    if (c == '\0') {
      (*array)[n][iv] = '\0';
      break;
    }

    switch (state) {
      case COCA_PARSE_START:
        if (c == COCA_SYM_ARRAY_START) state = COCA_PARSE_VALUE_START;
        else if (!isblank(c)) return COCA_PARSE_ERROR;
        break;
      case COCA_PARSE_VALUE_START:
        if (c == '"' || c == '\'') {
          quote = c;
          state = COCA_PARSE_QUOTE;
        }
        else if (isgraph(c)) {
          (*array)[n][iv++] = c;
          state = COCA_PARSE_VALUE;
        }
        break;
      case COCA_PARSE_VALUE:
        if (isblank(c)) {
          (*array)[n][iv] = '\0';
          state = COCA_PARSE_VALUE_END;
        }
        else if (c == COCA_SYM_ARRAY_SEP) {
          (*array)[n][iv] = '\0';
          if ((++n) >= *size) return COCA_PARSE_ERROR;
          iv = 0;
          state = COCA_PARSE_VALUE_START;
        }
        else if (c == COCA_SYM_ARRAY_END) {
          (*array)[n][iv] = '\0';
          state = COCA_PARSE_ARRAY_END;
        }
        else if (isgraph(c)) {
          (*array)[n][iv++] = c;
          if (iv >= COCA_MAX_VALUE_LEN) return COCA_PARSE_ERROR;
        }
        break;
      case COCA_PARSE_QUOTE:
        if (c == quote) {
          (*array)[n][iv] = '\0';
          quote = '\0';
          state = COCA_PARSE_VALUE_END;
        }
        else {
          (*array)[n][iv++] = c;
          if (iv >= COCA_MAX_VALUE_LEN) return COCA_PARSE_ERROR;
        }
        break;
      case COCA_PARSE_VALUE_END:
        if (c == COCA_SYM_ARRAY_SEP) {
          (*array)[n][iv] = '\0';
          if ((++n) >= *size) return COCA_PARSE_ERROR;
          iv = 0;
          state = COCA_PARSE_VALUE_START;
        }
        else if (c == COCA_SYM_ARRAY_END) {
          (*array)[n][iv] = '\0';
          state = COCA_PARSE_ARRAY_END;
        }
        break;
      case COCA_PARSE_ARRAY_END:
        if (c == COCA_SYM_COMMENT) {
          if (n != *size - 1) return COCA_ERR_OTHER;
          return COCA_PARSE_DONE;
        }
        break;
      default:
        return COCA_PARSE_ERROR;
    }
  }

  return COCA_PARSE_DONE;
}

int coca_get(const mystr value, const int dtype, void *var) {
  int i, n = 0;

  if (value == NULL || var == NULL) return COCA_ERR_ARG;
  if (value[0] == '\0') return COCA_ERR_NOT_SET;

  switch (dtype) {
    case COCA_DTYPE_BOOL:
      if (!strncmp(value, "1", COCA_MAX_VALUE_LEN) ||
          !strncmp(value, "T", COCA_MAX_VALUE_LEN) ||
          !strncmp(value, "true", COCA_MAX_VALUE_LEN) ||
          !strncmp(value, "TRUE", COCA_MAX_VALUE_LEN) ||
          !strncmp(value, "True", COCA_MAX_VALUE_LEN))
        *((bool *) var) = true;
      else if (!strncmp(value, "0", COCA_MAX_VALUE_LEN) ||
          !strncmp(value, "F", COCA_MAX_VALUE_LEN) ||
          !strncmp(value, "false", COCA_MAX_VALUE_LEN) ||
          !strncmp(value, "FALSE", COCA_MAX_VALUE_LEN) ||
          !strncmp(value, "False", COCA_MAX_VALUE_LEN))
        *((bool *) var) = false;
      else return COCA_ERR_DTYPE;
      break;
    case COCA_DTYPE_CHAR:
      *((char *) var) = value[0];
      if (value[1] != '\0') return COCA_ERR_DTYPE;
      break;
    case COCA_DTYPE_INT:
      if (sscanf(value, "%d%n", (int *) var, &n) != 1)
        return COCA_ERR_DTYPE;
      break;
    case COCA_DTYPE_LONG:
      if (sscanf(value, "%ld%n", (long *) var, &n) != 1)
        return COCA_ERR_DTYPE;
      break;
    case COCA_DTYPE_FLT:
      if (sscanf(value, "%f%n", (float *) var, &n) != 1)
        return COCA_ERR_DTYPE;
      break;
    case COCA_DTYPE_DBL:
      if (sscanf(value, "%lf%n", (double *) var, &n) != 1)
        return COCA_ERR_DTYPE;
      break;
    case COCA_DTYPE_STR:
      if (safe_strcpy(*((mystr *) var), value, COCA_MAX_VALUE_LEN))
        return COCA_ERR_STRING;
      break;
    default:
      return COCA_ERR_ARG;
  }

  if (n != 0) {
    for (i = n; i < COCA_MAX_VALUE_LEN; i++) {
      if (value[i] == '\0') break;
      if (!isspace(value[i])) return COCA_ERR_DTYPE;
    }
  }
  return 0;
}

int coca_get_var(coca *cfg, const char *name, const int dtype,
    void *var) {
  int ecode;
  coca *conf;

  if (cfg == NULL || name == NULL || var == NULL) return COCA_ERR_ARG;
  if (!isalpha(name[0]) && name[0] != '_') return COCA_ERR_ARG;

  for (conf = cfg; conf != NULL; conf = conf->next)
    if (!strncmp(name, conf->name, COCA_MAX_KEY_LEN)) break;
  if (conf == NULL) return COCA_ERR_NOT_FOUND;
  if (conf->value[0] == '\0') {
    if (dtype == COCA_DTYPE_BOOL) conf->value[0] = '1';
    else return COCA_ERR_NOT_SET;
  }

  if ((ecode = coca_get(conf->value, dtype, var))) return ecode;
  return 0;
}

int coca_get_arr(const coca *cfg, const char *name, const int dtype,
    void *arr, int *size) {
  int ecode, i, n, isarray = 0;
  char start;
  const coca *conf;
  mystr *element;

  if (cfg == NULL || name == NULL || arr == NULL) return COCA_ERR_ARG;
  if (!isalpha(name[0]) && name[0] != '_') return COCA_ERR_ARG;
  for (conf = cfg; conf != NULL; conf = conf->next)
    if (!strncmp(name, conf->name, COCA_MAX_KEY_LEN)) break;
  if (conf == NULL) return COCA_ERR_NOT_FOUND;
  if (conf->value[0] == '\0') return COCA_ERR_NOT_SET;

  for (i = 0; i < COCA_MAX_VALUE_LEN; i++) {
    start = conf->value[i];
    if (!isblank(start) || start == '\0') break;
  }
  if (start == COCA_SYM_ARRAY_START) isarray = 1;

  if (isarray) {
    if ((ecode = coca_parse_array(conf->value, &element, &n))) return ecode;
  }
  else {
    n = 1;
    element = malloc(sizeof(mystr));
    if (!element) return COCA_ERR_MEMORY;
    safe_strcpy(element[0], conf->value, COCA_MAX_VALUE_LEN);
  }
  *size = n;

  switch (dtype) {
    case COCA_DTYPE_BOOL:
      *((bool **) arr) = malloc(sizeof(bool) * n);
      if (!(*((bool **) arr))) {
        free(element);
        return COCA_ERR_MEMORY;
      }
      for (i = 0; i < n; i++) {
        if ((ecode = coca_get(element[i], dtype, *((bool **) arr) + i))) {
          free(element);
          return ecode;
        }
      }
      break;
    case COCA_DTYPE_CHAR:
      *((char **) arr) = malloc(sizeof(char) * n);
      if (!(*((char **) arr))) {
        free(element);
        return COCA_ERR_MEMORY;
      }
      for (i = 0; i < n; i++) {
        if ((ecode = coca_get(element[i], dtype, *((char **) arr) + i))) {
          free(element);
          return ecode;
        }
      }
      break;
    case COCA_DTYPE_INT:
      *((int **) arr) = malloc(sizeof(int) * n);
      if (!(*((int **) arr))) {
        free(element);
        return COCA_ERR_MEMORY;
      }
      for (i = 0; i < n; i++) {
        if ((ecode = coca_get(element[i], dtype, *((int **) arr) + i))) {
          free(element);
          return ecode;
        }
      }
      break;
    case COCA_DTYPE_LONG:
      *((long **) arr) = malloc(sizeof(long) * n);
      if (!(*((long **) arr))) {
        free(element);
        return COCA_ERR_MEMORY;
      }
      for (i = 0; i < n; i++) {
        if ((ecode = coca_get(element[i], dtype, *((long **) arr) + i))) {
          free(element);
          return ecode;
        }
      }
      break;
    case COCA_DTYPE_FLT:
      *((float **) arr) = malloc(sizeof(float) * n);
      if (!(*((float **) arr))) {
        free(element);
        return COCA_ERR_MEMORY;
      }
      for (i = 0; i < n; i++) {
        if ((ecode = coca_get(element[i], dtype, *((float **) arr) + i))) {
          free(element);
          return ecode;
        }
      }
      break;
    case COCA_DTYPE_DBL:
      *((double **) arr) = malloc(sizeof(double) * n);
      if (!(*((double **) arr))) {
        free(element);
        return COCA_ERR_MEMORY;
      }
      for (i = 0; i < n; i++) {
        if ((ecode = coca_get(element[i], dtype, *((double **) arr) + i))) {
          free(element);
          return ecode;
        }
      }
      break;
    case COCA_DTYPE_STR:
      *((mystr **) arr) = element;
      break;
    default:
      free(element);
      return COCA_ERR_ARG;
  }

  if (dtype != COCA_DTYPE_STR) free(element);
  return 0;
}

int coca_check_all(const coca *cfg) {
  const coca *conf;
  if (!cfg) return COCA_ERR_INIT;
  for (conf = cfg; conf != NULL; conf = conf->next)
    if (conf->src == COCA_SRC_NULL) return COCA_ERR_NOT_SET;
  return 0;
}

int safe_strcpy(char *dest, const char *src, const size_t num) {
  size_t i = 0;
  if (dest == NULL || src == NULL) return -1;
  if (dest > src && dest <= src + num) return -1;
  if (src > dest && src <= dest + num) return -1;
  while (i < num - 1 && src[i] != '\0') {
    dest[i] = src[i];
    ++i;
  }
  dest[i] = '\0';
  if (src[i] != '\0') return 1;
  return 0;
}

