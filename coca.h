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

#ifndef _COCA_H_
#define _COCA_H_

#include <stdio.h>
#include <stdbool.h>

/******************************************************************************
  Definitions for data types.
******************************************************************************/
#define COCA_DTYPE_BOOL         0
#define COCA_DTYPE_CHAR         1
#define COCA_DTYPE_INT          2
#define COCA_DTYPE_LONG         3
#define COCA_DTYPE_FLT          4
#define COCA_DTYPE_DBL          5
#define COCA_DTYPE_STR          6

/******************************************************************************
  Definitions for string lengths.
******************************************************************************/
#define COCA_MAX_LINE_LEN       2048
#define COCA_MAX_KEY_LEN        128
#define COCA_MAX_VALUE_LEN      2048

/******************************************************************************
  Definitions for the parser.
******************************************************************************/
#define COCA_SRC_NULL                   (-1)
#define COCA_PARSE_START                1
#define COCA_PARSE_KEYWORD              2
#define COCA_PARSE_EQUAL                3
#define COCA_PARSE_VALUE_START          4
#define COCA_PARSE_VALUE                5
#define COCA_PARSE_VALUE_END            6
#define COCA_PARSE_QUOTE                7
#define COCA_PARSE_ARRAY                8
#define COCA_PARSE_ARRAY_END            9
#define COCA_PARSE_ARRAY_QUOTE          10
#define COCA_PARSE_ARRAY_DONE           11

#define COCA_PARSE_DONE                 0
#define COCA_PARSE_ERROR                1
#define COCA_PARSE_WARN                 2
#define COCA_PARSE_DEBUG                3
#define COCA_PARSE_EMPTY                99

#define COCA_SYM_EQUAL                  '='
#define COCA_SYM_ARRAY_START            '['
#define COCA_SYM_ARRAY_END              ']'
#define COCA_SYM_ARRAY_SEP              ','
#define COCA_SYM_COMMENT                '#'

#define COCA_CMD_FLAG                   '-'
#define COCA_CMD_ASSIGN                 '='
#define COCA_IS_OPT(argv) (                                                   \
  argv[0] == COCA_CMD_FLAG && argv[1] &&                                      \
  ((isalpha(argv[1]) && (!argv[2] || argv[2] == COCA_CMD_ASSIGN)) ||          \
  (argv[1] == COCA_CMD_FLAG && (!argv[2] || (argv[2] && isgraph(argv[2])))))  \
  )

/******************************************************************************
  Definitions of error codes.
******************************************************************************/
#define COCA_ERR_MEMORY         101     /* Failed to allocate memory      */
#define COCA_ERR_FILE           102     /* Failed to read file            */
#define COCA_ERR_STRING         103     /* String length exceeding limits */
#define COCA_ERR_PARSE          104     /* Failed to parse line           */
#define COCA_ERR_INIT           105     /* Variable not initialised       */
#define COCA_ERR_EXIST          106     /* Key already exists             */
#define COCA_ERR_ARG            107     /* Incorrect argument             */
#define COCA_ERR_NOT_FOUND      108     /* Configuration not found        */
#define COCA_ERR_NOT_SET        109     /* Parameter not set              */
#define COCA_ERR_DTYPE          110     /* Invalid data type              */
#define COCA_ERR_OTHER          199     /* Unknown errors                 */

/******************************************************************************
  Definitions for printing messages.
******************************************************************************/
#define COCA_VERBOSE_WARN       1
#define COCA_VERBOSE_DEBUG      2

#define COCA_MSG_ERR(...)               \
  fprintf(stderr, "\x1B[31;1mError:\x1B[0m " __VA_ARGS__)
#define COCA_MSG_WARN(...)              \
  fprintf(stderr, "\x1B[35;1mWarning:\x1B[0m " __VA_ARGS__)
#define COCA_MSG_DEBUG(...)             \
  fprintf(stderr, "\x1B[33;1mDebug:\x1B[0m " __VA_ARGS__)

/******************************************************************************
  Definition of data types and structures.
******************************************************************************/
typedef char mystr[COCA_MAX_VALUE_LEN];

typedef struct coca_struct {
  struct coca_struct *next;             /* next coca node;            */
  int src;                              /* source of the value;         */
  int opt;                              /* short command line option;   */
  char long_opt[COCA_MAX_KEY_LEN];      /* long command line option;    */
  char name[COCA_MAX_KEY_LEN];          /* name of the param;           */
  mystr value;                          /* value of the param.          */
} coca;

typedef struct {
  int opt;                              /* short command line option;   */
  char *long_opt;                       /* long command line option;    */
  void (*func_ptr)(void *);             /* pointer to the function;     */
  void *func_arg;                       /* argument of the function.    */
} coca_funcs;

typedef struct {
  bool called;                          /* true after calling function; */
  int opt;                              /* validated short option;      */
  char *long_opt;                       /* validated long option.       */
} coca_func_valid;

/******************************************************************************
  Definition of functions.
******************************************************************************/

coca *coca_init(void);

void coca_destroy(coca *);

int coca_set_param(coca *, const char, const char*, const char *, const int);

int coca_check_funcs(coca *, const coca_funcs *, const int, coca_func_valid *,
    const int);

int coca_read_opts(coca *, const int, char * const *, const coca_funcs *,
    const int, int *, const int, const int);

int coca_read_file(coca *, const char *, const int, const int);

int coca_parse_line(const char *, char *, mystr);

int coca_parse_array(const mystr, mystr **, int *);

int coca_get(const mystr, const int, void *);

int coca_get_var(coca *, const char *, const int, void *);

int coca_get_arr(const coca *, const char *, const int, void *, int *);

int coca_check_all(const coca *);

int safe_strcpy(char *, const char *, const size_t);

#endif
