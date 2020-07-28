/*******************************************************
**                                                    **
**     Common definitions (macros) for linhalo        **
**     Author: Cheng Zhao <zhaocheng03@gmail.com>     **
**                                                    **
*******************************************************/

#ifndef _DEFINE_H_
#define _DEFINE_H_

/******************************************************************************
  Definitions for configurations.
******************************************************************************/
#define DEFAULT_CONF_FILE "linhalo.conf"      // Default configuration file.
#define DEFAULT_FORCE   "False"

/******************************************************************************
  Definitions of runtime constants
******************************************************************************/
#define CHUNK 10000 //1048576   // Size of block for reading ASCII files.
#define COMMENT '#'     // Comment symbol for input files.
#define TIMEOUT 20      // Maximum allowed number of failed inputs.
#define MAX_LEN_LINE 1024       // Maximum length of output lines.
#ifndef M_PI
#define M_PI    3.14159265358979323846
#endif

/******************************************************************************
  Definitions for the format of outputs.
******************************************************************************/
#define FMT_WARN "\n\x1B[35;1mWarning:\x1B[0m "    // Magenta "Warning".
#define FMT_ERR  "\n\x1B[31;1mError:\x1B[0m "      // Red "Error".
#define FMT_EXIT "\x1B[31;1mExit:\x1B[0m "         // Red "Exit".
#define FMT_DONE "\r\x1B[70C[\x1B[32;1mDONE\x1B[0m]\n"    // Green "DONE".
#define FMT_FAIL "\r\x1B[70C[\x1B[31;1mFAIL\x1B[0m]\n"    // Red "FAIL".
#define FMT_KEY(key)    "\x1B[36;1m" #key "\x1B[0m"       // Cyan keyword.

#define OFMT_REAL       "%.8g"          // Output format of real numbers.

/******************************************************************************
  Definitions for error codes.
******************************************************************************/
#define ERR_ARG         -1
#define ERR_MEM         -2
#define ERR_FILE        -3
#define ERR_INPUT       -4
#define ERR_RANGE       -5
#define ERR_STRING      -6
#define ERR_OTHER       -9

/******************************************************************************
  Definitions for small pieces of codes.
******************************************************************************/
#define P_ERR(...) fprintf(stderr, FMT_ERR __VA_ARGS__)
#define P_WRN(...) fprintf(stderr, FMT_WARN __VA_ARGS__)
#define P_EXT(...) fprintf(stderr, FMT_FAIL FMT_EXIT __VA_ARGS__)

#define CHECK_COCA(ecode, param)        {                               \
  if (ecode) {                                                          \
    P_ERR("failed to process configuration " FMT_KEY(param) ".\n");     \
    return ecode;                                                       \
  }                                                                     \
}

#define MESH_IDX(Ng,i,j,k) (i+j*Ng+(size_t)k*Ng*Ng)
#define SWAP(a,b,t)     {t=a; a=b; b=t;}

/******************************************************************************
  Definitions for data types.
******************************************************************************/

#endif
