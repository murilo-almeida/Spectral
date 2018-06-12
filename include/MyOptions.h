#ifndef _MyOptions
#define _MyOptions
// ****************************************************************************
// Para escoamentos de 2 fluidos imisciveis mixed DG
// ****************************************************************************
// Variaveis globais usadas por GeProb, DG_Prob e PhElem
// ****************************************************************************
//#define N_FIELDS 2
//#define N_VAR 2
//#define T_VAR N_VAR
//const int N_FIELDS = 2;
//const int N_VAR = 2;
#define GBNMAP_DG
#define _NEWMAT

#define _MAC_OS
#define HAVE_MPI

//#define ECHO_ON
#define PRINTF_ON

// Parametros para as matrizes locais dos elementos
# define MAXNFIELDS 2
# define MAXMODES  64
# define MAXQ  10
# define MAXNB 202
# define MAXNI 91
# define MAXNN 100
# define MAXNBL MAXNB*MAXNFIELDS
# define MAXNIL MAXNI*MAXNFIELDS

// Parametros para matrizes globais
# define MAXNL 30
# define MAXNG 30

#endif
