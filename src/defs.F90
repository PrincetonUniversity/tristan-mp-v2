#define STR_MAX   280
#define END_LINE  char(10)

#define CORR      1.025
#define TINYXYZ   1e-6

#define M_PI      3.141592653589793
#define CC        0.45
#define CCINV     2.222222222222222

#ifdef threeD
#define fldBoundZ           (-NGHOST) : ((this_meshblock%ptr%sz) - 1 + (NGHOST))
#else
#define fldBoundZ           (0) : (0)
#endif
