#ifndef SimTK_FORTRAN_LAPACK_H_
#define SimTK_FORTRAN_LAPACK_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
 * Authors: Jack Middleton                                                    *
 * Contributors: Michael Sherman, Christopher Bruns                           *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

/*
 * This header file contains const-correct function prototypes for C & C++ 
 * programs calling the legacy (Fortran) interface for BLAS and LAPACK
 * version 3. This header should work for almost any implementation of those
 * routines, whether our SimTKlapack, Apple's Accelerate lapack & blas, Intel's
 * Math Kernel Libraries (MKL), AMD's ACML, or Goto- or Atlas-generated 
 * libraries. This will also work with netlib's slow reference implementation,
 * and should work even using the f2c-translated version of that because
 * (although the code is in C) it still implements the Fortran calling
 * sequence.
 *
 * CAUTION: THIS INTERFACE USES 32 BIT INTEGERS. So even though addresses
 * are 8 bytes for a 64-bit lapack and blas, we're expecting all the integer
 * arguments to remain 4 bytes.
 *
 * Do not confuse this interface with the CBLAS and CLAPACK interfaces which 
 * are C-friendly wrappers around the legacy interface. Here we are dealing 
 * with direct calls to the legacy routines (which are Fortran-like) from C 
 * and C++ programs.
 * 
 * The basic rules for C programs calling Fortran-like routines with the 
 * convention we use (there are others) are:
 * 
 * 1) Function names are in lower case and have an underscore appended to the 
 *    name. For example, if a C program calls LAPACK's ZGEEV routine the call 
 *    would be:
 *       zgeev_(...).
 * 
 * 2) Fortran routines pass scalar arguments by reference. (except for 
 *    character string "length" arguments that are normally hidden from 
 *    FORTRAN programmers) Therefore a C program needs to pass pointers to 
 *    scalar arguments. C++ code can just pass the arguments; they will be 
 *    passed by reference automatically because of the declarations here.
 * 
 * 3) In Fortran 2-D arrays are stored in column major format meaning that
 *    the matrix    A = [ 1.0 2.0 ]
 *                      [ 3.0 4.0 ]
 *    declared as A(2,2) would be stored in memory as 1.0, 3.0, 2.0, 4.0.
 *    While a C 2-D array declared as a[2][2], would be stored in
 *    row-major order as 1.0, 2.0, 3.0, 4.0. Therefore C programs may need to
 *    transpose 2D arrays before calling the Fortran interface. Note that
 *    SimTK Matrix objects are normally stored using the Fortran convention,
 *    so their data is directly compatible with Lapack.
 * 
 * 4) The lengths of character strings need to be passed as additional 
 *    arguments which are added to the end of the parameter list. For example,
 *    LAPACK's ZGEEV  routine has two arguments which are character
 *    strings: JOBVL, JOBVR.
 * 
 *    ZGEEV(JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR,
 *                          WORK, LWORK, RWORK, INFO)
 * 
 *    A C program calling ZGEEV would need to add two additional arguments 
 *    at the end of the parameter list which contain the lengths of JOBVL, JOBVR
 *    arguments: 
 *    char* jobvl = "N";
 *    char* jobvr = "Vectors";
 *    int   len_jobvl = 1;
 *    int   len_jobvr = 7;
 *      .
 *      .
 *      .
 * 
 *    zgeev_(jobvl, jobvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, 
 *           work, &lwork, rwork, &info,
 *           len_jobvl, len_jobvr);
 *           ^^^^^^^^   ^^^^^^^^
 *           additional arguments
 * 
 *    In practice, only the first character is used for any Lapack option so 
 *    the length can always be passed as 1. Since these length arguments are 
 *    at the end, they can have defaults in C++ and are set to 1 below so 
 *    C++ programs do not need to be aware of the length arguments. But calls 
 *    from C will typically end with ",1,1,1)" or whatever.
 */

/*
 * We're going to define some temporary preprocessor macros here
 * (SimTK_C_ and SimTK_Z_) to represent complex types.
 * In C++ these will just be the built-in std::complex types. In C
 * we'll either use a type supplied by the including module, or we'll
 * declare complex types here if none are supplied. We assume the
 * binary representation is the same in all cases:
 * "float real,imag;" or "double real,imag;".
 * We define an assortment of temporary macros for other argument
 * passing situations.
 * We'll undefine these temporary macros at the end of this header.
 *
 * Yes, we know this is ugly with all the macros. Think of it as a "header
 * file generator" rather than a header file and it is more palatable. Our
 * goal is to capture all the argument semantics here and then generate the
 * right behavior.
 */

#ifdef __cplusplus

   /* This is C++, just use the built-in complex types. */
   #include <complex>
   #define SimTK_C_               std::complex<float>
   #define SimTK_Z_               std::complex<double>
   #define SimTK_S_INPUT_(s)      const float& s
   #define SimTK_D_INPUT_(d)      const double& d
   #define SimTK_I_INPUT_(i)      const int& i
   #define SimTK_C_INPUT_(c)      const std::complex<float>& c
   #define SimTK_Z_INPUT_(z)      const std::complex<double>& z

   #define SimTK_S_INOUT_(s)      float& s
   #define SimTK_D_INOUT_(d)      double& d
   #define SimTK_I_INOUT_(i)      int& i
   #define SimTK_C_INOUT_(c)      std::complex<float>& c
   #define SimTK_Z_INOUT_(z)      std::complex<double>& z

   #define SimTK_S_OUTPUT_(s)     float& s
   #define SimTK_D_OUTPUT_(d)     double& d
   #define SimTK_I_OUTPUT_(i)     int& i
   #define SimTK_C_OUTPUT_(c)     std::complex<float>& c
   #define SimTK_Z_OUTPUT_(z)     std::complex<double>& z

   #define SimTK_FDIM_(n)         const int& n        /* a dimension, e.g. N,M,lda */
   #define SimTK_FINC_(x)         const int& inc##x   /* increment, i.e. stride */
   #define SimTK_FOPT_(c)         const char& c       /* an option, passed as a single char */
   #define SimTK_CHAR_OUTPUT_(c)  char& c             /* returns a single char */
   #define SimTK_FLEN_(c)         int c##_len=1       /* dummy length parameter added by Fortran */
   #define SimTK_INFO_            int& info           /* returns error code */
#else

  /*
   * This is C, not C++.
   * Should check for 1999 standard C which has built-in SimTK_C_ type
   * here. For now we allow type override via preprocessor symbol;
   * users of 1999 C should provide these before including this file:
   *
   *   #define SimTK_C_FLOAT_COMPLEX_TYPE_TO_USE  float complex
   *   #define SimTK_C_DOUBLE_COMPLEX_TYPE_TO_USE double complex
   *
   */
  #ifdef SimTK_C_FLOAT_COMPLEX_TYPE_TO_USE
     #define SimTK_C_ SimTK_C_FLOAT_COMPLEX_TYPE_TO_USE
  #else
     typedef struct { float real, imag; } SimTK_float_complex;
     #define SimTK_C_ SimTK_float_complex
  #endif

  #ifdef SimTK_C_DOUBLE_COMPLEX_TYPE_TO_USE
     #define SimTK_Z_ SimTK_C_DOUBLE_COMPLEX_TYPE_TO_USE
  #else
     typedef struct { double real, imag; } SimTK_double_complex;
     #define SimTK_Z_ SimTK_double_complex
  #endif

  #define SimTK_S_INPUT_(s)      const float* s
  #define SimTK_D_INPUT_(d)      const double* d
  #define SimTK_I_INPUT_(i)      const int* i
  #define SimTK_C_INPUT_(c)      const SimTK_C_* c
  #define SimTK_Z_INPUT_(z)      const SimTK_Z_* z

  #define SimTK_S_INOUT_(s)      float* s
  #define SimTK_D_INOUT_(d)      double* d
  #define SimTK_I_INOUT_(i)      int* i
  #define SimTK_C_INOUT_(c)      SimTK_C_* c
  #define SimTK_Z_INOUT_(z)      SimTK_Z_* z

  #define SimTK_S_OUTPUT_(s)     float* s
  #define SimTK_D_OUTPUT_(d)     double* d
  #define SimTK_I_OUTPUT_(i)     int* i
  #define SimTK_C_OUTPUT_(c)     SimTK_C_* c
  #define SimTK_Z_OUTPUT_(z)     SimTK_Z_* z

  #define SimTK_FDIM_(n)         const int* n      /* a dimension, e.g. N,M,lda */
  #define SimTK_FINC_(x)         const int* inc##x /* increment, i.e. stride */
  #define SimTK_FOPT_(c)         const char* c     /* an option, passed as a single char */
  #define SimTK_CHAR_OUTPUT_(c)  char* c           /* returns a single char */
  #define SimTK_FLEN_(c)         int c##_len       /* dummy length parameter (must set to 1 in call) */
  #define SimTK_INFO_            int *info         /* returns error code */

#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * These are the standard routines provided by all SimTK libraries so that 
 * various information about the particulars of the library can be extracted 
 * from the binary.
 */
void SimTK_version_SimTKlapack(int*,int*,int*);
void SimTK_about_SimTKlapack(const char*, int, char*);

/*
 * These signatures define callouts to be made by some of the Lapack eigenvalue
 * routines for selecting eigenvalue subsets.
 */
typedef int (* SimTK_SELECT_2S)(SimTK_S_INPUT_(wr), SimTK_S_INPUT_(wi));
typedef int (* SimTK_SELECT_3F)(SimTK_S_INPUT_(ar), SimTK_S_INPUT_(ai), SimTK_S_INPUT_(b));
typedef int (* SimTK_SELECT_2D)(SimTK_D_INPUT_(wr), SimTK_D_INPUT_(wi));
typedef int (* SimTK_SELECT_3D)(SimTK_D_INPUT_(ar), SimTK_D_INPUT_(ai), SimTK_D_INPUT_(b));
typedef int (* SimTK_SELECT_C) (SimTK_C_INPUT_(w));
typedef int (* SimTK_SELECT_2C)(SimTK_C_INPUT_(a),  SimTK_C_INPUT_(b));
typedef int (* SimTK_SELECT_Z) (SimTK_Z_INPUT_(w));
typedef int (* SimTK_SELECT_2Z)(SimTK_Z_INPUT_(a),  SimTK_Z_INPUT_(b));

/*******************************************************************************
 * The BLAS routines. For documentation, see the LAPACK User's Guide, 3rd ed., *
 * Appendix C "Quick Reference Guide to the BLAS", pg. 180-4.                  *
 *******************************************************************************/

/*
 *  ****************
 *  * BLAS Level 1 *
 *  ****************
 *
 *  BLAS Level 1 functions (that is, value-returning methods).
 *
 *  TODO: The following functions return complex values. This is OK in C++ but
 *  what about C?
 *    cdotu_, zdotu_ (complex dot product without conjugation)
 *    cdotc_, zdotc_ (complex dot product with conjugation)
 * 
 *    SimTKlapack  1.2  was compiled with gfortran and the -ff2c flag on Mac 
 *    and Linux and g77 on Windows (gfortran still had issues on Windows). In 
 *    either case the four routines below return SimTK_C_ or SimTK_Z_ type in 
 *    C++ but for C programs they return a pointer to a SimTK_C_ or SimTK_Z_ 
 *    type as an extra parameter. 
 */
/*
#ifdef __cplusplus
SimTK_C_ cdotu_(SimTK_FDIM_(n), const SimTK_C_ *x, SimTK_FINC_(x), 
                                                const SimTK_C_ *y, SimTK_FINC_(y));
SimTK_Z_ zdotu_(SimTK_FDIM_(n), const SimTK_Z_ *x, SimTK_FINC_(x), 
                                                const SimTK_Z_ *y, SimTK_FINC_(y));
SimTK_C_ cdotc_(SimTK_FDIM_(n), const SimTK_C_ *x, SimTK_FINC_(x), 
                                                const SimTK_C_ *y, SimTK_FINC_(y));
SimTK_Z_ zdotc_(SimTK_FDIM_(n), const SimTK_Z_ *x, SimTK_FINC_(x), 
                                                const SimTK_Z_ *y, SimTK_FINC_(y));
#else
void cdotu_(SimTK_C_*, SimTK_FDIM_(n), const SimTK_C_ *x, SimTK_FINC_(x), 
                                                const SimTK_C_ *y, SimTK_FINC_(y));
void zdotu_(SimTK_Z_*, SimTK_FDIM_(n), const SimTK_Z_ *x, SimTK_FINC_(x), 
                                                const SimTK_Z_ *y, SimTK_FINC_(y));
void cdotc_(SimTK_C_*, SimTK_FDIM_(n), const SimTK_C_ *x, SimTK_FINC_(x), 
                                                const SimTK_C_ *y, SimTK_FINC_(y));
void zdotc_(SimTK_Z_*, SimTK_FDIM_(n), const SimTK_Z_ *x, SimTK_FINC_(x), 
                                                const SimTK_Z_ *y, SimTK_FINC_(y));
#endif
*/

float  sdot_  (SimTK_FDIM_(n), const float  *x, SimTK_FINC_(x), const float  *y, SimTK_FINC_(y));
double ddot_  (SimTK_FDIM_(n), const double *x, SimTK_FINC_(x), const double *y, SimTK_FINC_(y));
double dsdot_ (SimTK_FDIM_(n), const float  *x, SimTK_FINC_(x), const float  *y, SimTK_FINC_(y));
float  sdsdot_(SimTK_FDIM_(n), SimTK_S_INPUT_(alpha), 
                                      const float  *x, SimTK_FINC_(x), const float  *y, SimTK_FINC_(y));

/* Functions having prefixes S D SC DZ */
float snrm2_(SimTK_FDIM_(n), const float *x, SimTK_FINC_(x));
float sasum_(SimTK_FDIM_(n), const float *x, SimTK_FINC_(x));

double dnrm2_(SimTK_FDIM_(n), const double *x, SimTK_FINC_(x));
double dasum_(SimTK_FDIM_(n), const double *x, SimTK_FINC_(x));

float scnrm2_(SimTK_FDIM_(n), const SimTK_C_ *x, SimTK_FINC_(x));
float scasum_(SimTK_FDIM_(n), const SimTK_C_ *x, SimTK_FINC_(x));

double dznrm2_(SimTK_FDIM_(n), const SimTK_Z_ *x, SimTK_FINC_(x));
double dzasum_(SimTK_FDIM_(n), const SimTK_Z_ *x, SimTK_FINC_(x));

/* Int functions having standard 4 prefixes I(S D C Z) */
int isamax_(SimTK_FDIM_(n), const float    *x, SimTK_FINC_(x));
int idamax_(SimTK_FDIM_(n), const double   *x, SimTK_FINC_(x));
int icamax_(SimTK_FDIM_(n), const SimTK_C_ *x, SimTK_FINC_(x));
int izamax_(SimTK_FDIM_(n), const SimTK_Z_ *x, SimTK_FINC_(x));

/* BLAS Level 1 subroutines (that is, void methods). */

/* Routines with standard 4 prefixes _(s, d, c, z) */
void sswap_(SimTK_FDIM_(n), float    *x, SimTK_FINC_(x), float    *y, SimTK_FINC_(y));
void dswap_(SimTK_FDIM_(n), double   *x, SimTK_FINC_(x), double   *y, SimTK_FINC_(y));
void cswap_(SimTK_FDIM_(n), SimTK_C_ *x, SimTK_FINC_(x), SimTK_C_ *y, SimTK_FINC_(y));
void zswap_(SimTK_FDIM_(n), SimTK_Z_ *x, SimTK_FINC_(x), SimTK_Z_ *y, SimTK_FINC_(y));

/* assign y = x */
void scopy_(SimTK_FDIM_(n), const float    *x, SimTK_FINC_(x), float    *y, SimTK_FINC_(y));
void dcopy_(SimTK_FDIM_(n), const double   *x, SimTK_FINC_(x), double   *y, SimTK_FINC_(y));
void ccopy_(SimTK_FDIM_(n), const SimTK_C_ *x, SimTK_FINC_(x), SimTK_C_ *y, SimTK_FINC_(y));
void zcopy_(SimTK_FDIM_(n), const SimTK_Z_ *x, SimTK_FINC_(x), SimTK_Z_ *y, SimTK_FINC_(y));

/* y += ax */
void saxpy_(SimTK_FDIM_(n), SimTK_S_INPUT_(alpha), const float    *x, SimTK_FINC_(x), float    *y, SimTK_FINC_(y));
void daxpy_(SimTK_FDIM_(n), SimTK_D_INPUT_(alpha), const double   *x, SimTK_FINC_(x), double   *y, SimTK_FINC_(y));
void caxpy_(SimTK_FDIM_(n), SimTK_C_INPUT_(alpha), const SimTK_C_ *x, SimTK_FINC_(x), SimTK_C_ *y, SimTK_FINC_(y));
void zaxpy_(SimTK_FDIM_(n), SimTK_Z_INPUT_(alpha), const SimTK_Z_ *x, SimTK_FINC_(x), SimTK_Z_ *y, SimTK_FINC_(y));


/*
 * Routines with S and D prefix only
 */

/* a,b are in/out, c,s are output (all scalars) */
void srotg_(SimTK_S_OUTPUT_(a), SimTK_S_OUTPUT_(b), SimTK_S_OUTPUT_(c), SimTK_S_OUTPUT_(s));
void drotg_(SimTK_D_OUTPUT_(a), SimTK_D_OUTPUT_(b), SimTK_D_OUTPUT_(c), SimTK_D_OUTPUT_(s));

/* all parameters are in/out */
void srotmg_(SimTK_S_OUTPUT_(d1), SimTK_S_OUTPUT_(d2), SimTK_S_OUTPUT_(b1), SimTK_S_OUTPUT_(b2), float  P[5]);
void drotmg_(SimTK_D_OUTPUT_(d1), SimTK_D_OUTPUT_(d2), SimTK_D_OUTPUT_(b1), SimTK_D_OUTPUT_(b2), double P[5]);

void srot_(SimTK_FDIM_(n), float  *x, SimTK_FINC_(x), float  *y, SimTK_FINC_(y), SimTK_S_INPUT_(c), SimTK_S_INPUT_(s));
void drot_(SimTK_FDIM_(n), double *x, SimTK_FINC_(x), double *y, SimTK_FINC_(y), SimTK_D_INPUT_(c), SimTK_D_INPUT_(s));

void srotm_(SimTK_FDIM_(n), float  *x, SimTK_FINC_(x), float  *y, SimTK_FINC_(y), const float  P[5]);
void drotm_(SimTK_FDIM_(n), double *x, SimTK_FINC_(x), double *y, SimTK_FINC_(y), const double P[5]);


/*
 * Routines with S D C Z CS and ZD prefixes
 */
void sscal_ (SimTK_FDIM_(n), SimTK_S_INPUT_(alpha), float    *x, SimTK_FINC_(x));
void dscal_ (SimTK_FDIM_(n), SimTK_D_INPUT_(alpha), double   *x, SimTK_FINC_(x));
void cscal_ (SimTK_FDIM_(n), SimTK_C_INPUT_(alpha), SimTK_C_ *x, SimTK_FINC_(x));
void zscal_ (SimTK_FDIM_(n), SimTK_Z_INPUT_(alpha), SimTK_Z_ *x, SimTK_FINC_(x));
void csscal_(SimTK_FDIM_(n), SimTK_S_INPUT_(alpha), SimTK_C_ *x, SimTK_FINC_(x));
void zdscal_(SimTK_FDIM_(n), SimTK_D_INPUT_(alpha), SimTK_Z_ *x, SimTK_FINC_(x));

/*
 * Extra reference routines provided by ATLAS, but not mandated by the 
 * standard.
 */
void crotg_(SimTK_C_OUTPUT_(a), SimTK_C_INPUT_(b), SimTK_S_OUTPUT_(c), SimTK_C_OUTPUT_(s) );
void zrotg_(SimTK_Z_OUTPUT_(a), SimTK_Z_INPUT_(b), SimTK_D_OUTPUT_(c), SimTK_Z_OUTPUT_(s) );
void csrot_(SimTK_FDIM_(n), SimTK_C_ *x, SimTK_FINC_(x), SimTK_C_ *y, SimTK_FINC_(y),
                   SimTK_S_INPUT_(c),  SimTK_S_INPUT_(s));
void zdrot_(SimTK_FDIM_(n), SimTK_Z_ *x, SimTK_FINC_(x), SimTK_Z_ *y, SimTK_FINC_(y),
                   SimTK_D_INPUT_(c), SimTK_D_INPUT_(s));

/*
 *===========================================================================
 * Prototypes for level 2 BLAS
 *===========================================================================
 */

/*
 *Routines with standard 4 prefixes _(S, D, C, Z)
 */

/* y = alpha A x + beta y */
void sgemv_(SimTK_FOPT_(transA), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_S_INPUT_(alpha), const float *A, SimTK_FDIM_(lda), const float *x, SimTK_FINC_(x), SimTK_S_INPUT_(beta), float *Y, SimTK_FINC_(Y), SimTK_FLEN_(transA));
void sgbmv_(SimTK_FOPT_(transA), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_S_INPUT_(alpha), const float *A, SimTK_FDIM_(lda), const float *x, SimTK_FINC_(x), SimTK_S_INPUT_(beta), float *Y, SimTK_FINC_(Y), SimTK_FLEN_(transA));
void strmv_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(n), const float *A, SimTK_FDIM_(lda), float *x, SimTK_FINC_(x), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));
void stbmv_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(k), const float *A, SimTK_FDIM_(lda), float *x, SimTK_FINC_(x), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));
void stpmv_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(n), const float *Ap, float *x, SimTK_FINC_(x), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));
void strsv_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(n), const float *A, SimTK_FDIM_(lda), float *x, SimTK_FINC_(x), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));
void stbsv_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(k), const float *A, SimTK_FDIM_(lda), float *x, SimTK_FINC_(x), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));
void stpsv_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(n), const float *Ap, float *x, SimTK_FINC_(x), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));

void dgemv_(SimTK_FOPT_(transA), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_D_INPUT_(alpha), const double *A, SimTK_FDIM_(lda), const double *x, SimTK_FINC_(x), SimTK_D_INPUT_(beta), double *Y, SimTK_FINC_(Y), SimTK_FLEN_(transA));
void dgbmv_(SimTK_FOPT_(transA), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_D_INPUT_(alpha), const double *A, SimTK_FDIM_(lda), const double *x, SimTK_FINC_(x), SimTK_D_INPUT_(beta), double *Y, SimTK_FINC_(Y), SimTK_FLEN_(transA));
void dtrmv_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(n), const double *A, SimTK_FDIM_(lda), double *x, SimTK_FINC_(x), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));
void dtbmv_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(k), const double *A, SimTK_FDIM_(lda), double *x, SimTK_FINC_(x), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));
void dtpmv_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(n), const double *Ap, double *x, SimTK_FINC_(x), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));
void dtrsv_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(n), const double *A, SimTK_FDIM_(lda), double *x, SimTK_FINC_(x), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));
void dtbsv_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(k), const double *A, SimTK_FDIM_(lda), double *x, SimTK_FINC_(x), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));
void dtpsv_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(n), const double *Ap, double *x, SimTK_FINC_(x), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));

void cgemv_(SimTK_FOPT_(transA), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_C_INPUT_(alpha), const SimTK_C_ *A, SimTK_FDIM_(lda), const SimTK_C_ *x, SimTK_FINC_(x), SimTK_C_INPUT_(beta), SimTK_C_ *Y, SimTK_FINC_(Y), SimTK_FLEN_(transA));
void cgbmv_(SimTK_FOPT_(transA), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_C_INPUT_(alpha), const SimTK_C_ *A, SimTK_FDIM_(lda), const SimTK_C_ *x, SimTK_FINC_(x), SimTK_C_INPUT_(beta), SimTK_C_ *Y, SimTK_FINC_(Y), SimTK_FLEN_(transA));
void ctrmv_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(n), const SimTK_C_ *A, SimTK_FDIM_(lda), SimTK_C_*x, SimTK_FINC_(x), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));
void ctbmv_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_C_ *A, SimTK_FDIM_(lda), SimTK_C_*x, SimTK_FINC_(x), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));
void ctpmv_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(n), const SimTK_C_ *Ap, SimTK_C_*x, SimTK_FINC_(x), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));
void ctrsv_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(n), const SimTK_C_ *A, SimTK_FDIM_(lda), SimTK_C_*x, SimTK_FINC_(x), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));
void ctbsv_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_C_ *A, SimTK_FDIM_(lda), SimTK_C_*x, SimTK_FINC_(x), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));
void ctpsv_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(n), const SimTK_C_ *Ap, SimTK_C_*x, SimTK_FINC_(x), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));

void zgemv_(SimTK_FOPT_(transA), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_Z_INPUT_(alpha), const SimTK_Z_ *A, SimTK_FDIM_(lda), const SimTK_Z_ *x, SimTK_FINC_(x), SimTK_Z_INPUT_(beta), SimTK_Z_ *Y, SimTK_FINC_(Y), SimTK_FLEN_(transA));
void zgbmv_(SimTK_FOPT_(transA), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_Z_INPUT_(alpha), const SimTK_Z_ *A, SimTK_FDIM_(lda), const SimTK_Z_ *x, SimTK_FINC_(x), SimTK_Z_INPUT_(beta), SimTK_Z_ *Y, SimTK_FINC_(Y),
                  SimTK_FLEN_(transA));
void ztrmv_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(n), const SimTK_Z_ *A, SimTK_FDIM_(lda), SimTK_Z_*x, SimTK_FINC_(x), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));
void ztbmv_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_Z_ *A, SimTK_FDIM_(lda), SimTK_Z_*x, SimTK_FINC_(x), SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));
void ztpmv_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(n), const SimTK_Z_ *Ap, SimTK_Z_*x, SimTK_FINC_(x), SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));
void ztrsv_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(n), const SimTK_Z_ *A, SimTK_FDIM_(lda), SimTK_Z_*x, SimTK_FINC_(x), SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));
void ztbsv_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_Z_ *A, SimTK_FDIM_(lda), SimTK_Z_*x, SimTK_FINC_(x), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));
void ztpsv_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(n), const SimTK_Z_ *Ap, SimTK_Z_*x, SimTK_FINC_(x), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));

/*
 * Routines with S and D prefixes only
 */
/* y = alpha A x + beta y */
void ssymv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n),  SimTK_S_INPUT_(alpha), const float *A, SimTK_FDIM_(lda), const float *x, SimTK_FINC_(x),  SimTK_S_INPUT_(beta), float *y, SimTK_FINC_(y), SimTK_FLEN_(uplo));
void ssbmv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(k),  SimTK_S_INPUT_(alpha), const float *A, SimTK_FDIM_(lda), const float *x, SimTK_FINC_(x), SimTK_S_INPUT_(beta),  float *y, SimTK_FINC_(y), SimTK_FLEN_(uplo));
void sspmv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n),  SimTK_S_INPUT_(alpha), const float *Ap, const float *x, SimTK_FINC_(x),  SimTK_S_INPUT_(beta), float *y, SimTK_FINC_(y), SimTK_FLEN_(uplo));

void dsymv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_D_INPUT_(alpha), const double *A, SimTK_FDIM_(lda), const double *x, SimTK_FINC_(x), SimTK_D_INPUT_(beta), double *y, SimTK_FINC_(y), SimTK_FLEN_(uplo));
void dsbmv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_D_INPUT_(alpha), const double *A, SimTK_FDIM_(lda), const double *x, SimTK_FINC_(x), SimTK_D_INPUT_(beta), double *y, SimTK_FINC_(y), SimTK_FLEN_(uplo));
void dspmv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_D_INPUT_(alpha), const double *Ap, const double *x, SimTK_FINC_(x), SimTK_D_INPUT_(beta), double *y, SimTK_FINC_(y), SimTK_FLEN_(uplo));

/* x,y are const, A is in/out */
void sger_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_S_INPUT_(alpha), const float *x, SimTK_FINC_(x), const float *y, SimTK_FINC_(y), float *A, SimTK_FDIM_(lda));
void ssyr_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_S_INPUT_(alpha), const float *x, SimTK_FINC_(x), float *A, SimTK_FDIM_(lda), SimTK_FLEN_(uplo));
void sspr_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_S_INPUT_(alpha), const float *x, SimTK_FINC_(x), float *Ap, SimTK_FLEN_(uplo));
void ssyr2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_S_INPUT_(alpha), const float *x, SimTK_FINC_(x), const float *y, SimTK_FINC_(y), float *A, SimTK_FDIM_(lda), SimTK_FLEN_(uplo));
void sspr2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_S_INPUT_(alpha), const float *x, SimTK_FINC_(x), const float *y, SimTK_FINC_(y), float *A, SimTK_FLEN_(uplo));

void dger_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_D_INPUT_(alpha), const double *x, SimTK_FINC_(x), const double *y, SimTK_FINC_(y), double *A, SimTK_FDIM_(lda));
void dsyr_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_D_INPUT_(alpha), const double *x, SimTK_FINC_(x), double *A, SimTK_FDIM_(lda), SimTK_FLEN_(uplo));
void dspr_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_D_INPUT_(alpha), const double *x, SimTK_FINC_(x), double *Ap, SimTK_FLEN_(uplo));
void dsyr2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_D_INPUT_(alpha), const double *x, SimTK_FINC_(x), const double *y, SimTK_FINC_(y), double *A, SimTK_FDIM_(lda), SimTK_FLEN_(uplo));
void dspr2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_D_INPUT_(alpha), const double *x, SimTK_FINC_(x), const double *y, SimTK_FINC_(y), double *A, SimTK_FLEN_(uplo));

/*
 * Routines with C and Z prefixes only
 */
/* y = alpha A x + beta y */
void chemv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_INPUT_(alpha), const SimTK_C_ *A, SimTK_FDIM_(lda), const SimTK_C_ *x, SimTK_FINC_(x), SimTK_C_INPUT_(beta), SimTK_C_ *y, SimTK_FINC_(y), SimTK_FLEN_(uplo));
void chbmv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_C_INPUT_(alpha), const SimTK_C_ *A, SimTK_FDIM_(lda), const SimTK_C_ *x, SimTK_FINC_(x), SimTK_C_INPUT_(beta), SimTK_C_ *y, SimTK_FINC_(y), SimTK_FLEN_(uplo));
void chpmv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_INPUT_(alpha), const SimTK_C_ *Ap, const SimTK_C_ *x, SimTK_FINC_(x), SimTK_C_INPUT_(beta), SimTK_C_ *y, SimTK_FINC_(y), SimTK_FLEN_(uplo));

void zhemv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_INPUT_(alpha), const SimTK_Z_ *A, SimTK_FDIM_(lda), const SimTK_Z_ *x, SimTK_FINC_(x), SimTK_Z_INPUT_(beta), SimTK_Z_ *y, SimTK_FINC_(y), SimTK_FLEN_(uplo));
void zhbmv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_Z_INPUT_(alpha), const SimTK_Z_ *A, SimTK_FDIM_(lda), const SimTK_Z_ *x, SimTK_FINC_(x), SimTK_Z_INPUT_(beta), SimTK_Z_ *y, SimTK_FINC_(y), SimTK_FLEN_(uplo));
void zhpmv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_INPUT_(alpha), const SimTK_Z_ *Ap, const SimTK_Z_ *x, SimTK_FINC_(x), SimTK_Z_INPUT_(beta), SimTK_Z_ *y, SimTK_FINC_(y), SimTK_FLEN_(uplo));

/* x,y are const, A is in/out */
void cgeru_(SimTK_FDIM_(m),    SimTK_FDIM_(n), SimTK_C_INPUT_(alpha), const SimTK_C_ *x, SimTK_FINC_(x), const SimTK_C_ *y, SimTK_FINC_(y), SimTK_C_ *A, SimTK_FDIM_(lda));
void cgerc_(SimTK_FDIM_(m),    SimTK_FDIM_(n), SimTK_C_INPUT_(alpha), const SimTK_C_ *x, SimTK_FINC_(x), const SimTK_C_ *y, SimTK_FINC_(y), SimTK_C_ *A, SimTK_FDIM_(lda));
void cher_ (SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_S_INPUT_(alpha),    const SimTK_C_ *x, SimTK_FINC_(x), SimTK_C_ *A, SimTK_FDIM_(lda), SimTK_FLEN_(uplo));
void chpr_ (SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_S_INPUT_(alpha),    const SimTK_C_ *x, SimTK_FINC_(x), SimTK_C_ *A, SimTK_FLEN_(uplo));
void cher2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_INPUT_(alpha), const SimTK_C_ *x, SimTK_FINC_(x), const SimTK_C_ *y, SimTK_FINC_(y), SimTK_C_ *A, SimTK_FDIM_(lda), SimTK_FLEN_(uplo));
void chpr2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_INPUT_(alpha), const SimTK_C_ *x, SimTK_FINC_(x), const SimTK_C_ *y, SimTK_FINC_(y), SimTK_C_ *Ap, SimTK_FLEN_(uplo));

void zgeru_(SimTK_FDIM_(m),    SimTK_FDIM_(n), SimTK_Z_INPUT_(alpha), const SimTK_Z_ *x, SimTK_FINC_(x), const SimTK_Z_ *y, SimTK_FINC_(y), SimTK_Z_ *A, SimTK_FDIM_(lda));
void zgerc_(SimTK_FDIM_(m),    SimTK_FDIM_(n), SimTK_Z_INPUT_(alpha), const SimTK_Z_ *x, SimTK_FINC_(x), const SimTK_Z_ *y, SimTK_FINC_(y), SimTK_Z_ *A, SimTK_FDIM_(lda));
void zher_ (SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_D_INPUT_(alpha),   const SimTK_Z_ *x, SimTK_FINC_(x), SimTK_Z_ *A, SimTK_FDIM_(lda), SimTK_FLEN_(uplo));
void zhpr_ (SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_D_INPUT_(alpha),   const SimTK_Z_ *x, SimTK_FINC_(x), SimTK_Z_ *A, SimTK_FLEN_(uplo));
void zher2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_INPUT_(alpha), const SimTK_Z_ *x, SimTK_FINC_(x), const SimTK_Z_ *y, SimTK_FINC_(y), SimTK_Z_ *A, SimTK_FDIM_(lda), SimTK_FLEN_(uplo));
void zhpr2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_INPUT_(alpha), const SimTK_Z_ *x, SimTK_FINC_(x), const SimTK_Z_ *y, SimTK_FINC_(y), SimTK_Z_ *Ap, SimTK_FLEN_(uplo));

/*
 *===========================================================================
 * Prototypes for level 3 BLAS
 *===========================================================================
 */

/*
 * Routines with standard 4 prefixes _(S, D, C, Z)
 */
/* A, B are input, C in/out for gemm, symm, syrk, syr2k */
void sgemm_ (SimTK_FOPT_(transA), SimTK_FOPT_(transB), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_S_INPUT_(alpha), const float *A, SimTK_FDIM_(lda), const float *B, SimTK_FDIM_(ldb), SimTK_S_INPUT_(beta), float *C, SimTK_FDIM_(ldc), SimTK_FLEN_(transA), SimTK_FLEN_(transB));
void ssymm_ (SimTK_FOPT_(side), SimTK_FOPT_(uplo),   SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_S_INPUT_(alpha), const float *A, SimTK_FDIM_(lda), const float *B, SimTK_FDIM_(ldb), SimTK_S_INPUT_(beta), float *C, SimTK_FDIM_(ldc), SimTK_FLEN_(side), SimTK_FLEN_(uplo));
void ssyrk_ (SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_S_INPUT_(alpha), const float *A, SimTK_FDIM_(lda), SimTK_S_INPUT_(beta), float *C, SimTK_FDIM_(ldc), SimTK_FLEN_(uplo), SimTK_FLEN_(transA));
void ssyr2k_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_S_INPUT_(alpha), const float *A, SimTK_FDIM_(lda), const float *B, SimTK_FDIM_(ldb), SimTK_S_INPUT_(beta), float *C, SimTK_FDIM_(ldc), SimTK_FLEN_(uplo), SimTK_FLEN_(transA));

void dgemm_ (SimTK_FOPT_(transA), SimTK_FOPT_(transB), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_D_INPUT_(alpha), const double *A, SimTK_FDIM_(lda),
                    const double *B, SimTK_FDIM_(ldb), SimTK_D_INPUT_(beta), double *C, SimTK_FDIM_(ldc), SimTK_FLEN_(transA), SimTK_FLEN_(transB));
void dsymm_ (SimTK_FOPT_(side), SimTK_FOPT_(uplo),   SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_D_INPUT_(alpha), const double *A, SimTK_FDIM_(lda), const double *B, SimTK_FDIM_(ldb), SimTK_D_INPUT_(beta), double *C, SimTK_FDIM_(ldc), SimTK_FLEN_(side), SimTK_FLEN_(uplo));
void dsyrk_ (SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_D_INPUT_(alpha), const double *A, SimTK_FDIM_(lda), SimTK_D_INPUT_(beta), double *C, SimTK_FDIM_(ldc), SimTK_FLEN_(uplo), SimTK_FLEN_(transA));
void dsyr2k_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_D_INPUT_(alpha), const double *A, SimTK_FDIM_(lda), const double *B, SimTK_FDIM_(ldb), SimTK_D_INPUT_(beta), double *C, SimTK_FDIM_(ldc), SimTK_FLEN_(uplo), SimTK_FLEN_(transA));

void cgemm_ (SimTK_FOPT_(transA), SimTK_FOPT_(transB), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_C_INPUT_(alpha), const SimTK_C_ *A, SimTK_FDIM_(lda),
                    const SimTK_C_ *B, SimTK_FDIM_(ldb), SimTK_C_INPUT_(beta), SimTK_C_ *C, SimTK_FDIM_(ldc), SimTK_FLEN_(transA), SimTK_FLEN_(transB));
void csymm_ (SimTK_FOPT_(side), SimTK_FOPT_(uplo),   SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_C_INPUT_(alpha), const SimTK_C_ *A, SimTK_FDIM_(lda), const SimTK_C_ *B, SimTK_FDIM_(ldb), SimTK_C_INPUT_(beta), SimTK_C_ *C, SimTK_FDIM_(ldc), SimTK_FLEN_(side), SimTK_FLEN_(uplo));
void csyrk_ (SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_C_INPUT_(alpha), const SimTK_C_ *A, SimTK_FDIM_(lda), SimTK_C_INPUT_(beta), SimTK_C_ *C, SimTK_FDIM_(ldc), SimTK_FLEN_(uplo), SimTK_FLEN_(transA));
void csyr2k_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_C_INPUT_(alpha), const SimTK_C_ *A, SimTK_FDIM_(lda), const SimTK_C_ *B, SimTK_FDIM_(ldb), SimTK_C_INPUT_(beta), SimTK_C_ *C, SimTK_FDIM_(ldc), SimTK_FLEN_(uplo), SimTK_FLEN_(transA));

void zgemm_ (SimTK_FOPT_(transA), SimTK_FOPT_(transB), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_Z_INPUT_(alpha), const SimTK_Z_ *A, SimTK_FDIM_(lda),
                    const SimTK_Z_ *B, SimTK_FDIM_(ldb), SimTK_Z_INPUT_(beta), SimTK_Z_ *C, SimTK_FDIM_(ldc), SimTK_FLEN_(transA), SimTK_FLEN_(transB));
void zsymm_ (SimTK_FOPT_(side), SimTK_FOPT_(uplo),   SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_Z_INPUT_(alpha), const SimTK_Z_ *A, SimTK_FDIM_(lda), const SimTK_Z_ *B, SimTK_FDIM_(ldb), SimTK_Z_INPUT_(beta), SimTK_Z_ *C, SimTK_FDIM_(ldc), SimTK_FLEN_(side), SimTK_FLEN_(uplo));
void zsyrk_ (SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_Z_INPUT_(alpha), const SimTK_Z_ *A, SimTK_FDIM_(lda), SimTK_Z_INPUT_(beta), SimTK_Z_ *C, SimTK_FDIM_(ldc), SimTK_FLEN_(uplo), SimTK_FLEN_(trans));
void zsyr2k_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_Z_INPUT_(alpha), const SimTK_Z_ *A, SimTK_FDIM_(lda), const SimTK_Z_ *B, SimTK_FDIM_(ldb), SimTK_Z_INPUT_(beta), SimTK_Z_ *C, SimTK_FDIM_(ldc), SimTK_FLEN_(uplo), SimTK_FLEN_(trans));

/* A is input, B in/out for trmm and trsm */
void strmm_(SimTK_FOPT_(side), SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_S_INPUT_(alpha), const float *A, SimTK_FDIM_(lda), float *B, SimTK_FDIM_(ldb), SimTK_FLEN_(side), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));
void strsm_(SimTK_FOPT_(side), SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_S_INPUT_(alpha), const float *A, SimTK_FDIM_(lda), float *B, SimTK_FDIM_(ldb), SimTK_FLEN_(side), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));

void dtrmm_(SimTK_FOPT_(side), SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_D_INPUT_(alpha), const double *A, SimTK_FDIM_(lda), double *B, SimTK_FDIM_(ldb), SimTK_FLEN_(side), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));
void dtrsm_(SimTK_FOPT_(side), SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_D_INPUT_(alpha), const double *A, SimTK_FDIM_(lda), double *B, SimTK_FDIM_(ldb), SimTK_FLEN_(side), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));

void ctrmm_(SimTK_FOPT_(side), SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_C_INPUT_(alpha), const SimTK_C_ *A, SimTK_FDIM_(lda), SimTK_C_ *B, SimTK_FDIM_(ldb), SimTK_FLEN_(side), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));
void ctrsm_(SimTK_FOPT_(side), SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_C_INPUT_(alpha), const SimTK_C_ *A, SimTK_FDIM_(lda), SimTK_C_ *B, SimTK_FDIM_(ldb), SimTK_FLEN_(side), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));

void ztrmm_(SimTK_FOPT_(side), SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_Z_INPUT_(alpha), const SimTK_Z_ *A, SimTK_FDIM_(lda), SimTK_Z_ *B, SimTK_FDIM_(ldb), SimTK_FLEN_(side), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));
void ztrsm_(SimTK_FOPT_(side), SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FOPT_(diag), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_Z_INPUT_(alpha), const SimTK_Z_ *A, SimTK_FDIM_(lda), SimTK_Z_ *B, SimTK_FDIM_(ldb), SimTK_FLEN_(side), SimTK_FLEN_(uplo), SimTK_FLEN_(transA), SimTK_FLEN_(diag));

/*
 * Routines with prefixes C and Z only
 */
/* A, B are input, C in/out for hemm, herk, her2k */
void chemm_ (SimTK_FOPT_(side), SimTK_FOPT_(uplo),   SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_C_INPUT_(alpha), const SimTK_C_ *A, SimTK_FDIM_(lda), const SimTK_C_ *B, SimTK_FDIM_(ldb), SimTK_C_INPUT_(beta), SimTK_C_ *C, SimTK_FDIM_(ldc), SimTK_FLEN_(side), SimTK_FLEN_(uplo));
void cherk_ (SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_S_INPUT_(alpha),    const SimTK_C_ *A, SimTK_FDIM_(lda), SimTK_S_INPUT_(beta), SimTK_C_ *C, SimTK_FDIM_(ldc), SimTK_FLEN_(uplo), SimTK_FLEN_(trans));
void cher2k_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_C_INPUT_(alpha), const SimTK_C_ *A, SimTK_FDIM_(lda), const SimTK_C_ *B, SimTK_FDIM_(ldb), SimTK_S_INPUT_(beta), SimTK_C_ *C, SimTK_FDIM_(ldc), SimTK_FLEN_(uplo), SimTK_FLEN_(trans));

void zhemm_ (SimTK_FOPT_(side), SimTK_FOPT_(uplo),   SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_Z_INPUT_(alpha), const SimTK_Z_ *A, SimTK_FDIM_(lda), const SimTK_Z_ *B, SimTK_FDIM_(ldb), SimTK_Z_INPUT_(beta), SimTK_Z_ *C, SimTK_FDIM_(ldc), SimTK_FLEN_(side), SimTK_FLEN_(uplo));
void zherk_ (SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_D_INPUT_(alpha),   const SimTK_Z_ *A, SimTK_FDIM_(lda), SimTK_D_INPUT_(beta), SimTK_Z_ *C, SimTK_FDIM_(ldc), SimTK_FLEN_(uplo), SimTK_FLEN_(trans));
void zher2k_(SimTK_FOPT_(uplo), SimTK_FOPT_(transA), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_Z_INPUT_(alpha), const SimTK_Z_ *A, SimTK_FDIM_(lda), const SimTK_Z_ *B, SimTK_FDIM_(ldb), SimTK_D_INPUT_(beta), SimTK_Z_ *C, SimTK_FDIM_(ldc), SimTK_FLEN_(uplo), SimTK_FLEN_(trans));

/* END OF BLAS ROUTINES */






/********************************************************************************
 * The LAPACK routines. For documentation, see the LAPACK User's Guide, 3rd ed. *
 ********************************************************************************/

void cbdsqr_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(ncvt), SimTK_FDIM_(nru), SimTK_FDIM_(ncc), float *d__, float *e, SimTK_C_ *vt, SimTK_FDIM_(ldvt), SimTK_C_ *u, SimTK_FDIM_(ldu), SimTK_C_ *c__, SimTK_FDIM_(ldc), float *rwork, SimTK_INFO_, SimTK_FLEN_(uplo));
void cgbbrd_(SimTK_FOPT_(vect), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(ncc), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_C_ *ab, SimTK_FDIM_(ldab), float *d__, float *e, SimTK_C_ *q, SimTK_FDIM_(ldq), SimTK_C_ *pt, SimTK_FDIM_(ldpt), SimTK_C_ *c__, SimTK_FDIM_(ldc), SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(vect));

void cgbcon_(SimTK_FOPT_(norm), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), const SimTK_C_ *ab, SimTK_FDIM_(ldab), const int *ipiv, SimTK_S_INPUT_(anorm), SimTK_S_OUTPUT_(rcond), SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(norm));

void cgbequ_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), const SimTK_C_ *ab, SimTK_FDIM_(ldab), float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, SimTK_INFO_);

void cgbrfs_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_FDIM_(nrhs), const SimTK_C_ *ab, SimTK_FDIM_(ldab), const SimTK_C_ *afb, SimTK_FDIM_(ldafb), const int *ipiv, const SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *x, SimTK_FDIM_(ldx), float *ferr, float *berr, SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(trans));

void cgbsv_(SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_FDIM_(nrhs), SimTK_C_ *ab, SimTK_FDIM_(ldab), int *ipiv, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_INFO_);

void cgbsvx_(SimTK_FOPT_(fact), SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_FDIM_(nrhs), SimTK_C_ *ab, SimTK_FDIM_(ldab), SimTK_C_ *afb, SimTK_FDIM_(ldafb), int *ipiv, SimTK_CHAR_OUTPUT_(equed), float *r__, float *c__, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *x, SimTK_FDIM_(ldx), SimTK_S_OUTPUT_(rcond), float *ferr, float *berr, SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(trans), SimTK_FLEN_(equed));

void cgbtf2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_C_ *ab, SimTK_FDIM_(ldab), int *ipiv, SimTK_INFO_);

void cgbtrf_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_C_ *ab, SimTK_FDIM_(ldab), int *ipiv, SimTK_INFO_);

void cgbtrs_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_FDIM_(nrhs), const SimTK_C_ *ab, SimTK_FDIM_(ldab), const int *ipiv, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(trans));

void cgebak_(SimTK_FOPT_(job), SimTK_FOPT_(side), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), SimTK_S_INPUT_(scale), SimTK_I_INPUT_(m), SimTK_C_ *v, SimTK_FDIM_(ldv), SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(side));

void cgebal_(SimTK_FOPT_(job), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_I_OUTPUT_(ilo), SimTK_I_OUTPUT_(ihi), SimTK_S_OUTPUT_(scale), SimTK_INFO_, SimTK_FLEN_(job));

void cgebd2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), float *d__, float *e, SimTK_C_ *tauq, SimTK_C_ *taup, SimTK_C_ *work, SimTK_INFO_);

void cgebrd_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), float *d__, float *e, SimTK_C_ *tauq, SimTK_C_ *taup, SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void cgecon_(SimTK_FOPT_(norm), SimTK_FDIM_(n), const SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_S_INPUT_(anorm), SimTK_S_OUTPUT_(rcond), SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(norm));

void cgeequ_(SimTK_FDIM_(m), SimTK_FDIM_(n), const SimTK_C_ *a, SimTK_FDIM_(lda), float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, SimTK_INFO_);

void cgees_(SimTK_FOPT_(jobvs), SimTK_FOPT_(sort), SimTK_SELECT_C select, SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), int *sdim, SimTK_C_ *w, SimTK_C_ *vs, SimTK_FDIM_(ldvs), SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, int *bwork, SimTK_INFO_, SimTK_FLEN_(jobvs), SimTK_FLEN_(sort));

void cgeesx_(SimTK_FOPT_(jobvs), SimTK_FOPT_(sort), SimTK_SELECT_C select, SimTK_FOPT_(sense), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), int *sdim, SimTK_C_ *w, SimTK_C_ *vs, SimTK_FDIM_(ldvs), SimTK_S_OUTPUT_(rconde), SimTK_S_OUTPUT_(rcondv), SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, int *bwork, SimTK_INFO_, SimTK_FLEN_(jobvs), SimTK_FLEN_(sort), SimTK_FLEN_(sense));

void cgeev_(SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *w, SimTK_C_ *vl, SimTK_FDIM_(ldvl), SimTK_C_ *vr, SimTK_FDIM_(ldvr), SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, SimTK_INFO_, SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr));

void cgeevx_(SimTK_FOPT_(balanc), SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), SimTK_FOPT_(sense), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *w, SimTK_C_ *vl, SimTK_FDIM_(ldvl), SimTK_C_ *vr, SimTK_FDIM_(ldvr), SimTK_I_OUTPUT_(ilo), SimTK_I_OUTPUT_(ihi), SimTK_S_OUTPUT_(scale), float *abnrm, SimTK_S_OUTPUT_(rconde), SimTK_S_OUTPUT_(rcondv), SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, SimTK_INFO_, SimTK_FLEN_(balanc), SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr), SimTK_FLEN_(sense));

void cgegs_(SimTK_FOPT_(jobvsl), SimTK_FOPT_(jobvsr), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *alpha, SimTK_C_ *beta, SimTK_C_ *vsl, SimTK_FDIM_(ldvsl), SimTK_C_ *vsr, SimTK_FDIM_(ldvsr), SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, SimTK_INFO_, SimTK_FLEN_(jobvsl), SimTK_FLEN_(jobvsr));

void cgegv_(SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *alpha, SimTK_C_ *beta, SimTK_C_ *vl, SimTK_FDIM_(ldvl), SimTK_C_ *vr, SimTK_FDIM_(ldvr), SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, SimTK_INFO_, SimTK_FLEN_(jobvsl), SimTK_FLEN_(jobvsr));


void cgehd2_(SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_FDIM_(ihi), SimTK_C_ *a, SimTK_I_INPUT_(lda), SimTK_C_ *tau, SimTK_C_ *work, SimTK_INFO_);

void cgehrd_(SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *tau, SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void cgelq2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *tau, SimTK_C_ *work, SimTK_INFO_);

void cgelqf_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *tau, SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void cgels_(SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(trans));
void cgelss_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), float *s, SimTK_S_INPUT_(rcond), SimTK_I_OUTPUT_(rank), SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, SimTK_INFO_);

/* gelsx is deprecated; use gelsy instead */
void cgelsx_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), int *jpvt, SimTK_S_INPUT_(rcond), SimTK_I_OUTPUT_(rank), SimTK_C_ *work, float *rwork, SimTK_INFO_);
void cgelsy_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), int *jpvt, SimTK_S_INPUT_(rcond), SimTK_I_OUTPUT_(rank), SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, SimTK_INFO_);

void cgeql2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *tau, SimTK_C_ *work, SimTK_INFO_);

void cgeqlf_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *tau, SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void cgeqp3_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), int *jpvt, SimTK_C_ *tau, SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, SimTK_INFO_);

void cgeqpf_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), int *jpvt, SimTK_C_ *tau, SimTK_C_ *work, float *rwork, SimTK_INFO_);

void cgeqr2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *tau, SimTK_C_ *work, SimTK_INFO_);

void cgeqrf_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *tau, SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void cgerfs_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *af, SimTK_FDIM_(ldaf), const int *ipiv, const SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *x, SimTK_FDIM_(ldx), float *ferr, float *berr, SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(trans));

void cgerq2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *tau, SimTK_C_ *work, SimTK_INFO_);

void cgerqf_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *tau, SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void cgesc2_(SimTK_FDIM_(n), const SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *rhs, int *ipiv, int *jpiv, SimTK_S_OUTPUT_(scale));

void cgesv_(SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_C_ *a, SimTK_FDIM_(lda), int *ipiv, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_INFO_);

void cgesvx_(SimTK_FOPT_(fact), SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *af, SimTK_FDIM_(ldaf), int *ipiv, SimTK_CHAR_OUTPUT_(equed), float *r__, float *c__, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *x, SimTK_FDIM_(ldx), SimTK_S_OUTPUT_(rcond), float *ferr, float *berr, SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(trans), SimTK_FLEN_(equed));


void cgetc2_(SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), int *ipiv, int *jpiv, SimTK_INFO_);

void cgetf2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), int *ipiv, SimTK_INFO_);

void cgetrf_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), int *ipiv, SimTK_INFO_);

void cgetri_(SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), const int *ipiv, SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void cgetrs_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_C_ *a, SimTK_FDIM_(lda), const int *ipiv, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(trans));

void cggbak_(SimTK_FOPT_(job), SimTK_FOPT_(side), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), const float *lscale, const float *rscale, SimTK_FDIM_(m), SimTK_C_ *v, SimTK_FDIM_(ldv), SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(side));

void cggbal_(SimTK_FOPT_(job), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_I_OUTPUT_(ilo), SimTK_I_OUTPUT_(ihi), float *lscale, float *rscale, float *work, SimTK_INFO_, SimTK_FLEN_(job));

void cgges_(SimTK_FOPT_(jobvsl), SimTK_FOPT_(jobvsr), SimTK_FOPT_(sort), SimTK_SELECT_2C selctg, SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), int *sdim, SimTK_C_ *alpha, SimTK_C_ *beta, SimTK_C_ *vsl, SimTK_FDIM_(ldvsl), SimTK_C_ *vsr, SimTK_FDIM_(ldvsr), SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, int *bwork, SimTK_INFO_, SimTK_FLEN_(jobvsl), SimTK_FLEN_(jobvsr), SimTK_FLEN_(sort));

void cggesx_(SimTK_FOPT_(jobvsl), SimTK_FOPT_(jobvsr), SimTK_FOPT_(sort), SimTK_SELECT_2C selctg, SimTK_FOPT_(sense), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_I_OUTPUT_(sdim), SimTK_C_ *alpha, SimTK_C_ *beta, SimTK_C_ *vsl, SimTK_FDIM_(ldvsl), SimTK_C_ *vsr, SimTK_FDIM_(ldvsr), SimTK_S_OUTPUT_(rconde), SimTK_S_OUTPUT_(rcondv), SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, int *iwork, SimTK_FDIM_(liwork), int *bwork, SimTK_INFO_, SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr), SimTK_FLEN_(sort), SimTK_FLEN_(sense));

void cggev_(SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *alpha, SimTK_C_ *beta, SimTK_C_ *vl, SimTK_FDIM_(ldvl), SimTK_C_ *vr, SimTK_FDIM_(ldvr), SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, SimTK_INFO_, SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr));

void cggevx_(SimTK_FOPT_(balanc), SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), SimTK_FOPT_(sense), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *alpha, SimTK_C_ *beta, SimTK_C_ *vl, SimTK_FDIM_(ldvl), SimTK_C_ *vr, SimTK_FDIM_(ldvr), SimTK_I_OUTPUT_(ilo), SimTK_I_OUTPUT_(ihi), float *lscale, float *rscale, SimTK_S_OUTPUT_(abnrm), SimTK_S_OUTPUT_(bbnrm), SimTK_S_OUTPUT_(rconde), SimTK_S_OUTPUT_(rcondv), SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, int *iwork, int *bwork, SimTK_INFO_, SimTK_FLEN_(balanc), SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr), SimTK_FLEN_(sense));

void cggglm_(SimTK_FDIM_(n), SimTK_FDIM_(m),  SimTK_FDIM_(p), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *d__, SimTK_C_ *x, SimTK_C_ *y, SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void cgghrd_(SimTK_FOPT_(compq), SimTK_FOPT_(compz), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *q, SimTK_FDIM_(ldq), SimTK_C_ *z__, SimTK_FDIM_(ldz), SimTK_INFO_, SimTK_FLEN_(compq), SimTK_FLEN_(compz));

void cgglse_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(p), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *c__, SimTK_C_ *d__, SimTK_C_ *x, SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void cggqrf_(SimTK_FDIM_(n), SimTK_FDIM_(m), SimTK_FDIM_(p), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *taua, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *taub, SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void cggrqf_(SimTK_FDIM_(m), SimTK_FDIM_(p), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *taua, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *taub, SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void cggsvd_(SimTK_FOPT_(jobu), SimTK_FOPT_(jobv), SimTK_FOPT_(jobq), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(p), SimTK_I_OUTPUT_(k), SimTK_I_OUTPUT_(l), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), float *alpha, float *beta, SimTK_C_ *u, SimTK_FDIM_(ldu), SimTK_C_ *v, SimTK_FDIM_(ldv), SimTK_C_ *q, SimTK_FDIM_(ldq), SimTK_C_ *work, float *rwork, int *iwork, SimTK_INFO_, SimTK_FLEN_(jobu), SimTK_FLEN_(jobv), SimTK_FLEN_(jobq));

void cggsvp_(SimTK_FOPT_(jobu), SimTK_FOPT_(jobv), SimTK_FOPT_(jobq), SimTK_FDIM_(m), SimTK_FDIM_(p), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_S_INPUT_(tola), SimTK_S_INPUT_(tolb),  SimTK_I_OUTPUT_(k), SimTK_I_OUTPUT_(l), SimTK_C_ *u, SimTK_FDIM_(ldu), SimTK_C_ *v, SimTK_FDIM_(ldv), SimTK_C_ *q, SimTK_FDIM_(ldq), int *iwork, float *rwork, SimTK_C_ *tau, SimTK_C_ *work, SimTK_INFO_, SimTK_FLEN_(jobu), SimTK_FLEN_(jobv), SimTK_FLEN_(jobq));

void cgtcon_(SimTK_FOPT_(norm), SimTK_FDIM_(n), SimTK_C_ *dl, const SimTK_C_ *d__, const SimTK_C_ *du, const SimTK_C_ *du2, const int *ipiv, SimTK_S_INPUT_(anorm), SimTK_S_OUTPUT_(rcond), SimTK_C_ *work, SimTK_INFO_, SimTK_FLEN_(norm));

void cgtrfs_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_C_ *dl, const SimTK_C_ *d__, const SimTK_C_ *du, const SimTK_C_ *dlf, const SimTK_C_ *df, const SimTK_C_ *duf, const SimTK_C_ *du2, const int *ipiv, const SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *x, SimTK_FDIM_(ldx), float *ferr, float *berr, SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(trans));

void cgtsv_(SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_C_ *dl, SimTK_C_ *d__, SimTK_C_ *du, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_INFO_);

void cgtsvx_(SimTK_FOPT_(fact), SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_C_ *dl, const SimTK_C_ *d__, const SimTK_C_ *du, SimTK_C_ *dlf, SimTK_C_ *df, SimTK_C_ *duf, SimTK_C_ *du2, int *ipiv, const SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *x, SimTK_FDIM_(ldx), SimTK_S_OUTPUT_(rcond), float *ferr, float *berr, SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(trans));

void cgttrf_(SimTK_FDIM_(n), SimTK_C_ *dl, SimTK_C_ *d__, SimTK_C_ *du, SimTK_C_ *du2, int *ipiv, SimTK_INFO_);

void cgttrs_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_C_ *dl, const SimTK_C_ *d__, const SimTK_C_ *du, const SimTK_C_ *du2, const int *ipiv, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(trans));

void cgtts2_(int *itrans, SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_C_ *dl, const SimTK_C_ *d__, const SimTK_C_ *du, const SimTK_C_ *du2, const int *ipiv, SimTK_C_ *b, SimTK_FDIM_(ldb));

void chbev_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_C_ *ab, SimTK_FDIM_(ldab), float *w, SimTK_C_ *z__, SimTK_FDIM_(ldz), SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void chbevd_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_C_ *ab, SimTK_FDIM_(ldab), float *w, SimTK_C_ *z__, SimTK_FDIM_(ldz), SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, SimTK_FDIM_(lrwork), int *iwork, SimTK_FDIM_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void chbevx_(SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_C_ *ab, SimTK_FDIM_(ldab), SimTK_C_ *q, SimTK_FDIM_(ldq), SimTK_S_INPUT_(vl), SimTK_S_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_S_INPUT_(abstol), SimTK_I_OUTPUT_(m), float *w, SimTK_C_ *z__, SimTK_FDIM_(ldz), SimTK_C_ *work, float *rwork, int *iwork, int *ifail, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void chbgst_(SimTK_FOPT_(vect), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(ka), SimTK_FDIM_(kb), SimTK_C_ *ab, SimTK_FDIM_(ldab), const SimTK_C_ *bb, SimTK_FDIM_(ldbb), SimTK_C_ *x, SimTK_FDIM_(ldx), SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(vect), SimTK_FLEN_(uplo));

void chbgv_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(ka), SimTK_FDIM_(kb), SimTK_C_ *ab, SimTK_FDIM_(ldab), SimTK_C_ *bb, SimTK_FDIM_(ldbb), float *w, SimTK_C_ *z__, SimTK_FDIM_(ldz), SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void chbgvx_(SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(ka), SimTK_FDIM_(kb), SimTK_C_ *ab, SimTK_FDIM_(ldab), SimTK_C_ *bb, SimTK_FDIM_(ldbb), SimTK_C_ *q, SimTK_FDIM_(ldq), SimTK_S_INPUT_(vl), SimTK_S_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_S_INPUT_(abstol), SimTK_I_OUTPUT_(m), float *w, SimTK_C_ *z__, SimTK_FDIM_(ldz), SimTK_C_ *work, float *rwork, int *iwork, int *ifail, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void chbtrd_(SimTK_FOPT_(vect), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_C_ *ab, SimTK_FDIM_(ldab), float *d__, float *e, SimTK_C_ *q, SimTK_FDIM_(ldq), SimTK_C_ *work, SimTK_INFO_, SimTK_FLEN_(vect), SimTK_FLEN_(uplo));

void checon_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), const int *ipiv, SimTK_S_INPUT_(anorm), SimTK_S_OUTPUT_(rcond), SimTK_C_ *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void cheev_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), float *w, SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void cheevd_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), float *w, SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, SimTK_FDIM_(lrwork), int *iwork, SimTK_FDIM_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void cheevr_(SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_S_INPUT_(vl), SimTK_S_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_S_INPUT_(abstol),  SimTK_I_OUTPUT_(m), float *w, SimTK_C_ *z__, SimTK_FDIM_(ldz), int *isuppz, SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, SimTK_FDIM_(lrwork), int *iwork, SimTK_FDIM_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void cheevx_(SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_S_INPUT_(vl), SimTK_S_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_S_INPUT_(abstol),  SimTK_I_OUTPUT_(m), float *w, SimTK_C_ *z__, SimTK_FDIM_(ldz), SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, int *iwork, int *ifail, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void chegs2_(SimTK_I_INPUT_(itype), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void chegst_(SimTK_I_INPUT_(itype), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void chegv_(SimTK_I_INPUT_(itype), SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), float *w, SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void chegvd_(SimTK_I_INPUT_(itype), SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), float *w, SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, SimTK_FDIM_(lrwork), int *iwork, SimTK_FDIM_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void chegvx_(SimTK_I_INPUT_(itype), SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_S_INPUT_(vl), SimTK_S_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_S_INPUT_(abstol),  SimTK_I_OUTPUT_(m), float *w, SimTK_C_ *z__, SimTK_FDIM_(ldz), SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, int *iwork, int *ifail, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void cherfs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *af, SimTK_FDIM_(ldaf), const int *ipiv, const SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *x, SimTK_FDIM_(ldx), float *ferr, float *berr, SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void chesv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_C_ *a, SimTK_FDIM_(lda), const int *ipiv, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(uplo));

void chesvx_(SimTK_FOPT_(fact), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *af, SimTK_FDIM_(ldaf), int *ipiv, const SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *x, SimTK_FDIM_(ldx), SimTK_S_OUTPUT_(rcond), float *ferr, float *berr, SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(uplo));

void chetf2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), int *ipiv, SimTK_INFO_, SimTK_FLEN_(uplo));

void chetrd_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), float *d__, float *e, SimTK_C_ *tau, SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(uplo));

void chetrf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), int *ipiv, SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(uplo));

void chetri_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), const int *ipiv, SimTK_C_ *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void chetrs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_C_ *a, SimTK_FDIM_(lda), const int *ipiv, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void chgeqz_(SimTK_FOPT_(job), SimTK_FOPT_(compq), SimTK_FOPT_(compz), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *alpha, SimTK_C_ *beta, SimTK_C_ *q, SimTK_FDIM_(ldq), SimTK_C_ *z__, SimTK_FDIM_(ldz), SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(comq), SimTK_FLEN_(compz));

void chpcon_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const SimTK_C_ *ap, const int *ipiv, SimTK_S_INPUT_(anorm), SimTK_S_OUTPUT_(rcond), SimTK_C_ *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void chpev_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *ap, float *w, SimTK_C_ *z__, SimTK_FDIM_(ldz), SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void chpevd_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *ap, float *w, SimTK_C_ *z__, SimTK_FDIM_(ldz), SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, SimTK_FDIM_(lrwork), int *iwork, SimTK_FDIM_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void chpevx_(SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *ap, SimTK_S_INPUT_(vl), SimTK_S_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_S_INPUT_(abstol), SimTK_I_OUTPUT_(m), float *w, SimTK_C_ *z__, SimTK_FDIM_(ldz), SimTK_C_ *work, float *rwork, int *iwork, int *ifail, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void chpgst_(SimTK_I_INPUT_(itype), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *ap, const SimTK_C_ *bp, SimTK_INFO_, SimTK_FLEN_(uplo));

void chpgv_(SimTK_I_INPUT_(itype), SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *ap, SimTK_C_ *bp, float *w, SimTK_C_ *z__, SimTK_FDIM_(ldz), SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void chpgvd_(SimTK_I_INPUT_(itype), SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *ap, SimTK_C_ *bp, float *w, SimTK_C_ *z__, SimTK_FDIM_(ldz), SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, SimTK_FDIM_(lrwork), int *iwork, SimTK_FDIM_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void chpgvx_( SimTK_I_INPUT_(itype), SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *ap, SimTK_C_ *bp, SimTK_S_INPUT_(vl), SimTK_S_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_S_INPUT_(abstol), SimTK_I_OUTPUT_(m), float *w, SimTK_C_ *z__, SimTK_FDIM_(ldz), SimTK_C_ *work, float *rwork, int *iwork, int *ifail, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void chprfs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_C_ *ap, const SimTK_C_ *afp, const int *ipiv, const SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *x, SimTK_FDIM_(ldx), float *ferr, float *berr, SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void chpsv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_C_ *ap, int *ipiv, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void chpsvx_(SimTK_FOPT_(fact), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_C_ *ap, const SimTK_C_ *afp, const int *ipiv, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *x, SimTK_FDIM_(ldx), SimTK_S_OUTPUT_(rcond), float *ferr, float *berr, SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(uplo));

void chptrd_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *ap, float *d__, float *e, SimTK_C_ *tau, SimTK_INFO_, SimTK_FLEN_(uplo));

void chptrf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *ap, int *ipiv, SimTK_INFO_, SimTK_FLEN_(uplo));

void chptri_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *ap, const int *ipiv, SimTK_C_ *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void chptrs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_C_ *ap, const int *ipiv, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void chsein_(SimTK_FOPT_(side), SimTK_FOPT_(eigsrc), SimTK_FOPT_(initv), const int *select, SimTK_FDIM_(n), const SimTK_C_ *h__, SimTK_FDIM_(ldh), SimTK_C_ *w, SimTK_C_ *vl, SimTK_FDIM_(ldvl), SimTK_C_ *vr, SimTK_FDIM_(ldvr), SimTK_FDIM_(mm), SimTK_I_OUTPUT_(m), SimTK_C_ *work, float *rwork, int *ifaill, int *ifailr, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(eigsrc), SimTK_FLEN_(initv));

void chseqr_(SimTK_FOPT_(job), SimTK_FOPT_(compz), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), SimTK_C_ *h__, SimTK_FDIM_(ldh), SimTK_C_ *w, SimTK_C_ *z__, SimTK_FDIM_(ldz), SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(compz));

void clabrd_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(nb), SimTK_C_ *a, SimTK_FDIM_(lda), float *d__, float *e, SimTK_C_ *tauq, SimTK_C_ *taup, SimTK_C_ *x, SimTK_FDIM_(ldx), SimTK_C_ *y, SimTK_FDIM_(ldy));

void clacgv_(SimTK_FDIM_(n), SimTK_C_ *x, SimTK_FINC_(x));

void clacon_(SimTK_FDIM_(n), SimTK_C_ *v, SimTK_C_ *x, SimTK_S_OUTPUT_(est),  SimTK_I_OUTPUT_(kase) );

void clacn2_( SimTK_FDIM_(n), SimTK_C_ *v, SimTK_C_ *x, SimTK_S_OUTPUT_(est), SimTK_I_OUTPUT_(kase), int *isave   );

void clacp2_(SimTK_FOPT_(uplo), SimTK_FDIM_(m), SimTK_FDIM_(n), const float *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_FLEN_(uplo));

void clacpy_(SimTK_FOPT_(uplo), SimTK_FDIM_(m), SimTK_FDIM_(n), const SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_FLEN_(uplo));

void clacrm_(SimTK_FDIM_(m), SimTK_FDIM_(n), const SimTK_C_ *a, SimTK_FDIM_(lda), const float *b, SimTK_FDIM_(ldb), const SimTK_C_ *c__, SimTK_FDIM_(ldc), float *rwork);

void clacrt_(SimTK_FDIM_(n), SimTK_C_ *cx, SimTK_FINC_(x), SimTK_C_ *cy, SimTK_FINC_(y), const SimTK_C_ *c__, const SimTK_C_ *s);

void claed0_(SimTK_I_INPUT_(qsiz), SimTK_FDIM_(n), float *d__, float *e, SimTK_C_ *q, SimTK_FDIM_(ldq), SimTK_C_ *qstore, SimTK_FDIM_(ldqs), float *rwork, int *iwork, SimTK_INFO_);

void claed7_(SimTK_FDIM_(n), SimTK_I_INPUT_(cutpnt), SimTK_FDIM_(qsiz), SimTK_I_INPUT_(tlvls), SimTK_I_INPUT_(curlvl), SimTK_I_INPUT_(curpbm), float *d__, SimTK_C_ *q, SimTK_FDIM_(ldq), SimTK_S_INPUT_(rho), int *indxq, float *qstore, int *qptr, const int *prmptr, const int *perm, const int *givptr, const int *givcol, const float *givnum, SimTK_C_ *work, float *rwork, int *iwork, SimTK_INFO_);

void claed8_(SimTK_I_OUTPUT_(k), SimTK_FDIM_(n), SimTK_FDIM_(qsiz), SimTK_C_ *q, SimTK_FDIM_(ldq), float *d__, SimTK_S_INPUT_(rho), SimTK_I_INPUT_(cutpnt), float *z__, float *dlamda, SimTK_C_ *q2, SimTK_FDIM_(ldq2), float *w, int *indxp, int *indx, const int *indxq, int *perm, int *givptr, int *givcol, float *givnum, SimTK_INFO_);

void claein_(SimTK_I_INPUT_(rightv), SimTK_I_INPUT_(noinit), SimTK_FDIM_(n), const SimTK_C_ *h__, SimTK_FDIM_(ldh), SimTK_C_INPUT_(w), SimTK_C_ *v, SimTK_C_ *b, SimTK_FDIM_(ldb), float *rwork, SimTK_S_INPUT_(eps3), SimTK_S_INPUT_(smlnum), SimTK_INFO_);

void claesy_(SimTK_C_INPUT_(a), SimTK_C_INPUT_(b), SimTK_C_INPUT_(c__), SimTK_C_OUTPUT_(rt1), SimTK_C_OUTPUT_(rt2), SimTK_C_OUTPUT_(evscal), SimTK_C_OUTPUT_(cs1), SimTK_C_OUTPUT_(sn1) );

void claev2_(SimTK_C_INPUT_(a), SimTK_C_INPUT_(b), SimTK_C_INPUT_(c__), SimTK_S_OUTPUT_(rt1), SimTK_S_OUTPUT_(rt2), SimTK_S_OUTPUT_(cs1), SimTK_C_OUTPUT_(sn1) );

void clags2_(SimTK_I_INPUT_(upper), SimTK_S_INPUT_(a1), SimTK_C_INPUT_(a2), SimTK_S_INPUT_(a3), SimTK_S_INPUT_(b1), SimTK_C_INPUT_(b2), SimTK_S_INPUT_(b3), SimTK_S_OUTPUT_(csu), SimTK_C_OUTPUT_(snu), SimTK_S_OUTPUT_(csv), SimTK_C_OUTPUT_(snv), SimTK_S_OUTPUT_(csq), SimTK_C_OUTPUT_(snq) );

void clagtm_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_S_INPUT_(alpha), const SimTK_C_ *dl, const SimTK_C_ *d__, const SimTK_C_ *du, const SimTK_C_ *x, SimTK_FDIM_(ldx), SimTK_S_INPUT_(beta), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_FLEN_(trans));

void clag2z_( SimTK_I_INPUT_(m), SimTK_I_INPUT_(n), SimTK_C_ *sa, SimTK_FDIM_(ldsa), const SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_INFO_);

void clahef_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nb), SimTK_FDIM_(kb), SimTK_C_ *a, SimTK_FDIM_(lda), int *ipiv, SimTK_C_ *w, SimTK_FDIM_(ldw), SimTK_INFO_, SimTK_FLEN_(uplo));

void clahqr_(SimTK_I_INPUT_(wantt), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), SimTK_C_ *h__, SimTK_FDIM_(ldh), SimTK_C_ *w, SimTK_I_INPUT_(iloz), SimTK_I_INPUT_(ihiz), SimTK_C_ *z__, SimTK_FDIM_(ldz), SimTK_INFO_);

void clahrd_(SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_FDIM_(nb), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *tau, SimTK_C_ *t, SimTK_FDIM_(ldt), SimTK_C_ *y, SimTK_FDIM_(ldy));

void clahr2_(SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_FDIM_(nb), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *tau, SimTK_C_ *t, SimTK_FDIM_(ldt), SimTK_C_ *y, SimTK_FDIM_(ldy));

void claic1_(SimTK_I_INPUT_(job), SimTK_FDIM_(j), const SimTK_C_ *x, SimTK_S_INPUT_(sest), const SimTK_C_ *w, SimTK_C_INPUT_(gamma), SimTK_S_OUTPUT_(sestpr), SimTK_C_OUTPUT_(s), SimTK_C_OUTPUT_(c__) );

void clals0_(SimTK_I_INPUT_(icompq), SimTK_I_INPUT_(nl), SimTK_I_INPUT_(nr), SimTK_I_INPUT_(sqre), SimTK_FDIM_(nrhs), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *bx, SimTK_FDIM_(ldbx), int *perm, int *givptr, const int *givcol, SimTK_FDIM_(ldgcol), float *givnum, SimTK_FDIM_(ldgnum), float *poles, float *difl, float *difr, float *z__, SimTK_FDIM_(k), float *c__, float *s, float *rwork, SimTK_INFO_);

void clalsa_(SimTK_I_INPUT_(icompq), SimTK_I_INPUT_(smlsiz), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *bx, SimTK_FDIM_(ldbx), float *u, SimTK_FDIM_(ldu), float *vt, SimTK_FDIM_(k), float *difl, float *difr, float *z__, float *poles, int *givptr, const int *givcol, SimTK_FDIM_(ldgcol), int *perm, float *givnum, float *c__, float *s, float *rwork, int *iwork, SimTK_INFO_);

void clapll_(SimTK_FDIM_(n), SimTK_C_ *x, SimTK_FINC_(x), SimTK_C_ *y, SimTK_FINC_(y), SimTK_S_OUTPUT_(ssmin));

void clapmt_(SimTK_I_INPUT_(forwrd), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_C_ *x, SimTK_FDIM_(ldx), SimTK_FDIM_(k));

void claqgb_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_C_ *ab, SimTK_FDIM_(ldab), float *r__, float *c__, SimTK_S_OUTPUT_(rowcnd), SimTK_S_OUTPUT_(colcnd), SimTK_S_INPUT_(amax), SimTK_CHAR_OUTPUT_(equed), SimTK_FLEN_(equed));

void claqge_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), const float *r__, const float *c__, SimTK_S_INPUT_(rowcnd), SimTK_S_INPUT_(colcnd), SimTK_S_INPUT_(amax),  SimTK_CHAR_OUTPUT_(equed), SimTK_FLEN_(equed));

void claqhb_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_C_ *ab, SimTK_FDIM_(ldab), float *s, SimTK_S_INPUT_(scond), SimTK_S_INPUT_(amax), SimTK_CHAR_OUTPUT_(equed), SimTK_FLEN_(uplo), SimTK_FLEN_(equed));

void claqhe_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), const float *s, SimTK_S_INPUT_(scond), SimTK_S_INPUT_(amax), SimTK_CHAR_OUTPUT_(equed), SimTK_FLEN_(uplo), SimTK_FLEN_(equed));

void claqhp_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *ap, const float *s, SimTK_S_INPUT_(scond), SimTK_S_INPUT_(amax), SimTK_CHAR_OUTPUT_(equed), SimTK_FLEN_(uplo), SimTK_FLEN_(equed));

void claqp2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_I_INPUT_(offset), SimTK_C_ *a, SimTK_FDIM_(lda), int *jpvt, SimTK_C_ *tau, float *vn1, float *vn2, SimTK_C_ *work);

void claqps_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_I_INPUT_(offset), SimTK_I_INPUT_(nb), SimTK_I_OUTPUT_(kb), SimTK_C_ *a, SimTK_FDIM_(lda), int *jpvt, SimTK_C_ *tau, float *vn1, float *vn2, SimTK_C_ *auxv, SimTK_C_ *f, SimTK_FDIM_(ldf));

void claqr0_( SimTK_I_INPUT_(wantt), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), SimTK_C_ *h, SimTK_FDIM_(ldh), SimTK_C_ *w, SimTK_I_INPUT_(iloz), SimTK_I_INPUT_(ihiz), SimTK_C_ *z, SimTK_FDIM_(ldz), SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_ );

void claqr1_(SimTK_FDIM_(n), const SimTK_C_ *h, SimTK_FDIM_(ldh), SimTK_C_INPUT_(s1), SimTK_C_INPUT_(s2), SimTK_C_ *v );

void claqr2_( SimTK_I_INPUT_(wantt), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), SimTK_I_INPUT_(ktop), SimTK_I_INPUT_(kbot),  SimTK_I_INPUT_(nw), SimTK_C_ *h, SimTK_FDIM_(ldh), SimTK_I_INPUT_(iloz), SimTK_I_INPUT_(ihiz),  SimTK_C_ *z, SimTK_FDIM_(ldz), SimTK_I_OUTPUT_(ns), SimTK_I_OUTPUT_(nd), SimTK_C_ *sh, SimTK_C_ *v, SimTK_FDIM_(ldv), SimTK_I_INPUT_(nh), SimTK_C_ *t, SimTK_FDIM_(ldt), SimTK_I_INPUT_(nv), SimTK_C_ *wv, SimTK_I_INPUT_(ldwv), SimTK_C_ *work, SimTK_FDIM_(lwork) );

void claqr3_( SimTK_I_INPUT_(wantt), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), SimTK_I_INPUT_(ktop), SimTK_I_INPUT_(kbot),  SimTK_I_INPUT_(nw), SimTK_C_ *h, SimTK_FDIM_(ldh), SimTK_I_INPUT_(iloz), SimTK_I_INPUT_(ihiz),  SimTK_C_ *z, SimTK_FDIM_(ldz), SimTK_I_OUTPUT_(ns), SimTK_I_OUTPUT_(nd), SimTK_C_ *sh, SimTK_C_ *v, SimTK_FDIM_(ldv), SimTK_I_INPUT_(nh), SimTK_C_ *t, SimTK_FDIM_(ldt), SimTK_I_INPUT_(nv), SimTK_C_ *wv, SimTK_I_INPUT_(ldwv), SimTK_C_ *work, SimTK_FDIM_(lwork) );

void claqr4_( SimTK_I_INPUT_(wantt), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), SimTK_C_ *h, SimTK_FDIM_(ldh), SimTK_C_ *w, SimTK_I_INPUT_(iloz), SimTK_I_INPUT_(ihiz), SimTK_C_ *z, SimTK_FDIM_(ldz), SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_ );

void claqr5_( SimTK_I_INPUT_(wantt), SimTK_I_INPUT_(wantz), SimTK_I_INPUT_(kacc22), SimTK_FDIM_(n), SimTK_I_INPUT_(ktop), SimTK_I_INPUT_(kbot),  SimTK_I_INPUT_(nshifts), const SimTK_C_ *s, SimTK_C_ *h, SimTK_FDIM_(ldh), SimTK_I_INPUT_(iloz), SimTK_I_INPUT_(ihiz),  SimTK_C_ *z, SimTK_FDIM_(ldz), SimTK_C_ *v, SimTK_FDIM_(ldv), SimTK_C_ *u, SimTK_FDIM_(ldu), SimTK_I_INPUT_(nh), SimTK_C_ *wh, SimTK_FDIM_(ldwh), SimTK_I_INPUT_(nv), SimTK_C_ *wv, SimTK_FDIM_(ldwv)  );

void claqsb_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_C_ *ab, SimTK_FDIM_(ldab), float *s, SimTK_S_INPUT_(scond), SimTK_S_INPUT_(amax), SimTK_CHAR_OUTPUT_(equed), SimTK_FLEN_(uplo), SimTK_FLEN_(equed));

void claqsp_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *ap, const float *s, SimTK_S_INPUT_(scond), SimTK_S_INPUT_(amax), SimTK_CHAR_OUTPUT_(equed), SimTK_FLEN_(uplo), SimTK_FLEN_(equed));

void claqsy_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), const float *s, SimTK_S_INPUT_(scond), SimTK_S_INPUT_(amax), SimTK_CHAR_OUTPUT_(equed), SimTK_FLEN_(uplo), SimTK_FLEN_(equed));

void clar1v_(SimTK_FDIM_(n), SimTK_I_INPUT_(b1), SimTK_I_INPUT_(bn), SimTK_S_INPUT_(lambda), const float *d__, const float *l, const float *ld, const float *lld, SimTK_S_INPUT_(pivmin), SimTK_S_INPUT_(gaptol), SimTK_C_ *z__, SimTK_I_INPUT_(wantnc), SimTK_I_OUTPUT_(negcnt), SimTK_S_OUTPUT_(ztz), SimTK_S_OUTPUT_(mingma), SimTK_I_OUTPUT_(r__), int *isuppz, SimTK_S_OUTPUT_(nrminv), SimTK_S_OUTPUT_(resid), SimTK_S_OUTPUT_(rqcorr), float *work);

void clar2v_(SimTK_FDIM_(n), SimTK_C_ *x, SimTK_C_ *y, SimTK_C_ *z__, SimTK_FINC_(x), const float *c__, const SimTK_C_ *s, SimTK_FINC_(c));

void clarcm_(SimTK_FDIM_(m), SimTK_FDIM_(n), const float *a, SimTK_FDIM_(lda), const SimTK_C_ *b, SimTK_FDIM_(ldb), const SimTK_C_ *c__, SimTK_FDIM_(ldc), float *rwork);

void clarf_(SimTK_FOPT_(side), SimTK_FDIM_(m), SimTK_FDIM_(n), const SimTK_C_ *v, SimTK_FINC_(v), SimTK_C_INPUT_(tau), SimTK_C_ *c__, SimTK_FDIM_(ldc), SimTK_C_ *work, SimTK_FLEN_(side));

void clarfb_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FOPT_(direct), SimTK_FOPT_(storev), SimTK_I_OUTPUT_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_C_ *v, SimTK_FDIM_(ldv), const SimTK_C_ *t, SimTK_FDIM_(ldt), SimTK_C_ *c__, SimTK_FDIM_(ldc), SimTK_C_ *work, SimTK_FDIM_(ldwork), SimTK_FLEN_(side), SimTK_FLEN_(trans), SimTK_FLEN_(direct), SimTK_FLEN_(storev));

void clarfg_(SimTK_FDIM_(n), SimTK_C_OUTPUT_(alpha), SimTK_C_ *x, int *incx, SimTK_C_OUTPUT_(tau) );

void clarft_(SimTK_FOPT_(direct), SimTK_FOPT_(storev), SimTK_FDIM_(n), SimTK_I_INPUT_(k), SimTK_C_ *v, SimTK_FDIM_(ldv), const SimTK_C_ *tau, SimTK_C_ *t, SimTK_FDIM_(ldt), SimTK_FLEN_(direct), SimTK_FLEN_(storev));

void clarfx_(SimTK_FOPT_(side), SimTK_FDIM_(m), SimTK_FDIM_(n), const SimTK_C_ *v, const SimTK_C_ *tau, SimTK_C_ *c__, SimTK_FDIM_(ldc), SimTK_C_ *work, SimTK_FLEN_(side));

void clargv_(SimTK_FDIM_(n), SimTK_C_ *x, SimTK_FINC_(x), SimTK_C_ *y, SimTK_FINC_(y), float *c__, SimTK_FINC_(c));

void clarnv_(SimTK_I_INPUT_(idist), int *iseed, SimTK_FDIM_(n), SimTK_C_ *x);

void clarrv_(SimTK_FDIM_(n), SimTK_S_INPUT_(vl), SimTK_S_INPUT_(vu), float *d__, float *l, SimTK_S_INPUT_(pivmin), const int *isplit, SimTK_FDIM_(m), SimTK_I_INPUT_(dol), SimTK_I_INPUT_(dou), SimTK_S_INPUT_(minrgp), SimTK_S_INPUT_(rtol1), SimTK_S_INPUT_(rtol2), float *w, float *werr, float *wgap, const int *iblock, const int *indexw, const float *gers,  SimTK_C_ *z__, SimTK_FDIM_(ldz), int *isuppz, float *work, int *iwork, SimTK_INFO_);

void clartg_(const SimTK_C_ *f, const SimTK_C_ *g, float *cs, SimTK_C_ *sn, SimTK_C_ *r__);

void clartv_(SimTK_FDIM_(n), SimTK_C_ *x, SimTK_FINC_(x), SimTK_C_ *y, SimTK_FINC_(y), const float *c__, const SimTK_C_ *s, SimTK_FINC_(c));

void clarz_(SimTK_FOPT_(side), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_I_INPUT_(L), const SimTK_C_ *v, SimTK_FINC_(v), SimTK_C_INPUT_(tau), SimTK_C_ *c__, SimTK_FDIM_(ldc), SimTK_C_ *work, SimTK_FLEN_(side));

void clarzb_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FOPT_(direct), SimTK_FOPT_(storev), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_FDIM_(l), const SimTK_C_ *v, SimTK_FDIM_(ldv), const SimTK_C_ *t, SimTK_FDIM_(ldt), SimTK_C_ *c__, SimTK_FDIM_(ldc), SimTK_C_ *work, SimTK_FDIM_(ldwork), SimTK_FLEN_(side), SimTK_FLEN_(trans), SimTK_FLEN_(direct), SimTK_FLEN_(storev));

void clarzt_(SimTK_FOPT_(direct), SimTK_FOPT_(storev), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_C_ *v, SimTK_FDIM_(ldv), const SimTK_C_ *tau, SimTK_C_ *t, SimTK_FDIM_(ldt), SimTK_FLEN_(direct), SimTK_FLEN_(storev));

void clascl_(SimTK_FOPT_(type__), SimTK_FDIM_(kl), SimTK_FDIM_(ku), const float *cfrom, const float *cto, SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(type));

void claset_(SimTK_FOPT_(uplo), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_C_INPUT_(alpha), SimTK_C_INPUT_(beta),  SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_FLEN_(uplo));

void clasr_(SimTK_FOPT_(side), SimTK_FOPT_(pivot), SimTK_FOPT_(direct), SimTK_FDIM_(m), SimTK_FDIM_(n), const float *c__, const float *s, SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_FLEN_(side), SimTK_FLEN_(pivot), SimTK_FLEN_(direct));

void classq_(SimTK_FDIM_(n), const SimTK_C_ *x, SimTK_FINC_(x), SimTK_S_OUTPUT_(scale), float *sumsq);

void claswp_(SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_I_OUTPUT_(k1), SimTK_I_OUTPUT_(k2), const int *ipiv, SimTK_FINC_(x));

void clasyf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_I_INPUT_(nb), SimTK_I_INPUT_(kb), SimTK_C_ *a, SimTK_FDIM_(lda), int *ipiv, SimTK_C_ *w, SimTK_FDIM_(ldw), SimTK_INFO_, SimTK_FLEN_(uplo));




void clatbs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FOPT_(normin), SimTK_FDIM_(n), SimTK_I_INPUT_(kd), const SimTK_C_ *ab, SimTK_FDIM_(ldab), SimTK_C_ *x, SimTK_S_OUTPUT_(scale), float *cnorm, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag), SimTK_FLEN_(normin));

void clatdf_(SimTK_I_INPUT_(ijob), SimTK_FDIM_(n), const SimTK_C_ *z__, SimTK_FDIM_(ldz), SimTK_C_ *rhs, SimTK_S_OUTPUT_(rdsum), SimTK_S_OUTPUT_(rdscal), const int *ipiv, const int *jpiv);

void clatps_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FOPT_(normin), SimTK_FDIM_(n), const SimTK_C_ *ap, SimTK_C_ *x, SimTK_S_OUTPUT_(scale), float *cnorm, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag), SimTK_FLEN_(normin));

void clatrd_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nb), SimTK_C_ *a, SimTK_FDIM_(lda), float *e, SimTK_C_ *tau, SimTK_C_ *w, SimTK_FDIM_(ldw), SimTK_FLEN_(uplo));

void clatrs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FOPT_(normin), SimTK_FDIM_(n), const SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *x, SimTK_S_OUTPUT_(scale), float *cnorm, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag), SimTK_FLEN_(normin));

void clatrz_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(l), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *tau, SimTK_C_ *work);

void clatzm_(SimTK_FOPT_(side), SimTK_FDIM_(m), SimTK_FDIM_(n), const SimTK_C_ *v, SimTK_FINC_(v), SimTK_C_INPUT_(tau), SimTK_C_ *c1, SimTK_C_ *c2, SimTK_FDIM_(ldc), SimTK_C_ *work, SimTK_FLEN_(side));

void clauu2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo));

void clauum_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo));

void cpbcon_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), const SimTK_C_ *ab, SimTK_FDIM_(ldab), SimTK_S_INPUT_(anorm), SimTK_S_OUTPUT_(rcond), SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void cpbequ_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_I_INPUT_(kd), const SimTK_C_ *ab, SimTK_FDIM_(ldab), SimTK_S_OUTPUT_(s), SimTK_S_OUTPUT_(scond), SimTK_S_OUTPUT_(amax), SimTK_INFO_, SimTK_FLEN_(uplo));

void cpbrfs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_FDIM_(nrhs), const SimTK_C_ *ab, SimTK_FDIM_(ldab), const SimTK_C_ *afb, SimTK_FDIM_(ldafb), const SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *x, SimTK_FDIM_(ldx), float *ferr, float *berr, SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void cpbstf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_C_ *ab, SimTK_FDIM_(ldab), SimTK_INFO_, SimTK_FLEN_(uplo));

void cpbsv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_FDIM_(nrhs), SimTK_C_ *ab, SimTK_FDIM_(ldab), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));



void cpbsvx_(SimTK_FOPT_(fact), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_FDIM_(nrhs), SimTK_C_ *ab, SimTK_FDIM_(ldab), SimTK_C_ *afb, SimTK_FDIM_(ldafb), SimTK_CHAR_OUTPUT_(equed), float *s, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *x, SimTK_FDIM_(ldx), SimTK_S_INPUT_(rcond), float *ferr, float *berr, SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(uplo), SimTK_FLEN_(equed));

void cpbtf2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_C_ *ab, SimTK_FDIM_(ldab), SimTK_INFO_, SimTK_FLEN_(uplo));

void cpbtrf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_C_ *ab, SimTK_FDIM_(ldab), SimTK_INFO_, SimTK_FLEN_(uplo));

void cpbtrs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_FDIM_(nrhs), const SimTK_C_ *ab, SimTK_FDIM_(ldab), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void cpocon_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_S_INPUT_(anorm), SimTK_S_OUTPUT_(rcond), SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void cpoequ_(SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), float *s,  SimTK_S_OUTPUT_(scond), SimTK_S_OUTPUT_(amax), SimTK_INFO_);

void cporfs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *af, SimTK_FDIM_(ldaf), const SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *x, SimTK_FDIM_(ldx), float *ferr, float *berr, SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void cposv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void cposvx_(SimTK_FOPT_(fact), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *af, SimTK_FDIM_(ldaf), SimTK_CHAR_OUTPUT_(equed), float *s, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *x, SimTK_FDIM_(ldx), SimTK_S_OUTPUT_(rcond), float *ferr, float *berr, SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(uplo), SimTK_FLEN_(equed));

void cpotf2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo));

void cpotrf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo));

void cpotri_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo));

void cpotrs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void cppcon_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const SimTK_C_ *ap, SimTK_S_INPUT_(anorm), SimTK_S_OUTPUT_(rcond), SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void cppequ_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const SimTK_C_ *ap, float *s, SimTK_S_OUTPUT_(scond), SimTK_S_OUTPUT_(amax), SimTK_INFO_, SimTK_FLEN_(uplo));

void cpprfs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_C_ *ap, const SimTK_C_ *afp, const SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *x, SimTK_FDIM_(ldx), float *ferr, float *berr, SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void cppsv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_C_ *ap, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void cppsvx_(SimTK_FOPT_(fact), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_C_ *ap, SimTK_C_ *afp, SimTK_CHAR_OUTPUT_(equed), float *s, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *x, SimTK_FDIM_(ldx), SimTK_S_OUTPUT_(rcond), float *ferr, float *berr, SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(uplo), SimTK_FLEN_(equed));

void cpptrf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *ap, SimTK_INFO_, SimTK_FLEN_(uplo));

void cpptri_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *ap, SimTK_INFO_, SimTK_FLEN_(uplo));

void cpptrs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_C_ *ap, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void cptcon_(SimTK_FDIM_(n), const float *d__, const SimTK_C_ *e, SimTK_S_INPUT_(anorm), SimTK_S_OUTPUT_(rcond), float *rwork, SimTK_INFO_);

void cptrfs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const float *d__, const SimTK_C_ *e, const float *df, const SimTK_C_ *ef, const SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *x, SimTK_FDIM_(ldx), float *ferr, float *berr, SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void cptsv_(SimTK_FDIM_(n), SimTK_FDIM_(nrhs), float *d__, SimTK_C_ *e, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_INFO_);

void cptsvx_(SimTK_FOPT_(fact), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const float *d__, const SimTK_C_ *e, float *df, SimTK_C_ *ef, const SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *x, SimTK_FDIM_(ldx), SimTK_S_OUTPUT_(rcond), float *ferr, float *berr, SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(fact));

void cpttrf_(SimTK_FDIM_(n), float *d__, SimTK_C_ *e, SimTK_INFO_);

void cpttrs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const float *d__, const SimTK_C_ *e, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void cptts2_(SimTK_I_INPUT_(iuplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const float *d__, const SimTK_C_ *e, SimTK_C_ *b, SimTK_FDIM_(ldb));

void crot_(SimTK_FDIM_(n), SimTK_C_ *cx, SimTK_FINC_(x), SimTK_C_ *cy, SimTK_FINC_(y), const float *c__, const SimTK_C_ *s);

void cspcon_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const SimTK_C_ *ap, const int *ipiv, SimTK_S_INPUT_(anorm), SimTK_S_OUTPUT_(rcond), SimTK_C_ *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void cspmv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_INPUT_(alpha), const SimTK_C_ *ap, const SimTK_C_ *x, SimTK_FINC_(x), SimTK_C_INPUT_(beta), SimTK_C_ *y, int *incy, SimTK_FLEN_(uplo));

void cspr_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_INPUT_(alpha), const SimTK_C_ *x, SimTK_FINC_(x), SimTK_C_ *ap, SimTK_FLEN_(uplo));

void csprfs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_C_ *ap, const SimTK_C_ *afp, const int *ipiv, const SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *x, SimTK_FDIM_(ldx), float *ferr, float *berr, SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void cspsv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_C_ *ap, int *ipiv, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void cspsvx_(SimTK_FOPT_(fact), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_C_ *ap, SimTK_C_ *afp, int *ipiv, const SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *x, SimTK_FDIM_(ldx), SimTK_S_OUTPUT_(rcond), float *ferr, float *berr, SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(uplo));

void csptrf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *ap, int *ipiv, SimTK_INFO_, SimTK_FLEN_(uplo));

void csptri_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *ap, const int *ipiv, SimTK_C_ *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void csptrs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_C_ *ap, const int *ipiv, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void csrscl_(SimTK_FDIM_(n), SimTK_S_INPUT_(sa), SimTK_C_ *sx, SimTK_FINC_(x));

void cstedc_(SimTK_FOPT_(compz), SimTK_FDIM_(n), float *d__, float *e, SimTK_C_ *z__, SimTK_FDIM_(ldz), SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, int *lrwork, int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(compz));

void cstein_(SimTK_FDIM_(n), const float *d__, const float *e, SimTK_FDIM_(m), const float *w, const int *iblock, const int *isplit, SimTK_C_ *z__, SimTK_FDIM_(ldz), float *work, int *iwork, int *ifail, SimTK_INFO_);

void cstemr_( SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FDIM_(n), float *d, float *e, SimTK_S_INPUT_(vl), SimTK_S_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_I_OUTPUT_(m), float *w, SimTK_C_ *z, SimTK_FDIM_(ldz), SimTK_I_INPUT_(nzc), SimTK_I_OUTPUT_(isuppz), SimTK_I_OUTPUT_(tryrac), float *work, SimTK_FDIM_(lwork), int *iwork, SimTK_FDIM_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range) ); 

void csteqr_(SimTK_FOPT_(compz), SimTK_FDIM_(n), float *d__, float *e, SimTK_C_ *z__, SimTK_FDIM_(ldz), float *work, SimTK_INFO_, SimTK_FLEN_(compz));

void csycon_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const SimTK_C_ *a, SimTK_FDIM_(lda), const int *ipiv, SimTK_S_INPUT_(anorm), SimTK_S_OUTPUT_(rcond), SimTK_C_ *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void csymv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_INPUT_(alpha), const SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *x, SimTK_FINC_(x), SimTK_C_INPUT_(beta), SimTK_C_ *y, SimTK_FINC_(y), SimTK_FLEN_(uplo));

void csyr_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_INPUT_(alpha), const SimTK_C_ *x, SimTK_FINC_(x), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_FLEN_(uplo));

void csyrfs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *af, SimTK_FDIM_(ldaf), const int *ipiv, const SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *x, SimTK_FDIM_(ldx), float *ferr, float *berr, SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void csysv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_C_ *a, SimTK_FDIM_(lda), int *ipiv, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(uplo));

void csysvx_(SimTK_FOPT_(fact), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *af, SimTK_FDIM_(ldaf), int *ipiv, const SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *x, SimTK_FDIM_(ldx), SimTK_S_OUTPUT_(rcond), float *ferr, float *berr, SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(uplo));

void csytf2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), int *ipiv, SimTK_INFO_, SimTK_FLEN_(uplo));

void csytrf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), int *ipiv, SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(uplo));

void csytri_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), const int *ipiv, SimTK_C_ *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void csytrs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_C_ *a, SimTK_FDIM_(lda), const int *ipiv, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void ctbcon_(SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(kd), const SimTK_C_ *ab, SimTK_FDIM_(ldab), SimTK_S_OUTPUT_(rcond), SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(norm), SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

void ctbrfs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_FDIM_(nrhs), const SimTK_C_ *ab, SimTK_FDIM_(ldab), const SimTK_C_ *b, SimTK_FDIM_(ldb), const SimTK_C_ *x, SimTK_FDIM_(ldx), float *ferr, float *berr, SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));

void ctbtrs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_FDIM_(nrhs), const SimTK_C_ *ab, SimTK_FDIM_(ldab), const SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));

void ctgevc_(SimTK_FOPT_(side), SimTK_FOPT_(howmny), const int *select, SimTK_FDIM_(n), const SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *vl, SimTK_FDIM_(ldvl), SimTK_C_ *vr, SimTK_FDIM_(ldvr), SimTK_I_INPUT_(mm), SimTK_I_OUTPUT_(m), SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(howmny));




void ctgex2_(SimTK_I_INPUT_(wantq), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *q, SimTK_FDIM_(ldq), SimTK_C_ *z__, SimTK_FDIM_(ldz), SimTK_I_INPUT_(j1), SimTK_INFO_);

void ctgexc_(SimTK_I_INPUT_(wantq), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *q, SimTK_FDIM_(ldq), SimTK_C_ *z__, SimTK_FDIM_(ldz), SimTK_I_OUTPUT_(ifst), SimTK_I_OUTPUT_(ilst), SimTK_INFO_);

void ctgsen_(SimTK_I_INPUT_(ijob), SimTK_I_INPUT_(wantq), SimTK_I_INPUT_(wantz), const int *select, SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *alpha, SimTK_C_INPUT_(beta), SimTK_C_ *q, SimTK_FDIM_(ldq), SimTK_C_ *z__, SimTK_FDIM_(ldz), SimTK_I_OUTPUT_(m), SimTK_S_OUTPUT_(pl), SimTK_S_OUTPUT_(pr), float *dif, SimTK_C_ *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_);

void ctgsja_(SimTK_FOPT_(jobu), SimTK_FOPT_(jobv), SimTK_FOPT_(jobq), SimTK_FDIM_(m), SimTK_FDIM_(p), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_FDIM_(l), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_S_INPUT_(tola), SimTK_S_INPUT_(tolb), float *alpha, float *beta, SimTK_C_ *u, SimTK_FDIM_(ldu), SimTK_C_ *v, SimTK_FDIM_(ldv), SimTK_C_ *q, SimTK_FDIM_(ldq), SimTK_C_ *work, SimTK_I_OUTPUT_(ncycle), SimTK_INFO_, SimTK_FLEN_(jobu), SimTK_FLEN_(jobv), SimTK_FLEN_(jobq));

void ctgsna_(SimTK_FOPT_(job), SimTK_FOPT_(howmny), const int *select, SimTK_FDIM_(n), const SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *b, SimTK_FDIM_(ldb), const SimTK_C_ *vl, SimTK_FDIM_(ldvl), const SimTK_C_ *vr, SimTK_FDIM_(ldvr), float *s, float *dif, SimTK_FDIM_(mm), SimTK_I_OUTPUT_(m), SimTK_C_ *work, SimTK_FDIM_(lwork), int *iwork, SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(howmny));

void ctgsy2_(SimTK_FOPT_(trans), SimTK_I_INPUT_(ijob), SimTK_FDIM_(m), SimTK_FDIM_(n), const SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *c__, SimTK_FDIM_(ldc), const SimTK_C_ *d__, SimTK_FDIM_(ldd), const SimTK_C_ *e, SimTK_FDIM_(lde), SimTK_C_ *f, SimTK_FDIM_(ldf), SimTK_S_OUTPUT_(scale), SimTK_S_OUTPUT_(rdsum), SimTK_S_OUTPUT_(rdscal), SimTK_INFO_, SimTK_FLEN_(trans));

void ctgsyl_(SimTK_FOPT_(trans), SimTK_I_INPUT_(ijob), SimTK_FDIM_(m), SimTK_FDIM_(n), const SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *c__, SimTK_FDIM_(ldc), const SimTK_C_ *d__, SimTK_FDIM_(ldd), const SimTK_C_ *e, SimTK_FDIM_(lde), SimTK_C_ *f, SimTK_FDIM_(ldf), SimTK_S_OUTPUT_(scale), SimTK_S_OUTPUT_(dif), SimTK_C_ *work, SimTK_FDIM_(lwork), int *iwork, SimTK_INFO_, SimTK_FLEN_(trans));

void ctpcon_(SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), const SimTK_C_ *ap, SimTK_S_OUTPUT_(rcond), SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(norm), SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

void ctprfs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_C_ *ap, const SimTK_C_ *b, SimTK_FDIM_(ldb), const SimTK_C_ *x, SimTK_FDIM_(ldx), float *ferr, float *berr, SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));

void ctptri_(SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_C_ *ap, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

void ctptrs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_C_ *ap, SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));

void ctrcon_(SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), const SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_S_OUTPUT_(rcond), SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(norm), SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

void ctrevc_(SimTK_FOPT_(side), SimTK_FOPT_(howmny), const int *select, SimTK_FDIM_(n), SimTK_C_ *t, SimTK_FDIM_(ldt), SimTK_C_ *vl, SimTK_FDIM_(ldvl), SimTK_C_ *vr, SimTK_FDIM_(ldvr), SimTK_FDIM_(mm), SimTK_I_OUTPUT_(m), SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(howmny));

void ctrexc_(SimTK_FOPT_(compq), SimTK_FDIM_(n), SimTK_C_ *t, SimTK_FDIM_(ldt), SimTK_C_ *q, SimTK_FDIM_(ldq), SimTK_I_INPUT_(ifst), SimTK_I_INPUT_(ilst), SimTK_INFO_, SimTK_FLEN_(compq));


void ctrrfs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *b, SimTK_FDIM_(ldb), const SimTK_C_ *x, SimTK_FDIM_(ldx), float *ferr, float *berr, SimTK_C_ *work, float *rwork, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));

void ctrsen_(SimTK_FOPT_(job), SimTK_FOPT_(compq), const int *select, SimTK_FDIM_(n), SimTK_C_ *t, SimTK_FDIM_(ldt), SimTK_C_ *q, SimTK_FDIM_(ldq), SimTK_C_ *w, SimTK_I_OUTPUT_(m), SimTK_S_OUTPUT_(s), SimTK_S_OUTPUT_(sep), SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(compq));

void ctrsna_(SimTK_FOPT_(job), SimTK_FOPT_(howmny), const int *select, SimTK_FDIM_(n), const SimTK_C_ *t, SimTK_FDIM_(ldt), const SimTK_C_ *vl, SimTK_FDIM_(ldvl), const SimTK_C_ *vr, SimTK_FDIM_(ldvr), float *s, float *sep, SimTK_FDIM_(mm), SimTK_I_OUTPUT_(m), SimTK_C_ *work, SimTK_FDIM_(ldwork), float *rwork, SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(howmny));

void ctrsyl_(SimTK_FOPT_(trana), SimTK_FOPT_(tranb), SimTK_I_INPUT_(isgn), SimTK_FDIM_(m), SimTK_FDIM_(n), const SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_C_ *c__, SimTK_FDIM_(ldc), SimTK_S_OUTPUT_(scale), SimTK_INFO_, SimTK_FLEN_(trana), SimTK_FLEN_(tranb));

void ctrti2_(SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

void ctrtri_(SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

void ctrtrs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));

void ctzrqf_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *tau, SimTK_INFO_);

void ctzrzf_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *tau, SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void cung2l_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *tau, SimTK_C_ *work, SimTK_INFO_);

void cung2r_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *tau, SimTK_C_ *work, SimTK_INFO_);

void cungbr_(SimTK_FOPT_(vect), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *tau, SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(vect));

void cunghr_(SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *tau, SimTK_C_ *work, SimTK_FDIM_(lwork), int *info);

void cungl2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *tau, SimTK_C_ *work, SimTK_INFO_);

void cunglq_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *tau, SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void cungql_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *tau, SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void cungqr_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *tau, SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void cungr2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *tau, SimTK_C_ *work, SimTK_INFO_);

void cungrq_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *tau, SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void cungtr_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *tau, SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(uplo));

void cunm2l_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *tau, SimTK_C_ *c__, SimTK_FDIM_(ldc), SimTK_C_ *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void cunm2r_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *tau, SimTK_C_ *c__, SimTK_FDIM_(ldc), SimTK_C_ *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void cunmbr_(SimTK_FOPT_(vect), SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *tau, SimTK_C_ *c__, SimTK_FDIM_(ldc), SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(vect), SimTK_FLEN_(side), SimTK_FLEN_(trans));

void cunmhr_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), const SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *tau, SimTK_C_ *c__, SimTK_FDIM_(ldc), SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void cunml2_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *tau, SimTK_C_ *c__, SimTK_FDIM_(ldc), SimTK_C_ *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void cunmlq_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *tau, SimTK_C_ *c__, SimTK_FDIM_(ldc), SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void cunmql_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *tau, SimTK_C_ *c__, SimTK_FDIM_(ldc), SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void cunmqr_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *tau, SimTK_C_ *c__, SimTK_FDIM_(ldc), SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void cunmr2_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *tau, SimTK_C_ *c__, SimTK_FDIM_(ldc), SimTK_C_ *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void cunmr3_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_I_INPUT_(l), const SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *tau, SimTK_C_ *c__, SimTK_FDIM_(ldc), SimTK_C_ *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void cunmrq_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *tau, SimTK_C_ *c__, SimTK_FDIM_(ldc), SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void cunmrz_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_FDIM_(l), const SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *tau, SimTK_C_ *c__, SimTK_FDIM_(ldc), SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void cunmtr_(SimTK_FOPT_(side), SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), const SimTK_C_ *a, SimTK_FDIM_(lda), const SimTK_C_ *tau, SimTK_C_ *c__, SimTK_FDIM_(ldc), SimTK_C_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(uplo), SimTK_FLEN_(trans));

void cupgtr_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const SimTK_C_ *ap, const SimTK_C_ *tau, SimTK_C_ *q, SimTK_FDIM_(ldq), SimTK_C_ *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void cupmtr_(SimTK_FOPT_(side), SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), const SimTK_C_ *ap, const SimTK_C_ *tau, SimTK_C_ *c__, SimTK_FDIM_(ldc), SimTK_C_ *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(uplo), SimTK_FLEN_(trans));

void dbdsdc_(SimTK_FOPT_(uplo), SimTK_FOPT_(compq), SimTK_FDIM_(n), double *d__, double *e, double *u, SimTK_FDIM_(ldu), double *vt, SimTK_FDIM_(ldvt), double *q, int *iq, double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(compq));

void dbdsqr_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(ncvt), SimTK_FDIM_(nru), SimTK_FDIM_(ncc), double *d__, double *e, double *vt, SimTK_FDIM_(ldvt), double *u, SimTK_FDIM_(ldu), double *c__, SimTK_FDIM_(ldc), double *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void ddisna_(SimTK_FOPT_(job), SimTK_FDIM_(m), SimTK_FDIM_(n), const double *d__, double *sep, SimTK_INFO_, SimTK_FLEN_(job));

void dgbbrd_(SimTK_FOPT_(vect), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(ncc), SimTK_FDIM_(kl), SimTK_FDIM_(ku), double *ab, SimTK_FDIM_(ldab), double *d__, double *e, double *q, SimTK_FDIM_(ldq), double *pt, SimTK_FDIM_(ldpt), double *c__, SimTK_FDIM_(ldc), double *work, SimTK_INFO_, SimTK_FLEN_(vect));

void dgbcon_(SimTK_FOPT_(norm), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), double *ab, SimTK_FDIM_(ldab), const int *ipiv, SimTK_D_INPUT_(anorm), SimTK_D_OUTPUT_(rcond), double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(norm));

void dgbequ_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), double *ab, SimTK_FDIM_(ldab), double *r__, double *c__, double *rowcnd, double *colcnd, double *amax, SimTK_INFO_);

void dgbrfs_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_FDIM_(nrhs), double *ab, SimTK_FDIM_(ldab), double *afb, SimTK_FDIM_(ldafb), const int *ipiv, double *b, SimTK_FDIM_(ldb), double *x, SimTK_FDIM_(ldx), double *ferr, double *berr, double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(trans));

void dgbsv_(SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_FDIM_(nrhs), double *ab, SimTK_FDIM_(ldab), int *ipiv, double *b, SimTK_FDIM_(ldb), SimTK_INFO_);

void dgbsvx_(SimTK_FOPT_(fact), SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_FDIM_(nrhs), double *ab, SimTK_FDIM_(ldab), double *afb, SimTK_FDIM_(ldafb), int *ipiv, SimTK_CHAR_OUTPUT_(equed), double *r__, double *c__, double *b, SimTK_FDIM_(ldb), double *x, SimTK_FDIM_(ldx), SimTK_D_OUTPUT_(rcond), double *ferr, double *berr, double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(trans), SimTK_FLEN_(equed));

void dgbtf2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), double *ab, SimTK_FDIM_(ldab), int *ipiv, SimTK_INFO_);

void dgbtrf_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), double *ab, SimTK_FDIM_(ldab), int *ipiv, SimTK_INFO_);

void dgbtrs_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_FDIM_(nrhs), double *ab, SimTK_FDIM_(ldab), const int *ipiv, double *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(trans));

void dgebak_(SimTK_FOPT_(job), SimTK_FOPT_(side), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), SimTK_D_INPUT_(scale), SimTK_I_INPUT_(m), double *v, SimTK_FDIM_(ldv), SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(side));

void dgebal_(SimTK_FOPT_(job), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), SimTK_I_OUTPUT_(ilo), SimTK_I_OUTPUT_(ihi), SimTK_D_OUTPUT_(scale), SimTK_INFO_, SimTK_FLEN_(job));

void dgebd2_(SimTK_FDIM_(m), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *d__, double *e, double *tauq, double *taup, double *work, SimTK_INFO_);

void dgebrd_(SimTK_FDIM_(m), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *d__, double *e, double *tauq, double *taup, double *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void dgecon_(SimTK_FOPT_(norm), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), SimTK_D_INPUT_(anorm), SimTK_D_OUTPUT_(rcond), double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(norm));

void dgeequ_(SimTK_FDIM_(m), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *r__, double *c__, double *rowcnd, double *colcnd, double *amax, SimTK_INFO_);

void dgees_(SimTK_FOPT_(jobvs), SimTK_FOPT_(sort), SimTK_SELECT_2D select, SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), int *sdim, double *wr, double *wi, double *vs, SimTK_FDIM_(ldvs), double *work, SimTK_FDIM_(lwork), int *bwork, SimTK_INFO_, SimTK_FLEN_(jobvs), SimTK_FLEN_(sort));

void dgeesx_(SimTK_FOPT_(jobvs), SimTK_FOPT_(sort), SimTK_SELECT_2D select, char *sense, SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), int *sdim, double *wr, double *wi, double *vs, SimTK_FDIM_(ldvs), SimTK_D_OUTPUT_(rconde), SimTK_D_OUTPUT_(rcondv), double *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), int *bwork, SimTK_INFO_, SimTK_FLEN_(jobvs), SimTK_FLEN_(sort), SimTK_FLEN_(sense));

void dgeev_(SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *wr, double *wi, double *vl, SimTK_FDIM_(ldvl), double *vr, SimTK_FDIM_(ldvr), double *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr));

void dgeevx_(SimTK_FOPT_(balanc), SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), char *sense, SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *wr, double *wi, double *vl, SimTK_FDIM_(ldvl), double *vr, SimTK_FDIM_(ldvr), SimTK_I_OUTPUT_(ilo), SimTK_I_OUTPUT_(ihi), SimTK_D_OUTPUT_(scale), double *abnrm, SimTK_D_OUTPUT_(rconde), SimTK_D_OUTPUT_(rcondv), double *work, SimTK_FDIM_(lwork), int *iwork, SimTK_INFO_, SimTK_FLEN_(balanc), SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr), SimTK_FLEN_(sense));

void dgegs_(SimTK_FOPT_(jobvsl), SimTK_FOPT_(jobvsr), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), double *alphar, double *alphai, double *beta, double *vsl, SimTK_FDIM_(ldvsl), double *vsr, SimTK_FDIM_(ldvsr), double *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(jobvsl), SimTK_FLEN_(jobvsr));

void dgegv_(SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), double *alphar, double *alphai, double *beta, double *vl, SimTK_FDIM_(ldvl), double *vr, SimTK_FDIM_(ldvr), double *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(jobvsl), SimTK_FLEN_(jobvsr));

void dgehd2_(SimTK_FDIM_(n), SimTK_FDIM_(ilo), SimTK_I_INPUT_(ihi), double *a, SimTK_I_INPUT_(lda), double *tau, double *work, SimTK_INFO_);

void dgehrd_(SimTK_FDIM_(n), SimTK_FDIM_(ilo), SimTK_I_INPUT_(ihi), double *a, SimTK_I_INPUT_(lda), double *tau, double *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void dgelq2_(SimTK_FDIM_(m), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *tau, double *work, SimTK_INFO_);

void dgelqf_(SimTK_FDIM_(m), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *tau, double *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void dgels_(SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), double *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(trans));

void dgelsd_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), double *s, SimTK_D_INPUT_(rcond), SimTK_I_OUTPUT_(rank), double *work, SimTK_FDIM_(lwork), int *iwork, SimTK_INFO_);

void dgelss_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), double *s, SimTK_D_INPUT_(rcond), SimTK_I_OUTPUT_(rank), double *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void dgelsx_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), int *jpvt, SimTK_D_INPUT_(rcond), SimTK_I_OUTPUT_(rank), double *work, SimTK_INFO_);

void dgelsy_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), int *jpvt, SimTK_D_INPUT_(rcond), SimTK_I_OUTPUT_(rank), double *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void dgeql2_(SimTK_FDIM_(m), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *tau, double *work, SimTK_INFO_);

void dgeqlf_(SimTK_FDIM_(m), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *tau, double *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void dgeqp3_(SimTK_FDIM_(m), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), int *jpvt, double *tau, double *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void dgeqpf_(SimTK_FDIM_(m), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), int *jpvt, double *tau, double *work, SimTK_INFO_);

void dgeqr2_(SimTK_FDIM_(m), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *tau, double *work, SimTK_INFO_);

void dgeqrf_(SimTK_FDIM_(m), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *tau, double *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void dgerfs_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), double *a, SimTK_FDIM_(lda), double *af, SimTK_FDIM_(ldaf), const int *ipiv, double *b, SimTK_FDIM_(ldb), double *x, SimTK_FDIM_(ldx), double *ferr, double *berr, double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(trans));

void dgerq2_(SimTK_FDIM_(m), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *tau, double *work, SimTK_INFO_);

void dgerqf_(SimTK_FDIM_(m), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *tau, double *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void dgesc2_(SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *rhs, int *ipiv, int *jpiv, SimTK_D_OUTPUT_(scale));

void dgesdd_(SimTK_FOPT_(jobz), SimTK_FDIM_(m), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *s, double *u, SimTK_FDIM_(ldu), double *vt, SimTK_FDIM_(ldvt), double *work, SimTK_FDIM_(lwork), int *iwork, SimTK_INFO_, SimTK_FLEN_(jobz));

void dgesv_(SimTK_FDIM_(n), SimTK_FDIM_(nrhs), double *a, SimTK_FDIM_(lda), int *ipiv, double *b, SimTK_FDIM_(ldb), SimTK_INFO_);

void dgesvd_(SimTK_FOPT_(jobu), char *jobvt, SimTK_FDIM_(m), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *s, double *u, SimTK_FDIM_(ldu), double *vt, SimTK_FDIM_(ldvt), double *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(jobu), SimTK_FLEN_(jobvt));

void dgesvx_(SimTK_FOPT_(fact), SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), double *a, SimTK_FDIM_(lda), double *af, SimTK_FDIM_(ldaf), int *ipiv, SimTK_CHAR_OUTPUT_(equed), double *r__, double *c__, double *b, SimTK_FDIM_(ldb), double *x, SimTK_FDIM_(ldx), SimTK_D_OUTPUT_(rcond), double *ferr, double *berr, double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(trans), SimTK_FLEN_(equed));

void dgetc2_(SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), int *ipiv, int *jpiv, SimTK_INFO_);

void dgetf2_(SimTK_FDIM_(m), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), int *ipiv, SimTK_INFO_);

void dgetrf_(SimTK_FDIM_(m), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), int *ipiv, SimTK_INFO_);

void dgetri_(SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), const int *ipiv, double *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void dgetrs_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const double *a, SimTK_FDIM_(lda), const int *ipiv, double *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(trans));

void dggbak_(SimTK_FOPT_(job), SimTK_FOPT_(side), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), const double *lscale, const double *rscale, SimTK_FDIM_(m), double *v, SimTK_FDIM_(ldv), SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(side));

void dggbal_(SimTK_FOPT_(job), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), SimTK_I_OUTPUT_(ilo), SimTK_I_OUTPUT_(ihi), double *lscale, double *rscale, double *work, SimTK_INFO_, SimTK_FLEN_(job));

void dgges_(SimTK_FOPT_(jobvsl), SimTK_FOPT_(jobvsr), SimTK_FOPT_(sort), SimTK_SELECT_3D delctg, SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), int *sdim, double *alphar, double *alphai, SimTK_D_INPUT_(beta), double *vsl, SimTK_FDIM_(ldvsl), double *vsr, SimTK_FDIM_(ldvsr), double *work, SimTK_FDIM_(lwork), int *bwork, SimTK_INFO_, SimTK_FLEN_(jobvsl), SimTK_FLEN_(jobvsr), SimTK_FLEN_(sort));

void dggesx_(SimTK_FOPT_(jobvsl), SimTK_FOPT_(jobvsr), SimTK_FOPT_(sort), SimTK_SELECT_3D delctg, SimTK_FOPT_(sense), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), SimTK_I_OUTPUT_(sdim), double *alphar, double *alphai, double *beta, double *vsl, SimTK_FDIM_(ldvsl), double *vsr, SimTK_FDIM_(ldvsr), SimTK_D_OUTPUT_(rconde), SimTK_D_OUTPUT_(rcondv), double *work, SimTK_FDIM_(lwork), int *iwork, SimTK_FDIM_(liwork), int *bwork, SimTK_INFO_, SimTK_FLEN_(jobvsl), SimTK_FLEN_(jobvsr), SimTK_FLEN_(sort), SimTK_FLEN_(sense));

void dggev_(SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), double *alphar, double *alphai, double *beta, double *vl, SimTK_FDIM_(ldvl), double *vr, SimTK_FDIM_(ldvr), double *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr));

void dggevx_(SimTK_FOPT_(balanc), SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), SimTK_FOPT_(sense), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), double *alphar, double *alphai, double *beta, double *vl, SimTK_FDIM_(ldvl), double *vr, SimTK_FDIM_(ldvr), SimTK_I_OUTPUT_(ilo), SimTK_I_OUTPUT_(ihi), double *lscale, double *rscale, SimTK_D_OUTPUT_(abnrm), SimTK_D_OUTPUT_(bbnrm), SimTK_D_OUTPUT_(rconde), SimTK_D_OUTPUT_(rcondv), double *work, SimTK_FDIM_(lwork), int *iwork, int *bwork, SimTK_INFO_, SimTK_FLEN_(balanc), SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr), SimTK_FLEN_(sense));

void dggglm_(SimTK_FDIM_(n), SimTK_FDIM_(m), SimTK_FDIM_(p), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), double *d__, double *x, double *y, double *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void dgghrd_(SimTK_FOPT_(compq), SimTK_FOPT_(compz), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), double *q, SimTK_FDIM_(ldq), double *z__, SimTK_FDIM_(ldz), SimTK_INFO_, SimTK_FLEN_(compq), SimTK_FLEN_(compz));

void dgglse_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(p), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), double *c__, double *d__, double *x, double *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void dggqrf_(SimTK_FDIM_(n), SimTK_FDIM_(m), SimTK_FDIM_(p), double *a, SimTK_FDIM_(lda), double *taua, double *b, SimTK_FDIM_(ldb), double *taub, double *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void dggrqf_(SimTK_FDIM_(m), SimTK_FDIM_(p), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *taua, double *b, SimTK_FDIM_(ldb), double *taub, double *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void dggsvd_(SimTK_FOPT_(jobu), SimTK_FOPT_(jobv), SimTK_FOPT_(jobq), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(p), SimTK_I_OUTPUT_(k), SimTK_I_OUTPUT_(l), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), double *alpha, double *beta, double *u, SimTK_FDIM_(ldu), double *v, SimTK_FDIM_(ldv), double *q, SimTK_FDIM_(ldq), double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(jobu), SimTK_FLEN_(jobv), SimTK_FLEN_(jobq));

void dggsvp_(SimTK_FOPT_(jobu), SimTK_FOPT_(jobv), SimTK_FOPT_(jobq), SimTK_FDIM_(m), SimTK_FDIM_(p), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), SimTK_D_INPUT_(tola), SimTK_D_INPUT_(tolb),  SimTK_I_OUTPUT_(k), SimTK_I_OUTPUT_(l), double *u, SimTK_FDIM_(ldu), double *v, SimTK_FDIM_(ldv), double *q, SimTK_FDIM_(ldq), int *iwork, double *tau, double *work, SimTK_INFO_, SimTK_FLEN_(jobu), SimTK_FLEN_(jobv), SimTK_FLEN_(jobq));

void dgtcon_(SimTK_FOPT_(norm), SimTK_FDIM_(n), double *dl, double *d__, double *du, double *du2, const int *ipiv, SimTK_D_INPUT_(anorm), SimTK_D_OUTPUT_(rcond), double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(norm));

void dgtrfs_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), double *dl, double *d__, double *du, double *dlf, double *df, double *duf, double *du2, const int *ipiv, double *b, SimTK_FDIM_(ldb), double *x, SimTK_FDIM_(ldx), double *ferr, double *berr, double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(trans));

void dgtsv_(SimTK_FDIM_(n), SimTK_FDIM_(nrhs), double *dl, double *d__, double *du, double *b, SimTK_FDIM_(ldb), int *info);

void dgtsvx_(SimTK_FOPT_(fact), SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), double *dl, double *d__, double *du, double *dlf, double *df, double *duf, double *du2, int *ipiv, double *b, SimTK_FDIM_(ldb), double *x, SimTK_FDIM_(ldx), SimTK_D_OUTPUT_(rcond), double *ferr, double *berr, double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(trans));

void dgttrf_(SimTK_FDIM_(n), double *dl, double *d__, double *du, double *du2, int *ipiv, SimTK_INFO_);

void dgttrs_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), double *dl, double *d__, double *du, double *du2, const int *ipiv, double *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(trans));

void dgtts2_(int *itrans, SimTK_FDIM_(n), SimTK_FDIM_(nrhs), double *dl, double *d__, double *du, double *du2, const int *ipiv, double *b, SimTK_FDIM_(ldb));

void dhgeqz_(SimTK_FOPT_(job), SimTK_FOPT_(compq), SimTK_FOPT_(compz), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), double *alphar, double *alphai, double *beta, double *q, SimTK_FDIM_(ldq), double *z__, SimTK_FDIM_(ldz), double *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(compq), SimTK_FLEN_(compz));

void dhsein_(SimTK_FOPT_(side), SimTK_FOPT_(eigsrc), SimTK_FOPT_(initv), const int *select, SimTK_FDIM_(n), const double *h__, SimTK_FDIM_(ldh), double *wr, double *wi, double *vl, SimTK_FDIM_(ldvl), double *vr, SimTK_FDIM_(ldvr), SimTK_FDIM_(mm), SimTK_I_OUTPUT_(m), double *work, int *ifaill, int *ifailr, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(eigsrc), SimTK_FLEN_(initv));

double  dlamch_(SimTK_FOPT_(cmach), SimTK_FLEN_(cmach));

void dhseqr_(SimTK_FOPT_(job), SimTK_FOPT_(compz), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), double *h__, SimTK_FDIM_(ldh), double *wr, double *wi, double *z__, SimTK_FDIM_(ldz), double *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(comqz));

int disnan( SimTK_D_INPUT_(din) );

int dlaisnan( SimTK_D_INPUT_(din1), SimTK_D_INPUT_(din2) );

void dlabad_(SimTK_D_INOUT_(small1), SimTK_D_INOUT_(large));//small -> small1 conflict with definition when Qt is used

void dlabrd_(SimTK_FDIM_(m), SimTK_FDIM_(n),  SimTK_FDIM_(nb), double *a, SimTK_FDIM_(lda), double *d__, double *e, double *tauq, double *taup, double *x, SimTK_FDIM_(ldx), double *y, SimTK_FDIM_(ldy));

void dlacon_(SimTK_FDIM_(n), double *v, double *x, int *isgn, SimTK_D_OUTPUT_(est), SimTK_I_OUTPUT_(kase) );

void dlacn2_( SimTK_FDIM_(n), double *v, double *x, int *isgn, SimTK_D_INPUT_(est), SimTK_I_OUTPUT_(kase), int *isave );

void dlacpy_(SimTK_FOPT_(uplo), SimTK_FDIM_(m), SimTK_FDIM_(n), const double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), SimTK_FLEN_(uplo));

void dladiv_(SimTK_D_INPUT_(a), SimTK_D_INPUT_(b), SimTK_D_INPUT_(c__), SimTK_D_INPUT_(d__), double *p, double *q);

void dlae2_(SimTK_D_INPUT_(a), SimTK_D_INPUT_(b), SimTK_D_INPUT_(c__), double *rt1, double *rt2);

void dlaebz_(SimTK_I_INPUT_(ijob), SimTK_I_INPUT_(nitmax), SimTK_FDIM_(n), SimTK_I_INPUT_(mmax), SimTK_I_INPUT_(minp), SimTK_I_INPUT_(nbmin), SimTK_D_INPUT_(abstol), SimTK_D_INPUT_(reltol), SimTK_D_INPUT_(pivmin), const double *d__, const double *e, const double *e2, int *nval, double *ab, double *c__, SimTK_I_OUTPUT_(mout), int *nab, double *work, int *iwork, SimTK_INFO_);


void dlaed0_(SimTK_I_INPUT_(icompq), SimTK_I_INPUT_(qsiz), SimTK_FDIM_(n), double *d__, double *e, double *q, SimTK_FDIM_(ldq), double *qstore, SimTK_FDIM_(ldqs), double *work, int *iwork, SimTK_INFO_);

void dlaed1_(SimTK_FDIM_(n), double *d__, double *q, SimTK_FDIM_(ldq), int *indxq, SimTK_D_INPUT_(rho), SimTK_I_INPUT_(cutpnt), double *work, int *iwork, SimTK_INFO_);

void dlaed2_(SimTK_FDIM_(k), SimTK_FDIM_(n), SimTK_I_INPUT_(n1), double *d__, double *q, SimTK_FDIM_(ldq), int *indxq, double *rho, double *z__, double *dlamda, double *w, double *q2, int *indx, int *indxc, int *indxp, int *coltyp, SimTK_INFO_);

void dlaed3_(SimTK_FDIM_(k), SimTK_FDIM_(n),  SimTK_I_INPUT_(n1), double *d__, double *q, SimTK_FDIM_(ldq), SimTK_D_INPUT_(rho), double *dlamda, const double *q2, const int *indx, const int *ctot, double *w, double *s, SimTK_INFO_);

void dlaed4_(SimTK_FDIM_(n), SimTK_I_INPUT_(i__), const double *d__, const double *z__, double *delta, SimTK_D_INPUT_(rho), SimTK_D_OUTPUT_(dlam), SimTK_INFO_);

void dlaed5_(SimTK_I_INPUT_(i__), const double *d__, const double *z__, double *delta, SimTK_D_INPUT_(rho), SimTK_D_OUTPUT_(dlam));

void dlaed6_(SimTK_I_INPUT_(kniter), SimTK_I_INPUT_(orgati), SimTK_D_INPUT_(rho), const double *d__, const double *z__, SimTK_D_INPUT_(finit),  SimTK_D_INPUT_(tau), SimTK_INFO_);

void dlaed7_(SimTK_I_INPUT_(icompq), SimTK_FDIM_(n), SimTK_FDIM_(qsiz), SimTK_I_INPUT_(tlvls), SimTK_I_INPUT_(curlvl), SimTK_I_INPUT_(curpbm), double *d__, double *q, SimTK_FDIM_(ldq), int *indxq, SimTK_D_INPUT_(rho), SimTK_I_INPUT_(cutpnt), double *qstore, int *qptr, const int *prmptr, const int *perm, const int *givptr, const int *givcol, const double *givnum, double *work, int *iwork, SimTK_INFO_);

void dlaed8_(SimTK_I_INPUT_(icompq), SimTK_I_OUTPUT_(k), SimTK_FDIM_(n), SimTK_FDIM_(qsiz), double *d__, double *q, SimTK_FDIM_(ldq), const int *indxq, SimTK_D_OUTPUT_(rho), SimTK_I_INPUT_(cutpnt), double *z__, double *dlamda, double *q2, SimTK_FDIM_(ldq2), double *w, int *perm, int *givptr, int *givcol, double *givnum, int *indxp, int *indx, SimTK_INFO_);

void dlaed9_(SimTK_FDIM_(k), int *kstart, int *kstop, SimTK_FDIM_(n), double *d__, double *q, SimTK_FDIM_(ldq), double *rho, double *dlamda, double *w, double *s, SimTK_FDIM_(lds), SimTK_INFO_);

void dlaeda_(SimTK_FDIM_(n), SimTK_I_INPUT_(tlvls), SimTK_I_INPUT_(curlvl), SimTK_I_INPUT_(curpbm), const int *prmptr, const int *perm, const int *givptr, const int *givcol, const double *givnum, const double *q, const int *qptr, double *z__, double *ztemp, SimTK_INFO_);

void dlaein_(SimTK_I_INPUT_(rightv), SimTK_I_INPUT_(noinit), SimTK_FDIM_(n), const double *h__, SimTK_FDIM_(ldh), SimTK_D_INPUT_(wr), SimTK_D_INPUT_(wi), double *vr, double *vi, double *b, SimTK_FDIM_(ldb), double *work, SimTK_D_INPUT_(eps3), SimTK_D_INPUT_(smlnum), SimTK_D_INPUT_(bignum), SimTK_INFO_);

void dlaev2_(SimTK_D_INPUT_(a), SimTK_D_INPUT_(b), SimTK_D_INPUT_(c__), SimTK_D_OUTPUT_(rt1), SimTK_D_OUTPUT_(rt2), SimTK_D_OUTPUT_(cs1), SimTK_D_OUTPUT_(sn1));

void dlaexc_(SimTK_I_INPUT_(wantq), SimTK_FDIM_(n), double *t, SimTK_FDIM_(ldt), double *q, SimTK_FDIM_(ldq), SimTK_I_INPUT_(j1),  SimTK_I_INPUT_(n1), SimTK_I_INPUT_(n2), double *work, SimTK_INFO_);

void dlag2_(const double *a, SimTK_FDIM_(lda), const double *b, SimTK_FDIM_(ldb), SimTK_D_INPUT_(safmin), SimTK_D_OUTPUT_(scale1), SimTK_D_OUTPUT_(scale2), SimTK_D_OUTPUT_(wr1), SimTK_D_OUTPUT_(wr2), SimTK_D_OUTPUT_(wi) );

void dlags2_(SimTK_I_INPUT_(upper), SimTK_D_INPUT_(a1), SimTK_D_INPUT_(a2), SimTK_D_INPUT_(a3), SimTK_D_INPUT_(b1), SimTK_D_INPUT_(b2), SimTK_D_INPUT_(b3), SimTK_D_OUTPUT_(csu), SimTK_D_OUTPUT_(snu), SimTK_D_OUTPUT_(csv), SimTK_D_OUTPUT_(snv), SimTK_D_OUTPUT_(csq), SimTK_D_OUTPUT_(snq) );

void dlagtf_(SimTK_FDIM_(n), double *a, SimTK_D_INPUT_(lambda), double *b, double *c__, SimTK_D_INPUT_(tol), double *d__, int *in, SimTK_INFO_);

void dlagtm_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_D_INPUT_(alpha), const double *dl, const double *d__, const double *du, const double *x, SimTK_FDIM_(ldx), SimTK_D_INPUT_(beta), double *b, SimTK_FDIM_(ldb), SimTK_FLEN_(trans));

void dlag2s_( SimTK_I_INPUT_(m), SimTK_I_INPUT_(n),  const double *a, SimTK_FDIM_(lda), float *sa, SimTK_FDIM_(ldsa),  SimTK_INFO_);

void dlagts_(SimTK_I_INPUT_(job), SimTK_FDIM_(n), const double *a, const double *b, const double *c__, const double *d__, const int *in, double *y, SimTK_D_OUTPUT_(tol), SimTK_INFO_);

void dlagv2_(double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), double *alphar, double *alphai, double *beta, SimTK_D_OUTPUT_(csl), SimTK_D_OUTPUT_(snl), SimTK_D_OUTPUT_(csr), SimTK_D_OUTPUT_(snr));

void dlahqr_(SimTK_I_INPUT_(wantt), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), double *h__, SimTK_FDIM_(ldh), double *wr, double *wi, SimTK_I_INPUT_(iloz), SimTK_I_INPUT_(ihiz), double *z__, SimTK_FDIM_(ldz), SimTK_INFO_);

void dlahrd_(SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_FDIM_(nb), double *a, SimTK_FDIM_(lda), double *tau, double *t, SimTK_FDIM_(ldt), double *y, SimTK_FDIM_(ldy));

void dlahr2_(SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_FDIM_(nb), double *a, SimTK_FDIM_(lda), double *tau, double *t, SimTK_FDIM_(ldt), double *y, SimTK_FDIM_(ldy));

void dlaic1_(SimTK_I_INPUT_(job), SimTK_FDIM_(j), const double *x, SimTK_D_INPUT_(sest), const double *w, SimTK_D_INPUT_(gamma), SimTK_D_OUTPUT_(sestpr), SimTK_D_OUTPUT_(s), SimTK_D_OUTPUT_(c__) );

void dlaln2_(SimTK_I_INPUT_(ltrans), SimTK_I_INPUT_(na), SimTK_I_INPUT_(nw), SimTK_D_INPUT_(smin), SimTK_D_INPUT_(ca), const double *a, SimTK_FDIM_(lda), SimTK_D_INPUT_(d1), SimTK_D_INPUT_(d2), SimTK_D_INPUT_(b), SimTK_FDIM_(ldb), SimTK_D_INPUT_(wr), SimTK_D_INPUT_(wi), double *x, SimTK_FDIM_(ldx), SimTK_D_OUTPUT_(scale), double *xnorm, SimTK_INFO_);

void dlals0_(SimTK_I_INPUT_(icompq), SimTK_I_INPUT_(nl), SimTK_I_INPUT_(nr), SimTK_I_INPUT_(sqre), SimTK_FDIM_(nrhs), double *b, SimTK_FDIM_(ldb), double *bx, SimTK_FDIM_(ldbx), int *perm, int *givptr, const int *givcol, SimTK_FDIM_(ldgcol), double *givnum, SimTK_FDIM_(ldgnum), double *poles, double *difl, double *difr, double *z__, int *k, double *c__, double *s, double *work, SimTK_INFO_);

void dlalsa_(SimTK_I_INPUT_(icompq), SimTK_I_INPUT_(smlsiz), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), double *b, SimTK_FDIM_(ldb), double *bx, SimTK_FDIM_(ldbx), double *u, SimTK_FDIM_(ldu), double *vt, SimTK_FDIM_(k), double *difl, double *difr, double *z__, double *poles, int *givptr, const int *givcol, SimTK_FDIM_(ldgcol), int *perm, double *givnum, double *c__, double *s, double *work, int *iwork, SimTK_INFO_);

void dlalsd_(SimTK_FOPT_(uplo), SimTK_I_INPUT_(smlsiz), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), double *d__, double *e, double *b, SimTK_FDIM_(ldb), SimTK_D_INPUT_(rcond), SimTK_I_OUTPUT_(rank), double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void dlamc1_(int *beta, int *t, int *rnd, int *ieee1);

void dlamc2_(int *beta, int *t, int *rnd, double *eps, int *emin, double *rmin, int *emax, double *rmax);

void dlamc4_(int *emin, double *start, int *base);

void dlamc5_(int *beta, int *p, int *emin, int *ieee, int *emax, double *rmax);

void dlamrg_( SimTK_I_INPUT_(n1), SimTK_I_INPUT_(n2), double *a, int *dtrd1, int *dtrd2, int *index);

void dlanv2_(double *a, double *b, double *c__, double *d__, double *rt1r, double *rt1i, double *rt2r, double *rt2i, double *cs, double *sn);

void dlapll_(SimTK_FDIM_(n), double *x, SimTK_FINC_(x), double *y, SimTK_FINC_(y), SimTK_D_OUTPUT_(ssmin));

void dlapmt_(SimTK_I_INPUT_(forwrd), SimTK_FDIM_(m), SimTK_FDIM_(n), double *x, SimTK_FDIM_(ldx), SimTK_FDIM_(k));

void dlaqgb_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), double *ab, SimTK_FDIM_(ldab), double *r__, double *c__, SimTK_D_OUTPUT_(rowcnd), SimTK_D_OUTPUT_(colcnd), SimTK_D_INPUT_(amax), SimTK_CHAR_OUTPUT_(equed), SimTK_FLEN_(equed));

void dlaqge_(SimTK_FDIM_(m), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *r__, double *c__, SimTK_D_INPUT_(rowcnd), SimTK_D_INPUT_(colcnd), SimTK_D_INPUT_(amax),  SimTK_CHAR_OUTPUT_(equed), SimTK_FLEN_(equed));

void dlaqp2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_I_INPUT_(offset), double *a, SimTK_FDIM_(lda), int *jpvt, double *tau, double *vn1, double *vn2, double *work);

void dlaqps_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_I_INPUT_(offset), SimTK_I_INPUT_(nb), SimTK_I_OUTPUT_(kb), double *a, SimTK_FDIM_(lda), int *jpvt, double *tau, double *vn1, double *vn2, double *auxv, double *f, SimTK_FDIM_(ldf));

void dlaqr0_( SimTK_I_INPUT_(wantt), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), double *h, SimTK_FDIM_(ldh), double *wr, double *wi, SimTK_I_INPUT_(iloz), SimTK_I_INPUT_(ihiz), double *z, SimTK_FDIM_(ldz), double *work, SimTK_FDIM_(lwork), SimTK_INFO_ );

void dlaqr1_(SimTK_FDIM_(n), const double *h, SimTK_FDIM_(ldh), SimTK_D_INPUT_(sr1), SimTK_D_INPUT_(si1), SimTK_D_INPUT_(sr2), SimTK_D_INPUT_(si2), double *v );

void dlaqr2_( SimTK_I_INPUT_(wantt), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), SimTK_I_INPUT_(ktop), SimTK_I_INPUT_(kbot),  SimTK_I_INPUT_(nw), double *h, SimTK_FDIM_(ldh), SimTK_I_INPUT_(iloz), SimTK_I_INPUT_(ihiz),  double *z, SimTK_FDIM_(ldz), SimTK_I_OUTPUT_(ns), SimTK_I_OUTPUT_(nd), double *sr, double *si, double *v, SimTK_FDIM_(ldv), SimTK_I_INPUT_(nh), double *t, SimTK_FDIM_(ldt), SimTK_I_INPUT_(nv), double *wv, SimTK_I_INPUT_(ldwv), double *work, SimTK_FDIM_(lwork) );

void dlaqr3_( SimTK_I_INPUT_(wantt), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), SimTK_I_INPUT_(ktop), SimTK_I_INPUT_(kbot),  SimTK_I_INPUT_(nw), double *h, SimTK_FDIM_(ldh), SimTK_I_INPUT_(iloz), SimTK_I_INPUT_(ihiz),  double *z, SimTK_FDIM_(ldz), SimTK_I_OUTPUT_(ns), SimTK_I_OUTPUT_(nd), double *sr, double *si, double *v, SimTK_FDIM_(ldv), SimTK_I_INPUT_(nh), double *t, SimTK_FDIM_(ldt), SimTK_I_INPUT_(nv), double *wv, SimTK_I_INPUT_(ldwv), double *work, SimTK_FDIM_(lwork) );

void dlaqr4_( SimTK_I_INPUT_(wantt), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), double *h, SimTK_FDIM_(ldh), double *wr, double *wi, SimTK_I_INPUT_(iloz), SimTK_I_INPUT_(ihiz), double *z, SimTK_FDIM_(ldz), double *work, SimTK_FDIM_(lwork), SimTK_INFO_ );

void dlaqr5_( SimTK_I_INPUT_(wantt), SimTK_I_INPUT_(wantz), SimTK_I_INPUT_(kacc22), SimTK_FDIM_(n), SimTK_I_INPUT_(ktop), SimTK_I_INPUT_(kbot),  SimTK_I_INPUT_(nshifts), const double *sr, const double *si, double *h, SimTK_FDIM_(ldh), SimTK_I_INPUT_(iloz), SimTK_I_INPUT_(ihiz),  double *z, SimTK_FDIM_(ldz), double *v, SimTK_FDIM_(ldv), double *u, SimTK_FDIM_(ldu), SimTK_I_INPUT_(nh), double *wh, SimTK_FDIM_(ldwh), SimTK_I_INPUT_(nv), double *wv, SimTK_FDIM_(ldwv)  );

void dlaqsb_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), double *ab, SimTK_FDIM_(ldab), double *s, SimTK_D_INPUT_(scond), SimTK_D_INPUT_(amax), SimTK_CHAR_OUTPUT_(equed),  SimTK_FLEN_(uplo), SimTK_FLEN_(equedk));

void dlaqsp_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *ap, const double *s, SimTK_D_INPUT_(scond), SimTK_D_INPUT_(amax), SimTK_CHAR_OUTPUT_(equed), SimTK_FLEN_(uplo), SimTK_FLEN_(equed));

void dlaqsy_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), const double *s, SimTK_D_INPUT_(scond), SimTK_D_INPUT_(amax), SimTK_CHAR_OUTPUT_(equed), SimTK_FLEN_(uplo), SimTK_FLEN_(equed));

void dlaqtr_( SimTK_I_INPUT_(ltran), SimTK_I_INPUT_(lreal), SimTK_FDIM_(n), const double *t, SimTK_FDIM_(ldt), const double *b, const double *w, SimTK_D_OUTPUT_(scale), double *x, double *work, SimTK_INFO_);

void dlar1v_(SimTK_FDIM_(n), SimTK_I_INPUT_(b1), SimTK_I_INPUT_(bn), SimTK_D_INPUT_(lambda), const double *d__, const double *l, const double *ld, const double *lld, const double *pivmin, const double *gaptol, double *z__, SimTK_I_INPUT_(wantc), SimTK_I_OUTPUT_(negcnt), SimTK_D_OUTPUT_(ztz), SimTK_D_OUTPUT_(mingma), SimTK_I_OUTPUT_(r__), int *isuppz, SimTK_D_OUTPUT_(nrminv), SimTK_D_OUTPUT_(resid), SimTK_D_OUTPUT_(rqcorr), double *work);

void dlar2v_(SimTK_FDIM_(n), double *x, double *y, double *z__, SimTK_FINC_(x), const double *c__, const double *s, SimTK_FINC_(c));

void dlarf_(SimTK_FOPT_(side), SimTK_FDIM_(m), SimTK_FDIM_(n), double *v, SimTK_FINC_(v), double *tau, double *c__, SimTK_FDIM_(ldc), double *work, SimTK_FLEN_(side));

void dlarfb_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FOPT_(direct), SimTK_FOPT_(storev), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), const double *v, SimTK_FDIM_(ldv), const double *t, SimTK_FDIM_(ldt), double *c__, SimTK_FDIM_(ldc), double *work, SimTK_FDIM_(ldwork), SimTK_FLEN_(side), SimTK_FLEN_(trans), SimTK_FLEN_(direct), SimTK_FLEN_(storev));

void dlarfg_(SimTK_FDIM_(n), SimTK_D_OUTPUT_(alpha), double *x, SimTK_FINC_(x), SimTK_D_OUTPUT_(tau) );

void dlarft_(SimTK_FOPT_(direct), SimTK_FOPT_(storev), SimTK_FDIM_(n), SimTK_I_INPUT_(k), double *v, SimTK_FDIM_(ldv), const double *tau, double *t, SimTK_FDIM_(ldt), SimTK_FLEN_(direct), SimTK_FLEN_(storev));

void dlarfx_(SimTK_FOPT_(side), SimTK_FDIM_(m), SimTK_FDIM_(n), const double *v, const double *tau, double *c__, SimTK_FDIM_(ldc), double *work, SimTK_FLEN_(side));

void dlargv_(SimTK_FDIM_(n), double *x, SimTK_FINC_(x), double *y, SimTK_FINC_(y), double *c__, SimTK_FINC_(c));

void dlarnv_(SimTK_I_INPUT_(idist), int *iseed, SimTK_FDIM_(n), double *x);

void dlarra_(SimTK_FDIM_(n), const double *d, double *e, double *e2, SimTK_D_INPUT_(spltol), SimTK_D_INPUT_(tnrm), SimTK_I_OUTPUT_(nsplit), int *isplit, SimTK_INFO_);

void dlarrb_(SimTK_FDIM_(n), const double *d__, const double *lld, SimTK_I_INPUT_(ifirst), SimTK_I_INPUT_(ilast), SimTK_D_INPUT_(rtol1), SimTK_D_INPUT_(rtol2), SimTK_I_INPUT_(offset), double *w, double *wgap, double *werr, double *work, int *iwork, SimTK_D_INPUT_(pivmin), SimTK_D_INPUT_(spdiam),  SimTK_I_INPUT_(twist), SimTK_INFO_);


void dlarrc_(SimTK_FOPT_(jobt), SimTK_FDIM_(n), SimTK_D_INPUT_(vl), SimTK_D_INPUT_(vu), const double *d, const double *e,  SimTK_D_INPUT_(pivmin), SimTK_I_OUTPUT_(eigcnt), SimTK_I_OUTPUT_(lcnt), SimTK_I_OUTPUT_(rcnt), SimTK_INFO_, SimTK_FLEN_(jobt) );

void dlarrd_( SimTK_FOPT_(range), SimTK_FOPT_(order), SimTK_FDIM_(n), SimTK_D_INPUT_(vl), SimTK_D_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), const double *gers, SimTK_D_INPUT_(reltol), const double *d, const double *e, const double *e2, SimTK_D_INPUT_(pivmin), SimTK_I_INPUT_(nsplit), const int *isplit, SimTK_I_OUTPUT_(m), double *w, double *werr, SimTK_D_OUTPUT_(wl), SimTK_D_OUTPUT_(wu), int *iblock, int *indexw, double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(range), SimTK_FLEN_(order) );

void dlarre_(SimTK_FOPT_(range), SimTK_FDIM_(n), SimTK_D_OUTPUT_(vl), SimTK_D_OUTPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), double *d__, double *e, double *e2, SimTK_D_INPUT_(rtol1), SimTK_D_INPUT_(rtol2), SimTK_D_INPUT_(spltol), SimTK_I_OUTPUT_(nsplit),  int *isplit, SimTK_I_OUTPUT_(m), double *w, double *werr, double *wgap, int *iblock, int *indexw, double *gers, double *pivmin, double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(range) );


void dlarrf_(SimTK_FDIM_(n), double *d__, double *l, double *ld, SimTK_I_INPUT_(clstrt), SimTK_I_INPUT_(clend),  const double *w, double *wgap, double *werr, SimTK_D_INPUT_(spdiam), SimTK_D_INPUT_(clgapl), SimTK_D_INPUT_(clgapr), SimTK_D_INPUT_(pivmin), SimTK_D_OUTPUT_(sigma), double *dplus, double *lplus, double *work, SimTK_INFO_);


void dlarrj_(SimTK_FDIM_(n), const double *d, const double *e2, SimTK_I_INPUT_(ifirst), SimTK_I_INPUT_(ilast), SimTK_D_INPUT_(rtol), SimTK_I_INPUT_(offset), double *w, double *werr, double *work, int *iwork, SimTK_D_INPUT_(pivmin), SimTK_D_INPUT_(spdiam), SimTK_INFO_);

void dlarrk_( SimTK_FDIM_(n), SimTK_I_INPUT_(iw), SimTK_D_INPUT_(gl), SimTK_D_INPUT_(gu), const double *d, const double *e2, SimTK_D_INPUT_(pivmin),  SimTK_D_INPUT_(reltol), SimTK_D_OUTPUT_(w), SimTK_D_OUTPUT_(werr), SimTK_INFO_);

void dlarrr_( SimTK_FDIM_(n), const double *d__, double *e, SimTK_INFO_ );

void dlarrv_(SimTK_FDIM_(n), SimTK_D_INPUT_(vl), SimTK_D_INPUT_(vu), double *d__, double *l, SimTK_D_INPUT_(pivmin), const int *isplit, SimTK_FDIM_(m), SimTK_I_INPUT_(dol), SimTK_I_INPUT_(dou), SimTK_D_INPUT_(minrgp), SimTK_D_INPUT_(rtol1), SimTK_D_INPUT_(rtol2), double *w, double *werr, double *wgap, const int *iblock, const int *indexw, const double *gers,  double *z__, SimTK_FDIM_(ldz), int *isuppz, double *work, int *iwork, SimTK_INFO_);

void dlartg_(const double *f, const double *g, double *cs, double *sn, double *r__);

void dlartv_(SimTK_FDIM_(n), double *x, SimTK_FINC_(x), double *y, SimTK_FINC_(y), const double *c__, const double *s, int *incc);

void dlaruv_(int *iseed, SimTK_FDIM_(n), double *x);

void dlarz_(SimTK_FOPT_(side), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_I_INPUT_(L), const double *v, SimTK_FINC_(v), SimTK_D_INPUT_(tau), double *c__, SimTK_FDIM_(ldc), double *work, SimTK_FLEN_(side));

void dlarzb_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FOPT_(direct), SimTK_FOPT_(storev), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_FDIM_(l), const double *v, SimTK_FDIM_(ldv), const double *t, SimTK_FDIM_(ldt), double *c__, SimTK_FDIM_(ldc), double *work, SimTK_FDIM_(ldwork), SimTK_FLEN_(side), SimTK_FLEN_(trans), SimTK_FLEN_(direct), SimTK_FLEN_(storev));

void dlarzt_(SimTK_FOPT_(direct), SimTK_FOPT_(storev), SimTK_FDIM_(n),  SimTK_FDIM_(k), double *v, SimTK_FDIM_(ldv), const double *tau, double *t, SimTK_FDIM_(ldt), SimTK_FLEN_(direct), SimTK_FLEN_(storev));

void dlas2_(double *f, double *g, double *h__, SimTK_D_OUTPUT_(ssmin), double *ssmax);

void dlascl_(SimTK_FOPT_(type__), SimTK_FDIM_(kl), SimTK_FDIM_(ku), const double *cfrom, const double *cto, SimTK_FDIM_(m), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(type));

void dlasd0_(SimTK_FDIM_(n), SimTK_I_INPUT_(sqre), double *d__, double *e, double *u, SimTK_FDIM_(ldu), double *vt, SimTK_FDIM_(ldvt), SimTK_I_INPUT_(smlsiz), int *iwork, double *work, SimTK_INFO_);

void dlasd1_(SimTK_I_INPUT_(nl), SimTK_I_INPUT_(nr), SimTK_I_INPUT_(sqre), double *d__, SimTK_D_INPUT_(alpha), SimTK_D_INPUT_(beta), double *u, SimTK_FDIM_(ldu), double *vt, SimTK_FDIM_(ldvt), int *idxq, int *iwork, double *work, SimTK_INFO_);

void dlasd2_(SimTK_I_INPUT_(nl), SimTK_I_INPUT_(nr), SimTK_I_INPUT_(sqre), SimTK_I_OUTPUT_(k), double *d__, double *z__, SimTK_D_INPUT_(alpha), double *beta, double *u, SimTK_FDIM_(ldu), double *vt, SimTK_FDIM_(ldvt), double *dsigma, double *u2, SimTK_FDIM_(ldu2), double *vt2, SimTK_FDIM_(ldvt2), int *idxp, int *idx, int *idxc, int *idxq, int *coltyp, SimTK_INFO_);

void dlasd3_(SimTK_I_INPUT_(nl), SimTK_I_INPUT_(nr), SimTK_I_INPUT_(sqre), SimTK_I_INPUT_(k), double *d__, double *q, SimTK_FDIM_(ldq), const double *dsigma, const double *u, SimTK_FDIM_(ldu), const double *u2, SimTK_FDIM_(ldu2), const double *vt, SimTK_FDIM_(ldvt), const double *vt2, SimTK_FDIM_(ldvt2), const int *idxc, const int *ctot, const double *z__, SimTK_INFO_);

void dlasd4_(SimTK_FDIM_(n), SimTK_I_INPUT_(i__), const double *d__, const double *z__, double *delta, SimTK_D_INPUT_(rho), SimTK_D_OUTPUT_(sigma), double *work, SimTK_INFO_);

void dlasd5_(SimTK_I_INPUT_(i__), const double *d__, const double *z__, double *delta, SimTK_D_INPUT_(rho), SimTK_D_OUTPUT_(dsigma), double *work);

void dlasd6_(SimTK_I_INPUT_(icompq), SimTK_I_INPUT_(nl), SimTK_I_INPUT_(nr), SimTK_I_INPUT_(sqre), double *d__, double *vf, double *vl, SimTK_D_INPUT_(alpha), SimTK_D_INPUT_(beta), int *idxq, int *perm, SimTK_I_OUTPUT_(givptr), SimTK_I_OUTPUT_(givcol), SimTK_FDIM_(ldgcol), double *givnum, SimTK_FDIM_(ldgnum), double *poles, double *difl, double *difr, double *z__, SimTK_I_OUTPUT_(k), SimTK_D_OUTPUT_(c__), SimTK_D_OUTPUT_(s), double *work, int *iwork, SimTK_INFO_);

void dlasd7_(SimTK_I_INPUT_(icompq), SimTK_I_INPUT_(nl), SimTK_I_INPUT_(nr), SimTK_I_INPUT_(sqre), SimTK_I_OUTPUT_(k), double *d__, double *z__, double *zw, double *vf, double *vfw, double *vl, double *vlw, SimTK_D_INPUT_(alpha), SimTK_D_INPUT_(beta), double *dsigma, int *idx, int *idxp, const int *idxq, int *perm, SimTK_I_OUTPUT_(givptr), int *givcol, SimTK_FDIM_(LDGCOL), double *givnum, SimTK_FDIM_(ldgnum), SimTK_D_OUTPUT_(c__), SimTK_D_OUTPUT_(s), SimTK_INFO_);

void dlasd8_(SimTK_I_INPUT_(icompq), SimTK_FDIM_(k), double *d__, const double *z__, double *vf, double *vl, double *difl, double *difr, SimTK_FDIM_(lddifr), const double *dsigma, double *work, SimTK_INFO_);

void dlasd9_(SimTK_I_INPUT_(icompq), SimTK_FDIM_(ldu), SimTK_FDIM_(k), double *d__, const double *z__, double *vf, double *vl, double *difl, double *difr, double *dsigma, double *work, SimTK_INFO_);

void dlasda_(SimTK_I_INPUT_(icompq), SimTK_I_INPUT_(smlsiz), SimTK_FDIM_(n), SimTK_I_INPUT_(sqre), double *d__, double *e, double *u, SimTK_FDIM_(ldu), double *vt, SimTK_FDIM_(k), double *difl, double *difr, double *z__, double *poles, int *givptr, int *givcol, SimTK_FDIM_(ldgcol), int *perm, double *givnum, double *c__, double *s, double *work, int *iwork, SimTK_INFO_);

void dlasdq_(SimTK_FOPT_(uplo), SimTK_I_INPUT_(sqre), SimTK_FDIM_(n), SimTK_I_INPUT_(ncvt), SimTK_I_INPUT_(nru), SimTK_I_INPUT_(ncc), double *d__, double *e, double *vt, SimTK_FDIM_(ldvt), double *u, SimTK_FDIM_(ldu), double *c__, SimTK_FDIM_(ldc), double *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void dlasdt_(SimTK_FDIM_(n), SimTK_I_OUTPUT_(lvl), SimTK_I_OUTPUT_(nd), int *inode, int *ndiml, int *ndimr, SimTK_I_INPUT_(msub) );

void dlaset_(SimTK_FOPT_(uplo), SimTK_FDIM_(m), SimTK_FDIM_(n),  SimTK_D_INPUT_(alpha), SimTK_D_INPUT_(beta), double *a, SimTK_FDIM_(lda), SimTK_FLEN_(uplo));

void dlasq1_(SimTK_FDIM_(n), double *d__, double *e, double *work, SimTK_INFO_);

void dlasq2_(SimTK_FDIM_(n), double *z__, SimTK_INFO_);

void dlasq3_(SimTK_I_INPUT_(i0), SimTK_I_INPUT_(n0), const double *z__, SimTK_I_INPUT_(pp), SimTK_D_OUTPUT_(dmin__), SimTK_D_OUTPUT_(sigma), SimTK_D_OUTPUT_(desig), SimTK_D_INPUT_(qmax), SimTK_I_OUTPUT_(nfail), SimTK_I_OUTPUT_(iter), SimTK_I_OUTPUT_(ndiv), SimTK_I_INPUT_(ieee)) ; 

void dlasq4_(SimTK_I_INPUT_(i0), SimTK_I_INPUT_(n0), const double *z__, SimTK_I_INPUT_(pp), SimTK_I_INPUT_(n0in), SimTK_D_INPUT_(dmin),  SimTK_D_INPUT_(dmin1), SimTK_D_INPUT_(dmin2), SimTK_D_INPUT_(dn), SimTK_D_INPUT_(dn1), SimTK_D_INPUT_(dn2), SimTK_D_OUTPUT_(tau), SimTK_I_OUTPUT_(ttype) );

void dlasq5_(SimTK_I_INPUT_(i0), SimTK_I_INPUT_(n0), const double *z__, SimTK_I_INPUT_(pp), SimTK_D_INPUT_(tau), SimTK_D_OUTPUT_(dmin), SimTK_D_OUTPUT_(dmin1), SimTK_D_OUTPUT_(dmin2), SimTK_D_OUTPUT_(dn), SimTK_D_OUTPUT_(dnm1), SimTK_D_OUTPUT_(dnm2), SimTK_I_INPUT_(ieee) );

void dlasq6_(SimTK_I_INPUT_(i0), SimTK_I_INPUT_(n0), const double *z__, SimTK_I_INPUT_(pp), SimTK_D_OUTPUT_(dmin), SimTK_D_OUTPUT_(dmin1), SimTK_D_OUTPUT_(dmin2), SimTK_D_OUTPUT_(dn), SimTK_D_OUTPUT_(dnm1), SimTK_D_OUTPUT_(dnm2));

void dlasr_(SimTK_FOPT_(side), SimTK_FOPT_(pivot), SimTK_FOPT_(direct), SimTK_FDIM_(m), SimTK_FDIM_(n), const double *c__, const double *s, double *a, SimTK_FDIM_(lda), SimTK_FLEN_(side), SimTK_FLEN_(pivot), SimTK_FLEN_(direct));

void dlasrt_(SimTK_FOPT_(id), SimTK_FDIM_(n), double *d__, SimTK_INFO_, SimTK_FLEN_(id));

void dlassq_(SimTK_FDIM_(n), const double *x, SimTK_FINC_(x), SimTK_D_OUTPUT_(scale), double *sumsq);

void dlasv2_( SimTK_D_INPUT_(f), SimTK_D_INPUT_(g), SimTK_D_INPUT_(h__), SimTK_D_OUTPUT_(ssmin), SimTK_D_OUTPUT_(ssmax) , SimTK_D_OUTPUT_(snr) , SimTK_D_OUTPUT_(csr) , SimTK_D_OUTPUT_(snl) , SimTK_D_OUTPUT_(csl)) ;

void dlaswp_(SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), SimTK_I_OUTPUT_(k1), SimTK_I_OUTPUT_(k2), const int *ipiv, SimTK_FINC_(x));

void dlasy2_(SimTK_I_INPUT_(ltranl), SimTK_I_INPUT_(ltranr), SimTK_I_INPUT_(isgn),  SimTK_I_INPUT_(n1), SimTK_I_INPUT_(n2), const double *tl, SimTK_FDIM_(ldtl), const double *tr, SimTK_FDIM_(ldtr), const double *b, SimTK_FDIM_(ldb), SimTK_D_OUTPUT_(scale), double *x, SimTK_FDIM_(ldx), double *xnorm, SimTK_INFO_);

void dlasyf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_I_INPUT_(nb), SimTK_I_INPUT_(kb), double *a, SimTK_FDIM_(lda), int *ipiv, double *w, SimTK_FDIM_(ldw), SimTK_INFO_, SimTK_FLEN_(uplo));

void dlatbs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FOPT_(normin), SimTK_FDIM_(n), SimTK_I_INPUT_(kd), const double *ab, SimTK_FDIM_(ldab), double *x, SimTK_D_OUTPUT_(scale), double *cnorm, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag), SimTK_FLEN_(normin));

void dlatdf_(SimTK_I_INPUT_(ijob), SimTK_FDIM_(n), const double *z__, SimTK_FDIM_(ldz), double *rhs, SimTK_D_OUTPUT_(rdsum), SimTK_D_OUTPUT_(rdscal), const int *ipiv, const int *jpiv);

void dlatps_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FOPT_(normin), SimTK_FDIM_(n), const double *ap, double *x, SimTK_D_OUTPUT_(scale), double *cnorm, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag), SimTK_FLEN_(normin));

void dlatrd_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nb), double *a, SimTK_FDIM_(lda), double *e, double *tau, double *w, SimTK_FDIM_(ldw), SimTK_FLEN_(uplo));

void dlatrs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FOPT_(normin), SimTK_FDIM_(n), const double *a, SimTK_FDIM_(lda), double *x, SimTK_D_OUTPUT_(scale), double *cnorm, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag), SimTK_FLEN_(normin) );

void dlatrz_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(l), double *a, SimTK_FDIM_(lda), double *tau, double *work);

void dlatzm_(SimTK_FOPT_(side), SimTK_FDIM_(m), SimTK_FDIM_(n), const double *v, SimTK_FINC_(v), SimTK_D_INPUT_(tau), double *c1, double *c2, SimTK_FDIM_(ldc), double *work, SimTK_FLEN_(side));

void dlauu2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo));

void dlauum_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo));

void dlazq3( SimTK_I_INPUT_(io), SimTK_I_INPUT_(no), const double *z, SimTK_I_INPUT_(pp), SimTK_D_OUTPUT_(dmin), SimTK_D_OUTPUT_(sigma), SimTK_D_OUTPUT_(desig), SimTK_D_INPUT_(qmax), SimTK_I_OUTPUT_(nfail), SimTK_I_OUTPUT_(iter), SimTK_I_OUTPUT_(ndiv), SimTK_I_INPUT_(ieee),  SimTK_I_INPUT_(ttype), SimTK_D_OUTPUT_(dmin1), SimTK_D_OUTPUT_(dmin2), SimTK_D_OUTPUT_(dn), SimTK_D_OUTPUT_(dn1), SimTK_D_OUTPUT_(dn2), SimTK_D_OUTPUT_(tau) );

void dlazq4( SimTK_I_INPUT_(io), SimTK_I_INPUT_(no), const double *z, SimTK_I_INPUT_(pp), SimTK_I_INPUT_(noin),  SimTK_D_INPUT_(dmin), SimTK_D_INPUT_(dmin1), SimTK_D_INPUT_(dmin2), SimTK_D_INPUT_(dn), SimTK_D_INPUT_(dn1), SimTK_D_INPUT_(dn2), SimTK_D_OUTPUT_(tau), SimTK_I_OUTPUT_(ttype), SimTK_D_OUTPUT_(g) );

void dopgtr_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *ap, double *tau, double *q, SimTK_FDIM_(ldq), double *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void dopmtr_(SimTK_FOPT_(side), SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), double *ap, double *tau, double *c__, SimTK_FDIM_(ldc), double *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(uplo), SimTK_FLEN_(trans));

void dorg2l_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), double *a, SimTK_FDIM_(lda), double *tau, double *work, SimTK_INFO_);

void dorg2r_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), double *a, SimTK_FDIM_(lda), double *tau, double *work, SimTK_INFO_);

void dorgbr_(SimTK_FOPT_(vect), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), double *a, SimTK_FDIM_(lda), double *tau, double *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(vect));

void dorghr_(SimTK_FDIM_(n), SimTK_FDIM_(ilo), SimTK_FDIM_(ihi), double *a, SimTK_FDIM_(lda), double *tau, double *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void dorgl2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), double *a, SimTK_FDIM_(lda), double *tau, double *work, SimTK_INFO_);

void dorglq_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), double *a, SimTK_FDIM_(lda), double *tau, double *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void dorgql_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), double *a, SimTK_FDIM_(lda), double *tau, double *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void dorgqr_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), double *a, SimTK_FDIM_(lda), double *tau, double *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void dorgr2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), double *a, SimTK_FDIM_(lda), double *tau, double *work, SimTK_INFO_);

void dorgrq_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), double *a, SimTK_FDIM_(lda), double *tau, double *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void dorgtr_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *tau, double *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(uplo));

void dorm2l_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), double *a, SimTK_FDIM_(lda), double *tau, double *c__, SimTK_FDIM_(ldc), double *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void dorm2r_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), double *a, SimTK_FDIM_(lda), double *tau, double *c__, SimTK_FDIM_(ldc), double *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void dormbr_(SimTK_FOPT_(vect), SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), double *a, SimTK_FDIM_(lda), double *tau, double *c__, SimTK_FDIM_(ldc), double *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(vect), SimTK_FLEN_(side), SimTK_FLEN_(trans));

void dormhr_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(ilo), SimTK_FDIM_(ihi), double *a, SimTK_FDIM_(lda), double *tau, double *c__, SimTK_FDIM_(ldc), double *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void dorml2_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), double *a, SimTK_FDIM_(lda), double *tau, double *c__, SimTK_FDIM_(ldc), double *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void dormlq_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), double *a, SimTK_FDIM_(lda), double *tau, double *c__, SimTK_FDIM_(ldc), double *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void dormql_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), double *a, SimTK_FDIM_(lda), double *tau, double *c__, SimTK_FDIM_(ldc), double *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void dormqr_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), double *a, SimTK_FDIM_(lda), double *tau, double *c__, SimTK_FDIM_(ldc), double *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void dormr2_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), double *a, SimTK_FDIM_(lda), double *tau, double *c__, SimTK_FDIM_(ldc), double *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void dormr3_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), int *l, double *a, SimTK_FDIM_(lda), double *tau, double *c__, SimTK_FDIM_(ldc), double *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void dormrq_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), double *a, SimTK_FDIM_(lda), double *tau, double *c__, SimTK_FDIM_(ldc), double *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void dormrz_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_FDIM_(l), double *a, SimTK_FDIM_(lda), double *tau, double *c__, SimTK_FDIM_(ldc), double *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void dormtr_(SimTK_FOPT_(side), SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *tau, double *c__, SimTK_FDIM_(ldc), double *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(uplo), SimTK_FLEN_(trans));

void dpbcon_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), const double *ab, SimTK_FDIM_(ldab), SimTK_D_INPUT_(anorm), SimTK_D_OUTPUT_(rcond), double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void dpbequ_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_I_INPUT_(kd), const double *ab, SimTK_FDIM_(ldab), SimTK_D_OUTPUT_(s), SimTK_D_OUTPUT_(scond), SimTK_D_OUTPUT_(amax), SimTK_INFO_, SimTK_FLEN_(uplo));

void dpbrfs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_FDIM_(nrhs), const double *ab, SimTK_FDIM_(ldab), const double *afb, SimTK_FDIM_(ldafb), const double *b, SimTK_FDIM_(ldb), double *x, SimTK_FDIM_(ldx), double *ferr, double *berr, double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void dpbstf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), double *ab, SimTK_FDIM_(ldab), SimTK_INFO_, SimTK_FLEN_(uplo));

void dpbsv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_FDIM_(nrhs), double *ab, SimTK_FDIM_(ldab), double *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void dpbsvx_(SimTK_FOPT_(fact), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_FDIM_(nrhs), double *ab, SimTK_FDIM_(ldab), double *afb, SimTK_FDIM_(ldafb), SimTK_CHAR_OUTPUT_(equed), double *s, double *b, SimTK_FDIM_(ldb), double *x, SimTK_FDIM_(ldx), SimTK_D_OUTPUT_(rcond), double *ferr, double *berr, double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(uplo), SimTK_FLEN_(equed));

void dpbtf2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), double *ab, SimTK_FDIM_(ldab), SimTK_INFO_, SimTK_FLEN_(uplo));

void dpbtrf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), double *ab, SimTK_FDIM_(ldab), SimTK_INFO_, SimTK_FLEN_(uplo));

void dpbtrs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_FDIM_(nrhs), const double *ab, SimTK_FDIM_(ldab), double *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void dpocon_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const double *a, SimTK_FDIM_(lda), SimTK_D_INPUT_(anorm), SimTK_D_OUTPUT_(rcond), double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void dpoequ_(SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *s, SimTK_D_OUTPUT_(scond), SimTK_D_OUTPUT_(amax), SimTK_INFO_);

void dporfs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const double *a, SimTK_FDIM_(lda), const double *af, SimTK_FDIM_(ldaf), const double *b, SimTK_FDIM_(ldb), double *x, SimTK_FDIM_(ldx), double *ferr, double *berr, double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void dposv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void dposvx_(SimTK_FOPT_(fact), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), double *a, SimTK_FDIM_(lda), double *af, SimTK_FDIM_(ldaf), SimTK_CHAR_OUTPUT_(equed), double *s, double *b, SimTK_FDIM_(ldb), double *x, SimTK_FDIM_(ldx), SimTK_D_OUTPUT_(rcond), double *ferr, double *berr, double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(uplo), SimTK_FLEN_(equed));

void dpotf2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo));

void dpotrf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo));

void dpotri_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo));

void dpotrs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void dppcon_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const double *ap, SimTK_D_INPUT_(anorm), SimTK_D_OUTPUT_(rcond), double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void dppequ_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const double *ap, double *s, SimTK_D_OUTPUT_(scond), SimTK_D_OUTPUT_(amax), SimTK_INFO_, SimTK_FLEN_(uplo));

void dpprfs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const double *ap, const double *afp, const double *b, SimTK_FDIM_(ldb), double *x, SimTK_FDIM_(ldx), double *ferr, double *berr, double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void dppsv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), double *ap, double *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void dppsvx_(SimTK_FOPT_(fact), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), double *ap, double *afp, SimTK_CHAR_OUTPUT_(equed), double *s, double *b, SimTK_FDIM_(ldb), double *x, SimTK_FDIM_(ldx), SimTK_D_OUTPUT_(rcond), double *ferr, double *berr, double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(uplo), SimTK_FLEN_(equed));

void dpptrf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *ap, SimTK_INFO_, SimTK_FLEN_(uplo));

void dpptri_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *ap, SimTK_INFO_, SimTK_FLEN_(uplo));

void dpptrs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const double *ap, double *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void dptcon_(SimTK_FDIM_(n), const double *d__, const double *e, SimTK_D_INPUT_(anorm), SimTK_D_OUTPUT_(rcond), double *work, SimTK_INFO_);

void dpteqr_(SimTK_FOPT_(compz), SimTK_FDIM_(n), double *d__, double *e, double *z__, SimTK_FDIM_(ldz), double *work, SimTK_INFO_, SimTK_FLEN_(compz));

void dptrfs_(SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const double *d__, const double *e, const double *df, const double *ef, const double *b, SimTK_FDIM_(ldb), double *x, SimTK_FDIM_(ldx), double *ferr, double *berr, double *work, SimTK_INFO_);

void dptsv_(SimTK_FDIM_(n), SimTK_FDIM_(nrhs), double *d__, double *e, double *b, SimTK_FDIM_(ldb), SimTK_INFO_);

void dptsvx_(SimTK_FOPT_(fact), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const double *d__, const double *e, double *df, double *ef, const double *b, SimTK_FDIM_(ldb), double *x, SimTK_FDIM_(ldx), SimTK_D_OUTPUT_(rcond), double *ferr, double *berr, double *work, SimTK_INFO_, SimTK_FLEN_(fact));

void dpttrf_(SimTK_FDIM_(n), double *d__, double *e, SimTK_INFO_);

void dpttrs_(SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const double *d__, const double *e, double *b, SimTK_FDIM_(ldb), SimTK_INFO_);

void dptts2_(SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const double *d__, const double *e, double *b, SimTK_FDIM_(ldb));

void drscl_(SimTK_FDIM_(n), double *sa, double *sx, SimTK_FINC_(x));

void dsbev_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), double *ab, SimTK_FDIM_(ldab), double *w, double *z__, SimTK_FDIM_(ldz), double *work, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void dsbevd_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), double *ab, SimTK_FDIM_(ldab), double *w, double *z__, SimTK_FDIM_(ldz), double *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void dsbevx_(SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), double *ab, SimTK_FDIM_(ldab), double *q, SimTK_FDIM_(ldq), SimTK_D_INPUT_(vl), SimTK_D_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_D_INPUT_(abstol), SimTK_I_OUTPUT_(m), double *w, double *z__, SimTK_FDIM_(ldz), double *work, int *iwork, int *ifail, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void dsbgst_(SimTK_FOPT_(vect), SimTK_FOPT_(uplo), SimTK_FDIM_(n), int *ka, int *kb, double *ab, SimTK_FDIM_(ldab), double *bb, SimTK_FDIM_(ldbb), double *x, SimTK_FDIM_(ldx), double *work, SimTK_INFO_, SimTK_FLEN_(vect), SimTK_FLEN_(uplo));

void dsbgv_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), int *ka, int *kb, double *ab, SimTK_FDIM_(ldab), double *bb, SimTK_FDIM_(ldbb), double *w, double *z__, SimTK_FDIM_(ldz), double *work, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void dsbgvd_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), int *ka, int *kb, double *ab, SimTK_FDIM_(ldab), double *bb, SimTK_FDIM_(ldbb), double *w, double *z__, SimTK_FDIM_(ldz), double *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void dsbgvx_(SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), int *ka, int *kb, double *ab, SimTK_FDIM_(ldab), double *bb, SimTK_FDIM_(ldbb), double *q, SimTK_FDIM_(ldq), SimTK_D_INPUT_(vl), SimTK_D_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_D_INPUT_(abstol), SimTK_I_OUTPUT_(m), double *w, double *z__, SimTK_FDIM_(ldz), double *work, int *iwork, int *ifail, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void dsbtrd_(SimTK_FOPT_(vect), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), double *ab, SimTK_FDIM_(ldab), double *d__, double *e, double *q, SimTK_FDIM_(ldq), double *work, SimTK_INFO_, SimTK_FLEN_(vect), SimTK_FLEN_(uplo));

void dsgesv_( SimTK_FDIM_(n), SimTK_FDIM_(nrhs), double *a, SimTK_FDIM_(lda), int *ipiv, const double *b, SimTK_FDIM_(ldb), double *x,  SimTK_FDIM_(ldx), double *work, float *swork, SimTK_I_OUTPUT_(iter), SimTK_INFO_ );

void dspcon_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const double *ap, const int *ipiv, SimTK_D_INPUT_(anorm), SimTK_D_OUTPUT_(rcond), double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void dspev_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *ap, double *w, double *z__, SimTK_FDIM_(ldz), double *work, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void dspevd_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *ap, double *w, double *z__, SimTK_FDIM_(ldz), double *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void dspevx_(SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *ap, SimTK_D_INPUT_(vl), SimTK_D_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_D_INPUT_(abstol), SimTK_I_OUTPUT_(m), double *w, double *z__, SimTK_FDIM_(ldz), double *work, int *iwork, int *ifail, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void dspgst_(SimTK_I_INPUT_(itype), SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *ap, double *bp, SimTK_INFO_, SimTK_FLEN_(uplo));

void dspgv_(SimTK_I_INPUT_(itype), SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *ap, double *bp, double *w, double *z__, SimTK_FDIM_(ldz), double *work, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void dspgvd_(SimTK_I_INPUT_(itype), SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *ap, double *bp, double *w, double *z__, SimTK_FDIM_(ldz), double *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void dspgvx_(SimTK_I_INPUT_(itype), SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *ap, double *bp, SimTK_D_INPUT_(vl), SimTK_D_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_D_INPUT_(abstol), SimTK_I_OUTPUT_(m), double *w, double *z__, SimTK_FDIM_(ldz), double *work, int *iwork, int *ifail, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void dsprfs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const double *ap, const double *afp, const int *ipiv, const double *b, SimTK_FDIM_(ldb), double *x, SimTK_FDIM_(ldx), double *ferr, double *berr, double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void dspsv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), double *ap, int *ipiv, double *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void dspsvx_(SimTK_FOPT_(fact), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const double *ap, double *afp, int *ipiv, const double *b, SimTK_FDIM_(ldb), double *x, SimTK_FDIM_(ldx), SimTK_D_OUTPUT_(rcond), double *ferr, double *berr, double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(face), SimTK_FLEN_(uplo));

void dsptrd_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *ap, double *d__, double *e, double *tau, SimTK_INFO_, SimTK_FLEN_(uplo));

void dsptrf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *ap, int *ipiv, SimTK_INFO_, SimTK_FLEN_(uplo));

void dsptri_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *ap, const int *ipiv, double *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void dsptrs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const double *ap, const int *ipiv, double *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void dstebz_(SimTK_FOPT_(range), SimTK_FOPT_(order), SimTK_FDIM_(n), SimTK_D_INPUT_(vl), SimTK_D_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_D_INPUT_(abstol), const double *d__, const double *e, SimTK_FDIM_(m), SimTK_I_OUTPUT_(nsplit), double *w, int *iblock, int *isplit, double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(range), SimTK_FLEN_(order));

void dstedc_(SimTK_FOPT_(compz), SimTK_FDIM_(n), double *d__, double *e, double *z__, SimTK_FDIM_(ldz), double *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(compz));

void dstegr_(SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FDIM_(n), double *d__, double *e, SimTK_D_INPUT_(vl), SimTK_D_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_D_INPUT_(abstol), SimTK_I_OUTPUT_(m), double *w, double *z__, SimTK_FDIM_(ldz), int *isuppz, double *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range));

void dstein_(SimTK_FDIM_(n), const double *d__, const double *e, SimTK_FDIM_(m), const double *w, const int *iblock, const int *isplit, double *z__, SimTK_FDIM_(ldz), double *work, int *iwork, int *ifail, SimTK_INFO_);

void dstemr_( SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FDIM_(n), double *d, double *e, SimTK_D_INPUT_(vl), SimTK_D_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_I_OUTPUT_(m), double *w, double *z, SimTK_FDIM_(ldz), SimTK_I_INPUT_(nzc), SimTK_I_OUTPUT_(isuppz), SimTK_I_OUTPUT_(tryrac), double *work, SimTK_FDIM_(lwork), int *iwork, SimTK_FDIM_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range) ); 

void dsteqr_(SimTK_FOPT_(compz), SimTK_FDIM_(n), double *d__, double *e, double *z__, SimTK_FDIM_(ldz), double *work, SimTK_INFO_, SimTK_FLEN_(compz));

void dsterf_(SimTK_FDIM_(n), double *d__, double *e, SimTK_INFO_);

void dstev_(SimTK_FOPT_(jobz), SimTK_FDIM_(n), double *d__, double *e, double *z__, SimTK_FDIM_(ldz), double *work, SimTK_INFO_, SimTK_FLEN_(jobz));

void dstevd_(SimTK_FOPT_(jobz), SimTK_FDIM_(n), double *d__, double *e, double *z__, SimTK_FDIM_(ldz), double *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz));

void dstevr_(SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FDIM_(n), double *d__, double *e, SimTK_D_INPUT_(vl), SimTK_D_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_D_INPUT_(abstol), SimTK_I_OUTPUT_(m), double *w, double *z__, SimTK_FDIM_(ldz), int *isuppz, double *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range));

void dstevx_(SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FDIM_(n), double *d__, double *e, SimTK_D_INPUT_(vl), SimTK_D_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_D_INPUT_(abstol), SimTK_I_OUTPUT_(m), double *w, double *z__, SimTK_FDIM_(ldz), double *work, int *iwork, int *ifail, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range));

void dsycon_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const double *a, SimTK_FDIM_(lda), const int *ipiv, SimTK_D_INPUT_(anorm), SimTK_D_OUTPUT_(rcond), double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void dsyev_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *w, double *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void dsyevd_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *w, double *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void dsyevr_(SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), SimTK_D_INPUT_(vl), SimTK_D_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_D_INPUT_(abstol), SimTK_I_OUTPUT_(m), double *w, double *z__, SimTK_FDIM_(ldz), int *isuppz, double *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void dsyevx_(SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), SimTK_D_INPUT_(vl), SimTK_D_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_D_INPUT_(abstol), SimTK_I_OUTPUT_(m), double *w, double *z__, SimTK_FDIM_(ldz), double *work, SimTK_FDIM_(lwork), int *iwork, int *ifail, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void dsygs2_(SimTK_I_INPUT_(itype), SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void dsygst_(SimTK_I_INPUT_(itype), SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void dsygv_(SimTK_I_INPUT_(itype), SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), double *w, double *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void dsygvd_(SimTK_I_INPUT_(itype), SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), double *w, double *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void dsygvx_(SimTK_I_INPUT_(itype), SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), SimTK_D_INPUT_(vl), SimTK_D_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_D_INPUT_(abstol), SimTK_I_OUTPUT_(m), double *w, double *z__, SimTK_FDIM_(ldz), double *work, SimTK_FDIM_(lwork), int *iwork, int *ifail, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void dsyrfs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const double *a, SimTK_FDIM_(lda), const double *af, SimTK_FDIM_(ldaf), const int *ipiv, const double *b, SimTK_FDIM_(ldb), double *x, SimTK_FDIM_(ldx), double *ferr, double *berr, double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void dsysv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), double *a, SimTK_FDIM_(lda), int *ipiv, double *b, SimTK_FDIM_(ldb), double *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(uplo));

void dsysvx_(SimTK_FOPT_(fact), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const double *a, SimTK_FDIM_(lda), double *af, SimTK_FDIM_(ldaf), int *ipiv, const double *b, SimTK_FDIM_(ldb), double *x, SimTK_FDIM_(ldx), SimTK_D_OUTPUT_(rcond), double *ferr, double *berr, double *work, SimTK_FDIM_(lwork), int *iwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(uplo));

void dsytd2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *d__, double *e, double *tau, SimTK_INFO_, SimTK_FLEN_(uplo));

void dsytf2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), int *ipiv, SimTK_INFO_, SimTK_FLEN_(uplo));

void dsytrd_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *d__, double *e, double *tau, double *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(uplo));

void dsytrf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), int *ipiv, double *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(uplo));

void dsytri_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const double *a, SimTK_FDIM_(lda), const int *ipiv, double *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void dsytrs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), double *a, SimTK_FDIM_(lda), int *ipiv, double *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void dtbcon_(SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(kd), const double *ab, SimTK_FDIM_(ldab), SimTK_D_OUTPUT_(rcond), double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(norm), SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

void dtbrfs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_FDIM_(nrhs), const double *ab, SimTK_FDIM_(ldab), const double *b, SimTK_FDIM_(ldb), const double *x, SimTK_FDIM_(ldx), double *ferr, double *berr, double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));

void dtbtrs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_FDIM_(nrhs), const double *ab, SimTK_FDIM_(ldab), const double *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));

void dtgevc_(SimTK_FOPT_(side), SimTK_FOPT_(howmny), const int *select, SimTK_FDIM_(n), const double *a, SimTK_FDIM_(lda), const double *b, SimTK_FDIM_(ldb), double *vl, SimTK_FDIM_(ldvl), double *vr, SimTK_FDIM_(ldvr), SimTK_I_INPUT_(mm), SimTK_FDIM_(m), double *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(howmny));

void dtgex2_(SimTK_I_INPUT_(wantq), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), double *q, SimTK_FDIM_(ldq), double *z__, SimTK_FDIM_(ldz), SimTK_I_INPUT_(j1),  SimTK_I_INPUT_(n1), SimTK_I_INPUT_(n2), double *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void dtgexc_(SimTK_I_INPUT_(wantq), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), double *q, SimTK_FDIM_(ldq), double *z__, SimTK_FDIM_(ldz), SimTK_I_OUTPUT_(ifst), SimTK_I_OUTPUT_(ilst), double *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void dtgsen_(SimTK_I_INPUT_(ijob), SimTK_I_INPUT_(wantq), SimTK_I_INPUT_(wantz), const int *select, SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), double *alphar, double *alphai, double *beta, double *q, SimTK_FDIM_(ldq), double *z__, SimTK_FDIM_(ldz), SimTK_I_OUTPUT_(m),  SimTK_D_OUTPUT_(pl),  SimTK_D_OUTPUT_(pr), double *dif, double *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_);

void dtgsja_(SimTK_FOPT_(jobu), SimTK_FOPT_(jobv), SimTK_FOPT_(jobq), SimTK_FDIM_(m), SimTK_FDIM_(p), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_FDIM_(l), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), SimTK_D_INPUT_(tola), SimTK_D_INPUT_(tolb), double *alpha, double *beta, double *u, SimTK_FDIM_(ldu), double *v, SimTK_FDIM_(ldv), double *q, SimTK_FDIM_(ldq), double *work, SimTK_I_OUTPUT_(ncycle), SimTK_INFO_, SimTK_FLEN_(jobu), SimTK_FLEN_(jobv), SimTK_FLEN_(jobq));

void dtgsna_(SimTK_FOPT_(job), SimTK_FOPT_(howmny), const int *select, SimTK_FDIM_(n), const double *a, SimTK_FDIM_(lda), const double *b, SimTK_FDIM_(ldb), const double *vl, SimTK_FDIM_(ldvl), const double *vr, SimTK_FDIM_(ldvr), double *s, double *dif, SimTK_FDIM_(mm), SimTK_FDIM_(m), double *work, SimTK_FDIM_(lwork), int *iwork, SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(howmny));

void dtgsy2_(SimTK_FOPT_(trans), SimTK_I_INPUT_(ijob), SimTK_FDIM_(m), SimTK_FDIM_(n), const double *a, SimTK_FDIM_(lda), const double *b, SimTK_FDIM_(ldb), double *c__, SimTK_FDIM_(ldc), const double *d__, SimTK_FDIM_(ldd), const double *e, SimTK_FDIM_(lde), double *f, SimTK_FDIM_(ldf), SimTK_D_OUTPUT_(scale), SimTK_D_OUTPUT_(rdsum), SimTK_D_OUTPUT_(rdscal), int *iwork, int *pq, SimTK_INFO_, SimTK_FLEN_(trans));

void dtgsyl_(SimTK_FOPT_(trans), SimTK_I_INPUT_(ijob), SimTK_FDIM_(m), SimTK_FDIM_(n), const double *a, SimTK_FDIM_(lda), const double *b, SimTK_FDIM_(ldb), double *c__, SimTK_FDIM_(ldc), const double *d__, SimTK_FDIM_(ldd), const double *e, SimTK_FDIM_(lde), double *f, SimTK_FDIM_(ldf), SimTK_D_OUTPUT_(scale), double *dif, double *work, SimTK_FDIM_(lwork), int *iwork, SimTK_INFO_, SimTK_FLEN_(trans));

void dtpcon_(SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), const double *ap, SimTK_D_OUTPUT_(rcond), double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));

void dtprfs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const double *ap, const double *b, SimTK_FDIM_(ldb), const double *x, SimTK_FDIM_(ldx), double *ferr, double *berr, double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));

void dtptri_(SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), double *ap, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

void dtptrs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const double *ap, double *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));

void dtrcon_(SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), const double *a, SimTK_FDIM_(lda), SimTK_D_OUTPUT_(rcond), double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(norm), SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

void dtrevc_(SimTK_FOPT_(side), SimTK_FOPT_(howmny), int *select, SimTK_FDIM_(n), double *t, SimTK_FDIM_(ldt), double *vl, SimTK_FDIM_(ldvl), double *vr, SimTK_FDIM_(ldvr), SimTK_FDIM_(mm), SimTK_FDIM_(m), double *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(howmny));

void dtrexc_(SimTK_FOPT_(compq), SimTK_FDIM_(n), double *t, SimTK_FDIM_(ldt), double *q, SimTK_FDIM_(ldq), SimTK_I_OUTPUT_(ifst), SimTK_I_OUTPUT_(ilst), double *work, SimTK_INFO_, SimTK_FLEN_(compq));

void dtrrfs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const double *a, SimTK_FDIM_(lda), const double *b, SimTK_FDIM_(ldb), const double *x, SimTK_FDIM_(ldx), double *ferr, double *berr, double *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));

void dtrsen_(SimTK_FOPT_(job), SimTK_FOPT_(compq), const int *select, SimTK_FDIM_(n), double *t, SimTK_FDIM_(ldt), double *q, SimTK_FDIM_(ldq), double *wr, double *wi, SimTK_FDIM_(m), SimTK_D_OUTPUT_(s), SimTK_D_OUTPUT_(sep), double *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(compq));

void dtrsna_(SimTK_FOPT_(job), SimTK_FOPT_(howmny), const int *select, SimTK_FDIM_(n), double *t, SimTK_FDIM_(ldt), double *vl, SimTK_FDIM_(ldvl), double *vr, SimTK_FDIM_(ldvr), double *s, double *sep, SimTK_FDIM_(mm), SimTK_FDIM_(m), double *work, SimTK_FDIM_(ldwork), int *iwork, SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(howmny));

void dtrsyl_(SimTK_FOPT_(trana), SimTK_FOPT_(tranb), int *isgn, SimTK_FDIM_(m), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), double *c__, SimTK_FDIM_(ldc), SimTK_D_OUTPUT_(scale), SimTK_INFO_, SimTK_FLEN_(trana), SimTK_FLEN_(tranb));

void dtrti2_(SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

void dtrtri_(SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans));

void dtrtrs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const double *a, SimTK_FDIM_(lda), double *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));

void dtzrqf_(SimTK_FDIM_(m), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *tau, SimTK_INFO_);

void dtzrzf_(SimTK_FDIM_(m), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *tau, double *work, SimTK_FDIM_(lwork), SimTK_INFO_);

int icmax1_(SimTK_FDIM_(n), const SimTK_C_ *cx, SimTK_FINC_(x));

int ieeeck_(SimTK_FDIM_(ispec), SimTK_S_INPUT_(zero), SimTK_S_INPUT_(one));

int ilaenv_(SimTK_FDIM_(ispec), const char *name__, const char *opts, SimTK_FDIM_(n1), SimTK_FDIM_(n2), SimTK_FDIM_(n3), SimTK_FDIM_(n4), SimTK_FLEN_(name), SimTK_FLEN_(opts));

void ilaver_( SimTK_I_OUTPUT_(vers_major), SimTK_I_OUTPUT_(vers_minor), SimTK_I_OUTPUT_(vers_patch) );
int iparmq_( SimTK_I_INPUT_(ispec), const char *name, const char *opts, SimTK_I_INPUT_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), SimTK_I_INPUT_(lwork), SimTK_FLEN_(name), SimTK_FLEN_(opts));

int izmax1_(SimTK_FDIM_(n), const SimTK_Z_ *cx, SimTK_FINC_(x));

void sbdsdc_(SimTK_FOPT_(uplo), SimTK_FOPT_(compq), SimTK_FDIM_(n), float *d__, float *e, float *u, SimTK_FDIM_(ldu), float *vt, SimTK_FDIM_(ldvt), float *q, int *iq, float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(compq));

void sbdsqr_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(ncvt), SimTK_FDIM_(nru), SimTK_FDIM_(ncc), float *d__, float *e, float *vt, SimTK_FDIM_(ldvt), float *u, SimTK_FDIM_(ldu), float *c__, SimTK_FDIM_(ldc), float *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void sdisna_(SimTK_FOPT_(job), SimTK_FDIM_(m), SimTK_FDIM_(n), const float *d__, float *sep, SimTK_INFO_, SimTK_FLEN_(job));

void sgbbrd_(SimTK_FOPT_(vect), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(ncc), SimTK_FDIM_(kl), SimTK_FDIM_(ku), float *ab, SimTK_FDIM_(ldab), float *d__, float *e, float *q, SimTK_FDIM_(ldq), float *pt, SimTK_FDIM_(ldpt), float *c__, SimTK_FDIM_(ldc), float *work, SimTK_INFO_, SimTK_FLEN_(vect));

void sgbcon_(SimTK_FOPT_(norm), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), float *ab, SimTK_FDIM_(ldab), const int *ipiv, SimTK_S_INPUT_(anorm), SimTK_S_OUTPUT_(rcond), float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(norm));

void sgbequ_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), float *ab, SimTK_FDIM_(ldab), float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, SimTK_INFO_);

void sgbrfs_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_FDIM_(nrhs), float *ab, SimTK_FDIM_(ldab), float *afb, SimTK_FDIM_(ldafb), const int *ipiv, float *b, SimTK_FDIM_(ldb), float *x, SimTK_FDIM_(ldx), float *ferr, float *berr, float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(trans));

void sgbsv_(SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_FDIM_(nrhs), float *ab, SimTK_FDIM_(ldab), int *ipiv, float *b, SimTK_FDIM_(ldb), SimTK_INFO_);

void sgbsvx_(SimTK_FOPT_(fact), SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_FDIM_(nrhs), float *ab, SimTK_FDIM_(ldab), float *afb, SimTK_FDIM_(ldafb), int *ipiv, SimTK_CHAR_OUTPUT_(equed), float *r__, float *c__, float *b, SimTK_FDIM_(ldb), float *x, SimTK_FDIM_(ldx), SimTK_S_OUTPUT_(rcond), float *ferr, float *berr, float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(trans), SimTK_FLEN_(equed));

void sgbtf2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), float *ab, SimTK_FDIM_(ldab), int *ipiv, SimTK_INFO_);

void sgbtrf_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), float *ab, SimTK_FDIM_(ldab), int *ipiv, SimTK_INFO_);

void sgbtrs_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_FDIM_(nrhs), float *ab, SimTK_FDIM_(ldab), const int *ipiv, float *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(trans));

void sgebak_(SimTK_FOPT_(job), SimTK_FOPT_(side), SimTK_FDIM_(n),  SimTK_I_INPUT_(ilo),  SimTK_I_INPUT_(ihi),  SimTK_S_INPUT_(scale),  SimTK_I_INPUT_(m), float *v, SimTK_FDIM_(ldv), SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(side));

void sgebal_(SimTK_FOPT_(job), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), SimTK_I_OUTPUT_(ilo), SimTK_I_OUTPUT_(ihi), SimTK_S_OUTPUT_(scale), SimTK_INFO_, SimTK_FLEN_(job));

void sgebd2_(SimTK_FDIM_(m), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *d__, float *e, float *tauq, float *taup, float *work, SimTK_INFO_);

void sgebrd_(SimTK_FDIM_(m), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *d__, float *e, float *tauq, float *taup, float *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void sgecon_(SimTK_FOPT_(norm), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), SimTK_S_INPUT_(anorm), SimTK_S_OUTPUT_(rcond), float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(norm));

void sgeequ_(SimTK_FDIM_(m), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, int *info);

void sgees_(SimTK_FOPT_(jobvs), SimTK_FOPT_(sort), SimTK_SELECT_2S select, SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), int *sdim, float *wr, float *wi, float *vs, SimTK_FDIM_(ldvs), float *work, SimTK_FDIM_(lwork), int *bwork, SimTK_INFO_, SimTK_FLEN_(jobvs), SimTK_FLEN_(sort));

void sgeesx_(SimTK_FOPT_(jobvs), SimTK_FOPT_(sort), SimTK_SELECT_2S select, char *sense, SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), int *sdim, float *wr, float *wi, float *vs, SimTK_FDIM_(ldvs), SimTK_S_OUTPUT_(rconde), SimTK_S_OUTPUT_(rcondv), float *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), int *bwork, SimTK_INFO_, SimTK_FLEN_(jobvs), SimTK_FLEN_(sort), SimTK_FLEN_(sense));

void sgeev_(SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *wr, float *wi, float *vl, SimTK_FDIM_(ldvl), float *vr, SimTK_FDIM_(ldvr), float *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr));

void sgeevx_(SimTK_FOPT_(balanc), SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), char *sense, SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *wr, float *wi, float *vl, SimTK_FDIM_(ldvl), float *vr, SimTK_FDIM_(ldvr), SimTK_I_OUTPUT_(ilo), SimTK_I_OUTPUT_(ihi), SimTK_S_OUTPUT_(scale), float *abnrm, SimTK_S_OUTPUT_(rconde), SimTK_S_OUTPUT_(rcondv), float *work, SimTK_FDIM_(lwork), int *iwork, SimTK_INFO_, SimTK_FLEN_(balanc), SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr), SimTK_FLEN_(alpahr));

void sgegs_(SimTK_FOPT_(jobvsl), SimTK_FOPT_(jobvsr), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), float *alphar, float *alphai, float *beta, float *vsl, SimTK_FDIM_(ldvsl), float *vsr, SimTK_FDIM_(ldvsr), float *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvsr));

void sgegv_(SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), float *alphar, float *alphai, float *beta, float *vl, SimTK_FDIM_(ldvl), float *vr, SimTK_FDIM_(ldvr), float *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr));

void sgehd2_(SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), float *a, SimTK_FDIM_(lda), float *tau, float *work, SimTK_INFO_);

void sgehrd_(SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), float *a, SimTK_FDIM_(lda), float *tau, float *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void sgelq2_(SimTK_FDIM_(m), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *tau, float *work, SimTK_INFO_);

void sgelqf_(SimTK_FDIM_(m), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *tau, float *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void sgels_(SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), float *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(trans));

void sgelsd_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), float *s, SimTK_S_INPUT_(rcond), SimTK_I_OUTPUT_(rank), float *work, SimTK_FDIM_(lwork), int *iwork, SimTK_INFO_);

void sgelss_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), float *s, SimTK_S_INPUT_(rcond), SimTK_I_OUTPUT_(rank), float *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void sgelsx_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), int *jpvt, SimTK_S_INPUT_(rcond), SimTK_I_OUTPUT_(rank), float *work, SimTK_INFO_);

void sgelsy_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), int *jpvt, SimTK_S_INPUT_(rcond), SimTK_I_OUTPUT_(rank), float *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void sgeql2_(SimTK_FDIM_(m), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *tau, float *work, SimTK_INFO_);

void sgeqlf_(SimTK_FDIM_(m), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *tau, float *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void sgeqp3_(SimTK_FDIM_(m), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), int *jpvt, float *tau, float *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void sgeqpf_(SimTK_FDIM_(m), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), int *jpvt, float *tau, float *work, SimTK_INFO_);

void sgeqr2_(SimTK_FDIM_(m), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *tau, float *work, SimTK_INFO_);

void sgeqrf_(SimTK_FDIM_(m), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *tau, float *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void sgerfs_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), float *a, SimTK_FDIM_(lda), float *af, SimTK_FDIM_(ldaf), const int *ipiv, float *b, SimTK_FDIM_(ldb), float *x, SimTK_FDIM_(ldx), float *ferr, float *berr, float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(trans));

void sgerq2_(SimTK_FDIM_(m), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *tau, float *work, SimTK_INFO_);

void sgerqf_(SimTK_FDIM_(m), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *tau, float *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void sgesc2_(SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *rhs, int *ipiv, int *jpiv, SimTK_S_OUTPUT_(scale));

void sgesdd_(SimTK_FOPT_(jobz), SimTK_FDIM_(m), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *s, float *u, SimTK_FDIM_(ldu), float *vt, SimTK_FDIM_(ldvt), float *work, SimTK_FDIM_(lwork), int *iwork, SimTK_INFO_, SimTK_FLEN_(jobz));

void sgesv_(SimTK_FDIM_(n), SimTK_FDIM_(nrhs), float *a, SimTK_FDIM_(lda), int *ipiv, float *b, SimTK_FDIM_(ldb), SimTK_INFO_);

void sgesvd_(SimTK_FOPT_(jobu), char *jobvt, SimTK_FDIM_(m), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *s, float *u, SimTK_FDIM_(ldu), float *vt, SimTK_FDIM_(ldvt), float *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(jobu), SimTK_FLEN_(jobvt));

void sgesvx_(SimTK_FOPT_(fact), SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), float *a, SimTK_FDIM_(lda), float *af, SimTK_FDIM_(ldaf), int *ipiv, SimTK_CHAR_OUTPUT_(equed), float *r__, float *c__, float *b, SimTK_FDIM_(ldb), float *x, SimTK_FDIM_(ldx), SimTK_S_OUTPUT_(rcond), float *ferr, float *berr, float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(trans), SimTK_FLEN_(equed));

void sgetc2_(SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), int *ipiv, int *jpiv, SimTK_INFO_);

void sgetf2_(SimTK_FDIM_(m), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), int *ipiv, SimTK_INFO_);

void sgetrf_(SimTK_FDIM_(m), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), int *ipiv, SimTK_INFO_);

void sgetri_(SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), const int *ipiv, float *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void sgetrs_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const float *a, SimTK_FDIM_(lda), const int *ipiv, float *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(trans));

void sggbak_(SimTK_FOPT_(job), SimTK_FOPT_(side), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), const float *lscale, const float *rscale, SimTK_FDIM_(m), float *v, SimTK_FDIM_(ldv), SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(side));

void sggbal_(SimTK_FOPT_(job), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), SimTK_I_OUTPUT_(ilo), SimTK_I_OUTPUT_(ihi), float *lscale, float *rscale, float *work, SimTK_INFO_, SimTK_FLEN_(job));

void sgges_(SimTK_FOPT_(jobvsl), SimTK_FOPT_(jobvsr), SimTK_FOPT_(sort), SimTK_SELECT_3F selctg, SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), int *sdim, float *alphar, float *alphai, SimTK_S_INPUT_(eta), float *vsl, SimTK_FDIM_(ldvsl), float *vsr, SimTK_FDIM_(ldvsr), float *work, SimTK_FDIM_(lwork), int *bwork, SimTK_INFO_, SimTK_FLEN_(jobvsl), SimTK_FLEN_(jobvsr), SimTK_FLEN_(sort));

void sggesx_(SimTK_FOPT_(jobvsl), SimTK_FOPT_(jobvsr), SimTK_FOPT_(sort), SimTK_SELECT_3F selctg, SimTK_FOPT_(sense), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), int *sdim, float *alphar, float *alphai, float *beta, float *vsl, SimTK_FDIM_(ldvsl), float *vsr, SimTK_FDIM_(ldvsr), SimTK_S_OUTPUT_(rconde), SimTK_S_OUTPUT_(rcondv), float *work, SimTK_FDIM_(lwork), int *iwork, SimTK_FDIM_(liwork), int *bwork, SimTK_INFO_, SimTK_FLEN_(jobvsl), SimTK_FLEN_(jobvsr), SimTK_FLEN_(sort), SimTK_FLEN_(sense));

void sggev_(SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), float *alphar, float *alphai, float *beta, float *vl, SimTK_FDIM_(ldvl), float *vr, SimTK_FDIM_(ldvr), float *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(jobvsl), SimTK_FLEN_(jobvsr));

void sggevx_(SimTK_FOPT_(balanc), SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), SimTK_FOPT_(sense), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), float *alphar, float *alphai, float *beta, float *vl, SimTK_FDIM_(ldvl), float *vr, SimTK_FDIM_(ldvr), SimTK_I_OUTPUT_(ilo), SimTK_I_OUTPUT_(ihi), float *lscale, float *rscale, SimTK_S_OUTPUT_(abnrm), SimTK_S_OUTPUT_(bbnrm), SimTK_S_OUTPUT_(rconde), SimTK_S_OUTPUT_(rcondv), float *work, SimTK_FDIM_(lwork), int *iwork, int *bwork, SimTK_INFO_, SimTK_FLEN_(balanc), SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr), SimTK_FLEN_(alpahr));

void sggglm_(SimTK_FDIM_(n), SimTK_FDIM_(m), SimTK_FDIM_(p), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), float *d__, float *x, float *y, float *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void sgghrd_(SimTK_FOPT_(compq), SimTK_FOPT_(compz), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), float *q, SimTK_FDIM_(ldq), float *z__, SimTK_FDIM_(ldz), SimTK_INFO_, SimTK_FLEN_(compq), SimTK_FLEN_(compz));

void sgglse_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(p), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), float *c__, float *d__, float *x, float *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void sggqrf_(SimTK_FDIM_(n), SimTK_FDIM_(m), SimTK_FDIM_(p), float *a, SimTK_FDIM_(lda), float *taua, float *b, SimTK_FDIM_(ldb), float *taub, float *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void sggrqf_(SimTK_FDIM_(m), SimTK_FDIM_(p), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *taua, float *b, SimTK_FDIM_(ldb), float *taub, float *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void sggsvd_(SimTK_FOPT_(jobu), SimTK_FOPT_(jobv), SimTK_FOPT_(jobq), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(p), SimTK_I_OUTPUT_(k), SimTK_I_OUTPUT_(l), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), float *alpha, float *beta, float *u, SimTK_FDIM_(ldu), float *v, SimTK_FDIM_(ldv), float *q, SimTK_FDIM_(ldq), float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(jobu), SimTK_FLEN_(jobv), SimTK_FLEN_(jobq));

void sggsvp_(SimTK_FOPT_(jobu), SimTK_FOPT_(jobv), SimTK_FOPT_(jobq), SimTK_FDIM_(m), SimTK_FDIM_(p), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), SimTK_S_INPUT_(tola), SimTK_S_INPUT_(tolb),  SimTK_I_OUTPUT_(k), SimTK_I_OUTPUT_(l), float *u, SimTK_FDIM_(ldu), float *v, SimTK_FDIM_(ldv), float *q, SimTK_FDIM_(ldq), int *iwork, float *tau, float *work, SimTK_INFO_, SimTK_FLEN_(jobu), SimTK_FLEN_(jobv), SimTK_FLEN_(jobq));

void sgtcon_(SimTK_FOPT_(norm), SimTK_FDIM_(n), float *dl, float *d__, float *du, float *du2, const int *ipiv, SimTK_S_INPUT_(anorm), SimTK_S_OUTPUT_(rcond), float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(norm));

void sgtrfs_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), float *dl, float *d__, float *du, float *dlf, float *df, float *duf, float *du2, const int *ipiv, float *b, SimTK_FDIM_(ldb), float *x, SimTK_FDIM_(ldx), float *ferr, float *berr, float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(trans));

void sgtsv_(SimTK_FDIM_(n), SimTK_FDIM_(nrhs), float *dl, float *d__, float *du, float *b, SimTK_FDIM_(ldb), SimTK_INFO_);

void sgtsvx_(SimTK_FOPT_(fact), SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), float *dl, float *d__, float *du, float *dlf, float *df, float *duf, float *du2, int *ipiv, float *b, SimTK_FDIM_(ldb), float *x, SimTK_FDIM_(ldx), SimTK_S_OUTPUT_(rcond), float *ferr, float *berr, float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(trans));

void sgttrf_(SimTK_FDIM_(n), float *dl, float *d__, float *du, float *du2, int *ipiv, SimTK_INFO_);

void sgttrs_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), float *dl, float *d__, float *du, float *du2, const int *ipiv, float *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(trans));

void sgtts2_(int *itrans, SimTK_FDIM_(n), SimTK_FDIM_(nrhs), float *dl, float *d__, float *du, float *du2, const int *ipiv, float *b, SimTK_FDIM_(ldb));

void shgeqz_(SimTK_FOPT_(job), SimTK_FOPT_(compq), SimTK_FOPT_(compz), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), float *alphar, float *alphai, float *beta, float *q, SimTK_FDIM_(ldq), float *z__, SimTK_FDIM_(ldz), float *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(compq), SimTK_FLEN_(compz));

void shsein_(SimTK_FOPT_(side), SimTK_FOPT_(eigsrc), SimTK_FOPT_(initv), int *select, SimTK_FDIM_(n), const float *h__, SimTK_FDIM_(ldh), float *wr, float *wi, float *vl, SimTK_FDIM_(ldvl), float *vr, SimTK_FDIM_(ldvr), SimTK_FDIM_(mm), SimTK_I_OUTPUT_(m), float *work, int *ifaill, int *ifailr, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(eigsrc), SimTK_FLEN_(initv));

float  slamch_(SimTK_FOPT_(cmach), SimTK_FLEN_(cmach));

void shseqr_(SimTK_FOPT_(job), SimTK_FOPT_(compz), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), float *h__, SimTK_FDIM_(ldh), float *wr, float *wi, float *z__, SimTK_FDIM_(ldz), float *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(compz));

int sisnan( SimTK_S_INPUT_(sin) );

int slaisnan( SimTK_S_INPUT_(sin1), SimTK_S_INPUT_(sin2) );

void slabad_(SimTK_S_INOUT_(small1), SimTK_S_INOUT_(large));//small -> small1 conflict with definition when Qt is used

void slabrd_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(nb), float *a, SimTK_FDIM_(lda), float *d__, float *e, float *tauq, float *taup, float *x, SimTK_FDIM_(ldx), float *y, SimTK_FDIM_(ldy));

void slacon_(SimTK_FDIM_(n), float *v, float *x, int *isgn, SimTK_S_OUTPUT_(est), SimTK_I_OUTPUT_(kase) );

void slacn2_( SimTK_FDIM_(n), float *v, float *x, int *isgn, SimTK_S_INPUT_(est), SimTK_I_OUTPUT_(kase), int *isave );

void slacpy_(SimTK_FOPT_(uplo), SimTK_FDIM_(m), SimTK_FDIM_(n), const float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), SimTK_FLEN_(uplo));

void sladiv_(SimTK_S_INPUT_(a), SimTK_S_INPUT_(b), SimTK_S_INPUT_(c__), SimTK_S_INPUT_(d__), float *p, float *q);

void slae2_(SimTK_S_INPUT_(a), SimTK_S_INPUT_(b), SimTK_S_INPUT_(c__), float *rt1, float *rt2);

void slaebz_(SimTK_I_INPUT_(ijob), SimTK_I_INPUT_(nitmax), SimTK_FDIM_(n), SimTK_I_INPUT_(mmax), SimTK_I_INPUT_(minp), SimTK_I_INPUT_(nbmin), SimTK_S_INPUT_(abstol), SimTK_S_INPUT_(reltol), SimTK_S_INPUT_(pivmin), const float *d__, const float *e, const float *e2, int *nval, float *ab, float *c__, SimTK_I_OUTPUT_(mout), int *nab, float *work, int *iwork, SimTK_INFO_);

void slaed0_(SimTK_I_INPUT_(icompq), SimTK_I_INPUT_(qsiz), SimTK_FDIM_(n), float *d__, float *e, float *q, SimTK_FDIM_(ldq), float *qstore, SimTK_FDIM_(ldqs), float *work, int *iwork, SimTK_INFO_);

void slaed1_(SimTK_FDIM_(n), float *d__, float *q, SimTK_FDIM_(ldq), int *indxq, SimTK_S_INPUT_(rho), SimTK_I_INPUT_(cutpnt), float *work, int *iwork, SimTK_INFO_);

void slaed2_(SimTK_FDIM_(k), SimTK_FDIM_(n),  SimTK_I_INPUT_(n1), float *d__, float *q, SimTK_FDIM_(ldq), int *indxq, SimTK_S_INPUT_(rho), float *z__, float *dlamda, float *w, float *q2, int *indx, int *indxc, int *indxp, int *coltyp, SimTK_INFO_);

void slaed3_(SimTK_FDIM_(k), SimTK_FDIM_(n),  SimTK_I_INPUT_(n1), float *d__, float *q, SimTK_FDIM_(ldq), SimTK_S_INPUT_(rho), float *dlamda, const float *q2, const int *indx, const int *ctot, float *w, float *s, SimTK_INFO_);

void slaed4_(SimTK_FDIM_(n), SimTK_I_INPUT_(i__), const float *d__, const float *z__, float *delta, SimTK_S_INPUT_(rho), SimTK_S_OUTPUT_(dlam), SimTK_INFO_);

void slaed5_(SimTK_I_INPUT_(i__), const float *d__, const float *z__, float *delta, SimTK_S_INPUT_(rho), SimTK_S_OUTPUT_(dlam));

void slaed6_(SimTK_I_INPUT_(kniter), SimTK_I_INPUT_(orgati), SimTK_S_INPUT_(rho), const float *d__, const float *z__, SimTK_S_INPUT_(finit), SimTK_S_INPUT_(tau), SimTK_INFO_);

void slaed7_(SimTK_I_INPUT_(icompq), SimTK_FDIM_(n), SimTK_FDIM_(qsiz), SimTK_I_INPUT_(tlvls), SimTK_I_INPUT_(curlvl), SimTK_I_INPUT_(curpbm), float *d__, float *q, SimTK_FDIM_(ldq), int *indxq, SimTK_S_INPUT_(rho), SimTK_I_INPUT_(cutpnt), float *qstore, int *qptr, const int *prmptr, const int *perm, const int *givptr, const int *givcol, const float *givnum, float *work, int *iwork, SimTK_INFO_);

void slaed8_(SimTK_I_INPUT_(icompq), SimTK_I_OUTPUT_(k), SimTK_FDIM_(n), SimTK_FDIM_(qsiz), float *d__, float *q, SimTK_FDIM_(ldq), const int *indxq, SimTK_S_OUTPUT_(rho), SimTK_I_INPUT_(cutpnt), float *z__, float *dlamda, float *q2, SimTK_FDIM_(ldq2), float *w, int *perm, int *givptr, int *givcol, float *givnum, int *indxp, int *indx, SimTK_INFO_);

void slaed9_(SimTK_FDIM_(k), int *kstart, int *kstop, SimTK_FDIM_(n), float *d__, float *q, SimTK_FDIM_(ldq), float *rho, float *dlamda, float *w, float *s, SimTK_FDIM_(lds), SimTK_INFO_);

void slaeda_(SimTK_FDIM_(n), SimTK_I_INPUT_(tlvls), SimTK_I_INPUT_(curlvl), SimTK_I_INPUT_(curpbm), const int *prmptr, const int *perm, const int *givptr, const int *givcol, const float *givnum, const float *q, const int *qptr, float *z__, float *ztemp, SimTK_INFO_);

void slaein_( SimTK_I_INPUT_(rightv),  SimTK_I_INPUT_(noinit), SimTK_FDIM_(n), const float *h__, SimTK_FDIM_(ldh), SimTK_S_INPUT_(wr), SimTK_S_INPUT_(wi), float *vr, float *vi, float *b, SimTK_FDIM_(ldb), float *work, SimTK_S_INPUT_(eps3), SimTK_S_INPUT_(smlnum), SimTK_S_INPUT_(bignum), SimTK_INFO_);

void slaev2_(SimTK_S_INPUT_(a), SimTK_S_INPUT_(b), SimTK_S_INPUT_(c__), SimTK_S_OUTPUT_(rt1), SimTK_S_OUTPUT_(rt2), SimTK_S_OUTPUT_(cs1), SimTK_S_OUTPUT_(sn1));

void slaexc_(SimTK_I_INPUT_(wantq), SimTK_FDIM_(n), float *t, SimTK_FDIM_(ldt), float *q, SimTK_FDIM_(ldq), SimTK_I_INPUT_(j1),  SimTK_I_INPUT_(n1), SimTK_I_INPUT_(n2), float *work, SimTK_INFO_);

void slag2_(const float *a, SimTK_FDIM_(lda), const float *b, SimTK_FDIM_(ldb), SimTK_S_INPUT_(safmin), SimTK_S_OUTPUT_(scale1), SimTK_S_OUTPUT_(scale2), SimTK_S_OUTPUT_(wr1), SimTK_S_OUTPUT_(wr2), SimTK_S_OUTPUT_(wi)); 

void slags2_(SimTK_I_INPUT_(upper), SimTK_S_INPUT_(a1), SimTK_S_INPUT_(a2), SimTK_S_INPUT_(a3), SimTK_S_INPUT_(b1), SimTK_S_INPUT_(b2), SimTK_S_INPUT_(b3), SimTK_S_OUTPUT_(csu), SimTK_S_OUTPUT_(snu), SimTK_S_OUTPUT_(csv), SimTK_S_OUTPUT_(snv), SimTK_S_OUTPUT_(csq), SimTK_S_OUTPUT_(snq) );

void slagtf_(SimTK_FDIM_(n), float *a, SimTK_S_INPUT_(lambda), float *b, float *c__, SimTK_S_INPUT_(tol), float *d__, int *in, SimTK_INFO_);

void slagtm_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_S_INPUT_(alpha), const float *dl, const float *d__, const float *du, const float *x, SimTK_FDIM_(ldx), SimTK_S_INPUT_(beta), float *b, SimTK_FDIM_(ldb), SimTK_FLEN_(trans));

void slag2d_( SimTK_I_INPUT_(m), SimTK_I_INPUT_(n),  float *sa, SimTK_FDIM_(ldsa), const double *a, SimTK_FDIM_(lda),  SimTK_INFO_);

void slagts_(SimTK_I_INPUT_(job), SimTK_FDIM_(n), const float *a, const float *b, const float *c__, const float *d__, const int *in, float *y, SimTK_S_OUTPUT_(tol), SimTK_INFO_);

void slagv2_( float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), float *alphar, float *alphai, float *beta, SimTK_S_OUTPUT_(csl), SimTK_S_OUTPUT_(snl), SimTK_S_OUTPUT_(csr), SimTK_S_OUTPUT_(snr));

void slahqr_(SimTK_I_INPUT_(wantt), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), float *h__, SimTK_FDIM_(ldh), float *wr, float *wi, SimTK_I_INPUT_(iloz), SimTK_I_INPUT_(ihiz), float *z__, SimTK_FDIM_(ldz), SimTK_INFO_);

void slahrd_(SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_FDIM_(nb), float *a, SimTK_FDIM_(lda), float *tau, float *t, SimTK_FDIM_(ldt), float *y, SimTK_FDIM_(ldy));

void slahr2_(SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_FDIM_(nb), float *a, SimTK_FDIM_(lda), float *tau, float *t, SimTK_FDIM_(ldt), float *y, SimTK_FDIM_(ldy));

void slaic1_(SimTK_I_INPUT_(job), SimTK_FDIM_(j), const float *x, SimTK_S_INPUT_(sest), const float *w, SimTK_S_INPUT_(gamma), SimTK_S_OUTPUT_(sestpr), SimTK_S_OUTPUT_(s), SimTK_S_OUTPUT_(c__) ); 

void slaln2_(SimTK_I_INPUT_(ltrans), SimTK_I_INPUT_(na), SimTK_I_INPUT_(nw), SimTK_S_INPUT_(smin), SimTK_S_INPUT_(ca), const float *a, SimTK_FDIM_(lda), SimTK_S_INPUT_(d1), SimTK_S_INPUT_(d2), SimTK_S_INPUT_(b), SimTK_FDIM_(ldb), SimTK_S_INPUT_(wr), SimTK_S_INPUT_(wi), float *x, SimTK_FDIM_(ldx), SimTK_S_OUTPUT_(scale), float *xnorm, SimTK_INFO_);

void slals0_(SimTK_I_INPUT_(icompq), SimTK_I_INPUT_(nl), SimTK_I_INPUT_(nr), SimTK_I_INPUT_(sqre), SimTK_FDIM_(nrhs), float *b, SimTK_FDIM_(ldb), float *bx, SimTK_FDIM_(ldbx), int *perm, int *givptr, const int *givcol, SimTK_FDIM_(ldgcol), float *givnum, SimTK_FDIM_(ldgnum), float *poles, float *difl, float *difr, float *z__, SimTK_FDIM_(k), float *c__, float *s, float *work, SimTK_INFO_);

void slalsa_(SimTK_I_INPUT_(icompq), SimTK_I_INPUT_(smlsiz), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), float *b, SimTK_FDIM_(ldb), float *bx, SimTK_FDIM_(ldbx), float *u, SimTK_FDIM_(ldu), float *vt, SimTK_FDIM_(k), float *difl, float *difr, float *z__, float *poles, SimTK_FDIM_(givptr), const int *givcol, SimTK_FDIM_(ldgcol), int *perm, float *givnum, float *c__, float *s, float *work, int *iwork, SimTK_INFO_);

void slalsd_(SimTK_FOPT_(uplo), SimTK_I_INPUT_(smlsiz), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), float *d__, float *e, float *b, SimTK_FDIM_(ldb), SimTK_S_INPUT_(rcond), SimTK_I_OUTPUT_(rank), float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void slaneg_(  SimTK_FDIM_(n), SimTK_S_INPUT_(d), SimTK_S_INPUT_(lld), SimTK_S_INPUT_(sigma), SimTK_S_INPUT_(pivmin), SimTK_I_INPUT_(r) );

double slangb_( SimTK_FOPT_(norm), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), const float *ab, SimTK_FDIM_(ldab), float *work, SimTK_FLEN_(norm) );

double clangb_( SimTK_FOPT_(norm), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), const SimTK_C_ *ab, SimTK_FDIM_(ldab), float *work, SimTK_FLEN_(norm) );

void dlaneg_(  SimTK_FDIM_(n), SimTK_D_INPUT_(d), SimTK_D_INPUT_(lld), SimTK_D_INPUT_(sigma), SimTK_D_INPUT_(pivmin), SimTK_I_INPUT_(r) );

double dlangb_( SimTK_FOPT_(norm), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), const double *ab, SimTK_FDIM_(ldab), double *work, SimTK_FLEN_(norm) );

double zlangb_( SimTK_FOPT_(norm), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), const SimTK_Z_ *ab, SimTK_FDIM_(ldab), double *work, SimTK_FLEN_(norm) );

double slange_( SimTK_FOPT_(norm), SimTK_FDIM_(m), SimTK_FDIM_(n), const float *a,  SimTK_FDIM_(lda), float *work, SimTK_FLEN_(norm));
double clange_( SimTK_FOPT_(norm), SimTK_FDIM_(m), SimTK_FDIM_(n), const SimTK_C_ *a,  SimTK_FDIM_(lda), float *work, SimTK_FLEN_(norm));
double dlange_( SimTK_FOPT_(norm), SimTK_FDIM_(m), SimTK_FDIM_(n), const double *a,  SimTK_FDIM_(lda), double *work, SimTK_FLEN_(norm));
double zlange_( SimTK_FOPT_(norm), SimTK_FDIM_(m), SimTK_FDIM_(n), const SimTK_Z_ *a,  SimTK_FDIM_(lda), double *work, SimTK_FLEN_(norm));

double slansb_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(k), const float *ab,  SimTK_FDIM_(ldab), float *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo));

double dlansb_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(k), const double *ab,  SimTK_FDIM_(ldab), double *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo));

double clansb_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_C_ *ab,  SimTK_FDIM_(ldab), float *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo));

double zlansb_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_Z_ *ab,  SimTK_FDIM_(ldab), double *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo));


double clanhb_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_C_ *ab,  SimTK_FDIM_(ldab), float *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo));

double zlanhb_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_Z_ *ab,  SimTK_FDIM_(ldab), double *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo));

double slansp_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FDIM_(n), const float *ap,  float *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo));
double dlansp_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FDIM_(n), const double *ap,  double *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo));
double clansp_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FDIM_(n), const SimTK_C_ *ap,  float *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo));
double zlansp_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FDIM_(n), const SimTK_Z_ *ap,  double *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo));

double slanhs_( SimTK_FOPT_(norm), SimTK_FDIM_(n), const float *a,  SimTK_FDIM_(lda), float *work, SimTK_FLEN_(norm));
double dlanhs_( SimTK_FOPT_(norm), SimTK_FDIM_(n), const double *a,  SimTK_FDIM_(lda), double *work, SimTK_FLEN_(norm));
double clanhs_( SimTK_FOPT_(norm), SimTK_FDIM_(n), const SimTK_C_ *a,  SimTK_FDIM_(lda), float *work, SimTK_FLEN_(norm));
double zlanhs_( SimTK_FOPT_(norm), SimTK_FDIM_(n), const SimTK_Z_ *a,  SimTK_FDIM_(lda), double *work, SimTK_FLEN_(norm));

double slangt_( SimTK_FOPT_(norm), SimTK_FDIM_(n), const float *dl, const float *d,  const float *du, SimTK_FLEN_(norm));
double dlangt_( SimTK_FOPT_(norm), SimTK_FDIM_(n), const double *dl, const double *d,  const double *du, SimTK_FLEN_(norm));
double clangt_( SimTK_FOPT_(norm), SimTK_FDIM_(n), const SimTK_C_ *dl, const SimTK_C_ *d,  const SimTK_C_ *du, SimTK_FLEN_(norm));
double zlangt_( SimTK_FOPT_(norm), SimTK_FDIM_(n), const SimTK_Z_ *dl, const SimTK_Z_ *d,  const SimTK_Z_ *du, SimTK_FLEN_(norm));

double clanhp_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FDIM_(n), const SimTK_C_ *ap,  float *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo));
double zlanhp_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FDIM_(n), const SimTK_Z_ *ap,  double *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo));

double clanht_( SimTK_FOPT_(norm), SimTK_FDIM_(n), const float *d, const SimTK_C_ *e,  SimTK_FLEN_(norm) );
double zlanht_( SimTK_FOPT_(norm), SimTK_FDIM_(n), const double *d, const SimTK_Z_ *e,  SimTK_FLEN_(norm) );

double slansy_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FDIM_(n), const float *a,  SimTK_FDIM_(lda), float *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo));
double dlansy_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FDIM_(n), const double *a,  SimTK_FDIM_(lda), double *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo));
double clansy_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FDIM_(n), const SimTK_C_ *a,  SimTK_FDIM_(lda), float *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo));
double zlansy_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FDIM_(n), const SimTK_Z_ *a,  SimTK_FDIM_(lda), double *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo));

double clanhe_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FDIM_(n), const SimTK_C_ *a,  SimTK_FDIM_(lda), float *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo));
double zlanhe_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FDIM_(n), const SimTK_Z_ *a,  SimTK_FDIM_(lda), double *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo));

double slantb_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(k), const float *ab,  SimTK_FDIM_(ldab), float *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

double dlantb_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(k), const double *ab,  SimTK_FDIM_(ldab), double *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

double clantb_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_C_ *ab,  SimTK_FDIM_(ldab), float *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

double zlantb_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_Z_ *ab,  SimTK_FDIM_(ldab), double *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

double slantp_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), const float *ap, float *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

double dlantp_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), const double *ap, double *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

double clantp_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), const SimTK_C_ *ap, float *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

double zlantp_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), const SimTK_Z_ *ap, double *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

double slantr_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(m), SimTK_FDIM_(n), const float *a,  SimTK_FDIM_(lda), float *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

double dlantr_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(m), SimTK_FDIM_(n), const double *a,  SimTK_FDIM_(lda), double *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

double clantr_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(m), SimTK_FDIM_(n), const SimTK_C_ *a,  SimTK_FDIM_(lda), float *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

double zlantr_( SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(m), SimTK_FDIM_(n), const SimTK_Z_ *a,  SimTK_FDIM_(lda), double *work, SimTK_FLEN_(norm), SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

double slapy2_( SimTK_S_INPUT_(x), SimTK_S_INPUT_(y) );
double dlapy2_( SimTK_D_INPUT_(x), SimTK_D_INPUT_(y) );
double slapy3_( SimTK_S_INPUT_(x), SimTK_S_INPUT_(y), SimTK_S_INPUT_(z) );
double dlapy3_( SimTK_D_INPUT_(x), SimTK_D_INPUT_(y), SimTK_D_INPUT_(z) );

void chetd2_( SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), float *d, float *e, SimTK_C_ *tau, SimTK_INFO_, SimTK_FLEN_(uplo));

void zhetd2_( SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), double *d, double *e, SimTK_Z_ *tau, SimTK_INFO_, SimTK_FLEN_(uplo));

void slamc1_(int *beta, int *t, int *rnd, int *ieee1);

void slamc2_(int *beta, int *t, int *rnd, float *eps, int *emin, float *rmin, int *emax, float *rmax);

void slamc4_(int *emin, float *start, int *base);

void slamc5_(int *beta, int *p, int *emin, int *ieee, int *emax, float *rmax);

void slamrg_( SimTK_I_INPUT_(n1), SimTK_I_INPUT_(n2), const float *a, SimTK_I_INPUT_(strd1), SimTK_I_INPUT_(strd2), int *index);

void slanv2_(float *a, float *b, float *c__, float *d__, float *rt1r, float *rt1i, float *rt2r, float *rt2i, float *cs, float *sn);

void slapll_(SimTK_FDIM_(n), float *x, SimTK_FINC_(x), float *y, SimTK_FINC_(y), SimTK_S_OUTPUT_(ssmin));

void slapmt_(SimTK_I_INPUT_(forwrd), SimTK_FDIM_(m), SimTK_FDIM_(n), float *x, SimTK_FDIM_(ldx), SimTK_FDIM_(k));

void slaqgb_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), float *ab, SimTK_FDIM_(ldab), float *r__, float *c__, SimTK_S_OUTPUT_(rowcnd), SimTK_S_OUTPUT_(colcnd), SimTK_S_INPUT_(amax), SimTK_CHAR_OUTPUT_(equed), SimTK_FLEN_(equed));

void slaqge_(SimTK_FDIM_(m), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *r__, float *c__, SimTK_S_INPUT_(rowcnd), SimTK_S_INPUT_(colcnd), SimTK_S_INPUT_(amax),  SimTK_CHAR_OUTPUT_(equed), SimTK_FLEN_(equed));

void slaqp2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_I_INPUT_(offset), float *a, SimTK_FDIM_(lda), int *jpvt, float *tau, float *vn1, float *vn2, float *work);

void slaqps_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_I_INPUT_(offset), SimTK_I_INPUT_(nb), SimTK_I_OUTPUT_(kb), float *a, SimTK_FDIM_(lda), int *jpvt, float *tau, float *vn1, float *vn2, float *auxv, float *f, SimTK_FDIM_(ldf));

void slaqr0_( SimTK_I_INPUT_(wantt), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), float *h, SimTK_FDIM_(ldh), float *wr, float *wi, SimTK_I_INPUT_(iloz), SimTK_I_INPUT_(ihiz), float *z, SimTK_FDIM_(ldz), float *work, SimTK_FDIM_(lwork), SimTK_INFO_ );

void slaqr1_(SimTK_FDIM_(n), const float *h, SimTK_FDIM_(ldh), SimTK_S_INPUT_(sr1), SimTK_S_INPUT_(si1), SimTK_S_INPUT_(sr2), SimTK_S_INPUT_(si2), float *v );

void slaqr2_( SimTK_I_INPUT_(wantt), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), SimTK_I_INPUT_(ktop), SimTK_I_INPUT_(kbot),  SimTK_I_INPUT_(nw), float *h, SimTK_FDIM_(ldh), SimTK_I_INPUT_(iloz), SimTK_I_INPUT_(ihiz),  float *z, SimTK_FDIM_(ldz), SimTK_I_OUTPUT_(ns), SimTK_I_OUTPUT_(nd), float *sr, float *si, float *v, SimTK_FDIM_(ldv), SimTK_I_INPUT_(nh), float *t, SimTK_FDIM_(ldt), SimTK_I_INPUT_(nv), float *wv, SimTK_I_INPUT_(ldwv), float *work, SimTK_FDIM_(lwork) );

void slaqr3_( SimTK_I_INPUT_(wantt), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), SimTK_I_INPUT_(ktop), SimTK_I_INPUT_(kbot),  SimTK_I_INPUT_(nw), float *h, SimTK_FDIM_(ldh), SimTK_I_INPUT_(iloz), SimTK_I_INPUT_(ihiz),  float *z, SimTK_FDIM_(ldz), SimTK_I_OUTPUT_(ns), SimTK_I_OUTPUT_(nd), float *sr, float *si, float *v, SimTK_FDIM_(ldv), SimTK_I_INPUT_(nh), float *t, SimTK_FDIM_(ldt), SimTK_I_INPUT_(nv), float *wv, SimTK_I_INPUT_(ldwv), float *work, SimTK_FDIM_(lwork) );

void slaqr4_( SimTK_I_INPUT_(wantt), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), float *h, SimTK_FDIM_(ldh), float *wr, float *wi, SimTK_I_INPUT_(iloz), SimTK_I_INPUT_(ihiz), float *z, SimTK_FDIM_(ldz), float *work, SimTK_FDIM_(lwork), SimTK_INFO_ );

void slaqr5_( SimTK_I_INPUT_(wantt), SimTK_I_INPUT_(wantz), SimTK_I_INPUT_(kacc22), SimTK_FDIM_(n), SimTK_I_INPUT_(ktop), SimTK_I_INPUT_(kbot),  SimTK_I_INPUT_(nshifts), const float *sr, const float *si, float *h, SimTK_FDIM_(ldh), SimTK_I_INPUT_(iloz), SimTK_I_INPUT_(ihiz),  float *z, SimTK_FDIM_(ldz), float *v, SimTK_FDIM_(ldv), float *u, SimTK_FDIM_(ldu), SimTK_I_INPUT_(nh), float *wh, SimTK_FDIM_(ldwh), SimTK_I_INPUT_(nv), float *wv, SimTK_FDIM_(ldwv)  );

void slaqsb_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), float *ab, SimTK_FDIM_(ldab), float *s, SimTK_S_INPUT_(scond), SimTK_S_INPUT_(amax), SimTK_CHAR_OUTPUT_(equed), SimTK_FLEN_(uplo), SimTK_FLEN_(equed));

void slaqsp_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *ap, const float *s, SimTK_S_INPUT_(scond), SimTK_S_INPUT_(amax), SimTK_CHAR_OUTPUT_(equed), SimTK_FLEN_(uplo), SimTK_FLEN_(equed));

void slaqsy_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), const float *s, SimTK_S_INPUT_(scond), SimTK_S_INPUT_(amax), SimTK_CHAR_OUTPUT_(equed), SimTK_FLEN_(uplo), SimTK_FLEN_(equed));

void slaqtr_(SimTK_I_INPUT_(ltran), SimTK_I_INPUT_(lfloat), SimTK_FDIM_(n), const float *t, SimTK_FDIM_(ldt), const float *b, const float *w, SimTK_S_OUTPUT_(scale), float *x, float *work, SimTK_INFO_);

void slar1v_(SimTK_FDIM_(n), SimTK_I_INPUT_(b1), SimTK_I_INPUT_(bn), SimTK_S_INPUT_(lambda), const float *d__, const float *l, const float *ld, const float *lld, const float *pivmin, const float *gaptol, float *z__, SimTK_I_INPUT_(wantc), SimTK_I_OUTPUT_(negcnt), SimTK_S_OUTPUT_(ztz), SimTK_S_OUTPUT_(mingma), SimTK_I_OUTPUT_(r__), int *isuppz, SimTK_S_OUTPUT_(nrminv), SimTK_S_OUTPUT_(resid), SimTK_S_OUTPUT_(rqcorr), float *work);

void slar2v_(SimTK_FDIM_(n), float *x, float *y, float *z__, int *incx, const float *c__, const float *s, SimTK_FINC_(c));

void slarf_(SimTK_FOPT_(side), SimTK_FDIM_(m), SimTK_FDIM_(n), float *v, SimTK_FINC_(v), float *tau, float *c__, SimTK_FDIM_(ldc), float *work, SimTK_FLEN_(side));

void slarfb_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FOPT_(direct), SimTK_FOPT_(storev), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), const float *v, SimTK_FDIM_(ldv), const float *t, SimTK_FDIM_(ldt), float *c__, SimTK_FDIM_(ldc), float *work, SimTK_FDIM_(ldwork), SimTK_FLEN_(side), SimTK_FLEN_(trans), SimTK_FLEN_(direct), SimTK_FLEN_(storev) );

void slarfg_(SimTK_FDIM_(n), SimTK_S_OUTPUT_(alpha), float *x, SimTK_FINC_(x), SimTK_S_OUTPUT_(tau) );

void slarft_(SimTK_FOPT_(direct), SimTK_FOPT_(storev), SimTK_FDIM_(n), SimTK_I_INPUT_(k), float *v, SimTK_FDIM_(ldv), const float *tau, float *t, SimTK_FDIM_(ldt), SimTK_FLEN_(direct), SimTK_FLEN_(storev));

void slarfx_(SimTK_FOPT_(side), SimTK_FDIM_(m), SimTK_FDIM_(n), const float *v, const float *tau, float *c__, SimTK_FDIM_(ldc), float *work, SimTK_FLEN_(side));

void slargv_(SimTK_FDIM_(n), float *x, SimTK_FINC_(x), float *y, SimTK_FINC_(y), float *c__, SimTK_FINC_(c));

void slarnv_(SimTK_I_INPUT_(idist), int *iseed, SimTK_FDIM_(n), float *x);
 
void slarra_(SimTK_FDIM_(n), const float *d, float *e, float *e2, SimTK_S_INPUT_(spltol), SimTK_S_INPUT_(tnrm), SimTK_I_OUTPUT_(nsplit), int *isplit, SimTK_INFO_);

void slarrb_(SimTK_FDIM_(n), const float *d__, const float *lld, SimTK_I_INPUT_(ifirst), SimTK_I_INPUT_(ilast), SimTK_S_INPUT_(rtol1), SimTK_S_INPUT_(rtol2), SimTK_I_INPUT_(offset), float *w, float *wgap, float *werr, float *work, int *iwork, SimTK_S_INPUT_(pivmin), SimTK_S_INPUT_(spdiam),  SimTK_I_INPUT_(twist), SimTK_INFO_);

void slarrc_(SimTK_FOPT_(jobt), SimTK_FDIM_(n), SimTK_S_INPUT_(vl), SimTK_S_INPUT_(vu), const float *d, const float *e,  SimTK_S_INPUT_(pivmin), SimTK_I_OUTPUT_(eigcnt), SimTK_I_OUTPUT_(lcnt), SimTK_I_OUTPUT_(rcnt), SimTK_INFO_, SimTK_FLEN_(jobt) );

void slarrd_( SimTK_FOPT_(range), SimTK_FOPT_(order), SimTK_FDIM_(n), SimTK_S_INPUT_(vl), SimTK_S_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), const float *gers, SimTK_S_INPUT_(reltol), const float *d, const float *e, const float *e2, SimTK_S_INPUT_(pivmin), SimTK_I_INPUT_(nsplit), const int *isplit, SimTK_I_OUTPUT_(m), float *w, float *werr, SimTK_S_OUTPUT_(wl), SimTK_S_OUTPUT_(wu), int *iblock, int *indexw, float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(range), SimTK_FLEN_(order) );

void slarre_(SimTK_FOPT_(range), SimTK_FDIM_(n), SimTK_S_OUTPUT_(vl), SimTK_S_OUTPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), float *d__, float *e, float *e2, SimTK_S_INPUT_(rtol1), SimTK_S_INPUT_(rtol2), SimTK_S_INPUT_(spltol), SimTK_I_OUTPUT_(nsplit),  int *isplit, SimTK_I_OUTPUT_(m), float *w, float *werr, float *wgap, int *iblock, int *indexw, float *gers, float *pivmin, float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(range) );

void slarrf_(SimTK_FDIM_(n), float *d__, float *l, float *ld, SimTK_I_INPUT_(clstrt), SimTK_I_INPUT_(clend),  const float *w, float *wgap, float *werr, SimTK_S_INPUT_(spdiam), SimTK_S_INPUT_(clgapl), SimTK_S_INPUT_(clgapr), SimTK_S_INPUT_(pivmin), SimTK_S_OUTPUT_(sigma), float *dplus, float *lplus, float *work, SimTK_INFO_);

void slarrj_(SimTK_FDIM_(n), const float *d, const float *e2, SimTK_I_INPUT_(ifirst), SimTK_I_INPUT_(ilast), SimTK_S_INPUT_(rtol), SimTK_I_INPUT_(offset), float *w, float *werr, float *work, int *iwork, SimTK_S_INPUT_(pivmin), SimTK_S_INPUT_(spdiam), SimTK_INFO_);

void slarrk_( SimTK_FDIM_(n), SimTK_I_INPUT_(iw), SimTK_S_INPUT_(gl), SimTK_S_INPUT_(gu), const float *d, const float *e2, SimTK_S_INPUT_(pivmin),  SimTK_S_INPUT_(reltol), SimTK_S_OUTPUT_(w), SimTK_S_OUTPUT_(werr), SimTK_INFO_);

void slarrr_( SimTK_FDIM_(n), const float *d__, float *e, SimTK_INFO_ );

void slarrv_(SimTK_FDIM_(n), SimTK_S_INPUT_(vl), SimTK_S_INPUT_(vu), float *d__, float *l, SimTK_S_INPUT_(pivmin), const int *isplit, SimTK_FDIM_(m), SimTK_I_INPUT_(dol), SimTK_I_INPUT_(dou), SimTK_S_INPUT_(minrgp), SimTK_S_INPUT_(rtol1), SimTK_S_INPUT_(rtol2), float *w, float *werr, float *wgap, const int *iblock, const int *indexw, const float *gers,  float *z__, SimTK_FDIM_(ldz), int *isuppz, float *work, int *iwork, SimTK_INFO_);

void slartg_(const float *f, const float *g, float *cs, float *sn, float *r__);

void slartv_(SimTK_FDIM_(n), float *x, SimTK_FINC_(x), float *y, SimTK_FINC_(y), const float *c__, const float *s, SimTK_FINC_(c));

void slaruv_(int *iseed, SimTK_FDIM_(n), float *x);

void slarz_(SimTK_FOPT_(side), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_I_INPUT_(L), const float *v, SimTK_FINC_(v), SimTK_S_INPUT_(tau), float *c__, SimTK_FDIM_(ldc), float *work, SimTK_FLEN_(side));

void slarzb_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FOPT_(direct), SimTK_FOPT_(storev), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_FDIM_(l), const float *v, SimTK_FDIM_(ldv), const float *t, SimTK_FDIM_(ldt), float *c__, SimTK_FDIM_(ldc), float *work, SimTK_FDIM_(ldwork), SimTK_FLEN_(side), SimTK_FLEN_(trans), SimTK_FLEN_(direct), SimTK_FLEN_(storev));

void slarzt_(SimTK_FOPT_(direct), SimTK_FOPT_(storev), SimTK_FDIM_(n),  SimTK_FDIM_(k), float *v, SimTK_FDIM_(ldv), const float *tau, float *t, SimTK_FDIM_(ldt), SimTK_FLEN_(direct), SimTK_FLEN_(storev));

void slas2_(float *f, float *g, float *h__, SimTK_S_OUTPUT_(ssmin), float *ssmax);

void slascl_(SimTK_FOPT_(type__), SimTK_FDIM_(kl), SimTK_FDIM_(ku), const float *cfrom, const float *cto, SimTK_FDIM_(m), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(type));

void slasd0_(SimTK_FDIM_(n), SimTK_I_INPUT_(sqre), float *d__, float *e, float *u, SimTK_FDIM_(ldu), float *vt, SimTK_FDIM_(ldvt), SimTK_I_INPUT_(smlsiz), int *iwork, float *work, SimTK_INFO_);

void slasd1_(SimTK_I_INPUT_(nl), SimTK_I_INPUT_(nr), SimTK_I_INPUT_(sqre), float *d__, SimTK_S_INPUT_(alpha), SimTK_S_INPUT_(beta), float *u, SimTK_FDIM_(ldu), float *vt, SimTK_FDIM_(ldvt), int *idxq, int *iwork, float *work, SimTK_INFO_);

void slasd2_(SimTK_I_INPUT_(nl), SimTK_I_INPUT_(nr), SimTK_I_INPUT_(sqre), SimTK_I_OUTPUT_(k), float *d__, float *z__, SimTK_S_INPUT_(alpha), SimTK_S_INPUT_(beta), float *u, SimTK_FDIM_(ldu), float *vt, SimTK_FDIM_(ldvt), float *dsigma, float *u2, SimTK_FDIM_(ldu2), float *vt2, SimTK_FDIM_(ldvt2), int *idxp, int *idx, int *idxc, int *idxq, int *coltyp, SimTK_INFO_);

void slasd3_(SimTK_I_INPUT_(nl), SimTK_I_INPUT_(nr), SimTK_I_INPUT_(sqre), SimTK_I_INPUT_(k), float *d__, float *q, SimTK_FDIM_(ldq), const float *dsigma, const float *u, SimTK_FDIM_(ldu), const float *u2, SimTK_FDIM_(ldu2), const float *vt, SimTK_FDIM_(ldvt), const float *vt2, SimTK_FDIM_(ldvt2), const int *idxc, const int *ctot, const float *z__, SimTK_INFO_);


void slasd4_(SimTK_FDIM_(n), SimTK_I_INPUT_(i__), const float *d__, const float *z__, float *delta, SimTK_S_INPUT_(rho), SimTK_S_OUTPUT_(sigma), float *work, SimTK_INFO_);

void slasd5_( SimTK_I_INPUT_(i__), const float *d__, const float *z__, float *delta, SimTK_S_INPUT_(rho), SimTK_S_OUTPUT_(dsigma), float *work);

void slasd6_(SimTK_I_INPUT_(icompq), SimTK_I_INPUT_(nl), SimTK_I_INPUT_(nr), SimTK_I_INPUT_(sqre), float *d__, float *vf, float *vl, SimTK_S_INPUT_(alpha), SimTK_S_INPUT_(beta), int *idxq, int *perm, SimTK_I_OUTPUT_(givptr), SimTK_I_OUTPUT_(givcol), SimTK_FDIM_(ldgcol), float *givnum, SimTK_FDIM_(ldgnum), float *poles, float *difl, float *difr, float *z__, SimTK_I_OUTPUT_(k), SimTK_S_OUTPUT_(c__), SimTK_S_OUTPUT_(s), float *work, int *iwork, SimTK_INFO_);

void slasd7_( SimTK_I_INPUT_(icompq), SimTK_I_INPUT_(nl), SimTK_I_INPUT_(nr), SimTK_I_INPUT_(sqre), SimTK_I_OUTPUT_(k), float *d__, float *z__, float *zw, float *vf, float *vfw, float *vl, float *vlw, SimTK_S_INPUT_(alpha), SimTK_S_INPUT_(beta), float *dsigma, int *idx, int *idxp, const int *idxq, int *perm, SimTK_I_OUTPUT_(givptr), int *givcol, SimTK_FDIM_(LDGCOL), float *givnum, SimTK_FDIM_(ldgnum), SimTK_S_OUTPUT_(c__), SimTK_S_OUTPUT_(s), SimTK_INFO_);

void slasd8_(SimTK_I_INPUT_(icompq), SimTK_FDIM_(k), float *d__, const float *z__, float *vf, float *vl, float *difl, float *difr, SimTK_FDIM_(lddifr), const float *dsigma, float *work, SimTK_INFO_);

void slasd9_(SimTK_I_INPUT_(icompq), SimTK_FDIM_(ldu), SimTK_FDIM_(k), float *d__, float *z__, float *vf, float *vl, float *difl, float *difr, float *dsigma, float *work, SimTK_INFO_);

void slasda_(SimTK_I_INPUT_(icompq), SimTK_I_INPUT_(smlsiz), SimTK_FDIM_(n), SimTK_I_INPUT_(sqre), float *d__, float *e, float *u, SimTK_FDIM_(ldu), float *vt, SimTK_FDIM_(k), float *difl, float *difr, float *z__, float *poles, int *givptr, int *givcol, SimTK_FDIM_(ldgcol), int *perm, float *givnum, float *c__, float *s, float *work, int *iwork, SimTK_INFO_);

void slasdq_(SimTK_FOPT_(uplo), SimTK_I_INPUT_(sqre), SimTK_FDIM_(n), SimTK_I_INPUT_(ncvt), SimTK_I_INPUT_(nru), SimTK_I_INPUT_(ncc), float *d__, float *e, float *vt, SimTK_FDIM_(ldvt), float *u, SimTK_FDIM_(ldu), float *c__, SimTK_FDIM_(ldc), float *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void slasdt_( SimTK_FDIM_(n), SimTK_I_OUTPUT_(lvl), SimTK_I_OUTPUT_(nd), int *inode, int *ndiml, int *ndimr, SimTK_I_INPUT_(msub) );

void slaset_(SimTK_FOPT_(uplo), SimTK_FDIM_(m), SimTK_FDIM_(n),  SimTK_S_INPUT_(alpha), SimTK_S_INPUT_(beta), float *a, SimTK_FDIM_(lda), SimTK_FLEN_(uplo));

void slasq1_(SimTK_FDIM_(n), float *d__, float *e, float *work, SimTK_INFO_);

void slasq2_(SimTK_FDIM_(n), float *z__, SimTK_INFO_);

void slasq3_(SimTK_I_INPUT_(i0), SimTK_I_INPUT_(n0), const float *z__, SimTK_I_INPUT_(pp), SimTK_S_OUTPUT_(dmin__), SimTK_S_OUTPUT_(sigma), SimTK_S_OUTPUT_(desig), SimTK_S_INPUT_(qmax), SimTK_I_OUTPUT_(nfail), SimTK_I_OUTPUT_(iter), SimTK_I_OUTPUT_(ndiv), SimTK_I_INPUT_(ieee) );

void slasq4_( SimTK_I_INPUT_(i0), SimTK_I_INPUT_(n0), const float *z__, SimTK_I_INPUT_(pp), SimTK_I_INPUT_(noin), SimTK_S_INPUT_(dmin),  SimTK_S_INPUT_(dmin1), SimTK_S_INPUT_(dmin2), SimTK_S_INPUT_(dn), SimTK_S_INPUT_(dn1), SimTK_S_INPUT_(dn2), SimTK_S_OUTPUT_(tau), SimTK_I_OUTPUT_(ttype));

void slasq5_(SimTK_I_INPUT_(i0), SimTK_I_INPUT_(n0), const float *z__, SimTK_I_INPUT_(pp), SimTK_S_INPUT_(tau), SimTK_S_OUTPUT_(dmin), SimTK_S_OUTPUT_(dmin1), SimTK_S_OUTPUT_(dmin2), SimTK_S_OUTPUT_(dn), SimTK_S_OUTPUT_(dnm1), SimTK_S_OUTPUT_(dnm2), SimTK_I_INPUT_(ieee) );

void slasq6_(SimTK_I_INPUT_(i0), SimTK_I_INPUT_(n0), const float *z__, SimTK_I_INPUT_(pp), SimTK_S_OUTPUT_(dmin), SimTK_S_OUTPUT_(dmin1), SimTK_S_OUTPUT_(dmin2), SimTK_S_OUTPUT_(dn), SimTK_S_OUTPUT_(dnm1), SimTK_S_OUTPUT_(dnm2));

void slasr_(SimTK_FOPT_(side), SimTK_FOPT_(pivot), SimTK_FOPT_(direct), SimTK_FDIM_(m), SimTK_FDIM_(n), const float *c__, const float *s, float *a, SimTK_FDIM_(lda), SimTK_FLEN_(side), SimTK_FLEN_(pivot), SimTK_FLEN_(direct));

void slasrt_(SimTK_FOPT_(id), SimTK_FDIM_(n), float *d__, SimTK_INFO_, SimTK_FLEN_(id));

void slassq_(SimTK_FDIM_(n), const float *x, SimTK_FINC_(x), SimTK_S_OUTPUT_(scale), float *sumsq);

void slasv2_(SimTK_S_INPUT_(f), SimTK_S_INPUT_(g), SimTK_S_INPUT_(h__), SimTK_S_OUTPUT_(ssmin), SimTK_S_OUTPUT_(ssmax) , SimTK_S_OUTPUT_(snr) , SimTK_S_OUTPUT_(csr) , SimTK_S_OUTPUT_(snl) , SimTK_S_OUTPUT_(csl) );

void slaswp_(SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), SimTK_I_OUTPUT_(k1), SimTK_I_OUTPUT_(k2), const int *ipiv, SimTK_FINC_(x));

void slasy2_(SimTK_I_INPUT_(ltranl), SimTK_I_INPUT_(ltranr), SimTK_I_INPUT_(isgn),  SimTK_I_INPUT_(n1), SimTK_I_INPUT_(n2), const float *tl, SimTK_FDIM_(ldtl), const float *tr, SimTK_FDIM_(ldtr), const float *b, SimTK_FDIM_(ldb), SimTK_S_OUTPUT_(scale), float *x, SimTK_FDIM_(ldx), float *xnorm, SimTK_INFO_);


void slasyf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_I_INPUT_(nb), SimTK_I_INPUT_(kb), float *a, SimTK_FDIM_(lda), int *ipiv, float *w, SimTK_FDIM_(ldw), int *info, SimTK_FLEN_(uplo));

void slatbs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FOPT_(normin), SimTK_FDIM_(n), SimTK_I_INPUT_(kd), const float *ab, SimTK_FDIM_(ldab), float *x, SimTK_S_OUTPUT_(scale), float *cnorm, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag), SimTK_FLEN_(normin) );

void slatdf_(SimTK_I_INPUT_(ijob), SimTK_FDIM_(n), const float *z__, SimTK_FDIM_(ldz), float *rhs, SimTK_S_OUTPUT_(rdsum), SimTK_S_OUTPUT_(rdscal), const int *ipiv, const int *jpiv);

void slatps_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FOPT_(normin), SimTK_FDIM_(n), const float *ap, float *x, SimTK_S_OUTPUT_(scale), float *cnorm, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag), SimTK_FLEN_(normin) );

void slatrd_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nb), float *a, SimTK_FDIM_(lda), float *e, float *tau, float *w, SimTK_FDIM_(ldw), SimTK_FLEN_(uplo));

void slatrs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FOPT_(normin), SimTK_FDIM_(n), const float *a, SimTK_FDIM_(lda), float *x, SimTK_S_OUTPUT_(scale), float *cnorm, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag), SimTK_FLEN_(normin));

void slatrz_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(l), float *a, SimTK_FDIM_(lda), float *tau, float *work);

void slatzm_(SimTK_FOPT_(side), SimTK_FDIM_(m), SimTK_FDIM_(n), const float *v, SimTK_FINC_(v), SimTK_S_INPUT_(tau), float *c1, float *c2, SimTK_FDIM_(ldc), float *work, SimTK_FLEN_(side));

void slauu2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo));

void slauum_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo));

void slazq3( SimTK_I_INPUT_(io), SimTK_I_INPUT_(no), const float *z, SimTK_I_INPUT_(pp), SimTK_S_OUTPUT_(dmin), SimTK_S_OUTPUT_(sigma), SimTK_S_OUTPUT_(desig), SimTK_S_INPUT_(qmax), SimTK_I_OUTPUT_(nfail), SimTK_I_OUTPUT_(iter), SimTK_I_OUTPUT_(ndiv), SimTK_I_INPUT_(ieee),  SimTK_I_INPUT_(ttype), SimTK_S_OUTPUT_(dmin1), SimTK_S_OUTPUT_(dmin2), SimTK_S_OUTPUT_(dn), SimTK_S_OUTPUT_(dn1), SimTK_S_OUTPUT_(dn2), SimTK_S_OUTPUT_(tau) );

void slazq4( SimTK_I_INPUT_(io), SimTK_I_INPUT_(no), const float *z, SimTK_I_INPUT_(pp), SimTK_I_INPUT_(noin),  SimTK_S_INPUT_(dmin), SimTK_S_INPUT_(dmin1), SimTK_S_INPUT_(dmin2), SimTK_S_INPUT_(dn), SimTK_S_INPUT_(dn1), SimTK_S_INPUT_(dn2), SimTK_S_OUTPUT_(tau), SimTK_I_OUTPUT_(ttype), SimTK_S_OUTPUT_(g) );

void sopgtr_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *ap, float *tau, float *q, SimTK_FDIM_(ldq), float *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void sopmtr_(SimTK_FOPT_(side), SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), float *ap, float *tau, float *c__, SimTK_FDIM_(ldc), float *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(uplo), SimTK_FLEN_(trans));

void sorg2l_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), float *a, SimTK_FDIM_(lda), float *tau, float *work, SimTK_INFO_);

void sorg2r_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), float *a, SimTK_FDIM_(lda), float *tau, float *work, SimTK_INFO_);

void sorgbr_(SimTK_FOPT_(vect), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), float *a, SimTK_FDIM_(lda), float *tau, float *work, SimTK_FDIM_(lwork), int *info, SimTK_FLEN_(vect));

void sorghr_(SimTK_FDIM_(n), SimTK_FDIM_(ilo), SimTK_FDIM_(ihi), float *a, SimTK_FDIM_(lda), float *tau, float *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void sorgl2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), float *a, SimTK_FDIM_(lda), float *tau, float *work, SimTK_INFO_);

void sorglq_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), float *a, SimTK_FDIM_(lda), float *tau, float *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void sorgql_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), float *a, SimTK_FDIM_(lda), float *tau, float *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void sorgqr_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), float *a, SimTK_FDIM_(lda), float *tau, float *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void sorgr2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), float *a, SimTK_FDIM_(lda), float *tau, float *work, SimTK_INFO_);

void sorgrq_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), float *a, SimTK_FDIM_(lda), float *tau, float *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void sorgtr_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *tau, float *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(uplo));

void sorm2l_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), float *a, SimTK_FDIM_(lda), float *tau, float *c__, SimTK_FDIM_(ldc), float *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void sorm2r_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), float *a, SimTK_FDIM_(lda), float *tau, float *c__, SimTK_FDIM_(ldc), float *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void sormbr_(SimTK_FOPT_(vect), SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), float *a, SimTK_FDIM_(lda), float *tau, float *c__, SimTK_FDIM_(ldc), float *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(vect), SimTK_FLEN_(side), SimTK_FLEN_(trans));

void sormhr_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(ilo), SimTK_FDIM_(ihi), float *a, SimTK_FDIM_(lda), float *tau, float *c__, SimTK_FDIM_(ldc), float *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void sorml2_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), float *a, SimTK_FDIM_(lda), float *tau, float *c__, SimTK_FDIM_(ldc), float *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void sormlq_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), float *a, SimTK_FDIM_(lda), float *tau, float *c__, SimTK_FDIM_(ldc), float *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void sormql_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), float *a, SimTK_FDIM_(lda), float *tau, float *c__, SimTK_FDIM_(ldc), float *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void sormqr_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), float *a, SimTK_FDIM_(lda), float *tau, float *c__, SimTK_FDIM_(ldc), float *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void sormr2_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), float *a, SimTK_FDIM_(lda), float *tau, float *c__, SimTK_FDIM_(ldc), float *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void sormr3_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), int *l, float *a, SimTK_FDIM_(lda), float *tau, float *c__, SimTK_FDIM_(ldc), float *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void sormrq_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), float *a, SimTK_FDIM_(lda), float *tau, float *c__, SimTK_FDIM_(ldc), float *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void sormrz_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_FDIM_(l), float *a, SimTK_FDIM_(lda), float *tau, float *c__, SimTK_FDIM_(ldc), float *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void sormtr_(SimTK_FOPT_(side), SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *tau, float *c__, SimTK_FDIM_(ldc), float *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(uplo), SimTK_FLEN_(trans));

void spbcon_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), const float *ab, SimTK_FDIM_(ldab), SimTK_S_INPUT_(anorm), SimTK_S_OUTPUT_(rcond), float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void spbequ_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_I_INPUT_(kd), const float *ab, SimTK_FDIM_(ldab), SimTK_S_OUTPUT_(s), SimTK_S_OUTPUT_(scond), SimTK_S_OUTPUT_(amax), SimTK_INFO_, SimTK_FLEN_(uplo));

void spbrfs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_FDIM_(nrhs), const float *ab, SimTK_FDIM_(ldab), const float *afb, SimTK_FDIM_(ldafb), const float *b, SimTK_FDIM_(ldb), float *x, SimTK_FDIM_(ldx), float *ferr, float *berr, float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void spbstf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), float *ab, SimTK_FDIM_(ldab), SimTK_INFO_, SimTK_FLEN_(uplo));

void spbsv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_FDIM_(nrhs), float *ab, SimTK_FDIM_(ldab), float *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void spbsvx_(SimTK_FOPT_(fact), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_FDIM_(nrhs), float *ab, SimTK_FDIM_(ldab), float *afb, SimTK_FDIM_(ldafb), SimTK_CHAR_OUTPUT_(equed), float *s, float *b, SimTK_FDIM_(ldb), float *x, SimTK_FDIM_(ldx), SimTK_S_OUTPUT_(rcond), float *ferr, float *berr, float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(uplo), SimTK_FLEN_(equed));

void spbtf2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), float *ab, SimTK_FDIM_(ldab), SimTK_INFO_, SimTK_FLEN_(uplo) );

void spbtrf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), float *ab, SimTK_FDIM_(ldab), SimTK_INFO_, SimTK_FLEN_(uplo));

void spbtrs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_FDIM_(nrhs), const float *ab, SimTK_FDIM_(ldab), float *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void spocon_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const float *a, SimTK_FDIM_(lda), SimTK_S_INPUT_(anorm), SimTK_S_OUTPUT_(rcond), float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void spoequ_(SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *s, SimTK_S_OUTPUT_(scond), SimTK_S_OUTPUT_(amax), SimTK_INFO_);

void sporfs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const float *a, SimTK_FDIM_(lda), const float *af, SimTK_FDIM_(ldaf), const float *b, SimTK_FDIM_(ldb), float *x, SimTK_FDIM_(ldx), float *ferr, float *berr, float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void sposv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void sposvx_(SimTK_FOPT_(fact), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), float *a, SimTK_FDIM_(lda), float *af, SimTK_FDIM_(ldaf), SimTK_CHAR_OUTPUT_(equed), float *s, float *b, SimTK_FDIM_(ldb), float *x, SimTK_FDIM_(ldx), SimTK_S_OUTPUT_(rcond), float *ferr, float *berr, float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(uplo), SimTK_FLEN_(equed));

void spotf2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo));

void spotrf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo));

void spotri_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo));

void spotrs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void sppcon_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const float *ap, SimTK_S_INPUT_(anorm), SimTK_S_OUTPUT_(rcond), float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void sppequ_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const float *ap, float *s, SimTK_S_OUTPUT_(scond), SimTK_S_OUTPUT_(amax), SimTK_INFO_, SimTK_FLEN_(uplo));

void spprfs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const float *ap, const float *afp, const float *b, SimTK_FDIM_(ldb), float *x, SimTK_FDIM_(ldx), float *ferr, float *berr, float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void sppsv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), float *ap, float *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void sppsvx_(SimTK_FOPT_(fact), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), float *ap, float *afp, SimTK_CHAR_OUTPUT_(equed), float *s, float *b, SimTK_FDIM_(ldb), float *x, SimTK_FDIM_(ldx), SimTK_S_OUTPUT_(rcond), float *ferr, float *berr, float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(uplo), SimTK_FLEN_(equed));

void spptrf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *ap, SimTK_INFO_, SimTK_FLEN_(uplo));

void spptri_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *ap, SimTK_INFO_, SimTK_FLEN_(uplo));

void spptrs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const float *ap, float *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void sptcon_(SimTK_FDIM_(n), const float *d__, const float *e, SimTK_S_INPUT_(anorm), SimTK_S_OUTPUT_(rcond), float *work, SimTK_INFO_);

void spteqr_(SimTK_FOPT_(compz), SimTK_FDIM_(n), float *d__, float *e, float *z__, SimTK_FDIM_(ldz), float *work, SimTK_INFO_, SimTK_FLEN_(compz));

void sptrfs_(SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const float *d__, const float *e, const float *df, const float *ef, const float *b, SimTK_FDIM_(ldb), float *x, SimTK_FDIM_(ldx), float *ferr, float *berr, float *work, SimTK_INFO_);

void sptsv_(SimTK_FDIM_(n), SimTK_FDIM_(nrhs), float *d__, float *e, float *b, SimTK_FDIM_(ldb), SimTK_INFO_);

void sptsvx_(SimTK_FOPT_(fact), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const float *d__, const float *e, float *df, float *ef, const float *b, SimTK_FDIM_(ldb), float *x, SimTK_FDIM_(ldx), SimTK_S_OUTPUT_(rcond), float *ferr, float *berr, float *work, SimTK_INFO_, SimTK_FLEN_(fact));

void spttrf_(SimTK_FDIM_(n), float *d__, float *e, SimTK_INFO_);

void spttrs_(SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const float *d__, const float *e, float *b, SimTK_FDIM_(ldb), SimTK_INFO_);

void sptts2_(SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const float *d__, const float *e, float *b, SimTK_FDIM_(ldb));

void srscl_(SimTK_FDIM_(n), float *sa, float *sx, SimTK_FINC_(x));

void ssbev_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), float *ab, SimTK_FDIM_(ldab), float *w, float *z__, SimTK_FDIM_(ldz), float *work, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void ssbevd_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), float *ab, SimTK_FDIM_(ldab), float *w, float *z__, SimTK_FDIM_(ldz), float *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void ssbevx_(SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), float *ab, SimTK_FDIM_(ldab), float *q, SimTK_FDIM_(ldq), SimTK_S_INPUT_(vl), SimTK_S_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_S_INPUT_(abstol), SimTK_I_OUTPUT_(m), float *w, float *z__, SimTK_FDIM_(ldz), float *work, int *iwork, int *ifail, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void ssbgst_(SimTK_FOPT_(vect), SimTK_FOPT_(uplo), SimTK_FDIM_(n), int *ka, int *kb, float *ab, SimTK_FDIM_(ldab), float *bb, SimTK_FDIM_(ldbb), float *x, SimTK_FDIM_(ldx), float *work, SimTK_INFO_, SimTK_FLEN_(vect), SimTK_FLEN_(uplo));

void ssbgv_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), int *ka, int *kb, float *ab, SimTK_FDIM_(ldab), float *bb, SimTK_FDIM_(ldbb), float *w, float *z__, SimTK_FDIM_(ldz), float *work, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void ssbgvd_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), int *ka, int *kb, float *ab, SimTK_FDIM_(ldab), float *bb, SimTK_FDIM_(ldbb), float *w, float *z__, SimTK_FDIM_(ldz), float *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void ssbgvx_(SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), int *ka, int *kb, float *ab, SimTK_FDIM_(ldab), float *bb, SimTK_FDIM_(ldbb), float *q, SimTK_FDIM_(ldq), SimTK_S_INPUT_(vl), SimTK_S_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_S_INPUT_(abstol), SimTK_I_OUTPUT_(m), float *w, float *z__, SimTK_FDIM_(ldz), float *work, int *iwork, int *ifail, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void ssbtrd_(SimTK_FOPT_(vect), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), float *ab, SimTK_FDIM_(ldab), float *d__, float *e, float *q, SimTK_FDIM_(ldq), float *work, SimTK_INFO_, SimTK_FLEN_(vect), SimTK_FLEN_(uplo));

void sspcon_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const float *ap, const int *ipiv, SimTK_S_INPUT_(anorm), SimTK_S_OUTPUT_(rcond), float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void sspev_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *ap, float *w, float *z__, SimTK_FDIM_(ldz), float *work, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void sspevd_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *ap, float *w, float *z__, SimTK_FDIM_(ldz), float *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void sspevx_(SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *ap, SimTK_S_INPUT_(vl), SimTK_S_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_S_INPUT_(abstol), SimTK_I_OUTPUT_(m), float *w, float *z__, SimTK_FDIM_(ldz), float *work, int *iwork, int *ifail, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void sspgst_(SimTK_I_INPUT_(itype), SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *ap, float *bp, SimTK_INFO_, SimTK_FLEN_(uplo));

void sspgv_(SimTK_I_INPUT_(itype), SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *ap, float *bp, float *w, float *z__, SimTK_FDIM_(ldz), float *work, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void sspgvd_(SimTK_I_INPUT_(itype), SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *ap, float *bp, float *w, float *z__, SimTK_FDIM_(ldz), float *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void sspgvx_(SimTK_I_INPUT_(itype), SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *ap, float *bp, SimTK_S_INPUT_(vl), SimTK_S_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_S_INPUT_(abstol), SimTK_I_OUTPUT_(m), float *w, float *z__, SimTK_FDIM_(ldz), float *work, int *iwork, int *ifail, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void ssprfs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const float *ap, const float *afp, const int *ipiv, const float *b, SimTK_FDIM_(ldb), float *x, SimTK_FDIM_(ldx), float *ferr, float *berr, float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void sspsv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), float *ap, int *ipiv, float *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void sspsvx_(SimTK_FOPT_(fact), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const float *ap, float *afp, int *ipiv, const float *b, SimTK_FDIM_(ldb), float *x, SimTK_FDIM_(ldx), SimTK_S_OUTPUT_(rcond), float *ferr, float *berr, float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(uplo));

void ssptrd_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *ap, float *d__, float *e, float *tau, SimTK_INFO_, SimTK_FLEN_(uplo));

void ssptrf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *ap, int *ipiv, SimTK_INFO_, SimTK_FLEN_(uplo));

void ssptri_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *ap, const int *ipiv, float *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void ssptrs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const float *ap, const int *ipiv, float *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void sstebz_(SimTK_FOPT_(range), SimTK_FOPT_(order), SimTK_FDIM_(n), SimTK_S_INPUT_(vl), SimTK_S_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_S_INPUT_(abstol), const float *d__, const float *e, SimTK_I_OUTPUT_(m), SimTK_I_OUTPUT_(nsplit), float *w, int *iblock, int *isplit, float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(range), SimTK_FLEN_(order));

void sstedc_(SimTK_FOPT_(compz), SimTK_FDIM_(n), float *d__, float *e, float *z__, SimTK_FDIM_(ldz), float *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(compz));

void sstegr_(SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FDIM_(n), float *d__, float *e, SimTK_S_INPUT_(vl), SimTK_S_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_S_INPUT_(abstol), SimTK_I_OUTPUT_(m), float *w, float *z__, SimTK_FDIM_(ldz), int *isuppz, float *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range));

void sstein_(SimTK_FDIM_(n), const float *d__, const float *e, SimTK_FDIM_(m), const float *w, const int *iblock, const int *isplit, float *z__, SimTK_FDIM_(ldz), float *work, int *iwork, int *ifail, SimTK_INFO_);

void sstemr_( SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FDIM_(n), float *d, float *e, SimTK_S_INPUT_(vl), SimTK_S_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_I_OUTPUT_(m), float *w, float *z, SimTK_FDIM_(ldz), SimTK_I_INPUT_(nzc), SimTK_I_OUTPUT_(isuppz), SimTK_I_OUTPUT_(tryrac), float *work, SimTK_FDIM_(lwork), int *iwork, SimTK_FDIM_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range) ); 

void ssteqr_(SimTK_FOPT_(compz), SimTK_FDIM_(n), float *d__, float *e, float *z__, SimTK_FDIM_(ldz), float *work, SimTK_INFO_, SimTK_FLEN_(compz));

void ssterf_(SimTK_FDIM_(n), float *d__, float *e, SimTK_INFO_);

void sstev_(SimTK_FOPT_(jobz), SimTK_FDIM_(n), float *d__, float *e, float *z__, SimTK_FDIM_(ldz), float *work, SimTK_INFO_, SimTK_FLEN_(jobz));

void sstevd_(SimTK_FOPT_(jobz), SimTK_FDIM_(n), float *d__, float *e, float *z__, SimTK_FDIM_(ldz), float *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz));

void sstevr_(SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FDIM_(n), float *d__, float *e, SimTK_S_INPUT_(vl), SimTK_S_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_S_INPUT_(abstol), SimTK_I_OUTPUT_(m), float *w, float *z__, SimTK_FDIM_(ldz), int *isuppz, float *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range));

void sstevx_(SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FDIM_(n), float *d__, float *e, SimTK_S_INPUT_(vl), SimTK_S_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_S_INPUT_(abstol), SimTK_I_OUTPUT_(m), float *w, float *z__, SimTK_FDIM_(ldz), float *work, int *iwork, int *ifail, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range));

void ssycon_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const float *a, SimTK_FDIM_(lda), const int *ipiv, SimTK_S_INPUT_(anorm), SimTK_S_OUTPUT_(rcond), float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void ssyev_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *w, float *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void ssyevd_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *w, float *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void ssyevr_(SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), SimTK_S_INPUT_(vl), SimTK_S_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_S_INPUT_(abstol), SimTK_I_OUTPUT_(m), float *w, float *z__, SimTK_FDIM_(ldz), int *isuppz, float *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void ssyevx_(SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), SimTK_S_INPUT_(vl), SimTK_S_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_S_INPUT_(abstol), SimTK_I_OUTPUT_(m), float *w, float *z__, SimTK_FDIM_(ldz), float *work, SimTK_FDIM_(lwork), int *iwork, int *ifail, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void ssygs2_(SimTK_I_INPUT_(itype), SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void ssygst_(SimTK_I_INPUT_(itype), SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void ssygv_(SimTK_I_INPUT_(itype), SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), float *w, float *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void ssygvd_(SimTK_I_INPUT_(itype), SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), float *w, float *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void ssygvx_(SimTK_I_INPUT_(itype), SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), SimTK_S_INPUT_(vl), SimTK_S_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_S_INPUT_(abstol), SimTK_I_OUTPUT_(m), float *w, float *z__, SimTK_FDIM_(ldz), float *work, SimTK_FDIM_(lwork), int *iwork, int *ifail, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void ssyrfs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const float *a, SimTK_FDIM_(lda), const float *af, SimTK_FDIM_(ldaf), const int *ipiv, const float *b, SimTK_FDIM_(ldb), float *x, SimTK_FDIM_(ldx), float *ferr, float *berr, float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void ssysv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), float *a, SimTK_FDIM_(lda), int *ipiv, float *b, SimTK_FDIM_(ldb), float *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(uplo));

void ssysvx_(SimTK_FOPT_(fact), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const float *a, SimTK_FDIM_(lda), float *af, SimTK_FDIM_(ldaf), int *ipiv, const float *b, SimTK_FDIM_(ldb), float *x, SimTK_FDIM_(ldx), SimTK_S_OUTPUT_(rcond), float *ferr, float *berr, float *work, SimTK_FDIM_(lwork), int *iwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(uplo));

void ssytd2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *d__, float *e, float *tau, SimTK_INFO_, SimTK_FLEN_(uplo));

void ssytf2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), int *ipiv, SimTK_INFO_, SimTK_FLEN_(uplo));

void ssytrd_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *d__, float *e, float *tau, float *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(uplo));

void ssytrf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), int *ipiv, float *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(uplo));

void ssytri_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const float *a, SimTK_FDIM_(lda), const int *ipiv, float *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void ssytrs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), float *a, SimTK_FDIM_(lda), int *ipiv, float *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void stbcon_(SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(kd), const float *ab, SimTK_FDIM_(ldab), SimTK_S_OUTPUT_(rcond), float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(norm), SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

void stbrfs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_FDIM_(nrhs), const float *ab, SimTK_FDIM_(ldab), const float *b, SimTK_FDIM_(ldb), const float *x, SimTK_FDIM_(ldx), float *ferr, float *berr, float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));

void stbtrs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_FDIM_(nrhs), const float *ab, SimTK_FDIM_(ldab), const float *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));

void stgevc_(SimTK_FOPT_(side), SimTK_FOPT_(howmny), const int *select, SimTK_FDIM_(n), const float *a, SimTK_FDIM_(lda), const float *b, SimTK_FDIM_(ldb), float *vl, SimTK_FDIM_(ldvl), float *vr, SimTK_FDIM_(ldvr), SimTK_I_INPUT_(mm), SimTK_FDIM_(m), float *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(howmny));

void stgex2_(SimTK_I_INPUT_(wantq), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), float *q, SimTK_FDIM_(ldq), float *z__, SimTK_FDIM_(ldz), SimTK_I_INPUT_(j1),  SimTK_I_INPUT_(n1), SimTK_I_INPUT_(n2), float *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void stgexc_(SimTK_I_INPUT_(wantq), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), float *q, SimTK_FDIM_(ldq), float *z__, SimTK_FDIM_(ldz), SimTK_I_OUTPUT_(ifst), SimTK_I_OUTPUT_(ilst), float *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void stgsen_(SimTK_I_INPUT_(ijob), SimTK_I_INPUT_(wantq), SimTK_I_INPUT_(wantz), const int *select, SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), float *alphar, float *alphai, SimTK_S_INPUT_(beta), float *q, SimTK_FDIM_(ldq), float *z__, SimTK_FDIM_(ldz), SimTK_I_OUTPUT_(m), SimTK_S_OUTPUT_(pl),  SimTK_S_OUTPUT_(pr), float *dif, float *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_);

void stgsja_(SimTK_FOPT_(jobu), SimTK_FOPT_(jobv), SimTK_FOPT_(jobq), SimTK_FDIM_(m), SimTK_FDIM_(p), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_FDIM_(l), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), SimTK_S_INPUT_(tola), SimTK_S_INPUT_(tolb), float *alpha, float *beta, float *u, SimTK_FDIM_(ldu), float *v, SimTK_FDIM_(ldv), float *q, SimTK_FDIM_(ldq), float *work, SimTK_I_OUTPUT_(ncycle), SimTK_INFO_, SimTK_FLEN_(jobu), SimTK_FLEN_(jobv), SimTK_FLEN_(jobq));

void stgsna_(SimTK_FOPT_(job), SimTK_FOPT_(howmny), const int *select, SimTK_FDIM_(n), const float *a, SimTK_FDIM_(lda), const float *b, SimTK_FDIM_(ldb), const float *vl, SimTK_FDIM_(ldvl), const float *vr, SimTK_FDIM_(ldvr), float *s, float *dif, SimTK_FDIM_(mm), SimTK_FDIM_(m), float *work, SimTK_FDIM_(lwork), int *iwork, SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(howmny));

void stgsy2_(SimTK_FOPT_(trans), SimTK_I_INPUT_(ijob), SimTK_FDIM_(m), SimTK_FDIM_(n), const float *a, SimTK_FDIM_(lda), const float *b, SimTK_FDIM_(ldb), float *c__, SimTK_FDIM_(ldc), const float *d__, SimTK_FDIM_(ldd), const float *e, SimTK_FDIM_(lde), float *f, SimTK_FDIM_(ldf), SimTK_S_OUTPUT_(scale), SimTK_S_OUTPUT_(rdsum), SimTK_S_OUTPUT_(rdscal), int *iwork, int *pq, SimTK_INFO_, SimTK_FLEN_(trans));

void stgsyl_(SimTK_FOPT_(trans), SimTK_I_INPUT_(ijob), SimTK_FDIM_(m), SimTK_FDIM_(n), const float *a, SimTK_FDIM_(lda), const float *b, SimTK_FDIM_(ldb), float *c__, SimTK_FDIM_(ldc), const float *d__, SimTK_FDIM_(ldd), const float *e, SimTK_FDIM_(lde), float *f, SimTK_FDIM_(ldf), SimTK_S_OUTPUT_(scale), float *dif, float *work, SimTK_FDIM_(lwork), int *iwork, SimTK_INFO_, SimTK_FLEN_(trans));

void stpcon_(SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), const float *ap, SimTK_S_OUTPUT_(rcond), float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(norm), SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

void stprfs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const float *ap, const float *b, SimTK_FDIM_(ldb), const float *x, SimTK_FDIM_(ldx), float *ferr, float *berr, float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));

void stptri_(SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), float *ap, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

void stptrs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const float *ap, float *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));

void strcon_(SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), const float *a, SimTK_FDIM_(lda), SimTK_S_OUTPUT_(rcond), float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(norm), SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

void strevc_(SimTK_FOPT_(side), SimTK_FOPT_(howmny), int *select, SimTK_FDIM_(n), float *t, SimTK_FDIM_(ldt), float *vl, SimTK_FDIM_(ldvl), float *vr, SimTK_FDIM_(ldvr), SimTK_FDIM_(mm), SimTK_FDIM_(m), float *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(howmny));

void strexc_(SimTK_FOPT_(compq), SimTK_FDIM_(n), float *t, SimTK_FDIM_(ldt), float *q, SimTK_FDIM_(ldq), SimTK_I_OUTPUT_(ifst), SimTK_I_OUTPUT_(ilst), float *work, SimTK_INFO_, SimTK_FLEN_(compq));

void strrfs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const float *a, SimTK_FDIM_(lda), const float *b, SimTK_FDIM_(ldb), const float *x, SimTK_FDIM_(ldx), float *ferr, float *berr, float *work, int *iwork, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));

void strsen_(SimTK_FOPT_(job), SimTK_FOPT_(compq), const int *select, SimTK_FDIM_(n), float *t, SimTK_FDIM_(ldt), float *q, SimTK_FDIM_(ldq), float *wr, float *wi, SimTK_FDIM_(m), SimTK_S_OUTPUT_(s), SimTK_S_OUTPUT_(sep), float *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(compq));

void strsna_(SimTK_FOPT_(job), SimTK_FOPT_(howmny), const int *select, SimTK_FDIM_(n), float *t, SimTK_FDIM_(ldt), float *vl, SimTK_FDIM_(ldvl), float *vr, SimTK_FDIM_(ldvr), float *s, float *sep, SimTK_FDIM_(mm), SimTK_FDIM_(m), float *work, SimTK_FDIM_(ldwork), int *iwork, SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(howmny));

void strsyl_(SimTK_FOPT_(trana), SimTK_FOPT_(tranb), int *isgn, SimTK_FDIM_(m), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), float *c__, SimTK_FDIM_(ldc), SimTK_S_OUTPUT_(scale), SimTK_INFO_, SimTK_FLEN_(trana), SimTK_FLEN_(tranb));

void strti2_(SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

void strtri_(SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

void strtrs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const float *a, SimTK_FDIM_(lda), float *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));

void stzrqf_(SimTK_FDIM_(m), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *tau, SimTK_INFO_);

void stzrzf_(SimTK_FDIM_(m), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *tau, float *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void xerbla_(const char *srname, SimTK_INFO_, SimTK_FLEN_(srname));

void zbdsqr_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(ncvt), SimTK_FDIM_(nru), SimTK_FDIM_(ncc), double *d__, double *e, SimTK_Z_ *vt, SimTK_FDIM_(ldvt), SimTK_Z_ *u, SimTK_FDIM_(ldu), SimTK_Z_ *c__, SimTK_FDIM_(ldc), double *rwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void zcgesv_( SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_Z_ *a, SimTK_FDIM_(lda), int *ipiv, const SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *x,  SimTK_FDIM_(ldx), SimTK_Z_ *work, SimTK_C_ *swork, SimTK_I_OUTPUT_(iter), SimTK_INFO_ );

void zdrscl_(SimTK_FDIM_(n), SimTK_D_INPUT_(sa), SimTK_Z_ *sx, SimTK_FINC_(x));

void zgbbrd_(SimTK_FOPT_(vect), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(ncc), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_Z_ *ab, SimTK_FDIM_(ldab), double *d__, double *e, SimTK_Z_ *q, SimTK_FDIM_(ldq), SimTK_Z_ *pt, SimTK_FDIM_(ldpt), SimTK_Z_ *c__, SimTK_FDIM_(ldc), SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(vect));

void zgbcon_(SimTK_FOPT_(norm), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_Z_ *ab, SimTK_FDIM_(ldab), const int *ipiv, SimTK_D_INPUT_(anorm), SimTK_D_OUTPUT_(rcond), SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(norm));

void zgbequ_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_Z_ *ab, SimTK_FDIM_(ldab), double *r__, double *c__, double *rowcnd, double *colcnd, double *amax, SimTK_INFO_);

void zgbrfs_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_FDIM_(nrhs), SimTK_Z_ *ab, SimTK_FDIM_(ldab), SimTK_Z_ *afb, SimTK_FDIM_(ldafb), const int *ipiv, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *x, SimTK_FDIM_(ldx), double *ferr, double *berr, SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(trans));

void zgbsv_(SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_FDIM_(nrhs), SimTK_Z_ *ab, SimTK_FDIM_(ldab), int *ipiv, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_INFO_);

void zgbsvx_(SimTK_FOPT_(fact), SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_FDIM_(nrhs), SimTK_Z_ *ab, SimTK_FDIM_(ldab), SimTK_Z_ *afb, SimTK_FDIM_(ldafb), int *ipiv, SimTK_CHAR_OUTPUT_(equed), double *r__, double *c__, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *x, SimTK_FDIM_(ldx), SimTK_D_OUTPUT_(rcond), double *ferr, double *berr, SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(trans), SimTK_FLEN_(equed));

void zgbtf2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_Z_ *ab, SimTK_FDIM_(ldab), int *ipiv, SimTK_INFO_);

void zgbtrf_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_Z_ *ab, SimTK_FDIM_(ldab), int *ipiv, SimTK_INFO_);

void zgbtrs_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_FDIM_(nrhs), SimTK_Z_ *ab, SimTK_FDIM_(ldab), const int *ipiv, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(trans));

void zgebak_(SimTK_FOPT_(job), SimTK_FOPT_(side), SimTK_FDIM_(n),  SimTK_I_INPUT_(ilo),  SimTK_I_INPUT_(ihi),  SimTK_D_INPUT_(scale),  SimTK_I_INPUT_(m), SimTK_Z_ *v, SimTK_FDIM_(ldv), SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(side));

void zgebal_(SimTK_FOPT_(job), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_I_OUTPUT_(ilo), SimTK_I_OUTPUT_(ihi), SimTK_D_OUTPUT_(scale), SimTK_INFO_, SimTK_FLEN_(job));

void zgebd2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), double *d__, double *e, SimTK_Z_ *tauq, SimTK_Z_ *taup, SimTK_Z_ *work, SimTK_INFO_);

void zgebrd_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), double *d__, double *e, SimTK_Z_ *tauq, SimTK_Z_ *taup, SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void zgecon_(SimTK_FOPT_(norm), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_D_INPUT_(anorm), SimTK_D_OUTPUT_(rcond), SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(norm));

void zgeequ_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), double *r__, double *c__, double *rowcnd, double *colcnd, double *amax, SimTK_INFO_);

void zgees_(SimTK_FOPT_(jobvs), SimTK_FOPT_(sort), SimTK_SELECT_Z select, SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), int *sdim, SimTK_Z_ *w, SimTK_Z_ *vs, SimTK_FDIM_(ldvs), SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, int *bwork, SimTK_INFO_, SimTK_FLEN_(jobvs), SimTK_FLEN_(sort));

void zgeesx_(SimTK_FOPT_(jobvs), SimTK_FOPT_(sort), SimTK_SELECT_Z select, char *sense, SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), int *sdim, SimTK_Z_ *w, SimTK_Z_ *vs, SimTK_FDIM_(ldvs), SimTK_D_OUTPUT_(rconde), SimTK_D_OUTPUT_(rcondv), SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, int *bwork, SimTK_INFO_, SimTK_FLEN_(jobvs), SimTK_FLEN_(sort), SimTK_FLEN_(sense));

void zgeev_(SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *w, SimTK_Z_ *vl, SimTK_FDIM_(ldvl), SimTK_Z_ *vr, SimTK_FDIM_(ldvr), SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, SimTK_INFO_, SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr));

void zgeevx_(SimTK_FOPT_(balanc), SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), char *sense, SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *w, SimTK_Z_ *vl, SimTK_FDIM_(ldvl), SimTK_Z_ *vr, SimTK_FDIM_(ldvr), SimTK_I_OUTPUT_(ilo), SimTK_I_OUTPUT_(ihi), SimTK_D_OUTPUT_(scale), double *abnrm, SimTK_D_OUTPUT_(rconde), SimTK_D_OUTPUT_(rcondv), SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, SimTK_INFO_, SimTK_FLEN_(balanc), SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr), SimTK_FLEN_(sense));

void zgegs_(SimTK_FOPT_(jobvsl), SimTK_FOPT_(jobvsr), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *alpha, SimTK_Z_ *beta, SimTK_Z_ *vsl, SimTK_FDIM_(ldvsl), SimTK_Z_ *vsr, SimTK_FDIM_(ldvsr), SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, SimTK_INFO_, SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr));

void zgegv_(SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *alpha, SimTK_Z_ *beta, SimTK_Z_ *vl, SimTK_FDIM_(ldvl), SimTK_Z_ *vr, SimTK_FDIM_(ldvr), SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, SimTK_INFO_, SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr));

void zgehd2_(SimTK_FDIM_(n), SimTK_FDIM_(ilo), SimTK_I_INPUT_(ihi), SimTK_Z_ *a, SimTK_I_INPUT_(lda), SimTK_Z_ *tau, SimTK_Z_ *work, SimTK_INFO_);

void zgehrd_(SimTK_FDIM_(n), SimTK_FDIM_(ilo), SimTK_I_INPUT_(ihi), SimTK_Z_ *a, SimTK_I_INPUT_(lda), SimTK_Z_ *tau, SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void zgelq2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *tau, SimTK_Z_ *work, SimTK_INFO_);

void zgelqf_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *tau, SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void zgels_(SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(trans));
void zgelss_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), double *s, SimTK_D_INPUT_(rcond), SimTK_I_OUTPUT_(rank), SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, SimTK_INFO_);

/* gelsx is deprecated; use gelsy instead */
void zgelsx_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), int *jpvt, SimTK_D_INPUT_(rcond), SimTK_I_OUTPUT_(rank), SimTK_Z_ *work, double *rwork, SimTK_INFO_);
void zgelsy_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), int *jpvt, SimTK_D_INPUT_(rcond), SimTK_I_OUTPUT_(rank), SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, SimTK_INFO_);

void zgeql2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *tau, SimTK_Z_ *work, SimTK_INFO_);

void zgeqlf_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *tau, SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void zgeqp3_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), int *jpvt, SimTK_Z_ *tau, SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, SimTK_INFO_);

void zgeqpf_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), int *jpvt, SimTK_Z_ *tau, SimTK_Z_ *work, double *rwork, SimTK_INFO_);

void zgeqr2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *tau, SimTK_Z_ *work, SimTK_INFO_);

void zgeqrf_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *tau, SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void zgerfs_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *af, SimTK_FDIM_(ldaf), const int *ipiv, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *x, SimTK_FDIM_(ldx), double *ferr, double *berr, SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(trans));

void zgerq2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *tau, SimTK_Z_ *work, SimTK_INFO_);

void zgerqf_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *tau, SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void zgesc2_(SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *rhs, int *ipiv, int *jpiv, SimTK_D_OUTPUT_(scale));

void zgesv_(SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_Z_ *a, SimTK_FDIM_(lda), int *ipiv, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_INFO_);

void zgesvx_(SimTK_FOPT_(fact), SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *af, SimTK_FDIM_(ldaf), int *ipiv, SimTK_CHAR_OUTPUT_(equed), double *r__, double *c__, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *x, SimTK_FDIM_(ldx), SimTK_D_OUTPUT_(rcond), double *ferr, double *berr, SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(trans), SimTK_FLEN_(equed));

void zgetc2_(SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), int *ipiv, int *jpiv, SimTK_INFO_);

void zgetf2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), int *ipiv, SimTK_INFO_);

void zgetrf_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), int *ipiv, SimTK_INFO_);

void zgetri_(SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), const int *ipiv, SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void zgetrs_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_Z_ *a, SimTK_FDIM_(lda), const int *ipiv, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(trans));

void zggbak_(SimTK_FOPT_(job), SimTK_FOPT_(side), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), const double *lscale, const double *rscale, SimTK_FDIM_(m), SimTK_Z_ *v, SimTK_FDIM_(ldv), SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(side));

void zggbal_(SimTK_FOPT_(job), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_I_OUTPUT_(ilo), SimTK_I_OUTPUT_(ihi), double *lscale, double *rscale, double *work, SimTK_INFO_, SimTK_FLEN_(job));

void zgges_(SimTK_FOPT_(jobvsl), SimTK_FOPT_(jobvsr), SimTK_FOPT_(sort), SimTK_SELECT_2Z delctg, SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), int *sdim, SimTK_Z_ *alpha, SimTK_Z_ *beta, SimTK_Z_ *vsl, SimTK_FDIM_(ldvsl), SimTK_Z_ *vsr, SimTK_FDIM_(ldvsr), SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, int *bwork, SimTK_INFO_, SimTK_FLEN_(jobvsl), SimTK_FLEN_(jobvsr), SimTK_FLEN_(sort));

void zggesx_(SimTK_FOPT_(jobvsl), SimTK_FOPT_(jobvsr), SimTK_FOPT_(sort), SimTK_SELECT_2Z delctg, SimTK_FOPT_(sense), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_I_OUTPUT_(sdim), SimTK_Z_ *alpha, SimTK_Z_ *beta, SimTK_Z_ *vsl, SimTK_FDIM_(ldvsl), SimTK_Z_ *vsr, SimTK_FDIM_(ldvsr), SimTK_D_OUTPUT_(rconde), SimTK_D_OUTPUT_(rcondv), SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, int *iwork, SimTK_FDIM_(liwork), int *bwork, SimTK_INFO_, SimTK_FLEN_(jobvsl), SimTK_FLEN_(jobvsr), SimTK_FLEN_(sort), SimTK_FLEN_(sense));

void zggev_(SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *alpha, SimTK_Z_ *beta, SimTK_Z_ *vl, SimTK_FDIM_(ldvl), SimTK_Z_ *vr, SimTK_FDIM_(ldvr), SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, SimTK_INFO_, SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr));

void zggevx_(SimTK_FOPT_(balanc), SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), SimTK_FOPT_(sense), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *alpha, SimTK_Z_ *beta, SimTK_Z_ *vl, SimTK_FDIM_(ldvl), SimTK_Z_ *vr, SimTK_FDIM_(ldvr), SimTK_I_OUTPUT_(ilo), SimTK_I_OUTPUT_(ihi), double *lscale, double *rscale, SimTK_D_OUTPUT_(abnrm), SimTK_D_OUTPUT_(bbnrm), SimTK_D_OUTPUT_(rconde), SimTK_D_OUTPUT_(rcondv), SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, int *iwork, int *bwork, SimTK_INFO_, SimTK_FLEN_(blanc), SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr), SimTK_FLEN_(sense));

void zggglm_(SimTK_FDIM_(n), SimTK_FDIM_(m), SimTK_FDIM_(p), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *d__, SimTK_Z_ *x, SimTK_Z_ *y, SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void zgghrd_(SimTK_FOPT_(compq), SimTK_FOPT_(compz), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *q, SimTK_FDIM_(ldq), SimTK_Z_ *z__, SimTK_FDIM_(ldz), SimTK_INFO_, SimTK_FLEN_(compq), SimTK_FLEN_(compz));

void zgglse_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(p), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *c__, SimTK_Z_ *d__, SimTK_Z_ *x, SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void zggqrf_(SimTK_FDIM_(n), SimTK_FDIM_(m), SimTK_FDIM_(p), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *taua, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *taub, SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void zggrqf_(SimTK_FDIM_(m), SimTK_FDIM_(p), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *taua, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *taub, SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void zggsvd_(SimTK_FOPT_(jobu), SimTK_FOPT_(jobv), SimTK_FOPT_(jobq), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(p), SimTK_I_OUTPUT_(k), SimTK_I_OUTPUT_(l), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), double *alpha, double *beta, SimTK_Z_ *u, SimTK_FDIM_(ldu), SimTK_Z_ *v, SimTK_FDIM_(ldv), SimTK_Z_ *q, SimTK_FDIM_(ldq), SimTK_Z_ *work, double *rwork, int *iwork, SimTK_INFO_, SimTK_FLEN_(jobu), SimTK_FLEN_(jobv), SimTK_FLEN_(jobq));

void zggsvp_(SimTK_FOPT_(jobu), SimTK_FOPT_(jobv), SimTK_FOPT_(jobq), SimTK_FDIM_(m), SimTK_FDIM_(p), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_D_INPUT_(tola), SimTK_D_INPUT_(tolb),  SimTK_I_OUTPUT_(k), SimTK_I_OUTPUT_(l), SimTK_Z_ *u, SimTK_FDIM_(ldu), SimTK_Z_ *v, SimTK_FDIM_(ldv), SimTK_Z_ *q, SimTK_FDIM_(ldq), int *iwork, double *rwork, SimTK_Z_ *tau, SimTK_Z_ *work, SimTK_INFO_, SimTK_FLEN_(jobu), SimTK_FLEN_(jobv), SimTK_FLEN_(jobq));

void zgtcon_(SimTK_FOPT_(norm), SimTK_FDIM_(n), SimTK_Z_ *dl, SimTK_Z_ *d__, SimTK_Z_ *du, SimTK_Z_ *du2, const int *ipiv, SimTK_D_INPUT_(anorm),SimTK_D_OUTPUT_(rcond), SimTK_Z_ *work, SimTK_INFO_, SimTK_FLEN_(norm));

void zgtrfs_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_Z_ *dl, SimTK_Z_ *d__, SimTK_Z_ *du, SimTK_Z_ *dlf, SimTK_Z_ *df, SimTK_Z_ *duf, SimTK_Z_ *du2, const int *ipiv, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *x, SimTK_FDIM_(ldx), double *ferr, double *berr, SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(trans));

void zgtsv_(SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_Z_ *dl, SimTK_Z_ *d__, SimTK_Z_ *du, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_INFO_);

void zgtsvx_(SimTK_FOPT_(fact), SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_Z_ *dl, SimTK_Z_ *d__, SimTK_Z_ *du, SimTK_Z_ *dlf, SimTK_Z_ *df, SimTK_Z_ *duf, SimTK_Z_ *du2, int *ipiv, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *x, SimTK_FDIM_(ldx), SimTK_D_OUTPUT_(rcond), double *ferr, double *berr, SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(trans));

void zgttrf_(SimTK_FDIM_(n), SimTK_Z_ *dl, SimTK_Z_ *d__, SimTK_Z_ *du, SimTK_Z_ *du2, int *ipiv, SimTK_INFO_);

void zgttrs_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_Z_ *dl, SimTK_Z_ *d__, SimTK_Z_ *du, SimTK_Z_ *du2, const int *ipiv, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(trans));

void zgtts2_(int *itrans, SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_Z_ *dl, SimTK_Z_ *d__, SimTK_Z_ *du, SimTK_Z_ *du2, const int *ipiv, SimTK_Z_ *b, SimTK_FDIM_(ldb));

void zhbev_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_Z_ *ab, SimTK_FDIM_(ldab), double *w, SimTK_Z_ *z__, SimTK_FDIM_(ldz), SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void zhbevd_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_Z_ *ab, SimTK_FDIM_(ldab), double *w, SimTK_Z_ *z__, SimTK_FDIM_(ldz), SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, SimTK_FDIM_(lrwork), int *iwork, SimTK_FDIM_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void zhbevx_(SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_Z_ *ab, SimTK_FDIM_(ldab), SimTK_Z_ *q, SimTK_FDIM_(ldq), SimTK_D_INPUT_(vl), SimTK_D_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_D_INPUT_(abstol), SimTK_FDIM_(m), double *w, SimTK_Z_ *z__, SimTK_FDIM_(ldz), SimTK_Z_ *work, double *rwork, int *iwork, int *ifail, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void zhbgst_(SimTK_FOPT_(vect), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(ka), SimTK_FDIM_(kb), SimTK_Z_ *ab, SimTK_FDIM_(ldab), SimTK_Z_ *bb, SimTK_FDIM_(ldbb), SimTK_Z_ *x, SimTK_FDIM_(ldx), SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(vect), SimTK_FLEN_(uplo));

void zhbgv_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(ka), SimTK_FDIM_(kb), SimTK_Z_ *ab, SimTK_FDIM_(ldab), SimTK_Z_ *bb, SimTK_FDIM_(ldbb), double *w, SimTK_Z_ *z__, SimTK_FDIM_(ldz), SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void zhbgvx_(SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(ka), SimTK_FDIM_(kb), SimTK_Z_ *ab, SimTK_FDIM_(ldab), SimTK_Z_ *bb, SimTK_FDIM_(ldbb), SimTK_Z_ *q, SimTK_FDIM_(ldq), SimTK_D_INPUT_(vl), SimTK_D_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_D_INPUT_(abstol), SimTK_I_OUTPUT_(m), double *w, SimTK_Z_ *z__, SimTK_FDIM_(ldz), SimTK_Z_ *work, double *rwork, int *iwork, int *ifail, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void zhbtrd_(SimTK_FOPT_(vect), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_Z_ *ab, SimTK_FDIM_(ldab), double *d__, double *e, SimTK_Z_ *q, SimTK_FDIM_(ldq), SimTK_Z_ *work, SimTK_INFO_, SimTK_FLEN_(vect), SimTK_FLEN_(uplo));

void zhecon_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), const int *ipiv, SimTK_D_INPUT_(anorm), SimTK_D_OUTPUT_(rcond), SimTK_Z_ *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void zheev_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), double *w, SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void zheevd_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), double *w, SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, SimTK_FDIM_(lrwork), int *iwork, SimTK_FDIM_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void zheevr_(SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_D_INPUT_(vl), SimTK_D_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_D_INPUT_(abstol),  SimTK_I_OUTPUT_(m), double *w, SimTK_Z_ *z__, SimTK_FDIM_(ldz), int *isuppz, SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, SimTK_FDIM_(lrwork), int *iwork, SimTK_FDIM_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void zheevx_(SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_D_INPUT_(vl), SimTK_D_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_D_INPUT_(abstol),  SimTK_I_OUTPUT_(m), double *w, SimTK_Z_ *z__, SimTK_FDIM_(ldz), SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, int *iwork, int *ifail, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void zhegs2_(SimTK_I_INPUT_(itype), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void zhegst_(SimTK_I_INPUT_(itype), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void zhegv_(SimTK_I_INPUT_(itype), SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), double *w, SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void zhegvd_(SimTK_I_INPUT_(itype), SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), double *w, SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, int *lrwork, int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void zhegvx_(SimTK_I_INPUT_(itype), SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_D_INPUT_(vl), SimTK_D_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_D_INPUT_(abstol),  SimTK_I_OUTPUT_(m), double *w, SimTK_Z_ *z__, SimTK_FDIM_(ldz), SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, int *iwork, int *ifail, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void zherfs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *af, SimTK_FDIM_(ldaf), const int *ipiv, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *x, SimTK_FDIM_(ldx), double *ferr, double *berr, SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void zhesv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_Z_ *a, SimTK_FDIM_(lda), const int *ipiv, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(uplo));

void zhesvx_(SimTK_FOPT_(fact), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *af, SimTK_FDIM_(ldaf), int *ipiv, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *x, SimTK_FDIM_(ldx), SimTK_D_OUTPUT_(rcond), double *ferr, double *berr, SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(uplo));

void zhetf2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), int *ipiv, SimTK_INFO_, SimTK_FLEN_(uplo));

void zhetrd_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), double *d__, double *e, SimTK_Z_ *tau, SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(uplo));

void zhetrf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), int *ipiv, SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(uplo));

void zhetri_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), const int *ipiv, SimTK_Z_ *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void zhetrs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_Z_ *a, SimTK_FDIM_(lda), const int *ipiv, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void zhgeqz_(SimTK_FOPT_(job), SimTK_FOPT_(compq), SimTK_FOPT_(compz), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *alpha, SimTK_Z_ *beta, SimTK_Z_ *q, SimTK_FDIM_(ldq), SimTK_Z_ *z__, SimTK_FDIM_(ldz), SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(compq), SimTK_FLEN_(compz));

void zhpcon_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const SimTK_Z_ *ap, const int *ipiv, SimTK_D_INPUT_(anorm), SimTK_D_OUTPUT_(rcond), SimTK_Z_ *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void zhpev_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *ap, double *w, SimTK_Z_ *z__, SimTK_FDIM_(ldz), SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void zhpevd_(SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *ap, double *w, SimTK_Z_ *z__, SimTK_FDIM_(ldz), SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, SimTK_FDIM_(lrwork), int *iwork, SimTK_FDIM_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void zhpevx_(SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *ap, SimTK_D_INPUT_(vl), SimTK_D_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_D_INPUT_(abstol), SimTK_I_OUTPUT_(m), double *w, SimTK_Z_ *z__, SimTK_FDIM_(ldz), SimTK_Z_ *work, double *rwork, int *iwork, int *ifail, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void zhpgst_(SimTK_I_INPUT_(itype), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *ap, const SimTK_Z_ *bp, SimTK_INFO_, SimTK_FLEN_(uplo));

void zhpgv_(SimTK_I_INPUT_(itype), SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *ap, SimTK_Z_ *bp, double *w, SimTK_Z_ *z__, SimTK_FDIM_(ldz), SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void zhpgvd_(SimTK_I_INPUT_(itype), SimTK_FOPT_(jobz), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *ap, SimTK_Z_ *bp, double *w, SimTK_Z_ *z__, SimTK_FDIM_(ldz), SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, SimTK_FDIM_(lrwork), int *iwork, SimTK_FDIM_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(uplo));

void zhpgvx_(SimTK_I_INPUT_(itype), SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *ap, SimTK_Z_ *bp, SimTK_D_INPUT_(vl), SimTK_D_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_D_INPUT_(abstol), SimTK_I_OUTPUT_(m), double *w, SimTK_Z_ *z__, SimTK_FDIM_(ldz), SimTK_Z_ *work, double *rwork, int *iwork, int *ifail, SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range), SimTK_FLEN_(uplo));

void zhprfs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_Z_ *ap, const SimTK_Z_ *afp, const int *ipiv, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *x, SimTK_FDIM_(ldx), double *ferr, double *berr, SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void zhpsv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_Z_ *ap, int *ipiv, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void zhpsvx_(SimTK_FOPT_(fact), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_Z_ *ap, const SimTK_Z_ *afp, const int *ipiv, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *x, SimTK_FDIM_(ldx), SimTK_D_OUTPUT_(rcond), double *ferr, double *berr, SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(uplo));

void zhptrd_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *ap, double *d__, double *e, SimTK_Z_ *tau, SimTK_INFO_, SimTK_FLEN_(uplo));

void zhptrf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *ap, int *ipiv, SimTK_INFO_, SimTK_FLEN_(uplo));

void zhptri_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *ap, const int *ipiv, SimTK_Z_ *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void zhptrs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_Z_ *ap, const int *ipiv, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void zhsein_(SimTK_FOPT_(side), SimTK_FOPT_(eigsrc), SimTK_FOPT_(initv), const int *select, SimTK_FDIM_(n), const SimTK_Z_ *h__, SimTK_FDIM_(ldh), SimTK_Z_ *w, SimTK_Z_ *vl, SimTK_FDIM_(ldvl), SimTK_Z_ *vr, SimTK_FDIM_(ldvr), SimTK_FDIM_(mm), SimTK_I_OUTPUT_(m), SimTK_Z_ *work, double *rwork, int *ifaill, int *ifailr, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(eigsrc), SimTK_FLEN_(initv));

void zhseqr_(SimTK_FOPT_(job), SimTK_FOPT_(compz), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), SimTK_Z_ *h__, SimTK_FDIM_(ldh), SimTK_Z_ *w, SimTK_Z_ *z__, SimTK_FDIM_(ldz), SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(compz));

void zlabrd_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(nb), SimTK_Z_ *a, SimTK_FDIM_(lda), double *d__, double *e, SimTK_Z_ *tauq, SimTK_Z_ *taup, SimTK_Z_ *x, SimTK_FDIM_(ldx), SimTK_Z_ *y, SimTK_FDIM_(ldy));

void zlacgv_(SimTK_FDIM_(n), SimTK_Z_ *x, SimTK_FINC_(x));

void zlacon_(SimTK_FDIM_(n), SimTK_Z_ *v, SimTK_Z_ *x, SimTK_D_OUTPUT_(est), SimTK_I_OUTPUT_(kase) );

void zlacn2_( SimTK_FDIM_(n), SimTK_Z_ *v, SimTK_Z_ *x, SimTK_D_OUTPUT_(est), SimTK_I_OUTPUT_(kase), int *isave   );

void zlacp2_(SimTK_FOPT_(uplo), SimTK_FDIM_(m), SimTK_FDIM_(n), const double *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_FLEN_(uplo));

void zlacpy_(SimTK_FOPT_(uplo), SimTK_FDIM_(m), SimTK_FDIM_(n), const SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_FLEN_(uplo));

void zlacrm_(SimTK_FDIM_(m), SimTK_FDIM_(n), const SimTK_Z_ *a, SimTK_FDIM_(lda), const double *b, SimTK_FDIM_(ldb), const SimTK_Z_ *c__, SimTK_FDIM_(ldc), double *rwork);

void zlacrt_(SimTK_FDIM_(n), SimTK_Z_ *cx, SimTK_FINC_(x), SimTK_Z_ *cy, SimTK_FINC_(y), const SimTK_Z_ *c__, const SimTK_Z_ *s);

void zlaed0_(SimTK_I_INPUT_(qsiz), SimTK_FDIM_(n), double *d__, const double *e, SimTK_Z_ *q, SimTK_FDIM_(ldq), SimTK_Z_ *qstore, SimTK_FDIM_(ldqs), double *rwork, int *iwork, SimTK_INFO_);

void zlaed7_(SimTK_FDIM_(n), SimTK_I_INPUT_(cutpnt), SimTK_FDIM_(qsiz), SimTK_I_INPUT_(tlvls), SimTK_I_INPUT_(curlvl), SimTK_I_INPUT_(curpbm), double *d__, SimTK_Z_ *q, SimTK_FDIM_(ldq), SimTK_D_INPUT_(rho), int *indxq, double *qstore, int *qptr, const int *prmptr, const int *perm, const int *givptr, const int *givcol, const double *givnum, SimTK_Z_ *work, double *rwork, int *iwork, SimTK_INFO_);

void zlaed8_(SimTK_I_OUTPUT_(k), SimTK_FDIM_(n), SimTK_FDIM_(qsiz), SimTK_Z_ *q, SimTK_FDIM_(ldq), double *d__, SimTK_D_OUTPUT_(rho), SimTK_I_INPUT_(cutpnt), double *z__, double *dlamda, SimTK_Z_ *q2, SimTK_FDIM_(ldq2), double *w, int *indxp, int *indx, const int *indxq, int *perm, int *givptr, int *givcol, double *givnum, SimTK_INFO_);

void zlaein_(SimTK_I_INPUT_(rightv), SimTK_I_INPUT_(noinit), SimTK_FDIM_(n), const SimTK_Z_ *h__, SimTK_FDIM_(ldh), const SimTK_Z_ *w, SimTK_Z_ *v, SimTK_Z_ *b, SimTK_FDIM_(ldb), double *rwork, SimTK_D_INPUT_(eps3), SimTK_D_INPUT_(smlnum), SimTK_INFO_);

void zlaesy_(SimTK_Z_INPUT_(a), SimTK_Z_INPUT_(b), SimTK_Z_INPUT_(c__), SimTK_Z_OUTPUT_(rt1), SimTK_Z_OUTPUT_(rt2), SimTK_Z_OUTPUT_(evscal), SimTK_Z_OUTPUT_(cs1), SimTK_Z_OUTPUT_(sn1));

void zlaev2_(SimTK_Z_INPUT_(a), SimTK_Z_INPUT_(b), SimTK_Z_INPUT_(c__), SimTK_D_OUTPUT_(rt1), SimTK_D_OUTPUT_(rt2), SimTK_D_OUTPUT_(cs1), SimTK_Z_OUTPUT_(sn1));

void zlags2_(SimTK_I_INPUT_(upper), SimTK_D_INPUT_(a1), SimTK_Z_INPUT_(a2), SimTK_D_INPUT_(a3), SimTK_D_INPUT_(b1), SimTK_Z_INPUT_(b2), SimTK_D_INPUT_(b3), SimTK_D_OUTPUT_(csu), SimTK_Z_OUTPUT_(snu), SimTK_D_OUTPUT_(csv), SimTK_Z_OUTPUT_(snv), SimTK_D_OUTPUT_(csq), SimTK_Z_OUTPUT_(snq) );

void zlagtm_(SimTK_FOPT_(trans), SimTK_FDIM_(n), SimTK_FDIM_(nrhs),  SimTK_D_INPUT_(alpha), const SimTK_Z_ *dl, const SimTK_Z_ *d__, const SimTK_Z_ *du, const SimTK_Z_ *x, SimTK_FDIM_(ldx),  SimTK_D_INPUT_(beta), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_FLEN_(trans));

void zlag2c_( SimTK_I_INPUT_(m), SimTK_I_INPUT_(n),  const SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_C_ *sa, SimTK_FDIM_(ldsa),  SimTK_INFO_);

void zlahef_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nb), SimTK_FDIM_(kb), SimTK_Z_ *a, SimTK_FDIM_(lda), int *ipiv, SimTK_Z_ *w, SimTK_FDIM_(ldw), SimTK_INFO_, SimTK_FLEN_(uplo));

void zlahqr_(SimTK_I_INPUT_(wantt), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), SimTK_Z_ *h__, SimTK_FDIM_(ldh), SimTK_Z_ *w, SimTK_I_INPUT_(iloz), SimTK_I_INPUT_(ihiz), SimTK_Z_ *z__, SimTK_FDIM_(ldz), SimTK_INFO_);

void zlahrd_(SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_FDIM_(nb), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *tau, SimTK_Z_ *t, SimTK_FDIM_(ldt), SimTK_Z_ *y, SimTK_FDIM_(ldy));

void zlahr2_(SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_FDIM_(nb), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *tau, SimTK_Z_ *t, SimTK_FDIM_(ldt), SimTK_Z_ *y, SimTK_FDIM_(ldy));

void zlaic1_(SimTK_I_INPUT_(job), SimTK_FDIM_(j), const SimTK_Z_ *x, SimTK_D_INPUT_(sest), const SimTK_Z_ *w, SimTK_Z_INPUT_(gamma), SimTK_D_OUTPUT_(sestpr), SimTK_Z_OUTPUT_(s), SimTK_Z_OUTPUT_(c__));

void zlals0_(SimTK_I_INPUT_(icompq), SimTK_I_INPUT_(nl), SimTK_I_INPUT_(nr), SimTK_I_INPUT_(sqre), SimTK_FDIM_(nrhs), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *bx, SimTK_FDIM_(ldbx), int *perm, int *givptr, const int *givcol, SimTK_FDIM_(ldgcol), double *givnum, SimTK_FDIM_(ldgnum), double *poles, double *difl, double *difr, double *z__, SimTK_FDIM_(k), double *c__, double *s, double *rwork, SimTK_INFO_);

void zlalsa_(SimTK_I_INPUT_(icompq), SimTK_I_INPUT_(smlsiz), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *bx, SimTK_FDIM_(ldbx), double *u, SimTK_FDIM_(ldu), double *vt, int *k, double *difl, double *difr, double *z__, double *poles, int *givptr, const int *givcol, SimTK_FDIM_(ldgcol), int *perm, double *givnum, double *c__, double *s, double *rwork, int *iwork, SimTK_INFO_);

void zlapll_(SimTK_FDIM_(n), SimTK_Z_ *x, SimTK_FINC_(x), SimTK_Z_ *y, SimTK_FINC_(y), SimTK_D_OUTPUT_(ssmin));

void zlapmt_(SimTK_I_INPUT_(forwrd), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_Z_ *x, SimTK_FDIM_(ldx), SimTK_FDIM_(k));

void zlaqgb_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(kl), SimTK_FDIM_(ku), SimTK_Z_ *ab, SimTK_FDIM_(ldab), double *r__, double *c__, SimTK_D_OUTPUT_(rowcnd), SimTK_D_OUTPUT_(colcnd), SimTK_D_INPUT_(amax), SimTK_CHAR_OUTPUT_(equed),  SimTK_FLEN_(equed) );

void zlaqge_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), double *r__, double *c__, SimTK_D_INPUT_(rowcnd), SimTK_D_INPUT_(colcnd), SimTK_D_INPUT_(amax),  SimTK_CHAR_OUTPUT_(equed), SimTK_FLEN_(equed) );

void zlaqhb_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_Z_ *ab, SimTK_FDIM_(ldab), double *s, SimTK_D_INPUT_(scond), SimTK_D_INPUT_(amax), SimTK_CHAR_OUTPUT_(equed), SimTK_FLEN_(uplo), SimTK_FLEN_(equed) );

void zlaqhe_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), const double *s, SimTK_D_INPUT_(scond), SimTK_D_INPUT_(amax), SimTK_CHAR_OUTPUT_(equed), SimTK_FLEN_(uplo), SimTK_FLEN_(equed));

void zlaqhp_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *ap, const double *s, SimTK_D_INPUT_(scond), SimTK_D_INPUT_(amax), SimTK_CHAR_OUTPUT_(equed), SimTK_FLEN_(uplo), SimTK_FLEN_(equed));

void zlaqp2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_I_INPUT_(offset), SimTK_Z_ *a, SimTK_FDIM_(lda), int *jpvt, SimTK_Z_ *tau, double *vn1, double *vn2, SimTK_Z_ *work);

void zlaqps_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_I_INPUT_(offset), SimTK_I_INPUT_(nb), SimTK_I_OUTPUT_(kb), SimTK_Z_ *a, SimTK_FDIM_(lda), int *jpvt, SimTK_Z_ *tau, double *vn1, double *vn2, SimTK_Z_ *auxv, SimTK_Z_ *f, SimTK_FDIM_(ldf));

void zlaqr0_( SimTK_I_INPUT_(wantt), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), SimTK_Z_ *h, SimTK_FDIM_(ldh), SimTK_Z_ *w, SimTK_I_INPUT_(iloz), SimTK_I_INPUT_(ihiz), SimTK_Z_ *z, SimTK_FDIM_(ldz), SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_ );

void zlaqr1_(SimTK_FDIM_(n), const SimTK_Z_ *h, SimTK_FDIM_(ldh), SimTK_Z_INPUT_(s1), SimTK_Z_INPUT_(s2), SimTK_Z_ *v );

void zlaqr2_( SimTK_I_INPUT_(wantt), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), SimTK_I_INPUT_(ktop), SimTK_I_INPUT_(kbot),  SimTK_I_INPUT_(nw), SimTK_Z_ *h, SimTK_FDIM_(ldh), SimTK_I_INPUT_(iloz), SimTK_I_INPUT_(ihiz),  SimTK_Z_ *z, SimTK_FDIM_(ldz), SimTK_I_OUTPUT_(ns), SimTK_I_OUTPUT_(nd), SimTK_Z_ *sh, SimTK_Z_ *v, SimTK_FDIM_(ldv), SimTK_I_INPUT_(nh), SimTK_Z_ *t, SimTK_FDIM_(ldt), SimTK_I_INPUT_(nv), SimTK_Z_ *wv, SimTK_I_INPUT_(ldwv), SimTK_Z_ *work, SimTK_FDIM_(lwork) );

void zlaqr3_( SimTK_I_INPUT_(wantt), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), SimTK_I_INPUT_(ktop), SimTK_I_INPUT_(kbot),  SimTK_I_INPUT_(nw), SimTK_Z_ *h, SimTK_FDIM_(ldh), SimTK_I_INPUT_(iloz), SimTK_I_INPUT_(ihiz),  SimTK_Z_ *z, SimTK_FDIM_(ldz), SimTK_I_OUTPUT_(ns), SimTK_I_OUTPUT_(nd), SimTK_Z_ *sh, SimTK_Z_ *v, SimTK_FDIM_(ldv), SimTK_I_INPUT_(nh), SimTK_Z_ *t, SimTK_FDIM_(ldt), SimTK_I_INPUT_(nv), SimTK_Z_ *wv, SimTK_I_INPUT_(ldwv), SimTK_Z_ *work, SimTK_FDIM_(lwork) );

void zlaqr4_( SimTK_I_INPUT_(wantt), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), SimTK_Z_ *h, SimTK_FDIM_(ldh), SimTK_Z_ *w, SimTK_I_INPUT_(iloz), SimTK_I_INPUT_(ihiz), SimTK_Z_ *z, SimTK_FDIM_(ldz), SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_ );

void zlaqr5_( SimTK_I_INPUT_(wantt), SimTK_I_INPUT_(wantz), SimTK_I_INPUT_(kacc22), SimTK_FDIM_(n), SimTK_I_INPUT_(ktop), SimTK_I_INPUT_(kbot),  SimTK_I_INPUT_(nshifts), const SimTK_Z_ *s, SimTK_Z_ *h, SimTK_FDIM_(ldh), SimTK_I_INPUT_(iloz), SimTK_I_INPUT_(ihiz),  SimTK_Z_ *z, SimTK_FDIM_(ldz), SimTK_Z_ *v, SimTK_FDIM_(ldv), SimTK_Z_ *u, SimTK_FDIM_(ldu), SimTK_I_INPUT_(nh), SimTK_Z_ *wh, SimTK_FDIM_(ldwh), SimTK_I_INPUT_(nv), SimTK_Z_ *wv, SimTK_FDIM_(ldwv)  );

void zlaqsb_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_Z_ *ab, SimTK_FDIM_(ldab), double *s, SimTK_D_INPUT_(scond), SimTK_D_INPUT_(amax), SimTK_CHAR_OUTPUT_(equed), SimTK_FLEN_(uplo), SimTK_FLEN_(equed));

void zlaqsp_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *ap, const double *s, SimTK_D_INPUT_(scond), SimTK_D_INPUT_(amax), SimTK_CHAR_OUTPUT_(equed),  SimTK_FLEN_(uplo), SimTK_FLEN_(equed) );

void zlaqsy_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), const double *s, SimTK_D_INPUT_(scond), SimTK_D_INPUT_(amax), SimTK_CHAR_OUTPUT_(equed), SimTK_FLEN_(uplo), SimTK_FLEN_(equed) );

void zlar1v_(SimTK_FDIM_(n), SimTK_I_INPUT_(b1), SimTK_I_INPUT_(bn), SimTK_D_INPUT_(lambda), const double *d__, const double *l, const double *ld, const double *lld, SimTK_D_INPUT_(pivmin), SimTK_D_INPUT_(gaptol), SimTK_Z_ *z__, SimTK_I_INPUT_(wantnc), SimTK_I_OUTPUT_(negcnt), SimTK_D_OUTPUT_(ztz), SimTK_D_OUTPUT_(mingma), SimTK_I_OUTPUT_(r__), int *isuppz, SimTK_D_OUTPUT_(nrminv), SimTK_D_OUTPUT_(resid), SimTK_D_OUTPUT_(rqcorr), double *work);


void zlar2v_(SimTK_FDIM_(n), SimTK_Z_ *x, SimTK_Z_ *y, SimTK_Z_ *z__, SimTK_FINC_(x), const double *c__, const SimTK_Z_ *s, SimTK_FINC_(c));

void zlarcm_(SimTK_FDIM_(m), SimTK_FDIM_(n), const double *a, SimTK_FDIM_(lda), const SimTK_Z_ *b, SimTK_FDIM_(ldb), const SimTK_Z_ *c__, SimTK_FDIM_(ldc), double *rwork);

void zlarf_(SimTK_FOPT_(side), SimTK_FDIM_(m), SimTK_FDIM_(n), const SimTK_Z_ *v, SimTK_FINC_(v), SimTK_Z_INPUT_(tau), SimTK_Z_ *c__, SimTK_FDIM_(ldc), SimTK_Z_ *work, SimTK_FLEN_(side));

void zlarfb_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FOPT_(direct), SimTK_FOPT_(storev), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_Z_ *v, SimTK_FDIM_(ldv), const SimTK_Z_ *t, SimTK_FDIM_(ldt), SimTK_Z_ *c__, SimTK_FDIM_(ldc), SimTK_Z_ *work, SimTK_FDIM_(ldwork), SimTK_FLEN_(side), SimTK_FLEN_(trans), SimTK_FLEN_(direct), SimTK_FLEN_(storev));

void zlarfg_(SimTK_FDIM_(n), SimTK_Z_OUTPUT_(alpha), SimTK_Z_ *x, SimTK_FINC_(x), SimTK_Z_OUTPUT_(tau) );

void zlarft_(SimTK_FOPT_(direct), SimTK_FOPT_(storev), SimTK_FDIM_(n), SimTK_I_INPUT_(k), SimTK_Z_ *v, SimTK_FDIM_(ldv), const SimTK_Z_ *tau, SimTK_Z_ *t, SimTK_FDIM_(ldt), SimTK_FLEN_(direct), SimTK_FLEN_(storev));

void zlarfx_(SimTK_FOPT_(side), SimTK_FDIM_(m), SimTK_FDIM_(n), const SimTK_Z_ *v, const SimTK_Z_ *tau, SimTK_Z_ *c__, SimTK_FDIM_(ldc), SimTK_Z_ *work, SimTK_FLEN_(side));

void zlargv_(SimTK_FDIM_(n), SimTK_Z_ *x, SimTK_FINC_(x), SimTK_Z_ *y, SimTK_FINC_(y), double *c__, SimTK_FINC_(c));

void zlarnv_(SimTK_I_INPUT_(idist), int *iseed, SimTK_FDIM_(n), SimTK_Z_ *x);

void zlarrv_(SimTK_FDIM_(n), SimTK_D_INPUT_(vl), SimTK_D_INPUT_(vu), double *d__, double *l, SimTK_D_INPUT_(pivmin), const int *isplit, SimTK_FDIM_(m), SimTK_I_INPUT_(dol), SimTK_I_INPUT_(dou), SimTK_D_INPUT_(minrgp), SimTK_D_INPUT_(rtol1), SimTK_D_INPUT_(rtol2), double *w, double *werr, double *wgap, const int *iblock, const int *indexw, const double *gers,  SimTK_Z_ *z__, SimTK_FDIM_(ldz), int *isuppz, double *work, int *iwork, SimTK_INFO_);


void zlartg_(const SimTK_Z_ *f, const SimTK_Z_ *g, double *cs, SimTK_Z_ *sn, SimTK_Z_ *r__);

void zlartv_(SimTK_FDIM_(n), SimTK_Z_ *x, SimTK_FINC_(x), SimTK_Z_ *y, SimTK_FINC_(y), const double *c__, const SimTK_Z_ *s, SimTK_FINC_(c));

void zlarz_(SimTK_FOPT_(side), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_I_INPUT_(L), const SimTK_Z_ *v, SimTK_FINC_(v), SimTK_Z_INPUT_(tau), SimTK_Z_ *c__, SimTK_FDIM_(ldc), SimTK_Z_ *work, SimTK_FLEN_(side));

void zlarzb_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FOPT_(direct), SimTK_FOPT_(storev), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_FDIM_(l), const SimTK_Z_ *v, SimTK_FDIM_(ldv), SimTK_Z_ *t, SimTK_FDIM_(ldt), SimTK_Z_ *c__, SimTK_FDIM_(ldc), SimTK_Z_ *work, SimTK_FDIM_(ldwork), SimTK_FLEN_(side), SimTK_FLEN_(trans), SimTK_FLEN_(direct), SimTK_FLEN_(storev));

void zlarzt_(SimTK_FOPT_(direct), SimTK_FOPT_(storev), SimTK_FDIM_(n),  SimTK_FDIM_(k), SimTK_Z_ *v, SimTK_FDIM_(ldv), const SimTK_Z_ *tau, SimTK_Z_ *t, SimTK_FDIM_(ldt), SimTK_FLEN_(direct), SimTK_FLEN_(storev));

void zlascl_(SimTK_FOPT_(type__), SimTK_FDIM_(kl), SimTK_FDIM_(ku), const double *cfrom, const double *cto, SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(type));

void zlaset_(SimTK_FOPT_(uplo), SimTK_FDIM_(m), SimTK_FDIM_(n),  SimTK_Z_INPUT_(alpha), SimTK_Z_INPUT_(beta), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_FLEN_(uplo));

void zlasr_(SimTK_FOPT_(side), SimTK_FOPT_(pivot), SimTK_FOPT_(direct), SimTK_FDIM_(m), SimTK_FDIM_(n), const double *c__, const double *s, SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_FLEN_(side), SimTK_FLEN_(pivot), SimTK_FLEN_(direct));

void zlassq_(SimTK_FDIM_(n), const SimTK_Z_ *x, SimTK_FINC_(x), SimTK_D_OUTPUT_(scale), double *sumsq);

void zlaswp_(SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_I_OUTPUT_(k1), SimTK_I_OUTPUT_(k2), const int *ipiv, SimTK_FINC_(x));

void zlasyf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_I_INPUT_(nb), SimTK_I_INPUT_(kb), SimTK_Z_ *a, SimTK_FDIM_(lda), int *ipiv, SimTK_Z_ *w, SimTK_FDIM_(ldw), SimTK_INFO_, SimTK_FLEN_(uplo));

void zlatbs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FOPT_(normin), SimTK_FDIM_(n), SimTK_I_INPUT_(kd), const SimTK_Z_ *ab, SimTK_FDIM_(ldab), SimTK_Z_ *x, SimTK_D_OUTPUT_(scale), double *cnorm, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag), SimTK_FLEN_(normin) );

void zlatdf_(SimTK_I_INPUT_(ijob), SimTK_FDIM_(n), const SimTK_Z_ *z__, SimTK_FDIM_(ldz), SimTK_Z_ *rhs, SimTK_D_OUTPUT_(rdsum), SimTK_D_OUTPUT_(rdscal), const int *ipiv, const int *jpiv);

void zlatps_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FOPT_(normin), SimTK_FDIM_(n), const SimTK_Z_ *ap, SimTK_Z_ *x, SimTK_D_OUTPUT_(scale), double *cnorm, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag), SimTK_FLEN_(normin) );

void zlatrd_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nb), SimTK_Z_ *a, SimTK_FDIM_(lda), double *e, SimTK_Z_ *tau, SimTK_Z_ *w, SimTK_FDIM_(ldw), SimTK_FLEN_(uplo));

void zlatrs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FOPT_(normin), SimTK_FDIM_(n), const SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *x, SimTK_D_OUTPUT_(scale), double *cnorm, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag), SimTK_FLEN_(normin));

void zlatrz_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(l), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *tau, SimTK_Z_ *work);

void zlatzm_(SimTK_FOPT_(side), SimTK_FDIM_(m), SimTK_FDIM_(n), const SimTK_Z_ *v, SimTK_FINC_(v), SimTK_Z_INPUT_(tau), SimTK_Z_ *c1, SimTK_Z_ *c2, SimTK_FDIM_(ldc), SimTK_Z_ *work, SimTK_FLEN_(side));

void zlauu2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo));

void zlauum_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo));

void zpbcon_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), const SimTK_Z_ *ab, SimTK_FDIM_(ldab), SimTK_D_INPUT_(anorm), SimTK_D_OUTPUT_(rcond), SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void zpbequ_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_I_INPUT_(kd), const SimTK_Z_ *ab, SimTK_FDIM_(ldab), SimTK_D_OUTPUT_(s), SimTK_D_OUTPUT_(scond), SimTK_D_OUTPUT_(amax), SimTK_INFO_, SimTK_FLEN_(uplo));

void zpbrfs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_FDIM_(nrhs), const SimTK_Z_ *ab, SimTK_FDIM_(ldab), const SimTK_Z_ *afb, SimTK_FDIM_(ldafb), const SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *x, SimTK_FDIM_(ldx), double *ferr, double *berr, SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void zpbstf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_Z_ *ab, SimTK_FDIM_(ldab), SimTK_INFO_, SimTK_FLEN_(uplo));

void zpbsv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_FDIM_(nrhs), SimTK_Z_ *ab, SimTK_FDIM_(ldab), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void zpbsvx_(SimTK_FOPT_(fact), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_FDIM_(nrhs), SimTK_Z_ *ab, SimTK_FDIM_(ldab), SimTK_Z_ *afb, SimTK_FDIM_(ldafb), SimTK_CHAR_OUTPUT_(equed), double *s, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *x, SimTK_FDIM_(ldx), SimTK_D_OUTPUT_(rcond), double *ferr, double *berr, SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(uplo), SimTK_FLEN_(equed));

void zpbtf2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_Z_ *ab, SimTK_FDIM_(ldab), SimTK_INFO_, SimTK_FLEN_(uplo));

void zpbtrf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_Z_ *ab, SimTK_FDIM_(ldab), SimTK_INFO_, SimTK_FLEN_(uplo));

void zpbtrs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_FDIM_(nrhs), const SimTK_Z_ *ab, SimTK_FDIM_(ldab), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void zpocon_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_D_INPUT_(anorm), SimTK_D_OUTPUT_(rcond), SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void zpoequ_(SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), double *s, SimTK_D_OUTPUT_(scond), SimTK_D_OUTPUT_(amax), SimTK_INFO_);

void zporfs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *af, SimTK_FDIM_(ldaf), const SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *x, SimTK_FDIM_(ldx), double *ferr, double *berr, SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void zposv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void zposvx_(SimTK_FOPT_(fact), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *af, SimTK_FDIM_(ldaf), SimTK_CHAR_OUTPUT_(equed), double *s, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *x, SimTK_FDIM_(ldx), SimTK_D_OUTPUT_(rcond), double *ferr, double *berr, SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(uplo), SimTK_FLEN_(equed));

void zpotf2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo));

void zpotrf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo));

void zpotri_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo));

void zpotrs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void zppcon_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const SimTK_Z_ *ap, SimTK_D_INPUT_(anorm), SimTK_D_OUTPUT_(rcond), SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void zppequ_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const SimTK_Z_ *ap, double *s, SimTK_D_OUTPUT_(scond), SimTK_D_OUTPUT_(amax), SimTK_INFO_, SimTK_FLEN_(uplo));

void zpprfs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_Z_ *ap, const SimTK_Z_ *afp, const SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *x, SimTK_FDIM_(ldx), double *ferr, double *berr, SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void zppsv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_Z_ *ap, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(fact));

void zppsvx_(SimTK_FOPT_(fact), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_Z_ *ap, SimTK_Z_ *afp, SimTK_CHAR_OUTPUT_(equed), double *s, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *x, SimTK_FDIM_(ldx), SimTK_D_OUTPUT_(rcond), double *ferr, double *berr, SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(uplo), SimTK_FLEN_(equed));

void zpptrf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *ap, SimTK_INFO_, SimTK_FLEN_(uplo));

void zpptri_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *ap, SimTK_INFO_, SimTK_FLEN_(uplo));

void zpptrs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_Z_ *ap, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void zptcon_(SimTK_FDIM_(n), const double *d__, const SimTK_Z_ *e, SimTK_D_INPUT_(anorm), SimTK_D_OUTPUT_(rcond), double *rwork, SimTK_INFO_);

void zptrfs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const double *d__, const SimTK_Z_ *e, const double *df, const SimTK_Z_ *ef, const SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *x, SimTK_FDIM_(ldx), double *ferr, double *berr, SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void zptsv_(SimTK_FDIM_(n), SimTK_FDIM_(nrhs), double *d__, SimTK_Z_ *e, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_INFO_);

void zptsvx_(SimTK_FOPT_(fact), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const double *d__, const SimTK_Z_ *e, double *df, SimTK_Z_ *ef, const SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *x, SimTK_FDIM_(ldx), SimTK_D_OUTPUT_(rcond), double *ferr, double *berr, SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(fact));

void zpttrf_(SimTK_FDIM_(n), double *d__, SimTK_Z_ *e, SimTK_INFO_);

void zpttrs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const double *d__, const SimTK_Z_ *e, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void zptts2_(SimTK_I_INPUT_(iuplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const double *d__, const SimTK_Z_ *e, SimTK_Z_ *b, SimTK_FDIM_(ldb));

void zrot_(SimTK_FDIM_(n), SimTK_Z_ *cx, SimTK_FINC_(x), SimTK_Z_ *cy, SimTK_FINC_(y), const double *c__, const SimTK_Z_ *s);

void zspcon_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const SimTK_Z_ *ap, const int *ipiv, SimTK_D_INPUT_(anorm), SimTK_D_OUTPUT_(rcond), SimTK_Z_ *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void zspmv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_INPUT_(alpha), const SimTK_Z_ *ap, const SimTK_Z_ *x, SimTK_FINC_(x), SimTK_Z_ *beta, SimTK_Z_ *y, SimTK_FINC_(y), SimTK_FLEN_(uplo));

void zspr_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_INPUT_(alpha), const SimTK_Z_ *x, SimTK_FINC_(x), SimTK_Z_ *ap, SimTK_FLEN_(uplo));

void zsprfs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_Z_ *ap, const SimTK_Z_ *afp, const int *ipiv, const SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *x, SimTK_FDIM_(ldx), double *ferr, double *berr, SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void zspsv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_Z_ *ap, int *ipiv, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void zspsvx_(SimTK_FOPT_(fact), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_Z_ *ap, SimTK_Z_ *afp, int *ipiv, const SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *x, SimTK_FDIM_(ldx), SimTK_D_OUTPUT_(rcond), double *ferr, double *berr, SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(uplo));

void zsptrf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *ap, int *ipiv, SimTK_INFO_, SimTK_FLEN_(uplo));

void zsptri_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *ap, const int *ipiv, SimTK_Z_ *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void zsptrs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_Z_ *ap, const int *ipiv, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void zstedc_(SimTK_FOPT_(compz), SimTK_FDIM_(n), double *d__, double *e, SimTK_Z_ *z__, SimTK_FDIM_(ldz), SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, int *lrwork, int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_, SimTK_FLEN_(compz));

void zstein_(SimTK_FDIM_(n), const double *d__, const double *e, SimTK_FDIM_(m), const double *w, int *iblock, const int *isplit, const SimTK_Z_ *z__, SimTK_FDIM_(ldz), double *work, int *iwork, int *ifail, SimTK_INFO_);

void zstemr_( SimTK_FOPT_(jobz), SimTK_FOPT_(range), SimTK_FDIM_(n), double *d, double *e, SimTK_D_INPUT_(vl), SimTK_D_INPUT_(vu), SimTK_I_INPUT_(il), SimTK_I_INPUT_(iu), SimTK_I_OUTPUT_(m), double *w, SimTK_Z_ *z, SimTK_FDIM_(ldz), SimTK_I_INPUT_(nzc), SimTK_I_OUTPUT_(isuppz), SimTK_I_OUTPUT_(tryrac), double *work, SimTK_FDIM_(lwork), int *iwork, SimTK_FDIM_(liwork), SimTK_INFO_, SimTK_FLEN_(jobz), SimTK_FLEN_(range) ); 

void zsteqr_(SimTK_FOPT_(compz), SimTK_FDIM_(n), double *d__, double *e, SimTK_Z_ *z__, SimTK_FDIM_(ldz), double *work, SimTK_INFO_, SimTK_FLEN_(compz));

void zsycon_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const SimTK_Z_ *a, SimTK_FDIM_(lda), const int *ipiv, SimTK_D_INPUT_(anorm), SimTK_D_OUTPUT_(rcond), SimTK_Z_ *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void zsymv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_INPUT_(alpha), const SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *x, SimTK_FINC_(x), SimTK_Z_INPUT_(beta), SimTK_Z_ *y, SimTK_FINC_(y), SimTK_FLEN_(uplo));

void zsyr_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_INPUT_(alpha), SimTK_Z_ *x, SimTK_FINC_(x), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_FLEN_(uplo));

void zsyrfs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *af, SimTK_FDIM_(ldaf), const int *ipiv, const SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *x, SimTK_FDIM_(ldx), double *ferr, double *berr, SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(uplo));

void zsysv_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_Z_ *a, SimTK_FDIM_(lda), int *ipiv, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(uplo));

void zsysvx_(SimTK_FOPT_(fact), SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *af, SimTK_FDIM_(ldaf), int *ipiv, const SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *x, SimTK_FDIM_(ldx), SimTK_D_OUTPUT_(rcond), double *ferr, double *berr, SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, SimTK_INFO_, SimTK_FLEN_(fact), SimTK_FLEN_(uplo));

void zsytf2_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), int *ipiv, SimTK_INFO_, SimTK_FLEN_(uplo));

void zsytrf_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), int *ipiv, SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(uplo));

void zsytri_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const SimTK_Z_ *a, SimTK_FDIM_(lda), const int *ipiv, SimTK_Z_ *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void zsytrs_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), SimTK_Z_ *a, SimTK_FDIM_(lda), int *ipiv, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo));

void ztbcon_(SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(kd), const SimTK_Z_ *ab, SimTK_FDIM_(ldab), SimTK_D_OUTPUT_(rcond), SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(norm), SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

void ztbrfs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_FDIM_(nrhs), const SimTK_Z_ *ab, SimTK_FDIM_(ldab), const SimTK_Z_ *b, SimTK_FDIM_(ldb), const SimTK_Z_ *x, SimTK_FDIM_(ldx), double *ferr, double *berr, SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));

void ztbtrs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(kd), SimTK_FDIM_(nrhs), const SimTK_Z_ *ab, SimTK_FDIM_(ldab), const SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));

void ztgevc_(SimTK_FOPT_(side), SimTK_FOPT_(howmny),  const int *select, SimTK_FDIM_(n),  const SimTK_Z_ *a, SimTK_FDIM_(lda),  const SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *vl, SimTK_FDIM_(ldvl), SimTK_Z_ *vr, SimTK_FDIM_(ldvr), SimTK_I_INPUT_(mm), SimTK_FDIM_(m), SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(howmny));

void ztgex2_(SimTK_I_INPUT_(wantq), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *q, SimTK_FDIM_(ldq), SimTK_Z_ *z__, SimTK_FDIM_(ldz), SimTK_I_INPUT_(j1), SimTK_INFO_);

void ztgexc_(SimTK_I_INPUT_(wantq), SimTK_I_INPUT_(wantz), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *q, SimTK_FDIM_(ldq), SimTK_Z_ *z__, SimTK_FDIM_(ldz), SimTK_I_OUTPUT_(ifst), SimTK_I_OUTPUT_(ilst), SimTK_INFO_);

void ztgsen_(SimTK_I_INPUT_(ijob), SimTK_I_INPUT_(wantq), SimTK_I_INPUT_(wantz), const int *select, SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *alpha, SimTK_Z_ *beta, SimTK_Z_ *q, SimTK_FDIM_(ldq), SimTK_Z_ *z__, SimTK_FDIM_(ldz), SimTK_I_OUTPUT_(m), SimTK_D_OUTPUT_(pl),  SimTK_D_OUTPUT_(pr), double *dif, SimTK_Z_ *work, SimTK_FDIM_(lwork), int *iwork, SimTK_I_INPUT_(liwork), SimTK_INFO_);

void ztgsja_(SimTK_FOPT_(jobu), SimTK_FOPT_(jobv), SimTK_FOPT_(jobq), SimTK_FDIM_(m), SimTK_FDIM_(p), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_FDIM_(l), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_D_INPUT_(tola), SimTK_D_INPUT_(tolb), double *alpha, double *beta, SimTK_Z_ *u, SimTK_FDIM_(ldu), SimTK_Z_ *v, SimTK_FDIM_(ldv), SimTK_Z_ *q, SimTK_FDIM_(ldq), SimTK_Z_ *work, SimTK_I_OUTPUT_(ncycle), SimTK_INFO_, SimTK_FLEN_(jobu), SimTK_FLEN_(jobv), SimTK_FLEN_(jobq));

void ztgsna_(SimTK_FOPT_(job), SimTK_FOPT_(howmny), const int *select, SimTK_FDIM_(n), const SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *b, SimTK_FDIM_(ldb), const SimTK_Z_ *vl, SimTK_FDIM_(ldvl), const SimTK_Z_ *vr, SimTK_FDIM_(ldvr), double *s, double *dif, SimTK_FDIM_(mm), SimTK_FDIM_(m), SimTK_Z_ *work, SimTK_FDIM_(lwork), int *iwork, SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(howmny));

void ztgsy2_(SimTK_FOPT_(trans), SimTK_I_INPUT_(ijob), SimTK_FDIM_(m), SimTK_FDIM_(n), const SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *c__, SimTK_FDIM_(ldc), const SimTK_Z_ *d__, SimTK_FDIM_(ldd), const SimTK_Z_ *e, SimTK_FDIM_(lde), SimTK_Z_ *f, SimTK_FDIM_(ldf), SimTK_D_OUTPUT_(scale), SimTK_D_OUTPUT_(rdsum), SimTK_D_OUTPUT_(rdscal), SimTK_INFO_, SimTK_FLEN_(trans));

void ztgsyl_(SimTK_FOPT_(trans), SimTK_I_INPUT_(ijob), SimTK_FDIM_(m), SimTK_FDIM_(n), const SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *c__, SimTK_FDIM_(ldc), const SimTK_Z_ *d__, SimTK_FDIM_(ldd), const SimTK_Z_ *e, SimTK_FDIM_(lde), SimTK_Z_ *f, SimTK_FDIM_(ldf), SimTK_D_OUTPUT_(scale), double *dif, SimTK_Z_ *work, SimTK_FDIM_(lwork), int *iwork, SimTK_INFO_, SimTK_FLEN_(trans));

void ztpcon_(SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), const SimTK_Z_ *ap, SimTK_D_OUTPUT_(rcond), SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(norm), SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

void ztprfs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_Z_ *ap, const SimTK_Z_ *b, SimTK_FDIM_(ldb), const SimTK_Z_ *x, SimTK_FDIM_(ldx), double *ferr, double *berr, SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));

void ztptri_(SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_Z_ *ap, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

void ztptrs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_Z_ *ap, SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo),SimTK_FLEN_(trans), SimTK_FLEN_(diag));

void ztrcon_(SimTK_FOPT_(norm), SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), const SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_D_OUTPUT_(rcond), SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(norm), SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

void ztrevc_(SimTK_FOPT_(side), SimTK_FOPT_(howmny), const int *select, SimTK_FDIM_(n), SimTK_Z_ *t, SimTK_FDIM_(ldt), SimTK_Z_ *vl, SimTK_FDIM_(ldvl), SimTK_Z_ *vr, SimTK_FDIM_(ldvr), SimTK_FDIM_(mm), SimTK_FDIM_(m), SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(howmny));

void ztrexc_(SimTK_FOPT_(compq), SimTK_FDIM_(n), SimTK_Z_ *t, SimTK_FDIM_(ldt), SimTK_Z_ *q, SimTK_FDIM_(ldq), SimTK_I_INPUT_(ifst), SimTK_I_INPUT_(ilst), SimTK_INFO_, SimTK_FLEN_(compq));

void ztrrfs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *b, SimTK_FDIM_(ldb), const SimTK_Z_ *x, SimTK_FDIM_(ldx), double *ferr, double *berr, SimTK_Z_ *work, double *rwork, SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));

void ztrsen_(SimTK_FOPT_(job), SimTK_FOPT_(compq), const int *select, SimTK_FDIM_(n), SimTK_Z_ *t, SimTK_FDIM_(ldt), SimTK_Z_ *q, SimTK_FDIM_(ldq), SimTK_Z_ *w, SimTK_FDIM_(m), SimTK_D_OUTPUT_(s), SimTK_D_OUTPUT_(sep), SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(compq));

void ztrsna_(SimTK_FOPT_(job), SimTK_FOPT_(howmny), const int *select, SimTK_FDIM_(n), SimTK_Z_ *t, SimTK_FDIM_(ldt), SimTK_Z_ *vl, SimTK_FDIM_(ldvl), SimTK_Z_ *vr, SimTK_FDIM_(ldvr), double *s, double *sep, SimTK_FDIM_(mm), SimTK_FDIM_(m), SimTK_Z_ *work, SimTK_FDIM_(ldwork), double *rwork, SimTK_INFO_, SimTK_FLEN_(job), SimTK_FLEN_(howmny));

void ztrsyl_(char *tranA, char *tranB, int *isgn, SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_Z_ *c__, SimTK_FDIM_(ldc), SimTK_D_OUTPUT_(scale), SimTK_INFO_, SimTK_FLEN_(transA), SimTK_FLEN_(transB));

void ztrti2_(SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

void ztrtri_(SimTK_FOPT_(uplo), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(diag));

void ztrtrs_(SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FOPT_(diag), SimTK_FDIM_(n), SimTK_FDIM_(nrhs), const SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *b, SimTK_FDIM_(ldb), SimTK_INFO_, SimTK_FLEN_(uplo), SimTK_FLEN_(trans), SimTK_FLEN_(diag));

void ztzrqf_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *tau, SimTK_INFO_);

void ztzrzf_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *tau, SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void zung2l_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *tau, SimTK_Z_ *work, SimTK_INFO_);

void zung2r_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *tau, SimTK_Z_ *work, SimTK_INFO_);

void zungbr_(SimTK_FOPT_(vect), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *tau, SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(vect));

void zunghr_(SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *tau, SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void zungl2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *tau, SimTK_Z_ *work, SimTK_INFO_);

void zunglq_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *tau, SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void zungql_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *tau, SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void zungqr_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *tau, SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void zungr2_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *tau, SimTK_Z_ *work, SimTK_INFO_);

void zungrq_(SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *tau, SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_);

void zungtr_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *tau, SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(uplo));

void zunm2l_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *tau, SimTK_Z_ *c__, SimTK_FDIM_(ldc), SimTK_Z_ *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void zunm2r_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *tau, SimTK_Z_ *c__, SimTK_FDIM_(ldc), SimTK_Z_ *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void zunmbr_(SimTK_FOPT_(vect), SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *tau, SimTK_Z_ *c__, SimTK_FDIM_(ldc), SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(vect), SimTK_FLEN_(side), SimTK_FLEN_(trans));

void zunmhr_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_I_INPUT_(ilo), SimTK_I_INPUT_(ihi), const SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *tau, SimTK_Z_ *c__, SimTK_FDIM_(ldc), SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void zunml2_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *tau, SimTK_Z_ *c__, SimTK_FDIM_(ldc), SimTK_Z_ *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void zunmlq_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *tau, SimTK_Z_ *c__, SimTK_FDIM_(ldc), SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void zunmql_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *tau, SimTK_Z_ *c__, SimTK_FDIM_(ldc), SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void zunmqr_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *tau, SimTK_Z_ *c__, SimTK_FDIM_(ldc), SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void zunmr2_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *tau, SimTK_Z_ *c__, SimTK_FDIM_(ldc), SimTK_Z_ *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void zunmr3_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_I_INPUT_(l), const SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *tau, SimTK_Z_ *c__, SimTK_FDIM_(ldc), SimTK_Z_ *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(trans));

void zunmrq_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), const SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *tau, SimTK_Z_ *c__, SimTK_FDIM_(ldc), SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(uplo));

void zunmrz_(SimTK_FOPT_(side), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), SimTK_FDIM_(k), SimTK_FDIM_(l), const SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *tau, SimTK_Z_ *c__, SimTK_FDIM_(ldc), SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(uplo));

void zunmtr_(SimTK_FOPT_(side), SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), const SimTK_Z_ *a, SimTK_FDIM_(lda), const SimTK_Z_ *tau, SimTK_Z_ *c__, SimTK_FDIM_(ldc), SimTK_Z_ *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(uplo), SimTK_FLEN_(trans));

void zupgtr_(SimTK_FOPT_(uplo), SimTK_FDIM_(n), const SimTK_Z_ *ap, const SimTK_Z_ *tau, SimTK_Z_ *q, SimTK_FDIM_(ldq), SimTK_Z_ *work, SimTK_INFO_, SimTK_FLEN_(uplo));

void zupmtr_(SimTK_FOPT_(side), SimTK_FOPT_(uplo), SimTK_FOPT_(trans), SimTK_FDIM_(m), SimTK_FDIM_(n), const SimTK_Z_ *ap, const SimTK_Z_ *tau, SimTK_Z_ *c__, SimTK_FDIM_(ldc), SimTK_Z_ *work, SimTK_INFO_, SimTK_FLEN_(side), SimTK_FLEN_(uplo), SimTK_FLEN_(trans));

#ifdef __cplusplus
}   /* "C" */
#endif

#undef SimTK_C_
#undef SimTK_Z_
#undef SimTK_FDIM_
#undef SimTK_FOPT_
#undef SimTK_FLEN_
#undef SimTK_FINC_
#undef SimTK_INFO_

#endif /* SimTK_SIMTKLAPACK_H_ */

