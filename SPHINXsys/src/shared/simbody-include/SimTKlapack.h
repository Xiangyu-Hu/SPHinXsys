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
extern SimTK_C_ cdotu_(SimTK_FDIM_(n), const SimTK_C_ *x, SimTK_FINC_(x), 
                                                const SimTK_C_ *y, SimTK_FINC_(y));
extern SimTK_Z_ zdotu_(SimTK_FDIM_(n), const SimTK_Z_ *x, SimTK_FINC_(x), 
                                                const SimTK_Z_ *y, SimTK_FINC_(y));
extern SimTK_C_ cdotc_(SimTK_FDIM_(n), const SimTK_C_ *x, SimTK_FINC_(x), 
                                                const SimTK_C_ *y, SimTK_FINC_(y));
extern SimTK_Z_ zdotc_(SimTK_FDIM_(n), const SimTK_Z_ *x, SimTK_FINC_(x), 
                                                const SimTK_Z_ *y, SimTK_FINC_(y));
#else
extern void cdotu_(SimTK_C_*, SimTK_FDIM_(n), const SimTK_C_ *x, SimTK_FINC_(x), 
                                                const SimTK_C_ *y, SimTK_FINC_(y));
extern void zdotu_(SimTK_Z_*, SimTK_FDIM_(n), const SimTK_Z_ *x, SimTK_FINC_(x), 
                                                const SimTK_Z_ *y, SimTK_FINC_(y));
extern void cdotc_(SimTK_C_*, SimTK_FDIM_(n), const SimTK_C_ *x, SimTK_FINC_(x), 
                                                const SimTK_C_ *y, SimTK_FINC_(y));
extern void zdotc_(SimTK_Z_*, SimTK_FDIM_(n), const SimTK_Z_ *x, SimTK_FINC_(x), 
                                                const SimTK_Z_ *y, SimTK_FINC_(y));
#endif
*/

extern void cgeev_(SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *w, SimTK_C_ *vl, SimTK_FDIM_(ldvl), SimTK_C_ *vr, SimTK_FDIM_(ldvr), SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, SimTK_INFO_, SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr));

extern void cgeevx_(SimTK_FOPT_(balanc), SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), SimTK_FOPT_(sense), SimTK_FDIM_(n), SimTK_C_ *a, SimTK_FDIM_(lda), SimTK_C_ *w, SimTK_C_ *vl, SimTK_FDIM_(ldvl), SimTK_C_ *vr, SimTK_FDIM_(ldvr), SimTK_I_OUTPUT_(ilo), SimTK_I_OUTPUT_(ihi), SimTK_S_OUTPUT_(scale), float *abnrm, SimTK_S_OUTPUT_(rconde), SimTK_S_OUTPUT_(rcondv), SimTK_C_ *work, SimTK_FDIM_(lwork), float *rwork, SimTK_INFO_, SimTK_FLEN_(balanc), SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr), SimTK_FLEN_(sense));

extern void dgeev_(SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *wr, double *wi, double *vl, SimTK_FDIM_(ldvl), double *vr, SimTK_FDIM_(ldvr), double *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr));

extern void dgeevx_(SimTK_FOPT_(balanc), SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), char *sense, SimTK_FDIM_(n), double *a, SimTK_FDIM_(lda), double *wr, double *wi, double *vl, SimTK_FDIM_(ldvl), double *vr, SimTK_FDIM_(ldvr), SimTK_I_OUTPUT_(ilo), SimTK_I_OUTPUT_(ihi), SimTK_D_OUTPUT_(scale), double *abnrm, SimTK_D_OUTPUT_(rconde), SimTK_D_OUTPUT_(rcondv), double *work, SimTK_FDIM_(lwork), int *iwork, SimTK_INFO_, SimTK_FLEN_(balanc), SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr), SimTK_FLEN_(sense));

extern void sgeev_(SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *wr, float *wi, float *vl, SimTK_FDIM_(ldvl), float *vr, SimTK_FDIM_(ldvr), float *work, SimTK_FDIM_(lwork), SimTK_INFO_, SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr));

extern void sgeevx_(SimTK_FOPT_(balanc), SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), char *sense, SimTK_FDIM_(n), float *a, SimTK_FDIM_(lda), float *wr, float *wi, float *vl, SimTK_FDIM_(ldvl), float *vr, SimTK_FDIM_(ldvr), SimTK_I_OUTPUT_(ilo), SimTK_I_OUTPUT_(ihi), SimTK_S_OUTPUT_(scale), float *abnrm, SimTK_S_OUTPUT_(rconde), SimTK_S_OUTPUT_(rcondv), float *work, SimTK_FDIM_(lwork), int *iwork, SimTK_INFO_, SimTK_FLEN_(balanc), SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr), SimTK_FLEN_(alpahr));

extern void zgeev_(SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *w, SimTK_Z_ *vl, SimTK_FDIM_(ldvl), SimTK_Z_ *vr, SimTK_FDIM_(ldvr), SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, SimTK_INFO_, SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr));

extern void zgeevx_(SimTK_FOPT_(balanc), SimTK_FOPT_(jobvl), SimTK_FOPT_(jobvr), char *sense, SimTK_FDIM_(n), SimTK_Z_ *a, SimTK_FDIM_(lda), SimTK_Z_ *w, SimTK_Z_ *vl, SimTK_FDIM_(ldvl), SimTK_Z_ *vr, SimTK_FDIM_(ldvr), SimTK_I_OUTPUT_(ilo), SimTK_I_OUTPUT_(ihi), SimTK_D_OUTPUT_(scale), double *abnrm, SimTK_D_OUTPUT_(rconde), SimTK_D_OUTPUT_(rcondv), SimTK_Z_ *work, SimTK_FDIM_(lwork), double *rwork, SimTK_INFO_, SimTK_FLEN_(balanc), SimTK_FLEN_(jobvl), SimTK_FLEN_(jobvr), SimTK_FLEN_(sense));


#ifdef __cplusplus
}   /* extern "C" */
#endif

#undef SimTK_C_
#undef SimTK_Z_
#undef SimTK_FDIM_
#undef SimTK_FOPT_
#undef SimTK_FLEN_
#undef SimTK_FINC_
#undef SimTK_INFO_

#endif /* SimTK_SIMTKLAPACK_H_ */

