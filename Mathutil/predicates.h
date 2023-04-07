/***************************************************************************
**
**  		             tRIBS Version 1.0
**
**              TIN-based Real-time Integrated Basin Simulator
**		         Ralph M. Parsons Laboratory
**  		    Massachusetts Institute of Technology
**  
**
**		            Beta Release, 9/2001
**
**
**  predicates.h:   Header File for predicates.cpp
**
**  Functions called from member functions of tMesh to check for line
**  segment intersection since inexact arithmetic can lead to erroneous
**  answers. (see predicates.txt for more information)
**
***************************************************************************/

#ifndef PREDICATES_H
#define PREDICATES_H

//=========================================================================
//
//
//                  Section 1: predicates Include and Define Statements
//
//
//=========================================================================


#ifdef ALPHA_64
  #include <stdio.h> 
  #include <math.h> 
  #include <sys/time.h> 
#elif defined LINUX_32
  #include <cstdio> 
  #include <cmath> 
  #include <ctime> 
#elif defined WIN
  #include <stdio.h> 
  #include <math.h> 
  #include <sys/time.h> 
#else 
  #include <stdio.h> 
  #include <math.h> 
  #include <sys/time.h> 
#endif
 
#define INEXACT                         
#define REAL double                      // float or double 
 
// #define Absolute(a)  abs(a)              // defined in <macros.h>
#define Absolute(a)  ((a) >= 0.0 ? (a) : -(a)) 
 
#define Fast_Two_Sum_Tail(a, b, x, y) bvirt = x - a; y = b - bvirt 
 
#define Fast_Two_Sum(a, b, x, y) x = (REAL) (a + b); Fast_Two_Sum_Tail(a, b, x, y) 

#define Two_Sum_Tail(a, b, x, y) bvirt = (REAL) (x - a); avirt = x - bvirt; bround = b - bvirt; around = a - avirt; y = around + bround 
 
#define Two_Sum(a, b, x, y) x = (REAL) (a + b); Two_Sum_Tail(a, b, x, y) 
 
#define Two_Diff_Tail(a, b, x, y) bvirt = (REAL) (a - x); avirt = x + bvirt; bround = bvirt - b; around = a - avirt; y = around + bround 
 
#define Two_Diff(a, b, x, y) x = (REAL) (a - b); Two_Diff_Tail(a, b, x, y) 
 
#define Split(a, ahi, alo) c = (REAL) (splitter * a); abig = (REAL) (c - a); ahi = c - abig; alo = a - ahi 
 
#define Two_Product_Tail(a, b, x, y) Split(a, ahi, alo); Split(b, bhi, blo); err1 = x - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); y = (alo * blo) - err3 
 
#define Two_Product(a, b, x, y) x = (REAL) (a * b); Two_Product_Tail(a, b, x, y) 
 
#define Two_Product_Presplit(a, b, bhi, blo, x, y) x = (REAL) (a * b); Split(a, ahi, alo); err1 = x - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); y = (alo * blo) - err3 
 
#define Square_Tail(a, x, y) Split(a, ahi, alo); err1 = x - (ahi * ahi); err3 = err1 - ((ahi + ahi) * alo); y = (alo * alo) - err3 
 
#define Square(a, x, y) x = (REAL) (a * a); Square_Tail(a, x, y) 
 
#define Two_One_Sum(a1, a0, b, x2, x1, x0) Two_Sum(a0, b , _i, x0); Two_Sum(a1, _i, x2, x1) 
 
#define Two_One_Diff(a1, a0, b, x2, x1, x0) Two_Diff(a0, b , _i, x0); Two_Sum( a1, _i, x2, x1) 
 
#define Two_Two_Sum(a1, a0, b1, b0, x3, x2, x1, x0) Two_One_Sum(a1, a0, b0, _j, _0, x0); Two_One_Sum(_j, _0, b1, x3, x2, x1) 
 
#define Two_Two_Diff(a1, a0, b1, b0, x3, x2, x1, x0) Two_One_Diff(a1, a0, b0, _j, _0, x0); Two_One_Diff(_j, _0, b1, x3, x2, x1) 
 

//=========================================================================
//
//
//                  Section 2: predicates Class Declaration
//
//
//=========================================================================

class Predicates
{
public:
   Predicates(); // just calls exactinit()
   ~Predicates() {} // doesn't do anything
   
   // Constructor:
   void exactinit();
   
   // Basic "exact" arithmetic:
   int grow_expansion(int elen, REAL* e, REAL b, REAL* h);
   int grow_expansion_zeroelim(int elen, REAL* e, REAL b, REAL* h);
   int expansion_sum(int elen, REAL* e, int flen, REAL* f, REAL* h);
   int expansion_sum_zeroelim1(int elen, REAL* e, int flen, REAL* f,
                               REAL* h);
   int expansion_sum_zeroelim2(int elen, REAL* e, int flen, REAL* f,
                               REAL* h);
   int fast_expansion_sum(int elen, REAL* e, int flen, REAL* f, REAL* h);
   int fast_expansion_sum_zeroelim(int elen, REAL* e, int flen, REAL* f,
                                   REAL* h);
   int linear_expansion_sum(int elen, REAL* e, int flen, REAL* f, REAL* h);
   int linear_expansion_sum_zeroelim(int elen, REAL* e, int flen, REAL* f,
                                     REAL* h);
   int scale_expansion(int elen, REAL* e, REAL b, REAL* h);
   int scale_expansion_zeroelim(int elen, REAL* e, REAL b, REAL* h);
   int compress(int elen, REAL* e, REAL* h);
   REAL estimate( int elen, REAL* e );
   
   // Lancaster functions
   double DifferenceOfProductsOfDifferences( double a, double b,
                                             double c, double d,
                                             double e, double f,
                                             double g, double h );
   double AdaptDiffOfProdsOfDiffs( double a, double b,
                                   double c, double d,
                                   double e, double f,
                                   double g, double h,
                                   double sum );

   REAL orient2dfast(REAL *pa, REAL *pb, REAL *pc);
   REAL orient2dadapt(REAL *pa, REAL *pb, REAL *pc, REAL detsum);
   REAL orient2d(REAL *pa, REAL *pb, REAL *pc);
   REAL incirclefast(REAL *pa, REAL *pb, REAL *pc, REAL *pd);
   REAL incircleadapt(REAL *pa, REAL *pb, REAL *pc, REAL *pd,
                      REAL permanent);
   REAL incircle(REAL *pa, REAL *pb, REAL *pc, REAL *pd);
   
private:
   REAL splitter;     /* = 2^ceiling(p / 2) + 1.  Used to split floats in half. */ 
   REAL epsilon;                /* = 2^(-p).  Used to estimate roundoff errors. */ 
   REAL resulterrbound; 
   REAL ccwerrboundA, ccwerrboundB, ccwerrboundC; 
   REAL o3derrboundA, o3derrboundB, o3derrboundC; 
   REAL iccerrboundA, iccerrboundB, iccerrboundC; 
   REAL isperrboundA, isperrboundB, isperrboundC; 
};

#endif

//=========================================================================
//
//
//                           End of predicates.h 
//
//
//=========================================================================
