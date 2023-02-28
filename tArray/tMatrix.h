/***************************************************************************
**
**  		      tRIBS Distributed Hydrologic Model
**
**              TIN-based Real-time Integrated Basin Simulator
**		         Ralph M. Parsons Laboratory
**  		    Massachusetts Institute of Technology
**  
**  tMatrix.h: Header file for class tMatrix.
**
**  Class tMatrix implements matrices (2D arrays) as 1D arrays of tArray
**  objects. 
**
***************************************************************************/

#ifndef TMATRIX_H
#define TMATRIX_H

#include "tArray/tArray.h"

//=========================================================================
//
//
//                  Section 1: tMatrix Class Definition
//
//
//=========================================================================

template < class T >
class tMatrix
{
public:
    tMatrix();
    tMatrix( int nr, int nc );
    ~tMatrix();
    T &operator()( int row, int col );
   int getNumRows() {return nrows;}
   int getNumCols() {return ncols;}
    
private:
    tArray<T> *ptr;
    int nrows;
    int ncols;  
};

#endif

//=========================================================================
//
//
//                           End of tMatrix.h
//
//
//=========================================================================
