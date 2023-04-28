/***************************************************************************
**
**  		      tRIBS Distributed Hydrologic Model  
**
**              TIN-based Real-time Integrated Basin Simulator
**		         Ralph M. Parsons Laboratory
**  		    Massachusetts Institute of Technology
**  
**  tMatrix.cpp: Functions for tMatrix Class (See tMatrix.h)
**
***************************************************************************/

#include "src/tArray/tMatrix.h"

//=========================================================================
//
//
//                  Section 1: tMatrix Constructors/Destructors
//
//
//=========================================================================

// Default constructor: initializes an empty matrix

template < class T >tMatrix<T>::tMatrix(){
	ptr = 0;
	nrows = 0;
	ncols = 0;
}

// Constructor: Sets the size to nr by nc and sets all values to zero.

template < class T > tMatrix<T>::tMatrix( int nr, int nc ){
	int i;
	ptr = new tArray<T> [nr];
	for( i=0; i<nr; i++ )
		ptr[i].setSize(nc);
}

// Destructor: deletes the array of tArray objects 

template < class T >
tMatrix<T>::~tMatrix(){
	delete [] ptr;
}

//=========================================================================
//
//
//                  Section 2: tMatrix Overloaded Operators
//
//
//=========================================================================

// Overloaded () operator for referencing individual matrix entries

template < class T > T &tMatrix<T>::operator()( int row, int col ){
	return (ptr[row])[col];  
}

//=========================================================================
//
//
//                   	    End of tMatrix.cpp
//
//
//=========================================================================


