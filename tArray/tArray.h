/***************************************************************************
**
**  		      tRIBS Distributed Hydrologic Model 
**
**              TIN-based Real-time Integrated Basin Simulator
**		         Ralph M. Parsons Laboratory
**  		    Massachusetts Institute of Technology
**
**
**  tArray.h: Header file for tArray Class objects
**
**  A tArray is an object implementation of a one-dimensional array.
**  A template is used so that the array may be of any specified type.
**  Unlike regular arrays, tArray objects provide bounds checking,
**  memberwise equality/inequality comparison, and memberwise copy
**  operations. The size of the array is determined either by an 
**  argument passed to the constructor or by assignment of one array
**  to another. 
**
***************************************************************************/

#ifndef TARRAY_H
#define TARRAY_H

#include "Headers/tribs_os.h"

#ifdef ALPHA_64
  #include <iostream.h>
  #include <fstream.h>
  #include <assert.h>
#elif defined LINUX_32
  #include <iostream>
  #include <fstream>
  #include <cassert>
#elif defined WIN
  #include <iostream.h>
  #include <fstream.h>
  #include <assert.h>
#else 
  #include <iostream.h>
  #include <fstream.h>
  #include <assert.h>
#endif

using namespace std;

//=========================================================================
//
//
//                  Section 1: tArray Class Definition
//
//
//=========================================================================

template< class T > 
class tArray
{
public:
    friend class tFlowNet;
    tArray();                      	// default constructor
    tArray( int number ) :
      npts(number), avalue(0) // constructor initialize array size
    {
      assert( number > 0 );
      avalue = new T [npts];
      assert( avalue != 0 );
      for( int i=0; i<npts; i++ )
         avalue[i] = 0;
    }

    tArray( const tArray< T > & ); 	// copy constructor
    ~tArray()                     	// destructor
    { delete [] avalue;}

    const tArray< T > &operator=( const tArray< T > & ); // memberwise assignmt
    int operator==( const tArray< T > & ) const;    	// memberwise comparison
    int operator!=( const tArray< T > & ) const;     
    T &operator[]( int subscript )      // overloaded array index operator 
    {
      checkSubscript(subscript);
      return avalue[subscript];
    }
  
    const T &operator[]( int ) const;


    int getSize() const        // returns the number of elements in the array
    { return npts; } 

    void setSize( int );       // reinitializes and sets array size
    T *getArrayPtr();          
    void addValueTo( int , T );
    void checkSubscript(int) const;

protected:
    int npts;   		// size of array
    T * avalue; 		// an array
};

#endif

//=========================================================================
//
//
//                           End of tArray.h
//
//
//=========================================================================
