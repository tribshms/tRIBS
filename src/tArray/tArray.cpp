/***************************************************************************
**
**  		      tRIBS Distributed Hydrologic Model  
**
**              TIN-based Real-time Integrated Basin Simulator
**		         Ralph M. Parsons Laboratory
**  		    Massachusetts Institute of Technology
**
**
**  tArray.cpp: Functions for templated class tArray< T >
**
***************************************************************************/

#include "src/tArray/tArray.h"

//=========================================================================
//
//
//                  Section 1: tArray Constructors/Destructors
//
//
//=========================================================================

//Default constructor

template< class T > 
tArray< T >::
tArray() :
npts(1), avalue(0)
{
	avalue = new T [1];
	assert( avalue != 0 );
	avalue[0] = 0;
}

//Copy constructor

template< class T > 
tArray< T >::
tArray( const tArray< T > &original ) :
npts(original.npts), avalue(0)
{
	int i;
	assert( npts > 0 );
        avalue = new T[npts];
        assert( avalue != 0 );
        for( i = 0; i < npts; i++ )
            avalue[i] = original.avalue[i];
}

//=========================================================================
//
//
//                  Section 2: tArray Overloaded Operators
//
//
//=========================================================================

/**************************************************************************
**  Overloaded operators:
**
**    assignment, equality, inequality - memberwise operation
**    index - uses an assertion to check array bounds (assumed to be within
**           bounds at runtime) and returns value
**    left shift - sends array values to output stream, separated by
**                 spaces and with a carriage return after every 10 items
**
**************************************************************************/

//Overloaded assignment operator

template< class T >                                              
const tArray< T > &tArray< T >::operator=( const tArray< T > &right )
{
	assert( &right != 0 );
	int i;
	if( &right != this ){
		delete [] avalue; avalue = 0;
		npts = right.npts;
		if( npts>0 ){
			assert( right.avalue != 0 );
			avalue = new T [npts];
			for( i = 0; i < npts; i++ ){
				avalue[i] = right.avalue[i];
			}
		}
	}
	return *this;
}

//Overloaded equality operator:

template< class T >                                              
int tArray< T >::operator==( const tArray< T > &right ) const
{
	if( npts != right.npts ) return 0;
	int i;
	for( i = 0; i < npts; i++ )
		if( avalue[i] != right.avalue[i] )
			return 0;
	return 1;
}

//overloaded inequality operator:

template< class T >                                               
int tArray< T >::operator!=( const tArray< T > &right ) const
{
	if( npts != right.npts ) return 0;
	int i;
	for( i = 0; i < npts; i++ )
		if( avalue[i] != right.avalue[i] )
			return 1;
	return 0;
}

//checkSubscript function
template< class T >                                               //tArray
void tArray< T >::checkSubscript( int subscript ) const
{
	if( 0 > subscript || subscript >= npts ){
		cout<<"subscript is "<<subscript<<" npts is "<<npts<<endl<<flush;
		cout<<"bailing out of tArray[]"<<endl<<flush;
	}
	assert( 0 <= subscript && subscript < npts );
}

//Overloaded subscript operator 2
template< class T >                                               //tArray
const T &tArray< T >::operator[]( int subscript ) const
{
	checkSubscript(subscript);
	return avalue[subscript];
}

//Overloaded left shift operator

template< class T >                                           
ostream &operator<<( ostream &output, const tArray< T > &a )
{
	int i;
	
	for( i = 0; i < a.npts; i++ ){
		output << a.avalue[i] << " ";
		if( (i + 1) %10 == 0 ) output<<endl;
	}
	if( i % 10 != 0 ) output<<endl;
	return output;
}


//=========================================================================
//
//
//                  Section 3: tArray Get/Set Functions
//
//
//=========================================================================

//Returns a pointer to the head of the array 

template< class T > T *tArray< T >::getArrayPtr() {
	return avalue;
}

//Reinitializes and resized array

template< class T >void tArray<T>::setSize( int size )
{
	int i;  
	delete [] avalue;
	npts = size;
	avalue = new T [npts];
	assert( avalue!=0 && npts>=0 );
	for( i=0; i<npts; i++ ) avalue[i] = 0;
}

//Adds value to Array member

template< class T > void tArray< T >::addValueTo( int subscript, T value){
	avalue[subscript] = avalue[subscript] + value;
}


//=========================================================================
//
//
//                         End of tArray.cpp
//
//
//=========================================================================
