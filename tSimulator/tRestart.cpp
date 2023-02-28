/***************************************************************************
**
**  		     tRIBS Distributed Hydrologic Model
**
**              TIN-based Real-time Integrated Basin Simulator
**		         Ralph M. Parsons Laboratory
**  		    Massachusetts Institute of Technology
**           and Los Alamos National Laboratory
**
**
**  tRestart.cpp: Functions for class tRestart (see tRestart.h)
**
***************************************************************************/

#include <cassert>

#include "tSimulator/tRestart.h"

/*************************************************************************
**
** Constructor
**
*************************************************************************/

template< class tSubNode >
tRestart<tSubNode>::tRestart(
	tRunTimer* t,
	tMesh<tSubNode>* m,
	tKinemat* f,
	tWaterBalance* b,
	tHydroModel* h,
	tRainfall* r,
	tEvapoTrans* e,
	tIntercept* i,
   tSnowPack* s,
   tSnowIntercept* c)
{
  timer = t;
  mesh = m;
  flow = f;
  balance = b;
  hydro = h;
  rainfall = r;
  evap = e;
  intercept = i;
  snowpack = s;
  snowintercept = c;
}

/*************************************************************************
**
** Write restart information for all controlled objects
**
*************************************************************************/

template< class tSubNode >
void tRestart<tSubNode>::writeRestart(fstream & rStr)
{
  timer->writeRestart(rStr);
  flow->writeRestart(rStr);
  balance->writeRestart(rStr);
  hydro->writeRestart(rStr);
  rainfall->writeRestart(rStr);
  evap->writeRestart(rStr);
  intercept->writeRestart(rStr);
  mesh->writeRestart(rStr);
	// Giuseppe DEBUG Restart 2012 - START 
	// I have introduced an IF that checks whether
	// if the snow module is on. If not, the relative variables 
	// are not saved in the binary Restart files.
	if (snowpack->getSnowOpt() != 0){
		snowpack->writeRestart(rStr);
		snowintercept->writeRestart(rStr); // saving
	}// Giuseppe DEBUG Restart 2012 - END
		

}

/*************************************************************************
**
** Read restart information for all controlled objects
**
*************************************************************************/

template< class tSubNode >
void tRestart<tSubNode>::readRestart(fstream & rStr)
{
	timer->readRestart(rStr);
	flow->readRestart(rStr);
	balance->readRestart(rStr);
	hydro->readRestart(rStr);
	rainfall->readRestart(rStr);
	evap->readRestart(rStr);
	intercept->readRestart(rStr);
	mesh->readRestart(rStr);
	// Giuseppe DEBUG Restart 2012 - START 
	// I have introduced an IF that checks whether
	// if the snow module is on. This is to be consistent
	// with the writeRestart function.
	if (snowpack->getSnowOpt() != 0){	
		snowpack->readRestart(rStr);
		snowintercept->readRestart(rStr);
	}// Giuseppe DEBUG Restart 2012 - END
}

//=========================================================================
//
//
//                        End of tRestart.cpp
//
//
//=========================================================================
