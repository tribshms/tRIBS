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
**  tRestart.h: Header for tRestart class and objects
**
**  tRestart Class used in tRIBS for saving current state of a simulation
**  with the purpose of restarting at a later time
** 
***************************************************************************/

//=========================================================================
//
//
//                  Section 1: tRestart Include and Define Statements
//
//
//=========================================================================

#ifndef TRESTART_H
#define TRESTART_H

#include <iostream>
#include "tSimulator/tRunTimer.h"
#include "tFlowNet/tKinemat.h"
#include "tFlowNet/tReservoir.h" // JECR2015
#include "tFlowNet/tResData.h" // JECR2015
#include "tHydro/tSnowPack.h"
#include "tHydro/tSnowIntercept.h"

//=========================================================================
//
//
//                  Section 2: tRestart Class Definitions
//
//
//=========================================================================

template< class tSubNode >
class tRestart {

public:
  /// Constructor
  tRestart(
        tRunTimer* t,
        tMesh<tSubNode>* m,
        tKinemat* f,
        tWaterBalance* b,
        tHydroModel* h,
        tRainfall* r,
        tEvapoTrans* e,
        tIntercept* i,
        tSnowPack* s,
        tSnowIntercept* c);

  /// Destructor
  ~tRestart() {}

  /// Write a restart dump
  void writeRestart(fstream &);
  /// Read a restart dump
  void readRestart(fstream &);

private:
  tRunTimer*         timer;           //!< Run timer
  tMesh<tSubNode>*   mesh;            //!< Mesh
  tKinemat*          flow;            //!< Kinematic flow
  tWaterBalance*     balance;         //!< Water balance
  tHydroModel*       hydro;           //!< Hydro model
  tRainfall*         rainfall;        //!< Rainfall
  tEvapoTrans*       evap;            //!< Evapotranspiration
  tIntercept*        intercept;       //!< Intercept structure
  tSnowPack*         snowpack;        //!< Snow pack structure
  tSnowIntercept*    snowintercept;   //!< Snow intercept structure
};

#endif

//=========================================================================
//
//
//                          End of tRestart.h 
//
//
//=========================================================================
