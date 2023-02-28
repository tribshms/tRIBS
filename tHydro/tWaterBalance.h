/***************************************************************************
**
**                   tRIBS Distributed Hydrology Model
**
**              TIN-based Real-time Integrated Basin Simulator
**                       Ralph M. Parsons Laboratory
**                  Massachusetts Institute of Technology
**  
**
**  tWaterBalance.h:   Header file for tWaterBalance Class
**
**  Storage/computation class for the components of the water balance
**  in each node and for the entire watershed.
**
***************************************************************************/

#ifndef TWATERBALANCE_H 
#define TWATERBALANCE_H

//=========================================================================
//
//
//                  Section 1: tWaterBalance Include Statements
//
//
//=========================================================================

#include "Headers/Inclusions.h"

//=========================================================================
//
//
//                  Section 2: tWaterBalance Class Definitions
//
//
//=========================================================================

class tWaterBalance
{
public:
  tWaterBalance();
  tWaterBalance(SimulationControl*,tMesh<tCNode> *,tInputFile &);
  ~tWaterBalance();

  void initializeVariables();
  void SetWaterBalance(tInputFile &);
  void DeleteWaterBalance();
  void CanopyBalance();
  void UnSaturatedBalance();
  void SaturatedBalance();
  void BasinStorage(double);
  void Print(double *);
  void writeRestart(fstream &) const;
  void readRestart(fstream &);

protected:
  tMesh<tCNode> *gridPtr;      
  SimulationControl *simCtrl; 

  int metStep, unsStep, satStep, finalTime;
  double *BasinStorages;
};

#endif

//=========================================================================
//
//
//                       End of tWaterBalance.h
//
//
//=========================================================================
