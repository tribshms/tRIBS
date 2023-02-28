/***************************************************************************
**
**  		     tRIBS Distributed Hydrologic Model
**
**              TIN-based Real-time Integrated Basin Simulator
**		         Ralph M. Parsons Laboratory
**  		    Massachusetts Institute of Technology
**  
**
**  tResData.h:   Header file for class tResData
**
**  Functions provided to create and access Reservoir data.
**  Data members include elevation, discharge, storage and Node parameters.
**
***************************************************************************/

#ifndef RESDATA_H
#define RESDATA_H

#include "Headers/Inclusions.h"

using namespace std;
//=========================================================================
//
//
//                  Section 1: tResData Class Declaration
//
//
//=========================================================================

class tResData
{
 public:
  tResData();
  ~tResData();

  int getnumLines(char *);

  // Set/Get Functions for Reservoir properties:
  void setResType(int);
  void setResElev(double);
  void setResDischarge(double);
  void setResStorage(double);
  void setResEDS(double, int);
  void setResLines(int);
  void setRNum(int);

  int    getResType(int);
  double getResElev(int);
  double getResDischarge(int);
  double getResStorage(int);
  double getResEDS(int);
  int    getResLines();
  int    getRNum();

  // Set/Get Functions for Reservoir Polygon ID:
  void setResNodeID(int);
  void setResNodeType(int);
  void setInitial_H(double);
  void setResArraySize(int);

  int    getResNodeID();
  int    getResNodeType();
  double getInitial_H();

  int getRoutingStep();

  void setInflow(double);
  double getInflow(int);

  void setSTQnext(double, int);
  double getSTQnext(int);

 protected:

  int numLines;

  int setNum;
  int ResLines;
  int ResNodeID;
  int ResIDtype;
  double ResInElev;
  char *fileName;
  int rStep;

  int *ResType;
  double *rElev;
  double *rDischarge;
  double *rStorage;
  double *rInflow;
  double *rEDS;
  double *rSTQnext;
};

#endif

//=========================================================================
//
//
//                        End of tResData.h
//
//
//=========================================================================
