/***************************************************************************
**
**  		     tRIBS Distributed Hydrologic Model
**
**              TIN-based Real-time Integrated Basin Simulator
**		         Ralph M. Parsons Laboratory
**  		    Massachusetts Institute of Technology
**  
**
**  tInvariant.h:   Header File for tInvariant classes
**
**  tInvariant class is used for soil and land use attributes by defining
**  two Generic_Data and Data_Type classes.
**
***************************************************************************/

#ifndef  TINVARIANT_H
#define  TINVARIANT_H

#include "Headers/Inclusions.h"

class SoilType;
class LandType;   

//=========================================================================
//
//
//                  Section 1: tInvariant Generic Class Declarations
//
//
//=========================================================================

class GenericSoilData 
{
public:
  GenericSoilData(tMesh<tCNode> *, tInputFile *, tResample *);
  ~GenericSoilData();
  SoilType **SoilClass; 

  int numClass;    			// Number of soil classes
  int currClass;   			// Indicates current class 
  char soilGrid[kMaxNameSize];  	// Pathname to soil map
  char soilTable[kMaxNameSize]; 	// Pathname to soil reference table

  void setSoilPtr(int soilID);
  void printSoilPars();
  void SetSoilParameters(tMesh<tCNode>*, tResample*, tInputFile&, int); 
  double getSoilProp(int);
};

class GenericLandData 
{      
public:
  GenericLandData(tMesh<tCNode> *, tInputFile *, tResample *);
  ~GenericLandData();

  LandType **LandClass;   
  int numClass;      
  int currClass;     
  char landGrid[kMaxNameSize];  
  char landTable[kMaxNameSize]; 

  void setLandPtr(int landID);
  void SetLtypeParameters(tMesh<tCNode>*, tResample*, tInputFile&, int); 
  double getLandProp(int);
};


//=========================================================================
//
//
//                  Section 2: tInvariant Type Class Declarations
//
//
//=========================================================================

class SoilType 
{
public:
  SoilType();
  SoilType(double *, int);
  ~SoilType();

  int numProps;       	// Number of soil properties to store
  double *sProperty;    // Array containing soil properties

  void   setProperty(int, double);
  double getProperty(int);
};


class LandType 
{         
public:     
  LandType();
  LandType(double *, int);
  ~LandType();

  int numProps;         // Number of land properties to store
  double *lProperty;    // Array containing land properties
  
  void   setProperty(int, double);
  double getProperty(int);
};

#endif

//=========================================================================
//
//
//                         End of tInvariant.h 
//
//
//=========================================================================
