/***************************************************************************
**
**  		     tRIBS Distributed Hydrologic Model
**
**              TIN-based Real-time Integrated Basin Simulator
**		         Ralph M. Parsons Laboratory
**  		    Massachusetts Institute of Technology
**  
**
**  tPreProcess.h:  Header file for tPreProcess Class
**
**  This class is utilized for data pre-processing purposes which include
**  meteorological and rain gauge extraction from available sources. The
**  preprocessing class should be expandable to handle future data sources.
**
***************************************************************************/

#ifndef PREPROCESS_H
#define PREPROCESS_H

#include "Headers/tribs_os.h"
#include "Headers/Inclusions.h"

#ifdef ALPHA_64
#elif defined LINUX_32
#elif defined MAC
  #ifndef _WIN32
  #include <unistd.h>
  #endif
#elif defined WIN
#else 
#endif

//=========================================================================
//
//
//            Section 1: tPreProcess Class Declaration
//
//
//=========================================================================

class tPreProcess
{
 public:
  tPreProcess();
  tPreProcess(SimulationControl*, tInputFile &);
  ~tPreProcess();
  
  SimulationControl *simCtrl;

  void CheckInputFile(tInputFile &);
  void CheckFileExists(tInputFile &, char *, const char*);
  void CheckPathNameCorrect(tInputFile &, char*, const char*);

  int    IterReadItem(tInputFile &,    int, const char *);
  void   IterReadItem(tInputFile &, char *, const char *);
  double IterReadItem(tInputFile &, double, const char *);

  int convertData;

};

#endif

//=========================================================================
//
//
//                    End of tPreProcess.h
//
//
//=========================================================================
