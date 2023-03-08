/***************************************************************************
**
**  		     tRIBS Distributed Hydrologic Model
**
**              TIN-based Real-time Integrated Basin Simulator
**		         Ralph M. Parsons Laboratory
**  		    Massachusetts Institute of Technology
**  
**
**  tInputFile.h: Header for tInputFile class and objects
**
**  tInputFile Class used in tRIBS for inputing parameters and pathnames
**  using keywords in a *.in file read within various classes.
** 
***************************************************************************/

//=========================================================================
//
//
//                  Section 1: tInputFile Include and Define Statements
//
//
//=========================================================================

#ifndef TINPUTFILE_H
#define TINPUTFILE_H

#include "Headers/tribs_os.h"
#include "Headers/Definitions.h"
#include "Headers/Classes.h"

#ifdef ALPHA_64
  #include <iostream.h>
  #include <fstream.h>
  #include <assert.h>
  #include <string.h>
  #include <stdlib.h>
#elif defined LINUX_32
  #include <iostream>
  #include <fstream>
  #include <cassert>
  #include <string> 
  #include <cstdlib>
#elif defined WIN
  #include <iostream.h>
  #include <fstream.h>
  #include <assert.h>
  #include <string.h>
  #include <stdlib.h>
#else 
  #include <iostream.h>
  #include <fstream.h>
  #include <assert.h>
  #include <string.h>
  #include <stdlib.h>
#endif

using namespace std;

//=========================================================================
//
//
//                  Section 2: tInputFile Class Definitions
//
//
//=========================================================================

class tInputFile
{
public:
  tInputFile();
  tInputFile( const char * );	
  ~tInputFile();
  int    IsItemIn( const char * );
  int    ReadItem( const int &, const char * );
  long   ReadItem( const long &, const char * );
  double ReadItem( const double &, const char * );
  void   ReadItem( char *, const char * );
  void   CloseOldAndOpenNew( const char * ); 
  char*  GetInFileName() { return InFileName; } 

private:
  ifstream infile;
  char InFileName[kMaxNameSize]; 
};

#endif

//=========================================================================
//
//
//                          End of tInputFile.h 
//
//
//=========================================================================
