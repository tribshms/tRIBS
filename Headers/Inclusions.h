/***************************************************************************
**
**  		      tRIBS Distributed Hydrologic Model
**
**              TIN-based Real-time Integrated Basin Simulator
**		         Ralph M. Parsons Laboratory
**  		    Massachusetts Institute of Technology
**  
**
**  Inclusions.h: General tRIBS Include Statements
**
***************************************************************************/

#ifndef INCLUSIONS_H
#define INCLUSIONS_H

#include "Headers/tribs_os.h"

// INCLUDED LIBRARY HEADER FILES

#ifdef ALPHA_64
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <string.h>
#include <assert.h>
#include <memory.h>

#elif defined LINUX_32
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cassert>
#include <memory>

#elif defined MAC
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cassert>
#include <memory>

#elif defined WIN
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <string.h>
#include <assert.h>
#include <memory.h>

#else 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <string.h>
#include <assert.h>
#include <memory.h>

#endif

// INCLUDED FILES 

#include "Headers/Definitions.h"
#include "Headers/Classes.h"
#include "Headers/globalFns.h"
#include "tHydro/tHydroModel.h"
#include "tFlowNet/tFlowNet.h"
#include "tMesh/tMesh.h"
#include "tMeshList/tMeshList.h"
#include "tMeshElements/meshElements.h"
#include "tPtrList/tPtrList.h"
#include "tList/tList.h"
#include "tListInputData/tListInputData.h"
#include "tCNode/tCNode.h"
#include "Mathutil/mathutil.h"
#include "tArray/tArray.h"
#include "tArray/tMatrix.h"
#include "tInOut/tInputFile.h"
#include "tInOut/tOutput.h"
#include "tInOut/tOstream.h"
#include "tRasTin/tResample.h"
#include "tStorm/tStorm.h"
#include "tRasTin/tRainfall.h"
#include "tRasTin/tRainGauge.h"
#include "tRasTin/tVariant.h"
#include "tRasTin/tInvariant.h"
#include "tSimulator/tRunTimer.h"
#include "tSimulator/tControl.h"
#include "tHydro/tIntercept.h"
#include "tHydro/tEvapoTrans.h"
#include "tHydro/tHydroMet.h"
#include "tHydro/tHydroMetStoch.h"
#include "tHydro/tHydroMetConvert.h"
#include "Mathutil/predicates.h"
#include "Mathutil/geometry.h"
#include "tHydro/tWaterBalance.h"
#include "tSimulator/tPreProcess.h"

#endif

//=========================================================================
//
//
//                   End of Inclusions.h  
//
//
//=========================================================================
