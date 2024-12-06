#include "src/tHydro/tHydroModel.h"
#include "gtest/gtest.h"
#include <iostream>
#include "src/Headers/globalIO.h"
#include "src/Headers/Inclusions.h"
#include "src/Mathutil/predicates.h"
#include "src/tFlowNet/tKinemat.h"
#include "src/tFlowNet/tReservoir.h" // JECR2014
#include "src/tRasTin/tRainfall.h"
#include "src/tRasTin/tShelter.h" // SKY2008Snow from AJR2007
#include "src/tSimulator/tSimul.h"
#include "src/Headers/TemplDefinitions.h"
#include "src/tHydro/tSnowPack.h" // SKY2008Snow from AJR2007


#ifdef PARALLEL_TRIBS
#include "tGraph/tGraph.h"
#include "tParallel/tParallel.h"
#endif

using namespace std;

Predicates predicate;

// Create global output streams for serial/parallel
tOstream Cout(cout);
tOstream Cerr(cerr);


int argc = 1;
const char* input = "tests/TestInput/happy_jack/src/in_files/happy_jack.in";
char **argv = (char**)input;

// Check command-line arguments
	SimulationControl SimCtrl(argc, argv);

	
	tInputFile InputFile( SimCtrl.infile );
	
	// Preprocessing meteorological data
	tPreProcess PreProcessor( &SimCtrl, InputFile );
	
	// Timer, checks environmental variables 
	tRunTimer Timer( InputFile );
	
	// Creating Mesh 
	// Option 9 meshbuilder files can be handled like others in serial
	// so no special order has to be imposed as it is with the parallel simulation
	tMesh<tCNode> BasinMesh( &SimCtrl, InputFile );
	
	
	tKinemat Flow( &SimCtrl, &BasinMesh, InputFile, &Timer );

	
	tResample RsmplMaster( &SimCtrl, &BasinMesh ); 

	
	tShelter Shelter( &SimCtrl, &BasinMesh, InputFile );// SKY2008Snow from AJR2007

	
	tCOutput<tCNode> Output( &SimCtrl, &BasinMesh, InputFile, 
							 &RsmplMaster, &Timer );
	
	// Creating WaterBalance class
	tWaterBalance Balance( &SimCtrl, &BasinMesh, InputFile);
	
	
	tHydroModel Moisture( &SimCtrl, &BasinMesh, InputFile, 
						  &RsmplMaster, &Balance, &Timer );



TEST(HelloTest, BasicAssertions) {
  // Expect two strings not to be equal.
  EXPECT_STRNE("hello", "world");
  // Expect equality.
  EXPECT_EQ(7 * 6, 42);
}

TEST(CoolTest, BasicAsserions){
   Moisture.get_Lower_Moist(2.0,100.5);
   
}  