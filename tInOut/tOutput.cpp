/***************************************************************************
**
**  		     tRIBS Distributed Hydrologic Model
**
**              TIN-based Real-time Integrated Basin Simulator
**		         Ralph M. Parsons Laboratory
**  		    Massachusetts Institute of Technology
**  
**
**  tOutput.cpp: Functions for output objects for classes tOutput and 
**               tCOutput (see tOutput.h)
**
*************************************************************************/

#include "tInOut/tOutput.h"
#include "Headers/globalIO.h"
#include "Headers/Inclusions.h"

#ifdef PARALLEL_TRIBS
#include "tGraph/tGraph.h"
#include "tParallel/tParallel.h"
#endif

//=========================================================================
//
//
//                  Section 1: tOutput Constructors/Destructors
//
//
//=========================================================================

/*************************************************************************
**
**  tOutput::Constructor (1)
**
**  The constructor takes two arguments, a pointer to the grid mesh and
**  a reference to an open input file. It reads the base name for the
**  output files from the input file, and opens and initializes these.
**
*************************************************************************/
template< class tSubNode >
tOutput<tSubNode>::tOutput(SimulationControl *simCtrPtr,
						   tMesh<tSubNode> * gridPtr, tInputFile &infile,
						   tResample *resamp)
{
	assert(gridPtr > 0);
	g = gridPtr;  
	respPtr = resamp;
	simCtrl = simCtrPtr;
	
	Cout<<"\nOutput Files:"<<endl<<endl;
	infile.ReadItem(baseName, "OUTFILENAME" );    //basename of output
	infile.ReadItem(nodeFile, "NODEOUTPUTLIST");  //pathname of node output

	// Binary viz output
   vizOption = infile.ReadItem(vizOption, "OPTVIZ");
   if (this->vizOption > 0)
	    infile.ReadItem(vizName, "OUTVIZFILENAME" );
	
#ifdef PARALLEL_TRIBS
  // Nodes, edges, triangles, and Z files are only written
  // by the Master node
  if (tParallel::isMaster()) {
#endif

	char nodesext[10] = ".nodes";
	char edgesext[10] = ".edges";
	char trisext[10] = ".tri";
	char zofsext[10] = ".z";
	
	CreateAndOpenFile( &nodeofs, nodesext );
	CreateAndOpenFile( &edgofs,  edgesext );
	CreateAndOpenFile( &triofs,  trisext );
	CreateAndOpenFile( &zofs,    zofsext );
	
	nodeofs.setf( ios::fixed, ios::floatfield);
	zofs.setf   ( ios::fixed, ios::floatfield);
	
#ifdef PARALLEL_TRIBS
  }
#endif

	ReadNodeOutputList(); 
	CreateAndOpenPixel();
	dynvars = NULL;
}

/*************************************************************************
**
**  tOutput::Constructor (2)
**
**  The constructor takes two arguments, a pointer to the grid mesh and
**  a reference to an open input file. It reads the base name for the
**  output files from the input file, and opens and initializes these.
**
*************************************************************************/
template< class tSubNode >
tOutput<tSubNode>::tOutput( SimulationControl *simCtrPtr,
							tMesh<tSubNode> * gridPtr, 
                            tInputFile &infile, tResample *resamp,
                            tRunTimer * timptr )
{
	assert(gridPtr != nullptr); //Updated to new c++ standards -WR
	g = gridPtr;        
	timer = timptr;    
	respPtr = resamp;
	simCtrl = simCtrPtr;
	
	Cout<<"\nOutput Files:"<<endl<<endl;
	infile.ReadItem( baseName, "OUTFILENAME" );          
	infile.ReadItem( nodeFile, "NODEOUTPUTLIST");  
	
	// Binary viz output
   vizOption = infile.ReadItem(vizOption, "OPTVIZ");
   if (this->vizOption > 0)
	    infile.ReadItem(vizName, "OUTVIZFILENAME" );
	
#ifdef PARALLEL_TRIBS
  // Nodes, edges, triangles, and Z files are only written
  // by the Master node
  if (tParallel::isMaster()) {
#endif

	char nodesext[10] = ".nodes";
	char edgesext[10] = ".edges";
	char trisext[10] = ".tri";
	char zofsext[10] = ".z";
	
	CreateAndOpenFile( &nodeofs, nodesext);
	CreateAndOpenFile( &edgofs,  edgesext);
	CreateAndOpenFile( &triofs,  trisext);
	CreateAndOpenFile( &zofs,    zofsext);
	
	nodeofs.setf( ios::fixed, ios::floatfield);
	zofs.setf   ( ios::fixed, ios::floatfield);
	
#ifdef PARALLEL_TRIBS
  }
#endif

	ReadNodeOutputList(); 
	CreateAndOpenPixel();
	dynvars = NULL;
}

template< class tSubNode >
tOutput<tSubNode>::~tOutput()
{
	g     = NULL;
	timer = NULL;
	if (nodeList)
		delete [] nodeList;
	if (uzel)
		delete [] uzel;
	if (pixinfo)
		delete [] pixinfo;
	if (dynvars)
		delete [] dynvars;
	
	Cout<<"tOutput Object has been destroyed..."<<endl<<flush;
}

//=========================================================================
//
//
//                  Section 2: tOutput Open and Write Functions
//
//
//=========================================================================

/*************************************************************************
**
**  tOutput::CreateAndOpenFile
**
**  Opens the output file stream pointed by theOFStream, giving it the
**  name <baseName><extension>, and checks to make sure that the ofstream
**  is valid.
**
**  Input:  theOFStream -- ptr to an ofstream object
**          extension -- file name extension (e.g., ".nodes")
**  Output: theOFStream is initialized to create an open output file
**  Assumes: extension is a null-terminated string, and the length of
**           baseName plus extension doesn't exceed kMaxNameSize+6
**           (ie, the extension is expected to be <= 6 characters)
**
*************************************************************************/
template< class tSubNode >
void tOutput<tSubNode>::CreateAndOpenFile( ofstream *theOFStream,
                                           char *extension )
{
	char fullName[kMaxNameSize+6];
	
	strcpy( fullName, baseName );
	strcat( fullName, extension );
	
#ifdef PARALLEL_TRIBS
// Add processor extension if running in parallel
  char procex[10];
  sprintf( procex, ".%-d", tParallel::getMyProc());
  strcat(fullName, procex);
#endif

	theOFStream->open( fullName );
	
	if ( !theOFStream->good() )
		cerr << "File "<<fullName<<" not created.";
	
/*SMM
	Cout<<"Creating Output File: \t '"<<fullName<<"' "<<endl;
*/
	return;
}

template< class tSubNode >
void tOutput<tSubNode>::CreateAndOpenVizFile( ofstream *theOFStream,
                                           char *extension )
{ 
	char procex[10];
	char fullName[kMaxNameSize+6];
   
	strcpy( fullName, vizName );
	strcat( fullName, extension );
#ifdef PARALLEL_TRIBS
	sprintf( procex, ".%-d", tParallel::getMyProc());
	strcat(fullName, procex);
#else
	sprintf( procex, "");
	strcat(fullName, procex);
#endif 
  
	theOFStream->open( fullName, ios::binary );
  
	if ( !theOFStream->good() )
		cerr << "File "<<fullName<<" not created.";
  
	Cout<<"Creating Output File: \t '"<<fullName<<"' "<<endl;
	return;
} 

/*************************************************************************
**
**  tOutput::ReadNodeOutputList()
**
**  Opens and Reads the node list from a *.nol file whose structure is:
**
**  Number of Nodes
**  Node1 Node2 Node3 Node4 Node5 ...
**
*************************************************************************/
template< class tSubNode >
void tOutput<tSubNode>::ReadNodeOutputList() {
	
	ifstream readNOL(nodeFile);
	if (!readNOL) {
		Cout << "\nFile "<<nodeFile<<" not found..."<<endl;
		Cout<<"\tAttention: The specified file with node IDs does not "
			<<"exist.\n\t\t   No node output will be written."<<endl<<endl;
		numNodes = 0;
		return;
	}
	
	readNOL >> numNodes;
	nodeList = new int[numNodes];
	uzel = new tSubNode*[numNodes];
	pixinfo = new ofstream[numNodes];
	
#ifdef PARALLEL_TRIBS
  // Initialize to NULL
  for (int i = 0; i < numNodes; i++)
    uzel[i] = NULL;
#endif

	for (int i = 0; i < numNodes; i++) {
		readNOL >> nodeList[i]; 
	}
	
	readNOL.close();
	return;
}

/*************************************************************************
**
**  tOutput::CreateAndOpenPixel()
**
**  Write the header for the *.pixel output file.
** 
*************************************************************************/
template< class tSubNode >
void tOutput<tSubNode>::CreateAndOpenPixel()
{
	if ( nodeList ) {
		char pixelext[10] = ".pixel";
		char nodeNum[10], pixelnode[100];
		
      //SMM - Set interior nodes, added 08132008
      SetInteriorNode();

		for (int i = 0; i < numNodes; i++) {
#ifdef PARALLEL_TRIBS
         // Check if node is on this processor
         if ( (uzel[i] != NULL) && (nodeList[i] >= 0) ) {
#else
			if (nodeList[i] >= 0) {
#endif
				sprintf(nodeNum,"%d",nodeList[i]);
				strcpy(pixelnode, nodeNum);
				strcat(pixelnode, pixelext);		 
				
				CreateAndOpenFile( &pixinfo[i], pixelnode );
				
				if (simCtrl->Header_label == 'Y') {
					pixinfo[i]<<"1-NodeID  "
					<<"2-Time,hr  "
					<<"3-Nwt,mm  "
					<<"4-Nf,mm  "
					<<"5-Nt,mm  "
					<<"6-Mu,mm  "
					<<"7-Mi,mm  "
					<<"8-QpOut,mm/h  "
					<<"9-QpIn,mm/h  "
					<<"10-Trnsm,m2/h  "
					<<"11-GWflx,m3/h  "
					<<"12-Srf,mm  "
					<<"13-Rain,mm/h  "
					//<<"14-SoilMoist,d/l " 
					<<"14-SoilMoist,[]  " // SKYnGM2008LU
					//<<"15-RootMoist,d/l  " 
					<<"15-RootMoist,[]  " // SKYnGM2008LU
					<<"16-AirT,oC  "
					<<"17-DewT,oC  "
					<<"18-SurfT,oC  "
					<<"19-SoilT,oC  "
					<<"20-Press,Pa  "
					//<<"21-RelHum,d/l  " 
					<<"21-RelHum,[]  " // SKYnGM2008LU
					//<<"22-SkyCov,d/l  "
					<<"22-SkyCov,[]  " // SKYnGM2008LU
					<<"23-Wind,m/s  "
					<<"24-NetRad,W/m2  "
					<<"25-ShrtRadIn,W/m2  "
					<<"26-ShrtRadIn_dir,W/m2  "
					<<"27-ShrtRadIn_dif,W/m2  "
					<<"28-ShortAbsbVeg,W/m2  "
					<<"29-ShortAbsbSoi,W/m2  "
					<<"30-LngRadIn,W/m2  "
					<<"31-LngRadOut,W/m2,  "
					<<"32-PotEvp,mm/h  "
					<<"33-ActEvp,mm/h  "
					<<"34-EvpTtrs,mm/h  "
					<<"35-EvpWetCan,mm/h  "
					<<"36-EvpDryCan,mm/h  "
					<<"37-EvpSoil,mm/h  "
					<<"38-Gflux,W/m2  "
					<<"39-HFlux,W/m2  "
					<<"40-Lflux,W/m2  "
					<<"41-NetPrecip,mm/hr  "
					

					// SKY2008Snow from AJR2007
					<<"42-LiqWE,cm  "		//added by AJR 2007 @ NMT
					<<"43-IceWE,cm  "		//added by AJR 2007 @ NMT
					<<"44-SnWE,cm  "		//added by AJR 2007 @ NMT
					<<"45-SnSub,cm "		//added by CJC2020
					<<"46-SnEvap,cm "		//added by CJC2020
					<<"47-U,kJ/m2 "		//added by AJR 2007 @ NMT
					<<"48-RouteWE,cm  "		//added by AJR 2007 @ NMT
					<<"49-SnTemp,C  "		//added by AJR 2007 @ NMT
					<<"50-SurfAge,h "		//added by AJR 2007 @ NMT
					<<"51-DU,kJ/m2/etistep "		//added by AJR 2007 @ NMT
					<<"52-snLHF,kJ/m2/etistep "		//added by AJR 2007 @ NMT
					<<"53-snSHF,kJ/m2/etistep "		//added by AJR 2007 @ NMT
					<<"54-snGHF,kJ/m2/etistep "		//added by AJR 2007 @ NMT
					<<"55-snPHF,kJ/m2/etistep "		//added by AJR 2007 @ NMT
					<<"56-snRLout,kJ/m2/etistep "		//added by AJR 2007 @ NMT
					<<"57-snRLin,kJ/m2/etistep "		//added by AJR 2007 @ NMT
					<<"58-snRSin,kJ/m2/etistep "		//added by AJR 2007 @ NMT
					<<"59-Uerror,kJ/m2/etistep "		//added by AJR 2007 @ NMT
					<<"60-IntSWEq,cm "		//added by AJR 2007 @ NMT
					<<"61-IntSub,cm "		//added by AJR 2007 @ NMT
					<<"62-IntSnUnload,cm "		//added by AJR 2007 @ NMT

					// SKY2008Snow
					<<"63-CanStorage,mm  "
					<<"64-CumIntercept,mm  "
					<<"65-Interception,mm  "
					<<"66-Recharge,mm/hr  "
					<<"67-RunOn,mm  "
					<<"68-Srf_Hour,mm  "
					<<"69-Qstrm,m3/s  "
					<<"70-Hlevel,m"
					//<<"42-CanStorg,mm  "
					//<<"43-CumIntercept,mm  "
					//<<"44-Intercept,mm  "
					//<<"45-Recharge,mm/hr  "
					//<<"46-Runon,mm  "
					//<<"47-Srf_Hour,mm  "
					//<<"48-Qstrm,m3/s  "
					//<<"49-Hlev,m"

					// SKYnGM2008LU
					<<"71-CanStorParam,mm  "
					<<"72-IntercepCoeff,[]  "
					<<"73-ThroughFall,[]  "
					<<"74-CanFieldCap,mm  "
					<<"75-DrainCoeff,mm/hr  "
					<<"76-DrainExpPar,1/mm  "
					<<"77-LandUseAlb,[] "
					<<"78-VegHeight,m "
					<<"79-OptTransmCoeff,[]"
					<<"80-StomRes,s/m"
					<<"81-VegFraction,[] "
					<<"82-LeafAI,[] "

					<<"\n";
				}
				pixinfo[i].setf( ios::right, ios::adjustfield );
				pixinfo[i].setf( ios::fixed, ios::floatfield);  
			}
		}
	}
	return;
}

/*************************************************************************
**
**  tOutput::CreateAndOpenDynVar()
**
**  Opens a number of files for selected dynamic variables 
** 
*************************************************************************/
template< class tSubNode >
void tOutput<tSubNode>::CreateAndOpenDynVar()
{
	if (g->getNodeList()->getActiveSize() < 0) {
		
		// Allocate memory for ofstream objects
		//dynvars = new ofstream[19];
		dynvars = new ofstream[36]; // SKY2008Snow, AJR2008
		
		char Nwt[10] = "_Nwt";
		char Nf[10]  = "_Nf";
		char Nt[10]  = "_Nt";
		char Mu[10]  = "_Mu";
		char Mi[10]  = "_Mi";
		char QOut[10] = "_QOut";
		char QIn[10]  = "_QIn";
		char Ru[10]    = "_Ru";
		char Trsm[10]  = "_Trsm";
		char GWflx[10] = "_GWflx";
		char Srf[10]   = "_Srf";
		char SMS[10]   = "_SMS";
		char SMRT[10]  = "_SMRT";
		char ET[10]    = "_ET";
		char Edc[10]  = "_Edc";
		char Ewc[10]  = "_Ewc";
		char Eso[10]  = "_Eso";
		char NetP[10] = "_NetP";
		char Qbot[10] = "_Qbot";
	
		// SKY2008Snow from AJR2007
		char SnWE[10] = "_SnWE";	    //added by AJR 2007 @ NMT
		char SnTempC[10] = "_snTempC";  //added by AJR 2007 @ NMT
		char IceWE[10] = "_IcWE";	    //added by AJR 2007 @ NMT
		char LiqWE[10] = "_LiWE";	    //added by AJR 2007 @ NMT
		char DU[10] = "_DU";	    //added by AJR 2007 @ NMT
		char Upack[10] = "_UTOT";	    //added by AJR 2007 @ NMT
		char snLHF[10] = "_sLHF";	    //added by AJR 2007 @ NMT
		char snSHF[10] = "_sSHF";	    //added by AJR 2007 @ NMT
		char snGHF[10] = "_sGHF";	    //added by AJR 2007 @ NMT
		char snPHF[10] = "_sPHF";	    //added by AJR 2007 @ NMT
		char snRLo[10] = "_sRLo";	    //added by AJR 2007 @ NMT
		char snRLi[10] = "_sRLi";	    //added by AJR 2007 @ NMT
		char snRSi[10] = "_sRSi";	    //added by AJR 2007 @ NMT
		char Uerr[10] = "_Uerr";	    //added by AJR 2007 @ NMT
		char IntSn[10] = "_IntSn";	    //added by AJR 2007 @ NMT
		char IntSub[10] = "_IntSub";    //added by AJR 2007 @ NMT
		char IntUnl[10] = "_IntUnl";    //added by AJR 2007 @ NMT


		//for (int i = 0; i < 19; i++) {
		for (int i = 0; i < 36; i++) { // SKY2008Snow from AJR2007
			if (i == 0)
				CreateAndOpenFile( &dynvars[i], Nwt);
			else if (i == 1)
				CreateAndOpenFile( &dynvars[i], Nf);
			else if (i == 2)
				CreateAndOpenFile( &dynvars[i], Nt);
			else if (i == 3)
				CreateAndOpenFile( &dynvars[i], Mu);
			else if (i == 4)
				CreateAndOpenFile( &dynvars[i], Mi);
			else if (i == 5)
				CreateAndOpenFile( &dynvars[i], QOut);
			else if (i == 6)
				CreateAndOpenFile( &dynvars[i], QIn);
			else if (i == 7)
				CreateAndOpenFile( &dynvars[i], Ru);
			else if (i == 8)
				CreateAndOpenFile( &dynvars[i], Trsm);
			else if (i == 9)
				CreateAndOpenFile( &dynvars[i], GWflx);
			else if (i == 10)
				CreateAndOpenFile( &dynvars[i], Srf);
			else if (i == 11)
				CreateAndOpenFile( &dynvars[i], SMS);
			else if (i == 12)
				CreateAndOpenFile( &dynvars[i], SMRT);
			else if (i == 13)
				CreateAndOpenFile( &dynvars[i], ET);
			else if (i == 14)
				CreateAndOpenFile( &dynvars[i], Edc);
			else if (i == 15)
				CreateAndOpenFile( &dynvars[i], Ewc);
			else if (i == 16)
				CreateAndOpenFile( &dynvars[i], Eso);
			else if (i == 17)
				CreateAndOpenFile( &dynvars[i], NetP);
			else if (i == 18)

			// SKY2008Snow from AJR2007
				CreateAndOpenFile( &dynvars[i], SnWE);	    //added by AJR 2007 @ NMT
			else if (i == 19)
				CreateAndOpenFile( &dynvars[i], SnTempC);   //added by AJR 2007 @ NMT
			else if (i == 20)
				CreateAndOpenFile( &dynvars[i], IceWE);	    //added by AJR 2007 @ NMT
			else if (i == 21)
				CreateAndOpenFile( &dynvars[i], LiqWE);	    //added by AJR 2007 @ NMT
			else if (i == 22)
				CreateAndOpenFile( &dynvars[i], DU);	    //added by AJR 2007 @ NMT
			else if (i == 23 )
				CreateAndOpenFile( &dynvars[i], Upack);	    //added by AJR 2007 @ NMT
			else if (i == 24)
				CreateAndOpenFile( &dynvars[i], snLHF);	    //added by AJR 2007 @ NMT
			else if (i == 25)
				CreateAndOpenFile( &dynvars[i], snSHF);	    //added by AJR 2007 @ NMT
			else if (i == 26)
				CreateAndOpenFile( &dynvars[i], snGHF);	    //added by AJR 2007 @ NMT
			else if (i == 27)
				CreateAndOpenFile( &dynvars[i], snPHF);	    //added by AJR 2007 @ NMT
			else if ( i == 28)
				CreateAndOpenFile( &dynvars[i], snRLo);	    //added by AJR 2007 @ NMT
			else if ( i == 29 )
				CreateAndOpenFile( &dynvars[i], snRLi);	    //added by AJR 2007 @ NMT
			else if ( i == 30 )
				CreateAndOpenFile( &dynvars[i], snRSi);	    //added by AJR 2007 @ NMT
			else if (i == 31)
				CreateAndOpenFile( &dynvars[i], Uerr);	    //added by AJR 2007 @ NMT
			else if (i == 32)
				CreateAndOpenFile( &dynvars[i], IntSn);	    //added by AJR 2007 @ NMT
			else if (i == 33)
				CreateAndOpenFile( &dynvars[i], IntSub);    //added by AJR 2007 @ NMT
			else if (i == 34)
				CreateAndOpenFile( &dynvars[i], IntUnl);    //added by AJR 2007 @ NMT
			else if (i == 35)

				CreateAndOpenFile( &dynvars[i], Qbot);
			dynvars[i].setf( ios::right, ios::adjustfield );
			dynvars[i].setf( ios::fixed, ios::floatfield);
		}
	}
	return;
}

/*************************************************************************
**
**  tOutput::WriteOutput
**
**  This function writes information about the mesh to four files called
**  name.nodes, name.edges, name.tri, and name.z, where "name" is a
**  name that the user has specified in the input file and which is
**  stored in the data member baseName.
**
**  Input: time -- time of the current output time-slice
**  Output: the node, edge, and triangle ID numbers are modified so that
**          they are numbered according to their position on the list
**  Assumes: the four file ofstreams have been opened by the constructor
**           and are valid
**
*************************************************************************/
template< class tSubNode >
void tOutput<tSubNode>::WriteOutput( double time )
{
	tNode * cn;
	tEdge * ce;
	tTriangle * ct;

	tMeshListIter<tSubNode> niter( g->getNodeList() );
	tMeshListIter<tEdge>    eiter( g->getEdgeList() );
	tListIter<tTriangle>    titer( g->getTriList() );
	
	int nnodes = g->getNodeList()->getSize();
	int nedges = g->getEdgeList()->getSize();
	int ntri   = g->getTriList()->getSize();
	
#ifdef PARALLEL_TRIBS
	// If running parallel, sum sizes on all processors
	int nGlobalActiveNodes = g->getNodeList()->getGlobalActiveSize();
	int nGlobalActiveEdges = g->getEdgeList()->getGlobalActiveSize();
	int nActiveNodes = g->getNodeList()->getActiveSize();
	int nActiveEdges = g->getEdgeList()->getActiveSize();

	cout<<"Proc " << tParallel::getMyProc()
	    <<": tOutput Characteristics:"<<endl;
  
	cout<<"Proc " << tParallel::getMyProc()
	    <<": Number of nodes: \t\t"<<nnodes<<endl;
	cout<<"Proc " << tParallel::getMyProc()
	    <<": Number of edges: \t\t"<<nedges<<endl;
  
	cout<<"Proc " << tParallel::getMyProc()
	    <<": Number of active nodes: \t"<<nActiveNodes<<endl<<flush;
	cout<<"Proc " << tParallel::getMyProc()
	    <<": Number of active edges: \t"<<nActiveEdges<<endl<<flush;
  
	cout<<"Proc " << tParallel::getMyProc()
	    <<": Number of global active nodes: \t"<<nGlobalActiveNodes<<endl<<flush;
	cout<<"Proc " << tParallel::getMyProc()
	    <<": Number of global active edges: \t"<<nGlobalActiveEdges<<endl<<flush;

#else

	int nActiveNodes = g->getNodeList()->getActiveSize();
	int nActiveEdges = g->getEdgeList()->getActiveSize();

	cout<<"\ntOutput Characteristics:"<<endl<<endl;
	cout<<"Number of nodes: \t\t"<<nnodes<<endl;
	cout<<"Number of edges: \t\t"<<nedges<<endl;
	cout<<"Number of triangles: \t\t"<<ntri<<endl;
	cout<<"Number of active nodes: \t"
		<<nActiveNodes<<endl<<flush;
	cout<<"Number of active edges: \t"
		<<nActiveEdges<<endl<<flush;

#endif

	if (nnodes > 0 && nedges > 0 && ntri > 0) {
#ifdef PARALLEL_TRIBS
  // Master node ONLY writes these files
  if (tParallel::isMaster()) {
#endif

	// Formating floating point numbers in output
	nodeofs.setf(ios::fixed, ios::floatfield);
	zofs.setf   (ios::fixed, ios::floatfield);
	edgofs.setf (ios::fixed, ios::floatfield);
	triofs.setf (ios::fixed, ios::floatfield);
	
	// Write 'node' file and 'z' file 
	nodeofs<<" "<<time<<"\n"<<nnodes<<"\n";
	zofs   <<" "<<time<<"\n"<<nnodes<<"\n";
	
	for ( cn=niter.FirstP(); !(niter.AtEnd()); cn=niter.NextP() ) {
		nodeofs<<cn->getX()<<" "<<cn->getY()<<" "
		<<cn->getEdg()->getID()<<" "<<cn->getBoundaryFlag()<<"\n";
		zofs   <<cn->getZ()<<"\n";
	}
	
	// Write 'edge' file 
	edgofs<<" "<<time<<"\n"<<nedges<<"\n";
	for ( ce=eiter.FirstP(); !(eiter.AtEnd()); ce=eiter.NextP() )
		edgofs <<ce->getOriginPtrNC()->getID()<<" "
			<<ce->getDestinationPtrNC()->getID()<<" "
			<<ce->getCCWEdg()->getID()<<"\n";
	
	// Write 'triangle' file 
	int i;
	triofs<<" "<<time<<"\n"<<ntri<<"\n";
	for ( ct=titer.FirstP(); !(titer.AtEnd()); ct=titer.NextP() )  {
		for ( i=0; i<=2; i++ )
			triofs<<ct->pPtr(i)->getID()<<" ";
		for ( i=0; i<=2; i++ )  {
			if ( ct->tPtr(i) ) 
				triofs<<ct->tPtr(i)->getID()<<" ";
			else 
				triofs<<"-1 ";
		}
		triofs << ct->ePtr(0)->getID() << " " 
			<< ct->ePtr(1)->getID() << " " 
			<< ct->ePtr(2)->getID() << "\n";
	}
	
#ifdef PARALLEL_TRIBS
  }
#endif
	}

	// Set interior nodes
	SetInteriorNode();
	
	return;
}

/*************************************************************************
**
**  tOutput::SetInteriorNode()
**
**  Initializes pointers to basin interior outlets 
** 
*************************************************************************/
template< class tSubNode >
void tOutput<tSubNode>::SetInteriorNode()
{
	tSubNode * cnn;
	tMeshListIter<tSubNode> niter( g->getNodeList() );
	
	if (nodeList) {  // <--- If the node list in NOT empty
		for (int i=0; i < numNodes; i++) {
			if (nodeList[i] >= 0) {
#ifdef PARALLEL_TRIBS
           // The processor with the node of interest creates the file
           for( cnn=niter.FirstP(); niter.IsActive(); cnn=niter.NextP() ) {
#else
				for ( cnn=niter.FirstP(); !(niter.AtEnd()); cnn=niter.NextP() ) {
#endif
					if (cnn->getID() == nodeList[i]) {
						uzel[i] = cnn;   // <=== Defining node of interest 
						cout<<"\nNode of Interest ID: \t"<<nodeList[i]
							<<" has been set up..."<<endl<<flush;
					}
				}
			}
		}
	}
	return;
}

//=========================================================================
//
//
//                  Section 3: tOutput Void Write Functions
//
//
//=========================================================================
template< class tSubNode >
void tOutput<tSubNode>::WriteNodeData( double ) 
{
	cout<<"tOutput:WriteNodeData VOID function!"<<endl; 
}

template< class tSubNode >
void tOutput<tSubNode>::WriteNodeData( double , tResample *)
{
	cout<<"tOutput:WriteNodeData VOID function!"<<endl; 
}

template< class tSubNode >
void tOutput<tSubNode>::WriteGeometry( tResample *)
{
	cout<<"tOutput:WriteNodeData VOID function!"<<endl; 
}

template< class tSubNode >
void tOutput<tSubNode>::WriteDynamicVars( double )
{
	cout<<"tOutput:WriteDynamicVars VOID function!"<<endl; 
}

template< class tSubNode >
void tOutput<tSubNode>::WriteDynamicVarsBinary( double )
{
	cout<<"tOutput:WriteDynamicVarsBinary VOID function!"<<endl; 
}

template< class tSubNode >
void tOutput<tSubNode>::WriteDynamicVar( double time ) 
{
	cout<<"tOutput:WriteDynamicVar VOID function!"<<endl; 
}

template< class tSubNode >
void tOutput<tSubNode>::WritePixelInfo( double )
{
	cout<<"tOutput:WritePixelInfo VOID function!"<<endl; 
}

template< class tSubNode >
void tOutput<tSubNode>::end_simulation()
{
	for (int i = 0; i < numNodes; i++)
#ifdef PARALLEL_TRIBS
    // Check if node is on this processor
    if ( (uzel[i] != NULL) && (nodeList[i] >= 0) )
#endif
		pixinfo[i].close();
	return;  
}

//=========================================================================
//
//
//                  Section 4: tCOutput Derived Class (tRIBS Node)
//
//
//=========================================================================

/*************************************************************************
**
**  tCOutput::Constructors and Destructor
**
**  Constructor for the derived tCNode class. In addition to the 
**  inhereted constructor functions, it also opens a _voi file for
**  writing the voronoi mesh points as a ArcInfo generate file
**
*************************************************************************/
template< class tSubNode >
tCOutput<tSubNode>::tCOutput(SimulationControl *simCtrPtr, tMesh<tSubNode> *g, 
							 tInputFile &infile, tResample *resamp, 
							 tRunTimer *timptr) 
: tOutput<tSubNode>(simCtrPtr, g, infile, resamp, timptr)
{
	simCtrl = simCtrPtr;
	char nodeFileO[kMaxNameSize];
	char vorofsext[10] = "_voi";
	char widthsext[10] = "_width";
	char drarsext[16] = "_area";
	
	// Set up interior outlets and open .qout files
	infile.ReadItem( outletName, "OUTHYDROFILENAME");
	infile.ReadItem( nodeFileO, "OUTLETNODELIST" );
	ReadOutletNodeList(nodeFileO);

#ifdef PARALLEL_TRIBS
  // Change the order when running in parallel,
  // so outlets are set per processor
  SetInteriorOutlet();
  CreateAndOpenOutlet();
#else
	CreateAndOpenOutlet();
	SetInteriorOutlet();
#endif
	
	// Create files to output vertices of Voronoi
	// polygons and contributing area values
	this->CreateAndOpenFile( &vorofs, vorofsext);
	this->CreateAndOpenFile( &drareaofs, drarsext );
	this->CreateAndOpenFile( &widthsofs, widthsext );
	
	WriteNodeData( 0, resamp );
}

template< class tSubNode >
tCOutput<tSubNode>::tCOutput(SimulationControl *simCtrPtr, tMesh<tSubNode> *g, 
							 tInputFile &infile, tResample *resamp ) 
: tOutput<tSubNode>(simCtrPtr, g, infile, resamp, this->timptr)
{   
	char vorofsext[10] = "_voi";
	this->CreateAndOpenFile( &vorofs, vorofsext);
}

template< class tSubNode >
tCOutput<tSubNode>::~tCOutput() 
{

	// GMnSKY2008MLE to fix memory leaks
	if (numOutlets > 0) { 
		for (int j=0; j < numOutlets; j++) 
#ifdef PARALLEL_TRIBS
    if ( (Outlets[j] != NULL) && (OutletList[j] > 0) )
#endif
			outletinfo[j].close(); 
		delete [] OutletList; 
		delete [] Outlets; 
		delete [] outletinfo; 
	}       

    Cout<<"tCOutput Object has been destroyed..."<<endl<<flush;
}

/*************************************************************************
**
**  tCOutput::WriteNodeData()
**
**  Writes the voronoi vertices to a file *_voi in a format compatible
**  with ArcInfo generate files for input and transformation into
**  a polygon coverage. The output format should be readable by ArcInfo 
**  & Matlab
**
*************************************************************************/
template< class tSubNode >
void tCOutput<tSubNode>::WriteNodeData( double time )
{
	tSubNode *cn;
	tEdge * firstedg;
	tEdge * curedg; 
	tMeshListIter<tSubNode> ni( this->g->getNodeList() );
	tArray< double > xy, xy1, xy2, xy3, xy4;
	
	// Writing to a file voronoi vertices  of the current node 
	// The output format should be readable by ArcInfo & Matlab
	vorofs.setf( ios::fixed, ios::floatfield);
	if (time == 0) {
		cn = ni.FirstP();
		
		while (ni.IsActive()) {
			
			vorofs<<cn->getID()<<','<<cn->getX()<<','<<cn->getY()<<"\n";
			firstedg = cn->getFlowEdg();
			xy = firstedg->getRVtx();
			vorofs<<xy[0]<<','<<xy[1]<<"\n";
			curedg = firstedg->getCCWEdg();
			while (curedg!=firstedg) {
				xy = curedg->getRVtx();
				vorofs <<xy[0]<<','<<xy[1]<<"\n";
				curedg = curedg->getCCWEdg();
			}
			vorofs<<"END"<<"\n";
			cn = ni.NextP();
		}
		vorofs<<"END"<<"\n";
		
	}
	vorofs.close();
	return;
}

/*************************************************************************
**
**  tCOutput::WriteNodeData( double time, tResample *tresamp )
**
**  - Writes a file containing voronoi vertices: reads arrays
**    from the tResample object
**  
**  - Writes a file containing drainage areas for stream network
** 
*************************************************************************/
template< class tSubNode >
void tCOutput<tSubNode>::WriteNodeData( double time, tResample *tresamp )
{
	tSubNode *cn;
	tMeshListIter<tSubNode> ni( this->g->getNodeList() );
	
	int k, i;
	i = 0;
	
	// Writing to a file voronoi vertices  of the current node 
	// The output format should be readable by ArcInfo & Matlab
	vorofs.setf   (ios::fixed, ios::floatfield);
	drareaofs.setf(ios::fixed, ios::floatfield);
	widthsofs.setf(ios::fixed, ios::floatfield);
	
	if (time == 0) {
		cn = ni.FirstP();
		drareaofs<<"ID\tX\tY\tCArea"<<"\n";
		widthsofs<<"ID\tX\tY\tWidth\tEdgL\tSlp"<<"\n";
		
		while (ni.IsActive()) {
			vorofs<<cn->getID()<<','<<cn->getX()<<','<<cn->getY()<<"\n";
			for (k=0; k < tresamp->nPoints[i]; k++)
				vorofs<<tresamp->vXs[i][k]<<","<<tresamp->vYs[i][k]<<"\n";
			vorofs<<"END"<<"\n";
			
			if (cn->getBoundaryFlag() == 3) {
				drareaofs<<cn->getID()
				<<"\t"<<cn->getX()<<"\t"<<cn->getY()
				<<"\t"<<cn->getContrArea()<<"\n";
				widthsofs<<cn->getID()
					<<"\t"<<cn->getX()<<"\t"<<cn->getY()
					<<"\t"<<cn->getChannelWidth()
					<<"\t"<<cn->getFlowEdg()->getLength()
					<<"\t"<<cn->getFlowEdg()->getSlope()<<"\n";
			}
			i++;
			cn = ni.NextP();
		} 

		// Write geometry and static information for visualizer
		if (this->vizOption == 1)
			WriteGeometry(tresamp);

#ifdef PARALLEL_TRIBS
    // If the last reach is on this processor, report the final
    // outlet node
    if (tGraph::hasLastReach())
#endif

		if (cn->getBoundaryFlag() == 2) {
			drareaofs<<cn->getID()
			<<"\t"<<cn->getX()<<"\t"<<cn->getY()
			<<"\t"<<cn->getContrArea()<<"\n";
			widthsofs<<cn->getID()
				<<"\t"<<cn->getX()<<"\t"<<cn->getY()
				<<"\t"<<cn->getChannelWidth()<<"\t0.0\t0.0"<<"\n";
		}
		vorofs<<"END"<<"\n";
	}
	vorofs.close();
	drareaofs.close();
	widthsofs.close();
	return;
}

/*************************************************************************
**
**  tCOutput::WriteGeometry()
**
**  Writes the voronoi polygon information and static data suitable for
**  the tRIBS reader for EnSight or ParaView
**  After variables giving the number of polygons and number of points
**  making up those polygons comes all points in order by polygon
**  Next is the number of points per polygon which will all the polygons
**  to be created from the points array.  Finally the data is written for
**  each cell or polygon organized by data item across all cells.
**
*************************************************************************/
template< class tSubNode >
void tCOutput<tSubNode>::WriteGeometry(tResample* tresamp)
{
	tSubNode *cn;
	tMeshListIter<tSubNode> ni( this->g->getNodeList() );

	// Write the geometry data for visualizer using tResample because
	// boundary polygons have been clipped
	// Visualizer data is written in binary and organized by variable
	// across all cells rather than all data for a single cell
	ofstream geomofs;
	char geomext[16] = ".tribs";
	this->CreateAndOpenVizFile( &geomofs, geomext );

	// Count the total number of voronoi vertices making up the polygons
	int nActiveNodes = this->g->getNodeList()->getActiveSize();
	int nPoints = 0;
	for (int i = 0; i < nActiveNodes; i++)
		nPoints += tresamp->nPoints[i];

	// Write the counts of nodes (polygons) and points forming the polygons
	BinaryWrite(geomofs, nActiveNodes);
	BinaryWrite(geomofs, nPoints);

	// Write all the voronoi vertices which are points in the unstructured grid
	for (int i = 0; i < nActiveNodes; i++) {
		for (int k = 0; k < tresamp->nPoints[i]; k++) {
			BinaryWrite(geomofs, (float) tresamp->vXs[i][k]);
			BinaryWrite(geomofs, (float) tresamp->vYs[i][k]);
		}
	}

	// Write the edges per cell (polygon)
	for (int i = 0; i < nActiveNodes; i++)
		BinaryWrite(geomofs, tresamp->nPoints[i]);

	// Write the centroid points per cell
	for (cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP())
		BinaryWrite(geomofs, (float) cn->getX());
	for (cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP())
		BinaryWrite(geomofs, (float) cn->getY());

	// Centroid z is the elevation
	for (cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP())
		BinaryWrite(geomofs, (float) cn->getZ());

	// Write the boundary flag per cell
	for (cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP())
		BinaryWrite(geomofs, (float) cn->getBoundaryFlag());

	// WRite the reach number per cell
	for (cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP())
		BinaryWrite(geomofs, (float) cn->getReach());

	// Write the flood status per cell
	for (cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP())
		BinaryWrite(geomofs, (float) cn->getFloodStatus());

	// Write the flow width per cell
	for (cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP())
		BinaryWrite(geomofs, (float) cn->getFlowEdg()->getVEdgLen());

	// Write the flow length per cell
	for (cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP())
		BinaryWrite(geomofs, (float) cn->getFlowEdg()->getLength());

	// Write the contributing area per cell
	for (cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP())
		BinaryWrite(geomofs, (float) cn->getContrArea());

	// Write the slope per cell
	for (cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP())
		BinaryWrite(geomofs, (float) cn->getFlowEdg()->getSlope());

	geomofs.close();
}

/*************************************************************************
**
**  tCOutput::WritePixelInfo()
**
**  Writes the dynamic variables of node of interest to a *.pixel file
**  The output format should be readable by ArcInfo & Matlab 
**
*************************************************************************/
template< class tSubNode >
void tCOutput<tSubNode>::WritePixelInfo( double time )
{ 
	if (this->numNodes > 0) {
		int hour, minute;
		char extension[20];
		
		hour   = (int)floor(time);
		minute = (int)((time-hour)*100);
		sprintf(extension,"%04d.%02d", hour, minute);
		
		// Writing to a file dynamic variables of node of interest  
		// The output format should be readable by ArcInfo & Matlab 
		for (int i = 0; i < this->numNodes; i++) {
#ifdef PARALLEL_TRIBS
  // Doesn't need to be less than active size
      if ( (this->uzel[i] != NULL) && (this->nodeList[i] >= 0) ) {
#else
			if ( this->uzel[i] && this->nodeList[i] < this->g->getNodeList()->getActiveSize()) {
#endif
				this->pixinfo[i]<<setw(8)<<this->nodeList[i]
				<<setw(13)<<extension<<" "
				/* 3 */   <<setw(10)<<(this->uzel[i]->getNwtNew())<<" "
				
				<<setprecision(1)
				<<setw(6)<<this->uzel[i]->getNfNew()<<" "
				/* 5 */	  <<setw(6)<<this->uzel[i]->getNtNew()<<" "
				
				<<setprecision(5)
				<<setw(7)<<this->uzel[i]->getMuNew()<<" "
				<<setw(7)<<this->uzel[i]->getMiNew()<<"   "
				
				<<setprecision(7)
				<<setw(10)<<this->uzel[i]->getQpout()*1.E-6/this->uzel[i]->getVArea()<<"  "
				<<setw(10)<<this->uzel[i]->getQpin() *1.E-6/this->uzel[i]->getVArea()<<"   "   
				/* 10 */  <<setw(10)<<this->uzel[i]->getTransmiss()*1.E-6<<"    "
				<<setw(10)<<this->uzel[i]->getGwaterChng()*1.E-9<<"  "
				
				<<setprecision(5)
				<<setw(8) <<this->uzel[i]->getSrf()<<"  "
				<<setw(10)<<this->uzel[i]->getRain()<<"  "
				<<setw(10)<<this->uzel[i]->getSoilMoistureSC()<<"  "
				/* 15 */  <<setw(10)<<this->uzel[i]->getRootMoistureSC()<<" "
				
				<<setprecision(3)
				<<this->uzel[i]->getAirTemp()<<" "
				<<this->uzel[i]->getDewTemp()<<" "
				<<this->uzel[i]->getSurfTemp()<<" "
				<<this->uzel[i]->getSoilTemp()<<" "
				
				/* 20 */  <<this->uzel[i]->getAirPressure()<<" "
				<<this->uzel[i]->getRelHumid()<<" "
				<<this->uzel[i]->getSkyCover()<<" "
				<<this->uzel[i]->getWindSpeed()<<" "
				<<this->uzel[i]->getNetRad()<<" "
				
				/* 25 */  <<this->uzel[i]->getShortRadIn()<< " "
				<<this->uzel[i]->getShortRadIn_dir()<< " "
				<<this->uzel[i]->getShortRadIn_dif()<< " "
				<<this->uzel[i]->getShortAbsbVeg()<< " "
				<<this->uzel[i]->getShortAbsbSoi()<< " "
				
				/* 30 */  <<this->uzel[i]->getLongRadIn()<< " "
				<<this->uzel[i]->getLongRadOut()<< " "
				<<this->uzel[i]->getPotEvap()<<" "
				<<this->uzel[i]->getActEvap()<<" "
				<<this->uzel[i]->getEvapoTrans()<<" "
				
				/* 35 */  <<this->uzel[i]->getEvapWetCanopy()<<" "
				<<this->uzel[i]->getEvapDryCanopy()<<" "
				<<this->uzel[i]->getEvapSoil()<<" "
				<<this->uzel[i]->getGFlux()<<" "
				<<this->uzel[i]->getHFlux()<<" "
				/* 40 */  <<this->uzel[i]->getLFlux()<<" "
				
				<<this->uzel[i]->getNetPrecipitation()<<" "

				// SKY2008Snow from AJR2007
				<<this->uzel[i]->getLiqWE()<<" "	//added by AJR 2007 @ NMT
				<<this->uzel[i]->getIceWE()<<" "	//added by AJR 2007 @ NMT
				<<(this->uzel[i]->getLiqWE()+this->uzel[i]->getIceWE())<<" "    //added by AJR 2007 @ NMT
				/* 45 */ << this->uzel[i]->getSnSub()<<" "	// added by CJC2020
				<< this->uzel[i]->getSnEvap()<<" "	// added by CJC2020
				<<this->uzel[i]->getUnode()<<" "	//added by AJR 2007 @ NMT
				<<this->uzel[i]->getLiqRouted()<<" "	//added by AJR 2007 @ NMT
				<<this->uzel[i]->getSnTempC()<<" "	//added by AJR 2007 @ NMT
				/* 50 */ <<this->uzel[i]->getCrustAge()<<" "	//added by AJR 2007 @ NMT
				<<this->uzel[i]->getDU()<<" "		//added by AJR 2007 @ NMT
				<<this->uzel[i]->getSnLHF()<<" "	//added by AJR 2007 @ NMT
				<<this->uzel[i]->getSnSHF()<<" "	//added by AJR 2007 @ NMT
				<<this->uzel[i]->getSnGHF()<<" "	//added by AJR 2007 @ NMT
				/* 55 */ <<this->uzel[i]->getSnPHF()<<" "	//added by AJR 2007 @ NMT
				<<this->uzel[i]->getSnRLout()<<" "	//added by AJR 2007 @ NMT
				<<this->uzel[i]->getSnRLin()<<" "	//added by AJR 2007 @ NMT
				<<this->uzel[i]->getSnRSin()<<" "	//added by AJR 2007 @ NMT
				<<this->uzel[i]->getUerror()<<" "	//added by AJR 2007 @ NMT
				/* 60 */ <<this->uzel[i]->getIntSWE()<<" "	//added by AJR 2007 @ NMT
				<<this->uzel[i]->getIntSub()<<" "	//added by AJR 2007 @ NMT
				<<this->uzel[i]->getIntSnUnload()<<" " //added by AJR 2007 @ NMT
			
				<<this->uzel[i]->getCanStorage()<<" "
				<<this->uzel[i]->getCumIntercept()<<" "
				/* 65 */ <<this->uzel[i]->getInterceptLoss()<<" "
				<<this->uzel[i]->getRecharge()<<" "
				<<this->uzel[i]->getRunOn()<<" "
				<<this->uzel[i]->getSrf_Hr()<<" ";
				/* 69 */ 
				if (this->uzel[i]->getBoundaryFlag() == kStream)
					this->pixinfo[i]<<setw(10)<<this->uzel[i]->getQstrm()<<" "
						//<<setw(6)<<this->uzel[i]->getHlevel()<<endl<<flush;
						<<setw(6)<<this->uzel[i]->getHlevel()<<" ";// SKYnGM2008LU
				else
					//this->pixinfo[i]<<"0.0 0.0"<<endl<<flush;
					this->pixinfo[i]<<"0.0 0.0 "; // SKYnGM2008LU
				
				// SKYnGM2008LU
				this->pixinfo[i]<<setprecision(7)
				<<setw(10)<<this->uzel[i]->getCanStorParam()<<" "
				<< this->uzel[i]->getIntercepCoeff()<<" "
				<< this->uzel[i]->getThroughFall()<<" "
				<< this->uzel[i]->getCanFieldCap()<<" "
				/* 75 */ << this->uzel[i]->getDrainCoeff()<<" "
				<< this->uzel[i]->getDrainExpPar()<<" "
				<< this->uzel[i]->getLandUseAlb()<<" "
				<< this->uzel[i]->getVegHeight()<<" "
				<< this->uzel[i]->getOptTransmCoeff()<<" "
				/* 80 */ << this->uzel[i]->getStomRes()<<" "
				<< this->uzel[i]->getVegFraction()<<" "
				<< this->uzel[i]->getLeafAI()<<endl<< flush;

			}
		}
	}
	return;
}

/*************************************************************************
**
**  tCOutput::WriteDynamicVars()
**
**  Writes a file containing dynamic variables for all the nodes
**  The option is turned on with the '-R' option in the command line
**  
*************************************************************************/
template< class tSubNode >
void tCOutput<tSubNode>::WriteDynamicVars( double time )
{
	tSubNode *cn;
	tMeshListIter<tSubNode> ni( this->g->getNodeList() );
	
#ifdef PARALLEL_TRIBS
   int nActiveNodes = this->g->getNodeList()->getGlobalActiveSize();
#else
   int nActiveNodes = this->g->getNodeList()->getActiveSize();
#endif
	int nnodes       = this->g->getNodeList()->getSize();
	
	int hour, minute;
	char extension[20];
	
	hour   = (int)floor(time);
	minute = (int)floor((time-hour)*60);

	// SKY2008Snow from AJR2007
	if(simCtrl->Header_label == 'Y'){
		cout<<"\n\tHOUR = "<<hour<<"\tMINUTE = "<<minute<<"\n";
		//cout<<"\ttCOutput:     Time to write vars; nActiveNodes = "
		//   <<nActiveNodes<<";  nTotalNodes = "<<nnodes<<"\n";
	}

	sprintf(extension,".%04d_%02dd", hour, minute);
	this->CreateAndOpenFile( &arcofs, extension);  //Opens file for writing
	
	if (simCtrl->Header_label == 'Y') {
		arcofs<<"ID"<<','<<"Z"<<','<<"S"<<','<<"CAr"<<','<<"Nwt"<<','<<"Mu"<<','
		<<"Mi"<<','<<"Nf"<<','<<"Nt"<<','<<"Qpout"<<','<<"Qpin"<<','
		<<"Srf"<<','<<"Rain"<<','
		
		// SKY2008Snow from AJR2007
		<<"SWE"<<','//added by AJR 2007 @ NMT
		<<"ST"<<','<<"IWE"<<','<<"LWE"<<','<<"SnSub"<<','<<"SnEvap"<<','<<"DU"<<','<<"Upack"<<','//added by AJR 2007 @ NMT // Adjusted by CJC2020
		<<"sLHF"<<','<<"sSHF"<<','<<"sGHF"<<','<<"sPHF"<<','//added by AJR 2007 @ NMT
		<<"sRLo"<<','<<"sRLi"<<','<<"sRSi"<<','<<"Uerr"<<','//added by AJR 2007 @ NMT
		<<"IntSWE"<<','<<"IntSub"<<','<<"IntUnl"<<','//added by AJR 2007 @ NMT	
		
		<<"SoilMoist"<<','<<"RootMoist"<<','<<"CanStorage"<<','
		<<"ActEvp"<<','<<"EvpSoil"<<','<<"ET"<<','<<"GFlux"<<','
		<<"HFlux"<<','<<"LFlux"<<','<<"Qstrm"<<','<<"Hlev"<<','
		//<<"FlwVlc";
		<<"FlwVlc"<<','
		
		// SKYnGM2008LU
		<<"CanStorParam"<<','<<"IntercepCoeff"<<','<<"ThroughFall"<<','<<"CanFieldCap"<<','
		<<"DrainCoeff"<<','<<"DrainExpPar"<<','<<"LandUseAlb"<<','<<"VegHeight"<<','
		<<"OptTransmCoeff "<<','<<"StomRes"<<','<<"VegFraction "<<','<<"LeafAI";

		if ( time == 0 )
			arcofs<<','<<"SoilID"<<','<<"LUseID"<<endl<<flush;
		else
			arcofs<<"\n";    
	}
	
	cn = ni.FirstP();
	while (ni.IsActive()) {

		//cout << cn->getHlevel() << "\n";
		arcofs<<cn->getID()<<','
		<<setprecision(4)<<cn->getZ()<<','
		<<setprecision(5)<<fabs(atan(cn->getFlowEdg()->getSlope()))<<','
		<<setprecision(5)<<cn->getContrArea()<<','
		<<setprecision(5)<<cn->getNwtNew()<<','
		<<setprecision(5)<<cn->getMuNew()<<','
		<<setprecision(5)<<cn->getMiNew()<<','
		<<setprecision(5)<<cn->getNfNew()<<','
		<<setprecision(5)<<cn->getNtNew()<<','
		<<setprecision(5)
		<<cn->getQpout()*1.E-6/cn->getVArea()<<','
		<<cn->getQpin()*1.E-6/cn->getVArea()<<','
		<<setprecision(4)<<cn->getCumSrf()<<','
		<<setprecision(3)<<cn->getRain()<<','

		// SKY2008Snow from AJR2007
		<<setprecision(3)<<cn->getIceWE()+cn->getLiqWE()<<','//added by AJR 2007 @ NMT
		<<setprecision(3)<<cn->getSnTempC()<<','//added by AJR 2007 @ NMT
		<<setprecision(5)<<cn->getIceWE()<<','//added by AJR 2007 @ NMT
		<<setprecision(5)<<cn->getLiqWE()<<','//added by AJR 2007 @ NMT
		<<setprecision(7)<<cn->getCumSnSub()<<','//added by CJC2020
		<<setprecision(7)<<cn->getCumSnEvap()<<','//added by CJC2020
		<<setprecision(7)<<cn->getCumMelt()<<','//changed to cumulative melt CJC2021 <<setprecision(5)<<cn->getDU()<<','//added by AJR 2007 @ NMT
		<<setprecision(7)<<cn->getCumHrsSnow()<<','//changed to cumulative hours snow CJC2021	<<setprecision(5)<<cn->getUnode()<<','//added by AJR 2007 @ NMT
		<<setprecision(5)<<cn->getSnLHF()<<','//added by AJR 2007 @ NMT
		<<setprecision(5)<<cn->getSnSHF()<<','//added by AJR 2007 @ NMT
		<<setprecision(5)<<cn->getSnGHF()<<','//added by AJR 2007 @ NMT
		<<setprecision(5)<<cn->getSnPHF()<<','//added by AJR 2007 @ NMT
		<<setprecision(5)<<cn->getSnRLout()<<','//added by AJR 2007 @ NMT
		<<setprecision(5)<<cn->getSnRLin()<<','//added by AJR 2007 @ NMT
		<<setprecision(5)<<cn->getSnRSin()<<','//added by AJR 2007 @ NMT
		<<setprecision(5)<<cn->getUerror()<<','//added by AJR 2007 @ NMT
		<<setprecision(5)<<cn->getIntSWE()<<','//added by AJR 2007 @ NMT 
		<<setprecision(5)<<cn->getCumIntSub()<<','//added by AJR 2007 @ NMT			
		<<setprecision(7)<<cn->getCumIntUnl()<<','//changed to cumulative unload CJC2021	<<setprecision(5)<<cn->getIntSnUnload()<<','//added by AJR 2007 @ NMT

		<<setprecision(3)<<cn->getSoilMoistureSC()<<','
		<<setprecision(3)<<cn->getRootMoistureSC()<<','
		<<setprecision(3)<<cn->getCanStorage()<<','
		<<setprecision(3)<<cn->getActEvap()<<','
		<<setprecision(5)<<cn->getCumBarEvap()<<',' // change to cumulative outputs CJC2020 TODO review cumulative output added by josh
		<<setprecision(5)<<cn->getCumTotEvap()<<',' // change to cumulative outputs CJC2020
		<<setprecision(3)<<cn->getGFlux()<<','
		<<setprecision(3)<<cn->getHFlux()<<','
		<<setprecision(3)<<cn->getLFlux()<<','
		<<setprecision(3)<<cn->getQstrm()<<','
		<<setprecision(3)<<cn->getHlevel()<<','
		//<<setprecision(3)<<cn->getFlowVelocity();
		<<setprecision(3)<<cn->getFlowVelocity()<<','// SKYnGM2008LU

		// SKYnGM2008LU
		<<setprecision(5)<<cn->getCanStorParam()<<','
		<<setprecision(5)<<cn->getIntercepCoeff()<<','
		<<setprecision(5)<<cn->getThroughFall()<<','
		<<setprecision(5)<<cn->getCanFieldCap()<<','
		<<setprecision(5)<<cn->getDrainCoeff()<<','
		<<setprecision(5)<<cn->getDrainExpPar()<<','
		<<setprecision(5)<<cn->getLandUseAlb()<<','
		<<setprecision(5)<<cn->getVegHeight()<<','
		<<setprecision(5)<<cn->getOptTransmCoeff()<<','
		<<setprecision(5)<<cn->getStomRes()<<','
		<<setprecision(5)<<cn->getVegFraction()<<','
		<<setprecision(5)<<cn->getLeafAI();

		if ( time == 0 )
			arcofs<<','<<setprecision(0)<<cn->getSoilID()<<','
				<<setprecision(0)<<cn->getLandUse()<<endl;
		else 
			arcofs<<"\n";
		cn = ni.NextP();
	}
	arcofs.close();
	
	// Call another output function only at the beginning
	// and end of simulation
	if ( time == 0.0 || this->timer->IsFinished() )
		WriteIntegrVars( time );

	// Write binary information that varies with time step for visualizer
	if (this->vizOption == 1)
		WriteDynamicVarsBinary(time);
	
	return;
}

template< class tSubNode >
void tCOutput<tSubNode>::WriteDynamicVarsBinary( double time )
{
	tSubNode* cn;
	tMeshListIter<tCNode> niter(this->g->getNodeList());
	int hour = (int)floor(time);

	char extension[20];
	sprintf(extension, "_dyn.%04d", hour);

	ofstream ostr;
	this->CreateAndOpenVizFile(&ostr, extension);

	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getNwtNew());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getNfNew());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getNtNew());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getMuNew());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getQpout());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getQpin());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getGwaterChng());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getSrf());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getRain());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getSoilMoistureSC());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getRootMoistureSC());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getAirTemp());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getDewTemp());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getSurfTemp());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getSoilTemp());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getAirPressure());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getRelHumid());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getSkyCover());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getWindSpeed());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getNetRad());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getActEvap());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getEvapoTrans());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getEvapSoil());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getGFlux());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getHFlux());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getLFlux());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getNetPrecipitation());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getRecharge());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getQstrm());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getHlevel());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getCanStorage());
	for (cn = niter.FirstP(); niter.IsActive(); cn = niter.NextP())
		BinaryWrite(ostr, (float) cn->getFlowVelocity());
}

/*************************************************************************
**
**  WriteDynamicVar( double time )
**
**  Writes files containing various dynamic variables for all the nodes
**  The option is turned on with the '-R' option in the command line
**  
*************************************************************************/
template< class tSubNode >
void tCOutput<tSubNode>::WriteDynamicVar( double time )
{
	if (this->g->getNodeList()->getActiveSize() < 0) {
		
		int hour, minute;
		char extension[20];
		
		tSubNode *cn;
		tMeshListIter<tSubNode> ni( this->g->getNodeList() );
		
		hour   = (int)floor(time);
		minute = (int)(time-hour)*100;
		sprintf(extension,"%04d.%02d", hour, minute);
		
		//for (int i = 0; i < 19; i++) {
		for (int i = 0; i < 36; i++) {  // SKY2008Snow from AJR2007
			this->dynvars[i]<<extension<<" ";
			if (i == 0)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
					this->dynvars[i]<<setprecision(5)<<cn->getNwtNew()<<" ";
			else if (i == 1)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
					this->dynvars[i]<<setprecision(5)<<cn->getNfNew()<<" ";
			else if (i == 2)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
					this->dynvars[i]<<setprecision(5)<<cn->getNtNew()<<" ";
			else if (i == 3)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
					this->dynvars[i]<<setprecision(5)<<cn->getMuNew()<<" ";
			else if (i == 4)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
					this->dynvars[i]<<setprecision(5)<<cn->getMiNew()<<" ";
			else if (i == 5)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
					this->dynvars[i]<<setprecision(8)<<cn->getQpout()*1.E-6/cn->getVArea()<<" ";
			else if (i == 6)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
					this->dynvars[i]<<setprecision(8)<<cn->getQpin()*1.E-6/cn->getVArea()<<" ";
			else if (i == 7)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
					this->dynvars[i]<<setprecision(7)<<cn->getRuNew()<<" ";
			else if (i == 8)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
					this->dynvars[i]<<setprecision(7)<<cn->getTransmiss()*1.E-6<<" ";
			else if (i == 9)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
					this->dynvars[i]<<setprecision(7)<<cn->getGwaterChng()*1E-9<<" ";
			else if (i == 10)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
					this->dynvars[i]<<setprecision(7)<<cn->getSrf()<<" ";
			else if (i == 11)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
					this->dynvars[i]<<setprecision(5)<<cn->getSoilMoistureSC()<<" ";
			else if (i == 12)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
					this->dynvars[i]<<setprecision(5)<<cn->getRootMoisture()<<" ";
			else if (i == 13)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
					this->dynvars[i]<<setprecision(5)<<cn->getEvapoTrans()<<" ";
			else if (i == 14)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
					this->dynvars[i]<<setprecision(5)<<cn->getEvapDryCanopy()<<" ";
			else if (i == 15)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
					this->dynvars[i]<<setprecision(5)<<cn->getEvapWetCanopy()<<" ";
			else if (i == 16)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
					this->dynvars[i]<<setprecision(5)<<cn->getEvapSoil()<<" ";
			else if (i == 17)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
					this->dynvars[i]<<setprecision(5)<<cn->getNetPrecipitation();
			else if (i == 18)

			// SKY2008Snow from AJR2007
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
	  				this->dynvars[i]<<setprecision(5)<<cn->getIceWE()+cn->getLiqWE();//added by AJR 2007 @ NMT
      			else if (i == 19)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
	  				this->dynvars[i]<<setprecision(5)<<cn->getSnTempC();//added by AJR 2007 @ NMT
      			else if (i == 20)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
	  				this->dynvars[i]<<setprecision(5)<<cn->getIceWE();//added by AJR 2007 @ NMT
      			else if (i == 21)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
	  				this->dynvars[i]<<setprecision(5)<<cn->getLiqWE();//added by AJR 2007 @ NMT
      			else if (i == 22)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
	  				this->dynvars[i]<<setprecision(5)<<cn->getDU();//added by AJR 2007 @ NMT
      			else if (i == 23)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
	  				this->dynvars[i]<<setprecision(5)<<cn->getUnode();//added by AJR 2007 @ NMT
      			else if (i == 24)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
	  				this->dynvars[i]<<setprecision(5)<<cn->getSnLHF();//added by AJR 2007 @ NMT
      			else if (i == 25)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
	  				this->dynvars[i]<<setprecision(5)<<cn->getSnSHF();//added by AJR 2007 @ NMT
      			else if (i == 26)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
	  				this->dynvars[i]<<setprecision(5)<<cn->getSnGHF();//added by AJR 2007 @ NMT
      			else if (i == 27)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
					this->dynvars[i]<<setprecision(5)<<cn->getSnPHF();//added by AJR 2007 @ NMT
			else if (i == 28)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
	  				this->dynvars[i]<<setprecision(5)<<cn->getSnRLout();//added by AJR 2007 @ NMT
			else if (i == 29)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
					this->dynvars[i]<<setprecision(5)<<cn->getSnRLin();//added by AJR 2007 @ NMT
			else if (i == 30)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
					this->dynvars[i]<<setprecision(5)<<cn->getSnRSin();//added by AJR 2007 @ NMT
			else if (i == 31)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
					this->dynvars[i]<<setprecision(5)<<cn->getUerror();//added by AJR 2007 @ NMT
			else if (i == 32)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
					this->dynvars[i]<<setprecision(5)<<cn->getIntSWE();//added by AJR 2007 @ NMT
			else if (i == 33)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
					this->dynvars[i]<<setprecision(5)<<cn->getIntSub();//added by AJR 2007 @ NMT
			else if (i == 34)
				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
					this->dynvars[i]<<setprecision(5)<<cn->getIntSnUnload();//added by AJR 2007 @ NMT
			else if (i == 35)


				for ( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
					this->dynvars[i]<<setprecision(5)<<cn->getRiNew()<<" ";
			this->dynvars[i]<<"\n";
		}
	}
	return;
}

/*************************************************************************
**
**  WriteIntegrVars( double time )
**
**  Writes a file containing variables that contain integrated variables
**  The option is turned on with the '-R' option in the command line
**  
*************************************************************************/
template< class tSubNode >
void tCOutput<tSubNode>::WriteIntegrVars( double time )
{
	int hour, minute;
	int Occur, prec;
	char extension[20];
	double avRate, tmp1, tmp2;
	tSubNode *cn;
	tMeshListIter<tSubNode> ni( this->g->getNodeList() );
	
	hour   = (int)floor(time);
	minute = (int)floor((time-hour)*60);
	
	sprintf(extension,".%04d_%02di", hour, minute);
	this->CreateAndOpenFile(&intofs, extension);
	
	if (simCtrl->Header_label == 'Y') {
		intofs<<"ID"<<','<<"BndCd"<<','<<"Z"<<','<<"VAr"<<','<<"CAr"<<','<<"Curv"
		<<','<<"EdgL"<<','<<"Slp"
		<<','<<"FWidth"<<','<<"Aspect"

		// SKY2008Snow from AJR2007
		<<','<<"SV"<<','<<"LV"//added by AJR 2007 @ NMT

		<<','<<"AvSM"<<','<<"AvRtM"
		<<','<<"HOccr"<<','<<"HRt"
		<<','<<"SbOccr"<<','<<"SbRt"
		<<','<<"POccr"<<','<<"PRt"
		<<','<<"SatOccr"<<','<<"SatRt"
		<<','<<"SoiSatOccr"<<','<<"RchDsch"
		<<','<<"AvET"<<','<<"EvpFrct"

		// SKY2008Snow from AJR2007
		<<','<<"cLHF"<<','<<"cMelt"//added by AJR 2007 @ NMT
		<<','<<"cSHF"<<','<<"cPHF"//added by AJR 2007 @ NMT
		<<','<<"cRLIn"<<','<<"cRLo"//added by AJR 2007 @ NMT
		<<','<<"cRSIn"<<','<<"cGHF"//added by AJR 2007 @ NMT
		<<','<<"cUErr"<<','<<"cHrsSun"<<','<<"cHrsSnow"//added by AJR 2007 @ NMT
		<<','<<"persTime"<<','<<"peakWE"<<','<<"initTime"<<','<<"peakTime"
		<<','<<"cIntSub"<<','<<"cSnSub"<<','<<"cSnEvap"<<','<<"cIntUnl"//added by AJR 2007 @ NMT 

		// SKYnGM2008LU
		<<','<<"AvCanStorParam"<<','<<"AvIntercCoeff"
		<<','<<"AvTF"<<','<<"AvCanFieldCap"
		<<','<<"AvDrainCoeff"<<','<<"AvDrainExpPar"
		<<','<<"AvLUAlb"<<','<<"AvVegHeight"
		<<','<<"AvOTCoeff"<<','<<"AvStomRes"
		<<','<<"AvVegFract"<<','<<"AvLeafAI"

		<<"\n";
	}
	
	cn = ni.FirstP();
	while (ni.IsActive()) {
		
		intofs<<cn->getID()<<','<<cn->getBoundaryFlag()<<','
		<<setprecision(4)<<cn->getZ()<<','
		<<setprecision(7)<<cn->getVArea()<<','
		<<setprecision(7)<<cn->getContrArea()*1.E-6<<','
		<<setprecision(6)<<cn->getCurvature()<<','
		<<cn->getFlowEdg()->getLength()<<','
		<<cn->getFlowEdg()->getSlope()<<','
		<<cn->getFlowEdg()->getVEdgLen()<<','
		<<setprecision(4)<<cn->getAspect()<<',' //;

		// SKY2008Snow from AJR2007
		<<setprecision(7)<<cn->getSheltFact()<<','//added by AJR 2007 @ NMT
		<<cn->getLandFact()<<','//added by AJR 2007 @ NMT
		<<setprecision(4);
		
		tmp1 = floor(cn->getAvSoilMoisture())*1.E-4; 
		tmp2 = (cn->getAvSoilMoisture()-floor(cn->getAvSoilMoisture()))*1.E+1;
		intofs<<tmp1<<','<<tmp2<<',';
		
		// -----------
		Occur = (int)((cn->hsrfOccur-floor(cn->hsrfOccur))*1.E+6);
		if (Occur > 0) {
			avRate = floor(cn->hsrfOccur)/1000.0/Occur;
			prec = 3;
			if (avRate>100.0) 
				prec++;
		}
		else {
			avRate = 0.0;
			prec = 0;
		}
		intofs<<setprecision(6)<<Occur<<setprecision(prec)<<','<<avRate<<',';
		
		// -----------
		Occur = (int)((cn->sbsrfOccur-floor(cn->sbsrfOccur))*1.E+6);
		if (Occur > 0) {
			avRate = floor(cn->sbsrfOccur)/1000.0/Occur;
			prec = 3;
			if (avRate>100.) 
				prec++;
		}
		else {
			avRate = 0.0;
			prec = 0;
		}
		intofs<<setprecision(6)<<Occur<<setprecision(prec)<<','<<avRate<<',';
		
		// -----------
		Occur = (int)((cn->psrfOccur-floor(cn->psrfOccur))*1.E+6);
		if (Occur > 0) {
			avRate = floor(cn->psrfOccur)/1000./Occur;
			prec = 3;
			if (avRate>100.) 
				prec++;
		}
		else {
			avRate = 0.0;
			prec = 0;
		}
		intofs<<setprecision(6)<<Occur<<setprecision(prec)<<','<<avRate<<',';
		
		// -----------
		Occur = (int)((cn->satsrfOccur-floor(cn->satsrfOccur))*1.E+6);
		if (Occur > 0) {
			avRate = floor(cn->satsrfOccur)/1000./Occur;
			prec = 3;
			if (avRate>100.) 
				prec++;
		}
		else {
			avRate = 0.0;
			prec = 0;
		}
		intofs<<setprecision(6)<<Occur<<','
			<<setprecision(prec)<<avRate<<','
			<<setprecision(6)<<cn->satOccur<<','
			<<setprecision(3)<<cn->RechDisch<<','
			<<setprecision(4)<<cn->getAvET()<<','

			// SKY2008Snow from AJR2007
			//<<setprecision(4)<<cn->getAvEvapFract();
			<<setprecision(4)<<cn->getAvEvapFract()<<','
			<<setprecision(7)<<cn->getCumLHF()<<','//added by AJR 2007 @ NMT
			<<setprecision(7)<<cn->getCumMelt()<<','//added by AJR 2007 @ NMT
			<<setprecision(7)<<cn->getCumSHF()<<','//added by AJR 2007 @ NMT
			<<setprecision(7)<<cn->getCumPHF()<<','//added by AJR 2007 @ NMT
			<<setprecision(7)<<cn->getCumRLin()<<','//added by AJR 2007 @ NMT -- corrected from getCumRLout call, SKY2008Snow 
			<<setprecision(7)<<cn->getCumRLout()<<','//added by AJR 2007 @ NMT -- corrected from getCumRLin call, SKY2008Snow
			<<setprecision(7)<<cn->getCumRSin()<<','//added by AJR 2007 @ NMT
			<<setprecision(7)<<cn->getCumGHF()<<','//added by AJR 2007 @ NMT
			<<setprecision(7)<<cn->getCumUerror()<<','//added by AJR 2007 @ NMT
			<<setprecision(7)<<cn->getCumHrsSun()<<','//added by AJR 2007 @ NMT
			<<setprecision(7)<<cn->getCumHrsSnow()<<','//added by AJR 2007 @ NMT
			<<setprecision(7)<<cn->getPersTimeMax()<<',' //added by AJR 2007 @ NMT
			<<setprecision(7)<<cn->getPeakSWE()<<',' //added by AJR 2007 @ NMT
			<<setprecision(7)<<cn->getInitPackTime()<<',' //added by AJR 2007 @ NMT
			<<setprecision(7)<<cn->getPeakPackTime()<<',' //added by AJR 2007 @ NMT
			<<setprecision(7)<<cn->getCumIntSub()<<','//added by AJR 2007 @ NMT
			<<setprecision(7)<<cn->getCumSnSub()<<','//added by CJC2020
			<<setprecision(7)<<cn->getCumSnEvap()<<','//added by CJC2020			
			<<setprecision(7)<<cn->getCumIntUnl()<<','//added by AJR 2007 @ NMT

			// SKYnGM2008LU
			<<setprecision(7)<<cn->getAvCanStorParam()<<','
			<<setprecision(7)<<cn->getAvIntercepCoeff()<<','
			<<setprecision(7)<<cn->getAvThroughFall()<<','
			<<setprecision(7)<<cn->getAvCanFieldCap()<<','
			<<setprecision(7)<<cn->getAvDrainCoeff()<<','
			<<setprecision(7)<<cn->getAvDrainExpPar()<<','
			<<setprecision(7)<<cn->getAvLandUseAlb()<<',' 
			<<setprecision(7)<<cn->getAvVegHeight()<<',' 
			<<setprecision(7)<<cn->getAvOptTransmCoeff()<<',' 
			<<setprecision(7)<<cn->getAvStomRes()<<',' 
			<<setprecision(7)<<cn->getAvVegFraction()<<','
			<<setprecision(7)<<cn->getAvLeafAI();			

		intofs<<setprecision(6)<<"\n";
		
		cn = ni.NextP();
	}
	intofs.close();
	return;
}

/*************************************************************************
**
**  WriteOutletInfo( double time )
**
**  Writes a file containing simulated streamflow and stage values
**  The output format should be readable by ArcInfo & Matlab 
**  
*************************************************************************/
template< class tSubNode >
void tCOutput<tSubNode>::WriteOutletInfo( double time )
{
	if (numOutlets > 0) {
		int hour, minute;
		char extension[20];
		
		hour   = (int)floor(time); 
		minute = (int)((time-hour)*100);
		sprintf(extension,"%04d.%02d", hour, minute);
		
		for (int i = 0; i < numOutlets; i++) {
#ifdef PARALLEL_TRIBS
      // Node ID does not have to be less than the active size
      if ( (Outlets[i] != NULL) && (OutletList[i] > 0) ) {
#else
			if ( Outlets[i] && OutletList[i] < this->g->getNodeList()->getActiveSize()) {
#endif
				outletinfo[i]<<extension<<"\t"
				<<setw(10)<<Outlets[i]->getQstrm()<<"\t"
				<<setw(6) <<Outlets[i]->getHlevel()<<"\n";
			}
		}
	}
	return;
}

/*************************************************************************
**
**  tCOutput::ReadOutletList()
**
**  Opens and Reads the node list from a *.oul file whose structure is:
**
**  Number of Outlet Nodes
**  NodeID1 NodeID2 NodeID3 NodeID4 NodeID5 ...
**
*************************************************************************/
template< class tSubNode >
void tCOutput<tSubNode>::ReadOutletNodeList(char *nodeFileO)
{
	ifstream readOUL(nodeFileO);
	if (!readOUL) {
		Cout<<"\n>>>Outlet Node File "<<nodeFileO<<" not found..."<<endl;
		Cout<<">>>No output for interior nodes will be written..."<<endl<<endl;
		numOutlets = 0;
		return;
	}
	
	readOUL>>numOutlets;
	OutletList  = new int[numOutlets];
	Outlets     = new tSubNode*[numOutlets];
	outletinfo  = new ofstream[numOutlets];
	for (int i = 0; i < numOutlets; i++)
		readOUL>>OutletList[i]; 
	
#ifdef PARALLEL_TRIBS
  // Initialize Outlets to NULL, used to determine local Outlets
  for (int i = 0; i < numOutlets; i++)
    Outlets[i] = NULL;
#endif

	readOUL.close();
	return;
}

/*************************************************************************
**
**  tOutput::CreateAndOpenOutlet()
**
**  Write the header for the *.iout output file.
** 
*************************************************************************/
template< class tSubNode >
void tCOutput<tSubNode>::CreateAndOpenOutlet()
{
	if (numOutlets > 0)
	{
		if ( OutletList ) {
			char pixelext[10] = ".qout";
			char nodeNum[10], pixelnode[100];
			char fullName[kMaxNameSize+6];
		
			for (int i = 0; i < numOutlets; i++) {
#ifdef PARALLEL_TRIBS
       // Check if outlet node is on this processor
       if ( (Outlets[i] != NULL) && (OutletList[i] > 0) ) {
#else
				if (OutletList[i] > 0) {
#endif
					sprintf(nodeNum,"_%d", OutletList[i]);
					strcpy(pixelnode, nodeNum);
					strcat(pixelnode, pixelext);
				
					strcpy( fullName, outletName );
					strcat( fullName, pixelnode );
				
					outletinfo[i].open( fullName );
				
					if ( !outletinfo[i].good() )
						cerr<<"File "<<fullName<<"can not be created.";
					else
/*SMM
						Cout<<"Creating Output File: \t '"<<fullName<<"' "<<endl;
*/
					
					if (simCtrl->Header_label == 'Y') {
						outletinfo[i]<<"1-Time,hr\t "
						<<"2-Qstrm,m3/s\t"
						<<"3-Hlev,m"
						<<"\n";
					}
				}
			}
		}
	}
	return;
}

/*************************************************************************
**
**  tOutput::SetInteriorOutlet()
**
**  Initializes pointers to basin interior outlets 
** 
*************************************************************************/
template< class tSubNode >
void tCOutput<tSubNode>::SetInteriorOutlet()
{
	tSubNode * cnn;
	tMeshListIter<tSubNode> niter( this->g->getNodeList() );
	if ( OutletList ) {
		for (int i=0; i < numOutlets; i++) {
			if (OutletList[i] > 0) {
#ifdef PARALLEL_TRIBS
           // Each processor creates/writes outlet files for its points
           for (cnn=niter.FirstP(); niter.IsActive(); cnn=niter.NextP() ) {
#else
				for (cnn=niter.FirstP(); !(niter.AtEnd()); cnn=niter.NextP() ) {
#endif
					if (cnn->getID() == OutletList[i]) {
						Outlets[i] = cnn;   //Defining node of interest 
						Cout<<"\nOutlet of Interest ID: \t"<<OutletList[i]
							<<" has been set up..."<<endl<<flush;
					}
				}
			}
		}
	}
	return;
}

//=========================================================================
//
//
//                  Section 5: Update for New Run
//
//
//=========================================================================


/*************************************************************************
**
**  tCOutput::UpdateForNewRun( tInputFile &infile )
**
**  Used to update data members of tOutput when a new simulation run
**  is to be carried out (option -ON of the command line)
** 
**  Input: infile -- reference to an open input file, assumed valid
**
*************************************************************************/
template< class tSubNode >
void tCOutput<tSubNode>::UpdateForNewRun( tInputFile &infile )
{
	int j;
	char nodeFileO[kMaxNameSize];
	
	infile.ReadItem( this->baseName, "OUTFILENAME" ); 
	infile.ReadItem( outletName, "OUTHYDROFILENAME" );
	infile.ReadItem( this->nodeFile, "NODEOUTPUTLIST" );  
	infile.ReadItem( nodeFileO, "OUTLETNODELIST" );
	
	cout<<"\ntOutput basename: \t"<<this->baseName<<endl<<endl;
	
	this->nodeofs.close();
	this->edgofs.close();
	this->triofs.close();
	this->zofs.close();
	
	if (this->numNodes > 0) {
		for (j=0; j < this->numNodes; j++)
			this->pixinfo[j].close();
		delete [] this->nodeList;
		delete [] this->uzel;
		delete [] this->pixinfo;
	}
	this->numNodes = 0;
	
	if (numOutlets > 0) {
		for (j=0; j < numOutlets; j++)
			outletinfo[j].close();
		delete [] OutletList;
		delete [] Outlets;
		delete [] outletinfo;
	}
	numOutlets = 0;
	
	// Open pixel files
	this->ReadNodeOutputList(); 
	this->CreateAndOpenPixel();
	this->SetInteriorNode();
	
	// Open outlet files
	ReadOutletNodeList(nodeFileO);
	CreateAndOpenOutlet();
	SetInteriorOutlet();
	
	return;
}

//=========================================================================
//
//
//                       End of tOutput.cpp
//
//
//=========================================================================
