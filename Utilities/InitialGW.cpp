/***************************************************************************
**
**                           tRIBS Version 1.0
**
**              TIN-based Real-time Integrated Basin Simulator
**                       Ralph M. Parsons Laboratory
**                  Massachusetts Institute of Technology
**  
**
**                          Beta Release, 9/2001
**
**
**  InitialGW.cpp:   Utility program for tRIBS used for creating the initial
**                   groundwater elevation input to the model. It calculates
**                   the initial water table position in the basin area 
**                   based on an input steady state assumption.
**
**  Program compiled separately from tRIBS as:
**
**  UNIX% CC -o <executable> InitialGW.cpp
** 
***************************************************************************/


//=========================================================================
//
//
//                  Section 1: Include, Definitions
//                             Variable and Function Declarations
//
//
//=========================================================================


/***************************************************************************                        
** POINTERS USED HERE     -------------    ARC FLOW DIR ----------------- 
**                        I 4 I 3 I 2 I                 I 32 I 64 I 128 I           
**  0 = POINTS TO SELF    -------------                 ----------------- 
**      I.E. UNRESOLVED   I 5 I 0 I 1 I                 I 16 I  0 I   1 I
**  -1 = BOUNDARY PIXEL   -------------                 ----------------- 
**                        I 6 I 7 I 8 I                 I  8 I  4 I   2 I   
**                        -------------                 ----------------- 
***************************************************************************/

#include<iostream.h>
#include<iomanip.h>
#include<fstream.h>
#include<stdlib.h>
#include<string.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>

#define FACTOR 1 // Change factor used in division of GWT 

double AREA;    // km^2 - calculated value of basin area
double BASEF;   // BASEFLOW (m^3/s) - define steady-state profile
double BASEF_0; // BASEFLOW_0 (m^3/s) analysis of the recession
int    POROG;   // # of pixels used as a threshold for defining channels

double SAT, RESID, K0z,  F, POREIND;  //Soil parameters
double AR, UAR, PSIB, POROS, KS, CS;

int N = 0, M = 0;      // DIMENSIONS of the original array
double MIN = 99999;    // Minimum elevation of the basin is OUTLET
int MA = 3;            // Value by default for Moving Average  
int SMOO;              // Number of times smoothing applied
int DUMM;              // ArcView/Info dummy value for asc grids
double Zav, dx, dy, ddx, ddy, Zava, Zavb;
double Lambda = 0.0;    // Areal integral value of the topographic index
double Gamma = 0.0;     // Areal integral value of 

static int pi[]={0,0,-1,-1,-1,0,1,1,1};  // pointer for direction 
static int pj[]={0,1,1,0,-1,-1,-1,0,1};  // x pointer for direction
static int pi5[]={0,0,-1,-1,-1,0,1,1,1,0,-1,-2,-2,-2,-2,-2,-1,0,1,2,2,2,2,2,1};
static int pj5[]={0,1,1,0,-1,-1,-1,0,1,2,2,2,1,0,-1,-2,-2,-2,-2,-2,-1,0,1,2,2};

double get_cellEffectiveWidth(int);
double get_cellD8Width(int);
double getMoisture(double);
void smoothWTrelief(double**, double**, double**, double*);
void assignToZero(double**, int**, double**);
void zaglubiStream(double**, int**, double**, double);
void getWTmapSiva(double**, int**, double**, double**, double*, double*, int*);
void AdjustNegativeSlope(double **, int **, double *, int, int, double *);


//=========================================================================
//
//
//                  Section 2: Main Routine
//
//
//=========================================================================


int main(int argc, char *argv[]){
	
	char lineIn[300], head[20], namos[50];
	long int cnt = 0, cnt11 = 0;
	int nrows, ncols;
	int i, j, ii, jj, io, jo, dir;
	double slope, effwidth, Nwt;
	double summ = 0, summP = 0;
	double xllcorner, yllcorner;
	double cellsize, nodata;
	char *names[5];
	
	cout<<"\n-----------------------------------------------------------------"
		<<"-------";
	cout<<"\n\n\t\t tRIBS -- Version 1.0";
	cout<<"\n\t\t tRIBS Model: InitialGW Utility";
	cout<<"\n\t\t Ralph M. Parsons Laboratory";
	cout<<"\n\t\t Massachusetts Institute of Technology";
	cout<<"\n\n\t\t Release, 9/2001 \n\n";
	cout<<"-----------------------------------------------------------------"
		<<"-------"<<endl;
	
	cout<<"\nThis program implements a TOPMODEL approach of defining a steady-state GW\n"
		<<"surface based on the contributing area distribution in the basin and the\n"
		<<"local slope of the terrain surface."<<endl;
	
	// Command line arguments
	if(argc != 2) {
		cout<<"\nUsage: "<< argv[0] << "  *.gw\n" << endl;
		cout<<"Order of names in the file: " << endl;
		cout<<"\t 1. DEM file (*.asc)\n";
		cout<<"\t 2. Flow Accumulations file (*.asc)\n";
		cout<<"\t 3. Flow Directions file (*.asc)\n";
		cout<<"\t 4. Soil reclassification table (*.sdtt)\n";
		cout<<"\t 5. Basin area (km^2)\n";
		cout<<"\t 6. Baseflow of Reference Q(0) at outlet (m^3/sec)\n";
		cout<<"\t 7. Zero Baseflow Q_0 from recession analysis (m^3/sec)\n";
		cout<<"\t 8. Threshold value for stream network (# of pixels)\n";
		cout<<"\t 9. Size of MA smoothing window (3 or 5 pixels)\n";
		cout<<"\t 10. Number of smoothing iterations (2 or greater)\n";
		cout<<"\t 11. Average basin GWT depth (~ -ln(Q(0)/Q_0)/f)\n\n";
		exit(1);
	}
	
	ifstream Inp0(argv[1]);
	if (!Inp0) {
		cout<< "File "<<argv[1]<<" not found!"<<endl;
		exit(2);
	}
	
	cout<<"Input Parameters: \n"<<endl;
	for (i = 0; i < 4; i++) {
		Inp0 >> namos;
		names[i] = new char [strlen(namos) + 1];
		strcpy(names[i], namos);
		if(i == 0)
			cout<<"DEM Grid: \t\t\t"<<names[i]<<endl;
		else if(i==1)
			cout<<"Flow Accumulations Grid: \t"<<names[i]<<endl;
		else if(i==2)
			cout<<"Flow Directions Grid: \t\t"<<names[i]<<endl;
		else if(i==3)
			cout<<"Soil Reclassification Table: \t"<<names[i]<<endl;
	}
	
	//Read in parameters
	Inp0 >> AREA;
	Inp0 >> BASEF;
	Inp0 >> BASEF_0;
	Inp0 >> POROG;
	Inp0 >> MA;
	Inp0 >> SMOO;
	Inp0 >> Zav;
	
	cout<<"Basin Area: \t\t\t"<< AREA << endl;
	cout<<"Baseflow: \t\t\t" << BASEF << endl;
	cout<<"Baseflow Zero: \t\t\t" << BASEF_0 << endl;
	cout<<"Threshold pixels: \t\t" << POROG << endl;
	cout<<"Moving Average Window: \t\t" << MA << endl;
	cout<<"Average Depth to WT: \t\t"<<Zav<<endl;
	
	MA = MA * MA;
	
	//Read in the file names
	ifstream Inp1(names[0]);
	if (!Inp1) {
		cout << "\nFile '" << names[0] << "' not found!" << endl;
		cout << "Exiting Program..."<<endl;
		exit(2);
	}
	
	ifstream Inp2(names[1]);
	if (!Inp2) {
		cout << "\nFile '" << names[1] << "' not found!" << endl;
		cout << "Exiting Program..."<<endl;
		exit(2);
	}
	
	ifstream Inp3(names[2]);
	if (!Inp3) {
		cout << "\nFile '" << names[2] << "' not found!" << endl;
		cout << "Exiting Program..."<<endl;
		exit(2);
	}
	
	ifstream Inp4(names[3]);
	if (!Inp3) {
		cout << "\nFile '" << names[3] << "' not found!" << endl;
		cout << "Exiting Program..."<<endl;
		exit(2);
	}
	
	//Create Output files
	ofstream Otp2("_hillslope.hist");
	ofstream Otp3("_MeanStd.hist");
	ofstream Otp4("_occurence.asc");
	ofstream Otp5("_GWTabs.asc");
	ofstream Otp6("_GWTabs_smooth.asc");
	ofstream Otp ("_GWTdepth.asc");
	ofstream Otp7("_GWTdepth_smooth.asc");
	ofstream Otp8("_depth_instream.asc"); 
	
	//Read DEM File
	cout<<"\nReading DEM file header..."<<endl;
	
	Inp1 >> head >> ncols;
	Inp1 >> head >> nrows;
	Inp1 >> head >> xllcorner;
	Inp1 >> head >> yllcorner;
	Inp1 >> head >> cellsize;
	Inp1 >> head >> nodata;
	DUMM = (int)nodata;
	
	//Write header information to each output ASCII grid
	Otp << "ncols         "<<ncols<<endl;
	Otp << "nrows         "<<nrows<<endl;
	Otp << "xllcorner     "<<setiosflags(ios::fixed)<<xllcorner<<endl;
	Otp << "yllcorner     "<<setiosflags(ios::fixed)<<yllcorner<<endl;
	Otp << "cellsize      "<<resetiosflags(ios::fixed)<<cellsize<<endl;
	Otp << "NODATA_value  "<<nodata<<endl;
	
	Otp4 << "ncols        "<<ncols<<endl;
	Otp4 << "nrows        "<<nrows<<endl;
	Otp4 << "xllcorner    "<<setiosflags(ios::fixed)<<xllcorner<<endl;
	Otp4 << "yllcorner    "<<setiosflags(ios::fixed)<<yllcorner<<endl;
	Otp4 << "cellsize     "<<resetiosflags(ios::fixed)<<cellsize<<endl;
	Otp4 << "NODATA_value "<<nodata<<endl;
	
	Otp5 << "ncols        "<<ncols<<endl;
	Otp5 << "nrows        "<<nrows<<endl;
	Otp5 << "xllcorner    "<<setiosflags(ios::fixed)<<xllcorner<<endl;
	Otp5 << "yllcorner    "<<setiosflags(ios::fixed)<<yllcorner<<endl;
	Otp5 << "cellsize     "<<resetiosflags(ios::fixed)<<cellsize<<endl;
	Otp5 << "NODATA_value "<<nodata<<endl;
	
	Otp6 << "ncols        "<<ncols<<endl;
	Otp6 << "nrows        "<<nrows<<endl;
	Otp6 << "xllcorner    "<<setiosflags(ios::fixed)<<xllcorner<<endl;
	Otp6 << "yllcorner    "<<setiosflags(ios::fixed)<<yllcorner<<endl;
	Otp6 << "cellsize     "<<resetiosflags(ios::fixed)<<cellsize<<endl;
	Otp6 << "NODATA_value "<<nodata<<endl;
	
	Otp7 << "ncols        "<<ncols<<endl;
	Otp7 << "nrows        "<<nrows<<endl;
	Otp7 << "xllcorner    "<<setiosflags(ios::fixed)<<xllcorner<<endl;
	Otp7 << "yllcorner    "<<setiosflags(ios::fixed)<<yllcorner<<endl;
	Otp7 << "cellsize     "<<resetiosflags(ios::fixed)<<cellsize<<endl;
	Otp7 << "NODATA_value "<<nodata<<endl;
	
	Otp8 << "ncols        "<<ncols<<endl;
	Otp8 << "nrows        "<<nrows<<endl;
	Otp8 << "xllcorner    "<<setiosflags(ios::fixed)<<xllcorner<<endl;
	Otp8 << "yllcorner    "<<setiosflags(ios::fixed)<<yllcorner<<endl;
	Otp8 << "cellsize     "<<resetiosflags(ios::fixed)<<cellsize<<endl;
	Otp8 << "NODATA_value "<<nodata<<endl;
	
	N = nrows;
	M = ncols;
	dx = cellsize;
	dy = cellsize;
	
	//Read Flow Accumulations file
	cout<<"Reading Flow Accumulations file header..."<<endl;
	
	Inp2 >> head >> ncols;
	Inp2 >> head >> nrows;
	Inp2 >> head >> xllcorner;
	Inp2 >> head >> yllcorner;
	Inp2 >> head >> cellsize;
	Inp2 >> head >> nodata;
	
	//Read Flow Directions file
	cout<<"Reading Flow Directions file header..."<<endl;
	
	Inp3 >> head >> ncols;
	Inp3 >> head >> nrows;
	Inp3 >> head >> xllcorner;
	Inp3 >> head >> yllcorner;
	Inp3 >> head >> cellsize;
	Inp3 >> head >> nodata;
	
	//Read in soil parameters and output to screen
	cout<<"Reading Basin Averaged Soil Properties..."<<endl;
	int NS, NM, ID;
	
	Inp4 >> NS >> NM;
	Inp4 >> ID >> K0z >> SAT >> RESID >> POREIND >> PSIB;
	Inp4 >> F >> AR >> UAR >> POROS >> KS >> CS;
	
	cout <<"\nSoil types: \t\t"<<NS<<"\nSoil properties: \t"<<NM<<"\nID: \t\t\t"<<ID;
	cout <<"\nKOz: \t\t\t"<<K0z<<"\nSAT: \t\t\t"<<SAT<<"\nRESID: \t\t\t"<<RESID;
	cout <<"\nPOREIND: \t\t"<<POREIND<<"\nPSIB: \t\t\t"<<PSIB<<"\nF: \t\t\t"<<F;
	cout <<"\nAR: \t\t\t"<<AR<<"\nUAR: \t\t\t"<<UAR<<"\nPOROS: \t\t\t"<<POROS;
	cout <<"\nKS: \t\t\t"<<KS<<"\nCS: \t\t\t"<<CS;
	
	//Basin information
	AREA = AREA * 1.0e+6/(dx*dy);   // Total area (pixels)
	BASEF = BASEF*3.6*1.0e+12;      // Transform into mm/hr
	
	//Moving Average Windows
	//Window Size = 3
	double dist[9];
	for(i=1; i < 9; i++) {
		ddx = dx*(double)pi[i];  //dx = vertical dy = horizontal
		ddy = dy*(double)pj[i];
		dist[i] = sqrt(ddx*ddx + ddy*ddy);
	}
	
	//Window Size = 5
	double dist5[25];
	for(i=1; i < 25; i++) {
		ddx = dx*(double)pi5[i]; 
		ddy = dy*(double)pj5[i];
		dist5[i] = sqrt(ddx*ddx + ddy*ddy);
	}
	
	//Memory Allocation
	int    **area, **ptr;
	double *means, *stds;
	double **dem, **GWT, **GWTabs, **aIndex, **topoi;
	int *count, *ci; 
	
	area = new int* [N];  
	assert(area != 0);
	
	dem = new double* [N];
	assert(dem != 0);
	
	ptr = new int* [N]; 
	assert(ptr != 0);
	
	GWT = new double* [N]; 
	assert(GWT != 0); 
	
	GWTabs = new double* [N]; 
	assert(GWTabs!=0);
	
	aIndex = new double* [N];  //Stores topographic index
	assert(aIndex!=0);  
	
	topoi = new double* [POROG]; //Count slopes for areas<threshold
	assert(topoi!=0); 
	
	count = new int [POROG];   //Count pixels having area<threshold = i
	assert(count!=0);
	
	ci = new int [POROG];     //Count indices in the array topoi
	assert(ci!=0);   
	
	means = new double [POROG];  
	assert(means != 0);
	
	stds = new double [POROG]; 
	assert(stds != 0);
	
	for(i=0; i < POROG; i++) { 
		count[i] = 0;
		ci[i] = 0;
		means[i] = 0;
		stds[i] = 0;
	}
	
	for (i=0; i < N; i++){
		area[i] = new int[M];    
		assert(area[i] != 0);
		
		dem[i] = new double[M]; 
		assert(dem[i] != 0);
		
		ptr[i] = new int[M];     
		assert(ptr[i] != 0);
		
		GWT[i] = new double[M];  
		assert(GWT[i] != 0);
		
		GWTabs[i] = new double[M];  
		assert(GWTabs[i] != 0);
		
		aIndex[i] = new double[M];  
		assert(aIndex[i] != 0);
	}
	
	
	//  ------------------- DATA INPUT -------------------
	//  -------- DATA CHECK AND INITIAL SUMMATION --------
	
	//Data Input from ASCII grids
	cout<<"\n\nReading in grid data..."<<endl;
	for (i=0; i < N; i++) {
		for(j=0; j < M; j++) {
			Inp1 >> dem[i][j];
			Inp2 >> area[i][j];
			Inp3 >> dir;
			
			//Reassign direction values from Arc/Info to the values in range 1-8
			
			if (dir == 128)
				dir = 2;
			else if (dir == 64)
				dir = 3;
			else if (dir == 32)
				dir = 4;
			else if (dir == 16)
				dir = 5;
			else if (dir == 8)
				dir = 6;
			else if (dir == 4)
				dir = 7;
			else if (dir == 2)
				dir = 8;
			else if (dir < 0)
				dir = DUMM;
			
			ptr[i][j] = dir;
			
			//Modify Flow Accumulations to Fit 
			//Change from A = 0 (upstream) to A = 1
			
			if(area[i][j] < 0)
				area[i][j] = DUMM;
			else
				area[i][j] = area[i][j] + 1;
			
			
			// --- Checking for valid data in DEM file --- 
			// --- Creating histogram of hillslope cells --- 
			if (dem[i][j] != DUMM) {
				cnt++;
				if (area[i][j] <= POROG)   // Look at the pixels with contr.area<POROG
					count[area[i][j]-1]++;
			}
			// --- Finding the DEM outlet called MIN ---
			if (dem[i][j] < MIN && dem[i][j] > 0) {
				MIN = dem[i][j];
				io = i;
				jo = j;
			}
			
			// --- Cross-checking validity of input files ---
			// --- DEM and Flow Accumulations have equal valid pixels ---
			if((area[i][j]!=DUMM && dem[i][j]==DUMM) || (area[i][j]==DUMM && dem[i][j]!=DUMM)){
				cout<<"\nWarning: In row " << i <<" the Flow Accum mismatched with DEM"<<endl;
				cout<<"DEM value= "<<dem[i][j]<<"\tACCM value = "<<area[i][j]<<"\n";
				cout<<"ROW: "<< i <<"\tCOLUMN: "<<j<<"\n";
				cout<<"\nExiting Program..."<<endl<<endl;
				exit(2);
			}
			
			// --- DEM and Flow Directions have equal valid pixels --- 
			if((dem[i][j]==DUMM  && ptr[i][j]!=DUMM) || (dem[i][j]!=DUMM  && ptr[i][j]==DUMM)){
				cout<<"\nWarning: In row " << i <<" the Flow Dir mismatched with DEM"<<endl;
				cout<<"DEM value = "<<dem[i][j]<<"\tDIR value = "<<ptr[i][j]<<"\n";
				cout<<"ROW: "<< i <<"\tCOLUMN: "<<j<<"\n";
				cout<<"\nExiting Program..."<<endl<<endl;
				exit(2);
			}
			
			// --- Flow Directions have values in the range 1 to 8 --- 
			if((ptr[i][j]>8 || ptr[i][j]<1) && dem[i][j]!=DUMM){
				cout<<"\nWarning: In row " << i <<" the Flow Dir contains incorrect values"<<endl;
				cout<<"DEM value = "<<dem[i][j]<<"\tDIR value = "<<ptr[i][j]<<"\n";
				cout<<"ROW: "<< i <<"\tCOLUMN: "<<j<<"\n";
				cout<<"\nExiting Program..."<<endl<<endl;
				exit(2);
			}
		}
	}
	
	//Close input files
	Inp1.close();
	Inp2.close();
	Inp3.close();
	Inp4.close();
	
	cout <<"\n\n\t*************************************\n";
	cout << "\n\tTOTAL NUMBER OF FOUND NON-VOID pixels : "<<cnt<< endl;
	cout << "\n\tMINIMUM ELEVATION FOUND: "<< MIN << endl;
	cout << "\tROW = " << io << "  COLUMN = "<< jo << endl;
	cout << "\n#######  THIS IS CONSIDERED TO BE THE OUTLET #######";
	cout <<"\n\n\t*************************************\n\n";
	
	//  ------------------- HISTOGRAMM COMPUTATION -------------------
	
	// --- Hillslope histogram  --- 
	for (i=0; i < POROG; i++) {
		topoi[i] = new double[count[i]];
		assert(topoi[i] != 0);
		cnt11 += count[i];
		Otp2 << count[i] << endl;
	}
	Otp2<<"THE TOTAL # OF HILLSLOPE PIXELS = "<<cnt11<<endl;
	
	
	// --- Calculating TOPMOODEL statistics (mean and std) ---
	// --- for pixels with contributing area < POROG  --- 
	// --- (for each of contributing area value) ---
	
	cout<<"\nCalculating TOPMODEL statistics..."<<endl;
	double tempo;
	int k, l;
	for (i=0; i < N; i++) {
		for (j=0; j < M; j++) {
			if (dem[i][j] != DUMM) { 
				if (area[i][j] == 0) {
					area[i][j] = 1;  //Assumption
				}
				
				k = ptr[i][j];
				ii = i + pi[k];
				jj = j + pj[k];
				
				if (ptr[i][j]>0 && ii>-1 && jj>-1 && ii<N && jj<M) {
					
					slope = (double)(dem[i][j] - dem[ii][jj])/dist[k];
					
					// -----------------------------------------------
					// 1) [Slope < 0] 
					// ===> Check what we can do if the dem and/or 
					// ===> other coverages are corrupted
					if (slope < 0)
						AdjustNegativeSlope(dem, ptr, dist, i, j, &slope);
					
					// -----------------------------------------------
					// 2) [Slope >= 0]
					if (slope != 0.) {  // Only > 0. slopes 
						tempo = area[i][j]/slope;
						
						// Don't take [log() < 0] values
						if (area[i][j] <= POROG && tempo >= 1.0) {
							l = area[i][j];
							topoi[l-1][ci[l-1]] = log(tempo);
							means[l-1] += log(tempo);
							ci[l-1]++;   // ADD counting index
						}
					}
					else  // I.E. (slope == 0.)
						slope = 1e-5;  // <- Assumption
					
	}
				else
					slope = DUMM;
      }
    }
  }
	
	// ***************************************************************
	//Calculating the parameters of normal distribution: mean and std.
	
	cout<<"\nCalculating parameters of normal distribution..."<<endl;
	
	for (i=0; i < POROG; i++) {
		if (ci[i] > 0)
			means[i] /= ci[i];
		else
			means[i] = 0;
		
		// --- If number of values in a sample is larger than 5 --- 
		if (ci[i] >= 5) {
			for (j=0; j < ci[i]; j++) {
				stds[i] = stds[i] + (topoi[i][j] - means[i])*(topoi[i][j] - means[i]);
			}
			stds[i] = stds[i]/(ci[i] - 1);
			stds[i] = pow(stds[i], 0.5);
		}
		// --- If number of values in a sample is less than 5 --- 
		// --- ensure passing by default --- 
		else
			stds[i] = -1;
		
		// --- Print out mean and std histograms --- 
		Otp3<<means[i]<<"   "<<stds[i]<<endl;
	}
	
	// ***************************************************
	// --- Define where problematic pixels are located --- 
	// --- Print to Opt4 grid called "occurences.asc"  --- 
	cout<<"Defining problematic pixels..."<<endl;
	double gran1, gran2;
	double minn = 9999, maxx = 0;
	double lmaxx = 0;
	int cnt22 = 0;
	
	for (i=0; i < N; i++) {
		for (j=0; j < M; j++) {
			
			if (dem[i][j] > DUMM) { 
				
				k = ptr[i][j];
				ii = i + pi[k];
				jj = j + pj[k];
				
				// OUTLET pixel, assume a 5-degree slope ln[(L^2)/L] 
				if (i == io && j == jo) {  
					effwidth = get_cellD8Width(ptr[i][j]);
					// -- ln[(L^2)/L] --
					aIndex[i][j] = log(area[i][j]/0.0875*dx*dy/effwidth);
					Lambda += aIndex[i][j];
				}
				
				if (ptr[i][j]>0 && ii>-1 && jj>-1 && ii<N && jj<M) {
					
					slope = (double)(dem[i][j] - dem[ii][jj])/dist[k];
					
					if (slope < 0)
						AdjustNegativeSlope(dem, ptr, dist, i, j, &slope);
					
					l = area[i][j];
					
					// ----------------------------------------------
					// ========== CONSIDER HILLSOPE PIXELS ==========
					if (l <= POROG) {    // Only the ones < POROG
						
						// ############
						// To ensure that N < 5 values are taken as they are
						if (stds[l-1] > 0.) {
							gran1 = means[l-1] - 2*stds[l-1];   //(m - 2sigma)
							gran2 = means[l-1] + 2*stds[l-1];   //(m + 2sigma)
							
							// ---------------
							// I.e. non-zero slope -> Finite (area/slope) relation
							if (slope != 0.) {
								tempo = area[i][j]/slope;
								
								if (log(tempo)<gran1) {
									// *** Assign to the LOWER boundary of the conf. interval ***
									tempo = exp(gran1);
									Otp4 << l <<" "; 
								}
								else if (log(tempo)>gran2) {
									// *** Assign to the UPPER boundary of the conf. interval ***
									tempo = exp(gran2);  
									Otp4 << l <<" "; 		    
								}
								else
									// *** If within the specified range ***
									Otp4 << "0 "; 
							}
							
							// ---------------
							// I.e.  zero slope -> Infinite (area/slope) relation
							else {               
								// *** Assign to the UPPER boundary of the conf. interval ***
								tempo = exp(gran2);
								Otp4 << l <<" "; 
							}
						}   // closes if (stds > 0)
						
						// ############
						else {
							if (slope != 0.)
								tempo = area[i][j]/slope;
							else
								tempo = exp(means[l-1]);
							Otp4<<"0 ";		       
						}
					}   //less than POROG
					
					// ----------------------------------------------
					// =========== CONSIDER STREAM PIXELS ===========
					else {
						if (slope != 0.)
							tempo = area[i][j]/slope;
						else {
							// First order ASSUMPTION -- 5 degrees slope
							tempo = area[i][j]/0.0875;
						}
						Otp4 << "0 ";
					}
					
					// By NOW -tempo- should be all set ->> let's get 
					// Topographic index: ln (a_c/tg(b))     
					
					effwidth = get_cellD8Width(ptr[i][j]);
					aIndex[i][j] = log(tempo*dx*dy*1000/effwidth);  // ln[(L^2)/L] - ln[mm]
					Lambda += aIndex[i][j];
					
					if (log(tempo) > lmaxx)
						lmaxx = log(tempo);
				}
				else
					Otp4 << DUMM <<" ";
			}
			else {
				Otp4 << DUMM <<" ";
				Nwt = DUMM;
			}  
		} 
		Otp4 << endl;
	}
	
	// ***************************************************
	// --- Calculation of Topographic Index - Lambda --- 
	cout<<"\nCalculating the mean topographic index..."<<endl;
	Lambda /= cnt;
	cout <<"\n\n\t**************************************************\n";
	cout<<"\n\toooo AREAL INTEGRAL TOPOGRAPHIC INDEX: "<<Lambda<<" oooo"<<endl;
	
	Gamma = log(K0z*AR/F);
	cout<<"\n\toooo AREAL INTEGRAL VALUE OF Gamma: "<<Gamma<<" oooo\n"<<endl;
	cout <<"\t**************************************************\n\n";
	
	// ***************************************************
	// --- Call the getWTmapSiva function  --- 
	getWTmapSiva(dem, area, GWT, aIndex, &maxx, &minn, &cnt22);
	
	
	// ***************************************************
	// --- Calculate average depth to WT BEFORE smoothing --- 
	cout<<"Calculating the average depth before smoothing..."<<endl;
	Zav = 0.; 
	for (i=0; i < N; i++) {
		for (j=0; j < M; j++) {
			if (dem[i][j] == DUMM)
				Otp << DUMM << " ";
			else {
				Otp << GWT[i][j] <<" ";
				Zav += GWT[i][j];
			}
		}
		Otp  << endl;
	}
	Zav /= cnt;
	Zavb = Zav;
	
	// ***************************************************
	// --- Calculating the GW map in Absolute values --- 
	cout<<"Calculating the GW map in absolute values..."<<endl;
	for (i=0; i < N; i++) {
		for (j=0; j < M; j++) {
			if (dem[i][j] > DUMM) { 
				GWTabs[i][j] = dem[i][j] - GWT[i][j]/1000.;
				Otp5<<GWTabs[i][j]<<" ";
			}
			else {
				GWTabs[i][j] = DUMM;
				Otp5<<DUMM <<" ";
			}
		}
		Otp5 << endl;
	}
	
	
	// ***************************************************
	cout <<"\n\n\t*************************************\n";
	cout <<"\t--------- BEFORE SMOOTHING: ---------\n";
	cout <<"\n\tMIN value of GW defined:\t " << minn <<" (or '0')"<< endl;
	cout <<"\tMAX value of GW defined:\t "<< maxx << endl;
	cout <<"\n\tMAX value of ln(a/slope) found:\t "<< lmaxx << endl;
	cout <<"\n\tSum of pixels (area < POROG):           "<<cnt11;
	cout <<"\n\tSum of pixels (area > POROG & Nwt < 0): "<<cnt22;
	cout <<"\n\tTotal: "<<cnt22+cnt11<<" ---> left "<<AREA-cnt22-cnt11
		<<" DRY stream pixels"<<endl;
	cout <<"\n\t*************************************\n";
	
	// ***************************************************
	// --- Smooting the GW Table topography --- 
	for (i=0; i < SMOO; i++) 
		smoothWTrelief(dem, GWT, GWTabs, dist5);
	
	// ***************************************************
	// --- Printing the depth to WT Table in stream network --- 
	for (i=0; i < N; i++) {
		for(j=0; j < M; j++) {
			if(dem[i][j] > DUMM && area[i][j] > POROG)
				Otp8 << GWT[i][j] <<" ";
			else
				Otp8 << DUMM <<" ";
		}
		Otp8 << endl;
	}
	
	// ***************************************************
	// ---  Calculate average depth to WT AFTER smoothing --- 
	cout<<"\nCalculating the average depth after smoothing..."<<endl;
	Zav = 0.; 
	for (i=0; i < N; i++) {
		for(j=0; j < M; j++) {
			if (GWT[i][j] != DUMM) {
				if (GWT[i][j] >= 32000.) {
					GWTabs[i][j] += GWT[i][j]/1000. - 32.; // Meters
					GWT[i][j] = 32000.;
					cout <<"\nWarning: Pixel with Nwt > 32000mm\tNwt assigned to 32000 mm";
				}
				Otp6 << GWTabs[i][j] <<" ";
				Otp7 << GWT[i][j]/FACTOR <<" ";	  
				Zav += GWT[i][j];
			}
			else {
				Otp7 << GWT[i][j] <<" ";	
				Otp6 << GWTabs[i][j] <<" ";
			}  
		}
		Otp6 << endl;
		Otp7 << endl;
	}
	Zav /= cnt;
	Zava = Zav;
	
	// --- Printing summary results to screen --- 
	cout<<"\n\nSummary results:"<<endl;
	cout<<"----------------"<<endl;
	cout<<"Total Number of watershed pixels:              "<<cnt<<endl;
	cout<<"Areal Integral of Topographic Index:           "<<Lambda<<endl;
	cout<<"Areal Integral of Gamma:                       "<<Gamma<<endl;
	cout<<"Average Depth to Water Table before smoothing: "<<Zavb<<" mm"<<endl;
	cout<<"Average Depth to Water Table after smoothing:  "<<Zava<<" mm"<<endl;
	cout<<"\nMinimum value of GW defined (zero):            "<<minn<<" mm"<<endl;
	cout<<"Maximum value of GW defined:                   "<<maxx<<" mm"<<endl;
	cout<<"Maximum value of ln(a/slope):                  "<<lmaxx << endl;
	cout<<"\nSum of hillslope pixels (area < POROG):        "<<cnt11<< endl;
	cout<<"Sum of pixels (area > POROG & Nwt < 0):        "<<cnt22<< endl;
	cout<<"GW depth was divided by factor:                "<<FACTOR<<endl; 
	
	
	// --- Summarize file output --- 
	cout<<"\n\nFile output located in Input/waterTable/: "<<endl;
	cout<<"-----------------------------------------"<<endl;
	cout<<"Hillslope histogram:      _hillslope.hist"<<endl;
	cout<<"Mean and Std histogram:   _MeanStd.hist"<<endl;
	cout<<"Occurences Grid:          _occurence.asc"<<endl;
	cout<<"Absolute GWT Grid:        _GWTabs.asc"<<endl;
	cout<<"Smoothed Abs GWT Grid:    _GWTabs_smooth.asc"<<endl;
	cout<<"Depth GWT Grid:           _GWTdepth.asc" <<endl;
	cout<<"Smoothed Depth GWT Grid:  _GWTdepth_smooth.asc"<<endl;
	cout<<"GWT Depth in Stream Grid: _depth_instream.asc"<<endl; 
	cout<<"-----------------------------------------"<<endl<<endl;
	
	//Delete arrays
	delete [] count, ci, means, stds; 
	for(int ct=0;ct<N;ct++){
		delete area[ct];
		delete dem[ct];
		delete ptr[ct];
		delete GWT[ct];
		delete GWTabs[ct];
		delete aIndex[ct];
	}
	delete [] area, dem, ptr, GWT, GWTabs, aIndex;  
	delete [] topoi;
	
	return 0;
}

/***************************************************************************
**
** smoothWTrelief( )
**
***************************************************************************/

void smoothWTrelief(double **dem, double **GWT, double **GWTabs, double *dist5){
	
	int i, j, ii, jj, l;
	double ves = 0;
	double tempo;
	double Nwt;
	double minn = 9999, maxx = 0;
	int cnt = 0;
	
	cout <<"\n...Smoothing the GWT topography...";
	
	double **gwta;     //   Temporary array to store smoothed GW values
	gwta = new double* [N];
	assert(gwta != 0);
	
	for (i=0; i < N; i++) {
		
		gwta[i] = new double[M];
		assert(gwta[i] != 0);
		
		for(j=0; j < M; j++) {
			if(dem[i][j] != DUMM) { 
				tempo = GWTabs[i][j]; // Start from the current pixel
				ves = 1;
				
				for (l=1; l < MA; l++) {
					ii = i + pi5[l];
					jj = j + pj5[l];
					if (ii>-1 && jj>-1 && ii<N && jj<M && dem[ii][jj] != DUMM) {
						ves += dx/dist5[l];
						tempo += GWTabs[ii][jj]*dx/dist5[l];
					}
				}
				
				tempo /= ves;                     // absolute value of WT  
				Nwt = (dem[i][j] - tempo)*1000.;  // depth to WT in mm
				
				if (Nwt > maxx)
					maxx = Nwt;
				else if (Nwt < minn)
					minn = Nwt;
				
				// ---------------------------
				// Corrections... 
				if (Nwt < 0.) {        // cannot go above topography
					cnt++;	         
					Nwt   = GWT[i][j] ;     // assigned to where it was
					tempo = GWTabs[i][j];   // assigned to where it was
				}
				
				if (Nwt > maxx)
					maxx = Nwt;
				else if (Nwt < minn)
					minn = Nwt;
				
				GWTabs[i][j] = tempo;  // change to smoothed value
				GWT[i][j] = Nwt;       // change to smoothed value
				
			}
		}
	}
	
	for (i=0; i < N; i++)
		delete gwta[i];
	
	delete gwta;
	
	cout <<"\nIn total "<<cnt<<" pixels have been found with Nwt < 0";
	cout <<"\nMAX value of GW defined after smoothing: "<< maxx << endl;
	
	return;
}

/***************************************************************************
**
** zaglubiStream( )
**
***************************************************************************/

void zaglubiStream(double **dem, int **area, double **GWT, double ZZ)
{
	int i, j;
	for (i=0; i < N; i++) {
		for(j=0; j < M; j++) {
			if(dem[i][j] > DUMM && area[i][j] > POROG)
				GWT[i][j] += ZZ;   // goes deeper as much as ZZ
		}
	}
	return;
}

/***************************************************************************
**
** assignToZero( )
**
***************************************************************************/

void assignToZero(double **dem, int **area, double **GWT)
{
	int i, j;
	for (i=0; i < N; i++) {
		for(j=0; j < M; j++) {
			if(dem[i][j] > DUMM && area[i][j] > POROG)
				GWT[i][j] = 0.;
		}
	}
	return;
}

/***************************************************************************
**
** getMoisture()
**
***************************************************************************/

double getMoisture(double gw)
{
	double mc; 
	mc = (POROS - RESID)*POREIND/F*(1 - exp(-F*gw/POREIND)) + RESID*gw;  
	return mc;
}

/***************************************************************************
**
** get_cellEffectiveWidth(draindir)
**
** Arguments: int draindir
** Objective: Compute cell width in the effective direction for multiple
**     outflow directions
** Return value: cell width
** Algorithm:
**  switch on draiudir
**  get width
**
***************************************************************************/

double get_cellEffectiveWidth(int draindir){ 
	double Width;
	switch(draindir){  
		case 1: { Width = (0.5 * dy); break; }
		case 2: { Width = (0.354 * sqrt((dx*dx+dy*dy)/2.)); break; }
		case 3: { Width = (0.5 * dx); break; }
		case 4: { Width = (0.354 * sqrt((dx*dx+dy*dy)/2.)); break; }
		case 5: { Width = (0.5 * dy); break; }
		case 6: { Width = (0.354 * sqrt((dx*dx+dy*dy)/2.)); break; }
		case 7: { Width = (0.5 * dx); break; }
		case 8: { Width = (0.354 * sqrt((dx*dx+dy*dy)/2.)); break; }
		default: Width = (-1.0);
	}
	return Width*1000.;  // Return in mm
}


/***************************************************************************
**
** get_cellD8Width()
**
***************************************************************************/

double get_cellD8Width(int draindir){
	
	double Width;
	switch(draindir){  
		case 1: { Width = (dy); break; }
		case 2: { Width = (sqrt(dx*dx+dy*dy)); break; }
		case 3: { Width = (dx); break; }
		case 4: { Width = (sqrt(dx*dx+dy*dy)); break; }
		case 5: { Width = (dy); break; }
		case 6: { Width = (sqrt(dx*dx+dy*dy)); break; }
		case 7: { Width = (dx); break; }
		case 8: { Width = (sqrt(dx*dx+dy*dy)); break; }
		default: Width = (-1.0);
	}
	return Width;  // Return in meters
}

/***************************************************************************
**
** getWTmapSiva()
**
***************************************************************************/

void getWTmapSiva(double **dem, int **area, double **GWT, double **aIndex, 
				  double *maxx, double *minn, int *cnt22)
{
	int i,j;
	double Nwt;
	
	for (i=0; i < N; i++) {
		for (j=0; j < M; j++) {
			if (dem[i][j] != DUMM) { 
				
				Nwt = Zav - 1/F*( (aIndex[i][j] - Lambda) - (log(K0z*AR/F) - Gamma) ); 
				// Last two terms are supposed to calculate deviations when non-uniform
				
				if (Nwt > *maxx)
					*maxx = Nwt;
				else if (Nwt < *minn)
					*minn = Nwt;
				
				// CHECK FOR VALUE OF WATER TABLE POSITION
				if (Nwt < 0) {
					//cout <<"Nwt < 0: ln(a/slope) = " <<aIndex[i][j]
					//     <<";  Cont.AREA: "<<area[i][j]<<endl;
					Nwt = 0.;  // tempo in this case gets of the order 1e+6 - 1e+9
					
					if (area[i][j] < POROG) {
						//cout <<"\tNwt is negative for hillslope: I = "<<i<<" J = "<<j;
						//cout <<"\tCONTR. AREA = " << area[i][j] << endl;
					}
					else 
						*cnt22 = *cnt22+1;
				}
			}
			else {
				Nwt = DUMM; //Changed from DUMM to Zero since DUMM = -9999
			}
			GWT[i][j] = Nwt;
		}
	}
}


/***************************************************************************
**
**  AdjustNegativeSlope()
**
***************************************************************************/
void AdjustNegativeSlope(double **dem, int **ptr, double *dist, 
						 int i, int j, double *slope)
{
	int cntr = 0;
	int ii, jj, k;
	double tempo;
	
	k  = ptr[i][j];
	ii = i + pi[k];
	jj = j + pj[k];
	
	cout<<"\nNegative slope:  Row "<<i<<" Col "<<j
		<<"; dem0 = "<<dem[i][j]<<"; dem1 = "<<dem[ii][jj];
	
	for (int drainDir=1; drainDir < 9; drainDir++) {
        //===> Search for downslope pixel...
        ii = i + pi[drainDir];
        jj = i + pj[drainDir];
		
        if (ii>-1 && jj>-1 && ii<N && jj<M) {
			
			if (dem[i][j] >= dem[ii][jj] && dem[ii][jj] > 0) {
				tempo = (double)(dem[i][j]-dem[ii][jj])/dist[drainDir];
				if (tempo > (*slope)) {
					*slope = tempo;
					ptr[i][j] = drainDir;
				}
			}
		}
	}
	if (*slope > 0) {
		k  = ptr[i][j];
		ii = i + pi[k];
		jj = j + pj[k];
		cout<<"... adjusted to positive; dem1 = "<<dem[ii][jj];
	}
	else  {
		cout<<"... positive not found, assigned to 0";
		*slope = 0.; // <- this will be taken care of
	}
	
	return;
}

//=========================================================================
//
//
//                    End of InitialGW.cpp
//
//
//=========================================================================
