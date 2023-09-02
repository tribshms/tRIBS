/***************************************************************************
**
**  		     tRIBS Distributed Hydrologic Model
**
**              TIN-based Real-time Integrated Basin Simulator
**		         Ralph M. Parsons Laboratory
**  		    Massachusetts Institute of Technology
**
**
**  tSnowIntercept.cpp:   Function file for tSnowIntercept class 
**			    (see tSnowIntercept.h)
**
***************************************************************************/

#include "src/Headers/globalIO.h"

/****************************************************************************
**
**	      tSnowIntercept -- Constructor and Destructor
**
****************************************************************************/

//---------------------------------------------------------------------------
//
//		tSnowIntercept() Constructor and Destructor
//
//    These also necessarily call the tEvapoTrans constructors and destructor
//    as tSnowIntercept is derived from tEvapoTrans.
//
//---------------------------------------------------------------------------

tSnowIntercept::tSnowIntercept(){

}

tSnowIntercept::tSnowIntercept(SimulationControl *simCtrPtr, tMesh<tCNode> *gridRef, 
		       tInputFile &inFile, tRunTimer *tmrptr, tResample *resamp,
		       tHydroModel *hydro, tRainfall *storm)
: tEvapoTrans( simCtrPtr, gridRef, inFile, tmrptr, resamp, hydro, storm )

{

  SetSnowVariables(inFile);
  SetSnowInterceptVariables();

}

tSnowIntercept::~tSnowIntercept() {

  Cout << "tSnowIntercept Object has been destroyed..." << endl;
}

/*******************************************************************
**
**	    tSnowIntercept -- Initialization Routines
**
*******************************************************************/

//------------------------------------------------------------------
//
//		tSnowIntercept::SetSnowInterceptVariables()
//
//	  Initializes the basic variables needed for the Liston and
//	  Sturm (2006) implementation.
//
//------------------------------------------------------------------

void tSnowIntercept::SetSnowInterceptVariables()
{

  Qcs = Ce = I = psiS = 0.0;
  Imax = prec = LAI = 0.0;
  RH = D = rhoVap = Omega = 0.0;
  Sh = Nu = Re = 0.0;
  kc = 0.010; //-
  iceRad = 500e-6; //m
  Mwater = 18.01; //kg/kmol
  R = 8313; //J/kmol K
  RdryAir = 287; //J/kg K
  nu = 1.3e-5; //m^2/s
  KtAtm = 0.025; //J/msK
  esatIce = 0.0;
  beta = 0.9;
  acoefficient = 0.0;
  Lm = 0.0;
}
//---------------------------------------------------------------------------
//
//	tSnowIntercept::SetSnowVariables()
//
//	Auxiliary function used in robust constructor to initialize snow 
//	variables. This is the same as in tSnowPack w/ only a few minor
//	differences.
//
//---------------------------------------------------------------------------

void tSnowIntercept::SetSnowVariables(tInputFile &infile)
{

  //options	
  hillAlbedoOption = infile.ReadItem(hillAlbedoOption,"HILLALBOPT");

  //AJR2008, SKY2008Snow
  timeStepm = timeStep;  // Code originally assumed this was timesteps in hours but is minutes CJC 2022
  
  timeSteph = timeStepm/60;  // Code addition to calculate time step in hours CJC2022
  timeSteps = 60*timeStepm;

  //state variables
  liqWE = iceWE = snWE = 0.0;
  liqWEm = iceWEm = snWEm = 0.0;
  Utot = Usn = Uwat = 0.0;
  H = L = G = Prec = Rn = 0.0;
  liqTempC = iceTempC = snTempC = 0.0;
  liqTempK = iceTempK = snTempK = 0.0;
  snPrec = liqPrec = 0.0;
  snPrecm = liqPrecm = 0.0;
  snPrecmm = liqPrecmm = 0.0;
  vapPressSmb = vapPresskSPa = 0.0;
  rholiqcgs = 1.0;
  rhoicecgs = 0.94;
  rhosncgs = 0.1;
  rholiqkg = 1000.0;
  rhoicekg = 920.0;
  rhosnkg = 100.0;
  rhoAir = 1.3;
  crustAge = 0.0;

  //parameters
  cpicekJ = 2.1;
  cpwaterkJ = 4.190;
  cpairkJ = 1.01;
  latFreezekJ = 334;
  latVapkJ = 2470;
  latSubkJ = latFreezekJ + latVapkJ;
  albedo = 0.80;
  minutelyTimeStep = 0.0;
  
  //output variables
  snDepth = snDepthm = 0;
	  
  //conversions
  naughttokilo = 1e-3;
  kilotonaught = 1e3;
  cgsRHOtomks = 1e3;
  mksRHOtocgs = 1e-3;
  naughttocm = 100;
  cmtonaught = 0.01;
  ctom = 10;
  mtoc = 0.1;
  htos = 3600;

}

/****************************************************************************
**
**		    tSnowIntercept -- Calling Functions
**
****************************************************************************/




/****************************************************************************
**
**		  tSnowIntercept -- Conversion Functions
**
****************************************************************************/


//---------------------------------------------------------------------------
//
//			  tSnowIntercept::CtoK()
//	
//	Converts Celsius to Kelvin
//---------------------------------------------------------------------------

double tSnowIntercept::CtoK(double temperature) 
{

  return (temperature + 273.15);
  
}

//---------------------------------------------------------------------------
//
//			tSnowIntercept::KtoC()
//
//	Converts Kelvin to Celsius
//---------------------------------------------------------------------------

double tSnowIntercept::KtoC(double temperature)
{

  return (temperature - 273.15);
  
}

/****************************************************************************
**
**		tSnowIntercept -- EB Functions (Basic Caculations)
**
**	This set of functions is primarily taken from tSnowPack and is
**	included for further developments in snow interception models.
**
****************************************************************************/


double tSnowIntercept::inShortWaveSn()
{
  double Is, N, Iv, Isw, Ir;
  double v, t, cosi, scover;
  double RadGlobClr;
      
  double albedo(0.8);

  //Remaining variables in DirectDiffuse from v3 -- AJR2008, SKY2008Snow
  //  So, entire function changed to match inShortWave in tEvapoTrans
/*  double h0, m, pp0, Dh0ref, h0ref, drm, Tlinke;
  double TnTLK, Fdh0, A1p, A1, A2, A3;
  double pi = 4*atan(1.0);*/

  Ic=Is=Id=Ir=Ids=Ics=Isw=Iv=0.0;

  // Elevation, Slope and Aspect have been set before

  if (alphaD > 0.0) {

    DirectDiffuse(Tlinke, elevation);  // SKY2008Snow, AJR2007

    // Cloud cover information
    if (fabs(skyCover-9999.99) < 1.0E-3) {

	skyCover = compSkyCover();//ADDED BY RINEHART 2007 @ NMT
					// computes sky cover from relative
					// humidity and rain.
	scover = skyCover;

	//if (rain > 0.0) scover = 10.0;
	//else            scover = 1.0;
    }
    else 
	scover = skyCover;
    N = scover/10.0;
		
    // If observations (for a horizontal surface) exist - 
    // use them, at least in an approximate manner
    if (tsOption > 1 && !rainPtr->getoptStorm()) {
	RadGlobClr = (RadGlbObs/(1.0-0.65*pow(N,2.0)));
	Ic = Ic/(Ic*sinAlpha + Id)*RadGlobClr;
	Id = RadGlobClr - Ic*sinAlpha;
    }

    // 1) Slope aspect
    //Account for the aspect and slope of the element 
    //Estimate 'cosi' and compare it with the Sun position
    //  'cosi' = cos(i), where 'i' is the angle between 
    //  the sun beam and the normal to the slope surface
    //
    //Rinehart 2007 @ New Mexico Tech
    //
    //	We have incorporated sheltering options. Option 3 is
    //	no topographic shading. Option 0, the default, is
    //	local topographic shading. Option 1 is incorporation
    //	of horizon angles in calc of SV and LV. Option 2 is 
    //	the total integration of local and remote sheltering.
    //
    //	Here, if any sheltering is turned on, then we calculate
    //	the local controls of slope and aspect.
    //
    //	RMK: SLOPE AND ASPECT ARE CALCULATED FROM THE FLOW EDGE
    //	AND ARE IN RADIANS.
    
    if (shelterOption < 4) {//CHANGED IN 2008
      
      cosi = 1.0*(cos(slope)*sinAlpha + sin(slope)*cos(asin(sinAlpha))*cos(sunaz-aspect));

      if (cosi >= 0.0) {
	Ics = Ic*cosi;
      }
      else {
	Ics = 0.0;
      }
    }
    else {
	Ics = 1.0*Ic;
    }

    if ( (shelterOption == 2) || (shelterOption == 1) ) {//CHANGED IN 2008
	    
	Ics *= aboveHorizon(ID); //check to see if we can see the sun (aboveHorizon() in tEvapoTrans)
    }

    //2) Horizon factor for diffuse radiation?
    //Rinehart 2007 @ New Mexico Tech
    //
    //	See comment above about sheltering options.
    
    if ((shelterOption > 0)&&(shelterOption < 3)) {
      v = shelterFactorGlobal; //incorporate remote sheltering
    }
    else if (shelterOption == 0 || shelterOption == 3){ //CHANGED 2008
      v = 0.5*(1 + cos(slope)); //local sheltering
    }
    else {
      v = 1.0; // no sheltering
    }
    
    Ids = Id*v;

    // 3) Account for cloud cover
    Is = (1.0-0.65*pow(N,2))*(Ics + Ids);
    
    //Reflected from surrounded sites radiation
    //
    //Modified by Rinehart 2007 @ New Mexico Tech
    //

    if (hillAlbedoOption == 0) {
	hillalbedo = albedo; }
    else if (hillAlbedoOption == 1) {
	hillalbedo = coeffAl; 
    }
    else if (hillAlbedoOption == 2) {
	hillalbedo = (1-coeffV)*albedo + coeffV*coeffAl;//changed as do not see algorithm if no snow
    }

    if (shelterOption == 0) {
      //local
      Ir = hillalbedo*Is*(1-cos(slope))*0.5;  
      Is += Ir;
    }
    else if ( (shelterOption > 1)&&(shelterOption < 4) ) { //CHANGED IN 2008
      //remote    
      Ir = hillalbedo*Is*( 0.5*(1 + cos(slope)) - shelterFactorGlobal); //CHANGED IN 2008
      landRefGlobal = 0.5*(1 + cos(slope)) - shelterFactorGlobal;
      Is += Ir;

    }
    else { //CHANGED IN 2008
      Ir = 0.0;
      Is = Is;
    }

    // Account for albedo
    Isw = Iv*(1.0 - albedo);

  } //end -- alphaD > 0
  else {
    Ic=Is=N=Iv=Isw=Id=Ids=Ics=Ir=0.0;
  } // end -- alphaD <= 0
 
  return Isw;
}

/****************************************************************************
**
**		      tSnowIntercept -- I/O Functions
**
**	Functions that interact w/ the *.in file and get options.
**
****************************************************************************/

//---------------------------------------------------------------------------
//
//			tSnowIntercept::getSnowOpt()
//
//    Retrieves the snow option found during initialization.
//
//---------------------------------------------------------------------------


/***************************************************************************
**
** tSnowIntercept::writeRestart() Function
**
** Called from tSimulator during simulation loop
**        
***************************************************************************/
void tSnowIntercept::writeRestart(fstream & rStr) const
{


  BinaryWrite(rStr, timeSteph);
  BinaryWrite(rStr, timeSteps);
  BinaryWrite(rStr, timeStepm);
  BinaryWrite(rStr, minutelyTimeStep);

  BinaryWrite(rStr, liqWE);
  BinaryWrite(rStr, iceWE);
  BinaryWrite(rStr, snWE);
  BinaryWrite(rStr, liqRoute);
  BinaryWrite(rStr, liqWEm);
  BinaryWrite(rStr, iceWEm);
  BinaryWrite(rStr, snWEm);
  BinaryWrite(rStr, liqRoutem);
  BinaryWrite(rStr, Utot);
  BinaryWrite(rStr, Usn);
  BinaryWrite(rStr, Uwat);
  BinaryWrite(rStr, liqTempC);
  BinaryWrite(rStr, iceTempC);
  BinaryWrite(rStr, snTempC);
  BinaryWrite(rStr, liqTempK);
  BinaryWrite(rStr, iceTempK);
  BinaryWrite(rStr, snTempK);
  BinaryWrite(rStr, crustAge);
  BinaryWrite(rStr, albedo);
  BinaryWrite(rStr, hillalbedo);

  BinaryWrite(rStr, H);
  BinaryWrite(rStr, L);
  BinaryWrite(rStr, G);
  BinaryWrite(rStr, Prec);
  BinaryWrite(rStr, Rn);
  BinaryWrite(rStr, snPrec);
  BinaryWrite(rStr, liqPrec);
  BinaryWrite(rStr, snPrecm);
  BinaryWrite(rStr, liqPrecm);
  BinaryWrite(rStr, snPrecmm);
  BinaryWrite(rStr, liqPrecmm);
  BinaryWrite(rStr, vapPressSmb);
  BinaryWrite(rStr, vapPresskSPa);
  BinaryWrite(rStr, rholiqcgs);
  BinaryWrite(rStr, rhoicecgs);
  BinaryWrite(rStr, rhosncgs);
  BinaryWrite(rStr, rholiqkg);
  BinaryWrite(rStr, rhoicekg);
  BinaryWrite(rStr, rhosnkg);
  BinaryWrite(rStr, rhoAir);

  BinaryWrite(rStr, cpsnowkJ);
  BinaryWrite(rStr, cpicekJ);
  BinaryWrite(rStr, cpwaterkJ);
  BinaryWrite(rStr, cpairkJ);
  BinaryWrite(rStr, latFreezekJ);
  BinaryWrite(rStr, latVapkJ);
  BinaryWrite(rStr, latSubkJ);
  BinaryWrite(rStr, resFact);

  BinaryWrite(rStr, snDepth);
  BinaryWrite(rStr, snDepthm);

  BinaryWrite(rStr, naughttokilo);
  BinaryWrite(rStr, kilotonaught);
  BinaryWrite(rStr, cgsRHOtomks);
  BinaryWrite(rStr, mksRHOtocgs);
  BinaryWrite(rStr, naughttocm);
  BinaryWrite(rStr, cmtonaught);
  BinaryWrite(rStr, ctom);
  BinaryWrite(rStr, mtoc);
  BinaryWrite(rStr, htos);

  tEvapoTrans::writeRestart(rStr);
}

/***************************************************************************
**
** tSnowIntercept::readRestart() Function
**
***************************************************************************/
void tSnowIntercept::readRestart(fstream & rStr)
{


  BinaryRead(rStr, timeSteph);
  BinaryRead(rStr, timeSteps);
  BinaryRead(rStr, timeStepm);
  BinaryRead(rStr, minutelyTimeStep);

  BinaryRead(rStr, liqWE);
  BinaryRead(rStr, iceWE);
  BinaryRead(rStr, snWE);
  BinaryRead(rStr, liqRoute);
  BinaryRead(rStr, liqWEm);
  BinaryRead(rStr, iceWEm);
  BinaryRead(rStr, snWEm);
  BinaryRead(rStr, liqRoutem);
  BinaryRead(rStr, Utot);
  BinaryRead(rStr, Usn);
  BinaryRead(rStr, Uwat);
  BinaryRead(rStr, liqTempC);
  BinaryRead(rStr, iceTempC);
  BinaryRead(rStr, snTempC);
  BinaryRead(rStr, liqTempK);
  BinaryRead(rStr, iceTempK);
  BinaryRead(rStr, snTempK);
  BinaryRead(rStr, crustAge);
  BinaryRead(rStr, albedo);
  BinaryRead(rStr, hillalbedo);

  BinaryRead(rStr, H);
  BinaryRead(rStr, L);
  BinaryRead(rStr, G);
  BinaryRead(rStr, Prec);
  BinaryRead(rStr, Rn);
  BinaryRead(rStr, snPrec);
  BinaryRead(rStr, liqPrec);
  BinaryRead(rStr, snPrecm);
  BinaryRead(rStr, liqPrecm);
  BinaryRead(rStr, snPrecmm);
  BinaryRead(rStr, liqPrecmm);
  BinaryRead(rStr, vapPressSmb);
  BinaryRead(rStr, vapPresskSPa);
  BinaryRead(rStr, rholiqcgs);
  BinaryRead(rStr, rhoicecgs);
  BinaryRead(rStr, rhosncgs);
  BinaryRead(rStr, rholiqkg);
  BinaryRead(rStr, rhoicekg);
  BinaryRead(rStr, rhosnkg);
  BinaryRead(rStr, rhoAir);

  BinaryRead(rStr, cpsnowkJ);
  BinaryRead(rStr, cpicekJ);
  BinaryRead(rStr, cpwaterkJ);
  BinaryRead(rStr, cpairkJ);
  BinaryRead(rStr, latFreezekJ);
  BinaryRead(rStr, latVapkJ);
  BinaryRead(rStr, latSubkJ);
  BinaryRead(rStr, resFact);

  BinaryRead(rStr, snDepth);
  BinaryRead(rStr, snDepthm);

  BinaryRead(rStr, naughttokilo);
  BinaryRead(rStr, kilotonaught);
  BinaryRead(rStr, cgsRHOtomks);
  BinaryRead(rStr, mksRHOtocgs);
  BinaryRead(rStr, naughttocm);
  BinaryRead(rStr, cmtonaught);
  BinaryRead(rStr, ctom);
  BinaryRead(rStr, mtoc);
  BinaryRead(rStr, htos);

  tEvapoTrans::readRestart(rStr);
}

/*****************************************************************************
**
**		tSnowIntercept -- END OF TSNOWINTERCEPT.CPP
**
*****************************************************************************/
