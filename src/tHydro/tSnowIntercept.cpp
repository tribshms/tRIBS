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

#include "src/tHydro/tSnowIntercept.h"
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

  SetSnowVariables(inFile, hydro);
  SetSnowInterceptVariables(inFile, hydro);  

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

void tSnowIntercept::SetSnowInterceptVariables(tInputFile &infile, tHydroModel *hydro)
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

  return;
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

void tSnowIntercept::SetSnowVariables(tInputFile &infile, tHydroModel *hydro)
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
  
  return;
}

/****************************************************************************
**
**		    tSnowIntercept -- Calling Functions
**
****************************************************************************/

//---------------------------------------------------------------------------
//
//		tSnowIntercept::callSnowIntercept()
//
//    Calls the physical algorithms from tSnow::callSnowPack(). Some of 
//    tIntercept::callIntercept() is implemented for the case when there is
//    no snow.
//
//---------------------------------------------------------------------------

void tSnowIntercept::callSnowIntercept(tCNode *node, tIntercept *interceptModel)
{
  double CanStorage, ctos, evapWetCanopy;
  double subFrac, unlFrac, precip, Isnow, throughfall;// SKY2008Snow, AJR2008
  int count;
  count = 0;


  slope = fabs(atan(node->getFlowEdg()->getSlope()));
  aspect = node->getAspect();
  elevation = node->getZ();
  CanStorage = node->getCanStorage();


  SetSunVariables();

  // SKY2008Snow, AJR2008
  // Resample Meteorological Grids, if option
  if (metdataOption == 2){
	resampleGrids(timer);
  }

  //Derive the remote sheltering parameters for use in the computation
  //  of radiation. Use only when remote sheltering is turned on.

  //  Similar to tEvapoTrans::callEvapoPotential and tSnowPack::callSnowPack
  if ( (shelterOption > 0) && (shelterOption < 4) ) {//CHANGED IN 2008    
    for (int tempIndex = 0; tempIndex < 16; tempIndex++) {
      switch ( tempIndex ) {
	  case 0:
	    ha2700 = node->getHorAngle2700();
	    break;
	  case 1:
	    ha2925 = node->getHorAngle2925();
	    break;
	  case 2:
	    ha3150 = node->getHorAngle3150();
	    break;
	  case 3:
	    ha3375 = node->getHorAngle3375();
	    break;
	  case 4:
	    ha0000 = node->getHorAngle0000();
	    break;
	  case 5:
	    ha0225 = node->getHorAngle0225();
	    break;
	  case 6:
	    ha0450 = node->getHorAngle0450();
	    break;
	  case 7:
	    ha0675 = node->getHorAngle0675();
	    break;
	  case 8:
	    ha0900 = node->getHorAngle0900();
	    break;
	  case 9:
	    ha1125 = node->getHorAngle1125();
	    break;
	  case 10:
	    ha1350 = node->getHorAngle1350();
	    break;
	  case 11:
	    ha1575 = node->getHorAngle1575();
	    break;
	  case 12:
	    ha1800 = node->getHorAngle1800();
	    break;
	  case 13:
	    ha2025 = node->getHorAngle2025();
	    break;
	  case 14:
	    ha2250 = node->getHorAngle2250();
	    break;
	  case 15:
	    ha2475 = node->getHorAngle2475();
	    break;
	  default:
	    cout << "\nCheck tempInd -- did not exist or assign" << endl;
      }//end-switch
    }//end-for
    shelterFactorGlobal = node->getSheltFact();//remote sheltering
  }
  else if (shelterOption == 0) {
    shelterFactorGlobal = 0.5*(1 + cos(slope));//local sheltering
    node->setSheltFact(shelterFactorGlobal);
  }
  else {
    shelterFactorGlobal = 1.0;//no sheltering
  }

  //find meteorolgical conditions
  rHumidity = node->getRelHumid();
  airTemp = node->getAirTemp();
  airTempK = CtoK(airTemp);
  precip = node->getRain()*coeffV; //precip scaled by veg fraction
  windSpeed = node->getWindSpeed();
  //potEvap = node->getPotEvap();

  //set vegetation parameters from table
  setCoeffs(node);
  if (luOption == 1) {	
	newLUGridData(node);
  }
  LAI = coeffLAI;
  
  //reinitialize snow interception model
  Iold = rholiqkg*cmtonaught*(node->getIntSWE());
  Qcs = 0.0;
  Lm = 0.0;

  // Check snowpack conditions so surface temps are adequately set
  liqWE = node->getLiqWE(); //cm
  iceWE = node->getIceWE(); //cm
  snWE = liqWE + iceWE; //cm

  Tso = node->getSurfTemp() + 273.15;
  Tlo = node->getSoilTemp() + 273.15;


  //TODO WR-WB debug what happens to excess canopy storage?

  if ( (precip*snowFracCalc() < 1e-4) && (Iold < 1e-3) ) {
      //The below code block account for the case where there is no snow in canopy and it's not snowing
      // but could be raining In short this should accounts for the case of rain on snow, where the
      // net precip is then routed to SnowPack.
      I = Iold = Lm = Qcs = 0.0;

      // Below block of code added by WR 6/21/23 to set potEvap
      //Calculate the Potential and Actual Evaporation
      if(evapotransOption == 1){
          EvapPenmanMonteith(node); // call to get EvapPot, but energy balance for soil es
      }
      else if(evapotransOption == 2){
          EvapDeardorff(node); // SKY2008Snow
      }
      else if(evapotransOption == 3){
          EvapPriestlyTaylor(node); // SKY2008Snow
      }
      else if(evapotransOption == 4){
          EvapPan();
      }
      else{
          cout << "\nEvapotranspiration Option " << evapotransOption;
          cout <<" not valid." << endl;
          cout << "\tPlease use :" << endl;
          cout << "\t\t(1) for Penman-Monteith Method" << endl;
          cout << "\t\t(2) for Deardorff Method"<< endl;
          cout << "\t\t(3) for Priestly-Taylor Method" << endl;
          cout << "\t\t(4) for Pan Evaporation Measurements" << endl;
          cout << "Exiting Program...\n\n"<<endl;
          exit(1);
      }
      // Set ground element-scale fluxes to zero and snow temp--just so
      node->setNetRad(0.0);
      node->setGFlux(0.0);
      node->setHFlux(0.0);
      node->setLFlux(0.0);
      node->setLongRadOut(0.0);
      node->setSurfTemp(node->getSnTempC()); //assume surface temp == snow temp, note if snWE < 1e-4 snow intercept should not be called
      node->setSoilTemp(node->getSnTempC()); //this should be updated with n-layer snow model

      setToNode(node);

      ComputeETComponents(interceptModel, node, count, 1);

      //Set canopy snow components to 0
      node->setIntSWE(0);
      node->setIntSnUnload(0);
      node->setIntSub(0);
      node->setIntPrec(0);

  }//end -- no snow

  //snowing with or without snow in canopy
  else {
    
    albedo = 0.8;

    // Canopy storage in mm is same value as kg/m^2 when converted
    // add to I_old and reset to node canopy storage to zero
    if(CanStorage>1e-5){
        Iold += CanStorage; //CanStorage has been scaled by coeffV
        node->setCanStorage(0.0);
    }

    //precip in mm is same value when converted to kg/m^2

    //maximum mass of snow stored in canopy (kg/m^2)
    Imax = 4.4*LAI;
   
    //compute new intercepted snow (kg/m^2)
    Isnow = 0.7*(Imax - Iold)*(1 - exp(-precip/Imax));
    I = Iold + Isnow;

    //precip minus intercepted snow (i.e. throughfall)
    throughfall = precip - Isnow; //convert to mm while setting to node.

    //if there was old snow, sublimate and unload

    if (Iold > 0.0) {
      computeSub();
      computeUnload();
    }
    else {
      Qcs = 0.0;//sublimation term
      Lm = 0.0;//unloading term
    }

    //adjust amount of snow in canopy
    Iold = I;

    I += Qcs - Lm; //I == interception (kg), Qcs == sublimation (kg) (sign computed), Lm == unloading (computed positive) (kg)

    // SKY2008Snow based on AJR2008's recommendation starts here (water balance now preserved)
    if (I < 0.0) {

	    if (Qcs < 0.0) {
		    subFrac = fabs(Qcs)/(fabs(Qcs)+Lm);
		    unlFrac = Lm/(fabs(Qcs)+Lm);
	    }
	    else {
		    subFrac = 0.0;
		    unlFrac = 1.0;
            }
	    Qcs -= I*subFrac; 
	    Lm += I*unlFrac;
	    I = 0.0;
    }

    // SKY2008Snow based on AJR2008's recommendation ends here
      //set adjusted fluxes and states to node
      // because precip is now scaled by coeffV these values now only reflect fluxes and stores in the canopy
      node->setIntSWE( naughttocm*( 1/rholiqkg )*I );
      node->setIntSnUnload(naughttocm*( 1/rholiqkg )*Lm);
      node->setIntSub( naughttocm*( 1/rholiqkg )*Qcs);
      node->addIntSub( naughttocm*( 1/rholiqkg )*Qcs);
      node->addIntUnl( naughttocm*( 1/rholiqkg )*Lm );
      node->setIntPrec(Isnow*( 1/rholiqkg )*naughttocm);
      // Rate for the _ENTIRE_ cell:
      node->setNetPrecipitation(throughfall + (1-coeffV)*node->getRain());
      // note mm and kg/m^2 requires not conversion

      //set wet and dry evap to 0 when snow in canopy
      node->setEvapWetCanopy(0.0);
      node->setEvapDryCanopy(0.0);

  }//end -- snow exists
  count++;

  return;
}


/****************************************************************************
**
**		      tSnowIntercept -- Physical Routines
**
**	Functions that compute changes internal to the canopy for 
**	tSnowIntercept::callSnowIntercept. A loading function should probably
**	be implemented in order to fully modulate the algorithm.
**
****************************************************************************/

//---------------------------------------------------------------------------
//
//		tSnowIntercept::computeSub()
//
//    Uses Liston and Sturm (2006) implementation to compute the snow sublimated
//    from the canopy after the snow has existed in the canopy for more than a
//    singe time step. All notation is taken from Liston and Sturm (2006).
//
//----------------------------------------------------------------------------

void tSnowIntercept::computeSub()
{

  //compute incoming shortwave radiation
  inShortR = inShortWaveSn();// W

  //compute effective incident shortwave radiation on snow crystal
  Sp = 3.1416*pow(iceRad,2.0)*(1 - 0.8)*inShortR;//check units--check (W)

  //Find coefficient for changing windspeed
  acoefficient = beta*coeffLAI;

  //find windspeed
  if (windSpeed == 0.0) {
      windSpeed = 0.1;
  }
  windSpeed = windSpeed*exp(-acoefficient*0.4);

  //Calculate Reynolds number
  Re = 2*iceRad*windSpeed/nu;

  //Calculate Sherwood number
  Sh = 1.79 + 0.606*pow(Re,0.5);

  //Calculate Nusselt number
  Nu = Sh;
  
  //calculate saturated vapor pressure of at ice interface
  esatIce = 611.15*exp( 22.452*( airTempK - 273.16)/(airTempK - 0.61)); //check units--check

  //calculate density of vapor
  rhoVap = 0.622*esatIce/(RdryAir*airTempK);
  
  //compute vapor diffusivity
  D = 2.06e-5 * pow(airTempK/273,1.75);
  
  //Place holder in algorithm
  Omega = (1/(KtAtm*airTempK*Nu))*(1000*latSubkJ*Mwater/(R*airTempK) - 1);//check units--check
  
  //find change of mass of ice crystal with respect to time
  dmdt = (2*3.1416*iceRad*(rHumidity/100 - 1) - Sp*Omega) /
	    (1000*latSubkJ*Omega + (1/(D*rhoVap*Sh)));
  
  //relative sublimation from ice sphere
  psiS = dmdt/( (4/3)*3.1416*rhoicekg*iceRad*iceRad*iceRad );

  //canopy exposure coefficient
  Ce = kc*pow(I/Imax,-0.4);

  //compute total sublimated snow during timestep
  Qcs = Ce*I*psiS*timeSteps;

  //RMK: Qcs IS AN INTERNAL VARIABLE TO THE CLASS SO WE DO NOT NEED TO 
  //	 RETURN IT TO THE CALLING FUNCTION.

  return;
}
  
//---------------------------------------------------------------------------
//
//			tSnowIntercept::computeUnload()
//
//	Compute the amount of unloading during a timestep according to 
//	Liston and Sturm (2006). This is basically a degree day approach.
//
//----------------------------------------------------------------------------

void tSnowIntercept::computeUnload()
{

  //find if over critical temperature
  if (airTempK >= 273.16) {
    Lm = 5.8e-5*(airTempK - 273.16)*timeSteps;//unload
  }
  else {
    Lm = 0.0;//do not unload
  }

  //RMK: Lm IS AN INTERNAL VARIABLE AND DOES NOT NEED TO BE RETURNED TO THE
  //	 CALLING FUNCTION.
  
  return;

}


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

//-----------------------------------------------------------------------------
//
//		      tSnowIntercept::snowFracCalc()
//
//	Calculates the proportions of solid and liquid precipitation as a
//	function of air temperature. See same function in tSnowPack for further
//	documentation.
//
//	RMK:  NOT USED IN CURRENT INTERCEPTION IMPLEMENTATION
//
//-----------------------------------------------------------------------------

double tSnowIntercept::snowFracCalc()
{

  double snowfrac;

  double TMin(-1.1), TMax(3.3); //indices (Wigmosta et al. 1994)
  
  if ( airTemp <= TMin )
	  snowfrac = 1;//all solid
  if ( airTemp >= TMax )
	  snowfrac = 0;//all liquid
  if ( (airTemp >= TMin)&&(airTemp <= TMax) )
	  snowfrac = (TMax - airTemp)/(TMax - TMin);//mix
  
  return snowfrac;
}

//---------------------------------------------------------------------------
//
//			tSnowIntercept::latentHFCalc()
//
//    Currently, this is a dummy function, but could be used in a later 
//    implementation where this can actually be estimated or calculated
//    differently than the ubiquitous Pomeroy (1998) papers.
//
//---------------------------------------------------------------------------

double tSnowIntercept::latentHFCalc(double Kaero)
{

  double lhf;

  lhf = 0.0;
  return lhf;
}

//-----------------------------------------------------------------------------
//
//			tSnowIntercept::latHeatVapCalc()
//	
//	Calculates the latent heat of vaporization. Can be made more specific 
//	if needed in later model implementations.
//
//-----------------------------------------------------------------------------

double tSnowIntercept::latHeatVapCalc()
{

  double lhvap(2470);// kJ/kg

  return lhvap;
}

//-----------------------------------------------------------------------------
//
//			tSnowIntercept::latHeatFreezeCalc()
//
//	Calculates the latent heat of freezing. Can be made more specific
//	if need in later model implementations.
//			
//
//-----------------------------------------------------------------------------

double tSnowIntercept::latHeatFreezeCalc()
{

  double lhfreeze(334);// kJ/kg

  return lhfreeze;
}

//------------------------------------------------------------------------------
//
//		      tSnowIntercept::latHeatSubCalc()
//
//	  Calculates the latent heat of sublimation. Can be made more specific
//	  if needed in later model implementations.
//
//------------------------------------------------------------------------------

double tSnowIntercept::latHeatSubCalc()
{

  double lhsub(2470 + 334);// kJ/kg

  return lhsub;
}

//------------------------------------------------------------------------------
//
//			tSnowIntercept::heatCapVapCalc()
//
//	  Calculates the heat capacity of water vapor. Can be made more specific
//	  if needed in later implementations.
//
//------------------------------------------------------------------------------

double tSnowIntercept::heatCapAirCalc()
{
  double heatcapvap(1.01); //kJ/(K*kg)

  return heatcapvap;
}

//-------------------------------------------------------------------------------
//
//			  tSnowIntercept::heatCapSolCalc()
//	
//	Calculates the heat capacity of solid water. Can be made more specific
//	if needed in later implementations.
//
//-------------------------------------------------------------------------------

double tSnowIntercept::heatCapSolCalc()
{
  double heatcapsol(2.1); //kJ/(K*kg)

  return heatcapsol;
}

//--------------------------------------------------------------------------------
//
//			  tSnowIntercept::heatCapLiqCalc()
//			  
//	Calculates the heat capacity of liquid water. Can be made more specific if
//	needed in later implementations.
//
//--------------------------------------------------------------------------------

double tSnowIntercept::heatCapLiqCalc()
{
  double heatcapliq(4.19); // kJ/(K*kg)

  return heatcapliq;
}

//-------------------------------------------------------------------------
//
// tSnowIntercept::inShortWaveSn() Function
//
//
//  This code is a modified version of tEvapoTrans::inShortWave() and
//  tSnowPack::inShortWaveSn() functions. It does not adjust for overlying
//  vegetation.
//
//  It does, however, account for the different sheltering options and the
//  change in albedo in a canopy covered in snow.
//
//  
//
// Calculate the Clear Sky Incoming Direct Solar Radiation Intensity Ic 
//                                               for horizontal surface
//       Ic = Io*t^(1/sin(alpha))           Wilson and Gallant (Eq 4.7)
//     Transmission coefficient (t)
//        t = 0.65 + 0.00008*Z      Z = node elevation (meters) (Eq 4.8)
//     
//     Adjustment for Circumsolar Diffuse Radiation
//       Ic = Ic + Id*CIRC                    W&G Eq 4.17
//
// Calculate the Clear Sky Incoming Diffusive Radiation Intensity Id   
//                                            for horizontal surface 
//       Id = (0.271 - 0.294*t^(1/sin(alpha)))*Io
//     Adjustment for Circumsolar Diffuse Radiation
//       Id = Id - Id*CIRC                    W&G Eq 4.18
//       CIRC = {.07, .10, .12, .14, .16, .19, .23, .20, .16, .12, .08, .07}
//           for each month Jan - Dec.
// 
// Calculate the Direct Solar Radiation on Sloping Surface
//
//       Ics = Ic*cosi                        W&G Eq 4.20 - 4.24
//       cosi = A + B*cos(tau) + C*sin(tau)
//       A = sin(del)*sin(phi)*cos(beta) + sin(beta)*cos(asp)*cos(phi)
//       B = cos(del)*(cos(phi)*cos(beta) - sin(phi)*cos(asp)*sin(beta))
//       C = sin(beta)*cos(del)*sin(asp)
//
//       where beta = slope angle (radians), asp = aspect angle (radians)
//
// Calculate the Diffuse Solar Radiation on Sloping Surface
//
//       Ids = Id * v
//       v = sky view fraction ~ cos2(beta/2)         Moore et al (1993)
//
// Calculate Reflection Radiation Component
//
//       Ir = (Ic + Id)*(1 - v)*A      A = albedo     W&G Eq 4.26
//
// Total Incoming Solar Radiation on Sloping Surface (Clear Sky)
//
//       Ith = (Ics + Ids) + Ir
//
// Calculate the Cloudy Sky Incoming Solar Radiation Intensity Is
//       Is = (1-0.65*N^2)*(Ics + Ids) + Ir                              2.29
//     Fractional Sky Cover N 
//          N = skycover/10;
//
// Calculate the Effect of Vegetation Absorption on Incoming Solar Radiation
//       Iv = Kt*Is;
//     Optical Transmission Coefficient of Surface/Vegetation Kt
//     Landuse-derived coefficient coeffKt    Range (0-1)  bareground = 1;
//
// Calculate the Albedo Effect on Incoming Solar Radiation
//       Isw = Iv(1-A)                                                   2.32
//     Albedo Coefficient of Surface/Vegetation
//     Landuse-derived coefficient coeffA     Range (0-1) See Bras(1990) p.37
//
//----------------------------------------------------------------------------

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

int tSnowIntercept::getSnowOpt()
{
  return snowOption;
}

/***************************************************************************
**
** tSnowIntercept::writeRestart() Function
**
** Called from tSimulator during simulation loop
**        
***************************************************************************/
void tSnowIntercept::writeRestart(fstream & rStr) const
{  
  BinaryWrite(rStr, nID);
  BinaryWrite(rStr, hillAlbedoOption);
  BinaryWrite(rStr, Qcs);
  BinaryWrite(rStr, Ce);
  BinaryWrite(rStr, I);
  BinaryWrite(rStr, Iold);
  BinaryWrite(rStr, psiS);
  BinaryWrite(rStr, Imax);
  BinaryWrite(rStr, prec);
  BinaryWrite(rStr, LAI);
  BinaryWrite(rStr, kc);
  BinaryWrite(rStr, iceRad);
  BinaryWrite(rStr, dmdt);
  BinaryWrite(rStr, Omega);
  BinaryWrite(rStr, Sp);
  BinaryWrite(rStr, RH);
  BinaryWrite(rStr, D);
  BinaryWrite(rStr, rhoVap);
  BinaryWrite(rStr, Sh);
  BinaryWrite(rStr, Nu);
  BinaryWrite(rStr, Re); 
  BinaryWrite(rStr, KtAtm);
  BinaryWrite(rStr, Ta);
  BinaryWrite(rStr, Mwater);
  BinaryWrite(rStr, R);
  BinaryWrite(rStr, RdryAir);
  BinaryWrite(rStr, esatIce);
  BinaryWrite(rStr, nu);
  BinaryWrite(rStr, beta);
  BinaryWrite(rStr, acoefficient);
  BinaryWrite(rStr, Lm);
  BinaryWrite(rStr, airTempK);
  BinaryWrite(rStr, effPrecip);

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
  BinaryRead(rStr, nID);
  BinaryRead(rStr, hillAlbedoOption);
  BinaryRead(rStr, Qcs);
  BinaryRead(rStr, Ce);
  BinaryRead(rStr, I);
  BinaryRead(rStr, Iold);
  BinaryRead(rStr, psiS);
  BinaryRead(rStr, Imax);
  BinaryRead(rStr, prec);
  BinaryRead(rStr, LAI);
  BinaryRead(rStr, kc);
  BinaryRead(rStr, iceRad);
  BinaryRead(rStr, dmdt);
  BinaryRead(rStr, Omega);
  BinaryRead(rStr, Sp);
  BinaryRead(rStr, RH);
  BinaryRead(rStr, D);
  BinaryRead(rStr, rhoVap);
  BinaryRead(rStr, Sh);
  BinaryRead(rStr, Nu);
  BinaryRead(rStr, Re); 
  BinaryRead(rStr, KtAtm);
  BinaryRead(rStr, Ta);
  BinaryRead(rStr, Mwater);
  BinaryRead(rStr, R);
  BinaryRead(rStr, RdryAir);
  BinaryRead(rStr, esatIce);
  BinaryRead(rStr, nu);
  BinaryRead(rStr, beta);
  BinaryRead(rStr, acoefficient);
  BinaryRead(rStr, Lm);
  BinaryRead(rStr, airTempK);
  BinaryRead(rStr, effPrecip);

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
