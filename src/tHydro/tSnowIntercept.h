/***************************************************************************
**
**                   tRIBS Distributed Hydrology Model
**
**              TIN-based Real-time Integrated Basin Simulator
**                       Ralph M. Parsons Laboratory
**                  Massachusetts Institute of Technology
**  
**
**  tSnowIntercept.h:   Header file for tSnowIntercept Class
**
**  This class creates an object that handles the snow-in-canopy, i.e.
**  snow interception. It also is derived from tEvapoTrans, b/c of the
**  importance of radiation and the rest of the energy balance on the 
**  state of the canopy. The functional structure is almost identical to 
**  that of tSnowPack. It consists of the constructor/destructor, initialization
**  functions, calling function, interaction w/ the discretization, basic
**  physics, conversions, the physics necessary for an energy balance,
**  an energy balance function, and basic I/O.
**
**  The variables fall into a similar set of classifications. However, the
**  basic set-up is taken directly from Liston and Elder (2006). These are
**  the primary variables used in the current version of the code. Most of
**  the additional variables are for an implementation of Gelfan et al (2004),
**  which claims to have a canopy energy balance. In the future, it would
**  be a good exercise to actually compare the two implementations.
**
**  The object is called from tSnowPack::callSnowPack() only when there is
**  (a) interception is turned on, and (b) there is canopy. At this point,
**  the node information is retrieved and the physical algorithm is entered
**  into. If it is snowing out OR if there is snow in the canopy, then the 
**  snow interception algorithm is called. Otherwise,
**  tIntercept::callIntercept() is recreated, similar to the recreation of
**  tEvapoTrans::callEvapoPotential() and tEvapoTrans::callEvapoTrans() in
**  tSnowPack::callSnowPack().
**
**  If the snow interception is called, first, if their is snowfall, it is
**  intercepted in the canopy according to an adjusted version of the
**  Pomeroy (1998) algorithm. It has been adjusted to allow the sequential
**  intercpetion of snow in a "wet" (snow is already in) canopy. It is 
**  assumed that there is no unloading or sublimation during this part of
**  the period.
**
**  If the snow in the canopy is from a previos time step, then sublimation
**  and unloading can occur. The sublimation has been adapted from the
**  Pomeroy (1998) algorithm in the Liston and Sturm (2006) implementation.
**  Unloading is calculated using a linear temperature threshold, similar to
**  degree day approaches in melt modeling.
**
**  The effective rainfall is adjusted during the interception and unloading
**  portions of the algorithm and is called again in tSnowPack::callSnowPack().
**
**  11 July 2007 -- Rinehart @ New Mexico Tech
**
**			Pertinent References: Pomeroy et al (1998) -- this is a sequence of papers
**					      Gelfan et al (2004)
**					      Liston and Sturm (2007)
**
**************************************************************************/

#ifndef TSNOWINTERCEPT_H
#define TSNOWINTERCEPT_H

//=========================================================================
//
//		Section 1: tSnowIntercept Include Statements
//
//=========================================================================

#include "src/Headers/Inclusions.h"
#include "src/tHydro/tEvapoTrans.h"



//=========================================================================
//
//		Section 2: tSnowIntercept Class Definitions
//
//=========================================================================

class tSnowIntercept : public tEvapoTrans
{
public:
  tSnowIntercept();
  tSnowIntercept(SimulationControl *,tMesh<tCNode> *, tInputFile &, tRunTimer *,
		  tResample *, tHydroModel *, tRainfall *);
  ~tSnowIntercept();

  //initialization routine
  void SetSnowInterceptVariables(tInputFile &, tHydroModel *);
  void SetSnowVariables(tInputFile &, tHydroModel *);  

  //calling function
  void callSnowIntercept(tCNode *, tIntercept *);
  
  //general physical functions
  void computeSub();
  void computeUnload();
 
  //conversion functions
  double CtoK(double);
  double KtoC(double);
  
  //EB functions -- basic calculations
  double latentHFCalc(double);
  double snowFracCalc();
  double latHeatVapCalc();
  double latHeatFreezeCalc();
  double latHeatSubCalc();
  double heatCapAirCalc();
  double heatCapSolCalc();
  double heatCapLiqCalc();
  double inShortWaveSn();

  //EB function
  void snowEB(int);
  
  //communication functions
  int getSnowOpt();
  
  // Restart functions
  void writeRestart(fstream &) const;
  void readRestart(fstream &);

protected:

  //-------------------
  //intercept variables (from Liston and Elder 2006, section 3)
  int nID;
  int hillAlbedoOption;
  double Qcs, Ce, I, Iold, psiS; //
  double Imax, prec, LAI; //
  double kc, iceRad, dmdt; //
  double Omega, Sp, RH, D, rhoVap; //
  double Sh, Nu, Re; // 
  double KtAtm, Ta, Mwater, R; //
  double RdryAir, esatIce, nu, beta; //
  double acoefficient; //
  double Lm; // unloading
  double airTempK; //
  double effPrecip; //

  //-----------------
  //general variables

  //discretization
  double timeSteph, timeSteps, timeStepm;
  double minutelyTimeStep;
  
  //state variables (intrinsic)
  double liqWE, iceWE, snWE; //cm
  double liqRoute; //cm
  double liqWEm, iceWEm, snWEm; //m
  double liqRoutem; //m
  double Utot, Usn, Uwat; //internal energy (kJ/m^2), set to 0 at T=0 C
  double liqTempC, iceTempC, snTempC;
  double liqTempK, iceTempK, snTempK; //Kelvin
  double crustAge; //hrs
  double albedo;
  double hillalbedo;
  
  //energy/mass fluxes
  double H,L,G,Prec,Rn; //components of energy balance (kJ/(m^2)(s))  
  double snPrec, liqPrec; //cm
  double snPrecm, liqPrecm; //m
  double snPrecmm, liqPrecmm; //mm
  double vapPressSmb, vapPresskSPa; //vapor pressure (mb and Pa)
  double rholiqcgs, rhoicecgs, rhosncgs; //g/cm^3
  double rholiqkg, rhoicekg, rhosnkg; //kg/m^3
  double rhoAir; //kg/m^3
  
  //thermal properties
  double cpsnowkJ, cpicekJ, cpwaterkJ, cpairkJ; //heat capacity kJ/Kg
  double latFreezekJ, latVapkJ, latSubkJ; //latent heats kJ/g
  double resFact;

  //output variables
  double snDepth,snDepthm; //snow depths (cm,m)
  
  //conversion factors
  double naughttokilo, kilotonaught, cgsRHOtomks, mksRHOtocgs;
  double naughttocm, cmtonaught, ctom, mtoc;
  double htos;


};

#endif

//==========================================================================
//
//		End of tSnowIntercept.h
//
//==========================================================================
