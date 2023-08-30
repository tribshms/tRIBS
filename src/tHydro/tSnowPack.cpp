/***************************************************************************
**
**  		     tRIBS Distributed Hydrologic Model
**
**              TIN-based Real-time Integrated Basin Simulator
**		         Ralph M. Parsons Laboratory
**  		    Massachusetts Institute of Technology
**
**
**  tSnowPack.cpp:   Function file for tSnowPack class (see tSnowPack.h)
**
***************************************************************************/

#include "src/tHydro/tSnowPack.h"
#include "src/Headers/globalIO.h"

//===========================================================================
//
//		Section 1: Constructor and Initialization Routines
//
//===========================================================================

//---------------------------------------------------------------------------
//
//  			tSnowPack() Constructor and Destructor
//
//  Construct the tSnowPack object. This class inherits from tEvapoTrans but
//  does not require anymore inputed parameters than tEvapoTrans. Additional
//  parameters are initialized using SetSnowVariables and SetSnowPack variables.
//  The logic in having two set-functions is to allow the easy construction of
//  other snow-type classes that require similar information as all other snow
//  classes (e.g., thermal properties, energy balance information, SWE, ...).
//
//---------------------------------------------------------------------------

tSnowPack::tSnowPack() {

}

tSnowPack::tSnowPack(SimulationControl *simCtrPtr, tMesh<tCNode> *gridRef,
                     tInputFile &infile, tRunTimer *t, tResample *resamp,
                     tHydroModel *hydro, tRainfall *storm)
        : tEvapoTrans(simCtrPtr, gridRef, infile, t, resamp, hydro, storm) {

    gridPtr = gridRef;
    //timerET = t;
    timer = t; // SKY2008Snow
    simCtrl = simCtrPtr;

    //set variables
    SetSnowVariables(infile, hydro);
    SetSnowPackVariables(infile, hydro);

}

tSnowPack::~tSnowPack() {
    Cout << "tSnowPack Object has been destroyed..." << endl;
}

//---------------------------------------------------------------------------
//
//	tSnowPack::SetSnowPackVariables()
//
//	Auxiliary function used in robust constructor to initialize snow 
//	variables
//
//---------------------------------------------------------------------------

void tSnowPack::SetSnowPackVariables(tInputFile &infile, tHydroModel *hydro) {

    //parameters
    minSnTemp = infile.ReadItem(minSnTemp, "MINSNTEMP");
    snliqfrac = infile.ReadItem(snliqfrac, "SNLIQFRAC"); // Added by CJC 2020
    hillAlbedoOption = infile.ReadItem(hillAlbedoOption, "HILLALBOPT");
    densityAge = 0.0;
    ETAge = 0.0;
    compactParam = 0.3;
    rhoSnFreshkg = 100;
    snOnOff = 0.0;

    return;
}


//---------------------------------------------------------------------------
//
//	tSnowPack::SetSnowVariables()
//
//	Auxiliary function used in robust constructor to initialize snow 
//	variables. This function will remain more or less the same across
//	snow classes.
//
//---------------------------------------------------------------------------

void tSnowPack::SetSnowVariables(tInputFile &infile, tHydroModel *hydro) {

    //time steps
    timeStepm = infile.ReadItem(timeStepm, "METSTEP");
    timeSteph = timeStepm / 60;
    timeSteps = 60 * timeStepm;
    minutelyTimeStep = 0.0;

    //state variables
    liqWE = iceWE = snWE = 0.0;
    liqWEm = iceWEm = snWEm = 0.0;
    Utot = Usn = Uwat = 0.0;
    liqWatCont = 0.0;
    liqTempC = iceTempC = snTempC = 0.0;
    liqTempK = iceTempK = snTempK = 0.0;
    crustAge = 0.0;
    albedo = 0.88;
    canWE = 0.0;

    //fluxes
    H = L = G = Prec = Rn = 0.0;
    snPrec = liqPrec = 0.0;
    snPrecm = liqPrecm = 0.0;
    snPrecmm = liqPrecmm = 0.0;
    vapPressSmb = vapPresskSPa = 0.0;

    //density
    rholiqcgs = 1.0;
    rhoicecgs = 0.92;
    rhosncgs = 0.1;
    rholiqkg = 1000.0;
    rhoicekg = 920.0;
    rhosnkg = 100.0;
    rhoAir = 1.3;

    //thermal properties
    cpicekJ = 2.1;
    cpwaterkJ = 4.190;
    cpairkJ = 1.006;
    latFreezekJ = 334;
    latVapkJ = 2470;
    latSubkJ = latFreezekJ + latVapkJ;

    //output variables
    snDepth = snDepthm = 0.0;
    snOnOff = 0.0;
    peakSnWE = peakSnWEtemp = 0.0;
    persMax = persMaxtemp = 0.0;
    inittime = peaktime = 0.0;

    //conversions
    naughttokilo = 1e-3;
    kilotonaught = 1e3;
    cgsRHOtomks = 1e3;
    mksRHOtocgs = 1e-3;
    naughttocm = 100; // Used convert m to cm
    cmtonaught = 0.01; // Used convert cm to m
    ctom = 10;
    mtoc = 0.1;

    //canopy conditions

    return;
}

//---------------------------------------------------------------------------
//
//			tSnowPack::callSnowPack()
//
//  Calling function for snow pack dymamics. This function has the dual role
//  of representing both snow physics and evapotranspiration in the warm 
//  landscape. This is because we have to check to see if there is snow either
//  on the ground or coming down from the canopy or atmosphere before we can
//  tell if if we need to deal with snow or not.
//
//  The structure of this function is similar to that of 
//  tEvapoTrans::callEvapoPotential and tEvapoTrans::callEvapoTrans. It begins
//  by setting the sun variables and preparing the stochastic weather 
//  simulator if necessary. Then, the mesh list is initialized and the loop
//  through the node list begins.
//
//  At this point, the slope, aspect and elevation of every node is computed,
//  the appropriate sheltering algorithm is implemented, the appropriate met
//  and rainfall data is called, the land-surface is reinitialized, and the 
//  current state of the land-surface is found (snow or no-snow). If lapse
//  rates are implemented, then the air temperature and/or the rainfall is 
//  adjusted for elevation.
//
//  At this point, the actual physics begin. If the interception scheme is
//  on and there is canopy, then tSnowIntercept::callSnowIntercept() is
//  called and the snow-canopy interaction is accounted for.
//
//  Then existence of snow pack, snow precipitation or snow throughfall is 
//  checked for. If it does not exist, then tEvapoTrans::callEvapoPotential 
//  and tEvapoTrans::callEvapoTrans are recreated. RMK: THIS ASSUMES THAT 
//  ETISTEP AND METSTEP ARE THE SAME, WHICH HAS BEEN DONE IMPLICITLY THROUGH
//  MOST OF THE EXISTING MODEL.
//
//  If it does exist, or if there is a positive mass flux, we enter the actual
//  snow-physics portion of the code. Here, we begin by dealing with the 
//  atmospheric snow mass balance (adding/subtracting precip, throughfall,
//  turbulent latent heat flux). We have to make the distinction between a 
//  developed and undeveloped snow pack in order to accurately represent the
//  heat flux. If there is no snow on the ground, then there is, in a sense,
//  no precipitation heat flux. The pack's state is initialized to that of the
//  air temperature and precipitation phase masses. Otherwise, we deal with 
//  both the mass and energy precipitation fluxes, and initialize the energy
//  of the pack using the updated mass and the previous temperature.
//
//  RMK: IF SOMEONE REWORKS THIS PORTION OF THE MODEL, RESEQUENCING THE MB
//  AND EB WOULD BE A GOOD IDEA, IN ORDER TO ELIMINATE SOME OF THE ERROR THAT
//  WILL COME FROM THE CURRENT SEQUENCING. UNFORTUNATELY, I COULDN'T FIGURE
//  OUT A WAY TO HAVE THE EB OCCUR DURING THE ACCUMULATION PERIOD WITHOUT THIS
//  SEQUENCING.
//
//  We are now discussing the energy balance portion of the code. If there is
//  no snow left after the initial mass balance, then we do not actually 
//  compute the energy balance. This is to ensure stability of temperatures.
//
//  If there is snow left, then we initialize the internal energy of the pack,
//  compute the change in energy, add the change in energy, and compute the
//  new state of the model. If the internal energy is greater than 0J/m^2, then
//  there is liquid water in the pack. Otherwise, we compute the temperature
//  of the snow pack (which will be less than 0).
//
//  RMK: WE ASSUME A SINGLE EQUIVALENT TEMPERATURE IS REPRESENTATIVE, BUT WE
//  KNOW THAT IT IS NOT.
//
//  If there the amount of liquid mass is greater than 40% of the solid mass,
//  then the excess water is routed out of the pack and is treated as 
//  precipitation in tHydroModel::UnSaturatedZone(). That's it. It seems so
//  simple now that its done.
//
//  09 July 2007 -- AJR @ New Mexico Tech
//
//
//				    Pertinent references: Wigmosta et al (1994)
//							  Tarboton et al (?)
//							  Tuteja et al (1996)
//							  Marks et al (?)
//							  Anderson (1976)
//							  Jordan (1991)
//				    
//---------------------------------------------------------------------------
// SKY2008Snow, AJR2008
void tSnowPack::callSnowPack(tIntercept *Intercept, int flag, tSnowIntercept *SnIntercept)
{

    tCNode *cNode;
    tMeshListIter<tCNode> nodeIter(gridPtr->getNodeList());

    int count = 0;
    int cnt = 0;
    double EP = 0.0; //double tmp = 0.0;
    double SkyC = 0.0; //double tmpC = 0.0;
    double vegHeight;

    // SKY2008Snow, AJR2008
    //  metHour = metStep;
    //  etHour = etStep;

    if (simCtrl->Verbose_label == 'Y') {
        cout << "\nSnowPack Routine Call ..." << endl;
    }

    // Set time, sun, and meteorological variables -- AJR2008, SKY2008Snow
    SetEnvironment();

    // Resample Meteorological Grids, if option
    if (metdataOption == 2) {
        resampleGrids(timer);
    }

    if (luOption == 1) { // resampling Land Use grids done here, i.e., dynamic case
        if (AtFirstTimeStepLUFlag) {
            initialLUGridAssignment();
            AtFirstTimeStepLUFlag=0;
        }
        else {
            LUGridAssignment();
        }
    }

    // BEGIN LOOP THROUGH NODES
    cNode = nodeIter.FirstP();
    while (nodeIter.IsActive()) {
        // Set static land use tCnode members before being assigned dynamically further below.
        landPtr->setLandPtr(cNode->getLandUse());
        cNode->setCanStorParam(landPtr->getLandProp(1));
        cNode->setIntercepCoeff(landPtr->getLandProp(2));
        cNode->setThroughFall(landPtr->getLandProp(3));
        cNode->setCanFieldCap(landPtr->getLandProp(4));
        cNode->setDrainCoeff(landPtr->getLandProp(5));
        cNode->setDrainExpPar(landPtr->getLandProp(6));
        cNode->setLandUseAlb(landPtr->getLandProp(7));
        cNode->setVegHeight(landPtr->getLandProp(8));
        cNode->setOptTransmCoeff(landPtr->getLandProp(9));
        cNode->setStomRes(landPtr->getLandProp(10));
        cNode->setVegFraction(landPtr->getLandProp(11));
        cNode->setLeafAI(landPtr->getLandProp(12));
        cNode = nodeIter.NextP();

        if (luOption == 1) {
            if ( luInterpOption == 1) { // LU values linearly interpolated between 'previous' and 'until' values
                interpolateLUGrids(cNode);
            }
        }

        // Elapsed MET steps from the beginning, used for averaging dynamic LU grid values below over time for integ. output
        double te = (double)timer->getElapsedMETSteps(timer->getCurrentTime());
        integratedLUVars(cNode, te);

        //Get NodeID
        ID = cNode->getID();

        //Get Rainfall
        rain = cNode->getRain(); // get new rainfall

        //Set Elevation, Slope and Aspect
        slope = fabs(atan(cNode->getFlowEdg()->getSlope()));
        aspect = cNode->getAspect();
        elevation = cNode->getZ();

        snOnOff = 0.0;

        checkShelter(cNode);

        setCoeffs(cNode);
        if (luOption == 1) {
            newLUGridData(cNode);
        }

        //not in stochastic mode
        if (!rainPtr->getoptStorm()) {
            if (metdataOption == 1) {
                thisStation = assignedStation[count];
                newHydroMetData(hourlyTimeStep); //read in met data from station file -- inherited function
            } else if (metdataOption == 2) {
                //resampleGrids(timerET); // read in met grid data -- inherited function
                newHydroMetGridData(cNode); // set up and get appropriate data -- inherited function
            }

            // Set the observed values to the node:
            // they will be required by other function calls
            vPress = vaporPress(); //-- ADDED IN ORDER TO SET RH... CORRECTLY FOR SNOW
            cNode->setAirTemp(airTemp); // celsius
            cNode->setDewTemp(dewTemp);
            cNode->setRelHumid(rHumidity);
            cNode->setVapPressure(vPress);
            cNode->setSkyCover(skyCover);
            cNode->setWindSpeed(windSpeed);
            cNode->setAirPressure(atmPress);
            cNode->setShortRadIn(RadGlbObs);

            //Set Soil/Surface Temperature
            if (hourlyTimeStep == 0) {
                cNode->setSoilTemp(Tlo - 273.15);
                cNode->setSurfTemp(Tso - 273.15);
            }

        }

        if (Ioption == 0) {
            cNode->setNetPrecipitation(rain);
        }

        //Call Beta functions
        betaFunc(cNode); // inherited
        betaFuncT(cNode); // inherited

        //get the necessary information from tCNode for snow model
        getFrNodeSnP(cNode);

        // ensure routed liquid is reset
        cNode->setLiqRouted(0.0);

        snUnload = 0.0;
        canWE = cNode->getIntSWE();

        //No Snow on ground and canopy and not snowing
        if ((snWE <= 1e-4) && (rain * snowFracCalc() <= 5e-2) && rholiqkg * cmtonaught * (cNode->getIntSWE()) < 1e-3) {

            // Following block of code mirrors callEvapoPotential and callEvapoTrans in tEvapoTrans as no snow occurs at any
            // level of the system: i.e. snowpack, canopy, or snowing. —refactored by WR 6/21/23

            Tso = cNode->getSurfTemp() + 273.15;
            Tlo = cNode->getSoilTemp() + 273.15;

            //Calculate the Potential and Actual Evaporation
            if (evapotransOption == 1) {
                EvapPenmanMonteith(cNode); // SKY2008Snow
            } else if (evapotransOption == 2) {
                EvapDeardorff(cNode); // SKY2008Snow
            } else if (evapotransOption == 3) {
                EvapPriestlyTaylor(cNode); // SKY2008Snow
            } else if (evapotransOption == 4) {
                EvapPan();
            } else {
                cout << "\nEvapotranspiration Option " << evapotransOption;
                cout << " not valid." << endl;
                cout << "\tPlease use :" << endl;
                cout << "\t\t(1) for Penman-Monteith Method" << endl;
                cout << "\t\t(2) for Deardorff Method" << endl;
                cout << "\t\t(3) for Priestly-Taylor Method" << endl;
                cout << "\t\t(4) for Pan Evaporation Measurements" << endl;
                cout << "Exiting Program...\n\n" << endl;
                exit(1);
            }
            // Set
            setToNode(cNode);
            ComputeETComponents(Intercept, cNode, count, flag);

            snTempC = 0.0; // reinitialize snTemp
            snWE = 0.0; // reinitialize snWE
            snSub = 0.0; // No sublimation occurs CJC2020
            snEvap = 0.0; // No evaporation occurs CJC2020
            dUint = RLin = RLout = RSin = H = L = G = Prec = 0.0; //reinitialize energy terms
            ETAge = ETAge + timeStepm;
            liqRoute = 0.0;
            iceWE = 0.0;
            liqWE = 0.0;
            Utot = 0.0;
            Usn = 0.0;
            Uwat = 0.0;
            snTempC = 0.0;
            crustAge = 0.0;
            densityAge = 0.0;

            cNode->setIntSub(0.0);
            setToNodeSnP(cNode);
        }//end no-snow

        else //condtions include some combination of snowpack and snow in canopy, and snowing, raining, or no precip
        {

            // Implement interception schemes for snow, refactored WR 6/21/23
            if (Ioption && (Intercept->IsThereCanopy(cNode))) {
                SnIntercept->callSnowIntercept(cNode, Intercept);
                snUnload = cNode->getIntSnUnload(); //calculated in callSnowIntercept() units in cm
                snCanWE = cNode->getIntSWE();//units in cm
            } else {
                cNode->setNetPrecipitation(rain);
            }

            rain = cNode->getNetPrecipitation(); //units in mm
            rain += snUnload * ctom; // units in mm

            // Here rain is being set by net precipitation which is set from callSnowIntercept
            // and represents throughfall + unloaded snow  from snow interception (so scaled by coeffV) and precipitation
            // that falls on the no vegetated fraction of the cell (1-coeffV). refactored WR 6/21/23

            snDepthm = cmtonaught * snWE / 0.312; // 0.312 is value for bulk density of snow, Sturm et al. 2010
            //calculate current snow depth for use in the turbulent heat flux calculations and output.

            if (snWE < 1e-5) {

                //no precipitation heat flux, as it is totally accounted for in the snow pack energy
                //   initialization
                phfOnOff = 0.0;

                //account for veg height
                if (coeffH == 0) {
                    vegHeight = 0.1;
                } else {
                    vegHeight = coeffH;
                }

                //set the new density age
                densityAge = 0.0;

                //reinitialize crust age
                crustAge = 0.0;

                //change mass (volume) quantities to correct units (kJ, m, C, s)
                iceWE = iceWE * cmtonaught;
                liqWE = liqWE * cmtonaught;
                snWE = iceWE + liqWE;

                if (airTemp > 0.0) {
                    snTempC = 0.0;
                } else {
                    snTempC = airTemp;
                }

                //snowMB
                // evaporate liquid from ripe pack snWE +=
                // note evaporation/sublimation flux can be negative or positive, with latter representing condensation/deposition
                if (snTempC == 0.0) {

                    //liq WE update
                    snEvap = (1 - coeffV) * latentHFCalc(resFactCalc()) * timeSteps /
                             (rholiqkg * latVapkJ);//*naughttocm; // Calculate evaporation from snowpack in cm CJC2020
                    liqWE += cmtonaught * ((mtoc * rain) * (1 - snowFracCalc())) * timeSteps /
                             3600; // Removed snUnload term CJC2020

                    if (liqWE + snEvap <= 0) {
                        snEvap = liqWE;
                        liqWE = 0;
                    } else {
                        liqWE += snEvap;
                    }

                    //solid WE update
                    iceWE += cmtonaught * (mtoc * (rain * snowFracCalc())) * timeSteps / 3600;
                    snSub = 0.0; // No sublimation occurs CJC2020
                }

                    //sublimate solid from frozen pack
                else {

                    //liq WE update
                    liqWE += cmtonaught * (mtoc * (rain * (1 - snowFracCalc()))) * timeSteps /
                             3600; // Removed snUnload term CJC2020
                    snEvap = 0.0; // No evaporation occurs CJC2020

                    //ice WE update
                    snSub = (1 - coeffV) * latentHFCalc(resFactCalc()) * timeSteps /
                            (rholiqkg * latSubkJ);//*naughttocm; // Calculate sublimation from snowpack in cm CJC2020
                    iceWE += cmtonaught * (mtoc * (rain * snowFracCalc())) * timeSteps / 3600;

                    if (iceWE + snSub <= 0) {
                        snSub = iceWE;
                        iceWE = 0;
                    } else {
                        iceWE += snSub;
                    }
                }

                snWE = iceWE + liqWE;
                snSub *= naughttocm;
                snEvap *= naughttocm;
                //set other fluxes
                L = H = G = Prec = Utotold = 0.0;

                //initialize and record energy balance
                Utot = dUint = iceWE * rhoicekg * cpicekJ * snTempC // Changed to use rhoicekg CJC 2020
                               + liqWE * rholiqkg * latFreezekJ;

            } else {

                //account for precipitation heat flux
                phfOnOff = 1.0;

                //account for veg height
                if (coeffH == 0) {
                    vegHeight = 0.1;
                } else {
                    vegHeight = coeffH; // unused?
                }

                //find the new density age
                densityAge = (snWE * densityAge) / (snWE + mtoc * rain);

                //reset crust age if snowing out
                if (rain * snowFracCalc() > 1e-3) {
                    crustAge = 0.0;
                }

                //change mass (volume) quantities to correct units (kJs)
                iceWE = iceWE * cmtonaught;
                liqWE = liqWE * cmtonaught;
                snWE = iceWE + liqWE;

                //snowMB

                //ripe pack -- evaporate water
                if (snTempC == 0.0) {

                    //liq WE update
                    snEvap = (1 - coeffV) * latentHFCalc(resFactCalc()) * timeSteps /
                             (rholiqkg * latVapkJ);//*naughttocm; // Calculate evaporation from snowpack in cm CJC2020
                    liqWE += cmtonaught * ((mtoc * rain) * (1 - snowFracCalc())) * timeSteps /
                             3600; // Removed snUnload term CJC2020

                    if (liqWE + snEvap <= 0) {
                        snEvap = liqWE;
                        liqWE = 0;
                    } else {
                        liqWE += snEvap;
                    }

                    //ice WE update
                    iceWE += cmtonaught * (mtoc * (rain * snowFracCalc())) * timeSteps / 3600;
                    snSub = 0.0; // No sublimation occurs CJC2020
                }

                    //frozen pack -- sublimate water
                else {

                    //liq WE update
                    liqWE += cmtonaught * (mtoc * (rain * (1 - snowFracCalc()))) * timeSteps /
                             3600; // Removed snUnload term CJC2020
                    snEvap = 0.0; // No evaporation occurs CJC2020

                    //ice WE update
                    snSub = (1 - coeffV) * latentHFCalc(resFactCalc()) * timeSteps /
                            (rholiqkg * latSubkJ);//*naughttocm; // Calculate sublimation from snowpack in cm CJC2020
                    iceWE += cmtonaught * (mtoc * (rain * snowFracCalc())) * timeSteps / 3600;

                    if (iceWE + snSub <= 0) {
                        snSub = iceWE;
                        iceWE = 0;
                    } else {
                        iceWE += snSub;
                    }
                }

                snWE = iceWE + liqWE;
                snSub *= naughttocm;
                snEvap *= naughttocm;

                //snowEB
                ETAge = 0.0;

                L = latentHFCalc(resFactCalc());
                Prec = precipitationHFCalc();

                //if there is no snow left at this point, then bail out of
                //energy balance.
                if ((snWE <= 5e-6) || (snTempC < -800)) {
                    liqRoute = 0.0;
                    snWE = 0.0;
                    iceWE = 0.0;
                    liqWE = 0.0;
                    Utot = 0.0;
                    Usn = 0.0;
                    Uwat = 0.0;
                    snTempC = 0.0;
                    crustAge = 0.0;
                    densityAge = 0.0;
                } else {

                    //find initial state of energy
                    Utot = Utotold = iceWE * rhoicekg * cpicekJ * snTempC + liqWE * rholiqkg *
                                                                            latFreezekJ; // I am pretty sure this should be rhoicekg CJC 2020
                    //adjust albedo for age
                    albedo = agingAlbedo();

                    //calculate dU
                    snowEB(ID, cNode); // AJR2008, SKY2008Snow

                    //check for balance
                    Uerr = (Utot - Utotold) - dUint;

                    if (Utot < 0.0) { //frozen pack -- change temperature
                        Usn = Utot;// all energy in solid phase
                        Uwat = 0.0;// no energy in liquid phase
                        liqWatCont = 0.0;// no liquid content
                        liqWE = 0.0;// no liq WE
                        liqTempC = 0.0; // reset liq temp to default

                        iceWE = snWE;

                        //calculate sn temperature, modified THM 2012
                        if (iceWE < 0.1 && iceWE > 0) {
                            iceTempC = Usn / (cpicekJ * rhoicekg * 0.1);
                        } else {
                            iceTempC = Usn / (cpicekJ * rhoicekg * snWE);
                        }

                        //adjust to minimum snow temperature
                        //	RMK: THIS IS A KLUGE THAT IS NECESSARY B/C OF THE
                        //	ONE-LAYER ASSUMPTION. IT IS ALSO REQUIRED B/C OF THE
                        //	SIMPLISTIC WAY WE MODEL GROUND HEAT FLUX. SOMEONE
                        //	NEEDS TO INCORPORATE MULTIPLE LAYERS.
                        if (iceTempC <= minSnTemp) {
                            iceTempC = minSnTemp;
                        }

                        //set pack temperature to ice temperature
                        snTempC = iceTempC;

                    }//end -- frozen pack

                    else {//melt
                        Uwat = Utot;
                        Usn = 0.0;


                        liqWE = Uwat / (latFreezekJ * rholiqkg); //THM // The only thing THM change was = to +=

                        //make sure that there is enough SWE in the pack for the melt
                        if (liqWE >= snWE) {
                            liqWE = snWE; // this is here because the liqWE += term above
                        }                  //  can result in liqWE > snWE
                        //assign water equivalents
                        iceWE = snWE - liqWE;

                        //put in routing bucket
                        if (liqWE > snliqfrac * iceWE) { // Added snliqfrac by CJC2020
                            //there is enough water left over
                            if (liqWE != snWE) {
                                liqRoute = (liqWE - snliqfrac * iceWE); // Added snliqfrac by CJC2020
                                liqWE = liqWE - liqRoute;
                                snWE = liqWE + iceWE;
                            }
                                //there is no more pack
                            else {
                                liqRoute = snWE;
                                liqWE = 0.0;
                                iceWE = 0.0;
                                snWE = 0.0;
                            }
                        }
                        //set temperatures to 0 Celsius
                        snTempC = 0.0;
                        iceTempC = 0.0;
                        liqTempC = 0.0;
                    }//end -- melt
                }//end -- snow left after initial mass decrement
            }//end -- existing pack at beginning of time step

            //make sure that we still have snow
            if (snWE <= 5e-6) {
                liqRoute += snWE;
                snWE = 0.0;
                liqWE = 0.0;
                iceWE = 0.0;
                crustAge = 0.0;
                densityAge = 0.0;
                Utot = 0.0;
                Usn = 0.0;
                Uwat = 0.0;
            } else {
                crustAge += timeSteph;
                densityAge += timeSteph;
                snOnOff = 1.0;
            }

            //mass balance leaves >= 0 snow, then prepare for output in cm
            snWE = naughttocm * snWE;
            iceWE = naughttocm * iceWE;
            liqWE = naughttocm * liqWE;
            liqRoute = naughttocm * liqRoute;



            // Set ET variables equal to zero due to snowpack
            // ET variables are set to zero for canopy when snow in canopy, see tSnowIntercept
            cNode->setEvapSoil(0.0);
            cNode->setActEvap(0.0);

            setToNodeSnP(cNode);
            setToNode(cNode);

        }//end yes-snow


        // Estimate average Ep and cloudiness
        if (rainPtr->getoptStorm() && Io > 0.0) {
            potEvap = cNode->getPotEvap();
            SkyC += skyCover;
            cnt++;
        }

        //get the next node information for while loop
        cNode = nodeIter.NextP();
        count++;

    }//end while-nodes
    timeCount++;
    oldTimeStep = hourlyTimeStep;
    hourlyTimeStep++;

    // AJR2008, SKY2008Snow
    // Submit the basin average value
    if (rainPtr->getoptStorm()) {
        // SKY2008Snow -- Following bug corrected to account for SkyC division by cnt again in the next ComputeDailyEpCld call
        cnt ? SkyC = SkyC : SkyC = skyCover;

        // Get approximate EP from Priestley-Taylor method
        EP = ApproximateEP();

        // Submit values to climate simulator
        weatherSimul->ComputeDailyEpCld(EP, SkyC / cnt);

        // Assign the radiation variables to the 'tHydrometStoch'
        if (!count) {
            weatherSimul->setSunH(alphaD);
            weatherSimul->setSinH(sinAlpha);
            weatherSimul->setIo(Io);
            weatherSimul->setIdir(Ics);
            weatherSimul->setIdif(Ids);
            weatherSimul->OutputHydrometVars();
        }
    }
    return;
}

/****************************************************************************
**
**		      tSnowPack -- Interact w/ tCNode
**
**	This section of functions do the bulk on interacting w/ tCNode for
**	the tSnowPack class.
**
****************************************************************************/

//---------------------------------------------------------------------------
//
//			      tSnowPack::getFrNodeSnP()
//    
//	Get last time step's state variable information from tCNode object.
//
//---------------------------------------------------------------------------

void tSnowPack::getFrNodeSnP(tCNode *node) {

    liqWE = node->getLiqWE(); //cm
    iceWE = node->getIceWE(); //cm
    liqRoute = 0.0; //cm
    snWE = liqWE + iceWE; //cm

    //deal with the snow temperatures
    if (snWE > 1e-4) {

        snTempC = node->getSnTempC(); //Celsius
        iceTempC = snTempC; //Celsius
        liqTempC = 0.0; //default -- if this is 0.0, then the other should have been
        //		 set to 0.0 during the last time step.

    } else if (airTemp > 0.0) {
        //in case it is snowing out and the air temperature is greater than 0.0 (this is
        //actually a possibility). This also assumes that there is no pack.
        snTempC = 0.0;
        iceTempC = 0.0;
        liqTempC = 0.0;
    } else if (airTemp <= 0.0) {
        //if the air temperature is less than zero and there is no pack, assume
        //that the temperature of the phases are all the air temp.
        snTempC = airTemp;
        iceTempC = airTemp;
        liqTempC = airTemp;
    }

    crustAge = node->getCrustAge();
    densityAge = node->getDensityAge();
    ETAge = node->getEvapoTransAge();

    persMax = node->getPersTimeMax();
    persMaxtemp = node->getPersTime();

    peakSnWE = node->getPeakSWE();
    peakSnWEtemp = node->getPeakSWETemp();

    inittime = node->getInitPackTime();
    inittimeTemp = node->getInitPackTimeTemp();
    peaktime = node->getPeakPackTime();

    return;

}

//---------------------------------------------------------------------------
//
//			  tSnowPack::setToNodeSnP()
//
//	Sets the state variables, fluxes and outputs to the node.
//
//---------------------------------------------------------------------------

void tSnowPack::setToNodeSnP(tCNode *node) {

    //state variables
    node->setLiqWE(liqWE);
    node->setIceWE(iceWE);
    node->setSnTempC(snTempC);
    node->setCrustAge(crustAge);
    node->setDensityAge(densityAge);
    node->setEvapoTransAge(ETAge);
    node->setSnSub(snSub); // scaled by non-vegetated area
    node->setSnEvap(snEvap); // scaled by non-vegetated area

    //mass flux
    node->setLiqRouted(liqRoute);

    //energy fluxes, changes
    node->setDU(dUint);
    node->setUnode(Utot);
/*  node->setSnLHF(L*kilotonaught);
  node->setSnSHF(H*kilotonaught);
  node->setSnPHF(Prec*kilotonaught);
  node->setSnGHF(G*kilotonaught);
  node->setSnRLin(RLin*kilotonaught);
  node->setSnRLout(RLout*kilotonaught);
  node->setSnRSin(RSin*kilotonaught);*/

    //outputs
    node->setUerror(Uerr);

    //times and peaks
    if ((iceWE + liqWE) <= 1e-5) {
        persMaxtemp = 0.0;
        inittimeTemp = 0.0;
    }

    if (snWE > 1e-5) {
        persMaxtemp += 1.0;

        if (persMaxtemp > persMax) {
            persMax = persMaxtemp;

            if (persMaxtemp - 1 < 1)
                inittimeTemp = hourlyTimeStep;

        }
    }

    if (snWE >= peakSnWE) {
        peakSnWE = snWE;
        peaktime = hourlyTimeStep;
        inittime = inittimeTemp;
    }

    node->setPersTimeMax(persMax);
    node->setPersTime(persMaxtemp);
    node->setPeakSWE(peakSnWE);
    node->setPeakPackTime(double(peaktime));
    node->setInitPackTime(double(inittime));
    node->setInitPackTimeTemp(double(inittimeTemp));

    //cumulative outputs
    node->addLatHF(L);
    node->addSnSub(snSub); // cumulative snow sublimation CJC2020
    node->addSnEvap(snEvap); // cumulative snow evaporation CJC2020
    node->addMelt(liqRoute);
    node->addSHF(H * timeSteps);
    node->addPHF(Prec * timeSteps);
    node->addGHF(G * timeSteps);
    node->addRLin(RLin * timeSteps);
    node->addRLout(RLout * timeSteps);
    node->addRSin(RSin * timeSteps);
    node->addCumUerror(Uerr * timeSteps);
    node->addCumHrsSnow(snOnOff);

    //reset fluxes to zero
    L = H = Prec = G = RLin = RLout = RSin = dUint = 0.0;
    liqRoute = 0.0;

    return;
}


/****************************************************************************
**
**		  tSnowPack -- Internal Physical Routines
**
**	Deals with purely internal changes of the pack (not at surface)
**
****************************************************************************/

//---------------------------------------------------------------------------
//
//			  tSnowPack::densityFromAge()
//
//	This function should calculate density as a function of time. The
//	density should be in mks. (include a reference)
//
//						    (Tuteja et al 1996)
//
//---------------------------------------------------------------------------

double tSnowPack::densityFromAge() {

    double rhotemporary(400); //kg/m^3

    return rhotemporary;
}

/****************************************************************************
**
**		    tSnowPack -- EB Functions (basic calcs)
**
**	  This set of functions forms the basis of the energy balance 
**	  calculations later. This includes all of the heat fluxes (H, L, 
**	  P, G, Rn) and the ancillary information needed to calculate these.
**
****************************************************************************/

//---------------------------------------------------------------------------
//
//			  tSnowPack::latentHFCalc()
//
//	  Calculates the turbulent latent heat flux given a given aerodynamic
//	  resistivity.  There are two separate conditions. The first is 
//	  sublimation from dry snow and the second is evaporation of the 
//	  liquid phase. Other than that, it is straight forward (REF).
//
//	  RMK: DOES NOT TAKE INTO ACCOUNT ATMOSPHERIC INSTABILITY (DOES NOT
//	  INCORPORATE RICHARD'S NUMBER).
//
//							  (Dingman 2002)
//
//---------------------------------------------------------------------------

double tSnowPack::latentHFCalc(double Kaero) {

    double lhf;
    if (snTempC == 0.0)
        lhf = (latVapkJ * 0.622 * rhoAir * Kaero * (vPress - 6.111) / atmPress); //evaporation by THM 2012
    else
        lhf = (latSubkJ * 0.622 * rhoAir * Kaero * (vPress - 6.112 * exp((17.67 * snTempC) / (snTempC + 243.5))) /
               atmPress); //sublimation by THM 2012
    return lhf;
}

//----------------------------------------------------------------------------
//
//			    tSnowPack::sensibleHFCalc()
//
//	  Calculates the turbulent sensible heat from from the surface of the
//	  snow pack given a aerodynamic resistivity.
//
//	  RMK: DOES NOT ACCOUNT FOR ATMOSPHERIC INSTABILITY (DOES NOT
//	  INCORPORATE RICHARDS NUMBER).
//
//							    (Dingman 2002)
//
//----------------------------------------------------------------------------

double tSnowPack::sensibleHFCalc(double Kaero) {

    double shf;

    shf = (rhoAir * cpairkJ * Kaero * ((airTemp + 273.15) - (snTempC + 273.15)));
    return shf;
}

//-----------------------------------------------------------------------------
//
//			    tSnowPack::snowFracCalc()
//
//	  Partitions precipitation between solid and liquid phase based on 
//	  air temperature . This is linear interpolation b/t two
//	  temperatures. Ice can exist above freezing and liquid can exist below
//	  freezing, but there are certain minimun and maximum temperatures at 
//	  which the precipitation becomes all ice and all liquid, respectively.
//
//							(Wigmosta et al 1994)
//
//-----------------------------------------------------------------------------

double tSnowPack::snowFracCalc() {

    double snowfrac;

    double TMin(0), TMax(4.4); //indices (Wigmosta et al. 1994)—updated for CJC thesis (see table 11)

    if (airTemp <= TMin)
        snowfrac = 1; // all ice
    if (airTemp >= TMax)
        snowfrac = 0; // all liquid
    if ((airTemp >= TMin) && (airTemp <= TMax))
        snowfrac = (TMax - airTemp) / (TMax - TMin); //mixture


    return snowfrac;
}

//------------------------------------------------------------------------------
//
//			tSnowPack::precipitationHFCalc()
//			
//	  Finds the amount of energy advected into the pack from precipitation.
//	  It has two cases:
//
//		Ta > 0	  ==>	Tsol = 0,   Tliq = Ta
//		Ta <= 0	  ==>	Tsol = Ta,  Tliq = 0
//
//	  In other words, if the air temperature is greater than zero, the ice
//	  just froze and the liquid is at ambient temperature, but if the air
//	  temperature is less than 0, then the ice is at the ambient temperature
//	  and the liquid is freezing. This difference is displayed in the phase
//	  of the precipitation has the heat capacity and temperature incorportated
//	  into the calculation.
//
//------------------------------------------------------------------------------

double tSnowPack::precipitationHFCalc() {
    double phf;
    double frac;

    //  frac = snowFracCalc();
    snPrec = (snowFracCalc() * (rain + ctom * snUnload)) * mtoc; //convert from mm to cm
    liqPrec = ((1 - snowFracCalc()) * (rain + ctom * snUnload)) * mtoc; //convert from mm to cm

    if (airTemp > 0) {
        phf = (cmtonaught * snPrec * 0 * rholiqkg * cpicekJ +
               cmtonaught * liqPrec * (latFreezekJ + airTemp * rholiqkg * cpwaterkJ)) / 3600;

    } else if (airTemp <= 0) {
        phf = (cmtonaught * snPrec * airTemp * rholiqkg * cpicekJ +
               cmtonaught * liqPrec * latFreezekJ * rholiqkg) / 3600;

    }

    return phf;
}



//-----------------------------------------------------------------------------
//
//			    tSnowPack::latHeatVapCalc()
//
//	This function returns the latent heat of the vaporization in kJ/kg. It 
//	was constructed in order to ease incorporation of more precise functional
//	values of the latent heat of vaporization in a multilayered model.
//
//-----------------------------------------------------------------------------

double tSnowPack::latHeatVapCalc() {

    double lhvap(2470);

    return lhvap;
}

//-----------------------------------------------------------------------------
//
//			    tSnowPack::latHeatFreezeCalc()
//
//	This function returns the latent heat of the freezing in kJ/kg. It 
//	was constructed in order to ease incorporation of more precise functional
//	values of the latent heat of vaporization in a multilayered model.
//	
//-----------------------------------------------------------------------------

double tSnowPack::latHeatFreezeCalc() {

    double lhfreeze(334);

    return lhfreeze;
}

//------------------------------------------------------------------------------
//
//			      tSnowPack::latHeatSubCalc()
//
//
//	This function returns the latent heat of the sublimation in kJ/kg. It 
//	was constructed in order to ease incorporation of more precise functional
//	values of the latent heat of vaporization in a multilayered model.
//	
//------------------------------------------------------------------------------

double tSnowPack::latHeatSubCalc() {

    double lhsub(2470 + 334);

    return lhsub;
}

//------------------------------------------------------------------------------
//
//			    tSnowPack::heatCapVapCalc()
//
//	Returns the heat capacity of water vapor in kJ/kg. It was constructed
//	in order to allow easy incorporation of more precise values into a 
//	multilayered model.
//
//------------------------------------------------------------------------------

double tSnowPack::heatCapAirCalc() {
    double heatcapvap(1.01);

    return heatcapvap;
}

//-------------------------------------------------------------------------------
//
//			      tSnowPack::heatCapSolCalc()
//
//	Returns the heat capacity of solid water in kJ/kg. It was constructed
//	in order to allow easy incorporation of more precise values into a 
//	multilayered model.
//
//-------------------------------------------------------------------------------

double tSnowPack::heatCapSolCalc() {
    double heatcapsol(2.1);

    return heatcapsol;
}

//--------------------------------------------------------------------------------
//
//			      tSnowPack::heatCapLiqCalc()
//			      
//
//	Returns the heat capacity of liquid water in kJ/kg. It was constructed
//	in order to allow easy incorporation of more precise values into a 
//	multilayered model.
//
//--------------------------------------------------------------------------------

double tSnowPack::heatCapLiqCalc() {
    double heatcapliq(4.19);

    return heatcapliq;
}

//---------------------------------------------------------------------------------
//
//			    tSnowPack::vapPressSnowSurfCalc()
//
//	Returns the vapor pressure at the snow surface from the Clausisu-Clayperyon (sp?)
//	relationship. This assumes that the air is saturated at the snow surface.
//
//						    References: Dingman (2002)
//								Wigmosta et al (1994)
//								Anderson (1976)
//								    
//---------------------------------------------------------------------------------

double tSnowPack::vapPressSnowSurfCalc() {
    double vpresssnowsurf(0.0);

    // From Bras (1990)
    vpresssnowsurf = 6.112 * exp((17.67 * snTempC) / (snTempC + 243.5));

    return vpresssnowsurf;
}

//-----------------------------------------------------------------------------------
//
//				tSnowPack::agingAlbedo()
//
//	Returns the effective albedo of the surface for a given age of the snow
//	surface. This mainly is in response to crystal structure changes but also
//	heuristically deals with increased amounts of incorporated dust. It has 
//	two sets of curves. The first is during the accumulation period, assumed
//	to have predominantly dry snow; and the second is for during the melt period,
//	assumed to have wet snow.
//
//							      Wigmosta et al. (1994)
//							      CRREL Review of Albedo
//								  Modification.
//
//-----------------------------------------------------------------------------------

double tSnowPack::agingAlbedo() {
    double alb;

    if (liqWE < 1e-5)
        alb = 0.85 * pow(0.94, pow(crustAge / 24, 0.58)); // dry snow
    else
        alb = 0.85 * pow(0.82, pow(crustAge / 24, 0.46)); // wet snow

    return alb;
}

//-----------------------------------------------------------------------------------
//
//				tSnowPack::resFactCalc()
//
//	  Returns the effective aerodynamic resistance for use in the calculations
//	  of turbulent heat flux. This accounts for the change in effective 
//	  vegetation roughness height b/c of snow depth. The rest of the function is
//	  taken directly from tEvapoTrans.
//
//-----------------------------------------------------------------------------------

double tSnowPack::resFactCalc() {

    double rf, ra;
    double vonKarm = 0.41;
    double vegHeight, vegFrac, vegBare, windSpeedBare;
    double zm, zom, zov, d, rav, ras;

    if (coeffH == 0)
        vegHeight = 0.1;
    else
        vegHeight = coeffH; // vegHeight in meters

    //THM 2012 added for grassland
    if (coeffH < 1) {
        vegHeight = coeffH / 250;
        coeffV = 0.1;
    } else {
        vegHeight = coeffH;
    }
    if (vegHeight > snDepthm) {
        vegHeight = vegHeight - snDepthm;
    } else {
        vegHeight = 0.1; // aka height of snow
    }

    vegBare = 0.1; // height of bare soilc

    vegFrac = coeffV;

    if (windSpeed == 0.0 || fabs(windSpeed - 9999.99) < 1e-3) {
        windSpeedC = 0.01;    //Minimum wind speed (m/s)
    } else {
        windSpeedC = windSpeed;
    }

    // Compute below canopy windspeed at snow surface following equation Moreno et al. (2016) CJC 2020
    if (snDepthm < coeffH) {
        windSpeedC = windSpeedC * exp(-0.5 * coeffLAI * (1 - (snDepthm / coeffH)));
    }

    // Compute aerodynamic resistance for vegetation
    zm = 2.0 + vegHeight;
    zom = 0.13 * vegHeight;
    zov = 0.013 * vegHeight;
    d = 0.67 * vegHeight;
    rav = log((zm - d) / zom) * log((zm - d) / zov) / (windSpeedC * pow(vonKarm, 2));

    // Compute aerodynamic resistance for bare soil
    zm = 2.0 + vegBare;
    zom = 0.13 * vegBare;
    zov = 0.013 * vegBare;
    d = 0.67 * vegBare;

    ras = log((zm - d) / zom) * log((zm - d) / zov) / (windSpeedC * pow(vonKarm, 2));

    ra = (1 - vegFrac) * ras + vegFrac * rav;
    rf = 1 / ra; // Otherwise known as kaero

    return rf;
}

//-------------------------------------------------------------------------
//
// tSnowPack::inShortWaveSn() Function
//
//  This is modified, including comments, from tEvapoTrans::inShorWave().
//  Where we have modified things, it is noted.
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

double tSnowPack::inShortWaveSn(tCNode *cNode) {
    double Is, N, Iv, Isw, Ir;
    double v, t, cosi, scover;
    double RadGlobClr;

    //Remaining variables in DirectDiffuse from v3 -- AJR2008, SKY2008Snow
    //  So, entire function changed to match inShortWave in tEvapoTrans
/*  double h0, m, pp0, Dh0ref, h0ref, drm, Tlinke;
  double TnTLK, Fdh0, A1p, A1, A2, A3;
  double pi = 4*atan(1.0);*/

    Ic = Is = Id = Ir = Ids = Ics = Isw = Iv = 0.0;

    // Elevation, Slope and Aspect have been set before

    SunHour = 0.0; //Rinehart 2007 -- initialize whether we see sun or not to NO

    if (alphaD > 0.0) {

        elevation = cNode->getZ(); //SMM 10142008
        DirectDiffuse(Tlinke, elevation);  // SKY2008Snow, AJR2007

        // Cloud cover information
        if (fabs(skyCover - 9999.99) < 1.0E-3) {

            skyCover = compSkyCover();//ADDED BY RINEHART 2007 @ NMT
            // computes sky cover from relative
            // humidity and rain.
            scover = skyCover;

            //if (rain > 0.0) scover = 10.0;
            //else            scover = 1.0;
        } else
            scover = skyCover;
        N = scover / 10.0;

        // If observations (for a horizontal surface) exist -
        // use them, at least in an approximate manner
        if (tsOption > 1 && !rainPtr->getoptStorm()) {
            RadGlobClr = (RadGlbObs / (1.0 - 0.65 * pow(N, 2.0)));
            Ic = Ic / (Ic * sinAlpha + Id) * RadGlobClr;
            Id = RadGlobClr - Ic * sinAlpha;
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

            cosi = 1.0 * (cos(slope) * sinAlpha + sin(slope) * cos(asin(sinAlpha)) * cos(sunaz - aspect));

            if (cosi >= 0.0) {
                Ics = Ic * cosi;
                SunHour = 1.0; //YES SUN
            } else {
                Ics = 0.0;
                SunHour = 0.0; //NO SUN
            }
        } else {
            Ics = 1.0 * Ic;
            SunHour = 1.0;
        }

        if ((shelterOption == 2) || (shelterOption == 1)) {//CHANGED IN 2008

            Ics *= aboveHorizon(ID); //check to see if we can see the sun (aboveHorizon() in tEvapoTrans)
            SunHour *= aboveHorizon(ID);
        }

        //2) Horizon factor for diffuse radiation?
        //Rinehart 2007 @ New Mexico Tech
        //
        //	See comment above about sheltering options.

        if ((shelterOption > 0) && (shelterOption < 3)) {
            v = shelterFactorGlobal; //incorporate remote sheltering
        } else if (shelterOption == 0 || shelterOption == 3) { //CHANGED 2008
            v = 0.5 * (1 + cos(slope)); //local sheltering
        } else {
            v = 1.0; // no sheltering
        }

        Ids = Id * v;

        // 3) Account for cloud cover
        Is = (1.0 - 0.65 * pow(N, 2)) * (Ics + Ids);

        //Reflected from surrounded sites radiation
        //
        //Modified by Rinehart 2007 @ New Mexico Tech
        //
        if (hillAlbedoOption == 0) {
            hillalbedo = albedo;
        } else if (hillAlbedoOption == 1) {
            hillalbedo = coeffAl;
        } else if (hillAlbedoOption == 2) {

            if (snCanWE == 0)
                hillalbedo = coeffV * coeffAl + (1 - coeffV) * albedo;
            else
                hillalbedo = albedo;
        }

        if (shelterOption == 0) {
            //local
            Ir = hillalbedo * Is * (1 - cos(slope)) * 0.5;
            Is += Ir;
        } else if ((shelterOption > 1) && (shelterOption < 4)) { //CHANGED IN 2008
            //remote
            Ir = hillalbedo * Is * (0.5 * (1 + cos(slope)) - shelterFactorGlobal); //CHANGED IN 2008
            landRefGlobal = 0.5 * (1 + cos(slope)) - shelterFactorGlobal;
            Is += Ir;

        } else { //CHANGED IN 2008
            Ir = 0.0;
            Is = Is;
        }

        // Account for vegetation
        if ((evapotransOption == 1) || (snowOption)) {
            Iv = Is * coeffKt * coeffV + Is * (1.0 - coeffV);
            //Iv = Is*exp((coeffKt-1)*coeffLAI)*coeffV + Is*(1.0-coeffV); // Changed to use Beer-Lambert following Moreno et al. (2016) CJC 2020
        } else
            Iv = Is;

        // Account for albedo
        //Modified by Rinehart 2007 @ New Mexico Tech
        //	This is actually the main difference b/t the tEvapoTrans::inShortWave
        //	and this function. We no longer see the land surface and have
        //	calculated albedo as a function of surface age earlier in the
        //	algorithm.

        Isw = Iv * (1.0 - albedo);

    } //end -- alphaD > 0
    else {
        Ic = Is = N = Iv = Isw = Id = Ids = Ics = Ir = 0.0;
    } // end -- alphaD <= 0


    // Assign the radiation variables to the 'tHydrometStoch' for ID = 0
    if (rainPtr->getoptStorm() && (ID == 0)) {
        weatherSimul->setSunH(alphaD);
        weatherSimul->setIdir(Ics);
        weatherSimul->setIdir_vis(0.5 * Ics);
        weatherSimul->setIdir_nir(0.5 * Ics);
        weatherSimul->setIdif(Ids + Ir);//AJR2008, SKY2008Snow
        weatherSimul->setIdif_vis(0.5 * (Ids + Ir));//AJR2008, SKY2008Snow
        weatherSimul->setIdif_nir(0.5 * (Ids + Ir));//AJR2008, SKY2008Snow
        weatherSimul->OutputHydrometVars();
    }

    // Set shortwave variables to the node (partition is approximate)
    if (tsOption > 1 && !rainPtr->getoptStorm()) {
        cNode->setShortRadIn(RadGlbObs); //or set(Is), they must be equal
    } else {
        cNode->setShortRadIn(Isw);
    }
    cNode->setShortRadIn_dir(Ics * (1.0 - 0.65 * pow(N, 2.0)));
    cNode->setShortRadIn_dif((Ids + Ir) * (1.0 - 0.65 * pow(N, 2.0)));//AJR2008, SKY2008Snow

    return Isw;
}

//----------------------------------------------------------------------------
//
//			      tSnowPack::emmisSn()
//
//	Calculates the emmissivity of snow. Similar to the thermal properties
//	above, this is here for incorporation into a multilayer model later.
//
//----------------------------------------------------------------------------

double tSnowPack::emmisSn() {

    double emiss = 0.9;
    return emiss;
}

//-----------------------------------------------------------------------------
//
//			      tSnowPack::inLongWaveSn() 
//
// Function based on gray-body (Bras(1990) p44)
//
//        Rlin(t) = K(t)*Ea(t)*S*T(t)^4   (J/(m^2*s)) 
//
//        S      Stefan-Boltzmann Constant = 5.67e-8  J/(s*m^2*K^4)
//        K(t)   Cloud cover coefficient K(t) = (1+0.17*N^2) []
//        N(t)   Sky Cover/10
//        Ea(t)  Atmospheric Thermal Emissivity  Idso(1981)
//               Ea(t) = 0.74 + 0.0049*e(t)  []
//        T(t)   Air Temperature (K)  
//
// Assigns corrected SkyCover 
//
// This function was constructed when we were focusing on the controls of
// snow in canopy. Here, we could correct for the effect of 0 C temperatures
// in the canopy compared more variable ambient temperature. B/c of problems
// with the admittedly simple implementation, I have choosen not to 
// incorporate those physics.
//
//    Rinehart 2007 @ New Mexico Tech
//    
//---------------------------------------------------------------------------

double tSnowPack::inLongWaveSn() {
    double sigma, kCloud, airTempK, Ea, Rlin, scover, N;
    double v0;

    if (shelterOption < 3)
        v0 = shelterFactorGlobal;
    else
        v0 = 1;

    sigma = 5.67e-8;
    if (fabs(skyCover - 9999.99) < 1e-3) {
        skyCover = compSkyCover();
        scover = skyCover;
    } else
        scover = skyCover;

    skyCoverC = scover;

    N = scover / 10.0;
    kCloud = 1 + 0.17 * pow(N, 2);
    airTempK = airTemp + 273.15;
    Ea = 0.74 + 0.0049 * vaporPress();
    if (canWE <= 1e-3) // where we would account for no-snow in canopy
        Rlin = v0 * kCloud * Ea * sigma * pow(airTempK, 4);
    else // account for snow in canopy
        Rlin = v0 *
               ((1 - coeffV) * kCloud * Ea * sigma * pow(airTempK, 4) + coeffV * emmisSn() * sigma * pow(273.15, 4));

    return Rlin;
}

/************************************************************************************
**
**			tSnowPack -- Energy Balance Functions 
**
**	This set of functions forms the computational heart of the code. It calculates
**	the change in energy in the snow pack over a single time step.
**
**	This set of functions will need to be expanded for a multilayer model.
**	
**	Units are all kW/m^2
**
************************************************************************************/

//-----------------------------------------------------------------------------------
//
//				  tSnowPack::snowEB()
//
//	Calculates the change in energy over a single time step for the single-layer
//	model given all of the fluxes.
//
//-----------------------------------------------------------------------------------

void tSnowPack::snowEB(int nodeID, tCNode *node) {

    double sigma(5.67e-8);
    double v1;

    if (shelterOption > 0 && shelterOption < 3)
        v1 = shelterFactorGlobal;
    else
        v1 = 1;

    //convert temperature
    snTempK = CtoK(snTempC);

    //set up resistance
    resFact = resFactCalc();

    //turbulent heat fluxes
    H = sensibleHFCalc(resFact);
    L = latentHFCalc(resFact);

    //precipitation heat flux
    Prec = phfOnOff * precipitationHFCalc();

    //atmospheric heat flux
    RSin = naughttokilo * inShortWaveSn(node); // AJR2008, SKY2008Snow
    RLin = naughttokilo * inLongWave(node); // AJR2008, SKY2008Snow
    RLout = -naughttokilo * v1 * emmisSn() * sigma * pow(snTempK, 4.0);


    //set up for output
    inShortR = kilotonaught * RSin;
    inLongR = kilotonaught * RLin;
    outLongR = kilotonaught * RLout;

    // SKY2008Snow, AJR 2008
    //	Set the non-snow fluxes to zero
    node->setHFlux(0.0);
    node->setLFlux(0.0);
    node->setLongRadIn(0.0);
    node->setLongRadOut(outLongR);
    node->setShortAbsbVeg(0.0);
    node->setShortAbsbSoi(0.0);

    // Set the snow fluxes
    node->setSnLHF(L * kilotonaught);
    node->setSnSHF(H * kilotonaught);
    node->setSnPHF(Prec * kilotonaught);
    node->setSnGHF(G * kilotonaught);
    node->setSnRLin(RLin * kilotonaught);
    node->setSnRLout(RLout * kilotonaught);
    node->setSnRSin(RSin * kilotonaught);


    Rn = RSin + RLin + RLout;
    G = 0;

    //calculate total dU over given timestep
    dUint = (H + Rn + L + Prec + G) * timeSteps; // Changed *3600 to *timeSteps CJC2020

    //find new energy state of snow
    Utot += dUint;
    return;
}

/*****************************************************************************
**
**			  tSnowPack I/O Functions
**
**	This set of functions deal with the I/O of tSnowPack that is not already
**	derived from tEvapoTrans.
**
*****************************************************************************/

int tSnowPack::getSnowOpt() {
    return snowOption; //found in tEvapoTrans construction
}


/*****************************************************************************
**
**			  tSnowPack Conversion Functions
**
**	This set of functions deal with the I/O of tSnowPack that is not already
**	derived from tEvapoTrans.
**
*****************************************************************************/

//---------------------------------------------------------------------------
//
//				tSnowPack::CtoK()
//	
//	Convert Celsius to Kelvin
//	
//---------------------------------------------------------------------------

double tSnowPack::CtoK(double temperature) {

    return (temperature + 273.15);

}

//---------------------------------------------------------------------------
//
//				tSnowPack::KtoC()
//
//	Convert Kelvin to Celsius
//---------------------------------------------------------------------------

double tSnowPack::KtoC(double temperature) {

    return (temperature - 273.15);

}

/***************************************************************************
**
** tSnowPack::writeRestart() Function
**
** Called from tSimulator during simulation loop
**
***************************************************************************/
void tSnowPack::writeRestart(fstream &rStr) const {
    BinaryWrite(rStr, hillAlbedoOption);
    BinaryWrite(rStr, densityAge);
    BinaryWrite(rStr, rainTemp);
    BinaryWrite(rStr, ETAge);

    BinaryWrite(rStr, timeSteph);
    BinaryWrite(rStr, timeSteps);
    BinaryWrite(rStr, timeStepm);
    BinaryWrite(rStr, minutelyTimeStep);

    BinaryWrite(rStr, liqWE);
    BinaryWrite(rStr, iceWE);
    BinaryWrite(rStr, snWE);
    BinaryWrite(rStr, canWE);
    BinaryWrite(rStr, liqRoute);
    BinaryWrite(rStr, liqWEm);
    BinaryWrite(rStr, iceWEm);
    BinaryWrite(rStr, snWEm);
    BinaryWrite(rStr, liqRoutem);
    BinaryWrite(rStr, Utot);
    BinaryWrite(rStr, Usn);
    BinaryWrite(rStr, Uwat);
    BinaryWrite(rStr, Utotold);
    BinaryWrite(rStr, liqWatCont);
    BinaryWrite(rStr, liqTempC);
    BinaryWrite(rStr, iceTempC);
    BinaryWrite(rStr, snTempC);
    BinaryWrite(rStr, liqTempK);
    BinaryWrite(rStr, iceTempK);
    BinaryWrite(rStr, snTempK);
    BinaryWrite(rStr, crustAge);

    BinaryWrite(rStr, H);
    BinaryWrite(rStr, L);
    BinaryWrite(rStr, G);
    BinaryWrite(rStr, Prec);
    BinaryWrite(rStr, Rn);
    BinaryWrite(rStr, dUint);
    BinaryWrite(rStr, RLin);
    BinaryWrite(rStr, RLout);
    BinaryWrite(rStr, RSin);
    BinaryWrite(rStr, Uerr);
    BinaryWrite(rStr, snPrec);
    BinaryWrite(rStr, liqPrec);
    BinaryWrite(rStr, snPrecm);
    BinaryWrite(rStr, liqPrecm);
    BinaryWrite(rStr, snPrecmm);
    BinaryWrite(rStr, liqPrecmm);
    BinaryWrite(rStr, snUnload);
    BinaryWrite(rStr, snCanWE);
    BinaryWrite(rStr, vapPressSmb);
    BinaryWrite(rStr, vapPresskSPa);

    BinaryWrite(rStr, rholiqcgs);
    BinaryWrite(rStr, rhoicecgs);
    BinaryWrite(rStr, rhosncgs);
    BinaryWrite(rStr, rholiqkg);
    BinaryWrite(rStr, rhoicekg);
    BinaryWrite(rStr, rhosnkg);
    BinaryWrite(rStr, rhoAir);
    BinaryWrite(rStr, phfOnOff);

    BinaryWrite(rStr, cpsnowkJ);
    BinaryWrite(rStr, cpicekJ);
    BinaryWrite(rStr, cpwaterkJ);
    BinaryWrite(rStr, cpairkJ);
    BinaryWrite(rStr, latFreezekJ);
    BinaryWrite(rStr, latVapkJ);
    BinaryWrite(rStr, latSubkJ);

    BinaryWrite(rStr, resFact);
    BinaryWrite(rStr, albedo);
    BinaryWrite(rStr, hillalbedo);
    BinaryWrite(rStr, compactParam);
    BinaryWrite(rStr, rhoSnFreshkg);
    BinaryWrite(rStr, minSnTemp);

    BinaryWrite(rStr, snDepth);
    BinaryWrite(rStr, snDepthm);
    BinaryWrite(rStr, snOnOff);
    BinaryWrite(rStr, peakSnWE);
    BinaryWrite(rStr, peakSnWEtemp);
    BinaryWrite(rStr, persMax);
    BinaryWrite(rStr, persMaxtemp);
    BinaryWrite(rStr, inittime);
    BinaryWrite(rStr, inittimeTemp);
    BinaryWrite(rStr, peaktime);

    BinaryWrite(rStr, naughttokilo);
    BinaryWrite(rStr, kilotonaught);
    BinaryWrite(rStr, cgsRHOtomks);
    BinaryWrite(rStr, mksRHOtocgs);
    BinaryWrite(rStr, naughttocm);
    BinaryWrite(rStr, cmtonaught);
    BinaryWrite(rStr, ctom);
    BinaryWrite(rStr, mtoc);

    tEvapoTrans::writeRestart(rStr);
}

/***************************************************************************
**
** tSnowPack::readRestart() Function
**
***************************************************************************/
void tSnowPack::readRestart(fstream &rStr) {
    BinaryRead(rStr, hillAlbedoOption);
    BinaryRead(rStr, densityAge);
    BinaryRead(rStr, rainTemp);
    BinaryRead(rStr, ETAge);

    BinaryRead(rStr, timeSteph);
    BinaryRead(rStr, timeSteps);
    BinaryRead(rStr, timeStepm);
    BinaryRead(rStr, minutelyTimeStep);

    BinaryRead(rStr, liqWE);
    BinaryRead(rStr, iceWE);
    BinaryRead(rStr, snWE);
    BinaryRead(rStr, canWE);
    BinaryRead(rStr, liqRoute);
    BinaryRead(rStr, liqWEm);
    BinaryRead(rStr, iceWEm);
    BinaryRead(rStr, snWEm);
    BinaryRead(rStr, liqRoutem);
    BinaryRead(rStr, Utot);
    BinaryRead(rStr, Usn);
    BinaryRead(rStr, Uwat);
    BinaryRead(rStr, Utotold);
    BinaryRead(rStr, liqWatCont);
    BinaryRead(rStr, liqTempC);
    BinaryRead(rStr, iceTempC);
    BinaryRead(rStr, snTempC);
    BinaryRead(rStr, liqTempK);
    BinaryRead(rStr, iceTempK);
    BinaryRead(rStr, snTempK);
    BinaryRead(rStr, crustAge);

    BinaryRead(rStr, H);
    BinaryRead(rStr, L);
    BinaryRead(rStr, G);
    BinaryRead(rStr, Prec);
    BinaryRead(rStr, Rn);
    BinaryRead(rStr, dUint);
    BinaryRead(rStr, RLin);
    BinaryRead(rStr, RLout);
    BinaryRead(rStr, RSin);
    BinaryRead(rStr, Uerr);
    BinaryRead(rStr, snPrec);
    BinaryRead(rStr, liqPrec);
    BinaryRead(rStr, snPrecm);
    BinaryRead(rStr, liqPrecm);
    BinaryRead(rStr, snPrecmm);
    BinaryRead(rStr, liqPrecmm);
    BinaryRead(rStr, snUnload);
    BinaryRead(rStr, snCanWE);
    BinaryRead(rStr, vapPressSmb);
    BinaryRead(rStr, vapPresskSPa);

    BinaryRead(rStr, rholiqcgs);
    BinaryRead(rStr, rhoicecgs);
    BinaryRead(rStr, rhosncgs);
    BinaryRead(rStr, rholiqkg);
    BinaryRead(rStr, rhoicekg);
    BinaryRead(rStr, rhosnkg);
    BinaryRead(rStr, rhoAir);
    BinaryRead(rStr, phfOnOff);

    BinaryRead(rStr, cpsnowkJ);
    BinaryRead(rStr, cpicekJ);
    BinaryRead(rStr, cpwaterkJ);
    BinaryRead(rStr, cpairkJ);
    BinaryRead(rStr, latFreezekJ);
    BinaryRead(rStr, latVapkJ);
    BinaryRead(rStr, latSubkJ);

    BinaryRead(rStr, resFact);
    BinaryRead(rStr, albedo);
    BinaryRead(rStr, hillalbedo);
    BinaryRead(rStr, compactParam);
    BinaryRead(rStr, rhoSnFreshkg);
    BinaryRead(rStr, minSnTemp);

    BinaryRead(rStr, snDepth);
    BinaryRead(rStr, snDepthm);
    BinaryRead(rStr, snOnOff);
    BinaryRead(rStr, peakSnWE);
    BinaryRead(rStr, peakSnWEtemp);
    BinaryRead(rStr, persMax);
    BinaryRead(rStr, persMaxtemp);
    BinaryRead(rStr, inittime);
    BinaryRead(rStr, inittimeTemp);
    BinaryRead(rStr, peaktime);

    BinaryRead(rStr, naughttokilo);
    BinaryRead(rStr, kilotonaught);
    BinaryRead(rStr, cgsRHOtomks);
    BinaryRead(rStr, mksRHOtocgs);
    BinaryRead(rStr, naughttocm);
    BinaryRead(rStr, cmtonaught);
    BinaryRead(rStr, ctom);
    BinaryRead(rStr, mtoc);

    tEvapoTrans::readRestart(rStr);
}

void tSnowPack::checkShelter(tCNode *cNode) {
    if ((shelterOption > 0) && (shelterOption < 4)) {//CHANGED 2008
        int tempIndex = 0;
        for (tempIndex = 0; tempIndex < 16; tempIndex++) {
            switch (tempIndex) {
                case 0:
                    ha2700 = cNode->getHorAngle2700();
                    break;
                case 1:
                    ha2925 = cNode->getHorAngle2925();
                    break;
                case 2:
                    ha3150 = cNode->getHorAngle3150();
                    break;
                case 3:
                    ha3375 = cNode->getHorAngle3375();
                    break;
                case 4:
                    ha0000 = cNode->getHorAngle0000();
                    break;
                case 5:
                    ha0225 = cNode->getHorAngle0225();
                    break;
                case 6:
                    ha0450 = cNode->getHorAngle0450();
                    break;
                case 7:
                    ha0675 = cNode->getHorAngle0675();
                    break;
                case 8:
                    ha0900 = cNode->getHorAngle0900();
                    break;
                case 9:
                    ha1125 = cNode->getHorAngle1125();
                    break;
                case 10:
                    ha1350 = cNode->getHorAngle1350();
                    break;
                case 11:
                    ha1575 = cNode->getHorAngle1575();
                    break;
                case 12:
                    ha1800 = cNode->getHorAngle1800();
                    break;
                case 13:
                    ha2025 = cNode->getHorAngle2025();
                    break;
                case 14:
                    ha2250 = cNode->getHorAngle2250();
                    break;
                case 15:
                    ha2475 = cNode->getHorAngle2475();
                    break;
                default:
                    cout << "\nCheck tempInd -- did not exist or assign" << endl;
            }//end-switch
        }    //end-for

        shelterFactorGlobal = cNode->getSheltFact(); //computed in tShelter

    } else if (shelterOption == 0) { // local sheltering for factor only
        shelterFactorGlobal = 0.5 * (1 + cos(slope)); //computed here for output purposes
        cNode->setSheltFact(shelterFactorGlobal); // SKYnGM2008LU
    } else {
        shelterFactorGlobal = 1; //no sheltering
    }
    return;
}


/******************************************************************************
**
**		    tSnowPack -- END OF TSNOWPACK.CPP
**
******************************************************************************/
