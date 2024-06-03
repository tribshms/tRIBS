# This test confirms if the inputs for a provided model match the outputs.
# Currently this only setup for a single element case, but will need to be updated to address more complex multi-element
# models.

# pytest --binary /path/to/your/binary/executable --input /path/to/input/file

import numpy as np

def test_match_precip_input_to_output(setup_data):
    pixel, _, input_precip, m, _ = setup_data

    p_in = m.read_precip_station(input_precip)
    # Drop the last hour from precipitation data
    # this would automatically fail the test if not in place because last hour is not recorded.
    # This should be considered a bug.
    p_in.drop(p_in.index[-1], inplace=True)

    # Compare precipitation inputs and outputs
    pdiff = p_in.R.values - pixel.Rain_mm_h.values

    print("\nTesting Precip Match ")
    assert not np.any(pdiff), "Difference found in precipitation forcing"

def test_match_PA_input_to_output(setup_data):
    pixel, input_met, _, m, _ = setup_data
    # Read meteorological input
    met_in = m.read_met_station(input_met)
    # Drop the last hour from meteorological data
    met_in.drop(met_in.index[-1], inplace=True)

    tol = 0.01
    invar = 'PA'
    pvar = 'Press_Pa'

    var_dif = np.mean(met_in[invar].values - pixel[pvar].values)

    print("\nTesting Pressure Match ")
    assert var_dif <= tol, f"Difference found in variables {invar} and {pvar} at time {pixel.Time_hr[(met_in[invar] != pixel[pvar]).idxmax()]}"
def test_match_RH_input_to_output(setup_data):
    pixel, input_met, _, m, _ = setup_data
    # Read meteorological input
    met_in = m.read_met_station(input_met)
    # Drop the last hour from meteorological data
    met_in.drop(met_in.index[-1], inplace=True)

    tol = 0.01
    invar = 'RH'
    pvar = 'RelHum_[]'

    var_dif = np.mean(met_in[invar].values - pixel[pvar].values)

    print("\nTesting Relative Humidity Match ")
    assert var_dif <= tol, f"Difference found in variables {invar} and {pvar} at time {pixel.Time_hr[(met_in[invar] != pixel[pvar]).idxmax()]}"

def test_match_XC_input_to_output(setup_data):
    pixel, input_met, _, m, _ = setup_data
    # Read meteorological input
    met_in = m.read_met_station(input_met)
    # Drop the last hour from meteorological data
    met_in.drop(met_in.index[-1], inplace=True)

    tol = 0.01
    tribs_nodata = 9999.99
    invar = 'XC'
    pvar = 'SkyCov_[]'

    var_dif = np.mean(met_in[invar].values - pixel[pvar].values)

    print("\nTesting Sky Cover Match")

    # SkyCov_[] should all be less than 10 if estimated by tRIBS, which happens when supplied with nodata
    if all(met_in['XC'] == tribs_nodata):
        assert all(pixel[pvar]<=10.0)
    else:
        assert var_dif <= tol, f"Difference found in variables {invar} and {pvar} at time {pixel.Time_hr[(met_in[invar] != pixel[pvar]).idxmax()]}"
def test_match_US_input_to_output(setup_data):
    pixel, input_met, _, m, _ = setup_data
    # Read meteorological input
    met_in = m.read_met_station(input_met)
    # Drop the last hour from meteorological data
    met_in.drop(met_in.index[-1], inplace=True)

    min_windspeed = 0.1

    print("\nTesting Wind Speed Match ")
    assert all(pixel['Wind_m_s'][met_in['US'] != pixel['Wind_m_s']] == min_windspeed)

def test_match_TA_input_to_output(setup_data):
    pixel, input_met, _, m, _ = setup_data
    # Read meteorological input
    met_in = m.read_met_station(input_met)
    # Drop the last hour from meteorological data
    met_in.drop(met_in.index[-1], inplace=True)

    tol = 0.01
    invar = 'TA'
    pvar = 'AirT_oC'

    var_dif = np.mean(met_in[invar].values - pixel[pvar].values)

    print("\nTesting Air Temp. Match ")
    assert var_dif <= tol, f"Difference found in variables {invar} and {pvar} at time {pixel.Time_hr[(met_in[invar] != pixel[pvar]).idxmax()]}"

def test_match_IS_input_to_output(setup_data):
    pixel, input_met, _, m, _ = setup_data
    # Read meteorological input
    met_in = m.read_met_station(input_met)
    # Drop the last hour from meteorological data
    met_in.drop(met_in.index[-1], inplace=True)

    tol = 0.01
    invar = 'IS'
    pvar = 'ShrtRadIn_W_m2'

    var_dif = np.mean(met_in[invar].values - pixel[pvar].values)

    print("\nTesting Insolation Match ")
    assert var_dif <= tol, f"Difference found in variables {invar} and {pvar} at time {pixel.Time_hr[(met_in[invar] != pixel[pvar]).idxmax()]}"

def test_elem_water_balance(setup_data):

    # threshold for acceptable difference (or loss) of water between input precipitation and losses plus change in
    # storage. We want to achieve single digit differences or less, so this maybe further adjusted
    wb_threshold = 10 # mm/year

    print(f"\nTesting Conservation of Mass & Water Balance.\nThreshold = {wb_threshold} mm/yr")

    # read in results and assign relevant parameters
    pixel, _, _, _, r = setup_data
    porosity = r.int_spatial_vars.loc[r.int_spatial_vars.ID == 0, 'Porosity'].values[0]
    element_area = r.int_spatial_vars.loc[r.int_spatial_vars.ID == 0, 'VAr'].values[0]
    years = float(r.options['runtime']['value'])/(365.25*24)

    # change in storage
    unsat_mm = pixel['Mu_mm'].iloc[len(pixel)-1]-pixel['Mu_mm'].iloc[0] #mm
    sat_mm = (pixel['Nwt_mm'].iloc[0]-pixel['Nwt_mm'].iloc[len(pixel)-1]) * porosity #mm, ordering is reversed in cancelation of bedrock term
    can_swe_mm = 10*pixel['IntSWEq_cm'].iloc[len(pixel)-1]-pixel['IntSWEq_cm'].iloc[0] #mm
    can_mm = pixel['CanStorage_mm'].iloc[len(pixel)-1]-pixel['CanStorage_mm'].iloc[0] #mm
    swe_mm = 10*pixel['SnWE_cm'].iloc[len(pixel)-1]-pixel['SnWE_cm'].iloc[0] #mm

    # cumulative fluxes
    # assumed unit conversion for rates where each entry is scaled by a 1 hr time step then summed to a get per a area depth
    p_mm = pixel['Rain_mm_h'].sum()
    et_mm = pixel['EvpTtrs_mm_h'].sum() - (pixel['SnSub_cm'].sum() * 10 + pixel['SnEvap_cm'].sum() * 10 + pixel['IntSub_cm'].sum() * 10)
    qsurf_mm = pixel['Srf_Hour_mm'].sum()
    qunsat_mm =  np.sum(pixel['QpIn_mm_h'].values - pixel['QpOut_mm_h'].values)
    qsat_mm = pixel['GWflx_m3_h'].sum() / element_area * 1000 # conversion from m^3 to mm


    # water balance, ideally test should be ~= 0
    del_s = unsat_mm + sat_mm + can_mm + can_swe_mm
    net_loss = et_mm + qsurf_mm + qunsat_mm + qsat_mm
    net_flux = p_mm - net_loss

    # test avg. annual value
    test = (net_flux - del_s)/years
    print(f"Test metric (P - Loss - delS) = {test} mm/yr")

    assert test < wb_threshold

def test_model_skills(setup_data):
    pass