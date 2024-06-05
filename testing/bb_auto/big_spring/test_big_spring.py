def test_parallel_serial_match(setup_data):
    print("\nTesting congruence between parallel and serial simulations")
    parallel_result, serial_result, _ = setup_data

    diff_mask = parallel_result != serial_result
    different_columns = diff_mask.any(axis=0)
    result = different_columns.index[different_columns].tolist()

    if result:
        flag = True  # List is not empty
    else:
        flag = False  # List is empty

    assert not flag, f"The following variables do not match between the parallel and serial simulation: {result}"


def test_parallel_waterbalance(setup_data):
    mrf, _, r_par = setup_data

    print("\nTesting Big Spring (parallel) Water Balance")

    wb_threshold = 10  # mm/year

    porosity = r_par.int_spatial_vars.loc[r_par.int_spatial_vars.ID == 0, 'Porosity'].values[0]
    drainage_area = r_par.int_spatial_vars['VAr'].sum()
    years = float(r_par.options['runtime']['value']) / (365.25 * 24)

    # calculate changes in stores--note the assumption is that boundaries are closed except for streams
    # so sat and unsat fluxes = 0, also currently there is no average canopy storage variables so this portion of the
    # storage is ignored.

    unsat_mm = (mrf['MSMU'].iloc[len(mrf)-1] * mrf['MDGW'].iloc[len(mrf)-1] - mrf['MSMU'].iloc[0] * mrf['MDGW'].iloc[0]) * porosity
    sat_mm = (mrf['MDGW'].iloc[0] - mrf['MDGW'].iloc[0]) * porosity
    can_swe_mm = 10 * (mrf['AvInSn'].iloc[len(mrf)-1] - mrf['AvInSn'].iloc[0])
    swe_mm = 10 * (mrf['AvSWE'].iloc[len(mrf)-1] - mrf['AvSWE'].iloc[0])
    p_mm = mrf['MAP'].sum()
    et_mm = mrf['MET'].sum() - 10*(mrf['AvSnSub'].sum() + mrf['AvSnEvap'].sum() + mrf['AvInSu'].sum())
    q_mm = mrf['Srf'].sum() * 3600 * 1000 / drainage_area,

    del_s = unsat_mm + sat_mm + can_swe_mm + swe_mm
    loss = et_mm + q_mm
    test = (p_mm - del_s - loss)/years # mm/yr

    # test avg. annual value
    print(f"Test metric (P - Loss - delS) = {test} mm/yr")

    assert test < wb_threshold

