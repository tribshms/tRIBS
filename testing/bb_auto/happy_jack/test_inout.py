# This test confirms if the inputs for a provided model match the outputs.
# Currently this only setup for a single element case, but will need to be updated to address more complex multi-element
# models.

#pytest --binary /path/to/your/binary/executable --input /path/to/input/file

import sys
import os
import numpy as np
import pytest

extra_path = f"{os.environ['HOME']}/Documents/Repos/Forked/pytRIBS"
if extra_path not in sys.path:
    sys.path.append(extra_path)

from pytRIBS.classes import Model as model
from pytRIBS.classes import Results as results

@pytest.fixture
def setup_data():
    binary_path = '/Users/wr/Documents/Repos/Forked/tRIBS/cmake-build-serial/tRIBS'
    input_file = '/Users/wr/Desktop/tScenarioTesting/scenarios/one/scenario_one.in'
    input_precip = '/Users/wr/Desktop/tScenarioTesting/data/HJ_bench/data/HJ_PRECIP_2002-2018.mdf'
    input_met = '/Users/wr/Desktop/tScenarioTesting/data/HJ_bench/data/HJ_WEATHER_2002-2018_XC_mod.mdf'
    return binary_path, input_file, input_precip, input_met

def test_match_input_to_output(setup_data):
    binary_path, input_file, input_precip, input_met = setup_data
    try:
        print('Initializing model')
        m = model()
        m.read_input_file(input_file)

        # Add debugging statement before running the model
        print('Running tRIBS model')
        m.run(binary_path, input_file)
        # Add debugging statement after running the model
        print('Finished running tRIBS model')

        # Create instance of pytRIBS results class & read in results
        r = results(input_file)
        r.get_element_results()
        pixel = r.element[0]['pixel']

        # Read precipitation input
        p_in = m.read_precip_station(input_precip)
        # Drop the last hour from precipitation data
        # this would automatically fail the test if not in place because last hour is not recorded.
        # This should be considered a bug.
        p_in.drop(p_in.index[-1], inplace=True)

        # Read meteorological input
        met_in = m.read_met_station(input_met)
        # Drop the last hour from meteorological data
        met_in.drop(met_in.index[-1], inplace=True)

        # Compare precipitation inputs and outputs
        pdiff = p_in.R.values - pixel.Rain_mm_h.values
        assert not np.any(pdiff), "Difference found in precipitation forcing"

        # List of input variables and related variables from pixel file
        input_variables = ['PA', 'RH', 'XC', 'US', 'TA', 'IS']
        tribs_nodata = 9999.99
        pixel_variables = ['Press_Pa', 'RelHum_[]', 'SkyCov_[]', 'Wind_m_s', 'AirT_oC', 'ShrtRadIn_W_m2']

        # Compare meteorological inputs and outputs
        for invar, pvar in zip(input_variables, pixel_variables):
            if invar == 'XC' and any(met_in['XC'] == tribs_nodata):
                continue
            var_dif = np.not_equal(met_in[invar].values, pixel[pvar].values)

            assert not np.any(var_dif), f"Difference found in variables {invar} and {pvar} at time {pixel.Time_hr[(met_in[invar] != pixel[pvar]).idxmax()]}"

    except Exception as e:
        raise AssertionError(f"Test failed due to an unexpected error: {e}")












