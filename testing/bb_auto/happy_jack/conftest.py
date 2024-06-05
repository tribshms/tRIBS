import sys
import os
import pytest
import getpass
from datetime import datetime

extra_path = f"{os.environ['HOME']}/Documents/Repos/Forked/pytRIBS"
if extra_path not in sys.path:
    sys.path.append(extra_path)

from pytRIBS.classes import Model as model
from pytRIBS.classes import Results as results


@pytest.fixture(scope='module')
def setup_data():
    binary_path = '/Users/wr/Documents/Repos/Forked/tRIBS/cmake-build-serial/tRIBS'
    input_file = '/Users/wr/Desktop/tScenarioTesting/scenarios/one/scenario_one.in'
    input_precip = '/Users/wr/Desktop/tScenarioTesting/data/HJ_bench/data/HJ_PRECIP_2002-2018.mdf'
    input_met = '/Users/wr/Desktop/tScenarioTesting/data/HJ_bench/data/HJ_WEATHER_2002-2018_XC_mod.mdf'
    ref_data = '/Users/wr/Desktop/tScenarioTesting/data/HJ_bench/data/Snotel/bcqc_34.75000_-111.41000.txt'

    # Get the current username
    username = getpass.getuser()

    # Get the current date
    current_date = datetime.now()

    print(f'\ntRIBS Integrated Test: Happy Jack Benchmark\nUser: {username} \nDate: {current_date} \n')

    m = model()
    m.read_input_file(input_file)

    m.run(binary_path, input_file, verbose=False)  # set to true if model doesn't seem to be running.

    # Create instance of pytRIBS results class & read in results
    r = results(input_file)
    r.get_element_results()
    pixel = r.element[0]['pixel']

    yield pixel, input_met, input_precip, m, r, ref_data
