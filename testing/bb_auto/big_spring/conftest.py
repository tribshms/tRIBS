import sys
import os
import pytest
import getpass
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

extra_path = f"{os.environ['HOME']}/Documents/Repos/Forked/pytRIBS"
if extra_path not in sys.path:
    sys.path.append(extra_path)

from pytRIBS.classes import Model as model
from pytRIBS.classes import Results as results


@pytest.fixture(scope='module')
def setup_data():
    # setup conditional for refining or adding tests, for operational testing should always be set to True
    run_flag = True

    binary_par_path = '/Users/wr/Documents/Repos/Forked/tRIBS/cmake-build-parallel/tRIBSpar'
    binary_ser_path = '/Users/wr/Documents/Repos/Forked/tRIBS/cmake-build-serial/tRIBS'
    input_par_file = '/Users/wr/Desktop/big_spring/src/in_files/big_spring_par.in'
    input_ser_file = '/Users/wr/Desktop/big_spring/src/in_files/big_spring.in'

    if run_flag:
        # for future efficieny test, need to check data first
        #ref_data = '/Users/wr/Desktop/tScenarioTesting/data/HJ_bench/data/Snotel/bcqc_34.75000_-111.41000.txt'

        # Get the current username
        username = getpass.getuser()

        # Get the current date
        current_date = datetime.now()

        print(f'\ntRIBS Integrated Test: Big Spring Benchmark\nUser: {username} \nDate: {current_date} \n')

        # Use ThreadPoolExecutor to run tasks in parallel
        with ThreadPoolExecutor(max_workers=2) as executor:
            futures = {
                executor.submit(run_parallel_task, binary_par_path, input_par_file): 'parallel',
                executor.submit(run_serial_task, binary_ser_path, input_ser_file): 'serial'
            }

            results_dict = {}
            for future in as_completed(futures):
                task_name = futures[future]
                try:
                    mrf_result = future.result()
                    results_dict[task_name] = mrf_result
                    print(f'{task_name} simulation completed')
                except Exception as e:
                    print(f'{task_name} simulation generated an exception: {e}')

                # Access the results
        parallel_result = results_dict.get('parallel')
        serial_result = results_dict.get('serial')

    else:
        parallel_result = get_existing_results(input_par_file)
        serial_result = get_existing_results(input_ser_file)

    r_par = results(input_par_file)

    yield parallel_result, serial_result, r_par


def run_parallel_task(binary_path, input_file):
    m = model()
    m.read_input_file(input_file)
    m.run(binary_path, input_file, mpi_command='mpirun -n 3', verbose=False)
    r = results(input_file)
    r.get_mrf_results()
    return r.mrf['mrf']


def run_serial_task(binary_path, input_file):
    m = model()
    m.read_input_file(input_file)
    m.run(binary_path, input_file, verbose=False)
    r = results(input_file)
    r.get_mrf_results()
    return r.mrf['mrf']


def get_existing_results(input_file):
    r = results(input_file)
    r.get_mrf_results()
    return r.mrf['mrf']
