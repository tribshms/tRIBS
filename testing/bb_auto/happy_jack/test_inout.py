# This test confirms if the inputs for a provided model match the outputs.
# Currently this only setup for a single element case, but will need to be updated to address more complex multi-element
# models.

#pytest --binary /path/to/your/binary/executable --input /path/to/input/file


import sys
import os

import requests

extra_path = f"{os.environ['HOME']}/Documents/Repos/Forked/pytRIBS"
if extra_path not in sys.path:
    sys.path.append(extra_path)

from pytRIBS.classes import Model as model
from pytRIBS.classes import Results as results


def test_match_input_to_output(binary_path, input_file, input_precip, input_met):
    # create instance of pytRIBS model and set option based on input file
    m = model()
    m.read_input_file(input_file)

    # run tRIBS w/ inputs from input_file
    m.run(binary_path, input_file, verbose=False)

    # create instance of pytRIBS results class & read in results
    r = results(input_file)







