# conftest.py
import pytest


@pytest.fixture
def binary_path():
    return "/path/to/your/binary/executable"

@pytest.fixture
def input_file():
    return "/path/to/input/file"

@pytest.fixture
def input_precip():
    return "/path/to/input/file"

@pytest.fixture
def input_met():
    return "/path/to/input/file"