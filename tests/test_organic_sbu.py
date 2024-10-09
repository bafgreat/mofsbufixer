import pytest
from mofsbufixer.script import compile_sbu
from mofsbufixer.input_output import coords_library


@pytest.fixture
def ase_atoms():
    # Ensure the path to the test data is correct
    return coords_library.read_ase('tests/test_data/test_pyridine_organic_sbu.xyz')


def test_organic_sbu(ase_atoms):
    """
    Test the organic_sbu function to manipulate the SBU structure.
    """
    data = compile_sbu.manipulate_organic_sbu(ase_atoms)
    assert data[0] == '4-pyridin-4-ylpyridine'
    assert data[2] == 2
    assert len(data[1]) == 22
