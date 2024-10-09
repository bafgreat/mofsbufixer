import pytest
from mofsbufixer.script import compile_sbu
from mofsbufixer.input_output import coords_library


@pytest.fixture
def ase_atoms():
    # Ensure the path to the test data is correct
    return coords_library.read_ase('tests/test_data/ATOBAN_metal_sbu_1.xyz')


def test_metal_sbu(ase_atoms):
    """
    Test the organic_sbu function to manipulate the SBU structure.
    """
    data = compile_sbu.manipulate_metal_sbu(ase_atoms)
    assert len(data[0]) == 20
    assert data[1] == 6
