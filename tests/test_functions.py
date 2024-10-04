import pytest
from mofsbufixer.script import compile_sbu
from mofsbufixer.input_output import coords_library

iupac_name = "1,1'-biphenyl"

@pytest.fixture
def organic_sbu():
    # Ensure the path to the test data is correct
    return coords_library.read_ase('tests/test_data/ABADUG_organic_sbu_1.xyz')

@pytest.fixture
def metal_sbu():
    # Ensure the path to the test data is correct
    return coords_library.read_ase('tests/test_data/AVIHIY_metal_sbu_1.xyz')

def test_format_for_filename():
    '''
    Ensure that IUPAC names are properly formatted
    to be used as filenames.

    '''
    formatted_name = compile_sbu.format_for_filename(iupac_name)
    assert formatted_name == "1_1-biphenyl"


def test_read_ase(organic_sbu, metal_sbu):
    '''
    Test reading of SBU structures from file using read_ase.
    '''
    assert organic_sbu is not None, "Organic SBU could not be read."
    assert metal_sbu is not None, "Metal SBU could not be read."

def test_inchikey(organic_sbu, metal_sbu):
    '''
    Test InChIKey conversion using the PubChem PUG REST API.
    '''
    x_indices_org, organic_h_sbu = compile_sbu.sub_dummy_with_hydrogens(organic_sbu)
    x_indices_met, metal_h_sbu = compile_sbu.sub_dummy_with_hydrogens(metal_sbu)
    organic_inchi_key = compile_sbu.ase_to_inchi(organic_h_sbu )
    formatted_name, iupac_nanme = compile_sbu.inchikey_to_name(organic_inchi_key)

    assert len(x_indices_org) == 2
    assert organic_inchi_key == 'ZUOUZKKEUPVFJK-UHFFFAOYSA-N'
    assert formatted_name == '1_1-biphenyl'
    assert iupac_nanme == "1,1'-biphenyl"



    print (organic_inchi_key ,len(x_indices_org), len(x_indices_met))