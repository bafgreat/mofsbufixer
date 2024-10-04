from mofsbufixer.script import metal_sbu

def test_metal_sbu(filename):
    '''
    Test the metal_sbu script with a sample filename.
    '''
    metal_sbu.run(filename)


test_metal_sbu(filename)