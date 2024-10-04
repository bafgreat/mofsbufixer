from mofsbufixer.script import compile_sbu

def test_metal_sbu(filename):
    '''
    Test the metal_sbu script with a sample filename.
    '''
    compile_sbu.metal_sbu(filename)
meta_sbu_path = '../../SCM_SBU/XYZ/metal_sbu'
meta_organic_path = '../../SCM_SBU/XYZ/organic_sbu'

test_metal_sbu('../../SBU_BU_with_dummy/ADIQAH_organic_sbu_1.xyz', meta_sbu_path, meta_organic_path)

#test_metal_sbu('../../SBU_BU_with_dummy/AVIMIC_organic_sbu_1.xyz')
# test_metal_sbu('../../SBU_BU_with_dummy/ABADUG_organic_sbu_1.xyz')
# test_metal_sbu('../../SBU_BU_with_dummy/AVIHIY_metal_sbu_1.xyz')
# test_metal_sbu('../../SBU_BU_with_dummy/ABOZUO_metal_sbu_1.xyz')