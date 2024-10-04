from mofsbufixer.script import compile_sbu

def test_metal_sbu(filename):
    '''
    Test the metal_sbu script with a sample filename.
    '''
    compile_sbu.metal_sbu(filename)

test_metal_sbu('../../SBU_BU_with_dummy/')

#test_metal_sbu('../../SBU_BU_with_dummy/AVIMIC_organic_sbu_1.xyz')
# test_metal_sbu('../../SBU_BU_with_dummy/ABADUG_organic_sbu_1.xyz')
# test_metal_sbu('../../SBU_BU_with_dummy/AVIHIY_metal_sbu_1.xyz')
# test_metal_sbu('../../SBU_BU_with_dummy/ABOZUO_metal_sbu_1.xyz')