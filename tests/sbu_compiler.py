from mofsbufixer.script import compile_sbu
from mofsbufixer.input_output import filetyper
import glob

def test_metal_sbu(filename, meta_sbu_path, organic_sbu_path):
    '''
    Test the metal_sbu script with a sample filename.
    '''
    compile_sbu.sbu_compiler(filename,  meta_sbu_path, organic_sbu_path)
meta_sbu_path = '../../SCM_SBU/XYZ/metal_sbu'
organic_sbu_path = '../../SCM_SBU/XYZ/organic_sbu'

data1 = list(filetyper.get_contents('./Failed.txt'))
data2 = list(filetyper.get_contents('./mofstructure_error.txt'))
data3 = list(filetyper.get_contents('./mofstructure_error_connected_component.txt'))
data = list(set(data1+data2+data3))
data = [path.strip() for path in data]


# filenames = '../../SBU_BU_with_dummy/ABADUG_organic_sbu_1.xyz'

all_files = sorted(glob.glob('../../SBU_BU_with_dummy/*xyz'))
for filenames in all_files:
    if filenames not in data:
        print (f'processing file:{filenames}')
        test_metal_sbu(filenames, meta_sbu_path, organic_sbu_path)

# test_metal_sbu('../../SBU_BU_with_dummy/AVIMIC_organic_sbu_1.xyz', meta_sbu_path, organic_sbu_path)
# test_metal_sbu('../../SBU_BU_with_dummy/ABADUG_organic_sbu_1.xyz')
# test_metal_sbu('../../SBU_BU_with_dummy/AVIHIY_metal_sbu_1.xyz')
# test_metal_sbu('../../SBU_BU_with_dummy/ABOZUO_metal_sbu_1.xyz')