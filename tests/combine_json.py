import glob
from mofsbufixer.input_output import filetyper

all_json = sorted(glob.glob('../matrix_data2/*json'))

data_1 = filetyper.load_data(all_json[0])

for data in all_json[1:]:
    data_2 = filetyper.load_data(data)
    data_1.update(data_2)

filetyper.write_json(data_1, './adjacent_matrix.json')