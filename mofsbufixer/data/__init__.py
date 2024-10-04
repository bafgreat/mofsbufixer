import os
from mofsbufixer.input_output import filetyper
path_dir = os.path.dirname(__file__)
json_file = os.path.abspath(f'{path_dir}/inchi_to_name_and_smile.json')

inchi_name_data = filetyper.load_data(json_file )