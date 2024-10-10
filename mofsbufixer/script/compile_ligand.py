import re
import os
import shutil
import glob
import ase
import numpy as np
from mofstructure import mofdeconstructor
from mofsbufixer.input_output import coords_library, filetyper

def make_directory(dir_name):
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

def log_inchikey(data, outfile):
    """
    This function logs failures to a file called 'Failed.txt'.
    """
    with open(f"{outfile}.txt", "a") as file:
        file.write(f"{data}\n")


def ligand_compiler(all_filename, path_to_save, added_hydrogen):
    if os.path.exists('all_inchikey.txt'):
        all_inchi_keys = filetyper.get_contents('all_inchikey.txt')
        all_inchi_keys = [i.strip() for i in all_inchi_keys]

    else:
        all_inchi_keys = []
    for filename in all_filename:
        ase_atom = coords_library.read_ase(filename)
        inchikey = ase_atom.info['inchikey']
        all_dir = f'{path_to_save}/{inchikey}'
        added_hydrogen_dir = f'{added_hydrogen}/{inchikey}'
        if inchikey not in all_inchi_keys:
            print (inchikey)
            make_directory(all_dir)
            make_directory(added_hydrogen_dir)
            shutil.copy(filename, f'{all_dir}/{os.path.basename(filename)}')
            pybel_mol = mofdeconstructor.ase_2_pybel(ase_atom)
            pybel_mol.addh()
            full_path = f'{added_hydrogen_dir}/{os.path.basename(filename)}'
            h_index = full_path.rindex('.')
            add_h_name = full_path[:h_index] + ".sdf"
            pybel_mol.write("sdf", add_h_name, overwrite=True)

        else:

            shutil.copy(filename, f'{all_dir}/{os.path.basename(filename)}')
        all_inchi_keys.append(inchikey)
        log_inchikey(inchikey, 'all_inchikey')

    # print (inchikey)


