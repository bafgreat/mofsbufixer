import requests
import re
import ase
import numpy as np
from collections import Counter
from mofstructure import mofdeconstructor
from mofsbufixer.input_output import coords_library
from mofsbufixer.input_output import filetyper
from mofsbufixer.sym import symmetrize
from mofsbufixer.data import inchi_name_data


all_metals = mofdeconstructor.transition_metals()


def find_key_by_value(data, target):
    """
    Find the key in a dictionary where the target value is located within the nested lists.
    **Parameters**
        - data: Dictionary with lists of lists as values.
        - target: Value to search for within the nested lists.

    **Returns:**
        - Key in the dictionary where the target value is located within the nested lists.
    """
    for key, value in data.items():
        if key == target:
            chemical_names = value['name']
            print (chemical_names)
            return chemical_names
    return None


def format_for_filename(chemical_name):
    """
    Formats the chemical name to be used in a file name by:
    - Replacing brackets '(', ')', '[', ']' with underscores '_'
    - Keeping hyphens '-'
    - Removing any other special characters except alphanumeric and underscores

    **Parameters**
        - chemical_name: str
            The chemical name to be formatted.

    **Returns**
        - formatted_name:
            The formatted chemical name suitable for a file name.
    """
    # Replace brackets with underscores
    formatted_name = re.sub(r"[(){}\[\]]", '', chemical_name)
    # Replace commas, spaces, and apostrophes with underscores
    formatted_name = re.sub(r"[\s']", '', formatted_name)
    # Remove any remaining non-alphanumeric characters except hyphens and underscores
    formatted_name = re.sub(r"[^a-zA-Z0-9-_]", '_', formatted_name)
    return formatted_name


def find_metal(ase_atoms):
    """
    Find the metal(s) in an ASE Atom object.

    **Parameters**
        - ase_atoms: ASE Atoms object
            The ASE Atom object to find the metal(s) in.

    **returns**
        str
            The metal(s) found in the ASE Atom object, formatted as a hyphen-separated string.

    """
    metals = sorted(list(set([i.symbol for i in ase_atoms if i.symbol in all_metals])))
    return '_'.join(metals)


def inchikey_to_name(inchi_key: str) -> str:
    """
    Converts an InChIKey to a chemical name using the PubChem PUG REST API,
    and formats the name for use in a file name.

    **Parameters**
        - inchi_key: str
            The InChIKey to be converted to a chemical name.

    **Returns**
        str
            The formatted chemical name for a file, or None if no match is found.
    """
    # chemical_name = find_key_by_value(inchi_name_data, inchi_key)
    # if chemical_name is not None:
    #     formatted_name = format_for_filename(chemical_name)
    #     print (formatted_name, chemical_name )
    #     return formatted_name, chemical_name

    # else:
    pubchem_url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchi_key}/property/IUPACName/JSON'
    response = requests.get(pubchem_url)

    if response.status_code == 200:
        data = response.json()
        try:
            chemical_name = data['PropertyTable']['Properties'][0]['IUPACName']
            formatted_name = format_for_filename(chemical_name)
            return formatted_name, chemical_name
        except (KeyError, IndexError):
            return "Chemical name not found."
    else:
        return f"Error: {response.status_code}"


def ase_to_inchi(ase_atom):
    """
    This function converts an ASE Atom object to an InChI string.

    **Parameters**
        - ase_atom: ASE Atom object

    **Returns**
        str: The InChI string of the ASE Atom object.
    """
    _, _, inChiKey = mofdeconstructor.compute_openbabel_cheminformatic(ase_atom)
    return inChiKey


def add_dummy(ase_atom, indices):
    """
    This function removes atoms at the specified indices from an ASE Atom object
    and adds dummy atoms ('X') at the end of the structure at the same positions.

    **Parameters**
        - ase_atom: ASE Atom object
        - indices: list
            The indices at which to replace atoms with dummy atoms.

    **Returns**
        - ASE Atom object with the dummy atoms added at the end.
    """

    # Remove the atoms at the specified indices
    mask = [i for i in range(len(ase_atom)) if i not in indices]
    ase_without_dummy = ase_atom[mask]
    atom_neighbors, _ = mofdeconstructor.compute_ase_neighbour(ase_atom)
    connected_conponents = mofdeconstructor.connected_components(atom_neighbors)

    # Add dummy atoms ('X') at the end of the structure, using the positions of removed atoms
    dummy_atoms = ase.Atoms(['X'] * len(indices), positions=ase_atom.get_positions()[indices])

    # Merge the structure without the original atoms and the new dummy atoms
    final_structure = ase_without_dummy + dummy_atoms

    return final_structure, connected_conponents


def non_metal_formula(ase_atom, all_metals):
    """
    This function takes an ASE Atom object and returns the chemical formula
    of only the non-metal elements, sorted alphabetically by element symbol,
    excluding metals defined in all_metals.

    **Parameters**
        - ase_atom: ASE Atom object
            The atomic structure to analyze.
        - all_metals: list
            List of metal element symbols to exclude from the formula.

    **Returns**
        str: The chemical formula for the non-metal elements, sorted alphabetically.
    """
    # Get the atomic symbols of the atoms in the ASE Atom object
    atomic_symbols = ase_atom.get_chemical_symbols()

    # Filter out the symbols that are in the list of metals
    non_metal_symbols = [symbol for symbol in atomic_symbols if symbol not in all_metals]

    # Count the occurrences of each non-metal symbol
    symbol_counts = Counter(non_metal_symbols)

    # Sort the symbols alphabetically
    sorted_symbols = sorted(symbol_counts.items())

    # Create the chemical formula string
    formula = ''.join(f"{symbol}{count if count > 1 else ''}" for symbol, count in sorted_symbols)

    return formula

def paddlewheel_extension_form(ase_atom, all_metals):
    """
    This function takes an ASE Atom object and returns the chemical formula
    of only the non-metal elements, sorted alphabetically by element symbol,
    excluding metals defined in all_metals and excluding carbon (C), oxygen (O), and hydrogen (H).

    **Parameters**
        - ase_atom: ASE Atom object
            The atomic structure to analyze.
        - all_metals: list
            List of metal element symbols to exclude from the formula.

    **Returns**
        str: The chemical formula for the non-metal elements, sorted alphabetically.
    """
    # Get the atomic symbols of the atoms in the ASE Atom object
    atomic_symbols = ase_atom.get_chemical_symbols()

    # Filter out the symbols that are in the list of metals or are C, O, or H
    non_metal_symbols = [
        symbol for symbol in atomic_symbols
        if symbol not in all_metals and symbol not in ['C', 'O', 'H']
    ]

    # Count the occurrences of each non-metal symbol
    symbol_counts = Counter(non_metal_symbols)

    # Sort the symbols alphabetically
    sorted_symbols = sorted(symbol_counts.items())

    # Create the chemical formula string
    formula = ''.join(f"{symbol}{count if count > 1 else ''}" for symbol, count in sorted_symbols)

    return formula

def sub_dummy_with_hydrogens(ase_atom):
    """
    This function replaces dummy atoms ('X') in an ASE Atom object with hydrogen atoms ('H'),
    and returns the indices of the dummy atoms and the modified ASE Atom object.

    **Parameters**
        - ase_atom: ASE Atom object

    **Returns**
        - list: A list of indices of the dummy atoms that were replaced with hydrogen atoms,
        - ase_atom:and the modified ASE Atom object.
    """
    # Find the indices of dummy atoms ('X') in the ASE Atom object and replace them with hydrogen atoms ('H')
    # Also, update the symbol of the dummy atoms to 'H' in the ASE Atom object for consistency with other scripts and functions.
    x_indices = []
    for atoms in ase_atom:
        if atoms.symbol == 'X':
            x_indices.append(atoms.index)
            atoms.symbol = 'H'
    return x_indices, ase_atom


def log_failure(filename, outfile):
    """
    This function logs failures to a file called 'Failed.txt'.
    """
    with open(f"{outfile}.txt", "a") as file:
        file.write(f"{filename}\n")

def check_same_x_positions(ase_atom):
    """
    Check if all X positions are the same for atoms labeled as 'X' in the ASE Atoms object.

    **Parameters**
        - ase_atom: ASE Atoms object
            The atomic structure to analyze.

    **Returns**
        bool: True if all X positions are the same, False otherwise.
    """
    # Get atomic symbols and positions from ASE Atoms object
    symbols = ase_atom.get_chemical_symbols()
    positions = ase_atom.get_positions()

    # Extract X positions of atoms labeled 'X'
    x_positions = [pos[0] for symbol, pos in zip(symbols, positions) if symbol == 'X']

    # Check if all X positions are the same using numpy's isclose to handle floating-point precision
    return np.all(np.isclose(x_positions, x_positions[0]))


def sbu_compiler(filename, metal_sbu_path, organic_sbu_path):
    """
    This function reads a metal SBU file, extracts the coordinates and symmetry information,
    and returns them in the form of a dictionary.

    **Parameters**
        - filename: str
            The file path to the metal SBU file.

    **Returns**
        - dict: A dictionary containing the coordinates and symmetry information.
    """
    base_name = filename.split('/')[-1].split('.')[0].split('_')[0]
    ase_atom = coords_library.read_ase(filename)

    sbu_type = ase_atom.info.get('sbu_type', None)
    x_indices, ase_atom = sub_dummy_with_hydrogens(ase_atom)
    if len(x_indices) == 0:
        log_failure(filename, 'mofstructure_error')
        return
    else:
        try:
            sym_mol, point = symmetrize.symmetrize_molecule(ase_atom)
            new_atom, connected_conponents = add_dummy(sym_mol, x_indices)
            check_x = check_same_x_positions(new_atom)
            print (check_x)
            if len(connected_conponents) == 1 and not check_x:
                new_atom.info['refcode'] = base_name
                if sbu_type is not None and len(x_indices) < 25:
                    metal = find_metal(sym_mol)
                    if sbu_type != 'still checking!':
                        sbu_name = f'{metal}_{sbu_type}_X{len(x_indices)}.xyz'
                    elif 'paddlewheel' in sbu_type:
                        formular = paddlewheel_extension_form(ase_atom, all_metals)
                        sbu_name = f'{metal}_{sbu_type}_{formular}_X{len(x_indices)}.xyz'
                    else:
                        formular = non_metal_formula(ase_atom, all_metals)
                        sbu_name = f'{metal}_{formular}_X{len(x_indices)}.xyz'
                    new_atom.write(f'{metal_sbu_path}/{sbu_name}')
                else:
                    inchi_key = ase_to_inchi(sym_mol)
                    iupac_name = inchikey_to_name(inchi_key)
                    if len(iupac_name) == 2:
                        print (iupac_name[0], iupac_name[1])
                        new_atom.info['iupac_name'] = iupac_name[1]
                        sbu_name = f'{iupac_name[0]}_X{len(x_indices)}.xyz'
                        new_atom.write(f'{organic_sbu_path}/{sbu_name}')
            else:
                print(f"Error processing {filename}: More than one connected component")
                log_failure(filename, 'mofstructure_error_connected_component')
        except Exception as e:
            print(f"Error processing {filename}: {e}")
            log_failure(filename, 'Failed')
    return