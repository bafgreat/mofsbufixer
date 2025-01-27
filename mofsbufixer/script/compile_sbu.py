import requests
import re
import ase
import numpy as np
from collections import Counter
from ase.geometry import get_angles
from mofstructure import mofdeconstructor
from mofsbufixer.input_output import coords_library
from mofsbufixer.sym import symmetrize


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
    metals = [i for i in metals if i not in ['X']]
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
        if symbol not in all_metals and symbol not in ['C', 'O', 'H', 'X']
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


def remove_item_from_list(list1, list2):
    """
    This function removes items from one list that are also present in another list.

    **Parameters**
        - list1: list
        - list2: list

    **Returns**
        - new_list: list
        - containing the items that are not present in list2.
    """
    new_list = []
    for i in list1:
        if i not in list2:
            new_list.append(i)
    return new_list

def set_distance(ase_atom, index1, index2, distance):
    """
    Sets the distance between two atoms in an ASE Atoms object,
    only moving the atom at index2.

    **Parameters:**
        - ase_atom(Atoms): ASE Atoms object.
        - index1 (int): Index of the first atom (reference atom).
        - index2 (int): Index of the second atom (the atom to move).
        - distance (float): The target distance to set between the two atoms.

    **Returns:**
        - ase_atom: modified atom
        """
    # Get the positions of the two atoms
    pos1 = ase_atom.positions[index1]
    pos2 = ase_atom.positions[index2]

    # Calculate the current distance vector and its magnitude
    current_vector = pos2 - pos1
    current_distance = np.linalg.norm(current_vector)
    if current_distance == 0:
        current_distance = 1e-10

    # Normalize the vector to get the direction
    direction = current_vector / current_distance

    # Set the new position for the atom at index2
    new_pos2 = pos1 + direction * distance
    ase_atom.positions[index2] = new_pos2
    return ase_atom


def find_fourth_point(p1, p0, p2, distance=1.2):
    """
    Calculate a point P_new such that:
    - It lies in the same plane as the points P1, P0, and P2.
    - It is at a specified distance from P0 (default: 1.2).
    - It is in the direction opposite to the angle formed by P1, P0, and P2.

    **Parameters**:
        p1 (array-like): Coordinates of the first point (P1) as [x, y, z].
        p0 (array-like): Coordinates of the second point (P0) as [x, y, z].
        p2 (array-like): Coordinates of the third point (P2) as [x, y, z].
        distance (float): The distance from P0 to the new point (P_new). Default is 1.2.

    **Returns**:
        np.ndarray: Coordinates of the new point P_new as [x, y, z].
    """
    p1 = np.array(p1)
    p0 = np.array(p0)
    p2 = np.array(p2)

    # Step 1: Compute the normal vector to the plane defined by P1, P0, P2
    v1 = p1 - p0
    v2 = p2 - p0
    normal = np.cross(v1, v2)
    normal_unit = normal / np.linalg.norm(normal)

    # Step 2: Compute the bisector of the angle at P0 between P1 and P2
    v1_unit = v1 / np.linalg.norm(v1)
    v2_unit = v2 / np.linalg.norm(v2)
    bisector = v1_unit + v2_unit
    bisector_unit = bisector / np.linalg.norm(bisector)

    # Step 3: Reverse the bisector direction to point opposite to the angle
    opposite_dir = -bisector_unit

    # Step 4: Project the opposite direction onto the plane P1, P0, P2
    projection = opposite_dir - np.dot(opposite_dir, normal_unit) * normal_unit
    projected_dir_unit = projection / np.linalg.norm(projection)

    # Step 5: Calculate the new point P_new
    p_new = p0 + distance * projected_dir_unit

    return p_new

def find_fourth_point_previous(p1, p0, p2, distance=1.8):
    """
    Find a fourth point that lies in the plane of three points and is opposite
    the angle formed by vectors p1-p0 and p2-p0.

    Args:
        p1 (array-like): Coordinates of point 1.
        p0 (array-like): Coordinates of point 0 (origin point).
        p2 (array-like): Coordinates of point 2.
        distance (float): Distance from point 0 to the fourth point.

    Returns:
        numpy.ndarray: Coordinates of the fourth point.
    """
    p1 = np.array(p1)
    p0 = np.array(p0)
    p2 = np.array(p2)

    # Vectors defining the plane
    v1 = p1 - p0
    v2 = p2 - p0

    # Normal vector to the plane
    normal = np.cross(v1, v2)
    normal /= np.linalg.norm(normal)

    # Normalize v1 and v2
    v1_norm = v1 / np.linalg.norm(v1)
    v2_norm = v2 / np.linalg.norm(v2)

    # Bisector of the angle
    bisector = v1_norm + v2_norm
    bisector /= np.linalg.norm(bisector)

    # Reflect the bisector vector across the plane
    reflection = bisector - 2 * np.dot(bisector, normal) * normal

    # Negate the reflection to ensure it's on the opposite side of the angle
    reflection_opposite = -reflection

    # Scale the reflected vector to the desired distance and add to p0
    p3 = p0 + reflection_opposite * distance

    return p3


def find_correct_coordinats(sym_atom, atom_index, mapper):
    atom_neighbors, _ = mofdeconstructor.compute_ase_neighbour(sym_atom)
    new_positions = []
    done_index = []

    for i in atom_index:
        n_atom = mapper.get(i)
        connected_atoms = atom_neighbors[n_atom]
        connected_atoms = [x for x in connected_atoms if x not in atom_index]
        if len(connected_atoms):
            p1, p2 = connected_atoms
            coord_1  = sym_atom[p1].position
            coord_2  = sym_atom[p2].position
            coord_n = sym_atom[n_atom].position
            p_new = find_fourth_point(coord_1, coord_n, coord_2, 1.2)
            new_positions.append(p_new)
            done_index.append(i)

    if len(done_index) == 1:
        return ase.Atom(symbol='X' , position=new_positions[0])

    elif len(done_index) > 1:
        return ase.Atoms(['X'] * len(done_index), positions=new_positions)




def manipulate_organic_sbu(ase_atom, distance=1.2):
    """
    Manipulates organic secondary building units (SBUs) by identifying and transforming dummy atoms,
    symmetrizing the structure, and generating chemical identifiers. The function focuses on handling
    nitrogen-dummy atom interactions and adjusting distances between atoms.

    **Procedure:**
    1. Identify all dummy atoms (marked as 'X') in the atomic structure.
    2. Determine if any dummy atom is directly connected to a nitrogen atom (symbol 'N').
       - Nitrogen atoms in MOFs typically form dative covalent bonds with central metals,
         and in these cases, the dummy atom should be removed before symmetrization.
    3. For each dummy atom, identify the atoms connected to it based on proximity (distances).
       - Store the index of the atom directly connected to the dummy atom to help reassign
         atom distances later.
    4. Remove dummy atoms connected to nitrogen from the symmetrization process,
       while keeping track of their presence for later use.
    5. Convert the remaining dummy atoms into hydrogen atoms ('H') to facilitate symmetrization.
    6. Perform molecular symmetrization on the modified structure (with hydrogen atoms).
    7. Generate the IUPAC name of the symmetrized molecule using its InChIKey.
    8. After symmetrization, replace the hydrogen atoms back to dummy atoms ('X').
    9. Adjust the distances between atoms to the specified distance (default 1.2 Å),
       particularly between the dummy atoms and their mapped neighbors.
    10. Return both the IUPAC name and the final modified ASE Atoms object with symmetrized structure
        and adjusted atom distances.

    **Parameters:**
    - ase_atom: ASE Atoms object
        The atomic structure to manipulate. The structure should contain dummy atoms ('X')
        representing placeholders in the structure that will either be replaced or removed
        depending on their interaction with nitrogen atoms.
    - distance: float, optional (default=1.2)
        The target distance (in angstroms) between dummy atoms and their second closest neighboring atoms
        after manipulation. This is crucial for correctly placing the atoms during structure optimization.

    **Returns:**
    - iupac_name: str
        The IUPAC name of the symmetrized molecule, generated using its InChIKey.
    - edited_atom: ASE Atoms object
        The modified ASE Atoms object, with dummy atoms adjusted, symmetrized, and distances set accordingly.
    """
    # Define the indices of dummy atoms ('X') and hydrogen atoms ('H')
    x_indices = []
    nitrongen_x = []
    mapper = {}
    # com = ase_atom.get_center_of_mass()
    # atom_neighbors, _ = mofdeconstructor.compute_ase_neighbour(ase_atom)
    smi = ase_atom.info.get('smi')
    inchikey = ase_atom.info.get('inchikey')
    for atoms in ase_atom:
        if atoms.symbol == 'X':
            x_indices.append(atoms.index)
    for x in x_indices:
        distances = ase_atom.get_distances(x, range(len(ase_atom)))
        min_index = np.argmin(distances)
        distances[min_index] = np.inf
        second_min_index = np.argmin(distances)
        # neighbour = atom_neighbors[second_min_index]
        # neighbour = [i for i in neighbour if i!=x]
        # print (f"neighbour:{neighbour}")
        mapper[x] = second_min_index
        if  ase_atom[second_min_index].symbol in ['N']:
            nitrongen_x.append(x)
        # elif len(neighbour)==2 and ase_atom[second_min_index].symbol in ['O', 'S', 'Se']:
        #     nitrongen_x.append(x)

    if len(x_indices) > 0:
        no_n_indices = remove_item_from_list(x_indices, nitrongen_x)
        fixer_atom = ase_atom.copy()
        all_indices_with_no_n = [i for i in range(len(ase_atom)) if i not in nitrongen_x]
        for i in no_n_indices:
            fixer_atom[i].symbol = 'H'
        fixer_atom = fixer_atom[all_indices_with_no_n]
        sym_mol, _ = symmetrize.symmetrize_molecule(fixer_atom)
        # sym_mol.positions -= com
        inchi_key = ase_to_inchi(sym_mol)
        iupac_name = inchikey_to_name(inchi_key)
        # print ('This is mapper', mapper, nitrongen_x)
        nitrogen_x_atoms = find_correct_coordinats(sym_mol, nitrongen_x, mapper)
        if nitrogen_x_atoms is not None:
            edited_atom = sym_mol + nitrogen_x_atoms
        for i in no_n_indices:
            edited_atom[i].symbol = 'X'
        for i in mapper:
            atom_index = mapper[i]
            edited_atom = set_distance(edited_atom, atom_index, i, 1.2)
        if len(iupac_name) == 2:
            edited_atom.info['iupac_name'] = iupac_name[1]
            edited_atom.info['smile'] = smi
            edited_atom.info['inchikey'] = inchikey

            return iupac_name[0], edited_atom, len(x_indices)
    return None


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


def check_organic_in_filename(filename):
    """
    This function checks if the word 'organic' is present in the given filename.

    **Parameters:**
    - filename: str
        The path or name of the file to check.

    **Returns:**
    - bool: True if 'organic' is found in the filename, False otherwise.
    """
    # Use regex to search for 'organic' in the filename
    if re.search(r'organic', filename, re.IGNORECASE):
        return True
    else:
        return False


def check_metal_in_filename(filename):
    """
    This function checks if the word 'organic' is present in the given filename.

    **Parameters:**
    - filename: str
        The path or name of the file to check.

    **Returns:**
    - bool: True if 'organic' is found in the filename, False otherwise.
    """
    # Use regex to search for 'organic' in the filename
    if re.search(r'metal', filename, re.IGNORECASE):
        return True
    else:
        return False

def check_sbu_type(sbu_type, checker):
    """
    This function checks if the word 'organic' is present in the given filename.

    **Parameters:**
    - filename: str
        The path or name of the file to check.

    **Returns:**
    - bool: True if 'organic' is found in the filename, False otherwise.
    """
    # Use regex to search for 'organic' in the filename
    if re.search(rf'{checker}', sbu_type, re.IGNORECASE):
        return True
    else:
        return False


def manipulate_metal_sbu(ase_atom, distance=1.2):

    x_indices = []
    mapper = {}
    atom_neighbors, _ = mofdeconstructor.compute_ase_neighbour(ase_atom)

    for atoms in ase_atom:
        if atoms.symbol == 'X':
            x_indices.append(atoms.index)
    for x in x_indices:
        distances = ase_atom.get_distances(x, range(len(ase_atom)))
        min_index = np.argmin(distances)
        distances[min_index] = np.inf
        second_min_index = np.argmin(distances)
        mapper[x] = second_min_index

    if len(x_indices) > 0:
        fixer_atom = ase_atom.copy()
        for i in x_indices:
            fixer_atom[i].symbol = 'H'
        sym_mol, _ = symmetrize.symmetrize_molecule(fixer_atom)

        edited_atom = sym_mol
        for i in x_indices:
            edited_atom[i].symbol = 'X'
        mask = [i for i in range(len(ase_atom)) if i not in x_indices]
        sbu_atom = edited_atom[mask]
        graph, _ = mofdeconstructor.compute_ase_neighbour(sbu_atom)
        connected_conponents = mofdeconstructor.connected_components(graph)
        for i in mapper:
            atom_index = mapper[i]
            edited_atom = set_distance(edited_atom, atom_index, i, 1.2)
        if len(connected_conponents) == 1:
            return edited_atom, len(x_indices)

    return None


def remove_duplicate_atoms(atoms, tolerance=1e-5):
    """
    Removes duplicate atoms from an ASE Atoms object based on positions.

    Parameters:
    atoms (Atoms): ASE Atoms object from which to remove duplicates.
    tolerance (float): Tolerance for considering two positions identical. Defaults to 1e-5.

    Returns:
    Atoms: ASE Atoms object with duplicate positions removed.
    """
    unique_positions = []
    unique_indices = []

    positions = atoms.get_positions()
    for i, pos in enumerate(positions):
        if not any(np.allclose(pos, p, atol=tolerance) for p in unique_positions):
            unique_positions.append(pos)
            unique_indices.append(i)

    return atoms[unique_indices]


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
    ase_atom.pbc = False
    ase_atom.cell = None

    # Assuming these functions are defined elsewhere in your code
    is_organic = check_organic_in_filename(filename)
    is_metal = check_metal_in_filename(filename)

    if is_organic:
        try:
            organic_data = manipulate_organic_sbu(ase_atom, distance=1.2)
            if organic_data is not None:
                iupac_name, new_atom, length_x = organic_data
                sbu_name = f'{iupac_name}_X{length_x}.xyz'
                new_atom.write(f'{organic_sbu_path}/{sbu_name}')
        except Exception as e:
            print(f"Error occurred while manipulating organic SBU: {e}")
            log_failure(filename, 'organic_sbu_error')
    elif is_metal:
        try:
            sbu_type = ase_atom.info.get('sbu_type', None)
            is_rod = check_sbu_type(sbu_type, 'rodlike')

            if is_rod:
                return
            else:
                ase_atom = remove_duplicate_atoms(ase_atom, tolerance=1e-5)
                metal_data = manipulate_metal_sbu(ase_atom)
                is_paddlewheel = check_sbu_type(sbu_type, 'paddlewheel')

                if metal_data is not None:
                    edited_atom, length_x = metal_data
                    metal = find_metal(edited_atom)
                    if is_paddlewheel:
                        formular = paddlewheel_extension_form(edited_atom, metal)
                        if len(formular) == 0:
                            sbu_name = f'{metal}_{sbu_type}_X{length_x}.xyz'
                        else:
                            sbu_name = f'{metal}_{formular}_{sbu_type}_X{length_x}.xyz'
                    elif sbu_type != 'still checking!':
                        sbu_name = f'{metal}_{sbu_type}_X{length_x}.xyz'
                    else:
                        formular = non_metal_formula(edited_atom, metal)
                        sbu_name = f'{metal}_{formular}_X{length_x}.xyz'

                    # Save the edited structure
                    edited_atom.write(f'{metal_sbu_path}/{sbu_name}')

        except Exception as e:
            print(f"Error occurred while manipulating metal SBU: {e}")
            log_failure(filename, 'metal_sbu_error')





    # sbu_type = ase_atom.info.get('sbu_type', None)
    # x_indices, ase_atom = sub_dummy_with_hydrogens(ase_atom)
    # if len(x_indices) == 0:
    #     log_failure(filename, 'mofstructure_error')
    #     return
    # else:
    #     try:
    #         if

def sbu_compiler2(filename, metal_sbu_path, organic_sbu_path):
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
            overlap = mofdeconstructor.inter_atomic_distance_check(new_atom)
            print (check_x)
            if len(connected_conponents) == 1 and overlap and not check_x:
                new_atom.info['refcode'] = base_name
                if sbu_type is not None and len(x_indices) < 25:
                    metal = find_metal(sym_mol)
                    if sbu_type != 'still checking!':
                        sbu_name = f'{metal}_{sbu_type}_X{len(x_indices)}.xyz'
                    elif 'paddlewheel' in sbu_type:
                        formular = paddlewheel_extension_form(new_atom, all_metals)
                        sbu_name = f'{metal}_{sbu_type}_{formular}_X{len(x_indices)}.xyz'
                    else:
                        formular = non_metal_formula(new_atom, all_metals)
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



def add_atoms_to_molecule(molecule, atom_data):
    """
    Adds multiple atoms to specific atoms in a molecule.

    Parameters:
        molecule (pybel.Molecule): The Pybel molecule object.
        atom_data (list of tuples): List of tuples where each tuple contains:
            - atom_idx (int): The index of the atom to which the new atom will be attached (1-based).
            - new_atom_symbol (str): Symbol of the new atom to add (e.g., "H").
            - bond_order (int, optional): Bond order for the new bond (default is 1).

    Returns:
        pybel.Molecule: The modified molecule.
    """
    # Access the underlying OBMol object
    obmol = molecule.OBMol

    for data in atom_data:
        # Unpack the data
        atom_idx, new_atom_symbol, *bond_order = data
        bond_order = bond_order[0] if bond_order else 1  # Default bond order is 1

        # Get the target atom by index
        target_atom = obmol.GetAtom(atom_idx)
        if not target_atom:
            raise ValueError(f"Atom index {atom_idx} is invalid.")

        # Create a new atom
        new_atom = OBAtom()
        new_atom.SetAtomicNum(pybel.ob.GetAtomicNum(new_atom_symbol))
        obmol.AddAtom(new_atom)

        # Create a bond between the target atom and the new atom
        bond_success = obmol.AddBond(target_atom.GetIdx(), new_atom.GetIdx(), bond_order)
        if not bond_success:
            raise RuntimeError(f"Failed to add a bond between atom {atom_idx} and the new atom.")

    # Update the molecule and return it
    molecule.OBMol = obmol
    return pybel.Molecule(obmol)

# # Example usage
# mol = pybel.readstring("smi", "CCO")  # Example molecule (ethanol)

# # List of atoms to add: [(target_index, new_atom_symbol, bond_order)]
# atoms_to_add = [
#     (2, "H", 1),  # Add a hydrogen atom to the second atom (carbon)
#     (3, "Cl", 1)  # Add a chlorine atom to the third atom (oxygen)
# ]

# modified_mol = add_atoms_to_molecule(mol, atoms_to_add)

# # Output the modified molecule
# print(modified_mol.write("smi"))
