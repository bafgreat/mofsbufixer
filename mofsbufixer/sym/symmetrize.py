

from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, PointGroupAnalyzer


def symmetrize_molecule(ase_atom):
    '''
    Function to symmetrize a molecule. The function takes an
    ase atom object and converts it to pymatgen molecule object
    and then perform the symmetry transformation and returns both
    the symmetrized molecule and the point group.

    **parameter**
    ase_atom : ASE Atom object
    **returns**
    ase_atoms : ASE Atom object
    point_group : str

    '''
    molecule = AseAtomsAdaptor.get_molecule(ase_atom)
    symmetrizer = PointGroupAnalyzer(molecule)
    point_group = symmetrizer.get_pointgroup()
    symmetrized = symmetrizer.symmetrize_molecule()
    new_molecule = symmetrized['sym_mol']
    ase_atoms = AseAtomsAdaptor.get_atoms(new_molecule)
    return ase_atoms, point_group


def symmetrize_crystal(ase_atom):
    '''
    Function to symmetrize a crystal structure. The function takes an
    ASE Atom object representing a crystal structure, converts it to
    a pymatgen Structure object, performs symmetry analysis, and returns
    the symmetrized crystal structure and the space group.

    **parameter**
    ase_atom : ASE Atom object (crystal structure)

    **returns**
    ase_atoms : ASE Atom object (symmetrized crystal structure)
    space_group : str
    '''
    structure = AseAtomsAdaptor.get_structure(ase_atom)

    symmetrizer = SpacegroupAnalyzer(structure)
    space_group = symmetrizer.get_space_group_symbol()

    symmetrized_structure = symmetrizer.get_symmetrized_structure()

    ase_atoms = AseAtomsAdaptor.get_atoms(symmetrized_structure)

    return ase_atoms, space_group

