from __future__ import print_function
import numpy
from scm.autografs.utils.sbu import SBU 
from ase.io import iread
import os
import sys



if len(sys.argv) == 2:
    qcin = sys.argv[1]
else:
    print ('Incorrect filetype')
    sys.exit()

def xyz_to_inp(fxyz):
    for atoms in iread(fxyz):
        if "name" in atoms.info:
            name = atoms.info["name"]
        else:
            name = fxyz.split('.')[0]
        mol = SBU(name = name, atoms = atoms)
        with open(str(name)+'.inp', 'w') as of:
            of.write("Data: shape = {0}\n".format(mol.shape))
            of.write("Data: name  = {0}\n".format(mol.name))
            of.write("GEOMETRY CARTESIAN\n")
            for atom in mol.atoms:
                bondIndices = numpy.array(numpy.where(mol.bonds[atom.index]!=0)).ravel()
                bonds = [mol.bonds[atom.index, i] for i in bondIndices]
                sym = atom.symbol + " "*(5-len(atom.symbol))
                p0str = str(numpy.around(atom.position[0],decimals=6))
                p0 = p0str + " "*(10-len(p0str))
                p1str = str(numpy.around(atom.position[1],decimals=6))
                p1 = p1str + " "*(10-len(p1str))
                p2str = str(numpy.around(atom.position[2],decimals=6))
                p2 = p2str + " "*(10-len(p2str))
                of.write("{0}  {1}  {2}  {3}   MMTYPE={4}  QMMM=MM BOND={5}\n".format(sym,
                                                                                      p0,
                                                                                      p1,
                                                                                      p2,
                                                                                      mol.mmtypes[atom.index],
                                                                                      ":".join(["{0}/{1}".format(*b) for b in zip(bondIndices+1,bonds)])))
            of.write("END")




#def xyz_to_adf(fxyz):
#    from subprocess import call
#    import platform
#    for atoms in iread(fxyz):
#        if "name" in atoms.info:
#            name = atoms.info["name"]
#        else:
#            name = fxyz.split('.')[0]
#        mol = SBU(name = name, atoms = atoms)
#
#        os.environ["ADFBIN"] = os.environ["ADFHOME"]+'/bin'
#        os.environ["SCMLICENSE"]= os.environ["ADFHOME"]+'/license.txt'
#        adfprep = os.path.abspath(os.path.join(os.environ["ADFBIN"], "adfprep"))
#
#        if os.name == 'nt':
#            adfprep = "sh " + adfprep
#
#
#        args = [adfprep, '-t', 'UFF-GO']
#
#        bonds = numpy.argwhere(numpy.tril(mol.bonds)!=0)
#        args += ['-guibonds']
#        args += ['"{0}"'.format(' '.join(['{0} {1} {2}'.format(i+1, j+1, mol.bonds[i,j]) for i, j in bonds]))]
#        args += ['-m', fxyz, '-a', '"{0}.adf"'.format(name)]
#        args += ['-ufftypes', '"{0}"'.format(' '.join(mol.mmtypes))]
#        #generate the .adf file with correct bonds and regions
#        call(' '.join(args), shell=True)
#
#
##-----------------------------------------------
#
## directory to adfhome
os.environ["ADFHOME"] = '/Applications/ADF2019.103.app/Contents/Resources/adfhome'


xyz_to_inp(qcin)
#xyz_to_adf(qcin)

