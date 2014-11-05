#!/usr/bin/env python
__author__ = "Katra Kolsek and Samo Turk"
__copyright__ = "Copyright (C) 2014 by Katra Kolsek"
__license__ = "GPL"
__version__ = "0.3"
__status__ = "Development"

import pybel
from subprocess import call
from os import listdir

def anymol2vina(mols, conf, rec, gen3d):
    """
    converts mols to pdbq and docks with vina
    """
    for mol in mols:
        if str(mol.title)+".pdbqt" not in listdir('.'):
            if gen3d == True:
                mol.make3D(forcefield="MMFF94", steps=1000)
            mol.write("pdbqt", str(mol.title)+".pdbqt", overwrite=True)
            try:
                call(["vina", "--config", conf, "--receptor", rec, "--ligand", str(mol.title)+".pdbqt"])
            except:
                print "Error with: " + str(mol.title)
        else:
            mol.addh()
            mol1 = mol.write("inchi").rstrip()
            molread = pybel.readfile("pdbqt", filename=str(mol.title)+".pdbqt").next()
            molread.addh()
            mol2 = molread.write("inchi").rstrip()
            if mol1 != mol2:
                print str(mol.title) + " duplicated name, different mol, still docking, saving with _alt_ suffix"
                if gen3d == True:
                    mol.make3D(forcefield="MMFF94", steps=1000)
                mol.title = str(mol.title) + "_alt"
                mol.write("pdbqt", str(mol.title)+".pdbqt", overwrite=True)
                try:
                    call(["vina", "--config", conf, "--receptor", rec, "--ligand", str(mol.title)+".pdbqt"])
                except:
                    print "Error with: " + str(mol.title)
            else:
                print str(mol.title) + " duplicated.."

if __name__ == "__main__":
    """
    smi2vina smiles config_confeptor
    """
    from sys import argv
    if len(argv) == 5:
        molsin = argv[1]  #determine filetype
        ftypein = molsin.split(".")[-1]
        mols = pybel.readfile("%s" % ftypein, filename="%s" % molsin) #openfile
        conf = argv[2]
        rec = argv[3]
        if argv[4].lower() == "gen3d" or ftypein == "smi" or ftypein == "ism":
            print "3d coordinates will be generated."
            gen3d = True
        else:
            gen3d = False
        anymol2vina(mols, conf, rec, gen3d)
    elif len(argv) == 4:
        molsin = argv[1]  #determine filetype
        ftypein = molsin.split(".")[-1]
        mols = pybel.readfile("%s" % ftypein, filename="%s" % molsin) #openfile
        conf = argv[2]
        rec = argv[3]
        if ftypein == "smi" or ftypein == "ism":
            print "Mols are in smi or ism, generating 3d coordinates."
            gen3d = True
        else:
            gen3d = False
        anymol2vina(mols, conf, rec, gen3d)
    else:
        print "Usage: smi2vina mols config_receptor receptor.pdbqt gen3d."
        print "You can ommit gen3d if your mol file already has 3d coordinates."
        print "If mols are in smi or ism, 3d coordinates will allways be generated."
        print "File type is determined by extension."
