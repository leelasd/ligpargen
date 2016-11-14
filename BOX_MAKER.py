import argparse
import numpy as np
from LigParGenTools import *


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='Converter.py',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
	SCRIPT TO CREATE CUSTOM SOLVENT BOXES FOR 
	OPENMM AND NAMD FROM LIGPARGEN FILES
	Created on Mon Nov 14 15:10:05 2016
	@author: Leela S. Dodda leela.dodda@yale.edu
	@author: William L. Jorgensen Lab 
	
	if using PDB file 
	Usage: python BOX_MAKER.py -p UNK_647FFF.pdb -b 40
	
	REQUIREMENTS:
	Preferably Anaconda python with following modules
	argparse
	numpy
	"""
    )
    parser.add_argument(
        "-p", "--pdb", help="Submit PDB file from CHEMSPIDER or PubChem", type=str)
    parser.add_argument("-b", "--box_size", type=float,
                        help="SIZE of the CUBIC box in ANGSTROM")
    args = parser.parse_args()
    if (args.pdb is not None) and (args.box_size is not None):
        BOX_MAKER(args.pdb, BOX_SIZE=args.box_size)
    else:
        print('For Help: python BOX_MAKER.py -h')
