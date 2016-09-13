'''
Created on Aug 6, 2014

@author: oliwa

Project the calpha coordinates of a PDB file by a projection matrix and output the new PDB file with only calpahs
 
'''

import argparse
import sys
import os
import numpy as np
from prody import *
from prody.proteins.pdbfile import parsePDB, writePDB
from prody.measure.measure import calcCenter

def main():
    
    parser = argparse.ArgumentParser(description='Project the calpha coordinates of a PDB file by a projection matrix and output the new PDB file with only calpahs')
    parser.add_argument('pdbFile', help='The pdb file')
    parser.add_argument('P', help='the projection matrix')
    parser.add_argument('-outputName', help='name of the output file, default is beforeP_P_pdbFile and afterP_P_pdbFile')
     
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    
    assert os.path.isfile(args.pdbFile)
    assert os.path.isfile(args.P)  
    
    if args.outputName:
        outputName = args.outputName
    else:
        outputName = os.path.basename(args.P)+"_"+os.path.basename(args.pdbFile)
        
    # read files
    pdbFile = parsePDB(args.pdbFile)
    pdbFile_ca = pdbFile.select('calpha')
    
    P = np.loadtxt(args.P)

    writePDB("beforeP_"+outputName, pdbFile_ca)
    coord_shape = pdbFile_ca.getCoords().shape
    coords_P = P.dot(pdbFile_ca.getCoords().flatten())
    coords_P = coords_P.reshape(coord_shape)
    pdbFile_ca.setCoords(coords_P)
    writePDB("afterP_"+outputName, pdbFile_ca) 

    referenceSegment = "R"
    referenceOfComplex = pdbFile_ca.select('segment \"'+referenceSegment+'.\"')
    print "Center of fixed (receptor) frame is: ", calcCenter(referenceOfComplex.select('calpha'))
    
    print "calphas before projection written to: ", "beforeP_"+outputName
    print "calphas after projection written to: ", "afterP_"+outputName
    
if __name__ == '__main__':
    main()