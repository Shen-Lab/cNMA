'''
Created on Aug 6, 2014

@author: oliwa

Project the normal modes of a NMD file by a projection matrix and output the new NMD file
 
'''

import argparse
import sys
import os
import numpy as np

def projectNMDLine(line, P, normalize):
    """ Given a mode line from a NMD file and a projection matrix P, return the projected NMD mode line 
     
    Args:
        line: NMD mode line (with k coordinate elements)
        P: projection matrix (with k*k elements)
        
    Returns: projected NMD mode line by P
    """
    # find where the mode starts
    firstSpace = line.find(" ")
    secondSpace = line[firstSpace+1:].find(" ")
    thirdSpace = line[firstSpace+1+secondSpace+1:].find(" ")
    
    # extract mode
    mode = line[firstSpace+1+secondSpace+1+thirdSpace+1:].split()
    mode = np.array(mode, dtype='float')
    
    # test mode for center, 15:47 Aug 7th
    #modeM = mode.reshape(len(mode)/3, 3)
    #modeM = modeM[:736]
    #print "delta x ", modeM.T[:1].sum()
    #print "delta y ", modeM.T[1:2].sum()
    #print "delta z ", modeM.T[2:3].sum()
    
    # assert dimensionalities
    assert(len(mode) == P.shape[0])
    assert(len(mode) == P.shape[1])
    
    # multiply P times the mode
    p_mode = P.dot(mode)
    if normalize:
        p_mode =  p_mode/np.linalg.norm(p_mode)
    
    # recombine it to a string
    p_mode = ' '.join(str("{0:.3f}".format(e)) for e in p_mode.round(3))
    p_mode = line[:firstSpace+1+secondSpace+1+thirdSpace+1]+p_mode+"\n"
    
    return p_mode

def main():
    
    parser = argparse.ArgumentParser(description='Project the normal modes of a NMD file by a projection matrix and output the new NMD file')
    parser.add_argument('nmdFile', help='The nmd file')
    parser.add_argument('P', help='the projection matrix')
    parser.add_argument('-outputName', help='name of the output file, default is P_nmdFile')
    parser.add_argument('-normalize', action="store_true", help='normalize each mode after projection')
     
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    
    assert os.path.isfile(args.nmdFile)
    assert os.path.isfile(args.P)  
    
    if args.outputName:
        outputName = args.outputName
    else:
        outputName = os.path.basename(args.P)+"_"+os.path.basename(args.nmdFile)
        
    # read files
    nmdFile = open(args.nmdFile, "r")
    P = np.loadtxt(args.P)
    P_nmdFile = []
    
    lines = nmdFile.readlines()

    for line in lines:
        # filename
        if line.startswith("name"):
            P_nmdFile.append(line[:5]+outputName+"\n")

        if not line.startswith("name") and not line.startswith("mode"):
            P_nmdFile.append(line)

        if line.startswith("mode"):
            P_line = projectNMDLine(line, P, args.normalize)
            P_nmdFile.append(P_line)
            
    nmdFile.close()

    # write output
    f = open(outputName,'w')
    for item in P_nmdFile:
        f.write("%s" % item)
    f.close()
    print "output written to: ", outputName
    
if __name__ == '__main__':
    main()