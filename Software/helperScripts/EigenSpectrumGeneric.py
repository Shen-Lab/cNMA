'''
Created on Jul 9, 2014

@author: oliwa
'''

import sys
import glob
import os
from scriptutils import makeStringEndWith, mkdir_p
import argparse
import numpy as np
import traceback
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(description='Visualize eigenvalues and overlaps')
    parser.add_argument('resultsPath', help='Absolute path of the results folders')
    parser.add_argument('outputPath', help='outputPath')
    parser.add_argument('-title', help='title of the plot')
    parser.add_argument('-fileToLookFor', help='Specify the file with the overlap information')
    parser.add_argument('-modes', help='Specify how many modes to plot')
    parser.add_argument('-upperOverlapLimit', help='Upper overlap limit, force manually')
     
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()    
    
    if args.modes:
        modes = int(args.modes)
    else:
        modes = 10
        
    if args.title:
        title = args.title
    else:
        title = ""
        
    fileToLookFor = "singleModeOverlapsFromSuperset.txt"
    if args.fileToLookFor:
        fileToLookFor = args.fileToLookFor
        
    assert os.path.isdir(args.resultsPath)
    assert os.path.isdir(args.outputPath)


    all340proteinsPaths = glob.glob(args.resultsPath+"*/")
    
    dataToPlot = []
    proteins = []
    for proteinPath in sorted(all340proteinsPaths):
        proteinPath = makeStringEndWith(proteinPath, "/")
        try:
            overlap = np.loadtxt(proteinPath+fileToLookFor)
            overlap = overlap[:modes]
            overlap = abs(np.array(overlap))
            if args.upperOverlapLimit:
                for i in range(0, len(overlap)):
                    if overlap[i] > float(args.upperOverlapLimit):
                        overlap[i] = float(args.upperOverlapLimit)
            dataToPlot.append(overlap)
            protein = os.path.basename(os.path.normpath(proteinPath))
            proteins.append(protein)
        except IOError as err:
            print "IOError occurred, probably there is no such file at the path: ", err
            print traceback.format_exc()
            
    dataToPlot = np.array(dataToPlot)      
            
    fig, ax = plt.subplots(1)
    ax.set_yticklabels(proteins, minor=False)
    ax.xaxis.tick_top()    
    
    p = ax.pcolormesh(dataToPlot, cmap="gray")
    fig.colorbar(p)  
    
    # put the major ticks at the middle of each cell, notice "reverse" use of dimension
    ax.set_yticks(np.arange(dataToPlot.shape[0])+0.5, minor=False)
    ax.set_xticks(np.arange(dataToPlot.shape[1])+0.5, minor=False)   
    
    # want a more natural, table-like display (sorting)
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    
    ax.set_xticklabels(range(1, modes+1), minor=False)
    ax.set_yticklabels(proteins, minor=False) 
            
    if args.title:
        plt.title(args.title+"\n\n")
        
    # output
    outputPath = makeStringEndWith(args.outputPath, "/")+"eigenVis"
    mkdir_p(outputPath)
    plt.savefig(outputPath+'/eigenVis_'+title+'.eps', bbox_inches='tight')
    plt.savefig(outputPath+'/eigenVis_'+title+'.pdf', bbox_inches='tight') 
    #plt.show()
    # close and reset the plot 
    plt.clf()
    plt.cla()
    plt.close()         
            
if __name__ == '__main__':
    main()