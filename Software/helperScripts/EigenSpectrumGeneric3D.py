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
#import pylab
import matplotlib
matplotlib.use('Agg')

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from scriptutils import makeStringNotEndWith

def main():
    parser = argparse.ArgumentParser(description='Visualize eigenvalues and overlaps')
    parser.add_argument('resultsPath', help='Absolute path of the results folders')
    parser.add_argument('outputPath', help='outputPath')
    parser.add_argument('-title', help='title of the plot')
    parser.add_argument('-fileToLookFor_overlap', help='Specify the file with the overlap information')
    parser.add_argument('-fileToLookFor_differencesInRank', help='Specify the file with the differencesInRank information')    
    parser.add_argument('-modes', help='Specify how many modes to plot')
    parser.add_argument('-upperOverlapLimit', help='Upper overlap limit, force manually')
     
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()    
    
    if args.modes:
        modes = int(args.modes)
    else:
        modes = 4
        
    if args.title:
        title = args.title
    else:
        title = ""
        
    if args.outputPath:
        outputPath = args.outputPath
    else:
        outputPath = ""        
        
    fileToLookFor_overlap = "singleModeOverlapsFromSuperset.txt"
    fileToLookFor_differencesInRank = "differencesInRank.txt"
    if args.fileToLookFor_overlap:
        fileToLookFor_overlap = args.fileToLookFor
    if args.fileToLookFor_differencesInRank:
        fileToLookFor_differencesInRank = args.fileToLookFor_differencesInRank        
        
    assert os.path.isdir(args.resultsPath)
    assert os.path.isdir(args.outputPath)


    all340proteinsPaths = glob.glob(args.resultsPath+"*/")
    
    difficults = np.loadtxt("/home/oliwa/workspace/TNMA1/src/BenchmarkAssessmentsOfDifficulty/allinterfaceSuperposed/difficult.txt", dtype="string")
    difficults = set(difficults)    
    
    dataToPlot_overlaps = []
    dataToPlot_differencesInRank = []
    proteins = []
    counter = 0
    
    for proteinPath in sorted(all340proteinsPaths):
        proteinPath = makeStringEndWith(proteinPath, "/")
        
        protein = makeStringNotEndWith(os.path.basename(os.path.normpath(proteinPath)), "/")
        
        if protein not in difficults:
            continue
        counter += 1

        try:
            # load overlap
            overlap = np.loadtxt(proteinPath+fileToLookFor_overlap)
            overlap = overlap[:modes]
            overlap = abs(np.array(overlap))
            overlap = list(overlap)
            if args.upperOverlapLimit:
                for i in range(0, len(overlap)):
                    if overlap[i] > float(args.upperOverlapLimit):
                        overlap[i] = float(args.upperOverlapLimit)
            dataToPlot_overlaps.append(overlap)
            protein = os.path.basename(os.path.normpath(proteinPath))
            proteins.append(protein)
            # load ranking differences
            differenceInRank = np.loadtxt(proteinPath+fileToLookFor_differencesInRank, dtype="int")
            differenceInRank = list(differenceInRank)
            dataToPlot_differencesInRank.append(differenceInRank[:modes])

        except IOError as err:
            print "IOError occurred, probably there is no such file at the path: ", err
            print traceback.format_exc()           
            
            
    print proteins
            
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #x, y = np.random.rand(2, 100) * 4
    y = range(1, len(proteins)+1)
    x = range(1, modes+1)
    
    xpos, ypos = np.meshgrid(x, y)
    x = xpos.flatten()
    y = ypos.flatten()
    
    colors = []
    print "overlaps len: ", len(dataToPlot_overlaps)
    print "overlaps: ", dataToPlot_overlaps    
    
    dataToPlot_overlaps_flattened = np.array(dataToPlot_overlaps).flatten()
    maxOverlap = max(dataToPlot_overlaps_flattened)
    print "maxOverlap:", maxOverlap
    
    for element in dataToPlot_overlaps_flattened:
        colors.append(plt.cm.jet(element/maxOverlap))
        #print plt.cm.jet(element/maxOverlap)

    print "x", len(x)
    print "y", len(y)
    #print "colors", len(colors)
    
    print "dataToPlot_differencesInRank len: ",dataToPlot_differencesInRank
    dataToPlot_differencesInRank = np.array(dataToPlot_differencesInRank).flatten() + 0.0001
    print "dataToPlot_differencesInRank len: ", len(dataToPlot_differencesInRank.flatten())
    
    dx=np.ones(len(x))*0.5
    dy=dx
    
    p = ax.bar3d(x-0.25, y-0.25, np.zeros(len(x)), dx, dy, dataToPlot_differencesInRank, color=colors, zsort='average')
    
    ax.set_zlim([min(dataToPlot_differencesInRank), max(dataToPlot_differencesInRank)])
    #ax.set_title(title)
    
    # x label for the ascending modes
    #ax.set_xticklabels(range(1, modes+1), minor=False)
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    ax.set_xlabel("ascending lambda^R modes") 
    # y label for the proteins
    #ax.set_yticklabels(proteins, minor=False)
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    ax.set_ylabel("proteins")
# #     dataToPlot_overlaps = np.array(dataToPlot_overlaps)      
# #             
# #     fig, ax = plt.subplots(1)
# #     ax.set_yticklabels(proteins, minor=False)
# #     ax.xaxis.tick_top()    
# #     
# #     p = ax.pcolormesh(dataToPlot_overlaps, cmap="bone")
# #     fig.colorbar(p)  
# #     
# #     # put the major ticks at the middle of each cell, notice "reverse" use of dimension
# #     ax.set_yticks(np.arange(dataToPlot_overlaps.shape[0])+0.5, minor=False)
# #     ax.set_xticks(np.arange(dataToPlot_overlaps.shape[1])+0.5, minor=False)   
# #     
# #     # want a more natural, table-like display (sorting)
# #     ax.invert_yaxis()
# #     ax.xaxis.tick_top()
# #     
# #     ax.set_xticklabels(range(1, modes+1), minor=False)
# #     ax.set_yticklabels(proteins, minor=False) 
# #             
# #     if args.title:
# #         plt.title(args.title+"\n\n")
        
    # output
    #outputPath = makeStringEndWith(args.outputPath, "/")+"eigenVis"
    
    #mkdir_p(outputPath)
    plt.savefig(outputPath+'/eigenVis_'+title+'.eps', bbox_inches='tight')
    plt.savefig(outputPath+'/eigenVis_'+title+'.pdf', bbox_inches='tight') 
    #plt.show()
    # close and reset the plot 
    plt.clf()
    plt.cla()
    plt.close()         
    print "total proteins: ", counter
            
if __name__ == '__main__':
    main()