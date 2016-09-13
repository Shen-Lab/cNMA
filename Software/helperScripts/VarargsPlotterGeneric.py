'''
Created on June 29, 2014

@author: oliwa
'''

import argparse
import sys
from helperScripts.scriptutils import startModeIndexFromOne, mkdir_p
import numpy as np
import matplotlib
from numpy import dtype
matplotlib.use('Agg')
import matplotlib.pyplot as plt

class VarargsPlotterGeneric(object):
    '''
    Plots a variable number of arguments, the input is expected to given in multiples of three, always in the format: 
    y-axis, x-axis, plotname, for a variable number of experimental results
    '''

    def __init__(self):
        '''
        Constructor of SuccessFailurePlotter
        '''
        parser = argparse.ArgumentParser(description='Creates the success-failure plots based on output of the SuccessFailureAnalyzer')
        parser.add_argument("inputs", help="File input for the VarargsPlotter. The program expects inputs in multiples of three (3,6,9,...), corresponding to a files with the y-axis values, the x-axis values and a string with a description of the plot", nargs="+")
        parser.add_argument('-title', help='Title of the plot')
        parser.add_argument('-xLabel', help='xLabel of the plot')
        parser.add_argument('-yLabel', help='yLabel of the plot')
        parser.add_argument('-ticks', help='number of ticks (modes) to plot')
        parser.add_argument('-outputFile', help='output path and filename of the plots (provide no extensions), default is ./plotOutput.pdf|eps')
        if len(sys.argv)==1:
            parser.print_help()
            sys.exit(1)
        args = parser.parse_args()
        self.inputs = args.inputs
        assert len(self.inputs) % 2 == 0
        if args.title:
            self.title = args.title
        else:
            self.title = ""
        if args.xLabel:
            self.xLabel = args.xLabel
        else:
            self.xLabel = ""
        if args.yLabel:
            self.yLabel = args.yLabel
        else:
            self.yLabel = ""        
        if args.outputFile:
            self.outputFile = args.outputFile
        else:
            self.outputFile = "outputPlot" 
        if args.ticks:
            self.ticks = args.ticks
        else:
            self.ticks = None

    def mainPlot(self):
        # make the plot
        fontSize = 19
        lineWidth = 8
        legendSize = 10
        
        plt.title(self.title, size=fontSize)
        plt.tick_params(axis='both', which='major', labelsize=fontSize)
        plt.tick_params(axis='both', which='minor', labelsize=fontSize)
        plt.xlabel(self.xLabel, size=fontSize)
        plt.ylabel(self.yLabel, size=fontSize)      
        
        while self.inputs:
            yAxis = np.loadtxt(self.inputs.pop(0))
            description = self.inputs.pop(0)
            if description.startswith("lambdaC"):
                yAxis = yAxis[6:]
                
            if self.ticks is None:
                ticks = range(1,len(yAxis)+1)
            else:
                ticks = range(1,int(self.ticks)+1)
                yAxis = yAxis[:int(self.ticks)]
            plt.plot(ticks, yAxis, linewidth=lineWidth-2, label=description)

        plt.legend(loc='best', prop={'size':legendSize})
        #plt.show()
        plt.savefig(self.outputFile+'.eps', bbox_inches='tight')
        plt.savefig(self.outputFile+'.pdf', bbox_inches='tight')
        # close and reset the plot
        plt.clf()
        plt.cla()
        plt.close()
        
        
if __name__ == '__main__':
    varargsPlotter = VarargsPlotterGeneric()
    varargsPlotter.mainPlot()