'''
Created on Nov 8, 2013

@author: oliwa
'''

from pylab import *
from munkres import Munkres, print_matrix, make_cost_matrix

class Hungarian(object):
    ''' Provides the hungarian algorithm to be applied on a profit matrix with
    regards to protein chain matching.
    
    This class imports the Hungarian algorithm from the Munkres package. 
    Given all chains of two proteins (or two protein conformations), and the
    matching results provided by prody, it calculates a profit matrix 
    (matching between chains scored by the profit function of this class). 
    Then, it casts the optimal matching as the maximum weight perfect matching
    problem of a bisection of a graph and solves it.
    '''


    def __init__(self):
        '''
        Constructor.
        '''
        
    def getHungarianIndices(self, 
                            referenceChainsCount, 
                            mobileChainsCount, 
                            matches):
        """
        Run the hungarian algorithm from the Munkres module to determine 
        the optimal matches of chains.
            
            Args:
                referenceChainsCount: number of chains in reference
                mobileChainsCount: number of chains in mobile
                matches: the chain matches as determined by prody
            Returns:
                optimal matches based on the hungarian algorithm in the form
                of indices and the corresponding matchesMatrix for the indicies
        """
        profitStack = [None] * (referenceChainsCount*mobileChainsCount)
        matchesStack = [None] * (referenceChainsCount*mobileChainsCount)
        for element in range(0, len(matches)):
            profitStack[element] = self.hungarianProfit(
                                                    matches[element][2], 
                                                    matches[element][3])
            matchesStack[element] = matches[element]
        profitMatrix = np.zeros((referenceChainsCount,mobileChainsCount))
        matchesMatrix = np.zeros((referenceChainsCount,mobileChainsCount), dtype=object)
        for row in range(0, referenceChainsCount):
            for column in range(0, mobileChainsCount):
                profitMatrix[row][column] = profitStack.pop(0)
                matchesMatrix[row][column] = matchesStack.pop(0)
        profitMatrix = profitMatrix.tolist()
        cost_matrix = make_cost_matrix(
                                       profitMatrix, 
                                       lambda cost: 1000.0 - cost)
        m = Munkres()
        indices = m.compute(cost_matrix)
        return indices, matchesMatrix
            
    def hungarianProfit(self, seqid, overlap):
        """ Denotes the entry of a matching between chains in the profit matix 
        of the hungarian algorithm. This is the profit of an edge. The profit 
        is to be maximized.
        """
        return seqid + overlap