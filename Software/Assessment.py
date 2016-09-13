'''
Created on Nov 20, 2013

@author: oliwa
'''

#from pylab import *
import numpy as np, numpy
import matplotlib.pyplot as plt
from numpy.core.numerictypes import float64
from datetime import datetime
from prody.dynamics.functions import writeArray
import math
from prody.dynamics.compare import calcOverlap
from scipy.stats import pearsonr
from prody.dynamics.analysis import calcCollectivity
import sys
import itertools
from prody.dynamics.mode import Vector, Mode

class Assessment(object):
    ''' Assessment class with measures to quantify NMA.'''


    def __init__(self, utils, numberOfReferenceProteins):
        '''Constructor'''
        self.utils = utils
        self.experimentName = "experimentName"
        self.numberOfReferenceProteins = numberOfReferenceProteins
        self.sumHighestOverlap = 0.0
        self.sumHighestOverlapisWithin5LowestModes = 0
#         self.arrayCountBestMode = np.zeros(numberOfModes, dtype=int)
#         self.arrayCountHighestCollectivity = np.zeros(numberOfModes, dtype=int)
        self.TapproxOverlapSum = 0.0
#         self.arraySumofAbsOverlaps = np.zeros(numberOfModes, dtype=float64)
#         self.arraySumofAbsCollectivities = np.zeros(numberOfModes, dtype=float64)
#         self.arrayCumulSumofAbsOverlap = np.zeros(numberOfModes, dtype=float64)
#         self.arrayCumulSumofAbsCollectivities = np.zeros(numberOfModes, dtype=float64)
        self.listOfHighestOverlaps = []
        self.listOfTApproxOverlaps = []
        self.listofRMSDbefore = []
        self.listofRMSDInterfacebefore = []
        self.listofRMSDafterTapprox = []
        self.listofRMSDReductionsWhole = []
        self.listofRMSDReductionsInterface = []
        self.listofoverlapTApproxWhole = []
        self.listofoverlapTApproxInterface = []
        self.listofArraysofModeOverlaps = []
        self.listofArraysofModeCollectivities = []
        self.listofArraysofModeCorrelations = []
        self.listofArraysofRMSDReductionsWholePerc = []
        self.listofArraysofRMSDReductionsInterfacePerc = []
        self.listofArraysofModeOverlapsInterface = []
        self.listofArraysofModeCollectivitiesInterface = []
        self.listofArraysofModeCorrelationsInterface = []
        self.listofArraysofEigenvalues = []
        self.countHighestOverlapHighestCorrelationWhole = 0
        self.countHighestOverlapHighestCorrelationInterface = 0
        self.listOfNumberofModes = []
        ## for y
        self.listofHighestOverlapsPositions = []
        self.listOfHighestOverlapsInterface = []
        self.listofHighestOverlapsInterfacePositions = []
        self.listofWorkedOnTuples = []
        
    def calcDiscretizedCount01(self, myList):
        """ Discretize (in 0.1 steps) the real-valued counts from myList. 
        
            Args:
                myList: real-valued list
                
            Returns: arraydiscretizedCount: the discretized list
        """
        arraydiscretizedCount = np.zeros(10, dtype=int)
        for item in myList:
            assert(item >= -0.00000001 and item <= 1.00000001)
            if item < 0.1:
                arraydiscretizedCount[0] += 1
            elif item < 0.2:
                arraydiscretizedCount[1] += 1
            elif item < 0.3:
                arraydiscretizedCount[2] += 1
            elif item < 0.4:
                arraydiscretizedCount[3] += 1
            elif item < 0.5:
                arraydiscretizedCount[4] += 1
            elif item < 0.6:
                arraydiscretizedCount[5] += 1
            elif item < 0.7:
                arraydiscretizedCount[6] += 1
            elif item < 0.8:
                arraydiscretizedCount[7] += 1
            elif item < 0.9:
                arraydiscretizedCount[8] += 1
            else:
                arraydiscretizedCount[9] += 1
        return arraydiscretizedCount
                
    def calcDiscretizedCount005(self, myList):
        """ Discretize (in 0.05 steps) the real-valued counts from myList. 
        
            Args:
                myList: real-valued list
                
            Returns: arraydiscretizedCount: the discretized list
        """
        arraydiscretizedCount = np.zeros(20, dtype=int)
        for item in myList:
            assert(item >= -0.00000001 and item <= 1.00000001)
            if item < 0.05:
                arraydiscretizedCount[0] += 1
            elif item < 0.1:
                arraydiscretizedCount[1] += 1
            elif item < 0.15:
                arraydiscretizedCount[2] += 1
            elif item < 0.2:
                arraydiscretizedCount[3] += 1
            elif item < 0.25:
                arraydiscretizedCount[4] += 1
            elif item < 0.3:
                arraydiscretizedCount[5] += 1
            elif item < 0.35:
                arraydiscretizedCount[6] += 1
            elif item < 0.4:
                arraydiscretizedCount[7] += 1
            elif item < 0.45:
                arraydiscretizedCount[8] += 1
            elif item < 0.5:
                arraydiscretizedCount[9] += 1
            elif item < 0.55:
                arraydiscretizedCount[10] += 1
            elif item < 0.6:
                arraydiscretizedCount[11] += 1
            elif item < 0.65:
                arraydiscretizedCount[12] += 1
            elif item < 0.7:
                arraydiscretizedCount[13] += 1
            elif item < 0.75:
                arraydiscretizedCount[14] += 1
            elif item < 0.8:
                arraydiscretizedCount[15] += 1
            elif item < 0.85:
                arraydiscretizedCount[16] += 1
            elif item < 0.9:
                arraydiscretizedCount[17] += 1
            elif item < 0.95:
                arraydiscretizedCount[18] += 1
            else:
                arraydiscretizedCount[19] += 1
        return arraydiscretizedCount
    
#     def calcOverlap(self, ref_chain, mob_chain, mode):
#         """ Get overlap, which provides a measure of how well the displacement 
#         of atoms can be explained by the normal mode.
#         
#         This measure is defined in:
#         
#         Moal, Iain H., and Paul A. Bates. "SwarmDock and the Use of Normal 
#         Modes in Protein-Protein Docking." International Journal of 
#         Molecular Sciences 11, no. 10 (September 28, 2010): 3623-3648. 
#         doi:10.3390/ijms11103623. 
#         
#         Args:
#             ref_chain: unbound structure
#             rob_chain: bound structure
#             mode: normal mode
#             
#         Returns:
#             A measure of overlap between mode and displacement
#             
#         """
#         
#         if (ref_chain.numAtoms() != mob_chain.numAtoms()) or (ref_chain.numAtoms() != mode.numAtoms()):
#             raise ValueError("Cannot calculate overlap, the number of atoms in ref_chain, mob_chain and mode do not match")
#         
#         # get iterators, writePDB workaround to use the iterators from selections
#         ruIter = ref_chain.iterAtoms()
#         rbIter = mob_chain.iterAtoms()
#         
#         # calculate the nominator
#         nominator = 0.0
#         for i in range (0, mode.numAtoms()):
#             ru = ruIter.next().getCoords()
#             rb = rbIter.next().getCoords()
#             ai = mode.getArrayNx3()[i]
#             nominator += numpy.sum(ai * (rb-ru))
#         nominator = numpy.abs(nominator)
#             
#         # reset the iterators to the first atom
#         ruIter = ref_chain.iterAtoms()
#         rbIter = mob_chain.iterAtoms()
#         
#         # calculate the denominator
#         rbruSquaredSum = 0.0
#         aiSum = numpy.sum(numpy.power(mode.getArray(), 2))
#         for i in range (0, mode.numAtoms()):
#             ru = ruIter.next().getCoords()
#             rb = rbIter.next().getCoords()
#             rbruSquaredSum += numpy.sum(numpy.power(rb-ru, 2))
#         denominator = numpy.sqrt(aiSum * rbruSquaredSum)
#         
#         # calculate the overlap
#         try:
#             overlap = nominator/denominator
#             if math.isnan(overlap):
#                 overlap = 0.0
#         except ZeroDivisionError as e:
#             print "denominator is zero, divide by zero exception" + str(e)
#         return overlap
    
    def calcOverlapArray(self, ref_chain, mob_chain, anm):
        """ Calculate the overlap array by calling calcOverlap(...) on 
        every mode. 
        
        Args:
            ref_chain: unbound structure
            rob_chain: bound structure
            anm: the ANM with normal modes
            
        Returns: an array with the overlaps for each mode
        
        """
        if (ref_chain.numAtoms() != mob_chain.numAtoms()) or (ref_chain.numAtoms() != anm.numAtoms()):
            raise ValueError("Cannot calculate overlapArray, the number of atoms in ref_chain, mob_chain and anm do not match")
        
        overlapArray = numpy.array([], dtype=float64)
        for mode in anm:
            overlapArray = np.append(overlapArray, self.calcOverlap(ref_chain, mob_chain, mode))
        return overlapArray
    
    def calcOverlapArrayWithProDyOverlap(self, anm, defvec):
        """ Calculate the overlap array by calling ProDys calcOverlap(...) on 
        every mode. 
        
        Args:
            anm: the ANM with normal modes
            defvec: the true deformation vector
            
        Returns: an array with the overlaps for each mode
        
        """
        if anm.numAtoms() != defvec.numAtoms():
            raise ValueError("Cannot calculate overlapArray, the number of atoms in anm and the defvec do not match")
        
        overlapArray = numpy.array([], dtype=float64)
        for mode in anm:
            overlapArray = np.append(overlapArray, calcOverlap(anm[mode], defvec))
        return overlapArray
    
    def calcCorrelationArray(self, anm, defvec):
        """ Calculate the pearson correlation of the magnitudes of displacement for all modes
            against the deformation vector. 
            
            Args: 
                anm: the ANM with normal modes
                defvec: the true deformation vector
            
            Returns: an array with the correlations for each mode      
            """
        if anm.numAtoms() != defvec.numAtoms():
            raise ValueError("Cannot calculate correlationArray, the number of atoms in anm and the defvec do not match")
        
        correlationArray = []
        defvecMagnitudes = self.calcMagnitudeArray(defvec)
        
        for mode in anm:
            modeMagnitudes = self.calcMagnitudeArray(anm[mode])
            correlationArray.append(self.calcPearson2(modeMagnitudes, defvecMagnitudes))
        correlationArray = numpy.array(correlationArray, dtype=float64)
        return correlationArray
    
    def calcCumulArray(self, overlap, array=False):
        """This method is taken straight from prody and modified to work 
        directly on an array, so that the output of for instance 
        calcOverlapArray(...) can be given as input  here.
    
        For example, in case of given an overlap array this returns cumulative 
        overlap of the overlap array given in the argument.
        
        Args:
            overlap: the overlap array
            array: should it return the whole array (array == False) or the final end value (array == True)
        Return: The cumulative overlap (either array or final end value)
        """
    
        if array:
            return np.sqrt(np.power(overlap, 2).sum(axis=overlap.ndim-1))
        else:
            return np.sqrt(np.power(overlap, 2).cumsum(axis=overlap.ndim-1))
        
#     def calcCollectivity(self, mode):
#         """ Calculate the collectivity of a mode, this scalar measures the 
#         number of highly mobile atoms. 
#         
#         The calculation of collectivity from mode is implemented as defined in:
#         Dobbins, Sara E., Victor I. Lesk, and Michael J. E. Sternberg. 
#         "Insights into Protein Flexibility: The Relationship between Normal 
#         Modes and Conformational Change upon Protein-protein Docking." 
#         Proceedings of the National Academy of Sciences 105, no. 30 (July 29, 
#         2008): 10390-10395. doi:10.1073/pnas.0802496105.
# 
#         Args:
#             mode: the mode that the collectivity should be calculated on
#         Return: the collectivity of the mode
#         """
#         modeArray = mode.getArray()
#         modeSquared = numpy.power(modeArray, 2)
#         sumDisplacement = (-1) * numpy.sum(modeSquared * numpy.log(modeSquared))
#         collectivity = (1.0 / mode.numAtoms()) * numpy.exp(sumDisplacement)
#         return collectivity

    def calcProdyCollectivity2(self, mode, masses=None):
        """Return collectivity of the mode.  This function implements collectivity
        as defined in equation 5 of [BR95]_.  If *masses* are provided, they will
        be incorporated in the calculation.  Otherwise, atoms are assumed to have
        uniform masses.
    
        .. [BR95] Bruschweiler R. Collective protein dynamics and nuclear
           spin relaxation. *J Chem Phys* **1995** 102:3396-3403.
    
        :arg mode: mode or vector
        :type mode: :class:`.Mode` or :class:`.Vector`
    
        :arg masses: atomic masses
        :type masses: :class:`numpy.ndarray`"""
    
        if not isinstance(mode, (Mode, Vector)):
            raise TypeError('mode must be a Mode or Vector instance')
    
        is3d = mode.is3d()
        if masses is not None:
            if len(masses) != mode.numAtoms():
                raise ValueError('length of massesmust be equal to number of atoms')
            if is3d:
                u2in = (mode.getArrayNx3() ** 2).sum(1) / masses
        else:
            if is3d:
                u2in = (mode.getArrayNx3() ** 2).sum(1)
            else:
                u2in = (mode.getArrayNx3() ** 2)
        u2in = u2in * (1 / u2in.sum() ** 0.5)
        coll = np.exp(-(u2in * np.log(u2in)).sum()) / mode.numAtoms()
        return coll
    
    def calcProDyCollectivityVector(self, vector):
        is3d = vector.is3d()
        if is3d:
            u2in = (vector.getArrayNx3() ** 2).sum(1)
        else:
            u2in = (vector.getArrayNx3() ** 2)
        u2in = u2in * (1 / u2in.sum() ** 0.5)
        coll = np.exp(-(u2in * np.log(u2in)).sum()) / vector.numAtoms()
        return coll
    
    def calcCollectivityArray(self, anm):
        """ Calculate the collectivity array by calling self.calcProdyCollectivity2(...) on 
        every mode. Every mode is normalized prior to the collectivity calculation.
        
        Args:
            anm: the ANM with normal modes
            
        Returns: the array with the collectivities of each mode
        
        """
   
        collectivityArray = numpy.array([], dtype=float64)
        for mode in anm:
            modeArray = mode.getArray()
            modeVector = Vector(modeArray)
            modeNormalized = modeVector.getNormed()
            collectivityArray = np.append(collectivityArray, self.calcProdyCollectivity2(modeNormalized))
        return collectivityArray
    
    def calcCollectivityOfSingleMode(self, mode):
        modeArray = mode.getArray()
        modeVector = Vector(modeArray)
        modeNormalized = modeVector.getNormed()
        return self.calcProdyCollectivity2(modeNormalized)
    
    def calcProcentualReductionArray(self, RMSD_before, RMSD_after):
        """ Calculate the procentual reduction between the values in the arrays 
        RMSD_before and RMSD_after 
        
        Args:
            RMSD_before: an array with RMSD information 
                (before modifying coordinates of a protein 
                via for example a linear combination of normal modes)
            RMSD_after: an array with RMSD information (after modifying the 
            coordinates of a protein)
        
        Returns: the array with the procentual decreases, each between 0 to 1
        
        """
        assert len(RMSD_before) == len(RMSD_after)
        result = np.zeros(len(RMSD_before), dtype=float64)
        for i in range(0, len(RMSD_before)):
            result[i] = self.percentageDecrease(RMSD_before[i], RMSD_after[i])
        return result
    
    def calcProcentualReductionArray2(self, RMSD_before, reductionArray):
        result = []
        for element in reductionArray:
            result.append(self.percentageDecrease(RMSD_before, element))
        return np.array(result)
    
    def collectAssessmentMeasures(self, 
                                  myTuple,
                                  overlapArray, 
                                  collectivityArray, 
                                  TapproxOverlap, 
                                  RMSD_unbound_to_superposed_bound, 
                                  RMSD_interface, 
                                  RMSDReductionsWhole, 
                                  RMSDReductionsInterface,
                                  RMSD_after_Tapprox,
                                  overlapTApproxWhole,
                                  overlapTApproxInterface,
                                  correlationArray,
                                  stepPointsReductionWhole,
                                  stepPointsReductionInterface,
                                  overlapArrayInterface,
                                  collectivityArrayInterface,
                                  correlationArrayInterface,
                                  eigenvalues,
                                  numberOfCurrentModes):
        """ Collect and update various statistics, currently based on the 
        overlap, collectivity and the overlap of an approximated deformation 
        vector through a linear combination of normal modes with the real one.
        
        Args:
            overlapArray: Array with overlap information
            collectivityArray: Array with collectivity information
            TapproxOverlap: Highest absolute overap of a mode on a protein 
            RMSD_unbound_to_superposed_bound: RMSD before applying NMA coordinate transformations
            RMSD_interface: RMSD on the interface before applying NMA coordinate transformations
            RMSDReductionsWhole: RMSD reductions on the whole protein through an increasing number of modes 
            RMSDReductionsInterface: RMSD reductions on the interface through an increasing number of modes
            RMSD_after_Tapprox: RMSD after applying an linear combination of all modes on the whole protein 
        
        """
        # mode overlap
        if(np.abs(overlapArray).max() > 1):
            print "overlap over 1" + str(np.abs(overlapArray).max())
            sys.exit()
        if(np.abs(overlapArray).max() < 0):
            print "overlap smaller 0" + str(np.abs(overlapArray).max())
            sys.exit()
        if math.isnan(np.abs(overlapArray).max()):
            print "overlap nan" + str(np.abs(overlapArray).max())
            sys.exit()
        self.sumHighestOverlap += np.abs(overlapArray).max()
        self.listOfHighestOverlaps.append(np.abs(overlapArray).max())
        highestOverlapPos = np.argmax(np.abs(overlapArray))
        self.listofHighestOverlapsPositions.append(highestOverlapPos)
        if highestOverlapPos < 5:
            self.sumHighestOverlapisWithin5LowestModes += 1
        # best mode count
#         self.arrayCountBestMode[highestOverlapPos] += 1
        # best collectivity count
#         highestCollectivityPos = np.argmax(np.abs(collectivityArray))
#         self.arrayCountHighestCollectivity[highestCollectivityPos] += 1
        # Tapprox
#         self.TapproxOverlapSum += numpy.abs(TapproxOverlap)
#         self.listOfTApproxOverlaps.append(numpy.abs(TapproxOverlap))
        # RMSD
        self.listofRMSDbefore.append(RMSD_unbound_to_superposed_bound)
        self.listofRMSDInterfacebefore.append(RMSD_interface)
        
#         # pad RMSD reductions with last values
#         RMSDReductionsWhole = self.padArrayWithLastValue(RMSDReductionsWhole, maxModes)
#         RMSDReductionsInterface = self.padArrayWithLastValue(RMSDReductionsInterface, maxModes)
        
        self.listofRMSDReductionsWhole.append(RMSDReductionsWhole)
        self.listofRMSDReductionsInterface.append(RMSDReductionsInterface)     
        
        self.listofRMSDafterTapprox.append(RMSD_after_Tapprox)
        # absolute overlaps/collectivities
        ##self.arraySumofAbsOverlaps += np.abs(overlapArray)
        ##self.arraySumofAbsCollectivities += np.abs(collectivityArray)
        # overlap of TApprox and the true defvec
        self.listofoverlapTApproxWhole.append(overlapTApproxWhole)
        self.listofoverlapTApproxInterface.append(overlapTApproxInterface)
        ## collect more measures
        #print "before self.listofArraysofModeOverlaps: " + str(self.listofArraysofModeOverlaps)
        self.listofArraysofModeOverlaps.append(overlapArray)
        #print "after appending overlaparray: "+str(overlapArray)
        #print "self.listofArraysofModeOverlaps: ", self.listofArraysofModeOverlaps
        #sys.exit()
        self.listofArraysofModeOverlapsInterface.append(overlapArrayInterface)
        self.listofArraysofModeCollectivities.append(collectivityArray)
        self.listofArraysofModeCollectivitiesInterface.append(collectivityArrayInterface)
        self.listofArraysofModeCorrelations.append(correlationArray)
        self.listofArraysofModeCorrelationsInterface.append(correlationArrayInterface)
        self.listofArraysofRMSDReductionsWholePerc.append(self.calcProcentualReductionArray2(RMSD_unbound_to_superposed_bound, RMSDReductionsWhole))
        self.listofArraysofRMSDReductionsInterfacePerc.append(self.calcProcentualReductionArray2(RMSD_interface, RMSDReductionsInterface))
        self.listofArraysofEigenvalues.append(eigenvalues)
        #self.countOverlapCorrelationAgreement(overlapArray, correlationArray, overlapArrayInterface, correlationArrayInterface)
        
        # collect the number of modes
        self.listOfNumberofModes.append(numberOfCurrentModes)
        
        # collect highest overlap information
        self.listOfHighestOverlapsInterface.append(np.abs(overlapArrayInterface).max())
        self.listofHighestOverlapsInterfacePositions.append(np.argmax(np.abs(overlapArrayInterface)))
        self.listofWorkedOnTuples.append(myTuple)
        #print "myTuple: ", myTuple
        #sys.exit()
        
    def calcCumulAssessments(self):
        """ This calculates the assessment measures after the NMA investigations 
        for all tuples have been performed. """
        print 'start cumulative assessment '
#         self.arrayCumulSumofAbsOverlap = self.calcCumulArray(self.arraySumofAbsOverlaps)
#         self.arrayCumulSumofAbsCollectivities = self.calcCumulArray(self.arraySumofAbsCollectivities)
#         
#         self.arrayOfHighestOverlaps = np.array(self.listOfHighestOverlaps)
#         self.arrayCountofBestModePercentage = self.arrayCountBestMode.astype(float) / self.arrayCountBestMode.sum()
#         self.arrayOfTapproxOverlaps = np.array(self.listOfTApproxOverlaps)
#         self.procentualRMSDreductionTapprox = self.calcProcentualReductionArray(self.listofRMSDbefore, self.listofRMSDafterTapprox)
#         
#         # RMSD reduction results
#         self.RMSDReductionsWholeperMode = self.getArrayOfColumns(self.listofRMSDReductionsWhole)
#         self.RMSDReductionsInterfaceperMode = self.getArrayOfColumns(self.listofRMSDReductionsInterface)
#         
#         self.RMSDReductionsWholeperModeProcentual = self.getProcentualRMSDReductions(self.listofRMSDbefore, self.RMSDReductionsWholeperMode)
#         self.RMSDReductionsInterfaceperModeProcentual = self.getProcentualRMSDReductions(self.listofRMSDInterfacebefore, self.RMSDReductionsInterfaceperMode)
#         
#         self.listofMeansWhole, self.listofStdWhole = self.getMeansandStd(self.RMSDReductionsWholeperModeProcentual)
#         self.listofMeansInterface, self.listofStdInterface = self.getMeansandStd(self.RMSDReductionsInterfaceperModeProcentual)
#         
#         # Tapprox overlap results
#         self.listofoverlapTapproxWholeperMode = self.getArrayOfColumns(self.listofoverlapTApproxWhole)
#         self.listofoverlapTapproxWholeMeans, self.listofoverlapTapproxWholeStd = self.getMeansandStd(self.listofoverlapTapproxWholeperMode)

# #         # TODO Plot 1, RMSD reduction
# #         self.RMSDReductionsWholeperMode = self.getArrayOfColumns2(self.listofRMSDReductionsWhole)
# #         self.RMSDReductionsInterfaceperMode = self.getArrayOfColumns2(self.listofRMSDReductionsInterface)
# #         
# #         self.RMSDReductionsWholeperModeProcentual = self.getArrayOfColumns2(self.listofArraysofRMSDReductionsWholePerc)
# #         self.RMSDReductionsInterfaceperModeProcentual = self.getArrayOfColumns2(self.listofArraysofRMSDReductionsInterfacePerc)      
# #             #to plot
# #         self.listofMeansWhole, self.listofStdWhole = self.getMeansandStd(self.RMSDReductionsWholeperModeProcentual)
# #         self.listofMeansInterface, self.listofStdInterface = self.getMeansandStd(self.RMSDReductionsInterfaceperModeProcentual)
        
        print "number of modes: ", self.listOfNumberofModes
        print "before plot 2 calculations"
        # Plot 2, chance to include the best overlap 
        self.chanceToIncludeBestOverlapModeWhole, self.chanceToIncludeBestOverlapModeInterface = self.getChanceToIncludeBestMode(self.listofArraysofModeOverlaps, self.listofArraysofModeOverlapsInterface)
        print "after plot 2 calculations"
        # Plot 3, how good is the overlap, histogram of x axis overlap, y axis percentage of the best
        self.discretizedCountHighestOverlapWholePercentage, self.discretizedCountHighestOverlapInterfacePercentage = self.getPercentageDistribution(self.listofArraysofModeOverlaps, self.listofArraysofModeOverlapsInterface) 
        print "after plot 3 calculations"
        # Plot 4, chance to include the best collectivity
        self.chanceToIncludeBestCollectivityModeWhole, self.chanceToIncludeBestCollectivityModeInterface = self.getChanceToIncludeBestMode(self.listofArraysofModeCollectivities, self.listofArraysofModeCollectivitiesInterface)
        print "after plot 4 calculations"
        # Plot 5, how good is the collectivity, histogram of x axis collectivity, y axis percentage of the best
        self.discretizedCountHighestCollectivityWholePercentage, self.discretizedCountHighestCollectivityInterfacePercentage = self.getPercentageDistribution(self.listofArraysofModeCollectivities, self.listofArraysofModeCollectivitiesInterface)
        print "after plot 5 calculations"
        # Plot 6, chance to include the best correlation
        self.chanceToIncludeBestCorrelationModeWhole, self.chanceToIncludeBestCorrelationModeInterface = self.getChanceToIncludeBestMode(self.listofArraysofModeCorrelations, self.listofArraysofModeCorrelationsInterface)
        print "after plot 6 calculations"
        # Plot 7, how good is the correlation, histogram of x axis correlation, y axis percentage of the best
        self.discretizedCountHighestCorrelationWholePercentage, self.discretizedCountHighestCorrelationInterfacePercentage = self.getPercentageDistribution(self.listofArraysofModeCorrelations, self.listofArraysofModeCorrelationsInterface)
        print "after plot 7 calculations"
        
        # Plot 3c, how good is the overlap assessed by the collectivity, histogram of x axis collectivity[overlap maxpos] , y axis percentage of the best collectivity
        self.discretizedCountHighestOverlapAssessedAsCollectivityWholePercentage = self.getPercentageDistributionOfMeasureFromCollectivity(self.listofArraysofModeOverlaps, self.listofArraysofModeCollectivities)
        self.discretizedCountHighestOverlapAssessedAsCollectivityInterfacePercentage = self.getPercentageDistributionOfMeasureFromCollectivity(self.listofArraysofModeOverlapsInterface, self.listofArraysofModeCollectivitiesInterface)
        print "after plot 3c calculations"
        # Plot 7c, how good is the correlation assessed by the collectivity, histogram of x axis collectivity[correlation maxpos] , y axis percentage of the best collectivity
        self.discretizedCountHighestCorrelationAssessedAsCollectivityWholePercentage = self.getPercentageDistributionOfMeasureFromCollectivity(self.listofArraysofModeCorrelations, self.listofArraysofModeCollectivities)
        self.discretizedCountHighestCorrelationAssessedAsCollectivityInterfacePercentage = self.getPercentageDistributionOfMeasureFromCollectivity(self.listofArraysofModeCorrelationsInterface, self.listofArraysofModeCollectivitiesInterface)
        print "after plot 7c calculations"

        print 'end cumulative assessment '
    def countOverlapCorrelationAgreement(self, overlapArray, correlationArray, overlapArrayInterface, correlationArrayInterface):
        """ Increment agreement variables of highest overlap and correlation if maxpos is the same. """
        if self.isMaxAtSameIndex(overlapArray, correlationArray):
            self.countHighestOverlapHighestCorrelationWhole += 1
        if self.isMaxAtSameIndex(overlapArrayInterface, correlationArrayInterface):
            self.countHighestOverlapHighestCorrelationInterface += 1     

    def isMaxAtSameIndex(self, arr1, arr2):
        """ Returns true if the max of arr1 is at the same index as the max of arr2, else returns false. 
            Internally, this method uses the absolute values of the arrays provided. 
            This test assumes the arrays have the same length, they are calculated on the same protein with the same number
            of maximum modes
            """
        assert len(arr1) == len(arr2)
        if np.argmax(np.abs(arr1)) == np.argmax(np.abs(arr2)):
            return True
        else:
            return False

    def getChanceToIncludeBestMode(self, listOfArraysWhole, listOfArraysInterface):
        """ Return arrays with the probability (or chance) to obtain the mode that is associated with the
        argmax() in the listOfArrays method arguments. 
        For example, the chance to include the mode of best (or highest) absolute overlap, collectivity or absolute correlation.
        Internally, this method uses the absolute values of the arrays provided.
        
        Args:
            listOfArraysWhole: a list of np.arrays with information on the whole proteins
            listOfArraysInterface: a list of np.arrays with information (for example overlap) on the interfaces
            
        Returns: two arrays with the chance to include the best mode, first for the whole protein, the second for the interface      
        """
        indexCountsOfHighestWhole     = np.zeros(np.max(self.listOfNumberofModes), dtype=np.int64)
        indexCountsOfHighestInterface = np.zeros(np.max(self.listOfNumberofModes), dtype=np.int64)
        for arr in listOfArraysWhole:
            indexCountsOfHighestWhole[np.argmax(np.abs(arr))] += 1
        for arrInterface in listOfArraysInterface:
            indexCountsOfHighestInterface[np.argmax(np.abs(arrInterface))] += 1       
        chanceToIncludeBestModeWhole     = (indexCountsOfHighestWhole.astype(float) / indexCountsOfHighestWhole.sum()).cumsum()
        chanceToIncludeBestModeInterface = (indexCountsOfHighestInterface.astype(float) / indexCountsOfHighestInterface.sum()).cumsum()
        return chanceToIncludeBestModeWhole, chanceToIncludeBestModeInterface
    
    def getPercentageDistribution(self, listOfArrays, listOfArraysInterface):
        """ Returns a discretized percentage distribution (in 0.1 steps) of the max() in the listOfArrays method arguments. 
        For example, this can be plotted as a histogram, with the x axis as overlap, collectivity or correlation, and the
        y axis the percentage of the best (highest).
        Internally, this method uses the absolute values of the arrays provided.
        
        Args:
            listOfArrays: a list of np.arrays with information on the whole proteins
            listOfArraysInterface: a list of np.arrays with information (for example overlap) on the interfaces
            
        Returns: two arrays with discretized percentage distribution of the max(), first for the whole protein, the second for the interface    
        
        """
        highestWhole = []
        highestInterface = []
        for arr in listOfArrays:
            highestWhole.append(np.abs(arr).max())
        for arrInterface in listOfArraysInterface:
            highestInterface.append(np.abs(arrInterface).max())   
        highestWholeDiscretized = self.calcDiscretizedCount01(highestWhole)
        highestInterfaceDiscretized = self.calcDiscretizedCount01(highestInterface) 
        highestWholeDiscretizedPercentage = (highestWholeDiscretized.astype(np.float64) / highestWholeDiscretized.sum()) * 100
        highestInterfaceDiscretizedPercentage = (highestInterfaceDiscretized.astype(np.float64) / highestInterfaceDiscretized.sum()) * 100
        return highestWholeDiscretizedPercentage, highestInterfaceDiscretizedPercentage
    
    def getPercentageDistributionOfMeasureFromCollectivity(self, listOfArrays, listOfCollectivities):
        """ Returns a discretized percentage distribution (in 0.1 steps) of the values from listOfArrays based on 
        argmax(listOfCollectivities). For example, this can be plotted as a histogram, with the x axis as overlap, 
        collectivity or correlation assessed by the highest collectivity, and the
        y axis the percentage of the best (highest) collectivity.
        Internally, this method uses the absolute values of the arrays provided.
        """
        assert len(listOfArrays) == len(listOfCollectivities)
        highestWhole = []
        for arr, arrCollectivity in zip(listOfArrays, listOfCollectivities):
            assert arr.shape == arrCollectivity.shape
            highestCollectivityPosition = np.argmax(np.abs(arrCollectivity))
            highestWhole.append(np.abs(arr[highestCollectivityPosition]))
        highestWholeDiscretized = self.calcDiscretizedCount01(highestWhole)
        highestWholeDiscretizedPercentage = (highestWholeDiscretized.astype(np.float64) / highestWholeDiscretized.sum()) * 100
        return highestWholeDiscretizedPercentage
    
#     def getPercentageDistributionOfMeasureFromCollectivity(self, listOfArrays, listOfCollectivities):
#         """ Returns a discretized percentage distribution (in 0.1 steps) of the collectivity based on the max() in the listOfArrays method arguments. 
#         For example, this can be plotted as a histogram, with the x axis as overlap, collectivity or correlation, and the
#         y axis the percentage of the best (highest) collectivity.
#         Internally, this method uses the absolute values of the arrays provided.
#         
#         Args:
#             listOfArrays: a list of np.arrays with information on the whole proteins
#             listOfCollectivities: a list of np.arrays with collectivity information
#             
#         Returns: an arrays with the discretized percentage distribution of the collectivity based on the max()
#         """
#         assert len(listOfArrays) == len(listOfCollectivities)
#         highestWhole = []
#         for arr, arrCollectivity in zip(listOfArrays, listOfCollectivities):
#             assert arr.shape == arrCollectivity.shape
#             print "arr: ", arr
#             print "arrCollectivity: ", arrCollectivity
#             highestModePosition = np.argmax(np.abs(arr))
#             print "highestposition: ", highestModePosition
#             print "arrCollectivity[highestModePosition]: ", arrCollectivity[highestModePosition]
#             print "************"
#             highestWhole.append(arrCollectivity[highestModePosition])    
#         highestWholeDiscretized = self.calcDiscretizedCount01(highestWhole)
#         highestWholeDiscretizedPercentage = (highestWholeDiscretized.astype(np.float64) / highestWholeDiscretized.sum()) * 100
#         return highestWholeDiscretizedPercentage

    def getProcentualRMSDReductions(self, listofRMSDbefore, inputArray):
        arrayofColumns = []
        for i in range (0, len(inputArray)):
            procentualRow = self.calcProcentualReductionArray(listofRMSDbefore, inputArray[i])
            arrayofColumns.append(procentualRow)
        return arrayofColumns
        
    def padArrayWithLastValue(self, myarray, maxmodes):
        assert maxmodes >= len(myarray)
        if len(myarray) == maxmodes:
            return myarray
        else:
            lastValue = myarray[-1]
            output = np.zeros(maxmodes, dtype=np.float64)
            for i in range(0, len(myarray)):
                output[i] = myarray[i]
            if maxmodes > len(myarray):
                for j in range(len(myarray), maxmodes):
                    output[j] = lastValue
            return output
    
    def getArrayOfColumns(self, inputArray):
        arrayofColumns = []
        for row in range (0, len(inputArray[0])):
            columns = np.zeros(len(inputArray), dtype=float64)
            for column in range(0, len(inputArray)):
                columns[column] = inputArray[column][row]
            arrayofColumns.append(columns)
        return arrayofColumns
    
    def getArrayOfColumns2(self, l):
        listofColumns = []
        for x in itertools.izip_longest(*l, fillvalue=None):
            listofColumns.append(np.array(filter(self.Exists, x)))
        return listofColumns
    
    def Exists(self, it):
        return (it is not None)
    
    def getMeansandStd(self, listOfArrays):
        means = np.zeros(len(listOfArrays), dtype=float64)
        std = np.zeros(len(listOfArrays), dtype=float64)
        for i in range(0, len(listOfArrays)):
            means[i] = listOfArrays[i].mean()
            std[i] = listOfArrays[i].std()
        return means, std
    
    def getArrayOfLastElements(self, listOfArrays):
        lastElements = np.zeros(len(listOfArrays), dtype=float64)
        for i in range(0, len(listOfArrays)):
            lastElements[i] = listOfArrays[i][-1:]
        return lastElements
            
    def printStatistics(self):
        # overview
        print "self.numberOfReferenceProteins: " + str(self.numberOfReferenceProteins)
#         print "self.numberOfModes: " + str(self.numberOfModes)
#         # plot 1, count of best modes
#         print "self.arrayCountBestMode: \n" + str(self.arrayCountBestMode)
#         # plot 2, count of best collectivity
#         print "self.arrayCountHighestCollectivity: \n" +str(self.arrayCountHighestCollectivity)
#         # plot 3, how good is overlap 
#         print "self.calcDiscretizedCount01(self.listOfHighestOverlaps) : \n" +str(self.calcDiscretizedCount01(self.listOfHighestOverlaps))
#         print "self.calcDiscretizedCount005(self.listOfHighestOverlaps) : \n" +str(self.calcDiscretizedCount005(self.listOfHighestOverlaps))
#         # plot 4, how good is overlap, average and standard deviation of overlap
#         print "self.arrayOfHighestOverlaps.mean(): " + str(self.arrayOfHighestOverlaps.mean())
#         print "self.arrayOfHighestOverlaps.std(): " + str(self.arrayOfHighestOverlaps.std())
#         print "self.sumHighestOverlapisWithin5LowestModes: " + str(self.sumHighestOverlapisWithin5LowestModes) + " , percentage: " + str(float(self.sumHighestOverlapisWithin5LowestModes)/self.numberOfReferenceProteins)
        # plot 5, chance to include best
#         print "self.arrayCountofBestModePercentageCumsum: \n" + str(self.arrayCountofBestModePercentage.cumsum())
#         print "self.calcCumulArray(self.arrayCountofBestModePercentage): \n" + str(self.calcCumulArray(self.arrayCountofBestModePercentage))
#         # plot 6, Tapprox (or make it table)
#         print "self.arrayOfTapproxOverlaps.mean(): \n" + str(self.arrayOfTapproxOverlaps.mean())
#         print "self.arrayOfTapproxOverlaps.std(): \n" + str(self.arrayOfTapproxOverlaps.std())
#         # plot 7, RMSD reduction
#         print "self.procentualRMSDreductionTapprox.mean(): \n" + str(self.procentualRMSDreductionTapprox.mean())
#         print "self.procentualRMSDreductionTapprox.std(): \n" + str(self.procentualRMSDreductionTapprox.std())
#         # plot 8, RMSD per mode
#         print "self.listofMeansWhole: \n" + str(self.listofMeansWhole)
#         print "list of RMSD percentage reductions non decreasing: " + str(self.non_decreasing(self.listofMeansWhole.tolist()))
#         print "self.RMSDReductionsWholeperModeProcentualLastMode: \n" + str(self.RMSDReductionsWholeperModeProcentual[-1:])
#         print "self.listofStdWhole: \n" + str(self.listofStdWhole)
#         print "self.listofMeansInterface: \n" + str(self.listofMeansInterface)
#         print "self.listofMeansInterface RMSD percentage reductions non decreasing: " + str(self.non_decreasing(self.listofMeansInterface.tolist()))
#         print "self.listofStdInterface: \n" + str(self.listofStdInterface)
#         # self.listofoverlapTapproxWholeMeans
#         print "self.listofoverlapTapproxWholeMeans: \n" + str(self.listofoverlapTapproxWholeMeans)
#         print "self.listofoverlapTapproxWholeStdDev: \n" + str(self.listofoverlapTapproxWholeStd)
#         print "self.listofoverlapTapproxWholeMeans non decreasing: " + str(self.non_decreasing(self.listofoverlapTapproxWholeMeans.tolist()))
        
        ## new statistics
#         print "new self.listofRMSDbefore: \n", self.listofRMSDbefore
#         print "new self.listofRMSDInterfacebefore: \n", self.listofRMSDInterfacebefore
#         print "new self.listofRMSDReductionsWhole: \n", self.listofRMSDReductionsWhole
#         print "new self.listofRMSDafterTapprox: \n", self.listofRMSDafterTapprox
#         print "new self.listofRMSDReductionsInterface: \n", self.listofRMSDReductionsInterface
#         print "new self.listofoverlapTApproxWhole : \n", self.listofoverlapTApproxWhole
#         print "new self.listofoverlapTApproxInterface: \n", self.listofoverlapTApproxInterface
#         print "new self.listofArraysofModeOverlaps: \n", self.listofArraysofModeOverlaps
#         print "new self.listofArraysofModeCollectivities: \n", self.listofArraysofModeCollectivities
#         print "new self.listofArraysofModeCorrelations: \n", self.listofArraysofModeCorrelations
#         print "new self.RMSDReductionsWholeperMode: \n", self.RMSDReductionsWholeperMode
#         print "new self.RMSDReductionsWholeperModeProcentual: \n", self.RMSDReductionsWholeperModeProcentual
#         print "new self.RMSDReductionsInterfaceperModeProcentual: \n", self.RMSDReductionsInterfaceperModeProcentual
#         print "new self.listofMeansWhole: \n", self.listofMeansWhole
#         print "new self.listofStdWhole: \n", self.listofStdWhole
#         print "new self.listofMeansInterface: \n", self.listofMeansInterface
#         print "new self.listofStdInterface: \n", self.listofStdInterface
# #         print "list of RMSD percentage reductions non decreasing                  : " + str(self.non_decreasing(self.listofMeansWhole.tolist()))
# #         print "self.listofMeansInterface RMSD percentage reductions non decreasing: " + str(self.non_decreasing(self.listofMeansInterface.tolist()))
        print "number of modes: ", self.listOfNumberofModes

    def writeStatistics(self, maxModesWithAllTuples, maxModesWithAllTuplesInterface, argumentNumberOfTuples, whatAtomsToMatch, whichBenchmarkPartChoosen):
        # print out all values
        np.set_printoptions(threshold=np.nan)
        
        # Make the base path 
        path = self.utils.config.outputPath + self.experimentName + "/"
        self.utils.mkdir_p(path)
        
        # Write one liner statistics into a text file
        f = open(path+self.experimentName+"basics.txt", 'wt')
        f.write("#"+path+self.experimentName+"basics.txt : \n")
#         f.write("RMSDReductionsWholeNonDecreasing                           : " + str(self.non_decreasing(self.listofMeansWhole.tolist())) + "\n")
#         f.write("RMSDReductionsInterfacENonDecreasing                       : " + str(self.non_decreasing(self.listofMeansInterface.tolist())) + "\n")
        f.write("maxModesWithAllTuples                                      : " + str(maxModesWithAllTuples) + "\n")
        f.write("maxModesWithAllTuplesInterface                             : " + str(maxModesWithAllTuplesInterface) + "\n")
        f.write("self.countHighestOverlapHighestCorrelationWhole            : " + str(self.countHighestOverlapHighestCorrelationWhole) + "\n")
        f.write("self.countHighestOverlapHighestCorrelationInterface        : " + str(self.countHighestOverlapHighestCorrelationInterface) + "\n")
        f.write("self.numberOfReferenceProteins                                        : " + str(self.numberOfReferenceProteins) + "\n")
        f.write("argumentNumberOfTuples                                     : " + str(argumentNumberOfTuples) + "\n")
        f.write("self.countHighestOverlapHighestCorrelationWhole      /#t   : " + str(float(self.countHighestOverlapHighestCorrelationWhole) / float(self.numberOfReferenceProteins)) + "\n")
        f.write("self.countHighestOverlapHighestCorrelationInterface  /#t   : " + str(float(self.countHighestOverlapHighestCorrelationInterface) / float(self.numberOfReferenceProteins)) + "\n")
        f.write("whatAtomsToMatch                                           : " + str(whatAtomsToMatch) + "\n")
        f.write("whichBenchmarkPartChoosen                                  : " + str(whichBenchmarkPartChoosen) + "\n")
        f.close()
        
        # Output path for data in numpy format that can be plotted
        pathnp = path+"plotnp/"
        self.utils.mkdir_p(pathnp)
        
        ## Write the raw data to files
        self.writeToFile(pathnp, "self.numberOfReferenceProteins", self.numberOfReferenceProteins)
        
        # Plot 1 RMSD reductions whole and interface, mean and stddev
#         self.writeToFile(pathnp, "self.listofMeansWhole", self.listofMeansWhole)
#         self.writeToFile(pathnp, "self.listofStdWhole", self.listofStdWhole)
#         self.writeToFile(pathnp, "self.listofMeansInterface", self.listofMeansInterface)
#         self.writeToFile(pathnp, "self.listofStdInterface", self.listofStdInterface)
#         self.writeToFile(pathnp, "argumentNumberOfTuples", argumentNumberOfTuples)
        
        # Plot 2, chance to include the best overlap
        self.writeToFile(pathnp, "self.chanceToIncludeBestOverlapModeWhole", self.chanceToIncludeBestOverlapModeWhole)
        self.writeToFile(pathnp, "self.chanceToIncludeBestOverlapModeInterface", self.chanceToIncludeBestOverlapModeInterface)
        
        # Plot 3, how good is the overlap, histogram of x axis overlap, y axis percentage of the best
        self.writeToFile(pathnp, "self.discretizedCountHighestOverlapWholePercentage", self.discretizedCountHighestOverlapWholePercentage)
        self.writeToFile(pathnp, "self.discretizedCountHighestOverlapInterfacePercentage", self.discretizedCountHighestOverlapInterfacePercentage)
        
        # Plot 4, chance to include the best collectivity
        self.writeToFile(pathnp, "self.chanceToIncludeBestCollectivityModeWhole", self.chanceToIncludeBestCollectivityModeWhole)
        self.writeToFile(pathnp, "self.chanceToIncludeBestCollectivityModeInterface", self.chanceToIncludeBestCollectivityModeInterface)
        
        # Plot 5, how good is the collectivity, histogram of x axis collectivity, y axis percentage of the best
        self.writeToFile(pathnp, "self.discretizedCountHighestCollectivityWholePercentage", self.discretizedCountHighestCollectivityWholePercentage)
        self.writeToFile(pathnp, "self.discretizedCountHighestCollectivityInterfacePercentage", self.discretizedCountHighestCollectivityInterfacePercentage)
        
        # Plot 6, chance to include the best correlation
        self.writeToFile(pathnp, "self.chanceToIncludeBestCorrelationModeWhole", self.chanceToIncludeBestCorrelationModeWhole)
        self.writeToFile(pathnp, "self.chanceToIncludeBestCorrelationModeInterface", self.chanceToIncludeBestCorrelationModeInterface)
        
        # Plot 7, how good is the correlation, histogram of x axis correlation, y axis percentage of the best
        self.writeToFile(pathnp, "self.discretizedCountHighestCorrelationWholePercentage", self.discretizedCountHighestCorrelationWholePercentage)
        self.writeToFile(pathnp, "self.discretizedCountHighestCorrelationInterfacePercentage", self.discretizedCountHighestCorrelationInterfacePercentage)

        self.writeToFile(pathnp, "self.listOfNumberofModes", self.listOfNumberofModes)

#         # Plot 3b, how good is the overlap assessed by the collectivity, histogram of x axis collectivity[overlap maxpos] , y axis percentage of the best collectivity
        self.writeToFile(pathnp, "self.discretizedCountHighestOverlapAssessedAsCollectivityWholePercentage", self.discretizedCountHighestOverlapAssessedAsCollectivityWholePercentage)
        self.writeToFile(pathnp, "self.discretizedCountHighestOverlapAssessedAsCollectivityInterfacePercentage", self.discretizedCountHighestOverlapAssessedAsCollectivityInterfacePercentage) 

#         # Plot 7b, how good is the correlation assessed by the collectivity, histogram of x axis collectivity[correlation maxpos] , y axis percentage of the best collectivity
        self.writeToFile(pathnp, "self.discretizedCountHighestCorrelationAssessedAsCollectivityWholePercentage", self.discretizedCountHighestCorrelationAssessedAsCollectivityWholePercentage)
        self.writeToFile(pathnp, "self.discretizedCountHighestCorrelationAssessedAsCollectivityInterfacePercentage", self.discretizedCountHighestCorrelationAssessedAsCollectivityInterfacePercentage)


        # Output path for raw data, mostly list of arrays
        pathraw = path+"raw/"
        self.utils.mkdir_p(pathraw)
#         self.writeToFile(pathraw, "self.listofRMSDbefore", self.listofRMSDbefore)
#         self.writeToFile(pathraw, "self.listofRMSDInterfacebefore", self.listofRMSDInterfacebefore)
#         self.writeToFile(pathraw, "self.listofRMSDReductionsWhole", self.listofRMSDReductionsWhole)
#         self.writeToFile(pathraw, "self.listofRMSDReductionsInterface", self.listofRMSDReductionsInterface)
#         self.writeToFile(pathraw, "self.listofoverlapTApproxWhole", self.listofoverlapTApproxWhole)
#         self.writeToFile(pathraw, "self.listofoverlapTApproxInterface", self.listofoverlapTApproxInterface)
        self.writeToFile(pathraw, "self.listofArraysofModeOverlaps", self.listofArraysofModeOverlaps)
        self.writeToFile(pathraw, "self.listofArraysofModeOverlapsInterface", self.listofArraysofModeOverlapsInterface)
        self.writeToFile(pathraw, "self.listofArraysofModeCollectivities", self.listofArraysofModeCollectivities)
        self.writeToFile(pathraw, "self.listofArraysofModeCollectivitiesInterface", self.listofArraysofModeCollectivitiesInterface)
        self.writeToFile(pathraw, "self.listofArraysofModeCorrelations", self.listofArraysofModeCorrelations)
        self.writeToFile(pathraw, "self.listofArraysofModeCorrelationsInterface", self.listofArraysofModeCorrelationsInterface)
#         self.writeToFile(pathraw, "self.RMSDReductionsWholeperMode", self.RMSDReductionsWholeperMode)
#         self.writeToFile(pathraw, "self.RMSDReductionsInterfaceperMode", self.RMSDReductionsInterfaceperMode)
#         self.writeToFile(pathraw, "self.RMSDReductionsWholeperModeProcentual", self.RMSDReductionsWholeperModeProcentual)
#         self.writeToFile(pathraw, "self.RMSDReductionsInterfaceperModeProcentual", self.RMSDReductionsInterfaceperModeProcentual)
        self.writeToFile(pathraw, "self.listofArraysofEigenvalues", self.listofArraysofEigenvalues)
        
        # for list of highest overlap all tuples
        self.writeToFile(pathraw, "self.listofWorkedOnTuples", self.listofWorkedOnTuples)

        ## Write output using ProDys writeArray method that can be plotted
        pathpr = path+"plotprody/"
        self.utils.mkdir_p(pathpr)
        
        # Plot 1 RMSD reductions whole and interface, mean and stddev
#         writeArray(pathpr+self.experimentName+"listofMeansWhole.txt", self.listofMeansWhole, format='%f', delimiter=' ')
#         writeArray(pathpr+self.experimentName+"listofStdWhole.txt", self.listofStdWhole, format='%f', delimiter=' ')
#         writeArray(pathpr+self.experimentName+"listofMeansInterface.txt", self.listofMeansInterface, format='%f', delimiter=' ')
#         writeArray(pathpr+self.experimentName+"listofStdInterface.txt", self.listofStdInterface, format='%f', delimiter=' ')
        
        # Plot 2, chance to include the best overlap
        writeArray(pathpr+self.experimentName+"chanceToIncludeBestOverlapModeWhole.txt", self.chanceToIncludeBestOverlapModeWhole, format='%f', delimiter=' ')
        writeArray(pathpr+self.experimentName+"chanceToIncludeBestOverlapModeInterface.txt", self.chanceToIncludeBestOverlapModeInterface, format='%f', delimiter=' ')

        # Plot 3, how good is the overlap, histogram of x axis overlap, y axis percentage of the best
        writeArray(pathpr+self.experimentName+"discretizedCountHighestOverlapWholePercentage.txt", self.discretizedCountHighestOverlapWholePercentage, format='%f', delimiter=' ')
        writeArray(pathpr+self.experimentName+"discretizedCountHighestOverlapInterfacePercentage.txt", self.discretizedCountHighestOverlapInterfacePercentage, format='%f', delimiter=' ')

        # Plot 4, chance to include the best collectivity
        writeArray(pathpr+self.experimentName+"chanceToIncludeBestCollectivityModeWhole.txt", self.chanceToIncludeBestCollectivityModeWhole, format='%f', delimiter=' ')
        writeArray(pathpr+self.experimentName+"chanceToIncludeBestCollectivityModeInterface.txt", self.chanceToIncludeBestCollectivityModeInterface, format='%f', delimiter=' ')

        # Plot 5, how good is the collectivity, histogram of x axis collectivity, y axis percentage of the best
        writeArray(pathpr+self.experimentName+"discretizedCountHighestCollectivityWholePercentage.txt", self.discretizedCountHighestCollectivityWholePercentage, format='%f', delimiter=' ')
        writeArray(pathpr+self.experimentName+"discretizedCountHighestCollectivityInterfacePercentage.txt", self.discretizedCountHighestCollectivityInterfacePercentage, format='%f', delimiter=' ')

        # Plot 6, chance to include the best correlation
        writeArray(pathpr+self.experimentName+"chanceToIncludeBestCorrelationModeWhole.txt", self.chanceToIncludeBestCorrelationModeWhole, format='%f', delimiter=' ')
        writeArray(pathpr+self.experimentName+"chanceToIncludeBestCorrelationModeInterface.txt", self.chanceToIncludeBestCorrelationModeInterface, format='%f', delimiter=' ')

        # Plot 7, how good is the correlation, histogram of x axis correlation, y axis percentage of the best
        writeArray(pathpr+self.experimentName+"discretizedCountHighestCorrelationWholePercentage.txt", self.discretizedCountHighestCorrelationWholePercentage, format='%f', delimiter=' ')
        writeArray(pathpr+self.experimentName+"discretizedCountHighestCorrelationInterfacePercentage.txt", self.discretizedCountHighestCorrelationInterfacePercentage, format='%f', delimiter=' ')
        
#         # Statistical Tests, last results of lowest RMSD
#         writeArray(pathpr+self.experimentName+"self.RMSDReductionsWholeperMode[-1]", self.RMSDReductionsWholeperMode[-1], format='%f', delimiter=' ')
#         writeArray(pathpr+self.experimentName+"self.RMSDReductionsInterfaceperMode[-1]", self.RMSDReductionsInterfaceperMode[-1], format='%f', delimiter=' ')
#         writeArray(pathpr+self.experimentName+"self.RMSDReductionsWholeperModeProcentual[-1]", self.RMSDReductionsWholeperModeProcentual[-1], format='%f', delimiter=' ')
#         writeArray(pathpr+self.experimentName+"self.RMSDReductionsInterfaceperModeProcentual[-1]", self.RMSDReductionsInterfaceperModeProcentual[-1], format='%f', delimiter=' ')

        writeArray(pathpr+self.experimentName+"listOfNumberofModes", np.array(self.listOfNumberofModes), format='%d', delimiter=' ')

#         # Plot 3b, how good is the overlap assessed by the collectivity, histogram of x axis collectivity[overlap maxpos] , y axis percentage of the best collectivity
        writeArray(pathpr+self.experimentName+"discretizedCountHighestOverlapAssessedAsCollectivityWholePercentage", self.discretizedCountHighestOverlapAssessedAsCollectivityWholePercentage, format='%f', delimiter=' ')
        writeArray(pathpr+self.experimentName+"discretizedCountHighestOverlapAssessedAsCollectivityInterfacePercentage", self.discretizedCountHighestOverlapAssessedAsCollectivityInterfacePercentage, format='%f', delimiter=' ') 

        # Plot 7b, how good is the correlation assessed by the collectivity, histogram of x axis collectivity[correlation maxpos] , y axis percentage of the best collectivity
        writeArray(pathpr+self.experimentName+"discretizedCountHighestCorrelationAssessedAsCollectivityWholePercentage", self.discretizedCountHighestCorrelationAssessedAsCollectivityWholePercentage, format='%f', delimiter=' ')
        writeArray(pathpr+self.experimentName+"discretizedCountHighestCorrelationAssessedAsCollectivityInterfacePercentage", self.discretizedCountHighestCorrelationAssessedAsCollectivityInterfacePercentage, format='%f', delimiter=' ')

        # Plot for highest overlap and position
        writeArray(pathpr+self.experimentName+"self.listofHighestOverlapsPositions", np.array(self.listofHighestOverlapsPositions), format='%d', delimiter=' ')
        writeArray(pathpr+self.experimentName+"self.listOfHighestOverlaps",  np.array(self.listOfHighestOverlaps), format='%f', delimiter=' ')
        writeArray(pathpr+self.experimentName+"self.listofHighestOverlapsInterfacePositions",  np.array(self.listofHighestOverlapsInterfacePositions), format='%d', delimiter=' ')
        writeArray(pathpr+self.experimentName+"self.listOfHighestOverlapsInterface",  np.array(self.listOfHighestOverlapsInterface), format='%f', delimiter=' ')

#         f.write("\nself.numberOfModes: " + str(self.numberOfModes))
#         # plot 1
#         f.write("\nself.arrayCountBestMode: " + str(self.arrayCountBestMode))
#         # plot 2, count of best collectivity
#         f.write("\nself.arrayCountHighestCollectivity: " +str(self.arrayCountHighestCollectivity))
#         # plot 3, how good is overlap 
#         f.write("\nself.calcDiscretizedCount01(self.listOfHighestOverlaps): " +str(self.calcDiscretizedCount01(self.listOfHighestOverlaps)))
#         f.write("\nself.calcDiscretizedCount005(self.listOfHighestOverlaps): " +str(self.calcDiscretizedCount005(self.listOfHighestOverlaps)))
#         # plot 4, how good is overlap, average and standard deviation of overlap
#         f.write("\nself.arrayOfHighestOverlaps.mean(): " + str(self.arrayOfHighestOverlaps.mean()))
#         f.write("\nself.arrayOfHighestOverlaps.std(): " + str(self.arrayOfHighestOverlaps.std()))
#         f.write("\nself.sumHighestOverlapisWithin5LowestModes: " + str(self.sumHighestOverlapisWithin5LowestModes) + " , percentage: " + str(float(self.sumHighestOverlapisWithin5LowestModes)/self.numberOfReferenceProteins))
#         # plot 5, chance to include best
#         f.write("\nself.arrayCountofBestModePercentageCumsum: " + str(self.arrayCountofBestModePercentage.cumsum()))
#         f.write("\nself.calcCumulArray(self.arrayCountofBestModePercentage): " + str(self.calcCumulArray(self.arrayCountofBestModePercentage)))
#         # plot 6, Tapprox (or make it table)
#         f.write("\nself.arrayOfTapproxOverlaps.mean(): " + str(self.arrayOfTapproxOverlaps.mean()))
#         f.write("\nself.arrayOfTapproxOverlaps.std(): " + str(self.arrayOfTapproxOverlaps.std()))
#         # plot 7, RMSD reduction
#         f.write("\nself.procentualRMSDreductionTapprox.mean(): " + str(self.procentualRMSDreductionTapprox.mean()))
#         f.write("\nself.procentualRMSDreductionTapprox.std(): " + str(self.procentualRMSDreductionTapprox.std()))
#         # plot 8, RMSD per mode
#         f.write("\nself.listofMeansWhole: " + str(self.listofMeansWhole))
#         f.write("\nlist of RMSD percentage reductions non decreasing: " + str(self.non_decreasing(self.listofMeansWhole.tolist())))
#         f.write("\nself.RMSDReductionsWholeperModeProcentualLastMode: " + str(self.RMSDReductionsWholeperModeProcentual[-1:]))
#         f.write("\nself.listofStdWhole: " + str(self.listofStdWhole))
#         f.write("\nself.listofMeansInterface: " + str(self.listofMeansInterface))
#         f.write("\nself.listofMeansInterface non decreasing: " + str(self.non_decreasing(self.listofMeansInterface.tolist())))
#         f.write("\nself.listofStdInterface: " + str(self.listofStdInterface))
#         # self.listofoverlapTapproxWholeMeans
#         f.write("\nself.listofoverlapTapproxWholeMeans: " + str(self.listofoverlapTapproxWholeMeans))
#         f.write("\nself.listofoverlapTapproxWholeStdDev: " + str(self.listofoverlapTapproxWholeStd))
#         f.write("\nself.listofoverlapTapproxWholeMeans non decreasing: " + str(self.non_decreasing(self.listofoverlapTapproxWholeMeans.tolist())))
#        f.close()
#         writeArray(path+self.experimentName+"arrayCountBestMode.txt", self.arrayCountBestMode, format='%d', delimiter=' ')
#         writeArray(path+self.experimentName+"arrayCountHighestCollectivity.txt", self.arrayCountHighestCollectivity, format='%d', delimiter=' ')
#         writeArray(path+self.experimentName+"calcDiscretizedCount01(self.listOfHighestOverlaps).txt", self.calcDiscretizedCount01(self.listOfHighestOverlaps), format='%d', delimiter=' ')
#         writeArray(path+self.experimentName+"calcDiscretizedCount005(self.listOfHighestOverlaps).txt", self.calcDiscretizedCount005(self.listOfHighestOverlaps), format='%d', delimiter=' ')
#         writeArray(path+self.experimentName+"self.arrayCountofBestModePercentagecumsum.txt", self.arrayCountofBestModePercentage.cumsum(), format='%f', delimiter=' ')
#         writeArray(path+self.experimentName+"self.calcCumulArray(self.arrayCountofBestModePercentage).txt", self.calcCumulArray(self.arrayCountofBestModePercentage), format='%f', delimiter=' ')
#         writeArray(path+self.experimentName+"oldarrayAVGSumofAbsOverlaps.txt", (self.arraySumofAbsOverlaps/self.numberOfReferenceProteins), format='%f', delimiter=' ')
#         writeArray(path+self.experimentName+"oldarrayAVGSumofAbsCollectivities.txt", (self.arraySumofAbsCollectivities/self.numberOfReferenceProteins), format='%f', delimiter=' ')
#         writeArray(path+self.experimentName+"oldarrayAVGCumulSumofAbsOverlap.txt", (self.arrayCumulSumofAbsOverlap/self.numberOfReferenceProteins), format='%f', delimiter=' ')
#         writeArray(path+self.experimentName+"oldarrayAVGCumulSumofAbsCollectivities.txt", (self.arrayCumulSumofAbsCollectivities/self.numberOfReferenceProteins), format='%f', delimiter=' ')
#         writeArray(path+self.experimentName+"listofMeansWhole.txt", self.listofMeansWhole, format='%f', delimiter=' ')
#         writeArray(path+self.experimentName+"listofStdWhole.txt", self.listofStdWhole, format='%f', delimiter=' ')
#         writeArray(path+self.experimentName+"listofMeansInterface.txt", self.listofMeansInterface, format='%f', delimiter=' ')
#         writeArray(path+self.experimentName+"listofStdInterface.txt", self.listofStdInterface, format='%f', delimiter=' ')
#         writeArray(path+self.experimentName+"listofoverlapTapproxWholeMeans.txt", self.listofoverlapTapproxWholeMeans, format='%f', delimiter=' ')
#         writeArray(path+self.experimentName+"listofoverlapTapproxWholeStdDev.txt", self.listofoverlapTapproxWholeStd, format='%f', delimiter=' ')
        #self.showRMSDReductions(self.listofMeansWhole, self.listofoverlapTapproxWholeMeans)
        ## new statistics
#         writeArray(path+self.experimentName+"newself.listofRMSDbefore.txt", np.array(self.listofRMSDbefore), format='%f', delimiter=' ')
#         writeArray(path+self.experimentName+"newself.listofRMSDInterfacebefore.txt", np.array(self.listofRMSDInterfacebefore), format='%f', delimiter=' ')
#         self.writeArrayofArraystoFile(path+self.experimentName+"MYnewself.listofRMSDReductionsWhole.txt", np.array(self.listofRMSDReductionsWhole))
#         writeArray(path+self.experimentName+"newself.listofRMSDReductionsWhole.txt", np.array(self.listofRMSDReductionsWhole), format='%f', delimiter=' ')
#         writeArray(path+self.experimentName+"newself.listofRMSDafterTapprox.txt", np.array(self.listofRMSDafterTapprox), format='%f', delimiter=' ')
#         writeArray(path+self.experimentName+"newself.listofRMSDReductionsInterface.txt", np.array(self.listofRMSDReductionsInterface), format='%f', delimiter=' ')
#         writeArray(path+self.experimentName+"newself.listofoverlapTApproxWhole.txt", np.array(self.listofoverlapTApproxWhole), format='%f', delimiter=' ')
#         writeArray(path+self.experimentName+"newself.listofoverlapTApproxInterface.txt", np.array(self.listofoverlapTApproxInterface), format='%f', delimiter=' ')
#         writeArray(path+self.experimentName+"newself.listofArraysofModeOverlaps.txt", np.array(self.listofArraysofModeOverlaps), format='%f', delimiter=' ')
#         writeArray(path+self.experimentName+"newself.listofArraysofModeCollectivities.txt", np.array(self.listofArraysofModeCollectivities), format='%f', delimiter=' ')
#         writeArray(path+self.experimentName+"newself.listofArraysofModeCorrelations.txt", np.array(self.listofArraysofModeCorrelations), format='%f', delimiter=' ')
    
    def writeDirectOuput(self, pdbQueueItem, overlapArrayWhole, collectivityArrayWhole, RMSD_unbound_to_superposed_bound, RMSD_interface, RMSDReductionsWhole, RMSDReductionsInterface):
        np.set_printoptions(threshold=np.nan)
        path = self.utils.config.outputPath+self.experimentName+"/"+pdbQueueItem[0]+"/"
        self.utils.mkdir_p(path)
        np.savetxt(path+"pdbQueueItem.txt", [pdbQueueItem], fmt='%s', header='PDB name')
        np.savetxt(path+"overlapArrayWhole.txt", overlapArrayWhole, fmt='%f', header='overlapArrayWhole')
        np.savetxt(path+"collectivityArrayWhole.txt", collectivityArrayWhole, fmt='%f', header='collectivityArrayWhole')
        np.savetxt(path+"RMSD_unbound_to_superposed_bound.txt", [RMSD_unbound_to_superposed_bound], fmt='%f', header='RMSD_unbound_to_superposed_bound')
        np.savetxt(path+"RMSD_interface.txt", [RMSD_interface], fmt='%f', header='RMSD_interface')
        np.savetxt(path+"RMSDReductionsWhole.txt", RMSDReductionsWhole, fmt='%f', header='RMSDReductionsWhole')
        np.savetxt(path+"RMSDReductionsInterface.txt", RMSDReductionsInterface, fmt='%f', header='RMSDReductionsInterface')
    
    def appendArrayToFile(self, fileName, singleArray):
        with open(fileName, "a") as myfile:
            np.savetxt(myfile, singleArray, fmt='%f', delimiter=' ')
            
    def writeArrayofArraystoFile(self, fileName, wholeArray):
        for singleArray in wholeArray:
            self.appendArrayToFile(fileName, singleArray)
            
    def writeToFile(self, path, fileName, data):
        f = open(path+self.experimentName+fileName+"raw.txt", 'wt')
        f.write("#"+fileName+" : \n" + repr(data))
        f.close()
        
    def showRMSDReductions(self, RMSDreductions, overlapWhole):
        plt.plot(RMSDreductions, linewidth=2.5, label="RMSD reduction")
        plt.plot(overlapWhole, linewidth=2.5, label="overlap")
        plt.legend()
        plt.show()
    
    def percentageDecrease(self, x, y):
        """Return the percentage decrease from x towards y
        
            Args:
                x: original value
                y: changed value
            
            Returns: percentage decrease from x towards y.
                
        """
        return (float(x)-float(y))/float(x)
    
    def calcPearson(self, compareVector, defvec):
        """ Calculate the pearson correlation coefficent between the magnitudes of displacement in 
            compareVector (can be a mode or an approximated deformation vector) and defvec (which is the true
            deformation vector.
            
            Args:
                compareVector: a mode or approximated deformation vector
                defvec: the true deformation vector
                
            Returns: The pearson correlation coefficent between the magnitudes of displacement of the two vectors 
                compareVector and defvec (between -1 and 1)
            
        """
        A = compareVector
        R = defvec
        Aarray = A.getArray()
        Rarray = R.getArray()
        assert(Aarray.shape == Rarray.shape)
        AarrayNx3 = A.getArrayNx3()
        RarrayNx3 = R.getArrayNx3()
        
        # make magnitudes
        AaverageDisplacementArray = np.array([])
        AaverageDisplacement = 0.0
        RaverageDisplacementArray = np.array([])
        RaverageDisplacement = 0.0
        for i in range(0, len(AarrayNx3)):
           AaverageDisplacement = np.sqrt(np.power(AarrayNx3[i][0], 2) +  np.power(AarrayNx3[i][1], 2) +  np.power(AarrayNx3[i][2], 2))
           RaverageDisplacement = np.sqrt(np.power(RarrayNx3[i][0], 2) +  np.power(RarrayNx3[i][1], 2) +  np.power(RarrayNx3[i][2], 2))
           AaverageDisplacementArray = np.append(AaverageDisplacementArray, AaverageDisplacement)
           RaverageDisplacementArray = np.append(RaverageDisplacementArray, RaverageDisplacement)
        return pearsonr(AaverageDisplacementArray, RaverageDisplacementArray)[0]
    
    def calcPearson2(self, estimatedMagnitudes, trueMagnitudes):
        """ Calculate the pearson correlation coefficent between the magnitudes of displacement in 
            given by estimatedMagnitudes and trueMagnitudes.
            
            Args:
                estimatedMagnitudes: array of magnitudes of a mode or approximated deformation vector
                trueMagnitudes: array of magnitudes of the true deformation vector
                
            Returns: The pearson correlation coefficent between estimatedMagnitudes and trueMagnitudes,
                    the  value is (between -1 and 1)
            
        """
        assert (estimatedMagnitudes.shape == trueMagnitudes.shape)
        return pearsonr(estimatedMagnitudes, trueMagnitudes)[0]
    
    def calcMagnitudeArray(self, myVector):
        """ Calculate an array with the magnitudes of the 3D structure described by vector. 
        
            Args:
                myVector: a vector for 3D coordinates, for example a deformation vector or 
                          a normal mode vector
                          
            Returns: array of magnitudes of myVector
        """
        assert myVector.is3d()
        myVectorNx3 = myVector.getArrayNx3()
        vectorMagnitudes = np.linalg.norm(myVectorNx3, axis=1)
        vectorMagnitudes = np.array(vectorMagnitudes)
        return vectorMagnitudes
        
    
    def setExperimentName(self, whatAtomsToMatch, typeOfExperiment=None):
        if not typeOfExperiment:
            self.experimentName = "numberOfReferences" + str(self.numberOfReferenceProteins) + "matchedOn" + str(whatAtomsToMatch) + datetime.now().strftime("-%d-%m-%Y--%H:%M:%S")            
        else:
            self.experimentName = typeOfExperiment + "NumberOfReferences" + str(self.numberOfReferenceProteins) + "matchedOn" + str(whatAtomsToMatch) + datetime.now().strftime("-%d-%m-%Y--%H:%M:%S")
        
    def non_decreasing(self, L):
        threshold = 0.00001
        result = all(x<= (y+threshold) for x, y in zip(L, L[1:]))
        assert result == True
        return result

    def non_increasing(self, L):
        threshold = 0.00001
        result = all(x>= (y-threshold) for x, y in zip(L, L[1:]))
        assert result == True
        return result
