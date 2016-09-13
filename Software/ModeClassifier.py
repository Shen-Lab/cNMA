'''
Created on Feb 16, 2014

@author: not
'''

import numpy as np
import sys as sys
from prody.dynamics.mode import Vector
from prody.dynamics.anm import ANM
from prody.dynamics.functions import writeArray

class ModeClassifier(object):
    '''
    Classifies a mode and its eigenvalue from a bound complex towards receptor or ligand. 
    '''


    def __init__(self, utils):
        '''
        Constructor
        '''
        self.utils = utils
        
    def calcValidModes(self, anm_slc, anm_slc_counterpart, reference):
        eigenvalueClassifications = []
        measureCounter = 0
        
        assert anm_slc.numModes() == anm_slc_counterpart.numModes()
        assert all(anm_slc.getEigvals() == anm_slc_counterpart.getEigvals())
        modeArray = anm_slc.getArray().T
        modeArrayCounterpart = anm_slc_counterpart.getArray().T
           
        for i in range(0, anm_slc.numModes()):
            complexMode = np.concatenate((modeArray[i], modeArrayCounterpart[i]))
            complexModeMagnitude = np.linalg.norm(complexMode)
            
            # calculate contribution measure
            refContribution = np.power(np.linalg.norm(modeArray[i]/complexModeMagnitude), 2.0)
            counterpartContribution = np.power(np.linalg.norm(modeArrayCounterpart[i]/complexModeMagnitude), 2.0)
            measure = refContribution - counterpartContribution
            if measure > 0:
                measureCounter += 1
            
            # eigenvalue readjusting
            complexContribution = np.concatenate((modeArray[i], modeArrayCounterpart[i]))
            complexContribution = np.power(np.linalg.norm(complexContribution/complexModeMagnitude), 2.0)
            eigScaling = refContribution / complexContribution
            eigValue = anm_slc.getEigvals()[i] * eigScaling
            
            eigenvalueClassifications.append([[modeArray[i]], eigValue, measure]) 
        
        # sort from biggest to smallest measure
        eigenvalueClassifications.sort(key=lambda x: x[2], reverse=True)
        assert eigenvalueClassifications[0][2] <= 1.01
        
        print "eigenvalueClassifications[0]: ", eigenvalueClassifications[0]
        print "eigenvalueClassifications[-1]: ", eigenvalueClassifications[-1]
        # create the anm array and the list of corresponding eigenvalues
        selectedModes = []
        selectedEigenVals = []
        for j in range(0, reference.select('calpha').numAtoms()*3 -6):
            selectedModes.append(eigenvalueClassifications[j][0][0])
            selectedEigenVals.append(eigenvalueClassifications[j][1])
        
        print "measureCounter: ", measureCounter
        #print "selectedModes: ", selectedModes
        return selectedModes, selectedEigenVals

#     def eigenvalueDecider(matrix, topDimension, bottomDimension):
#         eigenvecs = np.linalg.eigh(matrix)[1].T
#         eigClassifications = []
#         
#         for eigenvec in eigenvecs:
#             eigenvec = getUnitVector(eigenvec)
#             topContribution = np.linalg.norm(eigenvec[0:topDimension])
#             bottomContribution = np.linalg.norm(eigenvec[topDimension:(topDimension+bottomDimension)])
#             #print eigenvec, "\n", topContribution, bottomContribution
#             measure = topContribution - bottomContribution
#             eigClassifications.append([[measure], [eigenvec]])
#         return eigClassifications

    def getRescaledANM(self, anmComplex, anmReference, typeOfProtein, normalize=True):
        """ Return an ANM with rescaled eigenvalues and eigenvectors from anmComplex to be used 
        on the reference """       
        if self.utils.isReceptor(typeOfProtein):
            typeOfProtein = "receptor"
        else:
            typeOfProtein = 'ligand'
        # init method variables
        eigenvalueClassifications = []
        if self.utils.config.calculateZeroEigvalModes == True:
            Mcomplex = anmComplex[0].getArray().T[6:]
            eigenvalsComplex = anmComplex[0].getEigvals()[6:]
            Mreference = anmReference[0].getArray().T[6:]
            eigenvalsReference = anmReference[0].getEigvals()[6:]
        else:
            Mcomplex = anmComplex[0].getArray().T
            eigenvalsComplex = anmComplex[0].getEigvals()
            Mreference = anmReference[0].getArray().T
            eigenvalsReference = anmReference[0].getEigvals()
        numOfreferenceCoords = len(Mreference[0])
        
#### testing variables for the rescaling code below, uncomment them to test the algorithm
#         Mcomplex = np.array([[0.10482848, 0.20965697, 0.31448545, 0.41931393, 0.52414242, 0.6289709], [0.54028161, 0.47274641, 0.60781681, 0.2026056, 0.1350704, 0.2363732]])
#         eigenvalsComplex = np.array([1.2, 2.0])
#         Mreference = np.array([[100, 200, 250], [20, 30, 40]])
#         eigenvalsReference = np.array([10, 20])
#         numOfreferenceCoords = len(Mreference[0])
#         normalize=False
        
        # rescale eigenvalues
        assert len(eigenvalsComplex) == len(Mcomplex)
        assert len(eigenvalsReference) == len(Mreference)
        for i in range (0, len(eigenvalsComplex)):
            assert np.isclose(np.linalg.norm(Mcomplex[i]), 1.0)
            if typeOfProtein == "receptor":
                eigenvectorReference = Mcomplex[i][0:numOfreferenceCoords]
            else:
                eigenvectorReference = Mcomplex[i][-numOfreferenceCoords:]
            eigenvalRescaled = eigenvalsComplex[i] / np.power(np.linalg.norm(eigenvectorReference), 2)
            # normalize eigenvectorReference before appending to the collection
            if normalize:
                eigenvectorReference = Vector(eigenvectorReference)
                eigenvectorReference = eigenvectorReference.getNormed()
                eigenvectorReference = eigenvectorReference.getArray()
            eigenvalueClassifications.append([[eigenvectorReference], eigenvalRescaled])
        
        # sort from smallest to biggest eigenvalue
        eigenvalueClassifications.sort(key=lambda x: x[1], reverse=False)

        # create the mode array and the eigenvalue array in a ProDy compatible shape
        selectedModes = []
        selectedEigenVals = []
        for i in range(0, len(Mreference)):
            selectedModes.append(eigenvalueClassifications[i][0][0].tolist())
            selectedEigenVals.append(eigenvalueClassifications[i][1])
        selectedModes = np.array(selectedModes)
        selectedEigenVals = np.array(selectedEigenVals)
        # create custom ANM and return it
        anmReference[0].setEigens(selectedModes.T, selectedEigenVals)
        return anmReference
    
    def getRescaledComplexANM(self, anmComplex, reference):
        """ Return the complex ANM with rescaled eigenvalues and eigenvectors based on lambda^R from the reference.
        Also output an array with the new ordering of modes, array indices indicating the current position after lambda^R, and the 
        integer element the original index of this mode.
        
        Args:
            anmComplex: ANM of the complex
            reference: protein, a part of the complex, on its dimensionalities the lambda^R scaling is defined
            
        Returns:
            Tuple: (abmComplex with lambda^R rescaling, new ordering of modes)
        """        
        # init method variables
        eigenvalueClassifications = []
        if self.utils.config.calculateZeroEigvalModes == True:
            Mcomplex = anmComplex[0].getArray().T[6:]
            eigenvalsComplex = anmComplex[0].getEigvals()[6:]
            MComplex_trivial = anmComplex[0].getArray().T[:6]
            eigenvalsComplex_trivial = anmComplex[0].getEigvals()[:6]
        else:
            Mcomplex = anmComplex[0].getArray().T
            eigenvalsComplex = anmComplex[0].getEigvals()
        # reference coordinates in the c-alpha hessian
        numOfreferenceCoords = reference.select('calpha').numAtoms()*3
        
# #### testing variables for the rescaling code below, uncomment them to test the algorithm
#         Mcomplex = np.array([[0.10482848, 0.20965697, 0.31448545, 0.41931393, 0.52414242, 0.6289709], [0.54028161, 0.47274641, 0.60781681, 0.2026056, 0.1350704, 0.2363732]])
#         eigenvalsComplex = np.array([1.2, 2.0])
#         Mreference = np.array([[100, 200, 250], [20, 30, 40]])
#         eigenvalsReference = np.array([10, 20])
#         numOfreferenceCoords = len(Mreference[0])
        
        # rescale eigenvalues
        assert len(eigenvalsComplex) == len(Mcomplex)
        for i in range (0, len(eigenvalsComplex)):
            assert(len(Mcomplex[i]) == anmComplex[0].numAtoms()*3)
            assert np.isclose(np.linalg.norm(Mcomplex[i]), 1.0)
            if self.utils.isReceptor(reference.getTitle()):
                eigenvectorReference = Mcomplex[i][0:numOfreferenceCoords]
            else:
                eigenvectorReference = Mcomplex[i][-numOfreferenceCoords:]
                
            eigenvalRescaled = eigenvalsComplex[i] / np.power(np.linalg.norm(eigenvectorReference), 2)

            eigenvalueClassifications.append([[Mcomplex[i]], eigenvalRescaled])
        
        # get the indices of the sorting
        # Example: 
        # [[[array([50])], 2],
        # [[array([51])], 1],
        # [[array([49])], 4],
        # [[array([449])], 3]]
        #
        # Out[151]: [1, 0, 3, 2]
        #
        indicesOfSorting = [j[0] for j in sorted(enumerate(eigenvalueClassifications), key=lambda x:x[1][1])]
        
        # sort from smallest to biggest eigenvalue
        eigenvalueClassifications.sort(key=lambda x: x[1], reverse=False)

        # create the mode array and the eigenvalue array in a ProDy compatible shape
        selectedModes = []
        selectedEigenVals = []
        for i in range(0, len(eigenvalsComplex)):
            selectedModes.append(eigenvalueClassifications[i][0][0].tolist())
            selectedEigenVals.append(eigenvalueClassifications[i][1])
        selectedModes = np.array(selectedModes)
        selectedEigenVals = np.array(selectedEigenVals)
        # if trivial modes have been calculated, put them again in front
        Marray_reranked = selectedModes.T
        if self.utils.config.calculateZeroEigvalModes == True:        
            Marray_reranked = np.hstack((MComplex_trivial.T, Marray_reranked))
            selectedEigenVals = np.hstack((eigenvalsComplex_trivial, selectedEigenVals))
        # create custom ANM and return it
        anmComplex[0].setEigens(Marray_reranked, selectedEigenVals)

        return anmComplex, indicesOfSorting
        