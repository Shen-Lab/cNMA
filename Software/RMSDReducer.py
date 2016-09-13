'''
Created on Jan 24, 2014

@author: oliwa
'''
from prody.measure.measure import calcDeformVector
import numpy as np
from prody.dynamics.compare import calcOverlap
from prody.dynamics.mode import Vector
from prody.measure.transform import calcRMSD
from scipy.sparse.linalg import cg
from timeout import timeout
from timeout import TimeoutError
from collections import OrderedDict

class RMSDReducer(object):
    '''
    The RMSDReducer contains method to reduce the RMSD between proteins. 
    '''
    
    def __init__(self, utils):
        '''
        Constructor
        '''
        self.utils = utils
        
    def setupMTMforBetas(self, anm):
        """ Calculate and return the dot product of all ANM modes transposed times 
        all ANM modes."""
        M = anm.getArray()
        Mtrans = M.T
        MTM = np.dot(Mtrans, M)
        return MTM
    
    
    def calcRMSDReductions(self, anm_slc, ref_chain, mob_chain, defvec):
        """ Calculate a list of RMSD reductions based increasing number of modes, that are 
        combined in a linear combination with betas. 
        
        Args:
            anm_slc: The sliced ANM, with the corresponding entries of the  eigenvectors 
            towards the matched atoms
            ref_chain: The overall matched chain atoms from the unbound structure
            mob_chain: The overall matched chain atoms from the bound structure
            defvec: the deformation vector
            
        Returns:
            RMSDReductions: The reduction list of obtained RMSD values 
        """
        RMSDReductions = []
        overlap = []
        MTM = self.setupMTMforBetas(anm_slc[0])
        betasListWhole = []
        stepPointsReduction = self.utils.getRMSDReductionStepPoints(10, 10, anm_slc[0].numModes())
        guard = 0
        for i in stepPointsReduction:
            if self.utils.config.stopRMSDReductionAt:
                if i > self.utils.config.stopRMSDReductionAt:
                    # temporary, to speedup other calculations
                    continue
                
#             elif RMSDReductions and (RMSDReductions[-1] == 1):
#                 # we already reached a RMSD Rreduction of 1.0
#                 betasListWhole.append(betasListWhole[-1])
#                 RMSDReductions.append(RMSDReductions[-1])
#                 overlap.append(overlap[-1])
#                 print "already reached RMSD = 1 at i:", i
#                 raw_input()
#                 continue
            
            if guard < self.utils.config.guard:
                # calculate betas
                try:
                    betas = self.obtainLstSqBetas(anm_slc[0][0:i+1], defvec, MTM, i, betasListWhole, anm_slc)
                except TimeoutError:
                    print "RMSD timeout at modes", i,"using previous betas"
#                     with open("RMSDtimeoutMAX"+self.utils.config.whatAtomsToMatch+".txt", "a") as myfile:
#                         myfile.write(referenceName+" RMSD timeout at modes " +str(i)+" using previous betas\n ")
                    betas = self.getInitialGuess(betasListWhole, i)
                Tapprox = np.dot(betas[0:i+1], anm_slc[0][0:i+1].getArray().T)
                TapproxVector = Vector(Tapprox, "Tapprox")
                        
                # apply Tapprox to a copy of the unbound structure and get the reduced RMSD
                ref_chain_copy = ref_chain.copy()
                ref_chain_copy.setCoords(ref_chain_copy.getCoords() + TapproxVector.getArrayNx3())
                RMSD_after_Tapprox = calcRMSD(ref_chain_copy, mob_chain)
            
                if RMSDReductions:
                    if RMSD_after_Tapprox < RMSDReductions[-1]:
                    # store betas and RMSD reduction results
                        betasListWhole.append(betas)
                        RMSDReductions.append(RMSD_after_Tapprox)
            
                        # calc overlap
                        currentOverlap = calcOverlap(TapproxVector, defvec)
                        if np.isnan(currentOverlap) or np.isinf(currentOverlap):
                            print "overlap has a numerical problem"
                            if overlap:
                                overlap.append(overlap[-1])
                            else:
                                currentOverlap = 0
                        overlap.append(currentOverlap)
                        guard = 0
                    else:
                        print "previous RMSD lower at ", i
                        # else the previous RMSD was actually lower, the beta calculation was not successful
                        guard += 1
                        betasListWhole.append(betasListWhole[-1])
                        RMSDReductions.append(RMSDReductions[-1])
                        overlap.append(overlap[-1])
                else:
                    # else it is the first RMSD reduction run, no need to compare against previous RMSD
                    # store betas and RMSD reduction results
                    betasListWhole.append(betas)
                    RMSDReductions.append(RMSD_after_Tapprox)
        
                    # calc overlap
                    currentOverlap = calcOverlap(TapproxVector, defvec)
                    if np.isnan(currentOverlap) or np.isinf(currentOverlap):
                        print "overlap has a numerical problem"
                        currentOverlap = 0
                    overlap.append(currentOverlap)
            else:
                # else guard is >= self.utils.config.guard, and the RMSD reduction should go preconceived
                # calculate betas
                try:
                    betas = self.obtainLstSqBetas(anm_slc[0][0:i+1], defvec, MTM, i, betasListWhole, anm_slc, preconceived=True)
                except TimeoutError:
                    print "RMSD timeout at modes", i, "using previous betas"
#                     with open("RMSDtimeoutMAX"+self.utils.config.whatAtomsToMatch+".txt", "a") as myfile:
#                         myfile.write(referenceName+" RMSD timeout at modes " +str(i)+" using previous betas\n ")
                    betas = self.getInitialGuess(betasListWhole, i)
                Tapprox = np.dot(betas[0:i+1], anm_slc[0][0:i+1].getArray().T)
                TapproxVector = Vector(Tapprox, "Tapprox")
                        
                # apply Tapprox to a copy of the unbound structure and get the reduced RMSD
                ref_chain_copy = ref_chain.copy()
                ref_chain_copy.setCoords(ref_chain_copy.getCoords() + TapproxVector.getArrayNx3())
                RMSD_after_Tapprox = calcRMSD(ref_chain_copy, mob_chain)
                
                if self.utils.isLessOrEqualThen(RMSD_after_Tapprox, RMSDReductions[-1]):
                # store betas and RMSD reduction results
                    betasListWhole.append(betas)
                    RMSDReductions.append(RMSD_after_Tapprox)
        
                    # calc overlap
                    currentOverlap = calcOverlap(TapproxVector, defvec)
                    if np.isnan(currentOverlap) or np.isinf(currentOverlap):
                        print "overlap has a numerical problem"
                        if overlap:
                            overlap.append(overlap[-1])
                        else:
                            currentOverlap = 0
                    overlap.append(currentOverlap)
                else:
                    # else the previous RMSD was actually lower, the beta calculation was not successful
                    betasListWhole.append(betasListWhole[-1])
                    RMSDReductions.append(RMSDReductions[-1])
                    overlap.append(overlap[-1])
                                
        # cast objects
        overlap = np.array(overlap, dtype=np.float64)
        RMSDReductions = np.array(RMSDReductions, dtype=np.float64)
         
        return RMSDReductions, overlap, stepPointsReduction
    
    def calcRMSDReductionsReverse(self, anm_slc, ref_chain, mob_chain, defvec, referenceName, filePrefix):
        """ Calculate a list of RMSD reductions based increasing number of modes, that are 
        combined in a linear combination with betas. RMSD change from mob_chain to ref_chain 
        
        Args:
            anm_slc: The sliced ANM, with the corresponding entries of the eigenvectors 
            towards the matched atoms
            ref_chain: The overall matched chain atoms from the unbound structure
            mob_chain: The overall matched chain atoms from the bound structure
            defvec: the deformation vector
            referenceName: the name of the reference
            
        Returns:
            RMSDReductions: The reduction list of obtained RMSD values 
        """
        print "anm_slc[0].getArray(): ", anm_slc[0][0:2].getArray().shape
        RMSDReductions = []
        overlap = []
        MTM = self.setupMTMforBetas(anm_slc[0])
        betasListWhole = []
        stepPointsReduction = self.utils.getRMSDReductionStepPoints(10, 10, anm_slc[0].numModes())
        guard = 0
        for i in stepPointsReduction:
            if self.utils.config.stopRMSDReductionAt:
                if i > self.utils.config.stopRMSDReductionAt:
                    # temporary, to speedup other calculations
                    continue
                
#             elif RMSDReductions and (RMSDReductions[-1] == 1):
#                 # we already reached a RMSD Rreduction of 1.0
#                 betasListWhole.append(betasListWhole[-1])
#                 RMSDReductions.append(RMSDReductions[-1])
#                 overlap.append(overlap[-1])
#                 print "already reached RMSD = 1 at i:", i
#                 raw_input()
#                 continue
            
            if guard < self.utils.config.guard:
                # calculate betas
                try:
                    betas = self.obtainLstSqBetas(anm_slc[0][0:i+1], defvec, MTM, i, betasListWhole, anm_slc)
                except TimeoutError:
                    print "RMSD timeout at modes", i,"using previous betas"
                    with open("RMSDtimeout"+filePrefix+self.utils.config.whatAtomsToMatch+".txt", "a") as myfile:
                        myfile.write(referenceName+" RMSD timeout at modes " +str(i)+" using previous betas\n ")
                    betas = self.getInitialGuess(betasListWhole, i)
                Tapprox = np.dot(betas[0:i+1], anm_slc[0][0:i+1].getArray().T)
                TapproxVector = Vector(Tapprox, "Tapprox")
                        
                # apply Tapprox to a copy of the bound structure and get the reduced RMSD
                mob_chain_copy = mob_chain.copy()
                mob_chain_copy.setCoords(mob_chain_copy.getCoords() + TapproxVector.getArrayNx3())
                RMSD_after_Tapprox = calcRMSD(mob_chain_copy, ref_chain)
            
                if RMSDReductions:
                    if RMSD_after_Tapprox < RMSDReductions[-1]:
                    # store betas and RMSD reduction results
                        betasListWhole.append(betas)
                        RMSDReductions.append(RMSD_after_Tapprox)
            
                        # calc overlap
                        currentOverlap = calcOverlap(TapproxVector, defvec)
                        if np.isnan(currentOverlap) or np.isinf(currentOverlap):
                            print "overlap has a numerical problem"
                            if overlap:
                                overlap.append(overlap[-1])
                            else:
                                currentOverlap = 0
                        overlap.append(currentOverlap)
                        guard = 0
                    else:
                        print "previous RMSD lower at ", i
                        # else the previous RMSD was actually lower, the beta calculation was not successful
                        guard += 1
                        betasListWhole.append(betasListWhole[-1])
                        RMSDReductions.append(RMSDReductions[-1])
                        overlap.append(overlap[-1])
                else:
                    # else it is the first RMSD reduction run, no need to compare against previous RMSD
                    # store betas and RMSD reduction results
                    betasListWhole.append(betas)
                    RMSDReductions.append(RMSD_after_Tapprox)
        
                    # calc overlap
                    currentOverlap = calcOverlap(TapproxVector, defvec)
                    if np.isnan(currentOverlap) or np.isinf(currentOverlap):
                        print "overlap has a numerical problem"
                        currentOverlap = 0
                    overlap.append(currentOverlap)
            else:
                # else guard is >= self.utils.config.guard, and the RMSD reduction should go preconceived
                # calculate betas
                try:
                    betas = self.obtainLstSqBetas(anm_slc[0][0:i+1], defvec, MTM, i, betasListWhole, anm_slc, preconceived=True)
                except TimeoutError:
                    print "RMSD timeout at modes", i, "using previous betas"
                    with open("RMSDtimeout"+filePrefix+self.utils.config.whatAtomsToMatch+".txt", "a") as myfile:
                        myfile.write(referenceName+" RMSD timeout at modes " +str(i)+" using previous betas\n ")
                    betas = self.getInitialGuess(betasListWhole, i)
                Tapprox = np.dot(betas[0:i+1], anm_slc[0][0:i+1].getArray().T)
                TapproxVector = Vector(Tapprox, "Tapprox")
                        
                # apply Tapprox to a copy of the bound structure and get the reduced RMSD
                mob_chain_copy = mob_chain.copy()
                mob_chain_copy.setCoords(mob_chain_copy.getCoords() + TapproxVector.getArrayNx3())
                RMSD_after_Tapprox = calcRMSD(mob_chain_copy, ref_chain)
                
                if self.utils.isLessOrEqualThen(RMSD_after_Tapprox, RMSDReductions[-1]):
                # store betas and RMSD reduction results
                    betasListWhole.append(betas)
                    RMSDReductions.append(RMSD_after_Tapprox)
        
                    # calc overlap
                    currentOverlap = calcOverlap(TapproxVector, defvec)
                    if np.isnan(currentOverlap) or np.isinf(currentOverlap):
                        print "overlap has a numerical problem"
                        if overlap:
                            overlap.append(overlap[-1])
                        else:
                            currentOverlap = 0
                    overlap.append(currentOverlap)
                else:
                    # else the previous RMSD was actually lower, the beta calculation was not successful
                    betasListWhole.append(betasListWhole[-1])
                    RMSDReductions.append(RMSDReductions[-1])
                    overlap.append(overlap[-1])
                                
        # cast objects
        overlap = np.array(overlap, dtype=np.float64)
        RMSDReductions = np.array(RMSDReductions, dtype=np.float64)
         
        return RMSDReductions, overlap, stepPointsReduction    
    
    def calcRMSDReductionsReverseGeneral(self, Marray, ref_chain, mob_chain, defvec, referenceName, filePrefix):
        """ Calculate a list of RMSD reductions based increasing number of modes, that are 
        combined in a linear combination with betas. RMSD change from mob_chain to ref_chain 
        
        Args:
            Marray: Array of normal modes, same shape as getArray from an ANM object
            ref_chain: The overall matched chain atoms from the unbound structure
            mob_chain: The overall matched chain atoms from the bound structure
            defvec: the deformation vector
            referenceName: the name of the reference
            
        Returns:
            RMSDReductions: The reduction list of obtained RMSD values 
        """
        #print "Marray: ", Marray[0:2]
        RMSDReductions = []
        overlap = []
        numModes = Marray.shape[1]
        #MTM = self.setupMTMforBetas(anm_slc[0])
        Mtrans = Marray.T
        MTM = np.dot(Mtrans, Marray)        
        betasListWhole = []
        stepPointsReduction = self.utils.getRMSDReductionStepPoints(10, 10, numModes, initialStep=1)
        print "stepPointsReduction: ", stepPointsReduction
        guard = 0
        for i in stepPointsReduction:
            if self.utils.config.stopRMSDReductionAt:
                if i > self.utils.config.stopRMSDReductionAt:
                    # temporary, to speedup other calculations
                    continue
                
#             elif RMSDReductions and (RMSDReductions[-1] == 1):
#                 # we already reached a RMSD Rreduction of 1.0
#                 betasListWhole.append(betasListWhole[-1])
#                 RMSDReductions.append(RMSDReductions[-1])
#                 overlap.append(overlap[-1])
#                 print "already reached RMSD = 1 at i:", i
#                 raw_input()
#                 continue
            
            if guard < self.utils.config.guard:
                # calculate betas
                try:
                    betas = self.obtainLstSqBetasGeneral(Marray.T[0:i+1].T, defvec, MTM, i, betasListWhole, numModes)
                except TimeoutError:
                    print "RMSD timeout at modes", i,"using previous betas"
                    with open("RMSDtimeoutgeneral"+filePrefix+self.utils.config.whatAtomsToMatch+".txt", "a") as myfile:
                        myfile.write(referenceName+" RMSD timeout at modes " +str(i)+" using previous betas\n ")
                    betas = self.getInitialGuess(betasListWhole, i)
                Tapprox = np.dot(betas[0:i+1], Marray.T[0:i+1])
                TapproxVector = Vector(Tapprox, "Tapprox")
                        
                # apply Tapprox to a copy of the bound structure and get the reduced RMSD
                mob_chain_copy = mob_chain.copy()
                mob_chain_copy.setCoords(mob_chain_copy.getCoords() + TapproxVector.getArrayNx3())
                RMSD_after_Tapprox = calcRMSD(mob_chain_copy, ref_chain)
            
                if RMSDReductions:
                    if RMSD_after_Tapprox < RMSDReductions[-1]:
                    # store betas and RMSD reduction results
                        betasListWhole.append(betas)
                        RMSDReductions.append(RMSD_after_Tapprox)
            
                        # calc overlap
                        currentOverlap = calcOverlap(TapproxVector, defvec)
                        if np.isnan(currentOverlap) or np.isinf(currentOverlap):
                            print "overlap has a numerical problem"
                            if overlap:
                                overlap.append(overlap[-1])
                            else:
                                currentOverlap = 0
                        overlap.append(currentOverlap)
                        guard = 0
                    else:
                        print "previous RMSD lower at ", i
                        # else the previous RMSD was actually lower, the beta calculation was not successful
                        guard += 1
                        betasListWhole.append(betasListWhole[-1])
                        RMSDReductions.append(RMSDReductions[-1])
                        overlap.append(overlap[-1])
                else:
                    # else it is the first RMSD reduction run, no need to compare against previous RMSD
                    # store betas and RMSD reduction results
                    betasListWhole.append(betas)
                    RMSDReductions.append(RMSD_after_Tapprox)
        
                    # calc overlap
                    currentOverlap = calcOverlap(TapproxVector, defvec)
                    if np.isnan(currentOverlap) or np.isinf(currentOverlap):
                        print "overlap has a numerical problem"
                        currentOverlap = 0
                    overlap.append(currentOverlap)
            else:
                # else guard is >= self.utils.config.guard, and the RMSD reduction should go preconceived
                # calculate betas
                try:
                    betas = self.obtainLstSqBetasGeneral(Marray.T[0:i+1].T, defvec, MTM, i, betasListWhole, numModes, preconceived=True)
                except TimeoutError:
                    print "RMSD timeout at modes", i, "using previous betas"
                    with open("RMSDtimeoutgeneral"+filePrefix+self.utils.config.whatAtomsToMatch+".txt", "a") as myfile:
                        myfile.write(referenceName+" RMSD timeout at modes " +str(i)+" using previous betas\n ")
                    betas = self.getInitialGuess(betasListWhole, i)
                Tapprox = np.dot(betas[0:i+1], Marray.T[0:i+1])
                TapproxVector = Vector(Tapprox, "Tapprox")
                        
                # apply Tapprox to a copy of the bound structure and get the reduced RMSD
                mob_chain_copy = mob_chain.copy()
                mob_chain_copy.setCoords(mob_chain_copy.getCoords() + TapproxVector.getArrayNx3())
                RMSD_after_Tapprox = calcRMSD(mob_chain_copy, ref_chain)
                
                if self.utils.isLessOrEqualThen(RMSD_after_Tapprox, RMSDReductions[-1]):
                # store betas and RMSD reduction results
                    betasListWhole.append(betas)
                    RMSDReductions.append(RMSD_after_Tapprox)
        
                    # calc overlap
                    currentOverlap = calcOverlap(TapproxVector, defvec)
                    if np.isnan(currentOverlap) or np.isinf(currentOverlap):
                        print "overlap has a numerical problem"
                        if overlap:
                            overlap.append(overlap[-1])
                        else:
                            currentOverlap = 0
                    overlap.append(currentOverlap)
                else:
                    # else the previous RMSD was actually lower, the beta calculation was not successful
                    betasListWhole.append(betasListWhole[-1])
                    RMSDReductions.append(RMSDReductions[-1])
                    overlap.append(overlap[-1])
                                
        # cast objects
        overlap = np.array(overlap, dtype=np.float64)
        RMSDReductions = np.array(RMSDReductions, dtype=np.float64)
         
        return RMSDReductions, overlap, stepPointsReduction
    
    def calcRMSDReductionsExpandingSet(self, Marray, ref_chain, mob_chain, defvec, stepPointsReduction, referenceName, filePrefix):
        """ Calculate a list of RMSD reductions based increasing number of modes, that are 
        combined in a linear combination with betas. RMSD change from mob_chain to ref_chain 
        
        Args:
            Marray: Array of normal modes, same shape as getArray from an ANM object
            ref_chain: The overall matched chain atoms from the unbound structure
            mob_chain: The overall matched chain atoms from the bound structure
            defvec: the deformation vector
            stepPointsReduction: list of number of modes to successively calculate the RMSD reductions on 
            referenceName: the name of the reference, for output debugging purposes
            filePrefix: file prefix, for output debugging purposes
            
        Returns:
            RMSDReductions: The reduction list of obtained RMSD values 
        """
        RMSDReductions = []
        L_RMSReductions = []
        overlap = []
        numModes = Marray.shape[1]
        Mtrans = Marray.T
        MTM = np.dot(Mtrans, Marray)
        stepPointsReduction = stepPointsReduction - 1 # reduce every value by one to have the index match the range 0 to n-1
        print stepPointsReduction
        betasListWhole = [[0] * stepPointsReduction[0]]
        deformationSnapshots = OrderedDict()
        deformationSnapshots["proteinFrom"] = mob_chain.copy()

        for i in stepPointsReduction:
            if self.utils.config.stopRMSDReductionAt:
                if i > self.utils.config.stopRMSDReductionAt or i > numModes:
                    # temporary, to speedup other calculations
                    continue
            
            # calculate betas
            try:
                betas = self.obtainLstSqBetasGeneralizedExpanding(Marray.T[0:i+1].T, defvec, MTM, i, betasListWhole, numModes)
            except TimeoutError:
                print "RMSD timeout at modes", i,"using previous betas"
                with open("RMSDtimeoutgeneral"+filePrefix+self.utils.config.whatAtomsToMatch+".txt", "a") as myfile:
                    myfile.write(referenceName+" RMSD timeout at modes " +str(i)+" using previous betas\n ")
                betas = self.getInitialGuessExpanding(betasListWhole, i, numModes)
            Tapprox = np.dot(betas[0:i+1], Marray.T[0:i+1])
            TapproxVector = Vector(Tapprox, "Tapprox")
                    
            # apply Tapprox to a copy of the bound structure and get the reduced RMSD
            mob_chain_copy = mob_chain.copy()
            mob_chain_copy.setCoords(mob_chain_copy.getCoords() + TapproxVector.getArrayNx3())
            RMSD_after_Tapprox = calcRMSD(mob_chain_copy, ref_chain)
            L_RMSD_after_Tapprox = self.getL_RMS(mob_chain_copy, ref_chain, self.utils.config.investigationsOn)
            deformationSnapshots[i] = mob_chain_copy.copy()
        
            if RMSDReductions:
                if RMSD_after_Tapprox < RMSDReductions[-1]:
                # store betas and RMSD reduction results
                    betasListWhole.append(betas)
                    RMSDReductions.append(RMSD_after_Tapprox)
        
                    # calc overlap
                    currentOverlap = calcOverlap(TapproxVector, defvec)
                    if np.isnan(currentOverlap) or np.isinf(currentOverlap):
                        print "overlap has a numerical problem"
                        if overlap:
                            overlap.append(overlap[-1])
                        else:
                            currentOverlap = 0
                    overlap.append(currentOverlap)
                else:
                    print "previous RMSD lower at ", i
                    # else the previous RMSD was actually lower, the beta calculation was not successful
                    betasListWhole.append(betasListWhole[-1])
                    RMSDReductions.append(RMSDReductions[-1])
                    overlap.append(overlap[-1])
            else:
                # else it is the first RMSD reduction run, store betas and RMSD reduction results
                initial_RMSD = calcRMSD(mob_chain, ref_chain)
                if RMSD_after_Tapprox < initial_RMSD:
                    RMSDReductions.append(RMSD_after_Tapprox)
                else:
                    RMSDReductions.append(initial_RMSD)
                    print "first mode did not lower RMSD"
                betasListWhole.append(betas)       
                # calc overlap
                currentOverlap = calcOverlap(TapproxVector, defvec)
                if np.isnan(currentOverlap) or np.isinf(currentOverlap):
                    print "overlap has a numerical problem"
                    currentOverlap = 0
                overlap.append(currentOverlap)
                
            if L_RMSReductions:
                if L_RMSD_after_Tapprox < L_RMSReductions[-1]:
                    L_RMSReductions.append(L_RMSD_after_Tapprox)
                else:
                    print "previous L_RMS lower at ", i
                    # else the previous LRMS was actually lower, the beta calculation was not successful
                    L_RMSReductions.append(L_RMSReductions[-1])
            else:
                # else it is the first L_RMSD reduction run, store L_RMS reduction results
                initial_L_RMS = self.getL_RMS(mob_chain, ref_chain, self.utils.config.investigationsOn)
                if L_RMSD_after_Tapprox < initial_L_RMS:
                    L_RMSReductions.append(L_RMSD_after_Tapprox)
                else:
                    L_RMSReductions.append(initial_L_RMS)
                    print "first mode did not lower L_RMS"
                                
        # cast objects
        overlap = np.array(overlap, dtype=np.float64)
        RMSDReductions = np.array(RMSDReductions, dtype=np.float64)
        L_RMSReductions = np.array(L_RMSReductions, dtype=np.float64)
        deformationSnapshots["proteinTo"] = ref_chain.copy()
         
        return RMSDReductions, overlap, stepPointsReduction, L_RMSReductions, deformationSnapshots
    
    def getL_RMS(self, proteinFrom, proteinTo, investigationsOn):
        """ Get the L_RMS of proteinFrom and proteinTo (they need to be chain matched). 
        
        Args: 
            proteinFrom: Deformed protein
            proteinTo: Target protein (target of the deformation vector)
            investigationsON: "Complex" or "Individual"
        
        Returns:
            L_RMS of proteinFrom and proteinTo
        """
        if investigationsOn == "Complex":
            proteinFromL = proteinFrom.select('segment \"L.\"')
            proteinToL = proteinTo.select('segment \"L.\"')
            return calcRMSD(proteinFromL, proteinToL)
        else:
            # else it is an investigation on individual proteins, L_RMS does not apply, 
            # return RMSD of individual proteins instead
            return calcRMSD(proteinFrom, proteinTo)
    
    def calcRMSDReductionFromTo(self, Marray, proteinFrom, proteinTo, defvec, previousBetas, previousOverlap, previousRMSD, referenceName, filePrefix):
        """ Calculate a list of RMSD reductions based increasing number of modes, that are 
        combined in a linear combination with betas. RMSD change from mob_chain to ref_chain 
         
        Args:
            Marray: Array of normal modes, same shape as getArray from an ANM object
            proteinFrom: The overall matched chains of the protein to deform towards proteinTo
            proteinTo: The overall matched chains of the protein which is being deformed towards
            previousBetas: The previous betas, serves as part of the initial guess for the fitter
            previousOverlap: The previous overlap
            previousRMSD: The previous reduced RMSD
            defvec: the deformation vector from proteinFrom to proteinTo
            referenceName: the name of the reference, for output debugging if the RMSD fitter timeouts
            filePrefix: filePrefix, for output debugging if the RMSD fitter timeouts
             
        Returns:
            RMSDReduction, overlap, betas
        """
        Mtrans = Marray.T
        MTM = np.dot(Mtrans, Marray)
        
        if len(previousBetas) == 0:
            previousBetas = [0]
        else:
            previousBetas = previousBetas[-1]
        if len(previousOverlap) == 0:
            previousOverlap = 0
        else:
            previousOverlap = previousOverlap[-1]
        if len(previousRMSD) == 0:
            previousRMSD = calcRMSD(proteinFrom, proteinTo)
        else:
            previousRMSD = previousRMSD[-1]
 
        try:
            betas = self.obtainLstSqBetasGeneralized2(Marray, defvec, MTM)
        except TimeoutError:
            print "RMSD timeout at modes", Marray.shape[1]," using previous betas"
            with open("RMSDtimeoutgeneral"+filePrefix+self.utils.config.whatAtomsToMatch+".txt", "a") as myfile:
                myfile.write(referenceName+" RMSD timeout at modes " +str(Marray.shape[1])+" using previous betas\n ")
            betas = self.getInitialGuess(previousBetas, Marray.shape[1])
        Tapprox = np.dot(betas, Marray.T)
        TapproxVector = Vector(Tapprox, "Tapprox")
                 
        # apply Tapprox to a copy of proteinFrom and get the RMSD towards proteinTo           
        proteinFrom_copy = proteinFrom.copy()
        proteinFrom_copy.setCoords(proteinFrom_copy.getCoords() + TapproxVector.getArrayNx3())                 
        RMSD_after_Tapprox = calcRMSD(proteinFrom_copy, proteinTo)

        # RMSD comparison
        if previousRMSD:
            if np.isnan(RMSD_after_Tapprox) or np.isinf(RMSD_after_Tapprox) or previousRMSD < RMSD_after_Tapprox:
                print "RMSD_after_Tapprox has a numerical problem, maybe the two structures are already too close or the mode vectors are problematic"
                RMSD_after_Tapprox = previousRMSD
        
        # calc overlap
        currentOverlap = calcOverlap(TapproxVector, defvec)
        if np.isnan(currentOverlap) or np.isinf(currentOverlap):
            print "overlap has a numerical problem, maybe the two structures are already too close or the mode vectors are problematic"
            if previousOverlap:
                currentOverlap = previousOverlap
            else:
                currentOverlap = 0
          
        return RMSD_after_Tapprox, currentOverlap, betas
    
    def RMSDReductionFixedset(self, Marray, proteinFrom, proteinTo, defvec, referenceName, filePrefix):
        """ One shot calculation for the RMSD reduction. 
         
        Args:
            Marray: Array of normal modes, same shape as getArray from an ANM object
            proteinFrom: The overall matched chains of the protein to deform towards proteinTo
            proteinTo: The overall matched chains of the protein which is being deformed towards
            defvec: the deformation vector from proteinFrom to proteinTo
            referenceName: the name of the reference, for output debugging if the RMSD fitter timeouts
            filePrefix: filePrefix, for output debugging if the RMSD fitter timeouts
             
        Returns:
            RMSDReduction
        """
        Mtrans = Marray.T
        MTM = np.dot(Mtrans, Marray)        
 
        try:
            betas = self.obtainLstSqBetasGeneralized2(Marray, defvec, MTM)
        except TimeoutError:
            print "RMSD timeout at modes", Marray.shape[1]," using previous betas"
            with open("RMSDtimeoutgeneral"+filePrefix+self.utils.config.whatAtomsToMatch+".txt", "a") as myfile:
                myfile.write(referenceName+" RMSD timeout at modes " +str(Marray.shape[1])+" using previous betas\n ")
            betas = self.getInitialGuess([0], Marray.shape[1])
        Tapprox = np.dot(betas, Marray.T)
        TapproxVector = Vector(Tapprox, "Tapprox")
                 
        # apply Tapprox to a copy of proteinFrom and get the RMSD towards proteinTo           
        proteinFrom_copy = proteinFrom.copy()
        proteinFrom_copy.setCoords(proteinFrom_copy.getCoords() + TapproxVector.getArrayNx3())                 
        RMSD_after_Tapprox = calcRMSD(proteinFrom_copy, proteinTo)
        
        # RMSD comparison
        if np.isnan(RMSD_after_Tapprox) or np.isinf(RMSD_after_Tapprox):
            print "RMSD_after_Tapprox has a numerical problem, maybe the two structures are already too close or the mode vectors are problematic, returning original RMSD"
            RMSD_after_Tapprox = calcRMSD(proteinFrom, proteinTo)
        
        # calc overlap
        currentOverlap = calcOverlap(TapproxVector, defvec)
        if np.isnan(currentOverlap) or np.isinf(currentOverlap):
            print "overlap has a numerical problem, maybe the two structures are already too close or the mode vectors are problematic, returning overlap 0"
            currentOverlap = 0
          
        return RMSD_after_Tapprox, currentOverlap, betas
    
    @timeout()
    def obtainLstSqBetas(self, anm, defvec, MTMfull, modesToConsider, listofPreviousBetas, anmTuple, preconceived=False):
        """ Obtain betas by a scipy optimizer fitting, the formula is given in :
    
            Moal, Iain H., and Paul A. Bates. "SwarmDock and the Use of Normal 
            Modes in Protein-Protein Docking." International Journal of 
            Molecular Sciences 11, no. 10 (September 28, 2010): 3623-3648. 
            doi:10.3390/ijms11103623.
            
            Args:
                anm: the ANM with modes
                defvec: the deformationvector
                MTMfull: dot product of the full ANM matrix inverse times 
                         the ANM matrix
                modesToConsider: up to how many modes the betas should be calculated
                listofPreviousBetas: the list of previously calculated betas
                anmTuple: anm tuple as generated by Prody
                preconceived: has guard from config been reached or not
                
            Returns:
                the beta coefficents
        """
        M = anm.getArray()
        #print "first M original: ", M

        Tdefvec = defvec.getArray()
        #print "shape(Tdefvec): ", np.shape(Tdefvec)
        #print "shape(M): ", np.shape(M)
        if len(M) != len(Tdefvec):
            raise ValueError("Cannot calculate betas, len(M) != len(Tdefvec)")
        Mtrans = M.T
        MTM = MTMfull[:modesToConsider+1,:modesToConsider+1]  # use pre-calculated MTM
        maximalIter = self.utils.config.maxIterBetas
        
        if modesToConsider < 1:
           #print "original MTM, np.dot(Mtrans, Tdefvec) ", MTM, np.dot(Mtrans, Tdefvec)
           betas, status = cg(MTM, np.dot(Mtrans, Tdefvec), maxiter=maximalIter)[0:2]
           print "modesToConsider, status: ", modesToConsider, status
        elif not preconceived:
            initialGuess = self.getInitialGuess(listofPreviousBetas, modesToConsider)
            betas, status = cg(MTM, np.dot(Mtrans, Tdefvec), x0=initialGuess, maxiter=maximalIter, tol=self.utils.config.precisionBetaFitting)[0:2]
#             betas, status = lsmr(MTM, np.dot(Mtrans, Tdefvec), atol=self.utils.config.precisionBetaFitting, btol=self.utils.config.precisionBetaFitting, conlim=1000000000.0, maxiter=maximalIter)[0:2]
            print "modesToConsider, status: ", modesToConsider, status
        else:
            # how many modes could be calculated on this structure
            nonTrivialModes = (anmTuple[1].select('calpha').numAtoms()*3) - 6
            initialGuess = self.getInitialGuess(listofPreviousBetas, modesToConsider)
            if modesToConsider > (nonTrivialModes+self.utils.config.goOverdetermined):
                if np.linalg.det(MTM) == 0.0 or np.linalg.det(MTM) == -0.0:
                    print "modesToConsider, nonTrivialModes, status: ", modesToConsider, nonTrivialModes, "det(MTM) == 0, skipped"
                    return initialGuess
            betas, status = cg(MTM, np.dot(Mtrans, Tdefvec), x0=initialGuess, maxiter=maximalIter, tol=self.utils.config.precisionBetaFitting)[0:2]        
#             betas, status = lsmr(MTM, np.dot(Mtrans, Tdefvec), atol=self.utils.config.precisionBetaFitting, btol=self.utils.config.precisionBetaFitting, conlim=1000000000.0, maxiter=maximalIter)[0:2]  
            if status != 0:
                print "modesToConsider, nonTrivialModes, status: ", modesToConsider, nonTrivialModes, "status == ", status, " skipped"
                return initialGuess
            print "modesToConsider, status: ", modesToConsider, status            
        return betas
    
    @timeout()
    def obtainLstSqBetasGeneral(self, anm, defvec, MTMfull, modesToConsider, listofPreviousBetas, maxModes, preconceived=False):
        """ Obtain betas by a scipy optimizer fitting, the formula is given in :
    
            Moal, Iain H., and Paul A. Bates. "SwarmDock and the Use of Normal 
            Modes in Protein-Protein Docking." International Journal of 
            Molecular Sciences 11, no. 10 (September 28, 2010): 3623-3648. 
            doi:10.3390/ijms11103623.
            
            Args:
                anm: the ANM with modes
                defvec: the deformationvector
                MTMfull: dot product of the full ANM matrix inverse times 
                         the ANM matrix
                modesToConsider: up to how many modes the betas should be calculated
                listofPreviousBetas: the list of previously calculated betas
                maxModes: the number of modes
                preconceived: has guard from config been reached or not
                
            Returns:
                the beta coefficents
        """
        M = anm

        Tdefvec = defvec.getArray()
        #print "shape(Tdefvec): ", np.shape(Tdefvec)
        #print "shape(M): ", np.shape(M)
        if len(M) != len(Tdefvec):
            print "len(M): ", M.shape
            print "len(Tdefvec): ", len(Tdefvec)
            raise ValueError("Cannot calculate betas, len(M) != len(Tdefvec)")
        Mtrans = M.T
        MTM = MTMfull[:modesToConsider+1,:modesToConsider+1]  # use pre-calculated MTM
        maximalIter = self.utils.config.maxIterBetas
        
        if modesToConsider < 1:
           #print "original MTM, np.dot(Mtrans, Tdefvec) ", MTM, np.dot(Mtrans, Tdefvec)
           betas, status = cg(MTM, np.dot(Mtrans, Tdefvec), maxiter=maximalIter)[0:2]
           print "modesToConsider, status: ", modesToConsider, status
        elif not preconceived:
            initialGuess = self.getInitialGuess(listofPreviousBetas, modesToConsider)
            betas, status = cg(MTM, np.dot(Mtrans, Tdefvec), x0=initialGuess, maxiter=maximalIter, tol=self.utils.config.precisionBetaFitting)[0:2]
#             betas, status = lsmr(MTM, np.dot(Mtrans, Tdefvec), atol=self.utils.config.precisionBetaFitting, btol=self.utils.config.precisionBetaFitting, conlim=1000000000.0, maxiter=maximalIter)[0:2]
            print "modesToConsider, status: ", modesToConsider, status
        else:
            # how many modes could be calculated on this structure
            nonTrivialModes = maxModes #(maxModes[1].select('calpha').numAtoms()*3) - 6
            initialGuess = self.getInitialGuess(listofPreviousBetas, modesToConsider)
            if modesToConsider > (nonTrivialModes+self.utils.config.goOverdetermined):
                if np.linalg.det(MTM) == 0.0 or np.linalg.det(MTM) == -0.0:
                    print "modesToConsider, nonTrivialModes, status: ", modesToConsider, nonTrivialModes, "det(MTM) == 0, skipped"
                    return initialGuess
            betas, status = cg(MTM, np.dot(Mtrans, Tdefvec), x0=initialGuess, maxiter=maximalIter, tol=self.utils.config.precisionBetaFitting)[0:2]
#             betas, status = lsmr(MTM, np.dot(Mtrans, Tdefvec), atol=self.utils.config.precisionBetaFitting, btol=self.utils.config.precisionBetaFitting, conlim=1000000000.0, maxiter=maximalIter)[0:2]           
            if status != 0:
                print "modesToConsider, nonTrivialModes, status: ", modesToConsider, nonTrivialModes, "status == ", status, " skipped"
                return initialGuess
            print "modesToConsider, status: ", modesToConsider, status            
        return betas
    
    @timeout()
    def obtainLstSqBetasGeneralizedExpanding(self, anm, defvec, MTMfull, modesToConsider, listofPreviousBetas, maxModes, preconceived=False):
        """ Obtain betas by a scipy optimizer fitting, the formula is given in :
    
            Moal, Iain H., and Paul A. Bates. "SwarmDock and the Use of Normal 
            Modes in Protein-Protein Docking." International Journal of 
            Molecular Sciences 11, no. 10 (September 28, 2010): 3623-3648. 
            doi:10.3390/ijms11103623.
            
            Args:
                anm: the ANM with modes
                defvec: the deformationvector
                MTMfull: dot product of the full ANM matrix inverse times 
                         the ANM matrix
                modesToConsider: up to how many modes the betas should be calculated, starting from 0 to n-1
                listofPreviousBetas: the list of previously calculated betas
                maxModes: the number of modes
                preconceived: has guard from config been reached or not
                
            Returns:
                the beta coefficents
        """
        M = anm

        Tdefvec = defvec.getArray()
        if len(M) != len(Tdefvec):
            print "len(M): ", M.shape
            print "len(Tdefvec): ", len(Tdefvec)
            raise ValueError("Cannot calculate betas, len(M) != len(Tdefvec)")
        Mtrans = M.T
        MTM = MTMfull[:modesToConsider+1,:modesToConsider+1]  # use pre-calculated MTM
        maximalIter = self.utils.config.maxIterBetas
        
        if modesToConsider < 1:
           #print "original MTM, np.dot(Mtrans, Tdefvec) ", MTM, np.dot(Mtrans, Tdefvec)
           betas, status = cg(MTM, np.dot(Mtrans, Tdefvec), maxiter=maximalIter)[0:2]
           print "modesToConsider, status: ", modesToConsider, status
        elif not preconceived:
            initialGuess = self.getInitialGuessExpanding(listofPreviousBetas, modesToConsider, maxModes)
            betas, status = cg(MTM, np.dot(Mtrans, Tdefvec), x0=initialGuess, maxiter=maximalIter, tol=self.utils.config.precisionBetaFitting)[0:2]
#             betas, status = lsmr(MTM, np.dot(Mtrans, Tdefvec), atol=self.utils.config.precisionBetaFitting, btol=self.utils.config.precisionBetaFitting, conlim=1000000000.0, maxiter=maximalIter)[0:2]
            print "modesToConsider, status: ", modesToConsider, status
        else:
            # how many modes could be calculated on this structure
            nonTrivialModes = maxModes #(maxModes[1].select('calpha').numAtoms()*3) - 6
            initialGuess = self.getInitialGuessExpanding(listofPreviousBetas, modesToConsider, maxModes)
            if modesToConsider > (nonTrivialModes+self.utils.config.goOverdetermined):
                if np.linalg.det(MTM) == 0.0 or np.linalg.det(MTM) == -0.0:
                    print "modesToConsider, nonTrivialModes, status: ", modesToConsider, nonTrivialModes, "det(MTM) == 0, skipped"
                    return initialGuess
            betas, status = cg(MTM, np.dot(Mtrans, Tdefvec), x0=initialGuess, maxiter=maximalIter, tol=self.utils.config.precisionBetaFitting)[0:2]
#             betas, status = lsmr(MTM, np.dot(Mtrans, Tdefvec), atol=self.utils.config.precisionBetaFitting, btol=self.utils.config.precisionBetaFitting, conlim=1000000000.0, maxiter=maximalIter)[0:2]           
            if status != 0:
                print "modesToConsider, nonTrivialModes, status: ", modesToConsider, nonTrivialModes, "status == ", status, " skipped"
                return initialGuess
            print "modesToConsider, status: ", modesToConsider, status            
        return betas    
    
    @timeout()
    def obtainLstSqBetasGeneralized2(self, M, defvec, MTM, previousBetas=None):
        """ Obtain betas by a scipy optimizer fitting, the formula is given in :
    
            Moal, Iain H., and Paul A. Bates. "SwarmDock and the Use of Normal 
            Modes in Protein-Protein Docking." International Journal of 
            Molecular Sciences 11, no. 10 (September 28, 2010): 3623-3648. 
            doi:10.3390/ijms11103623.
            
            Args:
                M: the modes array
                defvec: the deformation vector
                MTM: dot product of the ANM matrix inverse times the ANM matrix
                previousBetas: previously calculated betas
                
            Returns:
                the beta coefficents
        """

        Tdefvec = defvec.getArray()
        if len(M) != len(Tdefvec):
            print "len(M): ", M.shape
            print "len(Tdefvec): ", len(Tdefvec)
            raise ValueError("Cannot calculate betas, len(M) != len(Tdefvec)")
        Mtrans = M.T

        # the default maxiter is too low, increase the number
        maximalIter = self.utils.config.maxIterBetas
        
        if M.shape[1] == 1:
            betas, status = cg(MTM, np.dot(Mtrans, Tdefvec), maxiter=maximalIter)[0:2]
        else:
            if previousBetas is not None:
                initialGuess = self.expandInitialGuess(previousBetas, M.shape[1])
                betas, status = cg(MTM, np.dot(Mtrans, Tdefvec), x0=initialGuess, maxiter=maximalIter, tol=self.utils.config.precisionBetaFitting)[0:2]
            else:
                betas, status = cg(MTM, np.dot(Mtrans, Tdefvec), maxiter=maximalIter, tol=self.utils.config.precisionBetaFitting)[0:2]            
        print "modesToConsider, status: ", M.shape[1], status          
        return betas
    
#     def obtainLstSqBetasByCollectivity(self, M, defvec, MTMfull, modesToConsider, listofPreviousBetas, anmTuple, preconceived=False):
#         """ Obtain betas by a scipy optimizer fitting, the formula is given in :
#     
#             Moal, Iain H., and Paul A. Bates. "SwarmDock and the Use of Normal 
#             Modes in Protein-Protein Docking." International Journal of 
#             Molecular Sciences 11, no. 10 (September 28, 2010): 3623-3648. 
#             doi:10.3390/ijms11103623.
#             
#             Args:
#                 anm: the ANM with modes
#                 defvec: the deformationvector
#                 MTMfull: dot product of the full ANM matrix inverse times 
#                          the ANM matrix
#                 modesToConsider: up to how many modes the betas should be calculated
#                 
#             Returns:
#                 the beta coefficents
#         """
#         ### old
#         ### M = anm.getArray()
#         
#         Tdefvec = defvec.getArray()
#         #print "shape(Tdefvec): ", np.shape(Tdefvec)
#         #print "shape(M): ", np.shape(M)
#         if len(M) != len(Tdefvec):
#             raise ValueError("Cannot calculate betas, len(M) != len(Tdefvec)")
#         Mtrans = M.T
#         MTM = MTMfull[:modesToConsider+1,:modesToConsider+1]  # use pre-calculated MTM
#         maximalIter = self.utils.config.maxIterBetas
#         
#         if modesToConsider < 1:
#            print "using one column"
#            betas, status = cg(MTM, np.dot(Mtrans, Tdefvec), maxiter=maximalIter)[0:2]
#         elif not preconceived:
#             initialGuess = self.getInitialGuess(listofPreviousBetas, modesToConsider)
#             betas, status = cg(MTM, np.dot(Mtrans, Tdefvec), x0=initialGuess, maxiter=maximalIter, tol=self.utils.config.precisionBetaFitting)[0:2]
#             print "modesToConsider, status: ", modesToConsider, status
#         else:
#             # how many modes could be calculated on this structure
#             nonTrivialModes = (anmTuple[1].select('calpha').numAtoms()*3) - 6
#             initialGuess = self.getInitialGuess(listofPreviousBetas, modesToConsider)
#             if modesToConsider > (nonTrivialModes+self.utils.config.goOverdetermined):
#                 if np.linalg.det(MTM) == 0.0 or np.linalg.det(MTM) == -0.0:
#                     print "modesToConsider, nonTrivialModes, status: ", modesToConsider, nonTrivialModes, "det(MTM) == 0, skipped"
#                     return initialGuess
#             betas, status = cg(MTM, np.dot(Mtrans, Tdefvec), x0=initialGuess, maxiter=maximalIter, tol=self.utils.config.precisionBetaFitting)[0:2]
#             if status != 0:
#                 print "modesToConsider, nonTrivialModes, status: ", modesToConsider, nonTrivialModes, "status == ", status, " skipped"
#                 return initialGuess
#             print "modesToConsider, status: ", modesToConsider, status            
#         return betas 

    def getInitialGuessExpanding(self, listofPreviousBetas, modesToConsider, maxModesOverall):
        """ Create an initial guess vector, padded with 0.0 values to the correct length.
        
            Args:
                listofPreviousBetas: the list of previously calculated Betas
                modesToConsider: up to how many modes are given to get the betas
                
            Returns:
                The initial guess vector for the betas, padded with 0.0 to reach the 
                correct length
        """
        initialGuess = listofPreviousBetas[-1]
        initialGuess = np.append(initialGuess, [x*0.0 for x in range(len(initialGuess), modesToConsider+1)])
        if len(initialGuess) > maxModesOverall:
            initialGuess = initialGuess[:maxModesOverall]
        return initialGuess
    
    def getInitialGuess(self, listofPreviousBetas, modesToConsider):
        """ Create an initial guess vector, padded with 0.0 values to the correct length.
        
            Args:
                listofPreviousBetas: the list of previously calculated Betas
                modesToConsider: up to how many modes are given to get the betas
                
            Returns:
                The initial guess vector for the betas, padded with 0.0 to reach the 
                correct length
        """
        initialGuess = listofPreviousBetas[-1]
        initialGuess = np.append(initialGuess, [x*0.0 for x in range(len(initialGuess), modesToConsider+1)])
        return initialGuess
    
    def expandInitialGuess(self, listofPreviousBetas, modesToConsider):
        """ Create an initial guess vector, padded with 0.0 values to the correct length.
        
            Args:
                listofPreviousBetas: the list of previously calculated Betas
                modesToConsider: up to how many modes are given to get the betas
                
            Returns:
                The initial guess vector for the betas, padded with 0.0 to reach the 
                correct length
        """
        initialGuess = listofPreviousBetas
        initialGuess = np.append(initialGuess, [x*0.0 for x in range(len(initialGuess), modesToConsider)])
        return initialGuess
    
    def calcRMSDReductionsAidedByCollectivity(self, collectivity, highestN, excludeFirstK, anm_slc, ref_chain, mob_chain):
        indicesOfHighest = self.utils.getIndiciesofHighestN(np.abs(collectivity), highestN, excludeFirstK)
        M = self.getModeArrayBasedOnIndices(anm_slc[0], excludeFirstK, indicesOfHighest)

        defvec = calcDeformVector(ref_chain, mob_chain)
        RMSDReductions = []
        overlap = []
        Mtrans = M.T
        MTM = np.dot(Mtrans, M)
        betasListWhole = []
        stepPointsReduction = self.utils.getRMSDReductionStepPoints(10, 10, anm_slc[0].numModes())
        guard = 0
        for i in stepPointsReduction:
            if self.utils.config.stopRMSDReductionAt:
                if i > self.utils.config.stopRMSDReductionAt:
                    # temporary, to speedup other calculations
                    continue
            
            if guard < self.utils.config.guard:
                # calculate betas
                ## new Mmode instead of anm_slc and then [][]
                Mmode = self.getModeArrayKeepingFirstK(M, i)
                print "Mmode: ", np.shape(Mmode)
                betas = self.obtainLstSqBetasByCollectivity(Mmode, defvec, MTM, i, betasListWhole, anm_slc)
                Tapprox = np.dot(betas[0:i+1], anm_slc[0][0:i+1].getArray().T)
                TapproxVector = Vector(Tapprox, "Tapprox")
                        
                # apply Tapprox to a copy of the unbound structure and get the reduced RMSD
                ref_chain_copy = ref_chain.copy()
                ref_chain_copy.setCoords(ref_chain_copy.getCoords() + TapproxVector.getArrayNx3())
                RMSD_after_Tapprox = calcRMSD(ref_chain_copy, mob_chain)
            
                if RMSDReductions:
                    if RMSD_after_Tapprox < RMSDReductions[-1]:
                    # store betas and RMSD reduction results
                        betasListWhole.append(betas)
                        RMSDReductions.append(RMSD_after_Tapprox)
            
                        # calc overlap
                        currentOverlap = calcOverlap(TapproxVector, defvec)
                        if np.isnan(currentOverlap) or np.isinf(currentOverlap):
                            print "overlap has a numerical problem"
                            if overlap:
                                overlap.append(overlap[-1])
                            else:
                                currentOverlap = 0
                        overlap.append(currentOverlap)
                        guard = 0
                    else:
                        print "previous RMSD lower at ", i
                        # else the previous RMSD was actually lower, the beta calculation was not successful
                        guard += 1
                        betasListWhole.append(betasListWhole[-1])
                        RMSDReductions.append(RMSDReductions[-1])
                        overlap.append(overlap[-1])
                else:
                    # else it is the first RMSD reduction run, no need to compare against previous RMSD
                    # store betas and RMSD reduction results
                    betasListWhole.append(betas)
                    RMSDReductions.append(RMSD_after_Tapprox)
        
                    # calc overlap
                    currentOverlap = calcOverlap(TapproxVector, defvec)
                    if np.isnan(currentOverlap) or np.isinf(currentOverlap):
                        print "overlap has a numerical problem"
                        currentOverlap = 0
                    overlap.append(currentOverlap)
            else:
                # else guard is >= self.utils.config.guard, and the RMSD reduction should go preconceived
                # calculate betas
                Mmode = self.getModeArrayKeepingFirstK(M, i)
                
                betas = self.obtainLstSqBetasByCollectivity(Mmode, defvec, MTM, i, betasListWhole, anm_slc, preconceived=True)
                Tapprox = np.dot(betas[0:i+1], anm_slc[0][0:i+1].getArray().T)
                TapproxVector = Vector(Tapprox, "Tapprox")
                        
                # apply Tapprox to a copy of the unbound structure and get the reduced RMSD
                ref_chain_copy = ref_chain.copy()
                ref_chain_copy.setCoords(ref_chain_copy.getCoords() + TapproxVector.getArrayNx3())
                RMSD_after_Tapprox = calcRMSD(ref_chain_copy, mob_chain)
                
                if self.utils.isLessOrEqualThen(RMSD_after_Tapprox, RMSDReductions[-1]):
                # store betas and RMSD reduction results
                    betasListWhole.append(betas)
                    RMSDReductions.append(RMSD_after_Tapprox)
        
                    # calc overlap
                    currentOverlap = calcOverlap(TapproxVector, defvec)
                    if np.isnan(currentOverlap) or np.isinf(currentOverlap):
                        print "overlap has a numerical problem"
                        if overlap:
                            overlap.append(overlap[-1])
                        else:
                            currentOverlap = 0
                    overlap.append(currentOverlap)
                else:
                    # else the previous RMSD was actually lower, the beta calculation was not successful
                    betasListWhole.append(betasListWhole[-1])
                    RMSDReductions.append(RMSDReductions[-1])
                    overlap.append(overlap[-1])
                                
        # cast objects
        overlap = np.array(overlap, dtype=np.float64)
        RMSDReductions = np.array(RMSDReductions, dtype=np.float64)
         
        return RMSDReductions, overlap, stepPointsReduction
        
    def getModeArrayBasedOnIndices(self, anm_slc, excludeFirstK, indicesOfHighest):
        """ Create an array of np.arrays with the modes specified by the indices in excludeFirstK, 
        and the following modes as given by the indices in indicesOfHighest"""
        excludeFirstK = range(0, excludeFirstK)
        M = anm_slc[excludeFirstK[0]].getArray()
        #print "initial M: ", M
        for i in range(1, len(excludeFirstK)):
            M = np.dstack((M, anm_slc[excludeFirstK[i]].getArray()))
        #    print "first ",i," M: ", M
        for j in range(0, len(indicesOfHighest)):
            M = np.dstack((M, anm_slc[indicesOfHighest[j]].getArray()))
        #    print "highe ",j," M: ", M
        return M[0]
    
    def getModeArrayKeepingFirstK(self, arr, k):
        k += 1
        k = range(0, k)
        arrCopy = arr.copy()
        if len(k) == 1:
            Mbefore = np.array(np.dstack(arrCopy)[0][0])
            M = np.zeros((len(Mbefore), 1))
            #print "M: ", M
            for i in range(0, len(Mbefore)): 
                M[i] = Mbefore[i]
            return M
        elif len(arr[0]) == len(k):
            return arr
        else:
            M = np.dstack(arrCopy)[0][0]
            #print "first M in keep first k: ", M
            for i in range(1, len(k)):
                M = np.dstack((M, np.dstack(arrCopy)[0][i]))
                #print "M in keep first "+str(i)+": ", M
            return M[0]