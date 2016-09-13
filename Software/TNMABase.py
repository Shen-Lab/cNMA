'''
Created on Apr 18, 2014

@author: oliwa
'''

import os
import imp
import traceback
import sys
from prody.proteins.compare import matchTNMAChains
from Hungarian import Hungarian

class TNMABase(object):
    '''
    TNMABase is the TNMA base class
    '''


    def __init__(self):
        '''
        Constructor
        '''
        
    def importFromURI(self, uri, absl=False):
        """ Dynamic python module loading, adapted from http://stamat.wordpress.com/2013/06/30/dynamic-module-import-in-python/
        
        Args:
            uri: path to the module to be loaded
            
        Returns: a dynamically loaded python module
        
        """
        if not absl:
            uri = os.path.normpath(os.path.join(os.path.dirname(__file__), uri))
        path, fname = os.path.split(uri)
        mname, ext = os.path.splitext(fname)
            
        no_ext = os.path.join(path, mname)
             
        if os.path.exists(no_ext + '.pyc'):
            try:
                os.remove(no_ext + '.pyc')
            except:
                pass
        if os.path.exists(no_ext + '.py'):
            try:
                return imp.load_source(mname, no_ext + '.py')
            except:
                raise StandardError("StandardError occurred, the Configurations file could not be dynamically loaded")
            
    def getFullPathOfURI(self, uri):
        """ Get the full path of a file, adapted from http://stamat.wordpress.com/2013/06/30/dynamic-module-import-in-python/
        
        Args:
            uri: the provided path inclusive filename
            
        Returns:
            absolute (full) path of the file
        
        """
        return os.path.normpath(os.path.join(os.path.dirname(__file__), uri))
    
    def getOverallMatch(self, reference, mobile, subset):
        """
        Performs a matching of chains of the two elements reference and mobile returns
        the matches with all atoms as specified by subset.
        
        At first, a modified version of matchChains is called that only uses the 
        pairwise alignment matching of prody, but keeps the prody defaults of 
        minimum seqid and overlap/coverage. 
            - If the matches are a perfect bisection, this result is used 
            and returned.
            - Else, the hungarian algorithm is called to find the optimal 
            matches, and the result returned. 
        
        In case of the hungarian algorithm, the matchChains method has been 
        modified as follows the following addons:
            1. pairwise alignment is enforced (from Bio.pairwise2)
            2. pairwise alignment is the only matching algorithm, the prody 
                first choice of mapping based on residue numbers and type is 
                ignored
            3. minimum seqid and overlap criteria are set to 0.00001, the actual
                matching decision will be performed by the hungarian algorithm, 
                and pairwise alignment is only needed for the actual values of 
                seqid and overlap to create the cost matrix of the hungarian
                algorithm
                
        Remarks: prepareForHungarian needs to be set to True. Otherwise, the 
            ProDy matching sorts matched chains internally in decreasing order
            of sequence identity, but this order is sometimes not the order of 
            chains in the PDB file. 
        
            Args:
                reference: the unbound structure
                mobile: the bound structure
                subset: which matched atoms to return (calpha, bb, all ...)
            Returns:
                the overall match of chains from the given myTuple based on the
                    Bio.pairwise2 scores and possibly the hungarian algorithm
                
        """
        matches = matchTNMAChains(reference, 
                                       mobile, 
                                       prepareForHungarian = True,
                                       pwalign="True",
                                       subset=subset)
        # if the number of chains do not match, the behavior cannot be 
        # defined at this point
        assert reference.numChains() == mobile.numChains()
        
        if matches is None:
            return self.doHungarianMatching(reference, mobile, subset)
        elif not (reference.numChains() == mobile.numChains() == len(matches)):
            return self.doHungarianMatching(reference, mobile, subset)
        elif not self.isAOnetoOneMatch(matches):
            return self.doHungarianMatching(reference, mobile, subset)
        else:
            self.matches = matches
            # make overall match and return it
            noMatchYet = True
            for match in matches:
                ref_chain = match[0]
                mob_chain = match[1]
                if noMatchYet:
                    overallRefMatch = ref_chain
                    overallMobMatch = mob_chain
                    noMatchYet = False
                else: 
                    overallRefMatch += ref_chain
                    overallMobMatch += mob_chain
            if not noMatchYet:
                overallMatch = [overallRefMatch, overallMobMatch]
            else:
                overallMatch = [ref_chain, mob_chain]
            return overallMatch
            #return [matches[1][0]+matches[0][0], matches[1][1]+matches[0][1]]    
            
    def isAOnetoOneMatch(self, matches):
        """ Return False if matches does not have a one to one match for each 
        chain, else return True.
        
        It is assumed that len(matches) > 1 and that the number of matches
        equals the number of either chains.
        
        Args:
            matches: matches as returns by matchTNMAchains
        
        Returns: does matches contain a 1:1 perfect matching of the bisection 
        """
        
        baseSetUnbound = set(matches[0][0].getChids())
        baseSetBound = set(matches[0][1].getChids())
        assert len(baseSetUnbound) == 1, 'assert len(baseSetUnbound) == 1'
        assert len(baseSetBound) == 1, 'assert len(baseSetBound) == 1'
        for i in range(1, len(matches)):
            addonSetUnbound = set(matches[i][0].getChids())
            addonSetBound = set(matches[i][1].getChids())
            assert len(addonSetUnbound) == 1, 'assert len(addonSetUnbound) == 1'
            assert len(addonSetBound) == 1, 'assert len(addonSetBound) == 1'
            if len(baseSetUnbound.intersection(addonSetUnbound)) > 0:
                return False
            elif len(baseSetBound.intersection(addonSetBound)) > 0:
                return False
            elif len(baseSetUnbound.intersection(addonSetUnbound)) == 0:
                baseSetUnbound = baseSetUnbound.union(addonSetUnbound)
            elif len(baseSetBound.intersection(addonSetBound)) == 0:
                baseSetBound = baseSetBound.union(addonSetBound)
            else:
                print "**********\n\n\n set problem in isAOnetoOneMatch(...)"
                sys.exit()
        return True
            
    def doHungarianMatching(self, reference, mobile, subset):
        """ Do a chain matching with the help of the Hungarian Algorithm. 
        
        Args:
            reference: a structure (for instance protein) to be matched 
            mobile: another structure (of for instance the same protein in a different conformational state) to be matched
            subset: what atoms are considered for this matching (calpha, bb, all)
            
        Returns:
            object with the overall chain matchings  
        """
        print "Performing matching with the help of the Hungarian Algorithm."
        seqid = 0.00001
        overlap = 0.00001
        self.matches = matchTNMAChains(reference, 
                                       mobile, 
                                       prepareForHungarian = True,
                                       seqid=seqid, 
                                       overlap=overlap, 
                                       pwalign="True",
                                       subset=subset)
        hungarian = Hungarian()
        indices, matchesMatrix = hungarian.getHungarianIndices(
                                              reference.numChains(), 
                                              mobile.numChains(), 
                                              self.matches)
        noMatchYet = True
        for element in indices:
            ref_chain = (matchesMatrix[element[0]][element[1]])[0]
            mob_chain = (matchesMatrix[element[0]][element[1]])[1]
            if noMatchYet:
                overallRefMatch = ref_chain
                overallMobMatch = mob_chain
                noMatchYet = False
            else: 
                overallRefMatch += ref_chain
                overallMobMatch += mob_chain
        if not noMatchYet:
            overallMatch = [overallRefMatch, overallMobMatch]
        else:
            overallMatch = [ref_chain, mob_chain]
        return overallMatch
            
    def instantiateClassFromModule(self, moduleName, className):
        """ Instantiate the class "className" from the module "moduleName" and return it 
        
        Args:
            moduleName: the name of the module
            className: the name of the class to be instantiated
        
        Return:
            The class instantiated from the provided module
        """
        try:
            SomeClass = getattr(moduleName, className)
            obj = SomeClass()
            return obj
        except StandardError, e:
            print "StandardError occurred, the Configurations class could not be dynamically instantiated: ", e
            print traceback.format_exc()
            
    def equalityAssertionsOfComplexAndItsProteins(self, encounter, checkInterface=True):
        """ Assert/check equality of parsed proteins and their interfaces 
            
            Args:
                encounter: object with all parsed proteins
                checkInterface: if true, also perform the check on the interface
        """
        if self.utils.isReceptor(encounter.getReference().getTitle()):
            referenceSegment = "R"
            counterpartSegment = "L"
        else:
            referenceSegment = "L"
            counterpartSegment = "R"
        
        self.utils.checkEqualityOfProteins(encounter.getReference(), encounter.unboundComplexAligned.complex.select('segment "'+referenceSegment+'."'))
        self.utils.checkEqualityOfProteins(encounter.getRefChain(), encounter.getUnboundComplexAlignedChain().select('segment "'+referenceSegment+'."'))
        self.utils.checkEqualityOfProteins(encounter.getUnboundCounterpart(), encounter.unboundComplexAligned.complex.select('segment "'+counterpartSegment+'."'))
        self.utils.checkEqualityOfProteins(encounter.getUnboundCounterpartChain(), encounter.getUnboundComplexAlignedChain().select('segment "'+counterpartSegment+'."'))
        if checkInterface:
            self.utils.checkEqualityOfProteins(encounter.getRefChainInterface(), encounter.getUnboundComplexChainInterface().select('segment "'+referenceSegment+'."'))
            self.utils.checkEqualityOfProteins(encounter.getUnboundCounterpartChainInterface(), encounter.getUnboundComplexChainInterface().select('segment "'+counterpartSegment+'."'))
        
        self.utils.checkEqualityOfProteins(encounter.getMobile(), encounter.boundComplex.complex.select('segment "'+referenceSegment+'."'))
        self.utils.checkEqualityOfProteins(encounter.getMobChain(), encounter.getBoundComplexChain().select('segment "'+referenceSegment+'."'))
        self.utils.checkEqualityOfProteins(encounter.getBoundCounterpart(), encounter.boundComplex.complex.select('segment "'+counterpartSegment+'."'))
        self.utils.checkEqualityOfProteins(encounter.getBoundCounterpartChain(), encounter.getBoundComplexChain().select('segment "'+counterpartSegment+'."'))
        if checkInterface:
            self.utils.checkEqualityOfProteins(encounter.getMobChainInterface(), encounter.getBoundComplexChainInterface().select('segment "'+referenceSegment+'."'))
            self.utils.checkEqualityOfProteins(encounter.getBoundCounterpartChainInterface(), encounter.getBoundComplexChainInterface().select('segment "'+counterpartSegment+'."'))
            
    def equalityAssertionsOfComplexAndItsProteinsRoundCoords(self, encounter, checkInterface=True, roundTo=True):
        """ Assert/check equality of parsed proteins and their interfaces 
            
            Args:
                encounter: object with all parsed proteins
        """
        if self.utils.isReceptor(encounter.getReference().getTitle()):
            referenceSegment = "R"
            counterpartSegment = "L"
        else:
            referenceSegment = "L"
            counterpartSegment = "R"
        
        self.utils.checkEqualityOfProteins(encounter.getReference(), encounter.unboundComplexAligned.complex.select('segment "'+referenceSegment+'."'), roundTo)
        self.utils.checkEqualityOfProteins(encounter.getRefChain(), encounter.getUnboundComplexAlignedChain().select('segment "'+referenceSegment+'."'), roundTo)
        self.utils.checkEqualityOfProteins(encounter.getUnboundCounterpart(), encounter.unboundComplexAligned.complex.select('segment "'+counterpartSegment+'."'), roundTo)
        self.utils.checkEqualityOfProteins(encounter.getUnboundCounterpartChain(), encounter.getUnboundComplexAlignedChain().select('segment "'+counterpartSegment+'."'), roundTo)
        if checkInterface:
            self.utils.checkEqualityOfProteins(encounter.getRefChainInterface(), encounter.getUnboundComplexChainInterface().select('segment "'+referenceSegment+'."'), roundTo)
            self.utils.checkEqualityOfProteins(encounter.getUnboundCounterpartChainInterface(), encounter.getUnboundComplexChainInterface().select('segment "'+counterpartSegment+'."'), roundTo)
        
        self.utils.checkEqualityOfProteins(encounter.getMobile(), encounter.boundComplex.complex.select('segment "'+referenceSegment+'."'), roundTo)
        self.utils.checkEqualityOfProteins(encounter.getMobChain(), encounter.getBoundComplexChain().select('segment "'+referenceSegment+'."'), roundTo)
        self.utils.checkEqualityOfProteins(encounter.getBoundCounterpart(), encounter.boundComplex.complex.select('segment "'+counterpartSegment+'."'), roundTo)
        self.utils.checkEqualityOfProteins(encounter.getBoundCounterpartChain(), encounter.getBoundComplexChain().select('segment "'+counterpartSegment+'."'), roundTo)
        if checkInterface:
            self.utils.checkEqualityOfProteins(encounter.getMobChainInterface(), encounter.getBoundComplexChainInterface().select('segment "'+referenceSegment+'."'), roundTo)
            self.utils.checkEqualityOfProteins(encounter.getBoundCounterpartChainInterface(), encounter.getBoundComplexChainInterface().select('segment "'+counterpartSegment+'."'), roundTo)              