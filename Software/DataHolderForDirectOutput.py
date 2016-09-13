'''
Created on Apr 29, 2014

@author: oliwa
'''

class DataHolderForDirectOutput(object):
    '''
    Holds the direct results data from the NMAUnified that is to be outputted. 
    '''


    def __init__(self, protein1_A_name):
        '''
        Constructor
        Args:
            protein1_A_name: the name of the protein1_A_name
        '''
        self.protein1_A_name = protein1_A_name
        # whole references to the direct data
        RMSD_unbound_to_superposed_bound = None
        RMSDReductionsWhole = None
        overlapTApproxWhole = None
        stepPointsReductionWhole = None
        overlapArrayWhole = None
        cumulOverlapWholePrody = None 
        collectivityArrayWhole = None 
        correlationArrayWhole = None
        # interface references to the direct data
        RMSD_interface = None
        RMSDReductionsInterface = None
        overlapTApproxInterface = None
        stepPointsReductionInterface = None 
        overlapArrayInterface = None
        cumulOverlapInterfacePrody = None
        collectivityArrayInterface = None
        correlationArrayInterface = None
        # distance measures (I-rms, ligand-rms)
        counterpart_rms = None
        L_rms = None
        I_rms_before_align = None
        I_rms_after_align = None
        # L_RMS reduction
        L_RMSReductions = None
        L_RMSD_unbound_to_superposed_bound = None
        # pdbs
        reference = None
        mobile = None
        unboundCounterpart = None
        boundCountertpart = None
        unboundComplex = None
        boundComplex = None
        refChain = None
        refChainInterface = None
        mobChain = None
        mobChainInterface = None
        unboundCounterpartChain = None
        unboundCounterpartChainInterface = None
        boundCounterpartChain = None
        boundCounterpartChainInterface = None
        unboundComplexAlignedChain = None
        unboundComplexChainInterface = None
        boundComplexChain = None
        boundComplexChainInterface = None
        #overlap of Marray superset
        singleModeOverlapsFromSuperset = None
        deformationSnapshots = None
        #lambdaR for complex
        indicesOfLambdaRSorting = None
        
        
        