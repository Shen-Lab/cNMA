'''
Created on Jun 2, 2014

@author: oliwa
'''

import glob
import argparse
import sys
import os
import numpy as np
from scriptutils import mkdir_p, customCopytree, getFullPathOfURI, makeStringEndWith, runBashCommand
from prody.dynamics.functions import loadModel
from collections import OrderedDict
from prody.proteins.pdbfile import parsePDB, writePDB
from prody.measure.measure import calcDeformVector
from prody.dynamics.compare import calcOverlap
from prody.dynamics.mode import Vector
from prody.proteins.compare import matchTNMAChains
from Hungarian import Hungarian
from prody.dynamics.editing import extendModel, sliceModel
import matplotlib.pyplot as plt

def main():
    """ Visualizes normal modes memberships of an NMA experiment in a rectangle, where columns indicate the 
    membership (blue receptor, red ligand) and the tone of the color shows a linear relationship towards the overlap 
    towards the true deformation vector. """

    parser = argparse.ArgumentParser(description='Visualizes normal modes memberships of an NMA experiment in a rectangle, where columns indicate the membership (blue receptor, red ligand) and the tone of the color shows a linear relationship towards the overlap  towards the true deformation vector.')
    parser.add_argument('resultsPath', help='Path to the experimental resultsfolder')
    parser.add_argument('--extractNPZ', action="store_true", help='The NMA models are in a *.npz.tar.gz file. If this is not set, the input is expected in a *.npz subfolder in each experiment results folder')
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    resultsPath = makeStringEndWith(getFullPathOfURI(args.resultsPath), "/")
    folders = glob.glob(resultsPath+"*")

    if args.extractNPZ:
        extractNPZs(folders)
        
    getAllModeInformation(folders, resultsPath)
        
def getAllModeInformation(folders, resultsPath):
    fromMode = 12
    toMode = 42
    
    allModeMemberships = OrderedDict()
    allModeOverlaps = OrderedDict()
    
    for resultFolder in folders:
        # load models
        NPZfolder = glob.glob(resultFolder+"/*anms.npz")[0]
        NPZ_reference = glob.glob(NPZfolder+"/*reference_ANM.anm.npz")[0]
        NPZ_counterpart = glob.glob(NPZfolder+"/*anm_counterpart.anm.npz")[0]
        NPZ_complex = glob.glob(NPZfolder+"/*anm_complex.anm.npz")[0]
        anm_reference = loadModel(NPZ_reference)
        anm_counterpart = loadModel(NPZ_counterpart)
        anm_complex = loadModel(NPZ_complex)
        assert anm_reference.getArray().shape[0] + anm_counterpart.getArray().shape[0] == anm_complex.getArray().shape[0]
        # load resolution
        resolution = getResolution(folders)
        # load pdbs and deformation vector
        unboundComplex = parsePDB(glob.glob(resultFolder+"/pdbs/*ucomplex.pdb")[0])
        boundComplex = parsePDB(glob.glob(resultFolder+"/pdbs/*bcomplex.pdb")[0])        
        overallMatchComplex = getOverallMatch(unboundComplex, boundComplex, resolution)
        defvec = calcDeformVector(overallMatchComplex[0], overallMatchComplex[1])
        # get sliced ANM
        anm_complex_slc = getSlicedANM(unboundComplex, overallMatchComplex[0], anm_complex)
        # get mode memberships
        proteinTitle = os.path.basename(resultFolder)
        modeMemberships = getModeMemberships(anm_complex.getArray(), anm_reference.getArray().shape[0], anm_counterpart.getArray().shape[0])
        allModeMemberships[proteinTitle] = modeMemberships
        # get mode overlaps
        modeOverlaps = getModeOverlaps(anm_complex_slc[0].getArray(), defvec)
        allModeOverlaps[proteinTitle] = modeOverlaps
    
    visualizeModeMemberships(allModeMemberships, allModeOverlaps, fromMode, toMode, resultsPath)
    
def visualizeModeMemberships(allModeMemberships, allModeOverlaps, fromMode, toMode, resultsPath):
    dataMembership = None
    dataOverlap = None
    for k,v in allModeMemberships.items():
        print k, v[fromMode:toMode], allModeOverlaps[k][fromMode:toMode]
        if dataMembership is None:
            dataMembership = np.array(v[fromMode:toMode])
            dataOverlap = allModeOverlaps[k][fromMode:toMode]
        else:
            dataMembership = np.vstack((dataMembership, v[fromMode:toMode]))
            dataOverlap = np.vstack((dataOverlap, allModeOverlaps[k][fromMode:toMode]))
                                 
    print dataMembership
    print dataOverlap
    dataToPlot = combineMembershipAndOverlap(dataMembership, dataOverlap)
    dataToPlot = add1k1kVector(dataToPlot, 1.5, -1.5)
    print dataToPlot
        
    column_labels = allModeMemberships.keys()
    column_labels.append("1k1k_u")
    print column_labels
    row_labels = range(fromMode, toMode)
    fig, ax = plt.subplots()
    p = ax.pcolormesh(dataToPlot)
    cbar = fig.colorbar(p)
    cbar.set_ticks([-1.5, -0.5, 0.5, 1.5])
    cbar.set_ticklabels(['L bad', 'L good', 'R good', 'R bad'])  # put text labels on them
     
#     # put the major ticks at the middle of each cell, notice "reverse" use of dimension
    ax.set_yticks(np.arange(dataMembership.shape[0])+0.5, minor=False)
    ax.set_xticks(np.arange(dataMembership.shape[1])+0.5, minor=False)
  
    ax.set_xticklabels(row_labels, minor=False)
    ax.set_yticklabels(column_labels, minor=False)
    ax.xaxis.tick_top()
    #plt.show()
    plt.savefig(resultsPath+'/modeSpectrum.eps', bbox_inches='tight')
    plt.savefig(resultsPath+'/modeSpectrum.pdf', bbox_inches='tight')
    # close and reset the plot 
    plt.clf()
    plt.cla()
    plt.close()     
        
def combineMembershipAndOverlap(dataMembership, dataOverlap):
    assert dataMembership.shape == dataOverlap.shape
    plotData = np.empty([dataMembership.shape[0], dataMembership.shape[1]])
    for i in range(0, dataMembership.shape[0]):
        for j in range(0, dataMembership.shape[1]):
            if dataMembership[i][j] >= 0.0:
                plotData[i][j] = dataMembership[i][j] - np.abs(dataOverlap[i][j])
            else:
                plotData[i][j] = dataMembership[i][j] + np.abs(dataOverlap[i][j])
            print dataMembership[i][j], dataOverlap[i][j]
    return plotData

def add1k1kVector(plotData, receptorValue, ligandValue):
    k1kArray = np.empty(plotData.shape[1])
    for i in range(0, plotData.shape[1]):
        if i % 2 == 0:
            k1kArray[i] = receptorValue
        elif i % 2 == 1:
            k1kArray[i] = ligandValue
    plotData = np.vstack((plotData, k1kArray))
    return plotData
        
def getModeMemberships(modeArray, receptorDim, ligandDim):
    """ Returns the mode membership of a HC_0 model 
    
        Args: 
            modeArray: modes from a NMA in ProDy np array format (rows xyz of atoms, colums modes)
            receptorDim: dimension of the receptor atoms XYZs (top part of a columvector mode)
            ligandDim: dimension of the ligand atoms XYZs (bottom part of a columnvector mode)
        Returns:
            membership of modeArray (1.0 receptor, 0.0 ligand)
    """
    modeMemberships = []
    for mode in modeArray.T:
        assert np.allclose(mode[:receptorDim], 0.0) or np.allclose(mode[receptorDim:receptorDim+ligandDim], 0.0)
        if np.allclose(mode[:receptorDim], 0.0):
            modeMemberships.append(-1.5)
        else:
            modeMemberships.append(1.5)
    return modeMemberships

def getModeOverlaps(modeArray, defvec):
    """ Returns the mode overals of a modeArray towards a deformation vector  
        Args: 
            modeArray: modes from a NMA in ProDy np array format (rows xyz of atoms, colums modes)
            defvec: deformation vector
        Returns:
            list with overlaps of modes towards the deformation vector
    """
    modeOverlaps = []
    defvecArray = defvec.getArray()
    for mode in modeArray.T:
        assert mode.shape == defvecArray.shape
        modeOverlaps.append(calcOverlap(Vector(mode), defvec))
    print len(modeOverlaps)
    return modeOverlaps
    
def extractNPZs(folders):
    """ Extract the NPZs in each folder 
    
    Args:
        folders: path to the resultsfolders
        
    Result: extracted NPZs in each resultfolder
    """
    for resultFolder in folders:
        zipFileFullPath = glob.glob(resultFolder+"/*npz.tar.gz")[0]
        zipFile = os.path.basename(zipFileFullPath)
        extractDir = zipFile[0:zipFile.rfind(".tar")]
        runBashCommand("mkdir "+resultFolder+"/"+extractDir)
        runBashCommand("tar xfz "+zipFileFullPath+" -C "+resultFolder+"/"+extractDir)
        
def getResolution(folders):
    """ Get resolution (calpha, bb, all) of the experiment. It is assumed for speed purposed that all subfolders have 
    the same resolution """
    for resultFolder in folders:
        configFile = glob.glob(resultFolder+"/*.py")[0]
        whatAtomsToMatchLine = runBashCommand("grep self.whatAtomsToMatch "+configFile)
        return whatAtomsToMatchLine[whatAtomsToMatchLine.find("\"")+1:whatAtomsToMatchLine.rfind("\"")]
    
def getOverallMatch(reference, mobile, subset):
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
        return doHungarianMatching(reference, mobile, subset)
    elif not (reference.numChains() == mobile.numChains() == len(matches)):
        return doHungarianMatching(reference, mobile, subset)
    elif not isAOnetoOneMatch(matches):
        return doHungarianMatching(reference, mobile, subset)
    else:
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
        
def isAOnetoOneMatch( matches):
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

def doHungarianMatching(reference, mobile, subset):
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
    matches = matchTNMAChains(reference, 
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
                                          matches)
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

def getSlicedANM(reference, ref_chain, anm_reference):
    """ Get the sliced anm, given an already calculated anm_reference and matched chains """
    # Extend the anm_reference on all atoms
    anm_reference_extend = extendModel(anm_reference, reference.select('calpha'), reference, norm=True)
    # Then slice the anm_reference to the matched
    anm_reference_slc = sliceModel(anm_reference_extend[0], anm_reference_extend[1], ref_chain.getSelstr())
    # Normalize the slices anm
    anm_reference_slc = getNormalizedANM(anm_reference_slc)
    
    return anm_reference_slc

def normalizeM(M):
    """ Normalize a set of modes, which are the columnvectors in M. 
    
        Args:
            M: set of modes as columnvectors
            
        Returns: normalized (magnitude of each mode is 1) set of modes as columnvectors in M
    """        
    Mnormed = None
    if M.ndim == 1:
        modeVector = Vector(M)
        return modeVector.getNormed().getArray()
    else:
        for element in M.T:
            modeVector = Vector(element)
            modeNormalized = modeVector.getNormed()
            if Mnormed == None:
                Mnormed = modeNormalized.getArray()
            else:
                Mnormed = np.column_stack((Mnormed, modeNormalized.getArray()))
    return Mnormed

def getNormalizedANM(anm):
    """ Normalize the modes of the anm and return this anm object 
    
        Args:
            anm: the anm with modes calculated
            
        Returns: anm with normalized modes
    """
    M = normalizeM(anm[0].getArray())
    eigenvals = anm[0].getEigvals()
    anm[0].setEigens(M, eigenvals)
    return anm

if __name__ == '__main__':
    main()