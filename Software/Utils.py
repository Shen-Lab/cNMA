'''
Created on Nov 6, 2013

@author: oliwa
'''

import glob
import sys
import os as os
import errno as errno
#from pylab import *
import numpy as np, numpy
import StringIO
from os.path import basename, splitext
from prody.proteins.pdbfile import parsePDBStream, writePDB, parsePDB,\
    writePDBStream
from prody.dynamics.anm import calcANM, ANM
from prody.dynamics.functions import saveModel, loadModel, writeArray
from prody.dynamics.editing import extendModel, sliceModel
from prody.dynamics.mode import Vector
from itertools import izip
from os.path import split, normpath
from collections import OrderedDict
from prody.proteins.localpdb import iterPDBFilenames
import string as string
import traceback
from prody.dynamics.compare import calcOverlap
from scipy.integrate import trapz


class Utils(object):
    """Utilities class for TNMA.
    """
    
    def __init__(self, config):
        """Constructor. """
        self.config = config
        
    def getAllPDBs(self):
        """Return a list of all *.pdb files in the folder specified by the 
        class variable self.path.

        Returns:
            a list of all *.pdb files in self.path.
        """
        return glob.glob(self.path+"*.pdb")
        
    def _splitIntoLRBU(self):
        """Split the contents of allPDBfs into receptors/ligand (un)bound lists
        into the four variables: ligandBound, ligandUnbound, receptorBound, 
        receptorUnbound.

        Args:
            allPDBs: a list of all *.pdb files.

        """
        # empty the member lists
        self.ligandBound[:] = []
        self.ligandUnbound[:] = []
        self.receptorBound[:] = []
        self.receptorUnbound[:] = []
        # populate the empty lists with entries
        for entry in self.allPDBs:
            if entry.endswith('l_b.pdb'):
                self.ligandBound.append(entry)
            elif entry.endswith('l_u.pdb'):
                self.ligandUnbound.append(entry)
            elif entry.endswith('r_b.pdb'):
                self.receptorBound.append(entry)
            elif entry.endswith('r_u.pdb'):
                self.receptorUnbound.append(entry);

    def getUnboundBoundTuples(self, name):
        """Get a list of tuples of unbound and bound structures. The argument 
        name decides if receptor or ligand tuples will be returned.

        Args:
            name: receptor or ligand.
        Returns:
            a list of tuples of matched receptor or ligand unbound and bound 
            structures.

        """
        # populate ligandBound, ligandUnbound, receptorBound, receptorUnbound
        self._splitIntoLRBU()
        assert name != ""
        if name == 'receptor':
            assert self.receptorBound != []
            assert self.receptorUnbound != []
            assert len(self.receptorBound) == len(self.receptorUnbound)
            return zip(sorted(self.receptorUnbound), sorted(self.receptorBound))
        elif name == 'ligand':
            assert self.ligandBound != []
            assert self.ligandUnbound != []
            assert len(self.ligandBound) == len(self.ligandUnbound)
            return zip(sorted(self.ligandUnbound), sorted(self.ligandBound))
        else: 
            raise Exception("name is neither \"receptor\" nor \"ligand\" but \""
                            +name+"\".")
            
#     def writeExperimentResultsToFile(self, fileName, data):
#         """Write the experiment result contents of data into the file named 
#         fileName located at /home/oliwa/workspace/TNMA1/output/<fileName>
#         """
#         f = open("/home/oliwa/workspace/TNMA1/output/"+fileName, 'a')
#         f.write(data)
#         f.close()
#         
#     def writeToFile(self, fileName, data):
#         """ Write (append) data into file and append a \n. """
#         f = open("/home/oliwa/workspace/TNMA1/inspections/"+fileName, 'a')
#         f.write(data)
#         f.write("\n")
#         f.close()
    
    def espaceStringForLatex(self, inputString):
        """ Replaces _ with \_ """
        return inputString.replace("_", "\_")
        
    def filterPDB(self, path, exclude=None):
        """Parses the PDB file in path, excluding lines starting with exclude, 
        and returns it.
        
        parsePDBStream() from prody is called on a stringIO object that has 
        all lines of the file in path excluding the lines staring with the 
        string specified by exclude (for example HETATM).
        
        Args:
            path: path to the PDB file
            exclude: string based on which lines are excluded from being 
                parsed (for example HETATM) 
        
        Returns:
            The parsed and filtered (based on exclude) PDB object 
        """
        pdbFile = open(path, 'rt')
        stringIO = StringIO.StringIO()
        for line in pdbFile:
            if exclude is not None and line.startswith(exclude):
                continue
            # if it does not start with exclude, it is a valid line
            stringIO.write(line)
        stringIO.seek(0)
        fileName = splitext(basename(path))[0]
        filteredPDBfile = parsePDBStream(stringIO, title=fileName)
        return filteredPDBfile
    
    def filterPDBLines(self, path, isReceptor, receptorLines, exclude=None):
        """Parses the PDB file in path, excluding lines starting with exclude and returns it.
        
        parsePDBStream() from prody is called on a stringIO object that has 
        all lines of the file in path excluding the lines staring with the 
        string specified by exclude (for example HETATM).
        
        Args:
            path: path to the PDB file
            exclude: string based on which lines are excluded from being 
                parsed (for example HETATM) 
        
        Returns:
            The parsed and filtered (based on exclude and excludeLinesUntil) PDB object 
        """
        pdbFile = open(path, 'rt')
        stringIO = StringIO.StringIO()
        currentLine = 1
        if isReceptor:
            for line in pdbFile:
                if exclude is not None and line.startswith(exclude):
                    currentLine += 1
                    continue
                if currentLine <= receptorLines:
                    currentLine += 1
                    continue                    
                # if it does not start with exclude and the currentLine criterion is obeyed, it is a valid line
                stringIO.write(line)
        else:
            for line in pdbFile:
                if exclude is not None and line.startswith(exclude):
                    currentLine += 1
                    continue
                if currentLine > receptorLines:
                    currentLine += 1
                    continue 
                # if it does not start with exclude and the currentLine criterion is obeyed, it is a valid line
                stringIO.write(line)
                currentLine += 1            
        stringIO.seek(0)
        fileName = splitext(basename(path))[0]
        filteredPDBfile = parsePDBStream(stringIO, title=fileName)
        return filteredPDBfile    
            
    def getPDBfromSelection(self, selection):
        """ Write a selection to a PDB stream, parse and return it. """
        fileIO = StringIO.StringIO()
        writePDBStream(fileIO, selection)
        fileIO.seek(0)
        fileIOtitle = selection.getSelstr()+" from: "+str(selection.getAtomGroup())
        output = parsePDBStream(fileIO, title=fileIOtitle)
        return output
    
    def getUnifiedPDBfromPDBs(self, selection, selection2):
        """ Write two pdb objects to a PDB stream, parse and return it. """
        fileIO = StringIO.StringIO()
        writePDBStream(fileIO, selection)
        writePDBStream(fileIO, selection2)
        fileIO.seek(0)
        fileIOtitle = selection.getTitle() + " and " + selection2.getTitle()
        fileIOtitle = fileIOtitle.replace(" ", "_")        
        output = parsePDBStream(fileIO, title=fileIOtitle)
        return output
    
    def isReceptor(self, structure):
        "If name has r, it is receptor, if l, it is ligand, else error."
        if structure[-3] == 'r':
            return True
        elif structure[-3] == 'l':
            return False
        else: raise Exception(str(structure) +" at value "+ str(structure[-3]) + "is neither \"l\" nor \"r\".")
        
    def parseCorrenspondingBoundLigandPDB(self, receptor):
        receptorList = list(receptor.getTitle())
        receptorList[-3] = 'l'
        receptorList[-1] = 'b'
        ligandName = "".join(receptorList)
        return parsePDB("benchmark/"+ligandName+".pdb")

    def parseFilteredPDB(self, path, selection=None):
        """Parse the PDB object located in the path, excluding atoms defined as "HETATM"
        and optionally selecting the "protein" part of the PDB object.
            
            Args:
                path: the location of the PDB object
                selection: selection which atoms are to be considered, currently
                it can be set to None to use all atoms (excluding HETATM) or 
                "protein", where the ProDY selection("protein") will be performed
                on the read PDB
        """
        assert (selection == None) or (selection == "protein")
        exclude = "HETATM"
        if selection == None:
            return self.filterPDB(path, exclude=exclude)
        elif selection == "protein":
            pdbWithoutHetatm = self.filterPDB(path, exclude=exclude)
            filteredPdb = self.getPDBfromSelection(pdbWithoutHetatm.select('protein'))
            filteredPdb.setTitle(pdbWithoutHetatm.getTitle())
            return filteredPdb
        else:
            raise Exception("Exception in method parseFilteredPDB.")
        
    def parseFilteredPDB2c(self, path, isReceptor, receptorLines, selection=None):
        """Parse the PDB object located in the path, excluding atoms defined as "HETATM" and 
        excluding the unbound counterpart lines and optionally selecting the "protein" part of the PDB object.
            
            Args:
                path: the location of the PDB object
                selection: selection which atoms are to be considered, currently
                it can be set to None to use all atoms (excluding HETATM) or 
                "protein", where the ProDY selection("protein") will be performed
                on the read PDB
                cCase: which of the pdb files from 2c to load
                isReceptor: is the reference receptor or not
        """
        assert (selection == None) or (selection == "protein")
        exclude = "HETATM"
        if selection == None:
            return self.filterPDBLines(path, isReceptor, receptorLines, exclude=exclude)
        elif selection == "protein":
            pdbWithoutHetatm = self.filterPDBLines(path, isReceptor, receptorLines, exclude=exclude)
            filteredPdb = self.getPDBfromSelection(pdbWithoutHetatm.select('protein'))
            filteredPdb.setTitle(pdbWithoutHetatm.getTitle())
            return filteredPdb
        else:
            raise Exception("Exception in method parseFilteredPDB.")        
        
    def parseCorrenspondingCounterpart(self, structure, wantedConformation):
        if wantedConformation == "bound":
            if self.isReceptor(structure.getTitle()):
                receptorList = list(structure.getTitle())
                receptorList[5] = 'l'
                receptorList[7] = 'b'
                ligandName = "".join(receptorList)
                return self.parseFilteredPDB(self.config.pathToBenchmark40+ligandName+".pdb", selection=self.config.filterPDB)
            else:
                ligandList = list(structure.getTitle())
                ligandList[5] = 'r'
                ligandList[7] = 'b'
                receptorName = "".join(ligandList)
                return self.parseFilteredPDB(self.config.pathToBenchmark40+receptorName+".pdb", selection=self.config.filterPDB)
        elif wantedConformation == "unbound":
            if self.isReceptor(structure.getTitle()):
                receptorList = list(structure.getTitle())
                receptorList[5] = 'l'
                receptorList[7] = 'u'
                ligandName = "".join(receptorList)
                return self.parseFilteredPDB(self.config.pathToBenchmark40+ligandName+".pdb", selection=self.config.filterPDB)
            else:
                ligandList = list(structure.getTitle())
                ligandList[5] = 'r'
                ligandList[7] = 'u'
                receptorName = "".join(ligandList)
                return self.parseFilteredPDB(self.config.pathToBenchmark40+receptorName+".pdb", selection=self.config.filterPDB)
        else:
            raise Exception(str(structure) +" or "+ str(wantedConformation) + " unknown.")
        
    def parseCorrenspondingCounterpart2c(self, structure, wantedConformation, cCase, unfilteredReferenceAtomCount):
        if wantedConformation == "bound":
            if self.isReceptor(structure.getTitle()):
                receptorList = list(structure.getTitle())
                receptorList[5] = 'l'
                receptorList[7] = 'b'
                ligandName = "".join(receptorList)
                return self.parseFilteredPDB(self.config.pathToBenchmark40+ligandName+".pdb", selection=self.config.filterPDB)
            else:
                ligandList = list(structure.getTitle())
                ligandList[5] = 'r'
                ligandList[7] = 'b'
                receptorName = "".join(ligandList)
                return self.parseFilteredPDB(self.config.pathToBenchmark40+receptorName+".pdb", selection=self.config.filterPDB)
        elif wantedConformation == "unbound":
            # make path based on cCase
            if self.isReceptor(structure.getTitle()):
                receptorList = list(structure.getTitle())
                receptorList[5] = 'l'
                receptorList[7] = 'u'
                ligandName = "".join(receptorList)
                return self.parseFilteredPDB2c(self.config.pathTo2cFiles+ligandName[0:4]+"_"+str(cCase)+".pdb", self.isReceptor(structure.getTitle()), unfilteredReferenceAtomCount, selection=self.config.filterPDB)
            else:
                ligandList = list(structure.getTitle())
                ligandList[5] = 'r'
                ligandList[7] = 'u'
                receptorName = "".join(ligandList)
                return self.parseFilteredPDB2c(self.config.pathTo2cFiles+receptorName[0:4]+"_"+str(cCase)+".pdb", self.isReceptor(structure.getTitle()), unfilteredReferenceAtomCount, selection=self.config.filterPDB)
        else:
            raise Exception(str(structure) +" or "+ str(wantedConformation) + " unknown.")
        
    def getMatchingStructure(self, bound, interface, unbound):
        """Return the indices of atom elements in bound matching the same 
        atom elements in interface. 
        
        Args:
            bound: the bound structure that has elements in interface
            interface: the interface of the structure
            unbound: the unbound structure
            
        Returns: the list of indices of the elements in bound matching interface
        """
        
        interfaceIter = interface.iterAtoms()
        boundIter = bound.iterAtoms()
        interfaceAtom = interfaceIter.next()
        structureAtom = boundIter.next()
        structureIndices = []
        structureCounter = 0
        while(True):
            try:
                if (interfaceAtom == structureAtom):
                    #print bound[structureCounter].getResname() +" "+ unbound[structureCounter].getResname()
                    structureIndices.append(structureCounter)
                    #assert bound[structureCounter].getResname() == unbound[structureCounter].getResname()
                    assert bound[structureCounter].getName() == unbound[structureCounter].getName()
                    interfaceAtom = interfaceIter.next()
                else:
                    structureAtom = boundIter.next()
                    structureCounter += 1
            except StopIteration:
                break
        assert len(structureIndices) == interface.numAtoms()
        return structureIndices
    
    def getMatchingStructureSelections(self, bound, interface, unbound):
        """Return the indices of atom elements in bound matching the same 
        atom elements in interface. 
        
        Args:
            bound: the bound structure that has elements in interface
            interface: the interface of the structure
            unbound: the unbound structure
            
        Returns: the list of indices of the elements in bound matching interface
        """
        
        interfaceIter = interface.iterAtoms()
        boundIter = bound.iterAtoms()
        interfaceAtom = interfaceIter.next()
        structureAtom = boundIter.next()
        structureIndices = []
        structureCounter = 0
        while(True):
            try:
                if (interfaceAtom == structureAtom):
                    #print bound[structureCounter].getResname() +" "+ unbound[structureCounter].getResname()
                    structureIndices.append(structureCounter)
                    #assert bound[structureCounter].getResname() == unbound[structureCounter].getResname()
                    try:
                        assert bound[structureCounter].getName() == unbound[structureCounter].getName()
                    except TypeError:
                        assert zip(bound)[structureCounter][0].getName() == zip(unbound)[structureCounter][0].getName()
                    interfaceAtom = interfaceIter.next()
                else:
                    structureAtom = boundIter.next()
                    structureCounter += 1
            except StopIteration:
                break
        assert len(structureIndices) == interface.numAtoms()
        return structureIndices
    
    def getIndicesOFfSelection(self, structureSelection, subset):
        structureSelectionIter = structureSelection.iterAtoms()
        subsetIter = subset.iterAtoms()
        structureAtom = structureSelectionIter.next()
        subsetAtom = subsetIter.next()
        structureIndices = []
        structureCounter = 0
        while(True):
            try:
                if structureAtom == subsetAtom:
                    structureIndices.append(structureCounter)
                    subsetAtom = subsetIter.next()
                else:
                    structureAtom = structureSelectionIter.next()
                    structureCounter += 1
            except StopIteration:
                break
        assert len(structureIndices) == subset.numAtoms()
        return structureIndices
    
    def printAtomInfo(self, atom1):
        print atom1, atom1.getResname(), atom1.getCoords()
    
    def assertTwoAtomsAreEqual(self, atom1, atom2, useCoords=True, useResname=False):
        assert atom1.getName() == atom2.getName()
        if useResname:
            assert atom1.getResname() == atom2.getResname()
        if useCoords:
            assert np.allclose(atom1.getCoords(), atom2.getCoords())
    
    def getInterfaceVec(self, ref_chain_interface_indices, defvec, vectorName):
        """ Obtain the interface deformation vector. """
        defvec_interface = None
        for ref_chain_index in ref_chain_interface_indices:
            if defvec_interface == None:
                defvec_interface = defvec.getArrayNx3()[ref_chain_index]
                defvec_shape = defvec_interface.shape
            elif defvec_interface.shape == defvec_shape:
                defvec_interface = np.append([defvec_interface], [defvec.getArrayNx3()[ref_chain_index]], axis=0)
            else:
                defvec_interface = np.append(defvec_interface, [defvec.getArrayNx3()[ref_chain_index]], axis=0)        
        defvec_interfaceVector = Vector(defvec_interface.flatten(), vectorName)
        return defvec_interfaceVector
    
    def getSubsetOfSelection(self, selection, listWithIndices):
        """ Gets a subset of a selection that is defined by the indices in 
        the listWithIndices variable.
        
        Args:
            selection: a selection of atoms
            listWithIndices: a list with integers denoting the positions of atoms
            in the selection when the selection is viewed as a list
            
        Return:
            a subset of the selection as defined by the indices from 
            listWithIndices
        """
        subset = None
        for element in listWithIndices:
            if (subset == None):
                subset = selection[element]
            else:
                subset += selection[element]
        return subset
        
#     def getANMPath(self, reference, numberOfModes, selstr, whatAtomsToMatch, modified="", path="/home/oliwa/anmmodels/"):
#         if modified == "":
#             return path+reference.getTitle()+"_modes"+str(numberOfModes)+"_buildOn"+selstr+"_matchedOn"+whatAtomsToMatch
#         elif modified == "extended":
#             return path+"extended/"+reference.getTitle()+"_modes"+str(numberOfModes)+"_buildOn"+selstr+"_matchedOn"+whatAtomsToMatch+"_extended"
#         elif modified == "slicedback":
#             return path+"slicedback/"+reference.getTitle()+"_modes"+str(numberOfModes)+"_buildOn"+selstr+"_matchedOn"+whatAtomsToMatch+"_slicedback"
#         else:
#             raise Exception("the variable modified is not the empty string, extended or slicedback.")
#     
#     def doesANMExist(self, reference, numberOfModes, selstr, whatAtomsToMatch, modified="", path="/home/oliwa/anmmodels/"):
#         try:
#             with open(self.getANMPath(reference, numberOfModes, selstr, whatAtomsToMatch, modified, path)+".anm.npz"):
#                 return True
#         except IOError:
#             return False
        
#     def obtainANM(self, reference, ref_chain, numberOfModes, selstr='calpha', whatAtomsToMatch='calpha', modified="", forceRebuild=False):
#         # if the base model does not exist, it needs to be created along with the
#         # extended and slicedback models
#         if forceRebuild or not self.doesANMExist(reference, numberOfModes, selstr, whatAtomsToMatch, modified):
#             # Create the anm
#             anm = calcANM(reference, n_modes=numberOfModes, selstr=selstr)
#             # First extend the anm on all atoms, then slice it back to matched
#             anm_extend = extendModel(anm[0], anm[1], reference, norm=True)
#             anm_slc = sliceModel(anm_extend[0], anm_extend[1], ref_chain.getSelstr())
#             # Save the models
#             saveModel(anm[0], 
#                       filename=self.getANMPath(reference, numberOfModes, selstr, whatAtomsToMatch), 
#                       matrices=True)
#             saveModel(anm_extend[0],
#                       filename=self.getANMPath(reference, numberOfModes, selstr, whatAtomsToMatch, modified="extended"),
#                       matrices=True
#                       )
#             saveModel(anm_slc[0],
#                       filename=self.getANMPath(reference, numberOfModes, selstr, whatAtomsToMatch, modified="slicedback"),
#                       matrices=True
#                       )
#             print "created and saved models"
#             print "anm size      : " + str(anm[0].getArray().shape)
#             print "anm_ext size  : " + str(anm_extend[0].getArray().shape)
#             print "anm_slice size: " + str(anm_slc[0].getArray().shape)
#             print "ref_chain.getSelstr()" + ref_chain.getSelstr()
#             # Save the models"
#             return anm, anm_extend, anm_slc
#         else:
#             raise Exception("Problem with capturing the selection of saved models, do not use load models from files now.")
#             try:
#                 # load models
#                 anmModel = loadModel(self.getANMPath(reference, numberOfModes, selstr, whatAtomsToMatch)+".anm.npz")
#                 anm_extendModel = loadModel(self.getANMPath(reference, numberOfModes, selstr, whatAtomsToMatch, modified="extended")+".nma.npz")
#                 anm_slcModel = loadModel(self.getANMPath(reference, numberOfModes, selstr, whatAtomsToMatch, modified="slicedback")+".nma.npz")
#                 
#                 # store models selections
#                 anmModelSelection = reference.select(selstr)
#                 anm_extendModelSelection = reference
#                 anm_slcModelSelection = ref_chain
#                 
#                 # recombine models and selections as tuples
#                 anm = (anmModel, anmModelSelection)
#                 anm_extend = (anm_extendModel, anm_extendModelSelection)
#                 anm_slc = (anm_slcModel, anm_slcModelSelection)
#                 
#                 print "loaded models"
#                 print "anm size      : " + str(anm[0].getArray().shape)
#                 print "anm_ext size  : " + str(anm_extend[0].getArray().shape)
#                 print "anm_slice size: " + str(anm_slc[0].getArray().shape)
#                 print "ref_chain.getSelstr()" + ref_chain.getSelstr()
#                 return anm, anm_extend, anm_slc
#             except IOError as e:
#                 print "Error loading ANM models from disc: "+str(e)
                
    def getSlicedInterfaceANM(self, anm_ext, ref_chain_interface):
        anm_slc_interface = sliceModel(anm_ext[0], anm_ext[1], ref_chain_interface.getSelstr())
        return anm_slc_interface

    def checkElementsForOOB(self, myarray):
        """ Inspect myarray for Out of bounds (OOB) elements. 
        If an element is NaN or inf, exit"""
        for element in myarray:
            if np.isnan(element):
                print "element: " + str(element) + " is NaN, in " + str(myarray)
                sys.exit()
            elif np.isinf(element):
                print "element: " + str(element) + " is inf, in " + str(myarray)
                sys.exit()
                
    def doTwoChainsHaveSameAtoms(self, ref_chain, mob_chain):
        assert len(ref_chain) == len(mob_chain)
        for i in range(0, len(ref_chain)):
            if ref_chain[i].getName() != mob_chain[i].getName():
                print "ref_chain[i].getName() != mob_chain[i].getName()"
                print str(ref_chain[i]) +", " + str(mob_chain[i])
                print str(ref_chain[i].getName()) +", " + str(mob_chain[i].getName())
            if ref_chain[i].getResname() != mob_chain[i].getResname():
                print "ref_chain[i].getResname() != mob_chain[i].getResname()"
                print str(ref_chain[i]) +", " + str(mob_chain[i])
                print str(ref_chain[i].getResname()) +", " + str(mob_chain[i].getResname())
                
    def getRMSDReductionStepPoints(self, everyModeUntil, stepSize, maxModes, initialStep=1):
        """ Get a list with step points, at which the modes for the RMSD reduction and overlap
            calculation will be considered. For example getRMSDReductionStepPoints(10, 10, 30) 
            returns the list: [0,1,2,3,4,5,6,7,8,9,19,29]
            
            Args:
                everyModeUntil: Starting from 0, on up to how many lowest modes should the 
                calculation be explicitely done        
                stepSize: The stepsize after everyModeUntil is reached
                maxModes: The maximum number of modes        
            """
        assert maxModes >= stepSize
        assert maxModes >= everyModeUntil
        assert initialStep >= 1
        
        if initialStep == 1:
            modesToConsider = range(0, everyModeUntil)
        else:
            modesToConsider = range(1, everyModeUntil, 2)
        modesToConsider.extend([x+stepSize for x in range(everyModeUntil-1, maxModes-stepSize, stepSize)])
        
        assert modesToConsider[-1] <= maxModes
        return modesToConsider
    
#     def normalizeArray(self, arr):
#         """ Normalizes values of a list or a np.array between 0 and 1. """
#         try: 
#             arr = arr.astype(np.float64)
#         except AttributeError:
#             arr = np.array(arr, dtype=np.float64)
#         assert not np.isnan(arr).all()
#         assert not np.isinf(arr).all()
#         result = (arr - arr.min() ) / (arr.max() - arr.min())
#         print "normalization result: ", result
#         return result

    def mkdir_p(self, path):
        """ Reimplement mkdir -p, 
        from http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python """
        try:
            os.makedirs(path)
        except OSError as exc: 
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else: raise
    
    def pairwise(self, iterable):
        """ From: http://stackoverflow.com/questions/5389507/iterating-over-every-two-elements-in-a-list """
        "s -> (s0,s1), (s2,s3), (s4, s5), ..."
        a = iter(iterable)
        return izip(a, a)
    
    def makeListofTuples(self, allPDBs):
        listOfTuples = []
        for x,y in self.pairwise(allPDBs):
            listOfTuples.append((x, y))
        return listOfTuples
    
    ### new for restructoring
    def getFileExtension(self, path):
        """ Return the extension of path (i.e. return 'txt' for path='directory/test.txt')
        From: http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
        
        Args:
            path: the full path to the file
        
        Returns:
            the extension of the file in the path
        """
        return splitext(path)[1].split('.')[-1]

    def getUnboundPDBFile(self, path):
        """ Return an ordered dictionary containing a single PDB entry. 
        Adapted from findPDBFiles in prody.proteins.localpdb
        
        Args:
            path: path to the pdb file
            
        Returns:
            Ordered dictionary with a single entry of the pdb file
        """
        pdbs = OrderedDict()
        pdb = splitext(splitext(split(path)[1])[0])[0]
        if len(pdb) == 7 and pdb.startswith('pdb'):
            pdb = pdb[3:]
        # only include unbound proteins
        if pdb[7] == 'u':
            pdbs[pdb] = path
        else:
            raise Exception("path" + str(path)+ "does not have a pdb file following the template XYZA_[r|l]_[u].pdb")
        return pdbs
    
    def findUnboundPDBFiles(self, path, case=None, **kwargs):
        """Return an ordered dictionary that maps PDB filenames to file paths.  If *case*
        is specified (``'u[pper]'`` or ``'l[ower]'``), dictionary keys (filenames)
        will be modified accordingly.  If a PDB filename has :file:`pdb` prefix,
        it will be trimmed, for example ``'1mkp'`` will be mapped to file path
        :file:`./pdb1mkp.pdb.gz`).  If a file is present with multiple extensions,
        only one of them will be returned. See also :func:`.iterPDBFilenames`. Only files following the 
        template XYZA_[r|l]_[u].pdb are considered.
        
        Adapted from findPDBFiles in prody.proteins.localpdb 
        
        Args:
            path: the path to the PDB files
        
        Returns: ordered dictionary with the PDB files
        """
    
        case = str(case).lower()
        upper = lower = False
        if case.startswith('u'):
            upper = True
        elif case.startswith('l'):
            lower = True

        pdbs = OrderedDict()
        for fn in iterPDBFilenames(path, sort=True, reverse=False, **kwargs):
            fn = normpath(fn)
            pdb = splitext(splitext(split(fn)[1])[0])[0]
            if len(pdb) == 7 and pdb.startswith('pdb'):
                pdb = pdb[3:]
            if upper:
                pdbs[pdb.upper()] = fn
            elif lower:
                pdbs[pdb.lower()] = fn
            else:
                # only include unbound proteins
                if pdb[7] == 'u':
                    pdbs[pdb] = fn
        return pdbs
    
    def getPathToBoundProtein(self, pathUnbound):
        """ Get path to the bound protein. """
        s = list(pathUnbound)
        s[-5] = 'b'
        pathBound = "".join(s)
        return pathBound
        
    def getUnboundPDBFilesFromPDBNames(self, fileWithPDBNames):
        pdbNames = np.genfromtxt(fileWithPDBNames, dtype='str')
        pdbs = OrderedDict()
        for pdbName in pdbNames:
            pdbPath = self.config.pathToBenchmark40+pdbName+".pdb"
            pdbs[pdbName] = pdbPath
        return pdbs
    
    def getResultsDictFromPath(self, path, case=None, textFile=None, onlyPDBs=False, **kwargs):
        """Return an ordered dictionary that maps PDB result folders to file paths.  
        If *case* is specified (``'u[pper]'`` or ``'l[ower]'``), dictionary keys (filenames)
        will be modified accordingly.  If a PDB filename has :file:`pdb` prefix,
        it will be trimmed, for example ``'1mkp'`` will be mapped to file path
        :file:`./pdb1mkp.pdb.gz`).  If a file is present with multiple extensions,
        only one of them will be returned. See also :func:`.iterPDBFilenames`. Only files following the 
        template XYZA_[r|l]_[u].pdb are considered.
        
        Adapted from findPDBFiles in prody.proteins.localpdb 
        
        Args:
            path: the path to the PDB result folders
        
        Returns: ordered dictionary that maps PDB result folders to file paths
        """
        case = str(case).lower()
        upper = lower = False
        if case.startswith('u'):
            upper = True
        elif case.startswith('l'):
            lower = True
 
        pdbs = OrderedDict()
        folderPaths = sorted(glob.glob(path+"*"))
        for fn in folderPaths:
            fn = normpath(fn)
            pdb = splitext(splitext(split(fn)[1])[0])[0]
            if len(pdb) == 7 and pdb.startswith('pdb'):
                pdb = pdb[3:]
            if upper:
                pdbs[pdb.upper()] = fn
            elif lower:
                pdbs[pdb.lower()] = fn
            else:
                if not onlyPDBs:
                    # only include unbound proteins
                    if pdb[7] == 'u':
                        pdbs[pdb] = fn
                # if onlyPDBs argument, include all folders 
                else:
                    pdbs[pdb] = fn
        if textFile == None:
            return pdbs
        else:
            # now filter based on textfile entries
            folderNames = np.genfromtxt(textFile, dtype='str')
            # filter based on fixed l, u positions
            if not onlyPDBs:
                pdbsFiltered = OrderedDict()
                try:
                    folderNamesSet = set(folderNames)
                # only a single element in the list
                except TypeError:
                    folderNamesSet = set([folderNames.tolist()])
                for element in pdbs.keys():
                    if element[:8] in folderNamesSet:
                        pdbsFiltered[element] = pdbs[element]
            # filter just on the pdb itself (first 4 characters of the folder)
            else:
                pdbsFiltered = OrderedDict()
                try:
                    folderNamesSet = set(folderNames)
                # only a single element in the list
                except TypeError:
                    folderNamesSet = set([folderNames.tolist()])
                for element in pdbs.keys():
                    if element[:4] in folderNamesSet:
                        pdbsFiltered[element] = pdbs[element]
            return pdbsFiltered
        
    def getOrderedDictFromFiles(self, path, mask, textFile=None):
        """ Returns an alphabetically ordered dictionary of filenames and corresponding paths, 
            specified by mask and potentially filered by entries in textFile
        
        Args:
            path: the path of the folder with files
            extension: the mask of the files to consider (for example *.txt or list*)
            textFile: if only a subset of files with extensions should be returned, the 
                      textfile holds their full names
                      
        Returns: an alphabetically ordered dictionary of filenames and corresponding paths, 
            specified by mask and potentially filered by entries in textFile
            
        
        """
        path = self.assertPathEnding(path)
        fileDict = OrderedDict()
        folderPaths = sorted(glob.glob(path+mask))
        for filePath in folderPaths:
            filePath = normpath(filePath)
            fileName = splitext(splitext(split(filePath)[1])[0])[0]
            fileDict[fileName] = filePath
        if textFile == None:
            return fileDict
        else:
            # now filter based on textfile entries
            folderNames = np.genfromtxt(textFile, dtype='str')
            fileDictFiltered = OrderedDict(zip(folderNames, [fileDict[k] for k in folderNames]))
            return fileDictFiltered 
        
    def assertPathEnding(self, path):
        """ If path does not end with "/", make it so. 
            Then, return the path. 
        """
        if path[-1] != "/":
            return path+"/"
        else:
            return path
        
    def getCounterpartTitle(self, structure):
        assert structure.getTitle()[5] == 'r' or structure.getTitle()[5] == 'l'
        assert structure.getTitle()[7] == 'u'
        if self.isReceptor(structure.getTitle()):
            receptorList = list(structure.getTitle())
            receptorList[5] = 'l'
            receptorList[7] = 'u'
            ligandName = "".join(receptorList)
            return ligandName
        else:
            ligandList = list(structure.getTitle())
            ligandList[5] = 'r'
            ligandList[7] = 'u'
            receptorName = "".join(ligandList)
            return receptorName
        
    def tryToLoadData(self, folderPath, dataName, dataFormat=None):
        """ Try to load an results array. """
        folderPath = self.assertPathEnding(folderPath)
        try:
            if dataFormat == None:
                return np.loadtxt(folderPath+dataName)
            else:
                return np.loadtxt(folderPath+dataName, dtype=dataFormat)
        except IOError, err:
            print "IOError occurred, probably there is no such file at the path: ", err
            print traceback.format_exc()
            return None
        except Exception, err:
            print "Exception occurred when trying to load data, not an IOError: ", err
            print traceback.format_exc()
            return None
        
    def tryToLoadDatabyPath(self, fullPath, dataFormat=None):
        """ Try to load an results array. """
        try:
            if dataFormat == None:
                return np.loadtxt(fullPath)
            else:
                return np.loadtxt(fullPath, dtype=dataFormat)
        except IOError, err:
            print "IOError occurred, probably there is no such file at the path: ", err
            print traceback.format_exc()
            return None
        except Exception, err:
            print "Exception occurred when trying to load data, not an IOError: ", err
            print traceback.format_exc()
            return None
        
    def isLessOrEqualThen(self, x, y):
        """ Floating point comparision with a defined threshold, and check for overflows. """
        yWithThreshold = y+1e-15
        if x < yWithThreshold and not np.isinf(yWithThreshold) and not np.isnan(yWithThreshold):
            return True
        else:
            return False
        
    def getIndiciesofHighestN(self, arr, N, excludeFirstK):
        """ Return the N highest indicies of the array. This method makes a copy
        of the initial array to not modify the original."""
        arrCopy = arr.copy()
        for i in range(0, excludeFirstK):
            arrCopy[i] = 0.0
        return np.argsort(arrCopy)[::-1][:N]
    
    def addOffset(self, inputSelstr, offset):
        """ Adds the offset to each integer value in a selstr of the form:
        \"index 1 2 3 4 67 87\"
        or 
        \"index 1 to 20 40 to 215 216 to 432\"
        
        Args:
            inputSelStr: the input selection string
            offset: the offset to be added to each integer of inputSelStr
            
        Returns:
            the offsetted inputSelStr
        """
        inputSelstr = inputSelstr.split()
        output = []
        for i, x in enumerate(inputSelstr):
            try:
                output.append((int(x) + offset))
            except ValueError:
                output.append(x)
        return " ".join(map(str, output))
    
    def makeBlockZeroMatrix(self, HR, HL):
        HRpad = np.pad(HR, ((0,0),(0,HR.shape[1])), mode='constant')
        HLpad = np.pad(HL, ((0,0),(HL.shape[1],0)), mode='constant')
        modesV = np.vstack((HRpad, HLpad))
        m1 = modesV.T[0:HR.shape[1]]
        m2 = modesV.T[HR.shape[1]:HR.shape[1]+HL.shape[1]]
        final1 =  np.vstack((m1[0], m2[0]))
        for i in range(1, HR.shape[1]):
            final2 =  np.vstack((m1[i], m2[i]))
            final1 = np.vstack((final1, final2))
        return final1.T
    
    def get1k1kModes(self, M, Mtilde, encounter, assessment, diversity=True):
        """ Get 1k modes from M and 1k from Mtilde, with optional diversity for Mtilde"""
        # make matrix with the first k modes of M, they are always taken
        firstkModesOfM = M.T[0:self.config.selectKmodes].T
        # make matrix with the k modes of Mtilde
        if not diversity:
            firstkModesOfMtilde = Mtilde.T[0:self.config.selectKmodes].T
        else:
            firstkModesOfMtilde = self.chooseModesSelectivity(encounter, assessment, self.config.selectKmodes)
        # combine alternatingly modes from firstkModesOfMtilde and firstkModesOfMtilde
        Malternating = None
        for item1, item2 in zip(firstkModesOfM.T, firstkModesOfMtilde.T):
            if Malternating is None:
                Malternating = np.vstack((item1, item2))
            else:
                Malternating = np.vstack((Malternating, item1))
                Malternating = np.vstack((Malternating, item2))
        return Malternating.T
    
    def sliceComplexModestoMatchProtein(self, Mcomplex, reference, referenceSegment):
        """ Take the modes from the complex and, if reference is a receptor, 
            slice away the bottom k elements belonging to the ligand, if the 
            reference is a ligand, slice away the top m parts belonging to the 
            receptor.
            
            Args:
                Mcomplex: The mode array from the complex
                reference: the reference to which the modes will be sliced towards
                referenceSegment: R if reference is a receptor, L if reference is a ligand
            """
        if referenceSegment == "R":
            Marray = Mcomplex[:(reference.select('calpha').numAtoms()*3)]
        else:
            Marray = Mcomplex[-(reference.select('calpha').numAtoms()*3):]
        return Marray
    
    def chooseModesSelectivity(self, encounter, assessment, k):
        M = None
        collectivities = assessment.calcCollectivityArray(encounter.accessANMs()._anm_reference_slc_original[0][0:k])
        for i in range(0, encounter.accessANMs()._anm_reference_slc[0].numModes()):
            overlaps = calcOverlap(encounter.accessANMs()._anm_reference_original[0], encounter.accessANMs()._anm_reference[0][i])
            if ((-0.5 < overlaps) & (overlaps < 0.5)).sum():
                if assessment.calcCollectivityOfSingleMode(encounter.accessANMs()._anm_reference[0][i]) > collectivities.mean():
                    if M is None:
                        M = encounter.accessANMs()._anm_reference_slc[0][i].getArray()
                    else:
                        M = np.vstack((M, encounter.accessANMs()._anm_reference_slc[0][i].getArray()))
                        if M.shape[0] == k:
                            return M.T
    
    def testBlockMatrixMembership(self, HR, HL, HC, useRelativeError=True):
        """ Returns: does 2D array HC equal to [[HR, ?], [?, HL]] """
        print "HR, HL, HC shapes: ", HR.shape, HL.shape, HC.shape
        #print "HR[:HR.shape[0], :HR.shape[1]], HL[:HL.shape[0], :HL.shape[1]] shapes: ", HR[:HR.shape[0], :HR.shape[1]].shape, HL[:HL.shape[0], :HL.shape[1]].shape
        print "HC[:HR.shape[0], :HR.shape[1]].shape: ", HC[:HR.shape[0], :HR.shape[1]].shape
        print "HC[HR.shape[0]:HR.shape[0]+HL.shape[0], HR.shape[1]:HR.shape[1]+HL.shape[1]].shape: ", HC[HR.shape[0]:HR.shape[0]+HL.shape[0], HR.shape[1]:HR.shape[1]+HL.shape[1]].shape
        print "HR HC alltrue, allclose ", np.alltrue(HR == HC[:HR.shape[0], :HR.shape[1]]), np.allclose(HR, HC[:HR.shape[0], :HR.shape[1]])
        print "HL HC alltrue, allclose ", np.alltrue(HL == HC[HR.shape[0]:HR.shape[0]+HL.shape[0], HR.shape[1]:HR.shape[1]+HL.shape[1]]), np.allclose(HL, HC[HR.shape[0]:HR.shape[0]+HL.shape[0], HR.shape[1]:HR.shape[1]+HL.shape[1]])

        #writeArray("currentHR.txt", HR, format='%f')
        #writeArray("currentHC.txt", HC[:HR.shape[0], :HR.shape[1]], format='%f')
        print "written current HR and HC"
        significantDifferences = self.returnSignificantDifferences(HR, HC[:HR.shape[0], :HR.shape[1]], useRelativeError)
        print significantDifferences
        self.checkOnlyDiagonalEntriesDiffer(significantDifferences)
        print "only diagonal entries differ"
        return significantDifferences
    
    def testHessianSubMatrices(self, anmR, anmL, anmC):
        HR = anmR[0].getHessian()
        HL = anmL[0].getHessian()
        HC = anmC[0].getHessian()
        HRtilde = HC[:HR.shape[0], :HR.shape[1]]
        HLtilde = HC[HR.shape[0]:HR.shape[0]+HL.shape[0], HR.shape[1]:HR.shape[1]+HL.shape[1]]
        print "HR, HL, HC, HRtilde, HLtilde shapes: ", HR.shape, HL.shape, HC.shape, HRtilde.shape, HLtilde.shape
        print "HR HRtilde alltrue, allclose ", np.alltrue(HR == HRtilde), np.allclose(HR, HRtilde)
        print "HL HLtilde alltrue, allclose ", np.alltrue(HL == HLtilde), np.allclose(HL, HLtilde)
        anmHRtilde = ANM("anmHRtilde")
        anmHRtilde.setHessian(HRtilde)
        print "HRtilde calcmodes"
        anmHRtilde.calcModes(n_modes=anmR[0].numModes())
        anmHLtilde = ANM("anmHLtilde")
        anmHLtilde.setHessian(HLtilde)
        print "HLtilde calcmodes"
        anmHLtilde.calcModes(n_modes=anmL[0].numModes())
        
        print "\n Equality and Floating point np.allclose test"
        print "HR HRtilde eigenvecs alltrue, allclose ", np.alltrue(anmR[0].getEigvecs() == anmHRtilde.getEigvecs()), np.allclose(anmR[0].getEigvecs(), anmHRtilde.getEigvecs())
        print "HL HLtilde eigenvecs alltrue, allclose ", np.alltrue(anmL[0].getEigvecs() == anmHLtilde.getEigvecs()), np.allclose(anmL[0].getEigvecs(), anmHLtilde.getEigvecs())
        print "HR HRtilde eigenvals alltrue, allclose ", np.alltrue(anmR[0].getEigvals() == anmHRtilde.getEigvals()), np.allclose(anmR[0].getEigvals(), anmHRtilde.getEigvals())
        print "HL HLtilde eigenvals alltrue, allclose ", np.alltrue(anmL[0].getEigvals() == anmHLtilde.getEigvals()), np.allclose(anmL[0].getEigvals(), anmHLtilde.getEigvals())        
        #print "HR      eigenvecs: ", anmR[0].getEigvecs()[0]
        #print "HRtilde eigenvecs: ", anmHRtilde.getEigvecs()[0]
        #print "HL      eigenvecs: ", anmL[0].getEigvecs()
        #print "HLtilde eigenvecs: ", anmHLtilde.getEigvecs()
        
    def printEqualityOfTwoANMs(self, anmR, anmHRtilde):
        print "\nEquality and Floating point np.allclose test"
        print "H Htilde eigenvecs alltrue, allclose ", np.alltrue(anmR[0].getEigvecs() == anmHRtilde.getEigvecs()), np.allclose(anmR[0].getEigvecs(), anmHRtilde.getEigvecs())
        print "H Htilde eigenvals alltrue, allclose ", np.alltrue(anmR[0].getEigvals() == anmHRtilde.getEigvals()), np.allclose(anmR[0].getEigvals(), anmHRtilde.getEigvals())
        
    def returnSignificantDifferences(self, Mknown, Mdiffering, useRelativeError):
        epsilon = 0.00000001
        significantDifferences = []
        assert Mknown.shape == Mdiffering.shape
        deltaMatrix = Mknown-Mdiffering
        for i in range(0, deltaMatrix.shape[0]):
            for j in range(0, deltaMatrix.shape[1]):
                if useRelativeError:
                    if not np.isclose(deltaMatrix[i][j], 0.0):
                        absolutError = np.abs(deltaMatrix[i][j])
                        relativeError = absolutError / np.abs(Mknown[i][j])
                        if relativeError > epsilon:
                            significantDifferences.append([relativeError, i, j])
                else:
                    #if not np.isclose(deltaMatrix[i][j], 0.0, atol=0.0000000001, rtol=0.0000000001):
                    if deltaMatrix[i][j] != 0.0:
                        absolutError = np.abs(deltaMatrix[i][j])
                        significantDifferences.append([absolutError, i, j])
        print significantDifferences
        return significantDifferences
    
    def significantDifferencesToPymolResiduesString(self, significantDifferences, offset):
        atomsAffected = []
        for element in significantDifferences:
            atomsAffected.append(element[1]/3)
        atomsAffected = atomsAffected + offset 
        atomsAffected = sorted(list(set(atomsAffected)))
        pymolOutput = "select resi "+"+".join(map(str, atomsAffected))
        return pymolOutput
    
    def whichPatternsAreAffectedbySignificantDifferences(self, significantDifferences):
        dict = OrderedDict()
        for element in significantDifferences:
            atom = element[1]/3
            if not dict.get(atom):
                dict[atom] = [[element[1], element[2]]]
            else:
                dict[atom].append([[element[1], element[2]]])
        print "total number of shapes (affected atoms) : ", len(dict)
        print "the following atoms are in a shape not equal to length 9: "
        for k, v in dict.items():
            if len(v) != 9:
                print v
        print "end of pattern output"
    
    def checkOnlyDiagonalEntriesDiffer(self, results):
        for element in results:
            #print "abs(", element[1], " - ", element[2], ")", " = ", "abs(", element[1] - element[2], ")" 
            assert abs(element[1] - element[2]) <= 2
            #print element[1], " / 3 = ", (element[1] / 3) , " = ", element[2] , " / 3 = ", (element[2] / 3)
            assert (element[1] / 3) == (element[2] / 3)    
    
    def checkEqualityOfProteins(self, pro1, pro2, roundTo=None):
        assert pro1.numAtoms() == pro2.numAtoms()
        for item1, item2 in zip(pro1, pro2):
            assert item1.getResname() == item2.getResname()
            if roundTo == None:
                assert  np.alltrue(item1.getCoords() == item2.getCoords())
            else:
                item1roundedChoords = [round(x, 3) for x in item1.getCoords().tolist()] 
                item2roundedChoords = [round(x, 3) for x in item2.getCoords().tolist()] 
                assert  np.alltrue(item1roundedChoords == item2roundedChoords)
            assert item1.getName() == item2.getName()
            
    def printCoordsOfProteins(self, pro1, pro2, pro3):
#         assert pro1.numAtoms() == pro2.numAtoms()
#         assert pro2.numAtoms() == pro3.numAtoms()
        for item1, item2, item3 in zip(pro1, pro2, pro3):
            assert item1.getResname() == item2.getResname()
            assert item2.getResname() == item3.getResname()
            print item1.getCoords()
            print item2.getCoords()
            print item3.getCoords()
            print "\n"
            assert item1.getName() == item2.getName()
            assert item2.getName() == item3.getName()
            
    def returnEqualityOfProteins(self, pro1, pro2):
        for item1, item2 in zip(pro1, pro2):
            if item1.getResname() != item2.getResname():
                return False
            if np.alltrue(item1.getCoords() != item2.getCoords()):
                return False
            if item1.getName() != item2.getName():
                return False
        return True
    
    def returnEqualityOfProteinsWithoutCoordinates(self, pro1, pro2):
        if pro1.numAtoms() != pro2.numAtoms():
            return False
        for item1, item2 in zip(pro1, pro2):
            if item1.getResname() != item2.getResname():
                return False
            if item1.getName() != item2.getName():
                return False
        return True
    
    def getFilledChidProtein(self, pro1, pro2):
        """ Return an array of chids with every empty chid value replaced 
        by an uppercase letter not used in protein1 or protein2.  """
        # create an array of unused chids
        pro1ChidsSet = set(pro1.getChids())
        pro2ChidsSet = set(pro2.getChids())
        currentlyUsedChidSet = pro1ChidsSet.union(pro2ChidsSet)
        asciiSet = set(string.ascii_uppercase)
        freeChids = np.array(list(asciiSet.difference(currentlyUsedChidSet)))
        assert len(freeChids) > 0
        # fill in the chid
        chids = pro1.getChids()
        for i in range(0, len(chids)):
            if chids[i] == " ":
                chids[i] = freeChids[0]
        pro1.setChids(chids)
        return pro1
    
    def fill3DArrayWithValue(self, arr, val, pos):
        """ Fills the 3D array arr with the value val on each pos """
        for i in range(0, len(arr)):
            if i % 3 == pos:
                arr[i] = val
        return arr
    
    def createRx(self, coords):
        """ Create the rotation x vector, see PH'P formula from NMA book. """
        for i in range(0, len(coords)):
            oldX = coords[i][0]
            oldY = coords[i][1]
            oldZ = coords[i][2]
            coords[i][0] = 0.0
            coords[i][1] = -oldZ
            coords[i][2] = oldY
        return coords.flatten()
    
    def createRy(self, coords):
        """ Create the rotation y vector, see PH'P formula from NMA book. """
        for i in range(0, len(coords)):
            oldX = coords[i][0]
            oldY = coords[i][1]
            oldZ = coords[i][2]
            coords[i][0] = oldZ
            coords[i][1] = 0.0
            coords[i][2] = -oldX
        return coords.flatten()
    
    def createRz(self, coords):
        """ Create the rotation z vector, see PH'P formula from NMA book. """
        for i in range(0, len(coords)):
            oldX = coords[i][0]
            oldY = coords[i][1]
            oldZ = coords[i][2]
            coords[i][0] = -oldY
            coords[i][1] = oldX
            coords[i][2] = 0.0
        return coords.flatten()
    

    def independent_columns(self, A, tol = 1e-05):
        """
        Return an array composed of independent columns of A.
    
        Note the answer may not be unique; this function returns one of many
        possible answers.
    
        http://stackoverflow.com/q/13312498/190597 (user1812712)
        http://math.stackexchange.com/a/199132/1140 (Gerry Myerson)
        http://mail.scipy.org/pipermail/numpy-discussion/2008-November/038705.html
            (Anne Archibald)
    
        >>> A = np.array([(2,4,1,3),(-1,-2,1,0),(0,0,2,2),(3,6,2,5)])
        >>> independent_columns(A)
        np.array([[1, 4],
                  [2, 5],
                  [3, 6]])
        """
        Q, R = np.linalg.qr(A)
        independent = np.where(np.abs(R.diagonal()) > tol)[0]
        return A[:, independent]
    
    def matrixrank(self, A,tol=1e-8):
        """
        http://mail.scipy.org/pipermail/numpy-discussion/2008-February/031218.html
        """
        s = np.linalg.svd(A,compute_uv=0)
        return sum( np.where( s>tol, 1, 0 ) )
    
    def countOrthogonalColumns(self, M):
        tol = 0.0001
        M = M.T
        count = 0
        for i in range(0, M.shape[0]):
            for j in range (i+1, M.shape[0]):
                if np.abs(calcOverlap(Vector(M[i]), Vector(M[j]))) < tol:
                    count += 1
        return count
    
    def file_len(self, fname):
        """ From: http://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python"""
        with open(fname) as f:
            for i, l in enumerate(f):
                pass
        return i + 1
    
    def doesAtomInXatIExistInY(self, X, Y, i):
        atomCoords = X[i].getCoords()
        xString = str(atomCoords[0])
        yString = str(atomCoords[1])
        zString = str(atomCoords[2])
        return Y.select('x == '+xString+' and y == '+yString+' and z == '+zString) 
    
#     def doesAtomExistInY(self, atom, Y):
#         """ Based on the index, does atom exist in Y. Only usable if Y is a superset of where Y is from.
#         Useful to check if an atom from a selection is in the protein that the selection is from. """
#         resNum = str(atom.getResnum())
#         result = Y.select('resnum '+resNum)
#         return result    
        
    def doesAtomExistInY(self, atom, Y):
        """ Based on the index, does atom exist in Y. Only usable if Y is a superset of where Y is from.
        Useful to check if an atom from a selection is in the protein that the selection is from. """
        indexNum = str(atom.getIndex())
        result = Y.select('index '+indexNum)
        if result:
            assert result.getResnums()[0] == atom.getResnum()
            assert np.allclose(result.getCoords(), atom.getCoords())
        return result
    
    def normalized(self, a, axis=-1, order=2):
        """ http://stackoverflow.com/questions/21030391/how-to-normalize-array-numpy """
        l2 = np.atleast_1d(np.linalg.norm(a, order, axis))
        l2[l2==0] = 1
        return a / np.expand_dims(l2, axis)
    
    def returnMaxOfLastElementsInArrays(self, *arg):
        """ Find the Max from the last elements of the Arrays, where the Arrays can be of None type """
        maxSoFar = 0
        for arr in arg:
            if arr is not None:
                maxSoFar = max(maxSoFar, arr[-1])
        return maxSoFar
    
    def returnLengthOfLongestArray(self, *arg):
        """ Find the length of the longest array, where the Arrays can be of None type """
        maxSoFar = 0
        for arr in arg:
            if arr is not None:
                maxSoFar = max(maxSoFar, len(arr))
        return maxSoFar
    
    def getModeOverlaps(self, modeArray, defvec):
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
        return modeOverlaps
    
    def getStepsize(self, arr, startIndex, stopIndex):
        """ Get the stepsize of the array between startIndex and stopIndex. 
              Assert that this stepsize is obeyed between the indices. 
        
          Args:
            arr: array
            startIndex: the index from which the stepsize should be investigated
            stopIndex: the index until which the stepsize should be investigated
        
          Returns: stepsize
        
        """
        step = arr[startIndex+1] - arr[startIndex]
        for i in range(startIndex+1, stopIndex):
            assert arr[i] - arr[i-1] == step
        return step
    
    def getSmallestCommonValue(self, arr1, arr2):
        """ Returns the value of the smallest common element of two ascendingly sorted np.arrays arr1 and arr2 
        
            Args:
                arr1: first array 
                arr2: second array
                
            Returns: value of the smallest common element of arr1 and arr2
        """
        return min(set(arr1).intersection(set(arr2)))
    
    def getLargestCommonValue(self, arr1, arr2):
        """ Returns the value of the largest common element of two ascendingly sorted np.arrays arr1 and arr2 
        
            Args:
                arr1: first array 
                arr2: second array
                
            Returns: value of the largest common element of arr1 and arr2
        """
        return max(set(arr1).intersection(set(arr2)))
    
    def getQualityMeasureBetweenCurves(self, results, resultsStepPoints, optionalResults, optionalResultsStepPoints, endModes_measure=20, startModes_measure_large=30):
        """ Returns a quality measure given by the differences of the areas under the curves of results and optionalResults.
        The formula is measure = AUC(results) - AUC(optionalResults). 
        Two measures are returned, the measure from the first identical number modes until k_1 (default 20) modes, and the measure_large from
        30 modes on until the last identical mode number. The measure stands for the performance of the results model within the first
        20 modes, the measure_large for the performance of the results model beyond k_2 (default 30) modes. 
        
        A positive measure means that overall its AUC is bigger than the AUC of the compared model, which is desired.  
        
        The different x axis step sizes are taken into consideration.

        Args: 
            results: y values of the results
            resultsStepPoints: step points of results
            optionalResults: y values of optionalResults (such as canonical Results)
            optionalResultsStepPoints: step points of optionalResults
            endModes_measure: Until how many modes to calculate the measure
            startModes_measure_large: From how many modes on to calculate measure_large
            
        Returns: measure, measure_large, maxModes_used
        """
        # investigate only the first k modes or the whole curve
        if (resultsStepPoints.max() > startModes_measure_large) and (optionalResultsStepPoints.max() > startModes_measure_large):
            fullMeasureInvestigation = True
        else:
            fullMeasureInvestigation = False
        # starting x axis value for both measures
        startingXvalue = self.getSmallestCommonValue(resultsStepPoints, optionalResultsStepPoints)
        resultsStepPointsStartingIndex = np.where(resultsStepPoints==startingXvalue)[0][0]
        optionalResultsStepPointsStartingIndex = np.where(optionalResultsStepPoints==startingXvalue)[0][0]
        
        # end x axis 
        resultsStepPointsEndIndex = np.where(resultsStepPoints==endModes_measure)[0][0]
        optionalResultsStepPointsEndIndex = np.where(optionalResultsStepPoints==endModes_measure)[0][0]
        
        # stepsize for the first k modes
        resultsStepSize = self.getStepsize(resultsStepPoints, resultsStepPointsStartingIndex, resultsStepPointsEndIndex)
        optionalResultsStepSize = self.getStepsize(optionalResultsStepPoints, optionalResultsStepPointsStartingIndex, optionalResultsStepPointsEndIndex)
        
        # Calculate the areas by using the composite trapezoidal rule
        measure = trapz(results[resultsStepPointsStartingIndex:resultsStepPointsEndIndex], dx=resultsStepSize) - trapz(optionalResults[optionalResultsStepPointsStartingIndex:optionalResultsStepPointsEndIndex], dx=optionalResultsStepSize)
        measure_large = 0.0
        maxModesUsed = endModes_measure
        
        if fullMeasureInvestigation:
        # bigger stepsize after k modes
            # starting x axis value for both measures
            resultsStepPointsStartingIndex_large = np.where(resultsStepPoints==startModes_measure_large)[0][0]
            optionalResultsStepPointsStartingIndex_large = np.where(optionalResultsStepPoints==startModes_measure_large)[0][0]
            
            # end x axis 
            endingXvalue = self.getLargestCommonValue(resultsStepPoints, optionalResultsStepPoints)
            resultsStepPointsEndIndex_large = np.where(resultsStepPoints==endingXvalue)[0][0]
            optionalResultsStepPointsEndIndex_large = np.where(optionalResultsStepPoints==endingXvalue)[0][0]
            
            resultsStepSize_large = self.getStepsize(resultsStepPoints, resultsStepPointsStartingIndex_large, resultsStepPointsEndIndex_large)
            optionalResultsStepSize2_large = self.getStepsize(optionalResultsStepPoints, optionalResultsStepPointsStartingIndex_large, optionalResultsStepPointsEndIndex_large)
        # Calculate the areas by using the composite trapezoidal rule
            measure_large = trapz(results[resultsStepPointsStartingIndex_large:resultsStepPointsEndIndex_large], dx=resultsStepSize_large) - trapz(optionalResults[optionalResultsStepPointsStartingIndex_large:optionalResultsStepPointsEndIndex_large], dx=optionalResultsStepSize2_large)
            maxModesUsed = endingXvalue
        return measure, measure_large, maxModesUsed

    def getAUC(self, arr, stepPoints, startK, endK, normalize=False):
        """ Get the area under the curve (AUC) from the y values of arr and x values of stepPoints. 
        startK and endK indicate the starting and ending values of stepPoints for the AUC. The AUC can be 
        returned as normalized, where it is assumed that the optimal AUC curve goes between 0..1 on the y axis
        
        Args:
            arr: list of y values of the results
            stepPoints: list of x values of the results
            startK: starting point (such as starting number of modes) for the AUC, needs to be member of stepPoints
            endK: end point (such as ending number of modes) for the AUC, needs to be a member of stepPoints
            normalize: If True, divide the AUC by the optimal AUC curve, which is assumed to go between 0..1 on the y axis
            
        Return: AUC for arr, stepPoints between startK and endK
        """
        assert startK in stepPoints
        assert endK in stepPoints
        assert startK < endK
        assert len(arr) > 1
        assert len(stepPoints) > 1
        
        auc = 0.0
        startK_index = np.where(stepPoints==startK)[0][0]
        endK_index = np.where(stepPoints==endK)[0][0]
        for i in range(startK_index, endK_index):
            auc += trapz(arr[i:i+2], dx=stepPoints[i+1]-stepPoints[i])
        
        # normalization assumes that the optimal curve is between 0 and 1 on the y axis
        if normalize:
            auc = auc / float(stepPoints[endK_index] - stepPoints[startK_index])
        
        return auc
    
    def getAUCList(self, arr, stepPoints, startK, endK, normalize=False):
        """ Get a list of the areas under the curve (AUC) from the y values of arr and x values of stepPoints, successively 
        increasing between between startK and the values in stepPoints, until finally reaching endK. 
        They indicate the overall starting and ending values of stepPoints for the AUC. The AUC can be 
        returned as normalized, where it is assumed that the optimal AUC curve goes between 0..1 on the y axis
        
        Args:
            arr: list of y values of the results
            stepPoints: list of x values of the results
            startK: starting point (such as starting number of modes) for the AUC, needs to be member of stepPoints
            endK: end point (such as ending number of modes) for the AUC, needs to be a member of stepPoints
            normalize: If True, divide the AUC by the optimal AUC curve, which is assumed to go between 0..1 on the y axis
            
        Return: list of AUCs for arr, successively increasing the distance between startK and the element in stepPoints until endK
        """
        assert startK in stepPoints
        assert endK in stepPoints
        assert startK < endK
        assert len(arr) > 1
        assert len(stepPoints) > 1
        
        aucList = []
        startK_index = np.where(stepPoints==startK)[0][0]
        endK_index = np.where(stepPoints==endK)[0][0]
        for i in range(startK_index, endK_index):
            auc = trapz(arr[i:i+2], dx=stepPoints[i+1]-stepPoints[i])       
            aucList.append(auc)
            
        aucListCumulsum = np.array(aucList).cumsum()
        
        # normalization assumes that the optimal curve is between 0 and 1 on the y axis
        if normalize:
#             aucList_normalized = []
#             for i in range(0, len(aucListCumulsum)):
#                 aucList_normalized.append(aucListCumulsum[i] / (float(stepPoints[i+1]) - float(stepPoints[startK_index])))
#             aucListCumulsum = np.array(aucList_normalized)
            aucListCumulsum = aucListCumulsum / (endK_index - startK_index)
            
        return aucListCumulsum
        
    def writeResultToFile(self, path, arr, fileName, fileFormat, fileHeader):
        """ Write the assessment result to a file. """
        try: 
            np.savetxt(path+fileName, arr, fmt=fileFormat, header=fileHeader)
        except AttributeError, err:
            print "Exception AttributeError occurred for " + str(fileName)+" : ", err
            print traceback.format_exc()
        except Exception, err:
            print "Exception occurred for " + str(fileName)+" : ", err
            print traceback.format_exc()
