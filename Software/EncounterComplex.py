'''
Created on Jan 30, 2014

@author: oliwa
'''
import sys
from prody.proteins.pdbfile import writePDB
import numpy as np

class EncounterComplex(object):
    '''
    classdocs
    '''


    def __init__(self, receptorStructure, ligandStructure, complexType, utils):
        '''
        Constructor
        The complex types are:
        boundComplex:           bound complex from parsed PDB filed
        unboundComplex          (not implemented yet, for 2b)
        boundComplexAligned     (not implemented yet, for 2b)
        unboundComplexAlinged:  unbound complex from parsed PDBs and aligned to bound counterparts
        '''
        assert complexType == "boundComplex" or complexType == "unboundComplexAligned"
        self.utils = utils
        self.complexType = complexType       
        
        # complex of unmatched pdbs
        self.receptorStructure = self.setupSegments(receptorStructure, "R")
        self.ligandStructure = self.setupSegments(ligandStructure, "L")
        self.complex = self.setupComplex(self.receptorStructure, self.ligandStructure, utils)
        
    def setupSegments(self, protein, typeOfBindingProtein):
        """ For a protein, Write the segment of each atom as receptor or ligand followed by 
        its chainID in the way of: <R|L chainID> (no space between them) 
        Args:
            protein: the protein where the segments are to be written to
            typeOfBindingProtein: "R" for receptor or "L" for ligand
            
        Returns: the protein with their segments written
        """
        assert typeOfBindingProtein == "R" or typeOfBindingProtein == "L"      
        segmentNames = []
        for i in range(0, protein.numAtoms()):
            segmentNames.append(str(typeOfBindingProtein+protein[i].getChid()+"                                               "))
        try:
            protein.setSegnames(segmentNames)
        except Exception():
            raise Exception()
        return protein
    
    def setupComplex(self, receptorStructure, ligandStructure, utils):
        """ Setup the proteinComplex from receptor and ligand, then copy their coordinates, since 
        creating the proteinComplex involves writing/reading a PDB structure, which has only 
        three after comma digits per xyz value. """
        proteinComplex = utils.getUnifiedPDBfromPDBs(receptorStructure, ligandStructure)
        proteinComplex = self.copyCoords(receptorStructure, ligandStructure, proteinComplex)
        return proteinComplex
    
    def copyCoords(self, receptorStructure, ligandStructure, proteinComplex, checkIfStructuresAreClose=True):
        """ Copy the coordinates of protein to the same atoms of protein complex and return the protein complex. 
        
        Args:
            receptorStructure: the receptor structure, top part of the pdb file
            ligandStructure: the ligand structure, bottom part of the pdb file
            proteinComplex: the object representing parsed receptor and ligand
            checkIfStructuresAreClose: If true, assert that the coordinates of the same atoms in receptor/ligand and complex
                                        are the same after rounding. This is a useful check if the coordinates have to be equalized when they 
                                        are off by some floating point imprecision, and the check asserts that indeed for all atoms the right 
                                        coordinates are selected.
        Returns:
            proteinComplex with coordinates from receptorStructure and ligandStructure
        """
        assert receptorStructure.numAtoms() == proteinComplex.select('segment "R."').numAtoms()
        assert ligandStructure.numAtoms() == proteinComplex.select('segment "L."').numAtoms()
        newCoords = np.vstack((receptorStructure.getCoords(), ligandStructure.getCoords()))
        assert len(newCoords) == len(proteinComplex.getCoords())

        if checkIfStructuresAreClose:
            for i in range(0, receptorStructure.numAtoms()):
                item1roundedChoords = [round(x, 3) for x in receptorStructure[i].getCoords().tolist()] 
                item2roundedChoords = [round(x, 3) for x in proteinComplex[i].getCoords().tolist()] 
                assert np.alltrue(item1roundedChoords == item2roundedChoords)
            for j in range(0, ligandStructure.numAtoms()):
                item1roundedChoords = [round(x, 3) for x in ligandStructure[j].getCoords().tolist()] 
                item2roundedChoords = [round(x, 3) for x in proteinComplex[j+receptorStructure.numAtoms()].getCoords().tolist()]
                assert np.alltrue(item1roundedChoords == item2roundedChoords)            
        proteinComplex.setCoords(newCoords)
        return proteinComplex
    
    def copyCoordsFromComplex(self, proteinComplex, segment, individualProtein):
        """ Copy the relevant coordinates of the proteinComplex to the individualProtein. The individualProtein 
        is assumed to be a receptor or a ligand of this complex, as specified by the segment, with the same atoms, 
        this is asserted by this method. """
        assert self.utils.returnEqualityOfProteinsWithoutCoordinates(individualProtein, proteinComplex.select('segment \"'+segment+'.\"'))
        newCoords = proteinComplex.select('segment \"'+segment+'.\"').getCoords()
        individualProtein.setCoords(newCoords)
        assert np.alltrue(individualProtein.getCoords() == proteinComplex.select('segment \"'+segment+'.\"').getCoords())
