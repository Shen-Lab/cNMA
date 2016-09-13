'''
Created on Jan 17, 2014

@author: oliwa
'''

import traceback
from ANMs import ANMs
from EncounterComplex import EncounterComplex
from ResultsPrinter import ResultsPrinter
from prody.proteins.pdbfile import parsePDB, writePDB
import sys as sys
import numpy as np
from prody.measure.measure import calcDistance

class Encounter(object):
	"""
	This class holds the protein members of one encounter experiment and their 
	getter/setter with protection against getting unset members.
	"""

	def __init__(self, pdbQueueItem):
		"""
		Constructor, initialize the encounter with the PDB name and its path in the benchmark,
		specified in the pdbQueueItem tuple
		"""
		self._pdbQueueItem = pdbQueueItem
	
	def setReference(self, reference):
		""" Set the reference structure (unbound protein conformation)
		and create a backup member. """
		self._reference = reference
		self._referenceBackup = reference.copy()
		
	def updateReference(self, reference, title):
		self._reference = reference.copy()
		self._reference.setTitle(title)
		
	def setMobile(self, mobile):
		""" Set the mobile structure (bound protein conformation)
		and create a backup member. """
		self._mobile = mobile
		self._mobileBackup = mobile.copy()
		
	def updateMobile(self, mobile, title):
		self._mobile = mobile.copy()
		self._mobile.setTitle(title)
		
	def loadUnfilteredReceptor(self, utils):
		proteinName = self._pdbQueueItem[0]
		proteinName = list(proteinName)
		proteinName[5] = "r"
		proteinName = "".join(proteinName)
		path = utils.config.pathTo2cFiles+proteinName+".pdb.ms"
		self._unfilteredReferenceAtomCount = utils.file_len(path)
		self._unfilteredReference = parsePDB(path)
		
	def setRefChain(self, ref_chain):
		""" Set the matched chain of reference structure. """
		self._ref_chain = ref_chain
		self._ref_chainBackup = ref_chain.copy()
		
	def setMobChain(self, mob_chain):
		""" Set the matched chain of mobile structure. """
		self._mob_chain = mob_chain
		self._mob_chainBackup = mob_chain.copy()
		
	def setUnboundCounterpartChain(self, unboundCounterpartChain):
		""" Set the matched chain of the unbound counterpart. """
		self._unboundCounterpartChain = unboundCounterpartChain
		self._unboundCounterpartChainBackup = unboundCounterpartChain.copy()
		
	def setBoundCounterpartChain(self, boundCounterpartChain):
		""" Set the matched chain of the bound counterpart. """
		self._boundCounterpartChain = boundCounterpartChain
		self._boundCounterpartChainBackup = boundCounterpartChain.copy()
		
	def setUnboundComplexAlignedChain(self, unboundComplexAlignedChain):
		""" Set the matched chain of the unbound complex (which had its receptor and ligand transformed on their bound parts). """
		self._unboundComplexAlignedChain = unboundComplexAlignedChain
		self._unboundComplexAlignedChainBackup = unboundComplexAlignedChain.copy()
		
	def setBoundComplexChain(self, boundComplexChain):
		""" Set the matched chain of the bound complex. """
		self._boundComplexChain = boundComplexChain
		self._boundComplexChainBackup = boundComplexChain.copy()
			
	def setupBoundComplex(self, utils):
		if utils.isReceptor(self.getMobile().getTitle()):
			self.boundComplex = EncounterComplex(self._mobileBackup.copy(), self._boundCounterpartBackup.copy(), "boundComplex", utils)
		else:
			self.boundComplex = EncounterComplex(self._boundCounterpartBackup.copy(), self._mobileBackup.copy(), "boundComplex", utils)
			
	def setupUnboundComplexAligned(self, utils):
		if utils.isReceptor(self.getReference().getTitle()):
			self.unboundComplexAligned = EncounterComplex(self.getReference(), self.getUnboundCounterpart(), "unboundComplexAligned", utils)
		else:
			self.unboundComplexAligned = EncounterComplex(self.getUnboundCounterpart(), self.getReference(), "unboundComplexAligned", utils)
	
	def calcRefChainInterface(self, utils):
		""" Calculate the matched interface atoms (subset of the ref_chain) of the reference structure. """
		ref_chain_interface_indices = utils.getMatchingStructure(self.getMobChain(), self.getMobChainInterface(), self.getRefChain())
		self._ref_chain_interface = utils.getSubsetOfSelection(self.getRefChain(), ref_chain_interface_indices)
		
	def calcUnboundInterface(self, boundChain, boundChainInterface, unboundChain, utils):
		""" Get the interface on the unbound protein chain when given the bound protein chain and the known interface there. 
		Generalization of calcRefChainInterface """
		unbound_interface_indices = utils.getMatchingStructure(boundChain, boundChainInterface, unboundChain)
		return utils.getSubsetOfSelection(unboundChain, unbound_interface_indices) 
	
	def setUnboundCounterpartInterface(self, unboundCounterpart_interface):
		self._unboundCounterpartChain_interface = unboundCounterpart_interface
				 
	def calcMobChainInterface(self):
		""" Calculate the matched interface atoms (subset of the mob_chain) of the mobile structure. """
		# old way for interface
		#self._mob_chain_interface = self.getMobChain().select('same residue as within 6 of inhibitor', inhibitor=self.getBoundCounterpart())
		interfaceMobile = self.getMobile().select('same residue as within 6 of inhibitor', inhibitor=self.getBoundCounterpart())
		self._mob_chain_interface = self.getMobChain().select(interfaceMobile.getSelstr())
		
	def calcBoundCounterpartChainInterface(self):
		""" Calculate the matched interface atoms (subset of boundCounterpartChain) of the the bound counterpart"""
		# old way for interface
		#self._boundCounterpartChain_interface = self.getBoundCounterpartChain().select('same residue as within 6 of inhibitor', inhibitor=self.getMobile())
		interfaceCounterPart = self.getBoundCounterpart().select('same residue as within 6 of inhibitor', inhibitor=self.getMobile())
		self._boundCounterpartChain_interface = self.getBoundCounterpartChain().select(interfaceCounterPart.getSelstr())

	def calcBoundCounterpart(self, utils):
		""" If reference is XYZA_l_u, the boundCounterPart is XYZA_r_b. 
			If reference is XYZA_r_u, the boundCounterPart is XYZA_l_b
		"""
		self._boundCounterpart = utils.parseCorrenspondingCounterpart(self.getReference(), "bound")
		self._boundCounterpartBackup = self._boundCounterpart.copy()
		
	def updateBoundCounterpart(self, boundCounterpart, title):
#         print "before updating boundCounterpart: "
#         print "old boundCounterpart, new boundCounterpart"
#         for item1, item2 in zip(self._boundCounterpart, boundCounterpart):
#             print item1, item2
		self._boundCounterpart = boundCounterpart.copy()
		self._boundCounterpart.setTitle(title)
		
	def calcUnboundCounterpart(self, utils):
		""" If reference is XYZA_l_u, the boundCounterPart is XYZA_r_u. 
			If reference is XYZA_r_u, the boundCounterPart is XYZA_l_u
		"""
		self._unboundCounterpart = utils.parseCorrenspondingCounterpart(self.getReference(), "unbound")
		self._unboundCounterpartBackup = self._unboundCounterpart.copy()
		
	def updateUnboundCounterpart(self, unboundCounterpart, title):
		self._unboundCounterpart = unboundCounterpart.copy()
		self._unboundCounterpart.setTitle(title)
		
	def calcUnboundCounterpart2c(self, utils, cCase):
		""" Extract the unbound counterpart from the combined pdb file of 2c based on the cCase number
		"""
		self._unboundCounterpart = utils.parseCorrenspondingCounterpart2c(self.getReference(), "unbound", cCase, self._unfilteredReferenceAtomCount)
		self._unboundCounterpartBackup = self._unboundCounterpart.copy()
		
	def setBoundCounterpart(self, boundCounterpart):
		""" Set the boundCounterpart directly . """
		self._boundCounterpart = boundCounterpart
		self._boundCounterpartBackup = boundCounterpart.copy()
		
	def setUnboundCounterpart(self, unboundCounterpart):
		""" Set the unbound counterpart directly. """
		self._unboundCounterpart = unboundCounterpart
		self._unboundCounterpartBackup = unboundCounterpart.copy()
		
	def calcBoundComplexChainInterface(self):
		interfaceRB = self.boundComplex.complex.select('same residue as exwithin 6 of segment "L." ')
		interfaceLB = self.boundComplex.complex.select('same residue as exwithin 6 of segment "R." ')
		interfaceRBLB = interfaceRB + interfaceLB
		self._boundComplexChainInterface = self.getBoundComplexChain().select(interfaceRBLB.getSelstr())
		
	def setUnboundComplexChainInterface(self, unboundComplexChainInterface):
		self._unboundComplexChainInterface = unboundComplexChainInterface
		
	def setMobilePath(self, mobilePath):
		""" Set the path to the mobile structure. """
		self._mobilePath = mobilePath
		
	def initANMs(self, utils):
		""" Initialize AMNs, where all the ANMs are stored. """
		self._ANMs = ANMs(utils)
		
	def initResultsPrinter(self, pdbName=None):
		if pdbName==None:
			self.resultsPrinter = ResultsPrinter(self.getPdbQueueItem()[0])
		else: 
			self.resultsPrinter = ResultsPrinter(pdbName)
		
	def getReference(self):
		if self._reference == None:
			raise Exception('self._reference == None')
		return self._reference
	
	def getMobile(self):
		if self._mobile == None:
			raise Exception('self._mobile == None')
		return self._mobile
	
	def getPdbQueueItem(self):
		if self._pdbQueueItem == None:
			raise Exception('self._pdbQueueItem == None')
		return self._pdbQueueItem
	
	def getRefChain(self):
		if self._ref_chain == None:
			raise Exception('self._ref_chain == None')
		return self._ref_chain
	
	def getMobChain(self):
		if self._mob_chain == None:
			raise Exception('self._mob_chain == None')
		return self._mob_chain
	
	def getRefChainInterface(self):
		if self._ref_chain_interface == None:
			raise Exception('self._ref_chain_interface == None')
		return self._ref_chain_interface
	
	def getMobChainInterface(self):
		if self._mob_chain_interface == None:
			raise Exception('self._mob_chain_interface')
		return self._mob_chain_interface
	
	def getBoundCounterpart(self):
		if self._boundCounterpart == None:
			raise Exception('self._boundCounterpart == None')
		return self._boundCounterpart
	
	def getUnboundCounterpart(self):
		if self._unboundCounterpart == None:
			raise Exception('self._unboundCounterpart == None')
		return self._unboundCounterpart
	
	def getBoundCounterpartChainInterface(self):
		if self._boundCounterpartChain_interface == None:
			raise Exception('self._boundCounterpartChain_interface == None')
		return self._boundCounterpartChain_interface
	
	def getUnboundCounterpartChainInterface(self):
		if self._unboundCounterpartChain_interface == None:
			raise Exception('self._unboundCounterpartChain_interface == None')
		return self._unboundCounterpartChain_interface
	
	def getUnboundCounterpartChain(self):
		if self._unboundCounterpartChain == None:
			raise Exception('self._unboundCounterpartChain == None')
		return self._unboundCounterpartChain
	
	def getBoundCounterpartChain(self):
		if self._boundCounterpartChain == None:
			raise Exception('self._boundCounterpartChain == None')
		return self._boundCounterpartChain
	
	def getUnboundComplexAlignedChain(self):
		if self._unboundComplexAlignedChain == None:
			raise Exception('UnboundComplexAlignedChain == None')
		return self._unboundComplexAlignedChain
	
	def getBoundComplexChain(self):
		if self._boundComplexChain == None:
			raise Exception('BoundComplexChain == None')
		return self._boundComplexChain
	
	def getBoundComplexChainInterface(self):
		if self._boundComplexChainInterface == None:
			raise Exception('boundComplexChainInterface == None')
		return self._boundComplexChainInterface
	
	def getUnboundComplexChainInterface(self):
		if self._unboundComplexChainInterface == None:
			raise Exception('unboundComplexChainInterface == None')
		return self._unboundComplexChainInterface
	
	def getUnfilteredReferenceAtomCount(self):
		if self._unfilteredReferenceAtomCount == None:
			raise Exception('self._unfilteredReferenceAtomCount == None')
		return self._unfilteredReferenceAtomCount 
	
	def getMobilePath(self):
		if self._mobilePath == None:
			raise Exception('self._mobilePath == None')
		return self._mobilePath
	
	def accessANMs(self):
		if self._ANMs == None:
			raise Exception('self._ANMs == None')
		return self._ANMs
	
	def setReferenceSegment(self, referenceSegment):
		self._referenceSegment = referenceSegment
		
	def setCounterpartSegment(self, counterpartSegment):
		self._counterpartSegment = counterpartSegment
		
	def getReferenceSegment(self):
		if self._referenceSegment == None:
			raise Exception('self._referenceSegment == None')
		return self._referenceSegment
	
	def getCounterpartSegment(self):
		if self._counterpartSegment == None:
			raise Exception('self._counterpartSegment == None')
		return self._counterpartSegment
	
	def fixChids(self, utils, bound_provided):
		self.setReference(utils.getFilledChidProtein(self.getReference(), self.getUnboundCounterpart()))
		self.setUnboundCounterpart(utils.getFilledChidProtein(self.getUnboundCounterpart(), self.getReference()))
		if bound_provided == True:
			self.setMobile(utils.getFilledChidProtein(self.getMobile(), self.getBoundCounterpart()))
			self.setBoundCounterpart(utils.getFilledChidProtein(self.getBoundCounterpart(), self.getMobile()))
		
	def assert2cProteinsAreEqual(self, utils):
		proteinName = self._pdbQueueItem[0]
		if utils.isReceptor(self.getReference().getTitle()):
			proteinName = list(proteinName)
			proteinName[5] = "l"
			proteinName = "".join(proteinName)
		else:
			proteinName = list(proteinName)
			proteinName[5] = "r"
			proteinName = "".join(proteinName)
		path = utils.config.pathTo2cFiles+proteinName+".pdb.ms"
		counterpart = utils.parseFilteredPDB(path, utils.config.filterPDB)
		#print counterpart
		#print "self.getUnboundCounterpart().getTitle(): ", self.getUnboundCounterpart().getTitle()
		assert utils.returnEqualityOfProteinsWithoutCoordinates(self.getUnboundCounterpart(), counterpart)
		
	def printIntermolecularNeighbors(self, structure, counterpart, atomType, distance):
		""" Print the neighbors in counterpart, having a certain atom type and (based on distance in Angstrom) of atoms of the same type of
		atoms in structure."""
		structureElements = structure.select(atomType)
		counterpartElement = counterpart.select(atomType)
		elementCounter = 0
		for element in structureElements:
			elementCoords = element.getCoords()
			contacts = counterpartElement.select('within '+distance+' of somepoint', somepoint=elementCoords)
			if contacts:
				print str(elementCounter)+"'th calpha, hessian line number", str((elementCounter*3)+1)," contact calphas within", distance, "angstrom of", element, ": ", contacts 
				for singleContact in contacts:
					print "distance to "+str(singleContact)+": ", calcDistance(elementCoords, singleContact.getCoords())
#                 print "reference element coords: ", element.getCoords()
#                 print "coords of reference contact:", contacts.getCoords()
#                 sys.exit()
				#for neighbor in contacts:
				#    print neighbor
			elementCounter += 1
					
	def getIntermolecularNeighbors(self, structure, counterpart, structureAtomIndex, atomType, distance):
		""" Return the intermolecular contacts, situated in counterpart, of an atom of atomType at index structureAtomIndex from structure. """
		structureElements = structure.select(atomType)
		counterpartElement = counterpart.select(atomType)
		element = structureElements[structureAtomIndex]
		elementCoords = element.getCoords()
		contacts = counterpartElement.select('within '+distance+' of somepoint', somepoint=elementCoords)
		return contacts
	
	def getIntermolecularNeighborsOfAtom(self, atom1, counterpart, atomType, distance):
		""" Return the intermolecular contacts, situated in counterpart, of atom atom1 with atomType. """
		if atomType == 'calpha':
			assert atom1.getName() == 'CA'
		counterpartElements = counterpart.select(atomType)
		elementCoords = atom1.getCoords()
		contacts = counterpartElements.select('within '+distance+' of somepoint', somepoint=elementCoords)
		return contacts
	
	def storeOverlapTable(self, overlapTable):
		self.overlapTable = overlapTable
		
	def storeSubSpaceOverlaps(self, subspaceOverlaps, subspaceOverlapsRanges):
		self.subspaceOverlaps = subspaceOverlaps
		self.subspaceOverlapsRanges = subspaceOverlapsRanges
		
	def storeCovarianceOverlap(self, covarianceOverlaps, covarianceOverlapsRanges):
		self.covarianceOverlaps = covarianceOverlaps
		self.covarianceOverlapsRanges = covarianceOverlapsRanges
		
	def affirmProteinNames(self, protein1Type):
		""" Adhere Benchmark 4.0 naming standard for the parsed proteins """
		if protein1Type == "receptor":
			pro1 = "r"
			pro2 = "l"
		elif protein1Type == "ligand":
			pro1 = "l"
			pro2 = "r"
		self._pdbQueueItem = list(self._pdbQueueItem)              
		self._pdbQueueItem[0] = "pro1_"+pro1+"_u_"+self._pdbQueueItem[0]
		self._pdbQueueItem = tuple(self._pdbQueueItem)
		
		self._reference.setTitle("pro1_"+pro1+"_u_"+self._reference.getTitle())
		self._referenceBackup.setTitle("pro1_"+pro1+"_u_"+self._referenceBackup.getTitle())
		
		self._mobile.setTitle("pro1_"+pro1+"_b_"+self._mobile.getTitle())
		self._mobileBackup.setTitle("pro1_"+pro1+"_b_"+self._mobileBackup.getTitle())
		
		self._unboundCounterpart.setTitle("pro2_"+pro2+"_u_"+self._unboundCounterpart.getTitle())
		self._unboundCounterpartBackup.setTitle("pro2_"+pro2+"_u_"+self._unboundCounterpartBackup.getTitle())
		
		self._boundCounterpart.setTitle("pro2_"+pro2+"_b_"+self._boundCounterpart.getTitle())
		self._boundCounterpartBackup.setTitle("pro2_"+pro2+"_b_"+self._boundCounterpartBackup.getTitle())   
		
		