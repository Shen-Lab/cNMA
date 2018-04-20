'''
Created on Apr 27, 2014

@author: oliwa
'''

from TNMABase import TNMABase
from Assessment import Assessment
from Encounter import Encounter
from ModeClassifier import ModeClassifier
from RMSDReducer import RMSDReducer
from RMSDPrediction import RMSDPrediction
from Utils import Utils
from DataHolderForDirectOutput import DataHolderForDirectOutput
import argparse
import os
import sys
import numpy as np
from prody.measure.measure import calcDeformVector
from prody.measure.transform import calcRMSD, calcTransformation
from prody.dynamics.compare import calcCumulOverlap
# for profiling
import pstats
import cProfile
from prody.dynamics.sampling import sampleModes, traverseMode
from prody.dynamics.editing import extendModel, sliceModel
from prody.dynamics.functions import writeArray
from helperScripts.scriptutils import makeStringEndWith
from prody.proteins.pdbfile import writePDB
from prody.dynamics.anm import ANM
import traceback

class NMAUnified(TNMABase):
	'''
	NMAUnified inherits from TNMABase and performs a NMA for individual and complex proteins.
	May 6th: Complex Investigations implemented
	
	todo: Individual investigations, mode measures, later investigate complex <-> individual protein representation 
	'''

	def __init__(self, args):
		'''
		Constructor of NMAUnified. Dynamically loads the configuration file and creates Utils and Assessment helper objects.
		
		Args:
			args: argparse arguments
		'''       
		self.args = args
		if args.affirmProteinNames:
			self.affirmProteinNames = args.affirmProteinNames
			if (args.affirmProteinNames != "receptor") and (args.affirmProteinNames != "ligand"):
				raise Exception("Exception, with the affirmProteinNames parameter, a setting other than \"receptor\" or \"ligand\" was provided.")
		else:
			self.affirmProteinNames = None
		self.config = self.instantiateClassFromModule(self.importFromURI(args.config), 'Configurations')
		self.pathOfConfigFile = self.getFullPathOfURI(args.config)
		self.utils = Utils(self.config)
		self.assessment = Assessment(self.utils, 1)
		self.assessment.setExperimentName(self.config.whatAtomsToMatch, self.config.experimentNamePrefix)
	
	def setupEncounter(self, protein1_A, protein2_A=None, protein1_B=None, protein2_B=None, affirmProteinNames=None):
		""" Setup the encounter object which holds the proteins. 
		
		Args:
			protein1_A: Protein 1 in first conformational state
			protein2_A: Protein 2 in first conformational state
			protein1_B: Protein 1 in second conformational state
			protein2_B: Protein 2 in second conformational state
			
		Returns: 
			Encounter object    
		"""

		encounter = Encounter((os.path.splitext(os.path.basename(protein1_A))[0], protein1_A))
		encounter.setReference(self.utils.parseFilteredPDB(protein1_A, self.config.filterPDB))
	
		encounter.setUnboundCounterpart(self.utils.parseFilteredPDB(protein2_A, self.config.filterPDB))

		# use bound protein
		if self.bound_provided == True:
			encounter.setMobile(self.utils.parseFilteredPDB(protein1_B, self.config.filterPDB))
			encounter.setBoundCounterpart(self.utils.parseFilteredPDB(protein2_B, self.config.filterPDB))

		if affirmProteinNames:
			# adjust the titles of the proteins to match the Benchmark 4.0 standard
			encounter.affirmProteinNames(self.affirmProteinNames)
		encounter.fixChids(self.utils, self.bound_provided)
		if self.bound_provided == True:
			encounter.setupBoundComplex(self.utils) # still relies on naming, TODO
		encounter.setupUnboundComplexAligned(self.utils) # still relies on naming, TODO
		# set segment information, still relies on naming, TODO
		if self.utils.isReceptor(encounter.getReference().getTitle()):
			encounter.setReferenceSegment("R")
			encounter.setCounterpartSegment("L")
		else:
			encounter.setReferenceSegment("L")
			encounter.setCounterpartSegment("R")        
		return encounter
	
	def chainMatching(self, encounter):
		""" Match the chains of the proteins in the encounter object. 
		
		Args:
			encounter: the encounter object
		"""
		print "matching reference"
		overallMatch = self.getOverallMatch(encounter.getReference(), encounter.getMobile(), self.config.whatAtomsToMatch)

		encounter.setRefChain(overallMatch[0])
		encounter.setMobChain(overallMatch[1])
		print "overall match: \n" + str(encounter.getRefChain()) + "\n" + str(encounter.getMobChain())
		print "overallMatch: ", overallMatch        
		
		print "matching counterpart"
		overallMatchCounterpart = self.getOverallMatch(encounter.getUnboundCounterpart(), encounter.getBoundCounterpart(), self.config.whatAtomsToMatch)
		encounter.setUnboundCounterpartChain(overallMatchCounterpart[0])
		encounter.setBoundCounterpartChain(overallMatchCounterpart[1])
		print "overall match counterpart: \n" + str(encounter.getUnboundCounterpartChain()) + "\n" + str(encounter.getBoundCounterpartChain())
		print "overallMatchCounterpart: ", overallMatchCounterpart

		print "matching complex"
		overallMatchComplex = self.getOverallMatch(encounter.unboundComplexAligned.complex, encounter.boundComplex.complex, self.config.whatAtomsToMatch)
		encounter.setUnboundComplexAlignedChain(overallMatchComplex[0])
		encounter.setBoundComplexChain(overallMatchComplex[1])
		print "overall match complex : \n" + str(encounter.getUnboundComplexAlignedChain()) + "\n" + str(encounter.getBoundComplexChain())
		print "overallMatchComplex: ", overallMatchComplex              
	
	def alignProteins(self, encounter, alignStyle, alignWithWhole, dataHolder=None, storeUnboundComplex=False, storePDBs=False):
		""" Align (rigid-bodily transpose) proteins onto each other: 
		alpha - protein1 and protein2 independently superposed
		beta - protein1 superposed and protein2 "dragged along"
		complexOnComplex - unbound complex superposed on bound complex
		2bcomplex - first align alpha whole to create the complex, then align the complex again based on complexOnComplex for whole or interface
		2bindividual - first align alpha whole to create the complex, then align the complex again based on beta for whole or interface
		L-RMS - as beta, but always use the receptor
		
		Args:
			encounter: the encounter object
			alignStyle: alpha, beta or complexOnComplex 
			alignWithWhole: If true, use whole proteins as bases, if false use interfaces as bases of alignments
			dataHolder: holds the direct results
			storeUnboundComplex: store a copy of the moved, aligned unbound complex
			storePDBs: store pdbs from encounter into the dataholder
		"""
		assert alignStyle == "alpha" or alignStyle == "beta" or alignStyle == "complexOnComplex" or alignStyle == "2bcomplex" or alignStyle == "2bindividual" or alignStyle == "L-RMS"
		# pointer variable to receptor and ligands
		if self.utils.isReceptor(encounter.getReference().getTitle()):
			receptor = encounter.getReference()
			ligand = encounter.getUnboundCounterpart()
		else:
			receptor = encounter.getUnboundCounterpart()
			ligand = encounter.getReference()
		
		if alignWithWhole:
			# align protein1 and protein2 independently
			if alignStyle == "alpha":
				# align protein1
				t = calcTransformation(encounter.getRefChain(), encounter.getMobChain())
				t.apply(encounter.getReference())
				# align protein2
				t_counterpart = calcTransformation(encounter.getUnboundCounterpartChain(), encounter.getBoundCounterpartChain())
				t_counterpart.apply(encounter.getUnboundCounterpart())       
				# copy coordinates from the moved proteins 1 and 2 on the complex representation
				encounter.unboundComplexAligned.complex = encounter.unboundComplexAligned.copyCoords(receptor, ligand, encounter.unboundComplexAligned.complex, checkIfStructuresAreClose=False)
	
			# align protein 1 and rigid bodily "drag" protein 2 along
			elif alignStyle == "beta":
				t_newcomplex = calcTransformation(encounter.getUnboundComplexAlignedChain().select('segment \"'+encounter.getReferenceSegment()+'.\"'), 
												  encounter.getBoundComplexChain().select('segment \"'+encounter.getReferenceSegment()+'.\"'))
				t_newcomplex.apply(encounter.unboundComplexAligned.complex)
				# copy coordinates from the moved complex onto proteins 1 and 2 in their individual representation
				encounter.unboundComplexAligned.copyCoordsFromComplex(encounter.unboundComplexAligned.complex, encounter.getReferenceSegment(), encounter.getReference())
				encounter.unboundComplexAligned.copyCoordsFromComplex(encounter.unboundComplexAligned.complex, encounter.getCounterpartSegment(), encounter.getUnboundCounterpart())
			  
			# align protein 1 and rigid bodily "drag" protein 2 along, always aligning on receptor
			elif alignStyle == "L-RMS":
				t_newcomplex = calcTransformation(encounter.getUnboundComplexAlignedChain().select('segment \"R.\"'), 
												  encounter.getBoundComplexChain().select('segment \"R.\"'))
				t_newcomplex.apply(encounter.unboundComplexAligned.complex)
				# copy coordinates from the moved complex onto proteins 1 and 2 in their individual representation
				encounter.unboundComplexAligned.copyCoordsFromComplex(encounter.unboundComplexAligned.complex, encounter.getReferenceSegment(), encounter.getReference())
				encounter.unboundComplexAligned.copyCoordsFromComplex(encounter.unboundComplexAligned.complex, encounter.getCounterpartSegment(), encounter.getUnboundCounterpart())              
			  
			# align complex on complex
			elif alignStyle == "complexOnComplex":
				t_newcomplex = calcTransformation(encounter.getUnboundComplexAlignedChain(), encounter.getBoundComplexChain())
				t_newcomplex.apply(encounter.unboundComplexAligned.complex)
				# copy coordinates from the moved complex onto proteins 1 and 2 in their individual representation
				encounter.unboundComplexAligned.copyCoordsFromComplex(encounter.unboundComplexAligned.complex, encounter.getReferenceSegment(), encounter.getReference())
				encounter.unboundComplexAligned.copyCoordsFromComplex(encounter.unboundComplexAligned.complex, encounter.getCounterpartSegment(), encounter.getUnboundCounterpart())
					  
			# first align alpha whole to create the complex, then align the complex again based on complexOnComplex for whole or interface
			elif alignStyle == "2bcomplex":
				# align protein1
				t = calcTransformation(encounter.getRefChain(), encounter.getMobChain())
				t.apply(encounter.getReference())
				# align protein2
				t_counterpart = calcTransformation(encounter.getUnboundCounterpartChain(), encounter.getBoundCounterpartChain())
				t_counterpart.apply(encounter.getUnboundCounterpart())   
				# copy coordinates from the moved proteins 1 and 2 on the complex representation
				encounter.unboundComplexAligned.complex = encounter.unboundComplexAligned.copyCoords(receptor, ligand, encounter.unboundComplexAligned.complex, checkIfStructuresAreClose=False)
				# now align based on complexOnComplex
				t_newcomplex = calcTransformation(encounter.getUnboundComplexAlignedChain(), encounter.getBoundComplexChain())
				t_newcomplex.apply(encounter.unboundComplexAligned.complex)
				# copy coordinates from the moved complex onto proteins 1 and 2 in their individual representation
				encounter.unboundComplexAligned.copyCoordsFromComplex(encounter.unboundComplexAligned.complex, encounter.getReferenceSegment(), encounter.getReference())
				encounter.unboundComplexAligned.copyCoordsFromComplex(encounter.unboundComplexAligned.complex, encounter.getCounterpartSegment(), encounter.getUnboundCounterpart())                    
			
			# first align alpha whole to create the complex, then align the complex again based on beta of the reference segment for whole or interface
			elif alignStyle == "2bindividual":
				# align protein1
				t = calcTransformation(encounter.getRefChain(), encounter.getMobChain())
				t.apply(encounter.getReference())
				# align protein2
				t_counterpart = calcTransformation(encounter.getUnboundCounterpartChain(), encounter.getBoundCounterpartChain())
				t_counterpart.apply(encounter.getUnboundCounterpart())       
				# copy coordinates from the moved proteins 1 and 2 on the complex representation
				encounter.unboundComplexAligned.complex = encounter.unboundComplexAligned.copyCoords(receptor, ligand, encounter.unboundComplexAligned.complex, checkIfStructuresAreClose=False)                
				# now align beta 
				t_newcomplex = calcTransformation(encounter.getUnboundComplexAlignedChain().select('segment \"'+encounter.getReferenceSegment()+'.\"'), 
												  encounter.getBoundComplexChain().select('segment \"'+encounter.getReferenceSegment()+'.\"'))
				t_newcomplex.apply(encounter.unboundComplexAligned.complex)
				# copy coordinates from the moved complex onto proteins 1 and 2 in their individual representation
				encounter.unboundComplexAligned.copyCoordsFromComplex(encounter.unboundComplexAligned.complex, encounter.getReferenceSegment(), encounter.getReference())
				encounter.unboundComplexAligned.copyCoordsFromComplex(encounter.unboundComplexAligned.complex, encounter.getCounterpartSegment(), encounter.getUnboundCounterpart())
				
			else:
				raise Exception("No such alignment know.")
			
			if storeUnboundComplex:
				self.unboundComplexCopy = encounter.unboundComplexAligned.complex.copy()  
			
			# assert equality of individual proteins and the complexes they form
			self.equalityAssertionsOfComplexAndItsProteins(encounter, checkInterface=False)
		else:
			# align protein1 and protein2 independently
			if alignStyle == "alpha":
				# align protein1 interface
				t_interface = calcTransformation(encounter.getRefChainInterface(), encounter.getMobChainInterface())
				t_interface.apply(encounter.getReference())                
				# align protein2 interface
				t_counterpart = calcTransformation(encounter.getUnboundCounterpartChainInterface(), encounter.getBoundCounterpartChainInterface())
				t_counterpart.apply(encounter.getUnboundCounterpart())
				# copy coordinates from the moved proteins 1 and 2 on the complex representation
				encounter.unboundComplexAligned.complex = encounter.unboundComplexAligned.copyCoords(receptor, ligand, encounter.unboundComplexAligned.complex, checkIfStructuresAreClose=False)
			
			# align protein 1 interface part and rigid bodily "drag" protein 2 along
			elif alignStyle == "beta":
				t_newcomplex = calcTransformation(encounter.getUnboundComplexChainInterface().select('segment \"'+encounter.getReferenceSegment()+'.\"'), 
												  encounter.getBoundComplexChainInterface().select('segment \"'+encounter.getReferenceSegment()+'.\"'))
				t_newcomplex.apply(encounter.unboundComplexAligned.complex)
				# copy coordinates from the moved complex onto proteins 1 and 2 in their individual representation
				encounter.unboundComplexAligned.copyCoordsFromComplex(encounter.unboundComplexAligned.complex, encounter.getReferenceSegment(), encounter.getReference())
				encounter.unboundComplexAligned.copyCoordsFromComplex(encounter.unboundComplexAligned.complex, encounter.getCounterpartSegment(), encounter.getUnboundCounterpart())              
			  
			# align complex on complex
			elif alignStyle == "complexOnComplex":
				t_newcomplex = calcTransformation(encounter.getUnboundComplexChainInterface(), encounter.getBoundComplexChainInterface())
				t_newcomplex.apply(encounter.unboundComplexAligned.complex)
				# copy coordinates from the moved complex onto proteins 1 and 2 in their individual representation
				encounter.unboundComplexAligned.copyCoordsFromComplex(encounter.unboundComplexAligned.complex, encounter.getReferenceSegment(), encounter.getReference())
				encounter.unboundComplexAligned.copyCoordsFromComplex(encounter.unboundComplexAligned.complex, encounter.getCounterpartSegment(), encounter.getUnboundCounterpart())
						  
			# first align alpha whole to create the complex, then align the complex again based on complexOnComplex for whole or interface
			elif alignStyle == "2bcomplex":                          
				# align protein1
				t = calcTransformation(encounter.getRefChain(), encounter.getMobChain())
				t.apply(encounter.getReference())
				# align protein2
				t_counterpart = calcTransformation(encounter.getUnboundCounterpartChain(), encounter.getBoundCounterpartChain())
				t_counterpart.apply(encounter.getUnboundCounterpart())   
				# copy coordinates from the moved proteins 1 and 2 on the complex representation
				encounter.unboundComplexAligned.complex = encounter.unboundComplexAligned.copyCoords(receptor, ligand, encounter.unboundComplexAligned.complex, checkIfStructuresAreClose=False)
				# now align based on complexOnComplex                          
				t_newcomplex = calcTransformation(encounter.getUnboundComplexChainInterface(), encounter.getBoundComplexChainInterface())
				t_newcomplex.apply(encounter.unboundComplexAligned.complex)
				# copy coordinates from the moved complex onto proteins 1 and 2 in their individual representation
				encounter.unboundComplexAligned.copyCoordsFromComplex(encounter.unboundComplexAligned.complex, encounter.getReferenceSegment(), encounter.getReference())
				encounter.unboundComplexAligned.copyCoordsFromComplex(encounter.unboundComplexAligned.complex, encounter.getCounterpartSegment(), encounter.getUnboundCounterpart())                          
						  
			# first align alpha whole to create the complex, then align the complex again based on beta of the reference segment for whole or interface
			elif alignStyle == "2bindividual":                          
				# align protein1
				t = calcTransformation(encounter.getRefChain(), encounter.getMobChain())
				t.apply(encounter.getReference())
				# align protein2
				t_counterpart = calcTransformation(encounter.getUnboundCounterpartChain(), encounter.getBoundCounterpartChain())
				t_counterpart.apply(encounter.getUnboundCounterpart())       
				# copy coordinates from the moved proteins 1 and 2 on the complex representation
				encounter.unboundComplexAligned.complex = encounter.unboundComplexAligned.copyCoords(receptor, ligand, encounter.unboundComplexAligned.complex, checkIfStructuresAreClose=False)
				# now align based on beta
				t_newcomplex = calcTransformation(encounter.getUnboundComplexChainInterface().select('segment \"'+encounter.getReferenceSegment()+'.\"'), 
												  encounter.getBoundComplexChainInterface().select('segment \"'+encounter.getReferenceSegment()+'.\"'))
				t_newcomplex.apply(encounter.unboundComplexAligned.complex)
				# copy coordinates from the moved complex onto proteins 1 and 2 in their individual representation
				encounter.unboundComplexAligned.copyCoordsFromComplex(encounter.unboundComplexAligned.complex, encounter.getReferenceSegment(), encounter.getReference())
				encounter.unboundComplexAligned.copyCoordsFromComplex(encounter.unboundComplexAligned.complex, encounter.getCounterpartSegment(), encounter.getUnboundCounterpart())                 
						  
			else:
				raise Exception("No such alignment know.")                          
						  
			# save the unbound complex
			if storeUnboundComplex:
				self.unboundComplexCopy = encounter.unboundComplexAligned.complex.copy()
			
			# assert equality of individual proteins and the complexes they form
			self.equalityAssertionsOfComplexAndItsProteins(encounter, checkInterface=True)
			
		# save all pdbs
		if storePDBs and dataHolder is not None:
			self.savePDBs(encounter, dataHolder)
			
	def setupNMA(self, encounter, dataHolder):
		""" Setup and calculate the ANMs with parameters based on the config object.
		
		Args:
			encounter: the encounter object
		"""
		assert self.config.investigationsOn == "Individual" or self.config.investigationsOn == "Complex"
		encounter.initANMs(self.utils)
		
		if self.config.investigationsOn == "Individual" or self.config.investigationsOn == "Complex":
			maxModes = min((encounter.getReference().select('calpha').numAtoms()*3 -6), self.config.maxModesToCalculate)
			if self.bound_provided == True:
				encounter.accessANMs().calcANMsUnified(encounter.getReference(), 
												 encounter.getUnboundCounterpart(), 
												 encounter.unboundComplexAligned.complex,
												 maxModes,
												 encounter,
												 encounter.getRefChain(),
												 encounter.getUnboundCounterpartChain(),
												 encounter.getUnboundComplexAlignedChain(),
												 selstr='calpha',
												 whatAtomsToMatch=self.config.whatAtomsToMatch)
			else:
				encounter.accessANMs().calcANMsUnified(encounter.getReference(), 
												 encounter.getUnboundCounterpart(), 
												 encounter.unboundComplexAligned.complex,
												 maxModes,
												 encounter,
												 selstr='calpha', 
												 whatAtomsToMatch=self.config.whatAtomsToMatch)	
			# Mode classifier, rescale eigenvalues based on complex eigenvalues
			if self.config.rescaleEigenvalues:
				modeClassifier = ModeClassifier(self.utils)
				anm_reference_rescaledEigenvals = modeClassifier.getRescaledANM(encounter.accessANMs().getANMComplex(), encounter.accessANMs().getANMReference(), encounter.getReference().getTitle())
				if self.bound_provided == True:
					encounter.accessANMs().replaceReferenceANMs(anm_reference_rescaledEigenvals, encounter.getReference(), encounter.getRefChain())       
				else:
					encounter.accessANMs().replaceReferenceANMs(anm_reference_rescaledEigenvals, encounter.getReference())     
		
		# elif self.config.investigationsOn == "Complex":
		# 	assert self.config.complexRMSDreduction == "1k1k" or self.config.complexRMSDreduction == "2k" or self.config.complexRMSDreduction == "1k1k6"
		# 	assert self.config.whichCustomHC == "HC_0" or self.config.whichCustomHC == "HC_U1" or self.utils.config.whichCustomHC == "HC_06" or self.utils.config.whichCustomHC == "HC_U1_1k1k"

		# 	maxModes = min((encounter.getReference().select('calpha').numAtoms()*3 -6), (encounter.getUnboundCounterpart().select('calpha').numAtoms()*3 -6), self.config.maxModesToCalculate)
		# 	if self.bound_provided == True:
		# 		encounter.accessANMs().calcANMsUnified(encounter.getReference(), 
		# 										 encounter.getUnboundCounterpart(), 
		# 										 encounter.unboundComplexAligned.complex,
		# 										 maxModes,
		# 										 encounter,
		# 										 encounter.getRefChain(),
		# 										 encounter.getUnboundCounterpartChain(), 
		# 										 encounter.getUnboundComplexAlignedChain(),
		# 										 selstr='calpha', 
		# 										 whatAtomsToMatch=self.config.whatAtomsToMatch)
		# 	else:
		# 		encounter.accessANMs().calcANMsUnified(encounter.getReference(), 
		# 										 encounter.getUnboundCounterpart(), 
		# 										 encounter.unboundComplexAligned.complex,
		# 										 maxModes,
		# 										 encounter,
		# 										 selstr='calpha', 
		# 										 whatAtomsToMatch=self.config.whatAtomsToMatch)				
			
			# # Mode classifier, rescale eigenvalues 
			# if self.config.rescaleEigenvalues:     
			# 	print "re-ranking eigenvalues of HC"       
			# 	modeClassifier = ModeClassifier(self.utils)
			# 	anm_complex_rescaledEigenvals, indicesOfLambdaRSorting = modeClassifier.getRescaledComplexANM(encounter.accessANMs().getANMComplex(), encounter.getUnboundCounterpart())          
			# 	# replace the complex ANM
			# 	if self.bound_provided == True:
			# 		encounter.accessANMs().replaceComplexANMs(anm_complex_rescaledEigenvals, encounter.unboundComplexAligned.complex, encounter.getUnboundComplexAlignedChain())
			# 	else:
			# 		encounter.accessANMs().replaceComplexANMs(anm_complex_rescaledEigenvals, encounter.unboundComplexAligned.complex)

			# 	# put indices in the dataHolder
			# 	dataHolder.indicesOfLambdaRSorting = indicesOfLambdaRSorting
			
		return dataHolder
	
	def calculateMeasuresExpanding(self, encounter, dataHolder, proteinFrom, proteinTo, anmReference, anmCounterpart, anmComplex):
		""" Calculate the measures (RMSD reduction, cumul overlap, etc.) on the whole proteins. 
		
		Args:
			encounter: the encounter object
			dataHolder: holds the direct results
			
		Returns:
			dataHolder: dataHolder with obtained data
		"""
		RMSDreducer = RMSDReducer(self.utils)
		RMSD_unbound_to_superposed_bound = calcRMSD(proteinFrom, proteinTo)
		# If the RMSD is very close to zero, do not consider this pair  
		if RMSD_unbound_to_superposed_bound < self.config.floatingPointThreshold:
			print "RMSD_unbound_to_superposed_bound < self.config.floatingPointThreshold, skipping ", encounter._pdbQueueItem[0]
			sys.exit()            
		defvec = calcDeformVector(proteinFrom, proteinTo)      
		# stepPointsReduction as the list of numbers of non-trivial normal modes to be used for fitting
		stepPointsReduction = self.getStepPointsReduction(self.config.investigationsOn, encounter)
		
		if self.config.investigationsOn == "Individual":
			# RMSD reductions on the receptor using the modes from Marray
			if (self.config.calculateZeroEigvalModes == True) and (self.config.rescaleEigenvalues == False):
				Marray = anmReference[0].getArray().T[6:].T
			else:
				Marray = anmReference[0].getArray()

		elif self.config.investigationsOn == "Complex":           
						
			if self.config.complexRMSDreduction == "2k":
				# RMSD reductions on whole complex with 2k modes
				if self.config.calculateZeroEigvalModes == True:
					if self.config.whichCustomHC == "HC_0": # remove 12 modes
						Marray = anmComplex[0].getArray().T[12:].T
					elif self.config.whichCustomHC == "HC_U1": # remove 6 modes
						Marray = anmComplex[0].getArray().T[6:].T                      
					elif self.config.whichCustomHC == "HC_06": # HC_O and 6 trival modes of the ligand
						Marray = anmComplex[0].getArray().T[12:].T
						counterpartZeroEigenvalModes = anmCounterpart[0].getArray().T[:6].T
						counterpartZeroEigenvalModes = np.pad(counterpartZeroEigenvalModes, ((Marray.shape[0]-counterpartZeroEigenvalModes.shape[0],0),(0,0)), mode='constant')
						Marray = np.column_stack((counterpartZeroEigenvalModes, Marray))
						Marray = Marray.T[:-6].T # remove last 6 modes at the end of the array
					elif self.config.whichCustomHC == "HC_U1_1k1k": # calculate modes from HC, then split them into 1k 1k component vectors
						Marray = anmComplex[0].getArray().T[6:].T                      
						HRcomponent = Marray[:anmReference[0].getArray().shape[0]] 
						HLcomponent = Marray[anmReference[0].getArray().shape[0]:]
						assert HRcomponent.shape[0] == proteinFrom.select('segment "R."').numAtoms()*3
						assert HLcomponent.shape[0] == proteinFrom.select('segment "L."').numAtoms()*3
						assert HRcomponent.shape[0] + HLcomponent.shape[0] == Marray.shape[0]
						Marray = self.utils.makeBlockZeroMatrix(HRcomponent, HLcomponent)
						Marray = Marray.T[:anmComplex[0].numModes()-6].T # remove last modes in Marray to match the same number of modes as have been calculated                       
						# normalize the Marray
						Marray = encounter.accessANMs().normalizeM(Marray)
						# update the anm with these modes
						anmComplex[0].setEigens(Marray, anmComplex[0].getEigvals()[6:])
						if self.config.measuresOnWhole == True:
							encounter.accessANMs()._anm_complex_slc = (anmComplex[0], encounter.accessANMs()._anm_complex_slc[1])
						else:
							encounter.accessANMs()._anm_boundcomplex_slc_interface = (anmComplex[0], encounter.accessANMs()._anm_boundcomplex_slc_interface[1])
				else:
					Marray = anmComplex[0].getArray()           
			elif self.config.complexRMSDreduction == "1k1k":
				# RMSD reductions on whole complex used 1k modes from receptor & 1k modes from ligand (pairwise)
				if self.config.calculateZeroEigvalModes == True:
					Marray = self.utils.makeBlockZeroMatrix(anmReference[0].getArray().T[6:].T, anmCounterpart[0].getArray().T[6:].T)
				else:
					Marray = self.utils.makeBlockZeroMatrix(anmReference[0].getArray(), anmCounterpart[0].getArray())
			elif self.config.complexRMSDreduction == "1k1k6":
				assert self.config.calculateZeroEigvalModes == True
				Marray = self.utils.makeBlockZeroMatrix(anmReference[0].getArray().T[6:].T, anmCounterpart[0].getArray().T[6:].T)
				counterpartZeroEigenvalModes = anmCounterpart[0].getArray().T[:6].T
				counterpartZeroEigenvalModes = np.pad(counterpartZeroEigenvalModes, ((Marray.shape[0]-counterpartZeroEigenvalModes.shape[0],0),(0,0)), mode='constant')
				Marray = np.column_stack((counterpartZeroEigenvalModes, Marray))
				Marray = Marray.T[:-6].T # remove last 6 modes at the end of the array        

		RMSDReductions, overlapTApprox, stepPointsReduction, L_RMSReductions, deformationSnapshots = RMSDreducer.calcRMSDReductionsExpandingSet(Marray, 
																										   proteinTo,
																										   proteinFrom,  
																										   defvec,
																										   stepPointsReduction, 
																										   encounter._pdbQueueItem[0],
																										   self.config.experimentNamePrefix)            
		assert self.assessment.non_increasing(RMSDReductions.tolist())
		assert self.assessment.non_increasing(L_RMSReductions.tolist())
		print "RMSD before         : " + str(RMSD_unbound_to_superposed_bound)
		print "RMSDReductions shape, tolist: ", RMSDReductions.shape, RMSDReductions.tolist()
		
		# save results to the dataholder
		dataHolder.L_RMSReductions = L_RMSReductions
		dataHolder.L_RMSD_unbound_to_superposed_bound = RMSDreducer.getL_RMS(proteinFrom, proteinTo, self.config.investigationsOn)
		if self.config.measuresOnWhole == True:
			dataHolder.RMSD_unbound_to_superposed_bound = RMSD_unbound_to_superposed_bound
			dataHolder.RMSDReductionsWhole = RMSDReductions
			dataHolder.overlapTApproxWhole = overlapTApprox
			dataHolder.stepPointsReductionWhole = stepPointsReduction
			# capture overlap, collectivity, correlation
			if self.config.investigationsOn == "Individual":
				dataHolder.overlapArrayWhole = self.assessment.calcOverlapArrayWithProDyOverlap(anmReference[0], defvec)
				dataHolder.cumulOverlapWholePrody = calcCumulOverlap(anmReference[0], defvec)
				dataHolder.collectivityArrayWhole = self.assessment.calcCollectivityArray(anmReference[0])
				dataHolder.correlationArrayWhole =  self.assessment.calcCorrelationArray(anmReference[0], defvec)              
		else:
			dataHolder.RMSD_interface = RMSD_unbound_to_superposed_bound
			dataHolder.RMSDReductionsInterface = RMSDReductions
			dataHolder.overlapTApproxInterface = overlapTApprox
			dataHolder.stepPointsReductionInterface = stepPointsReduction
			# capture overlap, collectivity, correlation
			if self.config.investigationsOn == "Individual":
				dataHolder.overlapArrayInterface = self.assessment.calcOverlapArrayWithProDyOverlap(anmReference[0], defvec)
				dataHolder.cumulOverlapInterfacePrody = calcCumulOverlap(anmReference[0], defvec)
				dataHolder.collectivityArrayInterface = self.assessment.calcCollectivityArray(anmReference[0])
				dataHolder.correlationArrayInterface =  self.assessment.calcCorrelationArray(anmReference[0], defvec)   
				
		# save individual overlap of non-trivial modes:
		dataHolder.singleModeOverlapsFromSuperset = self.utils.getModeOverlaps(Marray, defvec)
		
		# save deformation snapshots, to show how the protein deformes with the application of modes
		dataHolder.deformationSnapshots = deformationSnapshots
			
		return dataHolder    

	def calculateMeasuresOuterLoop(self, encounter, dataHolder):
		""" Calculate the measures (RMSD reduction, cumul overlap, etc.), the creation of the modes array is performed
		in an outer loop prior to giving it to the fitter.
		
		Args:
			encounter: the encounter object
			dataHolder: holds the direct results
			
		Returns:
			dataHolder: dataHolder with obtained data
		"""
		if self.config.measuresOnWhole == True:
			if self.config.investigationsOn == "Complex": 
				print "whole"
				proteinFrom = encounter.getUnboundComplexAlignedChain()
				proteinTo = encounter.getBoundComplexChain()
				anmReference = encounter.accessANMs().getANMReference2a2kSlc()
				anmCounterPart = encounter.accessANMs().getANMCounterpart2a2kSlc()
		else:
			if self.config.investigationsOn == "Complex":
				print "interface"
				proteinFrom = encounter.getUnboundComplexChainInterface()
				proteinTo = encounter.getBoundComplexChainInterface()
				anmReference = encounter.accessANMs()._anm_reference_slc_interface
				anmCounterPart = encounter.accessANMs()._anm_counterpart_slc_interface
				
		RMSD_unbound_to_superposed_bound = calcRMSD(proteinFrom, proteinTo)
		# If the RMSD is very close to zero, do not consider this pair  
		if RMSD_unbound_to_superposed_bound < self.config.floatingPointThreshold:
			print "RMSD_unbound_to_superposed_bound < self.config.floatingPointThreshold, skipping ", encounter._pdbQueueItem[0]
			sys.exit()            
		defvec = calcDeformVector(proteinFrom, proteinTo)
			
		# stepPointsReduction as the list of numbers of non-trivial normal modes to be used for fitting
		stepPointsReduction = self.getStepPointsReduction(self.config.investigationsOn, encounter)
		print stepPointsReduction

		# create Marray, submit to fitter
		RMSDreducer = RMSDReducer(self.utils)
		RMSDReductions = []
		overlaps = []
		allBetas = []
		for element in stepPointsReduction:
			if self.config.complexRMSDreduction == "2k":
				# create Marray, submit to fitter
				print "2k"
			elif self.config.complexRMSDreduction == "1k1k":
				# create Marray
				element /= 2
				Marray = self.utils.makeBlockZeroMatrix(anmReference[0].getArray().T[6:6+element].T, anmCounterPart[0].getArray().T[6:6+element].T)      
				print Marray.shape
				# submit Marray to fitter             
				RMSD_after_Tapprox, currentOverlap, betas = RMSDreducer.calcRMSDReductionFromTo(Marray, proteinFrom, proteinTo, defvec, allBetas, overlaps, RMSDReductions, proteinFrom.getTitle(), "RMSDReductionFixedset")
			elif self.config.complexRMSDreduction == "1k1k6":
				element /= 2
				numRows = anmReference[0].getArray().shape[0] + anmCounterPart[0].getArray().shape[0]
				# zero eigenvalue modes of the ligand
				counterpartZeroEigenvalModes = anmCounterPart[0].getArray().T[:6].T
				counterpartZeroEigenvalModes = np.pad(counterpartZeroEigenvalModes, ((numRows-counterpartZeroEigenvalModes.shape[0],0),(0,0)), mode='constant')
				# create Marray, submit to fitter
				if element != 0:
					Marray = self.utils.makeBlockZeroMatrix(anmReference[0].getArray().T[6:6+element].T, anmCounterPart[0].getArray().T[6:6+element].T)
					Marray = np.column_stack((counterpartZeroEigenvalModes, Marray))   
				else:
					Marray = counterpartZeroEigenvalModes
				print Marray.shape
				# submit Marray to fitter
				RMSD_after_Tapprox, currentOverlap, betas = RMSDreducer.calcRMSDReductionFromTo(Marray, proteinFrom, proteinTo, defvec, allBetas, overlaps, RMSDReductions, proteinFrom.getTitle(), "RMSDReductionFixedset")
			# append the obtained data in lists
			RMSDReductions.append(RMSD_after_Tapprox)
			overlaps.append(currentOverlap)
			allBetas.append(betas) 
		print RMSDReductions
		print overlaps
		dataHolder.RMSDReductionsWhole = RMSDReductions
		dataHolder.overlapTApproxWhole = currentOverlap
		dataHolder.stepPointsReductionWhole = stepPointsReduction
		raw_input() 
		
	def getStepPointsReduction(self, reductionStyle, encounter):
		""" Return a list of number of non-trivial modes to be used for linear combinations.
		
		Args:
			reductionStyle: string describing the reductionstyle
			encounter: the encounter object
		
		Returns: stepPointsReduction, a list of number of modes to be used for linear combinations 
		"""
		if reductionStyle == "Complex":
			if self.config.complexRMSDreduction == "2k" and self.config.whichCustomHC != "HC_0" and self.config.whichCustomHC != "HC_06":
				maxModes = min(encounter.accessANMs()._anm_complex[0].numModes()-6, self.config.stopRMSDReductionAt)
				return np.concatenate((np.arange(1, 21, 1), np.arange(30, maxModes+10, 10)))       
			elif self.config.complexRMSDreduction == "1k1k" or (self.config.whichCustomHC == "HC_0" and self.config.complexRMSDreduction == "2k") or (self.config.whichCustomHC == "HC_06" and self.config.complexRMSDreduction == "2k"):
				maxModes = min(encounter.accessANMs()._anm_reference[0].numModes()-6 + encounter.accessANMs()._anm_counterpart[0].numModes()-6, self.config.stopRMSDReductionAt)
				return np.concatenate((np.arange(2, 22, 2), np.arange(30, maxModes+10, 10)))
			elif self.config.complexRMSDreduction == "1k1k6":
				maxModes = min(encounter.accessANMs()._anm_reference[0].numModes()-6 + encounter.accessANMs()._anm_counterpart[0].numModes()-6, self.config.stopRMSDReductionAt)
				return (np.concatenate((np.arange(0, 22, 2), np.arange(24, maxModes+10-6, 10))) +6)
		elif reductionStyle == "Individual":
			maxModes = min(encounter.accessANMs()._anm_reference[0].numModes()-6, self.config.stopRMSDReductionAt)
			return np.concatenate((np.arange(1, 21, 1), np.arange(30, maxModes+10, 10)))
		else: 
			raise ValueError("reductionStyle is not recognized in the method getStepPointsReduction")
		
	def setupInterface(self, encounter):
		""" Setup the interface of the proteins. 
		
		Args:
			encounter: the encounter object
		"""
		print "interface setup"
		# calculate interface independently on mobile and its bound counterpart
		encounter.calcMobChainInterface()
		encounter.calcRefChainInterface(self.utils)
		encounter.calcBoundCounterpartChainInterface()
		encounter.setUnboundCounterpartInterface(encounter.calcUnboundInterface(encounter.getBoundCounterpartChain(), 
																				encounter.getBoundCounterpartChainInterface(), 
																				encounter.getUnboundCounterpartChain(), 
																				self.utils))
		
		# calculate interface directly from boundComplex
		encounter.calcBoundComplexChainInterface()
		encounter.setUnboundComplexChainInterface(encounter.calcUnboundInterface(encounter.getBoundComplexChain(), 
																				 encounter.getBoundComplexChainInterface(), 
																				 encounter.getUnboundComplexAlignedChain(), 
																				 self.utils))
		
		# check equality of parsed proteins and their interfaces
		self.equalityAssertionsOfComplexAndItsProteins(encounter)
	
	def rotateM(self, M, rotationMatrix):
		""" Rotate a matrix M, consisting of columnwise normal modes, by a 3*3 rotationMatrix. The result is a matrix M
		with the normal modes rotated according to the rotationMatrix.
		
		Args:
			M: matrix with columnwise normal modes
			rotationMatrix: 3*3 rotationMatrix
			
		Returns:
			Mrotated: matrix with columnwise normal modes, which were rotated through the rotationMatrix
		"""
		assert rotationMatrix.shape == (3,3)
		#print "shape M: ", M.shape
		#print "shape M.T[:1]", M.T[:1].T.shape
		# construct a 3n*3n rotation matrix
		R = np.tile(rotationMatrix, (M.shape[0]/3, M.shape[0]/3))
		#print "shape R: ", R.shape
		# apply rotation
		Mrot = np.dot(R, M)
		writeArray("Mrot.txt", Mrot, format='%f')
		writeArray("R.txt", R, format='%f')
		writeArray("M.txt", M, format='%f')
		writeArray("MT[1]T.txt", M.T[:1].T, format='%f')
		return Mrot        
	
	def savePDBs(self, encounter, dataHolder):
		""" After aligning the proteins to the formations where the measures are collected, and setting up the chain matchings, interfaces 
		and the complex representations, this method stores a copy of them. 
		
		Args:
			encounter: the encounter object
			dataholder: holds the direct results
		"""
		dataHolder.reference = encounter.getReference().copy()
		dataHolder.mobile = encounter.getMobile().copy()
		dataHolder.unboundCounterpart = encounter.getUnboundCounterpart().copy()
		dataHolder.boundCountertpart = encounter.getBoundCounterpart().copy()
		dataHolder.unboundComplex = encounter.unboundComplexAligned.complex.copy()
		dataHolder.boundComplex = encounter.boundComplex.complex.copy()
		dataHolder.refChain = encounter.getRefChain().copy()
		dataHolder.refChainInterface = encounter.getRefChainInterface().copy()
		dataHolder.mobChain = encounter.getMobChain().copy()
		dataHolder.mobChainInterface = encounter.getMobChainInterface().copy()
		dataHolder.unboundCounterpartChain = encounter.getUnboundCounterpartChain().copy()
		dataHolder.unboundCounterpartChainInterface = encounter.getUnboundCounterpartChainInterface().copy()
		dataHolder.boundCounterpartChain = encounter.getBoundCounterpartChain().copy()
		dataHolder.boundCounterpartChainInterface = encounter.getBoundCounterpartChainInterface().copy()
		dataHolder.unboundComplexAlignedChain = encounter.getUnboundComplexAlignedChain().copy()
		dataHolder.unboundComplexChainInterface = encounter.getUnboundComplexChainInterface().copy()
		dataHolder.boundComplexChain = encounter.getBoundComplexChain().copy()
		dataHolder.boundComplexChainInterface = encounter.getBoundComplexChainInterface().copy()
				
	def calculateMeasuresDistance(self, encounter, dataHolder):
		""" Distance measures (I-rms, L-rms) following the idea of 
		"Docking and scoring protein complexes: CAPRI 3rd Edition" using backbone atoms of 
		the ligand (for L-RMS) or all interface (combined receptor and ligand interfaces) backbone atoms (for I-RMS)
		
		Args:
			encounter: the encounter object
			dataHolder: holds the direct results
			
		Returns:
			dataHolder: dataHolder with obtained data        
		"""
		# I-RMS
		self.alignProteins(encounter, "complexOnComplex", False, storeUnboundComplex=False)
		dataHolder.I_rms_after_align = calcRMSD(encounter.getUnboundComplexChainInterface(), encounter.getBoundComplexChainInterface())
		# L-RMS
		self.alignProteins(encounter, "L-RMS", True, storeUnboundComplex=False)
		if self.utils.isReceptor(encounter.getReference().getTitle()):
			dataHolder.L_rms = calcRMSD(encounter.getUnboundCounterpartChain(), encounter.getBoundCounterpartChain())
		else:
			dataHolder.L_rms = calcRMSD(encounter.getRefChain(), encounter.getMobChain())
		return dataHolder
	
	def outputResults(self, encounter, dataHolder, outputPath=None):
		""" Output the obtained results. 
		
		Args:
			encounter: the encounter object
			dataHolder: holds the direct results 
			
		Returns:
			resultsPath: absolute path to the results folder
		"""
		# store the direct results

		counterpartTitle = encounter.getUnboundCounterpart().getTitle()
		if self.config.investigationsOn == "Complex":
			encounter.initResultsPrinter(counterpartTitle[:4]+counterpartTitle[counterpartTitle.rfind("_"):])
		else:
			if self.utils.isReceptor(encounter.getReference().getTitle()):
				outputTitle = encounter.getUnboundCounterpart().getTitle()
				outputTitle = outputTitle.replace("_l_", "_r_")
				# print outputTitle
				# assert outputTitle[-3] == "r"
			else:
				outputTitle = encounter.getReference().getTitle()
			encounter.initResultsPrinter(outputTitle)
		if self.bound_provided == True:
			if self.config.measuresOnWhole == True:
				encounter.resultsPrinter.setRMSDReductionsWhole(dataHolder.RMSDReductionsWhole)
				encounter.resultsPrinter.setOverlapTApproxWhole(dataHolder.overlapTApproxWhole)
				encounter.resultsPrinter.setStepPointsReductionWhole(dataHolder.stepPointsReductionWhole)
				encounter.resultsPrinter.setRMSD_unbound_to_superposed_bound(dataHolder.RMSD_unbound_to_superposed_bound)
			if self.config.measuresOnWhole == False:
				encounter.resultsPrinter.setRMSDReductionsInterface(dataHolder.RMSDReductionsInterface)
				encounter.resultsPrinter.setOverlapTApproxInterface(dataHolder.overlapTApproxInterface)
				encounter.resultsPrinter.setStepPointsReductionInterface(dataHolder.stepPointsReductionInterface)
				encounter.resultsPrinter.setRMSD_interface(dataHolder.RMSD_interface)
			if self.config.investigationsOn == "Individual" or self.config.investigationsOn == "Complex":
				encounter.resultsPrinter.setEigenvectorsReference(encounter.accessANMs().getANMReference()[0].getArray().T)			
				encounter.resultsPrinter.setEigenvaluesReference(encounter.accessANMs().getANMReferenceSlc()[0].getEigvals())
				encounter.resultsPrinter.setEigenvaluesComplex(encounter.accessANMs().getANMComplexSlc()[0].getEigvals())
				encounter.resultsPrinter.setEigenvectorsComplex(encounter.accessANMs().getANMComplexSlc()[0].getArray().T)

				# eigenvectors output
			
			# encounter.resultsPrinter.setEigenvaluesReceptor1k(encounter.accessANMs()._anm_reference_slc[0].getEigvals())

			# encounter.resultsPrinter.setEigenvaluesLigand1k(encounter.accessANMs()._anm_counterpart[0].getEigvals())        
			encounter.resultsPrinter.setNumberOfModes(encounter.accessANMs().getANMReferenceSlc()[0].numModes())
			encounter.resultsPrinter.setNumberOfModesComplex(encounter.accessANMs().getANMComplexSlc()[0].numModes())
			encounter.resultsPrinter.setPathOfConfigFile(self.pathOfConfigFile)
			# measures
			if self.config.investigationsOn == "Individual":
				if self.config.measuresOnWhole == True:
					encounter.resultsPrinter.setOverlapArrayWhole(dataHolder.overlapArrayWhole)
					encounter.resultsPrinter.setCumulOverlapWholePrody(dataHolder.cumulOverlapWholePrody)
					encounter.resultsPrinter.setCollectivityArrayWhole(dataHolder.collectivityArrayWhole)
					encounter.resultsPrinter.setCorrelationArrayWhole(dataHolder.correlationArrayWhole)
				else:
					encounter.resultsPrinter.setOverlapArrayInterface(dataHolder.overlapArrayInterface)
					encounter.resultsPrinter.setCumulOverlapInterfacePrody(dataHolder.cumulOverlapInterfacePrody)
					encounter.resultsPrinter.setCollectivityArrayInterface(dataHolder.collectivityArrayInterface)
					encounter.resultsPrinter.setCorrelationArrayInterface(dataHolder.correlationArrayInterface)
			if self.config.calculateZeroEigvalModes == True:
				encounter.resultsPrinter.setZeroEigvecsProtein1(encounter.accessANMs()._anm_reference[0].getArray().T[:6])
				encounter.resultsPrinter.setZeroEigvecsProtein2(encounter.accessANMs()._anm_counterpart[0].getArray().T[:6])
				if self.config.whichCustomHC == "HC_0" or self.config.whichCustomHC == "HC_06":
					encounter.resultsPrinter.setZeroEigvecsComplex(encounter.accessANMs()._anm_complex[0].getArray().T[:12])
				else:
					encounter.resultsPrinter.setZeroEigvecsComplex(encounter.accessANMs()._anm_complex[0].getArray().T[:6])
			# dataHolder.indicesOfLambdaRSorting
			if self.config.investigationsOn == "Complex":
				if self.config.rescaleEigenvalues:
					encounter.resultsPrinter.setIndicesOfLambdaRSorting(dataHolder.indicesOfLambdaRSorting)
			# L-RMS reduction
			encounter.resultsPrinter.setL_RMSReductions(dataHolder.L_RMSReductions)
			encounter.resultsPrinter.setL_RMSD_unbound_to_superposed_bound(dataHolder.L_RMSD_unbound_to_superposed_bound)
			# L-RMS and IRMS
			encounter.resultsPrinter.setLRMS(dataHolder.L_rms)
			encounter.resultsPrinter.setIRMSafterAlign(dataHolder.I_rms_after_align)
			# Overlap from Marray superset
			encounter.resultsPrinter.setSingleModeOverlapsFromSuperset(dataHolder.singleModeOverlapsFromSuperset)
			# Deformation snapshots, set of pdbs
			encounter.resultsPrinter.setDeformationSnapshots(dataHolder.deformationSnapshots)
			# PDBs
			try:
				encounter.resultsPrinter.setReference(dataHolder.reference)
				encounter.resultsPrinter.setRefChain(dataHolder.refChain)
				encounter.resultsPrinter.setRefChainInterface(dataHolder.refChainInterface)
				encounter.resultsPrinter.setMobile(dataHolder.mobile)
				encounter.resultsPrinter.setMobChain(dataHolder.mobChain)
				encounter.resultsPrinter.setMobChainInterface(dataHolder.mobChainInterface)
				encounter.resultsPrinter.setUnboundCounterpart(dataHolder.unboundCounterpart)
				encounter.resultsPrinter.setUnboundCounterpartChain(dataHolder.unboundCounterpartChain)
				encounter.resultsPrinter.setUnboundCounterpartChainInterface(dataHolder.unboundCounterpartChainInterface)
				encounter.resultsPrinter.setBoundCounterpart(dataHolder.boundCountertpart)
				encounter.resultsPrinter.setBoundCounterpartChain(dataHolder.boundCounterpartChain)
				encounter.resultsPrinter.setBoundCounterpartChainInterface(dataHolder.boundCounterpartChainInterface)
				encounter.resultsPrinter.setUnboundComplexAligned(dataHolder.unboundComplex)
				encounter.resultsPrinter.setUnboundComplexAlignedChain(dataHolder.unboundComplexAlignedChain)
				encounter.resultsPrinter.setUnboundComplexChainInterface(dataHolder.unboundComplexChainInterface)
				encounter.resultsPrinter.setBoundComplex(dataHolder.boundComplex)
				encounter.resultsPrinter.setBoundComplexChain(dataHolder.boundComplexChain)
				encounter.resultsPrinter.setBoundComplexChainInterface(dataHolder.boundComplexChainInterface)
			except Exception:
				pass
		else:
			if self.config.investigationsOn == "Individual" or self.config.investigationsOn == "Complex":
				encounter.resultsPrinter.setEigenvectorsReference(encounter.accessANMs().getANMReference()[0].getArray().T)
				encounter.resultsPrinter.setEigenvaluesReference(encounter.accessANMs().getANMReferenceSlc()[0].getEigvals())
				encounter.resultsPrinter.setEigenvectorsComplex(encounter.accessANMs().getANMComplexSlc()[0].getArray().T)
				encounter.resultsPrinter.setEigenvaluesComplex(encounter.accessANMs().getANMComplexSlc()[0].getEigvals())

			encounter.resultsPrinter.setNumberOfModes(encounter.accessANMs().getANMReferenceSlc()[0].numModes())
			encounter.resultsPrinter.setPathOfConfigFile(self.pathOfConfigFile)
			if self.config.calculateZeroEigvalModes == True:
				encounter.resultsPrinter.setZeroEigvecsProtein1(encounter.accessANMs()._anm_reference[0].getArray().T[:6])
				encounter.resultsPrinter.setZeroEigvecsProtein2(encounter.accessANMs()._anm_counterpart[0].getArray().T[:6])
				if self.config.whichCustomHC == "HC_0" or self.config.whichCustomHC == "HC_06":
					encounter.resultsPrinter.setZeroEigvecsComplex(encounter.accessANMs()._anm_complex[0].getArray().T[:12])
				else:
					encounter.resultsPrinter.setZeroEigvecsComplex(encounter.accessANMs()._anm_complex[0].getArray().T[:6])
			
		# write the direct results 
		if outputPath:
			# if outputPath is ., make sure it refers to the directory where the program is located
			if outputPath == ".":
				outputPath = os.path.dirname(os.path.realpath(__file__))
			outputPath = makeStringEndWith(outputPath, "/")
			self.config.outputPath = outputPath     
		resultsPath = encounter.resultsPrinter.writeDirectResults(self.config.outputPath, self.config.experimentNamePrefix, self.utils, self.bound_provided, self.config.investigationsOn)
		encounter.resultsPrinter.writeNMDResults(self.config.outputPath, self.config.experimentNamePrefix, self.utils, encounter.accessANMs(), single = False, storeANMs=False)
		
		return resultsPath
		
	def mainCalculation(self, protein1_A, protein2_A=None, protein1_B=None, protein2_B=None, outputPath=None):
		""" Perform the NMA based on the provided proteins and the Configuration. 
		
		Args:
			protein1_A: Protein 1 in first conformational state
			protein2_A: Protein 2 in first conformational state
			protein1_B: Protein 1 in second conformational state
			protein2_B: Protein 2 in second conformational state
			outputPath: Optional argument, specifying the output path, ignoring the output path from the Configuration 
			
		Returns:
			resultsPath where the results have been written to inside the specified output folder
			
		Results: 
			Writes experimental results to specified output folder
			
		Example scenarios:
			2a    A=bound, B=unbound
			2b    A=unbound, B=bound
			2c    A=docked, B=bound, k=1..k, loop over the A conformations outside of this method
		
		Note that for now it is assumed that when doing conformational change calculations from whole unbound to bound complexes,
		protein_1 is the receptor and protein_2 the ligand. This can be generalized in the future. 
		"""
		print self.pathOfConfigFile

		# if only 1 protein provided with two conformations, do canonical NMA from unbound to bound
		if protein2_A == None:
			print "NNA with 2 proteins present, canonical NMA" 
			print "Set Configurations to canonical settings for canonical NMA"
			assert os.path.isfile(protein1_A)
			print protein1_A
			self.bound_provided = False	
			encounter = Encounter((os.path.splitext(os.path.basename(protein1_A))[0], protein1_A))
			encounter.setReference(self.utils.parseFilteredPDB(protein1_A, self.config.filterPDB))	
			maxModes = min((encounter.getReference().select('all').numAtoms()*3 -6), self.config.maxModesToCalculate)
			encounter.initANMs(self.utils)
			encounter.accessANMs()._calcANMsUnified(encounter.getReference(), maxModes, 'calpha', 'calpha', True)
			outputTitle = encounter.getReference().getTitle()
			encounter.initResultsPrinter(outputTitle)
			encounter.resultsPrinter.setEigenvectorsReference(encounter.accessANMs().getANMReference()[0].getArray().T[6:])
			encounter.resultsPrinter.setEigenvaluesReference(encounter.accessANMs().getANMReferenceSlc()[0].getEigvals()[6:])
			encounter.resultsPrinter.setNumberOfModes(encounter.accessANMs().getANMReferenceSlc()[0].numModes())
			encounter.resultsPrinter.setPathOfConfigFile(self.pathOfConfigFile)
			# encounter.resultsPrinter.setRMSD_unbound_to_superposed_bound(calcRMSD(encounter.getReference()))		
			# write the direct results 
			if outputPath:
				# if outputPath is ., make sure it refers to the directory where the program is located
				if outputPath == ".":
					outputPath = os.path.dirname(os.path.realpath(__file__))
				outputPath = makeStringEndWith(outputPath, "/")
				self.config.outputPath = outputPath          
			resultsPath = encounter.resultsPrinter.writeDirectResults(self.config.outputPath, self.config.experimentNamePrefix, self.utils, self.bound_provided, self.config.investigationsOn)
			encounter.resultsPrinter.writeNMDResults(self.config.outputPath, self.config.experimentNamePrefix, self.utils, encounter.accessANMs(), single = True, storeANMs=False)
			
			outputTitle = encounter.getReference().getTitle()
		# 	RMSDprediction = RMSDPrediction(protein1_A)
		# 	RMSDPredicted = RMSDprediction.prediction(resultsPath, RMSDprediction.modelinput(), RMSDprediction.calculation(RMSDprediction.eigeninput(resultsPath), outputTitle))
			

		# 	# sampling

		# 	proteinFrom = encounter.getReference()
		# 	anmReferenceTemp = encounter.accessANMs().getANMReference()
		# 	anmReference = extendModel(anmReferenceTemp[0], anmReferenceTemp[1], proteinFrom, norm=True)
		# 	ensem = sampleModes(anmReference[0][6:107], proteinFrom, 100, RMSDPredicted)
		# 	# RMSDtostart = []
		# 	# RMSDtobound = []
		# 	# print type(proteinFrom)
		# 	# print np.unique(proteinFrom.select('calpha').getIndices())
		# 	ensem.setAtoms(proteinFrom)
  #  # 			for sample_index in range (0, 100):
  #  # 				RMSDtostart.append(calcRMSD(ensem.getConformation(sample_index).getCoords(), proteinFrom))
  #  # 				if self.bound_provided == True:
		# # 			RMSDtobound.append(calcRMSD(ensem.getConformation(sample_index).getCoords(), proteinTo))
		# 	# proteinFromTemp = proteinFrom.copy()	
		# 	# proteinFromTemp.addCoordset(ensem)
		# 	encounter.resultsPrinter.setEnsemble(ensem)
		# 	# encounter.resultsPrinter.setRMSDtostart(RMSDtostart)
		# 	# encounter.resultsPrinter.setRMSDtobound(RMSDtobound)

		# 	encounter.resultsPrinter.setRMSDPrediction(RMSDPredicted)
		# 	encounter.resultsPrinter.writeRMSDResults(self.config.outputPath, self.config.experimentNamePrefix, self.utils)
		# 	encounter.resultsPrinter.writeSampleResults(self.config.outputPath, self.config.experimentNamePrefix, self.utils)


			print "Results have been written to: ", resultsPath
			return resultsPath

		else:
			# if only 2 unbound protein provided, do NMA without post analysis
			if protein1_B == None and protein2_B == None:
				assert os.path.isfile(protein1_A)
				assert os.path.isfile(protein2_A)
				print "NMA without post analysis"
				self.bound_provided = False
			else:
				assert os.path.isfile(protein1_A)
				assert os.path.isfile(protein1_B)  
				assert os.path.isfile(protein2_A)
				assert os.path.isfile(protein2_B)
				print "NMA with 4 proteins present"
				self.bound_provided = True

			# load the protein structures into the encounter object and match their chains
			encounter = self.setupEncounter(protein1_A, protein2_A, protein1_B, protein2_B, self.affirmProteinNames)
			dataHolder = DataHolderForDirectOutput(encounter._pdbQueueItem)
			if self.bound_provided == True:
				self.chainMatching(encounter)
				# setup the interfaces, even for investigations on whole proteins this is needed to obtain for instance I-RMS values
				self.setupInterface(encounter)
				# align the proteins
				self.alignProteins(encounter, self.config.align, self.config.measuresOnWhole, dataHolder, storeUnboundComplex=True, storePDBs=True)
			# create the NMAs
			dataHolder = self.setupNMA(encounter, dataHolder) 
			#  if interface investigations will be performed, setup the NMA for the interfaces
			if self.config.measuresOnWhole == False:     
				encounter.accessANMs().calcInterfaceANMsUnified(encounter.getReference(), 
																	encounter.getUnboundCounterpart(),
																	encounter.unboundComplexAligned.complex,
																	encounter.getRefChainInterface(),
																	encounter.getUnboundCounterpartChainInterface(),
																	encounter.getUnboundComplexChainInterface())
			if self.bound_provided == True:
				# data collection on whole
				if self.config.measuresOnWhole == True:
					if self.config.investigationsOn == "Complex":
						dataHolder = self.calculateMeasuresExpanding(encounter, 
																		  dataHolder, 
																		  encounter.getUnboundComplexAlignedChain(), 
																		  encounter.getBoundComplexChain(),
																		  encounter.accessANMs().getANMReference2a2kSlc(),
																		  encounter.accessANMs().getANMCounterpart2a2kSlc(),
																		  encounter.accessANMs()._anm_complex_slc)
					elif self.config.investigationsOn == "Individual":
						dataHolder = self.calculateMeasuresExpanding(encounter, 
																	 dataHolder, 
																	 encounter.getRefChain(), 
																	 encounter.getMobChain(), 
																	 encounter.accessANMs()._anm_reference_slc,
																	 None, 
																	 encounter.accessANMs()._anm_complex_slc)
				# data collection on the interface
				if self.config.measuresOnWhole == False:     
					if self.config.investigationsOn == "Complex":
						dataHolder = self.calculateMeasuresExpanding(encounter, 
																		  dataHolder,
																		  encounter.getUnboundComplexChainInterface(),
																		  encounter.getBoundComplexChainInterface(),
																		  encounter.accessANMs()._anm_reference_slc_interface,
																		  encounter.accessANMs()._anm_counterpart_slc_interface,
																		  encounter.accessANMs()._anm_boundcomplex_slc_interface)
					elif self.config.investigationsOn == "Individual":
						dataHolder = self.calculateMeasuresExpanding(encounter, 
																	 dataHolder, 
																	 encounter.getRefChainInterface(), 
																	 encounter.getMobChainInterface(), 
																	 encounter.accessANMs()._anm_reference_slc_interface,
																	 None,
																	 encounter.accessANMs()._anm_boundcomplex_slc_interface)                      
				# data collection I-rms, L-rms
				dataHolder = self.calculateMeasuresDistance(encounter, dataHolder)
				if self.config.investigationsOn == "Complex":
					encounter.resultsPrinter.setC_RMSD_unbound_to_superposed_bound(calcRMSD(encounter.getUnboundComplexAlignedChain(), encounter.getBoundComplexChain()))
					encounter.resultsPrinter.writeComplexRMSD(self.config.outputPath, self.config.experimentNamePrefix, self.utils, self.bound_provided)
			

		
			# output
			resultsPath = self.outputResults(encounter, dataHolder, outputPath)
			outputTitle = encounter.getReference().getTitle()
			RMSDprediction = RMSDPrediction(protein1_A, protein2_A, self.config.investigationsOn)
			RMSDPredicted = RMSDprediction.prediction(resultsPath, RMSDprediction.modelinput(), RMSDprediction.calculation(RMSDprediction.eigeninput(resultsPath), outputTitle))
			

			# sampling
			if self.config.investigationsOn == "Individual":
				proteinFrom = encounter.getReference()
				anmReferenceTemp = encounter.accessANMs().getANMReference()
				anmReference = extendModel(anmReferenceTemp[0], anmReferenceTemp[1], proteinFrom, norm=True)

				ensem = sampleModes(anmReference[0][6:106], proteinFrom, 100, RMSDPredicted)
			
			if self.config.investigationsOn == "Complex":

				proteinFrom = encounter.unboundComplexAligned.complex
				anmReferenceTemp = encounter.accessANMs().getANMComplexSlc()
				anmReference = extendModel(anmReferenceTemp[0], anmReferenceTemp[1], proteinFrom, norm=True)
				ensem = sampleModes(anmReference[0][6:106], proteinFrom, 100, RMSDPredicted)


			# RMSDtostart = []
			# RMSDtobound = []
			# print type(proteinFrom)
			# print np.unique(proteinFrom.select('calpha').getIndices())
			ensem.setAtoms(proteinFrom)
   # 			for sample_index in range (0, 100):
   # 				RMSDtostart.append(calcRMSD(ensem.getConformation(sample_index).getCoords(), proteinFrom))
   # 				if self.bound_provided == True:
		# 			RMSDtobound.append(calcRMSD(ensem.getConformation(sample_index).getCoords(), proteinTo))
			# proteinFromTemp = proteinFrom.copy()	
			# proteinFromTemp.addCoordset(ensem)
			encounter.resultsPrinter.setEnsemble(ensem)
			# encounter.resultsPrinter.setRMSDtostart(RMSDtostart)
			# encounter.resultsPrinter.setRMSDtobound(RMSDtobound)

			encounter.resultsPrinter.setRMSDPrediction(RMSDPredicted)
			encounter.resultsPrinter.writeRMSDResults(self.config.outputPath, self.config.experimentNamePrefix, self.utils)
			encounter.resultsPrinter.writeSampleResults(self.config.outputPath, self.config.experimentNamePrefix, self.utils)


			print "Results have been written to: ", resultsPath
			return resultsPath

			
if  __name__ =='__main__':
	# parse the arguments
	parser = argparse.ArgumentParser(description='NMAUnified inherits from TNMABase and performs a NMA for individual and complex proteins.',
	epilog="Provide a configurations class and 4 proteins, Result: NMA on the proteins with the results written in an output folder, parameters specified in the configuration.")
	parser.add_argument('config', help='Configurations class')
	parser.add_argument('protein1_A', help='Protein 1 in first conformational state')
	parser.add_argument('protein2_A', nargs='?', help='Protein 2 in first conformational state')
	parser.add_argument('protein1_B', nargs='?', help='Protein 1 in second conformational state')
	parser.add_argument('protein2_B', nargs='?', help='Protein 2 in second conformational state')
	parser.add_argument('--profiler', action="store_true", help='Run the program with a profiler attached to it, that writes runtime information after a successful run."')
	parser.add_argument('--outputPath', help='Directly specify the output path of the program, ignoring the output path in the Configuration')
	parser.add_argument('--affirmProteinNames', help='If the arguments for protein1_A, protein2_A, protein1_B or protein2_B do not follow the xxxx_r/l_l/b.pdb naming scheme of the Protein-Protein Docking Benchmark 4.0, provide "receptor" to tell the program that protein1 is a receptor (and protein2 therefore a ligand), or provide "ligand" to state the opposite. The titles of the proteins will be adjusted accordingly.')
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = parser.parse_args()    
	# instantiate the NMAUnified class
	nmaCalculation = NMAUnified(args)
	# run the program, with optionally a profiler attached
	if args.profiler:
		profiler = cProfile.Profile()
		resultsPath = profiler.runcall(nmaCalculation.mainCalculation, nmaCalculation.args.protein1_A, nmaCalculation.args.protein2_A, nmaCalculation.args.protein1_B, nmaCalculation.args.protein2_B, nmaCalculation.args.outputPath)
		profiler.dump_stats(resultsPath+'profiler.pro') # for visualization        
		# write human-readable output
		stream = open(resultsPath+'profiler.txt', 'w')
		try:
			p = pstats.Stats(resultsPath+'profiler.pro', stream=stream)
			p.sort_stats('cumulative').print_stats(15)
		finally:
			stream.close()
		# bz2 the profiler binary file
		try:
			os.chdir(resultsPath)
			shellCommand = "tar -jcvf profiler.pro.tar.bz2 profiler.pro"
			os.system(shellCommand)
			shellCommand = "rm -rf profiler.pro"
			os.system(shellCommand)            
			# change the working path back to the directory of the program
			os.chdir(os.path.dirname(os.path.abspath(__file__)))
		except Exception, err:
			print "Exception occurred when gzipping and removing profiler.pro: ", err
			print traceback.format_exc()
		# to visualize, see visualizePstats.sh in helperScrips or run: 
		# $ run gprof2dot -f pstats profiler.pro | dot -Tpng -o gprof2dot_output.png
		# OR with the RunSnakeRun visualizer
		# $ runsnake profiler.pro
	else:
		nmaCalculation.mainCalculation(nmaCalculation.args.protein1_A, nmaCalculation.args.protein2_A, nmaCalculation.args.protein1_B, nmaCalculation.args.protein2_B, nmaCalculation.args.outputPath)


	print "NMA calculations and RMSD prediction finished."
			
