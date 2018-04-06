'''
Created on Jan 17, 2014

@author: oliwa
'''
import sys as sys
import numpy as np
from prody.dynamics.anm import calcANM, ANM
from prody.dynamics.editing import extendModel, sliceModel
from prody.dynamics.functions import saveModel, loadModel, writeArray
from prody.proteins.pdbfile import writePDB, parsePDB
from prody.dynamics.mode import Vector
from prody.measure.measure import calcCenter, calcDistance
from prody.dynamics.compare import calcOverlap, calcCumulOverlap,\
	calcSubspaceOverlap, calcCovOverlap, printOverlapTable, getOverlapTable
from prody.apps.prody_apps.prody_contacts import prody_contacts
import traceback
from prody.dynamics.nmdfile import writeNMD
import scipy as sp

class ANMs(object):
	"""
	This class holds all the ANMs for an encounter.
	"""

	def __init__(self, utils):
		"""
		Constructor
		"""
		self.utils = utils
		
	def createSlcSelectionString(self, reference, isBoundComplex, ref_chain, referenceTitle):
		""" Under the assumption that is reflected in the Benchmark 4.0 that the receptor atoms are set before the
		ligand atoms (spacially in the PDB file), if the current protein under investigation is a ligand, 
		an offset is added to the selection string to match the atoms of the ligand from the complex. """
		if isBoundComplex and not self.utils.isReceptor(referenceTitle):
			print "adding offset"
			return self.utils.addOffset(ref_chain.getSelstr(), reference.select('segment "R."').numAtoms())
		else:
			print "using original selstr"
			return ref_chain.getSelstr()
		
	def calcANMs(self, reference, ref_chain, numberOfModes, encounter, selstr='calpha', whatAtomsToMatch='calpha', modified="", forceRebuild=False, isBoundComplex=False):
		# if the base model does not exist, it needs to be created along with the
		# extended and slicedback models
		if forceRebuild or not self.doesANMExist(reference, numberOfModes, selstr, whatAtomsToMatch, modified):
			# Create the anm
			anm = calcANM(reference, n_modes=numberOfModes, selstr=selstr)
			# First extend the anm on all atoms
			anm_extend = extendModel(anm[0], anm[1], reference, norm=True)
			# Then slice it back to matched
			selectionAtoms = self.createSlcSelectionString(reference, isBoundComplex, ref_chain, encounter.getReference().getTitle())
			anm_slc = sliceModel(anm_extend[0], anm_extend[1], selectionAtoms)
			# If isBoundComplex, slice one anm back to its overall matched chains
			if isBoundComplex:
				selectionAtomsCounterpart = self.createSlcSelectionString(reference, isBoundComplex, encounter.getBoundCounterpartChain(), encounter.getUnboundCounterpart().getTitle())
				anm_slc_counterpart= sliceModel(anm_extend[0], anm_extend[1], selectionAtomsCounterpart)
				
			# Save the models
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
			print "created and saved models"
#             print "reference, it is the complex: ", reference.select('calpha and segment "R."').numAtoms()
#             print "to slice on, it is the mob_chain: ", ref_chain.numAtoms()
			print "anm hessian               : " + str(anm[0].getHessian().shape)
			print "number of calpha          : " + str(reference.select('calpha').numAtoms())
			print "anm size                  : " + str(anm[0].getArray().shape)
			print "anm_ext size              : " + str(anm_extend[0].getArray().shape)
			print "anm_slice size            : " + str(anm_slc[0].getArray().shape)
			print "selectionAtoms            : " + selectionAtoms
			if isBoundComplex:
				print "anm slice counterpart size: " + str(anm_slc_counterpart[0].getArray().shape)
				print "selectionAtoms counterpart: " + selectionAtomsCounterpart
			# Save the models"
			self._anm = anm
			self._anm_extend = anm_extend
			self._anm_slc = anm_slc
			if isBoundComplex:
				self._anm_slc_counterpart = anm_slc_counterpart
		else:
			#raise Exception("Problem with capturing the selection of saved models, do not use load models from files now.")
			try:
				# load models
				anmModel = loadModel(self.getANMPath(reference, numberOfModes, selstr, whatAtomsToMatch)+".anm.npz")
				anm_extendModel = loadModel(self.getANMPath(reference, numberOfModes, selstr, whatAtomsToMatch, modified="extended")+".nma.npz")
				anm_slcModel = loadModel(self.getANMPath(reference, numberOfModes, selstr, whatAtomsToMatch, modified="slicedback")+".nma.npz")
				
				# store models selections
				anmModelSelection = reference.select(selstr)
				anm_extendModelSelection = reference
				selectionAtoms = self.createSlcSelectionString(reference, isBoundComplex, ref_chain)
				anm_slcModelSelection = reference.select(selectionAtoms)
				
				# recombine models and selections as tuples
				anm = (anmModel, anmModelSelection)
				anm_extend = (anm_extendModel, anm_extendModelSelection)
				anm_slc = (anm_slcModel, anm_slcModelSelection)
				
				print "loaded models"
				print "anm size      : " + str(anm[0].getArray().shape)
				print "anm_ext size  : " + str(anm_extend[0].getArray().shape)
				print "anm_slice size: " + str(anm_slc[0].getArray().shape)
				print "selectionAtoms: " + selectionAtoms
				self._anm = anm
				self._anm_extend = anm_extend
				self._anm_slc = anm_slc
			except IOError as e:
				print "Error loading ANM models from disc: "+str(e)
				
	def calcANMsForPart2a2k(self, reference, counterpart, proteinComplex, ref_chain, counterpart_chain, chain_complex, numberOfModes, selstr='calpha', whatAtomsToMatch='calpha'):
			# Create the anm of reference, counterpart and proteinComplex)
#             print "reference, counterpart, proteinComplex, chain_complex (calphas, calphas*3-6)          : ", (reference.select('calpha').numAtoms(), reference.select('calpha').numAtoms()*3 -6), (counterpart.select('calpha').numAtoms(), counterpart.select('calpha').numAtoms()*3-6), (proteinComplex.select('calpha').numAtoms(), proteinComplex.select('calpha').numAtoms()*3-6), (chain_complex.select('calpha').numAtoms(), chain_complex.select('calpha').numAtoms()*3 -6)
#             print "anm_reference, anm_counterpart, anm_complex hessian shapes                     : ", anm_reference[0].getHessian().shape, anm_counterpart[0].getHessian().shape, anm_complex[0].getHessian().shape
#             print "anm_reference, anm_counterpart, anm_complex, anm_complex_slc getArray() shapes : ", anm_reference[0].getArray().shape, anm_counterpart[0].getArray().shape, anm_complex[0].getArray().shape, anm_complex_slc[0].getArray().shape
			self._anm_reference, self._anm_reference_slc = self._calcANMsUnified(reference, ref_chain, numberOfModes/2, selstr, whatAtomsToMatch)
			self._anm_counterpart, self._anm_counterpart_slc = self._calcANMsUnified(counterpart, counterpart_chain, numberOfModes/2, selstr, whatAtomsToMatch)
			
#             print "15 ang contact before moving atoms:", proteinComplex.select('same residue as exwithin 15 of segment "L." ').numAtoms()
#             self._moveSegment(proteinComplex, "L", 30)
#             if proteinComplex.select('same residue as exwithin 15 of segment "L." ') != None:
#                 print "15 ang contact after moving atoms: ", proteinComplex.select('same residue as exwithin 15 of segment "L." ').numAtoms()
#             else:
#                 print "15 ang contact after moving atoms: 0"
			self._anm_complex, self._anm_complex_slc = self._calcANMsUnified(proteinComplex, chain_complex, numberOfModes, selstr, whatAtomsToMatch)
			
			#self.utils.testHessianSubMatrices(self._anm_reference, self._anm_counterpart, self._anm_complex)
			# check blockmatrix differences and pymol output
#             useRelError = True
			#significantDifferences = self.utils.testBlockMatrixMembership(self._anm_reference[0].getHessian(), self._anm_counterpart[0].getHessian(), self._anm_complex[0].getHessian(), useRelativeError=useRelError)
			#self.utils.whichPatternsAreAffectedbySignificantDifferences(significantDifferences)
#             assert reference.getResnums()[0] == proteinComplex.getResnums()[0]
			#print self.utils.significantDifferencesToPymolResiduesString(significantDifferences, reference.getResnums()[0])
			print "anm_reference_slc, anm_counterpart_slc, anm_complex_slc getArray() shapes : ", self._anm_reference_slc[0].getArray().shape, self._anm_counterpart_slc[0].getArray().shape, self._anm_complex_slc[0].getArray().shape

	def calcANMsUnified(self, reference, counterpart, proteinComplex, numberOfModes, encounter, ref_chain = None, counterpart_chain = None, chain_complex = None, selstr='calpha', whatAtomsToMatch='calpha',):
		""" Calculate the ANMs for the NMA. If examinations on the complex, it is assumed (for now) that the reference protein is the receptor. """
		if (ref_chain == None) and (counterpart_chain == None) and (chain_complex == None):
			self.bound_provided = False
		else:
			self.bound_provided = True
		if self.utils.config.investigationsOn == "Individual":
			assert self.utils.config.whichCustomHIndividual == "HC_subvector" or self.utils.config.whichCustomHIndividual == "submatrix" or self.utils.config.whichCustomHIndividual == "canonical"
			numberOfModesComplex = min((proteinComplex.select('calpha').numAtoms()*3 - 6), self.utils.config.maxModesToCalculate)
			if ref_chain != None:
				self._anm_reference, self._anm_reference_slc = self._calcANMsUnified(reference, numberOfModes, selstr, whatAtomsToMatch, ref_chain)
			else:
				self._anm_reference, self._anm_reference_slc = self._calcANMsUnified(reference, numberOfModes, selstr, whatAtomsToMatch)

			self._anm_counterpart = calcANM(counterpart, n_modes = numberOfModes, selstr = selstr, zeros = True)
			if chain_complex != None:
				self._anm_complex, self._anm_complex_slc = self._calcANMsUnified(proteinComplex, numberOfModesComplex, selstr, whatAtomsToMatch, chain_complex)
			else:
				self._anm_complex, self._anm_complex_slc = self._calcANMsUnified(proteinComplex, numberOfModesComplex, selstr, whatAtomsToMatch)
		elif self.utils.config.investigationsOn == "Complex":
			numberOfModesComplex = numberOfModes*2
			self._anm_reference, self._anm_reference_slc = self._calcANMsUnified(reference, numberOfModes, selstr, whatAtomsToMatch, ref_chain)
			self._anm_counterpart, self._anm_counterpart_slc = self._calcANMsUnified(counterpart, numberOfModes, selstr, whatAtomsToMatch, counterpart_chain)
			self._anm_complex, self._anm_complex_slc = self._calcANMsUnified(proteinComplex, numberOfModesComplex, selstr, whatAtomsToMatch, chain_complex)

		print "anm_reference anm_counterpart, anm_complex getArray() shapes : ", self._anm_reference[0].getArray().shape, self._anm_counterpart[0].getArray().shape, self._anm_complex[0].getArray().shape
		print "anm_reference_slc, anm_complex_slc getArray() shapes : ", self._anm_reference_slc[0].getArray().shape, self._anm_complex_slc[0].getArray().shape
		
		# create custom H via U1
		if self.utils.config.customH:
			HC = self._anm_complex[0].getHessian()
			if self.utils.isReceptor(reference.getTitle()):
				HR = self._anm_reference[0].getHessian()
				HL = self._anm_counterpart[0].getHessian()
			else:
				HR = self._anm_counterpart[0].getHessian()
				HL = self._anm_reference[0].getHessian()
			HRtilde = HC[:HR.shape[0], :HR.shape[1]]
			HLtilde = HC[HR.shape[0]:HR.shape[0]+HL.shape[0], HR.shape[1]:HR.shape[1]+HL.shape[1]]    
			assert HR.shape == HRtilde.shape
			assert HL.shape == HLtilde.shape
		
			# for now assert that reference is always the receptor
			if self.utils.config.investigationsOn == "Complex":
				assert self.utils.isReceptor(reference.getTitle())
			HCcustomBuild = np.zeros((HC.shape[0], HC.shape[1]))
		
			if self.utils.isReceptor(reference.getTitle()):
				if self.utils.config.whichCustomHC == "HC_U1" or self.utils.config.whichCustomHC == "HC_U1_1k1k":
					HRtildeH_ANew, interCalphaIndicesHR = self.calcCustomH_ANew(HR.copy(), encounter.getReference(), encounter.getUnboundCounterpart(), encounter, "C_u", "r_ij", True, selstr)
					HLtildeH_ANew, interCalphaIndicesHL = self.calcCustomH_ANew(HL.copy(), encounter.getUnboundCounterpart(), encounter.getReference(), encounter, "C_u", "r_ij", False, selstr)
					HRL_new = self.calcCustomH_ANew_IJ(encounter.getReference(), encounter.getUnboundCounterpart(), encounter, False, "r_ij", True, selstr)
				elif self.utils.config.whichCustomHC == "HC_0" or self.utils.config.whichCustomHC == "HC_06":
					HRtildeH_ANew = HR.copy()
					HLtildeH_ANew = HL.copy()
					HRL_new = np.zeros(((reference.select('calpha').numAtoms()*3), (counterpart.select('calpha').numAtoms()*3) ))
					interCalphaIndicesHR = None
					interCalphaIndicesHL = None                    
				print "reference is receptor, shapes of HRtilde, HLtilde, HRL: ", HRtildeH_ANew.shape, HLtildeH_ANew.shape, HRL_new.shape
			else:
				if self.utils.config.whichCustomHC == "HC_U1":
					HRtildeH_ANew, interCalphaIndicesHR = self.calcCustomH_ANew(HR.copy(), encounter.getUnboundCounterpart(), encounter.getReference(), encounter, "C_u", "r_ij", False, selstr)
					HLtildeH_ANew, interCalphaIndicesHL = self.calcCustomH_ANew(HL.copy(), encounter.getReference(), encounter.getUnboundCounterpart(), encounter, "C_u", "r_ij", True, selstr)
					HRL_new = self.calcCustomH_ANew_IJ(encounter.getUnboundCounterpart(), encounter.getReference(), encounter, False, "r_ij", False, selstr)                                               
				print "reference is ligand, shapes of HLtilde, HRtilde, HRL: ", HLtildeH_ANew.shape, HRtildeH_ANew.shape, HRL_new.shape                 
			#  put the new HRtilde and HLtilde inside HC
			HCcustomBuild[:HR.shape[0], :HR.shape[1]] = HRtildeH_ANew
			HCcustomBuild[HR.shape[0]:HR.shape[0]+HL.shape[0], HR.shape[1]:HR.shape[1]+HL.shape[1]] = HLtildeH_ANew
			HCcustomBuild[0:HR.shape[0], HR.shape[1]:HR.shape[1]+HL.shape[1]] = HRL_new
			HCcustomBuild[HR.shape[0]:HR.shape[0]+HL.shape[0], 0:HR.shape[1]] = HRL_new.T
			
			# optional assertion to test if HCcustomBuild equals the original HC if k = 1 and d = 15 (default ProDy settings)
			if (self.utils.config.whichCustomHC == "HC_U1" and self.utils.config.customHRdistance == 15 and self.utils.config.customForceConstant == 1.0):
				# assert np.allclose(HC, HCcustomBuild) # assert this if k = 1, A = 15
				print "not asserting HCcustomBuild equals original HC with k1 A15"
			# Projection
			#     def projectHessian(self, hessian, reference, proteinComplex, referenceSegment, projectionStyle, projectOnlyReferencePartOfHC=False, interCalphaIndices=None):
			if self.utils.config.projectHessian:
				if self.utils.config.investigationsOn == "Individual":
					if self.utils.isReceptor(reference.getTitle()):
						if self.utils.config.whichCustomHC == "HC_U1":
							if self.utils.config.projectionStyle == "full" or self.utils.config.projectionStyle == "intra":
								if self.utils.config.whichCustomHIndividual == "HC_subvector":
									HCcustomBuild = self.projectHessian(HCcustomBuild.copy(), reference, proteinComplex, "R", self.utils.config.projectionStyle, True, interCalphaIndicesHR)
									#HCcustomBuild = self.projectHessian(HCcustomBuild.copy(), proteinComplex, proteinComplex, '', self.utils.config.projectionStyle, False, interCalphaIndicesHR) 
								elif self.utils.config.whichCustomHIndividual == "submatrix":
									HRtildeH_ANew = self.projectHessian(HRtildeH_ANew.copy(), reference, proteinComplex, "R", self.utils.config.projectionStyle, False, interCalphaIndicesHR)
							elif self.utils.config.projectionStyle == "fixedDomainFrame":
								HCcustomBuild = self.transformHessianToFixedDomainFrame(HCcustomBuild.copy(), reference, proteinComplex, "R", self.utils.config.projectionStyle)
					# else reference is the ligand
					else:
						if self.utils.config.whichCustomHC == "HC_U1":
							if self.utils.config.projectionStyle == "full" or self.utils.config.projectionStyle == "intra":
								if self.utils.config.whichCustomHIndividual == "HC_subvector":
									HCcustomBuild = self.projectHessian(HCcustomBuild.copy(), reference, proteinComplex, "L", self.utils.config.projectionStyle, True, interCalphaIndicesHL)
									#HCcustomBuild = self.projectHessian(HCcustomBuild.copy(), proteinComplex, proteinComplex, '', self.utils.config.projectionStyle, False, interCalphaIndicesHL) 
								elif self.utils.config.whichCustomHIndividual == "submatrix":
									HLtildeH_ANew = self.projectHessian(HLtildeH_ANew.copy(), reference, proteinComplex, "L", self.utils.config.projectionStyle, False, interCalphaIndicesHL)
							elif self.utils.config.projectionStyle == "fixedDomainFrame":
								HCcustomBuild = self.transformHessianToFixedDomainFrame(HCcustomBuild.copy(), reference, proteinComplex, "L", self.utils.config.projectionStyle)
				elif self.utils.config.investigationsOn == "Complex":
					# project out the rigid body motions of the receptor. if the goal is to project the whole complex, do: HCcustomBuild = self.projectHessian(HCcustomBuild, proteinComplex, proteinComplex, '') 
					if self.utils.config.projectionStyle == "full" or self.utils.config.projectionStyle == "intra":
						HCcustomBuild = self.projectHessian(HCcustomBuild.copy(), reference, proteinComplex, "R", self.utils.config.projectionStyle, True, interCalphaIndicesHR)     
					elif self.utils.config.projectionStyle == "fullComplex":
						HCcustomBuild = self.projectHessian(HCcustomBuild.copy(), proteinComplex, proteinComplex, '', self.utils.config.projectionStyle)
					elif self.utils.config.projectionStyle == "fixedDomainFrame":
						HCcustomBuild = self.transformHessianToFixedDomainFrame(HCcustomBuild.copy(), reference, proteinComplex, "R", self.utils.config.projectionStyle)
					else:
						raise Exception('unknown projection style')
		
			if self.utils.config.investigationsOn == "Complex" or self.utils.config.whichCustomHIndividual == "HC_subvector":
				# Create the custom complex ANM 
				self._anm_complex_tilde = ANM(self._anm_complex[0].getTitle()+"_"+self.utils.config.whichCustomHC)
				self._anm_complex_tilde.setHessian(HCcustomBuild)
				if self.utils.config.calculateZeroEigvalModes:
					if self.utils.config.whichCustomHC == "HC_0" or self.utils.config.whichCustomHC == "HC_06":
						numberOfModesComplex += 6
					self._anm_complex_tilde.calcModes(n_modes=numberOfModesComplex, zeros=True)
				else:
					self._anm_complex_tilde.calcModes(n_modes=numberOfModesComplex)
								  
				# Extend the self._anm_reference_tilde on all atoms
				anm_complex_tilde_extend = extendModel(self._anm_complex_tilde, self._anm_complex[1], proteinComplex, norm=True)        
				# Then slice the anm_complex to the matched atoms
				self._anm_complex_tilde_slc = sliceModel(anm_complex_tilde_extend[0], anm_complex_tilde_extend[1], selstr)
				# Normalize the modes of the sliced anm
				self._anm_complex_tilde_slc = self.getNormalizedANM(self._anm_complex_tilde_slc)
				# Replace the complex anm and the complex_slc anm with the modified ANMs
				print "Replacing ANM H with ANM Htilde for the complex"
				self._anm_complex = (self._anm_complex_tilde, self._anm_complex[1])
				self._anm_complex_slc = self._anm_complex_tilde_slc
		
			# modify HR to have the sliced part of HC_tilde
			if self.utils.config.investigationsOn == "Individual":
				if self.utils.config.whichCustomHIndividual == "HC_subvector":
					Marray = self.utils.sliceComplexModestoMatchProtein(self._anm_complex[0].getArray(), reference, encounter.getReferenceSegment())
					self._anm_reference_tilde = ANM(self._anm_reference[0].getTitle()+"_"+self.utils.config.whichCustomHC)
					self._anm_reference_tilde.setEigens(Marray, self._anm_complex[0].getEigvals())
					self._anm_reference_tilde = (self._anm_reference_tilde, self._anm_reference[1])
					self._anm_reference_tilde = self.getNormalizedANM(self._anm_reference_tilde)
				# submatrix, take the new HRtilde/HLtilde, re-calculate its modes and replace the previous ANM 
				elif self.utils.config.whichCustomHIndividual == "submatrix":
					if self.utils.isReceptor(reference.getTitle()):
						submatrix = HRtildeH_ANew
					else:
						submatrix = HLtildeH_ANew
					self._anm_reference_tilde = ANM(self._anm_reference[0].getTitle()+"_"+self.utils.config.whichCustomHC)
					self._anm_reference_tilde.setHessian(submatrix)
					if self.utils.config.calculateZeroEigvalModes:
						self._anm_reference_tilde.calcModes(n_modes=numberOfModes, zeros=True)
					else:
						self._anm_reference_tilde.calcModes(n_modes=numberOfModes)
					self._anm_reference_tilde = (self._anm_reference_tilde, self._anm_reference[1])   
											   
				# Extend the self._anm_reference_tilde on all atoms
				anm_reference_tilde_extend = extendModel(self._anm_reference_tilde[0], self._anm_reference[1], reference, norm=True)        
				# Then slice the anm_reference to the matched
				self._anm_reference_tilde_slc = sliceModel(anm_reference_tilde_extend[0], anm_reference_tilde_extend[1], selstr)
				self._anm_reference_tilde_slc = self.getNormalizedANM(self._anm_reference_tilde_slc)
										 
				# Replace reference and reference_slc with the modified ANMs
				print "Replacing ANM H with ANM Htilde for the reference"
				self._anm_reference = self._anm_reference_tilde
				self._anm_reference_slc = self._anm_reference_tilde_slc

	def calcANMsForPart2b2k(self, reference, counterpart, proteinComplex, ref_chain, counterpart_chain, chain_complex, numberOfModes, encounter, selstr='calpha', whatAtomsToMatch='calpha'):
			""" Unbound complex to bound complex NMA, it is assumed that the reference is the receptor and is the first object in the complex pdb file 
				This method creates self.* NMA objects
			
				Args: 
					reference: the receptor protein
					counterpart: the ligand protein
					proteinComplex: the protein complex
					ref_chain: the matched part of the reference
					counterpart_chain: the matched part of the counterpart
					chain_complex: the matched part on the complex
					numberOfModes: the 2k number of modes
					encounter: object aggregating proteins
					selstr: the selection string for the NMA, course grained is calpha
			"""
			# Create the anm of reference, counterpart and proteinComplex)
			self._anm_reference, self._anm_reference_slc = self._calcANMsUnified(reference, ref_chain, numberOfModes/2, selstr, whatAtomsToMatch)
			self._anm_counterpart, self._anm_counterpart_slc = self._calcANMsUnified(counterpart, counterpart_chain, numberOfModes/2, selstr, whatAtomsToMatch)
			self._anm_complex, self._anm_complex_slc = self._calcANMsUnified(proteinComplex, chain_complex, numberOfModes, selstr, whatAtomsToMatch)
			
			print "anm_reference anm_counterpart, anm_complex getArray() shapes : ", self._anm_reference[0].getArray().shape, self._anm_counterpart[0].getArray().shape, self._anm_complex[0].getArray().shape
			print "anm_reference_slc, anm_counterpart_slc, anm_complex_slc getArray() shapes : ", self._anm_reference_slc[0].getArray().shape, self._anm_counterpart_slc[0].getArray().shape, self._anm_complex_slc[0].getArray().shape

			
			# modify the hessians
			if self.utils.config.customH:
				HC = self._anm_complex[0].getHessian()
				if self.utils.isReceptor(reference.getTitle()):
					HR = self._anm_reference[0].getHessian()
					HL = self._anm_counterpart[0].getHessian()
				else:
					HR = self._anm_counterpart[0].getHessian()
					HL = self._anm_reference[0].getHessian()
				HRtilde = HC[:HR.shape[0], :HR.shape[1]]
				HLtilde = HC[HR.shape[0]:HR.shape[0]+HL.shape[0], HR.shape[1]:HR.shape[1]+HL.shape[1]]    
				assert HR.shape == HRtilde.shape
				assert HL.shape == HLtilde.shape
 
				# for now assert that reference is always the receptor, in case of complex investigation
				assert self.utils.isReceptor(reference.getTitle())
				HCcustomBuild = np.zeros((HC.shape[0], HC.shape[1]))

				if self.utils.config.whichCustomHC == "HC_U1":
					# create the complex hessian with interactions on the off diagonal using U1
					print "HC_U1"
					HRtildeH_ANew = self.calcCustomH_ANew(HR.copy(), encounter.getReference(), encounter.getUnboundCounterpart(), encounter, "C_u", "r_ij", True, selstr)
					HLtildeH_ANew = self.calcCustomH_ANew(HL.copy(), encounter.getUnboundCounterpart(), encounter.getReference(), encounter, "C_u", "r_ij", False, selstr)
					HRL_new = self.calcCustomH_ANew_IJ(encounter.getReference(), encounter.getUnboundCounterpart(), encounter, False, "r_ij", True, selstr)
				elif self.utils.config.whichCustomHC == "HC_0" or self.utils.config.whichCustomHC == "HC_06":
					# create the hessian by just using canonical HR and HL and offmatrices zero
					print "HC_0 or HC_06"
					HRtildeH_ANew = HR.copy()
					HLtildeH_ANew = HL.copy()
					HRL_new = np.zeros(((reference.select('calpha').numAtoms()*3), (counterpart.select('calpha').numAtoms()*3) ))
				print "reference is receptor, shapes of HRtilde, HLtilde, HRL: ", HRtildeH_ANew.shape, HLtildeH_ANew.shape, HRL_new.shape            
				print "finished projecting H, anm_reference_tilde calc modes"
				#  put the new HRtilde and HLtilde inside HC
				HCcustomBuild[:HR.shape[0], :HR.shape[1]] = HRtildeH_ANew
				HCcustomBuild[HR.shape[0]:HR.shape[0]+HL.shape[0], HR.shape[1]:HR.shape[1]+HL.shape[1]] = HLtildeH_ANew
				HCcustomBuild[0:HR.shape[0], HR.shape[1]:HR.shape[1]+HL.shape[1]] = HRL_new
				HCcustomBuild[HR.shape[0]:HR.shape[0]+HL.shape[0], 0:HR.shape[1]] = HRL_new.T
				
				#if self.utils.config.whichCustomHC == "HC_U1":
				#    assert np.allclose(HC, HCcustomBuild) # assert this if k = 1, A = 15
				#    print "asserted HC with k1 A 15"
				if self.utils.config.projectHessian:
					HCcustomBuild = self.projectHessian(HCcustomBuild, proteinComplex, proteinComplex, '') 
				#  make HC anm
				self._anm_complex_tilde = ANM(self._anm_complex[0].getTitle()+"_"+self.utils.config.whichCustomHC)
				self._anm_complex_tilde.setHessian(HCcustomBuild)
				self._anm_complex_tilde.calcModes(n_modes=numberOfModes)
				 
				# Extend the self._anm_reference_tilde on all atoms
				anm_complex_tilde_extend = extendModel(self._anm_complex_tilde, self._anm_complex[1], proteinComplex, norm=True)        
				# Then slice the anm_complex to the matched atoms
				self._anm_complex_tilde_slc = sliceModel(anm_complex_tilde_extend[0], anm_complex_tilde_extend[1], chain_complex.getSelstr())
				# Replace the complex anm and the complex_slc anm with the modified ANMs
				print "Replacing ANM H with ANM Htilde for the complex"
				self._anm_complex = (self._anm_complex_tilde, self._anm_complex[1])
				self._anm_complex_slc = self._anm_complex_tilde_slc

	def calcANMsForPart2b(self, reference, counterpart, proteinComplex, ref_chain, counterpart_chain, chain_complex, numberOfModes, encounter, selstr='calpha', whatAtomsToMatch='calpha'):
			""" Create the ANMs of the reference, counterpart and complex objects. If set in config, project the hessian matrix of the reference 
			to ensure 6 zero eigenvalue modes, see formula 8.27 from the book "A practical introduction to the simulation of molecular dynamics", Field. """
			
			self._anm_reference, self._anm_reference_slc = self._calcANMsUnified(reference, ref_chain, numberOfModes, selstr, whatAtomsToMatch)
			self._anm_counterpart = calcANM(counterpart, selstr=selstr)
#             self._moveSegment(proteinComplex, "L", 50)
			numberOfModesComplex = min((proteinComplex.select('calpha').numAtoms()*3 - 6), self.utils.config.maxModesToCalculate)
			self._anm_complex, self._anm_complex_slc = self._calcANMsUnified(proteinComplex, chain_complex, numberOfModesComplex, selstr, whatAtomsToMatch)
			
			# project hessian matrix
			if self.utils.config.projectHessian:
				HC = self._anm_complex[0].getHessian()
				if self.utils.isReceptor(reference.getTitle()):
					HR = self._anm_reference[0].getHessian()
					HL = self._anm_counterpart[0].getHessian()
				else:
					HR = self._anm_counterpart[0].getHessian()
					HL = self._anm_reference[0].getHessian()
				HRtilde = HC[:HR.shape[0], :HR.shape[1]]
				HLtilde = HC[HR.shape[0]:HR.shape[0]+HL.shape[0], HR.shape[1]:HR.shape[1]+HL.shape[1]]    
				assert HR.shape == HRtilde.shape
				assert HL.shape == HLtilde.shape

				##
				#writeArray("HRtildefromHC.txt", HRtilde, format='%f')
				#writeArray("HLtildefromHC.txt", HLtilde, format='%f')
				##
				# Create the tilde ANM
				self._anm_reference_tilde = ANM(self._anm_reference[0].getTitle()+"_tilde")
				# Here the PH'P treatment for the hessian matrix from the normal modes book by Field
				if self.utils.isReceptor(reference.getTitle()):
					if self.utils.config.modifyHDelta:
						print "modifying HR with deltaHR"
						HRtilde = self.addscaledHdelta(HR, HRtilde, self.utils.config.deltamultiplicatorForH)
				# if using terms with true bound structure second derivation parts r_{ij}-r_{ij}^{2}
					if self.utils.config.customHR_A:
						#writeArray("originalHR.txt", self._anm_reference[0].getHessian(), format='%f')
						HRtilde = self.calcCustomH_A_NeighborsBound(self._anm_reference[0].getHessian(), encounter, selstr)
						#writeArray("customHRtilde.txt", HRtilde, format='%f')
					print "reference is receptor, shape of HRtilde: ", HRtilde.shape
					HRtilde = self.projectHessian(HRtilde, reference, proteinComplex, encounter.getReferenceSegment()) 
					self._anm_reference_tilde.setHessian(HRtilde)
				else:
					if self.utils.config.modifyHDelta:
						print "modifying HL with deltaHL"
						HLtilde = self.addscaledHdelta(HL, HLtilde, self.utils.config.deltamultiplicatorForH)
				# if using terms with true bound structure second derivation parts r_{ij}-r_{ij}^{2}
					if self.utils.config.customHR_A:
						#writeArray("originalHL.txt", self._anm_reference[0].getHessian(), format='%f')
						HLtilde = self.calcCustomH_A_NeighborsBound(self._anm_reference[0].getHessian(), encounter, selstr)
						#writeArray("customHLtilde.txt", HLtilde, format='%f')
					print "reference is ligand, shape of HLtilde: ", HLtilde.shape
					HLtilde = self.projectHessian(HLtilde, reference, proteinComplex, encounter.getReferenceSegment())
					self._anm_reference_tilde.setHessian(HLtilde)
				print "finished projecting H, anm_reference_tilde calc modes"
				#  testing of projected eigenvals
				self._anm_reference_tilde.calcModes(n_modes=numberOfModes)
				#print "HR      eigenvals: ", self._anm_reference[0].getEigvals()[0:10]
				#print "HRtilde eigenvals: ", self._anm_reference_tilde.getEigvals()[0:10]

				# Extend the self._anm_reference_tilde on all atoms
				anm_reference_tilde_extend = extendModel(self._anm_reference_tilde, self._anm_reference[1], reference, norm=True)        
				# Then slice the anm_reference to the matched
				self._anm_reference_tilde_slc = sliceModel(anm_reference_tilde_extend[0], anm_reference_tilde_extend[1], ref_chain.getSelstr())
				# Replace reference and reference_slc with the modified ANMs
				print "Replacing ANM H with ANM Htilde for the reference"
				self._anm_reference = (self._anm_reference_tilde, self._anm_reference[1])
				self._anm_reference_slc = self._anm_reference_tilde_slc
				if self.utils.config.HR1kHRtilde1k:
					self._anm_reference_original, self._anm_reference_slc_original = self._calcANMsUnified(reference, ref_chain, numberOfModes, selstr, whatAtomsToMatch)

	def calcANMsForPart2bIndividualProtein_U1(self, reference, counterpart, proteinComplex, ref_chain, counterpart_chain, chain_complex, numberOfModes, encounter, selstr='calpha', whatAtomsToMatch='calpha'):
			""" Create the ANMs of the reference, counterpart and complex objects. If set in config, project the hessian matrix of the reference 
			to ensure 6 zero eigenvalue modes, see formula 8.27 from the book "A practical introduction to the simulation of molecular dynamics", Field. """
			
			self._anm_reference, self._anm_reference_slc = self._calcANMsUnified(reference, ref_chain, numberOfModes, selstr, whatAtomsToMatch)
			self._anm_counterpart = calcANM(counterpart, selstr=selstr)
#             self._moveSegment(proteinComplex, "L", 50)
			numberOfModesComplex = min((proteinComplex.select('calpha').numAtoms()*3 - 6), self.utils.config.maxModesToCalculate)
			self._anm_complex, self._anm_complex_slc = self._calcANMsUnified(proteinComplex, chain_complex, numberOfModesComplex, selstr, whatAtomsToMatch)
			###
			print "anm_reference anm_counterpart, anm_complex getArray() shapes : ", self._anm_reference[0].getArray().shape, self._anm_counterpart[0].getArray().shape, self._anm_complex[0].getArray().shape
			print "anm_reference_slc, anm_complex_slc getArray() shapes : ", self._anm_reference_slc[0].getArray().shape, self._anm_complex_slc[0].getArray().shape
			# create custom H via U1
			if self.utils.config.customH:
				HC = self._anm_complex[0].getHessian()
				if self.utils.isReceptor(reference.getTitle()):
					HR = self._anm_reference[0].getHessian()
					HL = self._anm_counterpart[0].getHessian()
				else:
					HR = self._anm_counterpart[0].getHessian()
					HL = self._anm_reference[0].getHessian()
				HRtilde = HC[:HR.shape[0], :HR.shape[1]]
				HLtilde = HC[HR.shape[0]:HR.shape[0]+HL.shape[0], HR.shape[1]:HR.shape[1]+HL.shape[1]]    
				assert HR.shape == HRtilde.shape
				assert HL.shape == HLtilde.shape
 
				# for now assert that reference is always the receptor
				HCcustomBuild = np.zeros((HC.shape[0], HC.shape[1]))

				if self.utils.isReceptor(reference.getTitle()):
					if self.utils.config.customHR_A:
						#HR, referenceStructure, neighborStructure, encounter, neighborhoodFrom, equilibriumAt, workOnReceptor=True, selstr='calpha'
						HRtildeH_ANew = self.calcCustomH_ANew(HR.copy(), encounter.getReference(), encounter.getUnboundCounterpart(), encounter, "C_u", "r_ij", True, selstr)
						HLtildeH_ANew = self.calcCustomH_ANew(HL.copy(), encounter.getUnboundCounterpart(), encounter.getReference(), encounter, "C_u", "r_ij", False, selstr)
						HRL_new = self.calcCustomH_ANew_IJ(encounter.getReference(), encounter.getUnboundCounterpart(), encounter, False, "r_ij", True, selstr)
					print "reference is receptor, shapes of HRtilde, HLtilde, HRL: ", HRtildeH_ANew.shape, HLtildeH_ANew.shape, HRL_new.shape
				else:
					if self.utils.config.customHR_A:
						HRtildeH_ANew = self.calcCustomH_ANew(HR.copy(), encounter.getUnboundCounterpart(), encounter.getReference(), encounter, "C_u", "r_ij", False, selstr)
						HLtildeH_ANew = self.calcCustomH_ANew(HL.copy(), encounter.getReference(), encounter.getUnboundCounterpart(), encounter, "C_u", "r_ij", True, selstr)
						HRL_new = self.calcCustomH_ANew_IJ(encounter.getUnboundCounterpart(), encounter.getReference(), encounter, False, "r_ij", False, selstr)                                               
					print "reference is ligand, shapes of HLtilde, HRtilde, HRL: ", HLtildeH_ANew.shape, HRtildeH_ANew.shape, HRL_new.shape                 
				print "finished projecting H, anm_reference_tilde calc modes"
				#  put the new HRtilde and HLtilde inside HC
				HCcustomBuild[:HR.shape[0], :HR.shape[1]] = HRtildeH_ANew
				HCcustomBuild[HR.shape[0]:HR.shape[0]+HL.shape[0], HR.shape[1]:HR.shape[1]+HL.shape[1]] = HLtildeH_ANew
				HCcustomBuild[0:HR.shape[0], HR.shape[1]:HR.shape[1]+HL.shape[1]] = HRL_new
				HCcustomBuild[HR.shape[0]:HR.shape[0]+HL.shape[0], 0:HR.shape[1]] = HRL_new.T
				
				#assert np.allclose(HC, HCcustomBuild)
				#sys.exit()
				
				# Project the reference part in the HCcustomBuild matrix
				if self.utils.isReceptor(reference.getTitle()):
					if self.utils.config.customHR_A:
						HCcustomBuildprojected = self.projectHessian(HCcustomBuild.copy(), reference, proteinComplex, "R", True)
				else:
					if self.utils.config.customHR_A:
						HCcustomBuildprojected = self.projectHessian(HCcustomBuild.copy(), reference, proteinComplex, "L", True)

				# Create the custom complex ANM 
				self._anm_complex_tilde = ANM(self._anm_complex[0].getTitle()+"_tilde")
				self._anm_complex_tilde.setHessian(HCcustomBuildprojected)
				if self.utils.config.enforceAllModesAfterProjection:
					self._anm_complex_tilde.calcModes(n_modes=numberOfModes, zeros=True)
				else:
					self._anm_complex_tilde.calcModes(n_modes=numberOfModes)

				# Extend the self._anm_reference_tilde on all atoms
				anm_complex_tilde_extend = extendModel(self._anm_complex_tilde, self._anm_complex[1], proteinComplex, norm=True)        
				# Then slice the anm_complex to the matched atoms
				self._anm_complex_tilde_slc = sliceModel(anm_complex_tilde_extend[0], anm_complex_tilde_extend[1], chain_complex.getSelstr())
				# Replace the complex anm and the complex_slc anm with the modified ANMs
				print "Replacing ANM H with ANM Htilde for the complex"
				self._anm_complex = (self._anm_complex_tilde, self._anm_complex[1])
				self._anm_complex_slc = self._anm_complex_tilde_slc

				# Create custom anm for reference
				if self.utils.config.enforceAllModesAfterProjection:
					Marray = self.utils.sliceComplexModestoMatchProtein(self._anm_complex[0].getArray()[:,6:], reference, encounter.getReferenceSegment())
					self._anm_reference_tilde = ANM(self._anm_reference[0].getTitle()+"_tilde")
					self._anm_reference_tilde.setEigens(Marray, self._anm_complex[0].getEigvals()[6:])
				else:
					Marray = self.utils.sliceComplexModestoMatchProtein(self._anm_complex[0].getArray(), reference, encounter.getReferenceSegment())
					self._anm_reference_tilde = ANM(self._anm_reference[0].getTitle()+"_tilde")
					self._anm_reference_tilde.setEigens(Marray, self._anm_complex[0].getEigvals())

				# Extend the self._anm_reference_tilde on all atoms
				anm_reference_tilde_extend = extendModel(self._anm_reference_tilde, self._anm_reference[1], reference, norm=True)        
				# Then slice the anm_reference to the matched
				self._anm_reference_tilde_slc = sliceModel(anm_reference_tilde_extend[0], anm_reference_tilde_extend[1], ref_chain.getSelstr())
				
				#
				# try modes comparison
#                 ranges = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70]
#                 
#                 try:
#                     subspaceOverlaps = []
#                     for val in ranges:
#                         subspaceOverlaps.append(calcSubspaceOverlap(self._anm_reference[0][0:val], self._anm_reference_tilde[0:val]))
#                     encounter.storeSubSpaceOverlaps(subspaceOverlaps, ranges)
#                 except Exception:
#                     sys.exc_clear()
#                 
#                 try:
#                     MarrayNormed = self.utils.normalized(Marray.copy(), axis=0)
#                     anm_reference_tilde_normed = ANM(self._anm_reference[0].getTitle()+"_tildenormed")
#                     anm_reference_tilde_normed.setEigens(MarrayNormed, self._anm_complex[0].getEigvals())                
#                     covarianceOverlaps = []
#                     for val in ranges:
#                         covarianceOverlaps.append(calcCovOverlap(self._anm_reference[0][0:val], anm_reference_tilde_normed[0:val]))
#                     encounter.storeCovarianceOverlap(covarianceOverlaps, ranges)                
#                 except Exception, err:
#                     #sys.exc_clear()
#                     print "Exception covarianceoverlap occurred: ", err
#                     print traceback.format_exc()
#                 
#                 try:
#                     overlapTable = getOverlapTable(self._anm_reference[0], self._anm_reference_tilde)
#                     encounter.storeOverlapTable(overlapTable)
#                 except Exception:
#                     sys.exc_clear()
				#
				# Replace reference and reference_slc with the modified ANMs
				print "Replacing ANM H with ANM Htilde for the reference"
				self._anm_reference = (self._anm_reference_tilde, self._anm_reference[1])
				self._anm_reference_slc = self._anm_reference_tilde_slc

	def _calcANMsUnified(self, reference, numberOfModes, selstr='calpha', whatAtomsToMatch='calpha', direct_call = None, ref_chain = None):
			# Create the anm of the reference
			#writePDB(reference.getTitle()+"forANMmoved.pdb", reference)
			if self.utils.config.calculateZeroEigvalModes == True:
				anm_reference = calcANM(reference, n_modes=numberOfModes, selstr=selstr, zeros=True)
			else:
				anm_reference = calcANM(reference, n_modes=numberOfModes, selstr=selstr)
			# Extend the anm_reference on all atoms
			anm_reference_extend = extendModel(anm_reference[0], anm_reference[1], reference, norm=True)
			# Then slice the anm_reference to the matched
			if direct_call == None:
				if self.bound_provided == True:
					anm_reference_slc = sliceModel(anm_reference_extend[0], anm_reference_extend[1], ref_chain.getSelstr())
				else:
					anm_reference_slc = sliceModel(anm_reference_extend[0], anm_reference_extend[1], selstr)
			else:
				anm_reference_slc = sliceModel(anm_reference_extend[0], anm_reference_extend[1], selstr)
			# Normalize the slices anm
			anm_reference_slc = self.getNormalizedANM(anm_reference_slc)
			if direct_call == True:
				self._anm_reference = anm_reference
				self._anm_reference_slc = anm_reference_slc
			else:
				return anm_reference, anm_reference_slc
		
	def getNormalizedANM(self, anm):
		""" Normalize the modes of the anm and return this anm object 
		
			Args:
				anm: the anm with modes calculated
				
			Returns: anm with normalized modes
		"""
		M = self.normalizeM(anm[0].getArray())
		eigenvals = anm[0].getEigvals()
		anm[0].setEigens(M, eigenvals)
		return anm
	
	def _moveSegment(self, reference, segment, angstrom):
		""" Move all atoms x,y,z, belonging to the segment the number in angstrom """
		print "15 ang contact before moving atoms:", reference.select('same residue as exwithin 15 of segment "L." ').numAtoms()        
		ref_select = reference.select('segment \"'+segment+'.\"')
		ref_select.setCoords(ref_select.getCoords()+angstrom)
		if reference.select('same residue as exwithin 15 of segment "L." ') != None:
			print "15 ang contact after moving atoms: ", reference.select('same residue as exwithin 15 of segment "L." ').numAtoms()
		else:
			print "15 ang contact after moving atoms: 0"
			
	def replaceReferenceANMs(self, anm_new, reference, ref_chain = None):
		""" Replace the anm of reference with anm_new and normalize along the way. 
		
		Args:
			anm_new: the new ANM
			reference: the protein the ANM was created on
			ref_chain: the matched chains of reference
		Result:
			replaced self._anm_reference and self._anm_reference_slc based on normalized anm_new
		"""
		self._anm_reference = anm_new

		self._anm_reference = self.getNormalizedANM(self._anm_reference)
		# Extend the self._anm_reference_tilde on all atoms
		anm_reference_extend = extendModel(self._anm_reference[0], self._anm_reference[1], reference, norm=True)        
		# Then slice the anm_reference to the matched
		if ref_chain != None:
			self._anm_reference_slc = sliceModel(anm_reference_extend[0], anm_reference_extend[1], ref_chain.getSelstr())
		else:
			self._anm_reference_slc = sliceModel(anm_reference_extend[0], anm_reference_extend[1], 'calpha')
		self._anm_reference_slc = self.getNormalizedANM(self._anm_reference_slc)
		
	def replaceComplexANMs(self, anm_new, proteinComplex, complex_chain = None):
		""" Replace the anm of the complex with anm_new and normalize along the way. 
		
		Args:
			anm_new: the new ANM
			proteinComplex: the complex that the ANM was created on
			complex_chain: the matched chains of the complex
		Result:
			replaced self._anm_complex and self._anm_complex_slc based on normalized anm_new
		"""
		self._anm_complex = anm_new
		self._anm_complex = self.getNormalizedANM(self._anm_complex)
		# Extend the self.self._anm_complex_tilde on all atoms
		anm_complex_extend = extendModel(self._anm_complex[0], self._anm_complex[1], proteinComplex, norm=True)        
		# Then slice the anm_reference to the matched
		if complex_chain != None:
			self._anm_complex_slc = sliceModel(anm_complex_extend[0], anm_complex_extend[1], complex_chain.getSelstr())
		else:
			self._anm_complex_slc = sliceModel(anm_complex_extend[0], anm_complex_extend[1], complex_chain.getSelstr())			
		self._anm_complex_slc = self.getNormalizedANM(self._anm_complex_slc)        
		
	def calcANMSlcInterface(self, ref_chain_interface, reference, titleOfReferenceSingleProtein, isBoundComplex=False):
		self._anm_slc_interface = self.getSlicedInterfaceANM(self.getANMExtend(), ref_chain_interface, reference, titleOfReferenceSingleProtein, isBoundComplex)
				
	def getSlicedInterfaceANM(self, anm_ext, ref_chain_interface, reference, titleOfReferenceSingleProtein, isBoundComplex=False):
		selectionAtoms = self.createSlcSelectionString(reference, isBoundComplex, ref_chain_interface, titleOfReferenceSingleProtein)
		anm_slc_interface = sliceModel(anm_ext[0], anm_ext[1], selectionAtoms)
		return anm_slc_interface
				
	def calcInterfaceANMsforPart2a2k(self, encounter):
		self._anm_reference_slc_interface = self._slicedInterfaceANMs(self._anm_reference, encounter.getMobile(), encounter.getMobChainInterface())
		self._anm_counterpart_slc_interface = self._slicedInterfaceANMs(self._anm_counterpart, encounter.getBoundCounterpart(), encounter.getBoundCounterpartChainInterface())
		self._anm_boundcomplex_slc_interface = self._slicedInterfaceANMs(self._anm_complex, encounter.boundComplex.complex , encounter.getBoundComplexChainInterface())
		
		assert (self._anm_reference_slc_interface[1].numAtoms() 
				+ self._anm_counterpart_slc_interface[1].numAtoms() 
				== self._anm_boundcomplex_slc_interface[1].numAtoms())

		for i in range(0, self._anm_reference_slc_interface[1].numAtoms()):
			assert self._anm_reference_slc_interface[1][i].getResname() == self._anm_boundcomplex_slc_interface[1][i].getResname()
			assert np.alltrue(self._anm_reference_slc_interface[1][i].getCoords() == self._anm_boundcomplex_slc_interface[1][i].getCoords())
			assert self._anm_reference_slc_interface[1][i].getName() == self._anm_boundcomplex_slc_interface[1][i].getName()
			
		offsetAtoms = self._anm_reference_slc_interface[1].numAtoms()
		for i in range(0, self._anm_counterpart_slc_interface[1].numAtoms()):
			j = i + offsetAtoms
			assert self._anm_counterpart_slc_interface[1][i].getResname() == self._anm_boundcomplex_slc_interface[1][j].getResname()
			assert np.alltrue(self._anm_counterpart_slc_interface[1][i].getCoords() == self._anm_boundcomplex_slc_interface[1][j].getCoords())
			assert self._anm_counterpart_slc_interface[1][i].getName() == self._anm_boundcomplex_slc_interface[1][j].getName()
			
	def calcInterfaceANMsUnified(self, reference, counterpart, proteinComplex, ref_chain_interface, counterpart_chain_interface, complex_chain_interface):
		""" Calculate (slice) the ANMs according to the interfaces on prot1, prot2 and their complex representation.
		
		Args:
			reference: prot1
			counterpart: prot2
			proteinComplex: prot1 and prot2 as one parsed object
			ref_chain_interface: interface of prot1
			counterpart_chain_interface: interface of prot2
			complex_chain_interface: interface of the proteinComplex
		
		"""
		self._anm_reference_slc_interface = self._slicedInterfaceANMs(self._anm_reference, reference, ref_chain_interface)
		self._anm_counterpart_slc_interface = self._slicedInterfaceANMs(self._anm_counterpart, counterpart, counterpart_chain_interface)
		self._anm_boundcomplex_slc_interface = self._slicedInterfaceANMs(self._anm_complex, proteinComplex, complex_chain_interface)
		# normalize modes
		self._anm_reference_slc_interface = self.getNormalizedANM(self._anm_reference_slc_interface)
		self._anm_counterpart_slc_interface = self.getNormalizedANM(self._anm_counterpart_slc_interface)
		self._anm_boundcomplex_slc_interface = self.getNormalizedANM(self._anm_boundcomplex_slc_interface)
		
		assert (self._anm_reference_slc_interface[1].numAtoms() 
				+ self._anm_counterpart_slc_interface[1].numAtoms() 
				== self._anm_boundcomplex_slc_interface[1].numAtoms())
		
		assertANMAtomEquality = False
		if assertANMAtomEquality:
			if self.utils.isReceptor(reference.getTitle()):
				for i in range(0, self._anm_reference_slc_interface[1].numAtoms()):
		#             print i, self._anm_reference_slc_interface[1][i].getCoords(), self._anm_boundcomplex_slc_interface[1][i].getCoords()
					assert self._anm_reference_slc_interface[1][i].getResname() == self._anm_boundcomplex_slc_interface[1][i].getResname()
					assert np.alltrue(self._anm_reference_slc_interface[1][i].getCoords() == self._anm_boundcomplex_slc_interface[1][i].getCoords())
	#                 item1roundedChoords = [round(x, 3) for x in self._anm_reference_slc_interface[1][i].getCoords().tolist()] 
	#                 item2roundedChoords = [round(x, 3) for x in self._anm_boundcomplex_slc_interface[1][i].getCoords().tolist()] 
	#                 assert  np.alltrue(item1roundedChoords == item2roundedChoords)   
					assert self._anm_reference_slc_interface[1][i].getName() == self._anm_boundcomplex_slc_interface[1][i].getName()
					
				offsetAtoms = self._anm_reference_slc_interface[1].numAtoms()
				for i in range(0, self._anm_counterpart_slc_interface[1].numAtoms()):
					j = i + offsetAtoms
					assert self._anm_counterpart_slc_interface[1][i].getResname() == self._anm_boundcomplex_slc_interface[1][j].getResname()
					assert np.alltrue(self._anm_counterpart_slc_interface[1][i].getCoords() == self._anm_boundcomplex_slc_interface[1][j].getCoords())
	#                 item1roundedChoords = [round(x, 3) for x in self._anm_counterpart_slc_interface[1][i].getCoords().tolist()] 
	#                 item2roundedChoords = [round(x, 3) for x in self._anm_boundcomplex_slc_interface[1][j].getCoords().tolist()] 
	#                 assert  np.alltrue(item1roundedChoords == item2roundedChoords)             
					assert self._anm_counterpart_slc_interface[1][i].getName() == self._anm_boundcomplex_slc_interface[1][j].getName()
	
			else:
				offsetAtoms = self._anm_counterpart_slc_interface[1].numAtoms()
				for i in range(0, self._anm_reference_slc_interface[1].numAtoms()):
					j = i + offsetAtoms
		#             print i, self._anm_reference_slc_interface[1][i].getCoords(), self._anm_boundcomplex_slc_interface[1][i].getCoords()
					assert self._anm_reference_slc_interface[1][i].getResname() == self._anm_boundcomplex_slc_interface[1][j].getResname()
					assert np.alltrue(self._anm_reference_slc_interface[1][i].getCoords() == self._anm_boundcomplex_slc_interface[1][j].getCoords())
	#                 item1roundedChoords = [round(x, 3) for x in self._anm_reference_slc_interface[1][i].getCoords().tolist()] 
	#                 item2roundedChoords = [round(x, 3) for x in self._anm_boundcomplex_slc_interface[1][j].getCoords().tolist()] 
	#                 assert  np.alltrue(item1roundedChoords == item2roundedChoords)   
					assert self._anm_reference_slc_interface[1][i].getName() == self._anm_boundcomplex_slc_interface[1][j].getName()
					 
				for i in range(0, self._anm_counterpart_slc_interface[1].numAtoms()):
					assert self._anm_counterpart_slc_interface[1][i].getResname() == self._anm_boundcomplex_slc_interface[1][i].getResname()
					assert np.alltrue(self._anm_counterpart_slc_interface[1][i].getCoords() == self._anm_boundcomplex_slc_interface[1][i].getCoords())
	#                 item1roundedChoords = [round(x, 3) for x in self._anm_counterpart_slc_interface[1][i].getCoords().tolist()] 
	#                 item2roundedChoords = [round(x, 3) for x in self._anm_boundcomplex_slc_interface[1][i].getCoords().tolist()] 
	#                 assert  np.alltrue(item1roundedChoords == item2roundedChoords)             
					assert self._anm_counterpart_slc_interface[1][i].getName() == self._anm_boundcomplex_slc_interface[1][i].getName()
		
	def _slicedInterfaceANMs(self, anm, reference, interface):
		""" Slice an anm to match the provided interface. 
		
		Args:
			anm: the anm to be sliced
			reference: the protein that the anm is based upon, necessary for extention of the model first
			interface: the interface of the protein
		"""
		anm_ext = extendModel(anm[0], anm[1], reference, norm=True)
		anm_slc = sliceModel(anm_ext[0], anm_ext[1], interface.getSelstr())
		anm_slc = self.getNormalizedANM(anm_slc)
		return anm_slc
	
	def getANM(self):
		""" Get the ANM calculated on the reference (default) calpha atoms. """
		if self._anm == None:
			raise Exception('self._anm == None')
		return self._anm
	
	def getANMExtend(self):
		""" Get the ANM extended to the whole reference (all atoms). """
		if self._anm_extend == None:
			raise Exception('self._anm == None')
		return self._anm_extend
	
	def getANMSlc(self):
		""" Get the sliced back ANM to match all atoms in the ref_chain."""
		if self._anm_slc == None:
			raise Exception('self._anm_slc == None')
		return self._anm_slc
	
	def getANMSlcCounterpart(self):
		""" Get the sliced back ANM to match all atoms in the counterpart chain(s) """
		if self._anm_slc_counterpart == None:
			raise Exception('self._anm_slc == None')
		return self._anm_slc_counterpart
	
	def getANMSlcInterface(self):
		""" Get the sliced back ANM to match all atoms in the ref_chain_interface. """
		if self._anm_slc_interface == None:
			raise Exception('self._anm_slc_interface == None')
		return self._anm_slc_interface
				
	def getANMComplexSlc(self):
		""" Get the sliced back ANM to match all atoms in the chain_complex. """
		if self._anm_complex_slc == None:
			raise Exception('self._anm_complex_slc == None')
		return self._anm_complex_slc
	
	def getANMReference2a2kSlc(self):
		""" Get the sliced back self._anm_reference_slc ANM to match all atoms in the reference variable. """
		if self._anm_reference_slc == None:
			raise Exception('self._anm_reference_slc == None')
		return self._anm_reference_slc
	
	def getANMCounterpart2a2kSlc(self):
		""" Get the sliced back self._anm_counterpart_slc ANM to match all atoms in the counterpart variable. """
		if self._anm_counterpart_slc == None:
			raise Exception('self._anm_counterpart_slc == None')
		return self._anm_counterpart_slc
	  
	def getANMReference(self):
		if self._anm_reference == None:
			raise Exception('self._anm_reference == None')
		return self._anm_reference
	
	def getANMReferenceSlc(self):
		if self._anm_reference_slc == None:
			raise Exception('self._anm_reference_slc == None')
		return self._anm_reference_slc
	
	def getANMCounterpart(self):
		if self._anm_counterpart == None:
			raise Exception('self._anm_counterpart == None')
		return self._anm_counterpart

	def getANMComplex(self):
		if self._anm_complex == None:
			raise Exception('self._anm_complex == None')
		return self._anm_complex

	def getANMReferenceSlcInterface(self):
		if self._anm_reference_slc_interface == None:
			raise Exception('self._anm_reference_slc_interface == None')
		return self._anm_reference_slc_interface
	
	def getANMCounterpartSlcInterface(self):
		if self._anm_counterpart_slc_interface == None:
			raise Exception('self._anm_counterpart_slc_interface == None')
		return self._anm_counterpart_slc_interface
	
	def getANMComplexSlcInterface(self):
		if self._anm_boundcomplex_slc_interface == None:
			raise Exception('self._anm_boundcomplex_slc_interface == None')
		return self._anm_boundcomplex_slc_interface
				
	def getANMPath(self, reference, numberOfModes, selstr, whatAtomsToMatch, modified=""):
		path = self.utils.config.anmPath
		prefix = reference.getTitle()
		prefix = prefix.replace(" ", "_")
		if modified == "":
			return path+prefix+"_modes"+str(numberOfModes)+"_buildOn"+selstr+"_matchedOn"+whatAtomsToMatch
		elif modified == "extended":
			return path+"extended/"+prefix+"_modes"+str(numberOfModes)+"_buildOn"+selstr+"_matchedOn"+whatAtomsToMatch+"_extended"
		elif modified == "slicedback":
			return path+"slicedback/"+prefix+"_modes"+str(numberOfModes)+"_buildOn"+selstr+"_matchedOn"+whatAtomsToMatch+"_slicedback"
		else:
			raise Exception("the variable modified is not the empty string, extended or slicedback.")
	
	def doesANMExist(self, reference, numberOfModes, selstr, whatAtomsToMatch, modified=""):
		path = self.utils.config.anmPath
		try:
			with open(self.getANMPath(reference, numberOfModes, selstr, whatAtomsToMatch, modified)+".anm.npz"):
				return True
		except IOError:
			return False

	def projectHessian(self, hessian, reference, proteinComplex, referenceSegment, projectionStyle, projectOnlyReferencePartOfHC=False, interCalphaIndices=None):
		""" Return the PH'P hessian which has 6 zero eigenvalues according to the formula 8.27
		from the book "A practical introduction to the simulation of molecular dynamics", Field. 
		However, here it is made sure that the assumed basis is orthonormal via np.linalg.qr applied 
		on the six vectors discussed in this book.
		
		Args:
			hessian: the hessian to be projected
			reference: the protein the hessian or HRtilde/HLtilde of the hessian was created on
			proteinComplex: the whole protein that reference is part of
			referenceSegment: if reference is receptor, provide "R", else it needs to be ligand, provide "L"
			projectionStyle: project away from "full" (intra+inter) or "intra" (intra) or "fullComplex"
			pojectOnlyReferencePartOfHC: if true, the hessian was created on reference, if false, HRtilde or HLtilde 
										 of the hessian were created on the reference
			interCalphaIndices: list of calphas indices that have intermolecular interactions  
			
		Returns: projected hessian with 6 external degrees of freedom (rotation and translation) removed
		"""
		assert projectionStyle == "full" or projectionStyle == "intra" or projectionStyle == "fullComplex"
		normalize = True
		numAtoms = reference.select('calpha').numAtoms()
		numCoords = numAtoms*3
		centerOfCoords = calcCenter(reference.select('calpha')) 

		assert np.alltrue(centerOfCoords == calcCenter(proteinComplex.select('segment \"'+referenceSegment+'.\"').select('calpha')))
		print "before projection symmetry ==, allclose: ", np.all(hessian-hessian.T==0), np.allclose(hessian, hessian.T)
	   
		if projectOnlyReferencePartOfHC:
			numComplexAtoms = proteinComplex.select('calpha').numAtoms()
			numComplexCoords = numComplexAtoms*3
			numCounterpartCoords = numComplexCoords - numCoords
			
			if referenceSegment == "R":
				assert numCounterpartCoords == proteinComplex.select('segment \"L.\"').select('calpha').numAtoms() * 3
			else:
				assert numCounterpartCoords == proteinComplex.select('segment \"R.\"').select('calpha').numAtoms() * 3
			
			# Create null vector with length of the counterpart calphas
			counterPartNullVector = np.zeros(numCounterpartCoords)   
	   
		# Create I
		I = np.identity(numCoords)
		# Create the three translation vectors Tx, Ty, Tz
		Tx = np.zeros(numCoords)
		Tx = self.utils.fill3DArrayWithValue(Tx, 1.0, 0)
		Ty = np.zeros(numCoords)
		Ty = self.utils.fill3DArrayWithValue(Ty, 1.0, 1)
		Tz = np.zeros(numCoords)
		Tz = self.utils.fill3DArrayWithValue(Tz, 1.0, 2)
		# Create the three rotation vectors Rx, Ry, Rz
		coordsCopy = reference.select('calpha').getCoords().copy()
		Rx = self.utils.createRx(coordsCopy)
		coordsCopy2 = reference.select('calpha').getCoords().copy()
		Ry = self.utils.createRy(coordsCopy2)
		coordsCopy3 = reference.select('calpha').getCoords().copy()
		Rz = self.utils.createRz(coordsCopy3)
		
		# remove inter atoms from projection
		if projectionStyle == "intra":
			Tx = self.removeInterAtoms(Tx, interCalphaIndices)
			Ty = self.removeInterAtoms(Ty, interCalphaIndices)
			Tz = self.removeInterAtoms(Tz, interCalphaIndices)
			Rx = self.removeInterAtoms(Rx, interCalphaIndices)
			Ry = self.removeInterAtoms(Ry, interCalphaIndices)
			Rz = self.removeInterAtoms(Rz, interCalphaIndices)
   
		if projectOnlyReferencePartOfHC:
			# overwrite previous I
			I = np.identity(numComplexCoords)
			# extend (with the nullvector) the rotational and translational vectors to the dimension of the complex
			if referenceSegment == "R":
				Tx = np.concatenate((Tx, counterPartNullVector))
				Ty = np.concatenate((Ty, counterPartNullVector))
				Tz = np.concatenate((Tz, counterPartNullVector))
				Rx = np.concatenate((Rx, counterPartNullVector))
				Ry = np.concatenate((Ry, counterPartNullVector))
				Rz = np.concatenate((Rz, counterPartNullVector))
			else:
				Tx = np.concatenate((counterPartNullVector, Tx))
				Ty = np.concatenate((counterPartNullVector, Tz))
				Tz = np.concatenate((counterPartNullVector, Tz))
				Rx = np.concatenate((counterPartNullVector, Rx))
				Ry = np.concatenate((counterPartNullVector, Ry))
				Rz = np.concatenate((counterPartNullVector, Rz))                
   
		# Normalize translation vectors and apply rotational fix
		if normalize:
			Tx = Vector(Tx)
			#Tx = self.subtractCenterOfCoords(Tx, centerOfCoords[0], 0.0, 0.0)
			Tx = Tx.getNormed().getArray()
			
			Ty = Vector(Ty)
			#Ty = self.subtractCenterOfCoords(Ty, 0.0, centerOfCoords[1], 0.0)
			Ty = Ty.getNormed().getArray()
			
			Tz = Vector(Tz)
			#Tz = self.subtractCenterOfCoords(Tz, 0.0, 0.0, centerOfCoords[2])
			Tz = Tz.getNormed().getArray() 
						
			Rx = Vector(Rx)
			#Rx = self.subtractCenterOfCoords(Rx, 0.0, centerOfCoords[2], centerOfCoords[1])
			Rx = Rx.getNormed().getArray()
			
			Ry = Vector(Ry)
			#Ry = self.subtractCenterOfCoords(Ry, centerOfCoords[2], 0.0, centerOfCoords[0])
			Ry = Ry.getNormed().getArray()  
			
			Rz = Vector(Rz)
			#Rz = self.subtractCenterOfCoords(Rz, centerOfCoords[1], centerOfCoords[0], 0.0)
			Rz = Rz.getNormed().getArray()
			
		# Create P
		#P = I - np.outer(Rx, Rx) - np.outer(Ry, Ry) - np.outer(Rz, Rz) - np.outer(Tx, Tx) - np.outer(Ty, Ty) - np.outer(Tz, Tz)
		### corres P = I - P
		#print "independent columns P: ", self.utils.independent_columns(P).shape
		#print "matrix rank P: ", self.utils.matrixrank(P)
		#print "independent columns I-P: ", self.utils.independent_columns(I-P).shape
		#print "matrix rank I-P: ", self.utils.matrixrank(I-P)
		#print "np matrix rank I-P : ", np.linalg.matrix_rank(I-P)
		#print "np matrix as matrix rank I-P : ", np.linalg.matrix_rank(np.matrix(I-P))
		
		assumedBasis = np.array([Tx, Ty, Tz, Rx, Ry, Rz]).T
		MyQ, MyR = np.linalg.qr(assumedBasis)
		#print "MyQ.shape: ", MyQ.shape
		Rx = MyQ.T[0]
		Ry = MyQ.T[1]
		Rz = MyQ.T[2]
		Tx = MyQ.T[3]
		Ty = MyQ.T[4]
		Tz = MyQ.T[5]
		###
		print "before full projection"
		###
		P = I - np.outer(Rx, Rx) - np.outer(Ry, Ry) - np.outer(Rz, Rz) - np.outer(Tx, Tx) - np.outer(Ty, Ty) - np.outer(Tz, Tz)
		#print "assumedBasis : \n", assumedBasis.round(4)
		#print "basis after QR: \n", np.array([Tx, Ty, Tz, Rx, Ry, Rz]).T.round(4)
		#writeArray("assumedBasis.txt", assumedBasis.round(4), format="%f")
		#writeArray("basis_after_QR.txt", np.array([Tx, Ty, Tz, Rx, Ry, Rz]).T.round(4), format="%f")

		###
		#print "P", P
#         print "P.shape", P.shape
#         print "symmetric P: ", np.allclose(P, P.T)
#         print "complex calphas * 3: ", proteinComplex.select('calpha').numAtoms() * 3
#         print "rank of P projection", projectionStyle, ": ", np.linalg.matrix_rank(np.matrix(P)) 
#         P_eigenvals, P_eigenvecs = np.linalg.eigh(P)
#         print "number of P_eigenvals: ", len(P_eigenvals)
#         #print "P_eigenvals: ", P_eigenvals
#         print "number of P_eigenvecs: ", len(P_eigenvecs)
#         #print "P_eigenvecs: ", P_eigenvecs       
#         #writeArray("helperScripts/"+proteinComplex.getTitle()+"_P_"+projectionStyle+".txt", P, format='%10.7f')
#         #writeArray("P_eigenvals"+projectionStyle+".txt", P_eigenvals, format='%10.7f')
#         #writeArray("P_eigenvecs"+projectionStyle+".txt", P_eigenvecs, format='%10.7f')
#         
#         P_times_Peigenvecs = P.dot(P_eigenvecs)
#         P_times_Peigenvecs_T = P.dot(P_eigenvecs).T
#         P_orthonormalityTest = P_times_Peigenvecs_T.dot(P_times_Peigenvecs)
#         #writeArray("P_orthonormalityTest"+projectionStyle+".txt", P_orthonormalityTest, format='%10.7f')
#         # does this P_orthonormalityTest equal the identity matrix or part of it?
#         print "P_orthonormalityTest: ", np.allclose(P_orthonormalityTest, np.identity(len(P_eigenvecs)))
#         print "P_orthonormalityTest w/o upper 6x6: ", np.allclose(P_orthonormalityTest[6:,6:], np.identity(len(P_eigenvecs)-6))
#         zeroM = np.zeros((len(P_eigenvecs), len(P_eigenvecs)))
#         zeroM[6:,6:] = P_orthonormalityTest[6:,6:]
#         print "P_orthonormalityTest except lower n-6,n-6 zero: ", np.allclose(P_orthonormalityTest, zeroM)
		
#         proteinComplex_ca = proteinComplex.select('calpha')
#         writePDB("complex_allatoms.pdb", proteinComplex)
#         writePDB("complex_before_Ptimes.pdb", proteinComplex_ca)
#         coord_shape = proteinComplex_ca.getCoords().shape
#         coords_P = P.dot(proteinComplex_ca.getCoords().flatten())
#         coords_P = coords_P.reshape(coord_shape)
#         proteinComplex_ca.setCoords(coords_P)
#         writePDB("complex_after_Ptimes"+projectionStyle+".pdb", proteinComplex_ca)
		#raw_input()
		###
#         Q, R = np.linalg.qr(P, mode="complete")
#         print "independent columns Q: ", self.utils.independent_columns(Q).shape
#         print "matrix rank Q: ", self.utils.matrixrank(Q)
#         print "matrix np rank Q: ", np.linalg.matrix_rank(Q)," ", np.linalg.matrix_rank(np.matrix(Q))
#         print "log of determinant of Q: ", np.linalg.slogdet(Q)
		### corres Q = I - Q
		#P = I-Q
		# Apply PH'H, np.dot is matrix multiplication for 2D arrays
		#print "count orthogonal columns: ", self.utils.countOrthogonalColumns(I-P)
		Hprime = np.dot(P.T, hessian)
		Hprime = np.dot(Hprime, P)
		# Return the projected hessian
		#print "after projection symmetry ==, allclose: ", np.all(Hprime-Hprime.T==0), np.allclose(Hprime, Hprime.T)
		#print "H: ", hessian
		#print "Hprime: ", Hprime
		return Hprime
		
	
	def projectHessian_test2timesQR(self, hessian, reference, proteinComplex, referenceSegment, projectionStyle, projectOnlyReferencePartOfHC=False, interCalphaIndices=None):
		""" Return the PH'P hessian which has 6 zero eigenvalues according to the formula 8.27
		from the book "A practical introduction to the simulation of molecular dynamics", Field. 
		However, here it is made sure that the assumed basis is orthonormal via np.linalg.qr applied 
		on the six vectors discussed in this book.
		
		Args:
			hessian: the hessian to be projected
			reference: the protein the hessian or HRtilde/HLtilde of the hessian was created on
			proteinComplex: the whole protein that reference is part of
			referenceSegment: if reference is receptor, provide "R", else it needs to be ligand, provide "L"
			projectionStyle: project away from "full" (intra+inter) or "intra" (intra) or "fullComplex"
			pojectOnlyReferencePartOfHC: if true, the hessian was created on reference, if false, HRtilde or HLtilde 
										 of the hessian were created on the reference
			interCalphaIndices: list of calphas indices that have intermolecular interactions  
			
		Returns: projected hessian with 6 external degrees of freedom (rotation and translation) removed
		"""
		assert projectionStyle == "full"
		normalize = True
		numAtoms = reference.select('calpha').numAtoms()
		numCoords = numAtoms*3
		centerOfCoords = calcCenter(reference.select('calpha')) 

		assert np.alltrue(centerOfCoords == calcCenter(proteinComplex.select('segment \"'+referenceSegment+'.\"').select('calpha')))
		print "before projection symmetry ==, allclose: ", np.all(hessian-hessian.T==0), np.allclose(hessian, hessian.T)
	   
		numComplexAtoms = proteinComplex.select('calpha').numAtoms()
		numComplexCoords = numComplexAtoms*3
		numCounterpartCoords = numComplexCoords - numCoords
			
		if referenceSegment == "R":
			assert numCounterpartCoords == proteinComplex.select('segment \"L.\"').select('calpha').numAtoms() * 3
		else:
			assert numCounterpartCoords == proteinComplex.select('segment \"R.\"').select('calpha').numAtoms() * 3
			
		# Create null vector with length of the counterpart calphas
		counterPartNullVector = np.zeros(numCounterpartCoords)   
	   
		# Create I
		I = np.identity(numComplexCoords)
		# Create the three translation vectors Tx, Ty, Tz
		Tx = np.zeros(numComplexCoords)
		Tx = self.utils.fill3DArrayWithValue(Tx, 1.0, 0)
		Ty = np.zeros(numComplexCoords)
		Ty = self.utils.fill3DArrayWithValue(Ty, 1.0, 1)
		Tz = np.zeros(numComplexCoords)
		Tz = self.utils.fill3DArrayWithValue(Tz, 1.0, 2)
		# Create the three rotation vectors Rx, Ry, Rz
		coordsCopy = proteinComplex.select('calpha').getCoords().copy()
		Rx = self.utils.createRx(coordsCopy)
		coordsCopy2 = proteinComplex.select('calpha').getCoords().copy()
		Ry = self.utils.createRy(coordsCopy2)
		coordsCopy3 = proteinComplex.select('calpha').getCoords().copy()
		Rz = self.utils.createRz(coordsCopy3)
   
#         if projectOnlyReferencePartOfHC:
#             # overwrite previous I
#             I = np.identity(numComplexCoords)
#             # extend (with the nullvector) the rotational and translational vectors to the dimension of the complex
#             if referenceSegment == "R":
#                 Tx = np.concatenate((Tx, counterPartNullVector))
#                 Ty = np.concatenate((Ty, counterPartNullVector))
#                 Tz = np.concatenate((Tz, counterPartNullVector))
#                 Rx = np.concatenate((Rx, counterPartNullVector))
#                 Ry = np.concatenate((Ry, counterPartNullVector))
#                 Rz = np.concatenate((Rz, counterPartNullVector))
#             else:
#                 Tx = np.concatenate((counterPartNullVector, Tx))
#                 Ty = np.concatenate((counterPartNullVector, Tz))
#                 Tz = np.concatenate((counterPartNullVector, Tz))
#                 Rx = np.concatenate((counterPartNullVector, Rx))
#                 Ry = np.concatenate((counterPartNullVector, Ry))
#                 Rz = np.concatenate((counterPartNullVector, Rz))                
   
		# Normalize translation vectors and apply rotational fix
		if normalize:
			Tx = Vector(Tx)
			#Tx = self.subtractCenterOfCoords(Tx, centerOfCoords[0], 0.0, 0.0)
			Tx = Tx.getNormed().getArray()
			
			Ty = Vector(Ty)
			#Ty = self.subtractCenterOfCoords(Ty, 0.0, centerOfCoords[1], 0.0)
			Ty = Ty.getNormed().getArray()
			
			Tz = Vector(Tz)
			#Tz = self.subtractCenterOfCoords(Tz, 0.0, 0.0, centerOfCoords[2])
			Tz = Tz.getNormed().getArray() 
						
			Rx = Vector(Rx)
			#Rx = self.subtractCenterOfCoords(Rx, 0.0, centerOfCoords[2], centerOfCoords[1])
			Rx = Rx.getNormed().getArray()
			
			Ry = Vector(Ry)
			#Ry = self.subtractCenterOfCoords(Ry, centerOfCoords[2], 0.0, centerOfCoords[0])
			Ry = Ry.getNormed().getArray()  
			
			Rz = Vector(Rz)
			#Rz = self.subtractCenterOfCoords(Rz, centerOfCoords[1], centerOfCoords[0], 0.0)
			Rz = Rz.getNormed().getArray()
		
		assumedBasis = np.array([Tx, Ty, Tz, Rx, Ry, Rz]).T
		MyQ, MyR = np.linalg.qr(assumedBasis, mode='full')
		Rx = MyQ.T[0]
		Ry = MyQ.T[1]
		Rz = MyQ.T[2]
		Tx = MyQ.T[3]
		Ty = MyQ.T[4]
		Tz = MyQ.T[5]

		Rx = Rx[:numCoords]
		Ry = Ry[:numCoords]
		Rz = Rz[:numCoords]
		Tx = Tx[:numCoords]
		Ty = Ty[:numCoords]
		Tz = Tz[:numCoords]
		
		assumedBasis = np.array([Tx, Ty, Tz, Rx, Ry, Rz]).T
		MyQ, MyR = np.linalg.qr(assumedBasis, mode='full')
		Rx = MyQ.T[0]
		Ry = MyQ.T[1]
		Rz = MyQ.T[2]
		Tx = MyQ.T[3]
		Ty = MyQ.T[4]
		Tz = MyQ.T[5]             
		
		print "len(Rx): ", len(Rx)
		
		Tx = np.concatenate((Tx, counterPartNullVector))
		Ty = np.concatenate((Ty, counterPartNullVector))
		Tz = np.concatenate((Tz, counterPartNullVector))
		Rx = np.concatenate((Rx, counterPartNullVector))
		Ry = np.concatenate((Ry, counterPartNullVector))
		Rz = np.concatenate((Rz, counterPartNullVector))        
		
		print "Pr test"
		raw_input()
		
		P = I - np.outer(Rx, Rx) - np.outer(Ry, Ry) - np.outer(Rz, Rz) - np.outer(Tx, Tx) - np.outer(Ty, Ty) - np.outer(Tz, Tz)
		#print "assumedBasis : \n", assumedBasis.round(4)
		#print "basis after QR: \n", np.array([Tx, Ty, Tz, Rx, Ry, Rz]).T.round(4)
		#writeArray("assumedBasis.txt", assumedBasis.round(4), format="%f")
		#writeArray("basis_after_QR.txt", np.array([Tx, Ty, Tz, Rx, Ry, Rz]).T.round(4), format="%f")

		###
		print "P", P
		print "P.shape", P.shape
		print "symmetric P: ", np.allclose(P, P.T)
		print "complex calphas * 3: ", proteinComplex.select('calpha').numAtoms() * 3
		print "rank of P projection", projectionStyle, ": ", np.linalg.matrix_rank(np.matrix(P)) 
		P_eigenvals, P_eigenvecs = np.linalg.eigh(P)
		print "number of P_eigenvals: ", len(P_eigenvals)
		#print "P_eigenvals: ", P_eigenvals
		print "number of P_eigenvecs: ", len(P_eigenvecs)
		#print "P_eigenvecs: ", P_eigenvecs       
		writeArray("helperScripts/"+proteinComplex.getTitle()+"_P_"+projectionStyle+".txt", P, format='%10.7f')
		#writeArray("P_eigenvals"+projectionStyle+".txt", P_eigenvals, format='%10.7f')
		#writeArray("P_eigenvecs"+projectionStyle+".txt", P_eigenvecs, format='%10.7f')
		
		P_times_Peigenvecs = P.dot(P_eigenvecs)
		P_times_Peigenvecs_T = P.dot(P_eigenvecs).T
		P_orthonormalityTest = P_times_Peigenvecs_T.dot(P_times_Peigenvecs)
		#writeArray("P_orthonormalityTest"+projectionStyle+".txt", P_orthonormalityTest, format='%10.7f')
		# does this P_orthonormalityTest equal the identity matrix or part of it?
		print "P_orthonormalityTest: ", np.allclose(P_orthonormalityTest, np.identity(len(P_eigenvecs)))
		print "P_orthonormalityTest w/o upper 6x6: ", np.allclose(P_orthonormalityTest[6:,6:], np.identity(len(P_eigenvecs)-6))
		zeroM = np.zeros((len(P_eigenvecs), len(P_eigenvecs)))
		zeroM[6:,6:] = P_orthonormalityTest[6:,6:]
		print "P_orthonormalityTest except lower n-6,n-6 zero: ", np.allclose(P_orthonormalityTest, zeroM)
		
#         proteinComplex_ca = proteinComplex.select('calpha')
#         writePDB("complex_allatoms.pdb", proteinComplex)
#         writePDB("complex_before_Ptimes.pdb", proteinComplex_ca)
#         coord_shape = proteinComplex_ca.getCoords().shape
#         coords_P = P.dot(proteinComplex_ca.getCoords().flatten())
#         coords_P = coords_P.reshape(coord_shape)
#         proteinComplex_ca.setCoords(coords_P)
#         writePDB("complex_after_Ptimes"+projectionStyle+".pdb", proteinComplex_ca)
		raw_input()
		###
#         Q, R = np.linalg.qr(P, mode="complete")
#         print "independent columns Q: ", self.utils.independent_columns(Q).shape
#         print "matrix rank Q: ", self.utils.matrixrank(Q)
#         print "matrix np rank Q: ", np.linalg.matrix_rank(Q)," ", np.linalg.matrix_rank(np.matrix(Q))
#         print "log of determinant of Q: ", np.linalg.slogdet(Q)
		### corres Q = I - Q
		#P = I-Q
		# Apply PH'H, np.dot is matrix multiplication for 2D arrays
		#print "count orthogonal columns: ", self.utils.countOrthogonalColumns(I-P)
		Hprime = np.dot(P.T, hessian)
		Hprime = np.dot(Hprime, P)
		# Return the projected hessian
		#print "after projection symmetry ==, allclose: ", np.all(Hprime-Hprime.T==0), np.allclose(Hprime, Hprime.T)
		#print "H: ", hessian
		#print "Hprime: ", Hprime
		return Hprime    
	
	def transformHessianToFixedDomainFrame(self, hessian, reference, proteinComplex, referenceSegment, projectionStyle):
		""" Application of formula 20 from:
		Fuchigami, Sotaro, Satoshi Omori, Mitsunori Ikeguchi, and Akinori Kidera. 
		"Normal Mode Analysis of Protein Dynamics in a Non-Eckart Frame." 
		The Journal of Chemical Physics 132, no. 10 (March 11, 2010): 104109. doi:10.1063/1.3352566. 
		"""
		
		numAtoms = reference.select('calpha').numAtoms()
		numCoords = numAtoms*3
		centerOfCoords = calcCenter(reference.select('calpha')) 

		#assert np.alltrue(centerOfCoords == calcCenter(proteinComplex.select('segment \"'+referenceSegment+'.\"').select('calpha')))

		numComplexAtoms = proteinComplex.select('calpha').numAtoms()
		numComplexCoords = numComplexAtoms*3
		numCounterpartCoords = numComplexCoords - numCoords
		
		if referenceSegment == "R":
			# create the P matrix, receptor is fixed domain
			P = np.zeros((numComplexCoords, numComplexCoords))
			P[:numCoords, :numCoords] = np.identity(numCoords)            
			assert numCounterpartCoords == proteinComplex.select('segment \"L.\"').select('calpha').numAtoms() * 3
		else:
			# create the P matrix, ligand is fixed domain
			P = np.zeros((numComplexCoords, numComplexCoords))
			numCoords_receptor = proteinComplex.select('segment \"R.\"').select('calpha').numAtoms() * 3
			P[numCoords_receptor:, numCoords_receptor:] = np.identity(proteinComplex.select('segment \"L.\"').select('calpha').numAtoms() * 3)        
			assert numCounterpartCoords == proteinComplex.select('segment \"R.\"').select('calpha').numAtoms() * 3

		# create rigid body motion eigenvectors out_values
		out_vals, out_vectors = sp.linalg.eigh(hessian)
		
		# sort the eigenvalues and eigenvectors ascendingly, this is not asserted by the eigh return, see
		# http://stackoverflow.com/questions/8092920/sort-eigenvalues-and-associated-eigenvectors-after-using-numpy-linalg-eig-in-pyt
		idx = out_vals.argsort()   
		out_vals = out_vals[idx]
		out_vectors = out_vectors[:,idx]
		
		# take the first six eigenvalues and eigenvectors
		out_vals = out_vals[:6]
		out_vectors = out_vectors.T[:6].T   
		#print "P.shape: ", P.shape              
		#print "out_vectors.shape: ", out_vectors.shape
		
		# create the transformation matrix
		inv = (out_vectors.T.dot(P)).dot(out_vectors)
		
		inv = np.linalg.inv(inv)
		secondTerm = ((out_vectors.dot(inv)).dot(out_vectors.T)).dot(P)

		U = np.identity(numComplexCoords) - secondTerm
		print "calculated transformation matrix U"
		
		#writeArray("hessianbeforeU.txt", hessian, format='%10.7f')
		Hprime = np.dot(U, hessian)
		Hprime = np.dot(Hprime, U.T)
		#writeArray(proteinComplex.getTitle()+"U.txt", U, format='%10.7f')
		
		#writeArray("hessianafterU.txt", Hprime, format='%10.7f')
		print "obtained Hprime with a fixed domain frame"
		return Hprime    
	
	def subtractCenterOfCoords(self, vector, xElement, yElement, zElement):
		""" Subtract from a vector having a [i][3] dim array elementwise the center of coords and return the result. """
		coordsNx3 = vector.getArrayNx3()
		subtractArray = np.array([xElement, yElement, zElement])
		coordsNx3 = coordsNx3 - subtractArray
		resultVector = Vector(coordsNx3.flatten())
		return resultVector
	
	def addscaledHdelta(self, HR, HRtilde, deltaHRmultiplicator):
		assert HR.shape == HRtilde.shape
		deltaHR = HRtilde - HR
		deltaHR = deltaHR * deltaHRmultiplicator
		return (HR + deltaHR)
	
	def calcCustomH_ANew(self, HR, referenceStructure, neighborStructure, encounter, neighborhoodFrom, equilibriumAt, workOnReceptor=True, selstr='calpha'):
		""" Modifies the hessian HR or HL by adding additonal terms for intramolecular contacts.
		
			Args:
				HR: The original HR as calculated by prody
				referenceStructure: structure to take calphas from, the hessian HR belongs to it or to its superset if I is a chain
				neighborStructure: structure to apply the neighborhood calculations on
				encounter: object with all encounter information
				neighborhoodFrom: is the neighborhood calculated from the unbound complex C_u or the bound complex C_b
				equilibriumAt: is the equilibrium set to r_ij or r_ij_b
				workonReceptor: is the Hessian and the referenceStructure receptor or ligand
				selstr: atomType of the course grained ANM (by default calpha)
		"""
		assert equilibriumAt == "r_ij" or equilibriumAt == "r_ij_b"
		assert neighborhoodFrom == "C_u" or neighborhoodFrom == "C_b"
		
		if workOnReceptor:
			reference = encounter.getReference()
			if self.bound_provided == True:
				refchain = encounter.getRefChain()
				mobile = encounter.getMobile()
				mobChain = encounter.getMobChain()
				boundCounterpart = encounter.getBoundCounterpart()
				boundCounterpartChain = encounter.getBoundCounterpartChain()
				unboundCounterpartChain = encounter.getUnboundCounterpartChain()
		else:  
			reference = encounter.getUnboundCounterpart()
			if self.bound_provided == True:
				refchain = encounter.getUnboundCounterpartChain()
				mobile = encounter.getBoundCounterpart()
				mobChain = encounter.getBoundCounterpartChain()
				boundCounterpart = encounter.getMobile()
				boundCounterpartChain = encounter.getMobChain()
				unboundCounterpartChain = encounter.getRefChain()
				
		neighborStructureCalpha = neighborStructure.select('calpha')
		contactsCounter = 0
		interCalphaIndices = []
				
		for idx, element in enumerate(referenceStructure.select('calpha')):
			contactsOfI = encounter.getIntermolecularNeighborsOfAtom(element, neighborStructure, selstr, str(self.utils.config.customHRdistance))
			# if element has contacts in the neighborStructure, the hessian needs an update in the 3*3 matrix on the diagonal of this element atom
			if contactsOfI:
				contactsCounter += contactsOfI.numAtoms()
				interCalphaIndices.append(idx)
				print "intermolecular contacts: ", contactsOfI.numAtoms()
				contacts_counterpartChainIndices = self.utils.getMatchingStructureSelections(neighborStructureCalpha, contactsOfI, neighborStructureCalpha)
				assert len(contactsOfI) == len(contacts_counterpartChainIndices)            
				# access each element contact to create the deltaTerm
				overallTerm = np.zeros((3,3))
				for elementcontact, contacts_counterpartChainIndex in zip(contactsOfI, contacts_counterpartChainIndices):
					if neighborhoodFrom == "C_u":
						r_ij = calcDistance(element, elementcontact)
						if equilibriumAt == "r_ij":
							r_ij_b = r_ij
						#if element is not in matched reference or contact is not in matched counterpart: r_ij_b = r_ij
						elif not(element in refchain.select('calpha')) or not(elementcontact in unboundCounterpartChain.select('calpha')):
							r_ij_b = r_ij
						else:
							elementPositionInChain = encounter.accessANMs().getCalphaPosition(element, refchain.select('calpha'))
							contactPositionInChain = encounter.accessANMs().getCalphaPosition(elementcontact, unboundCounterpartChain.select('calpha'))
							r_ij_b = calcDistance(mobChain.select('calpha')[elementPositionInChain], boundCounterpartChain.select('calpha')[contactPositionInChain])                        
							self.utils.assertTwoAtomsAreEqual(mobChain.select('calpha')[elementPositionInChain], element, useCoords=False, useResname=True)
							self.utils.assertTwoAtomsAreEqual(refchain.select('calpha')[elementPositionInChain], element, useCoords=False, useResname=True)
							self.utils.assertTwoAtomsAreEqual(boundCounterpartChain.select('calpha')[contactPositionInChain], elementcontact, useCoords=False, useResname=True)                       
							self.utils.assertTwoAtomsAreEqual(unboundCounterpartChain.select('calpha')[contactPositionInChain], elementcontact, useCoords=False, useResname=True)
						# make the 3*3 hessian term for this contact (excluding gamma, gamma is multiplied at the end to the sum)  
						deltaTerm = self.make3By3HessianTerm(element, elementcontact, r_ij, r_ij_b)
						#print element, elementcontact, " r_ij, rij_b: ", r_ij, r_ij_b
						overallTerm += deltaTerm
					else:
						if equilibriumAt == "r_ij_b":
							r_ij_b = calcDistance(element, elementcontact)
							elementPositionInChain = encounter.accessANMs().getCalphaPosition(element, mobChain.select('calpha'))
							contactPositionInChain = encounter.accessANMs().getCalphaPosition(elementcontact, boundCounterpartChain.select('calpha'))
							r_ij = calcDistance(refchain.select('calpha')[elementPositionInChain], unboundCounterpartChain.select('calpha')[contactPositionInChain])                        
							self.utils.assertTwoAtomsAreEqual(mobChain.select('calpha')[elementPositionInChain], element, useCoords=False, useResname=True)
							self.utils.assertTwoAtomsAreEqual(refchain.select('calpha')[elementPositionInChain], element, useCoords=False, useResname=True)
							self.utils.assertTwoAtomsAreEqual(boundCounterpartChain.select('calpha')[contactPositionInChain], elementcontact, useCoords=False, useResname=True)                       
							self.utils.assertTwoAtomsAreEqual(unboundCounterpartChain.select('calpha')[contactPositionInChain], elementcontact, useCoords=False, useResname=True)
						else:
							elementPositionInChain = encounter.accessANMs().getCalphaPosition(element, mobChain.select('calpha'))
							contactPositionInChain = encounter.accessANMs().getCalphaPosition(elementcontact, boundCounterpartChain.select('calpha'))
							r_ij = calcDistance(refchain.select('calpha')[elementPositionInChain], unboundCounterpartChain.select('calpha')[contactPositionInChain])                        
							r_ij_b = r_ij
							self.utils.assertTwoAtomsAreEqual(mobChain.select('calpha')[elementPositionInChain], element, useCoords=False, useResname=True)
							self.utils.assertTwoAtomsAreEqual(refchain.select('calpha')[elementPositionInChain], element, useCoords=False, useResname=True)
							self.utils.assertTwoAtomsAreEqual(boundCounterpartChain.select('calpha')[contactPositionInChain], elementcontact, useCoords=False, useResname=True)                       
							self.utils.assertTwoAtomsAreEqual(unboundCounterpartChain.select('calpha')[contactPositionInChain], elementcontact, useCoords=False, useResname=True)
						# make the 3*3 hessian term for this contact (excluding gamma, gamma is multiplied at the end to the sum)  
						deltaTerm = self.make3By3HessianTerm(refchain.select('calpha')[elementPositionInChain], unboundCounterpartChain.select('calpha')[contactPositionInChain], r_ij, r_ij_b)
						#print refchain.select('calpha')[elementPositionInChain], unboundCounterpartChain.select('calpha')[contactPositionInChain], " r_ij, rij_b: ", r_ij, r_ij_b
						overallTerm += deltaTerm
				# multiply the overallTerm with the spring constant gamma
				overallTerm = overallTerm * self.utils.config.customForceConstant
				# add the overallterm to the hessian matrix
				if neighborhoodFrom == "C_b":
					elementPosition = encounter.accessANMs().getCalphaPosition(refchain.select('calpha')[elementPositionInChain], reference.select('calpha'))
				else:
					elementPosition = encounter.accessANMs().getCalphaPosition(element, reference.select('calpha'))
				HR = self.add3By3MatrixtoHessian(overallTerm, HR, elementPosition*3)
		print "added custom terms to hessian"
		print "total intermolecular contacts: ", contactsCounter
		return HR, interCalphaIndices    
	
	def calcCustomH_ANew_IJ(self, referenceStructure, neighborStructure, encounter, areStructuresChains, equilibriumAt, workOnReceptor=True, selstr='calpha'):
		""" Creates the HRL matrix made through intramolecular contacts.
		
			Args:
				referenceStructure: structure to take calphas from, the hessian HR belongs to it or to its superset if I is a chain
				neighborStructure: structure to apply the neighborhood calculations on
				encounter: object with all encounter information
				areStructuresChains: boolean to describe if the structures are chains (subsets)
				equilibriumAt: is the equilibrium set to r_ij or r_ij_b
				workonReceptor: is the Hessian and the referenceStructure receptor or ligand
				selstr: atomType of the course grained ANM (by default calpha)
		"""
		assert equilibriumAt == "r_ij" or equilibriumAt == "r_ij_b"
		
		if workOnReceptor:
			if areStructuresChains:
				if self.bound_provided == True:				
					mobile = encounter.getMobChain()
					boundCounterpart = encounter.getBoundCounterpartChain()
				else:
					pass
			else:
				reference = encounter.getReference()
				unboundCounterpart = encounter.getUnboundCounterpart()
				if self.bound_provided == True:
					refchain = encounter.getRefChain()
					mobile = encounter.getMobile()
					mobChain = encounter.getMobChain()
					boundCounterpart = encounter.getBoundCounterpart()
					boundCounterpartChain = encounter.getBoundCounterpartChain()
					unboundCounterpartChain = encounter.getUnboundCounterpartChain()
		else:
			if areStructuresChains:
				if self.bound_provided == True:	
					mobile = encounter.getBoundCounterpartChain()
					boundCounterpart = encounter.getMobChain()
				else:
					pass
			else:
				reference = encounter.getUnboundCounterpart()
				unboundCounterpart = encounter.getReference()
				if self.bound_provided == True:
					refchain = encounter.getUnboundCounterpartChain()
					mobile = encounter.getBoundCounterpart()
					mobChain = encounter.getBoundCounterpartChain()
					boundCounterpart = encounter.getMobile()
					boundCounterpartChain = encounter.getMobChain()
					unboundCounterpartChain = encounter.getRefChain()
				
		neighborStructureCalpha = neighborStructure.select('calpha')
		offDiagonalHessianMatrix = np.zeros(((reference.select('calpha').numAtoms()*3), (unboundCounterpart.select('calpha').numAtoms()*3) ))
		contactsCounter = 0
				
		for idx, element in enumerate(referenceStructure.select('calpha')):
			contactsOfI = encounter.getIntermolecularNeighborsOfAtom(element, neighborStructure, selstr, str(self.utils.config.customHRdistance))
			# if element has contacts in the neighborStructure, the hessian needs an update in the 3*3 matrix on the diagonal of this element atom
			if contactsOfI:             
				print "intermolecular contacts: ", contactsOfI.numAtoms()
				contactsCounter += contactsOfI.numAtoms()
#               print "contact at i, refChainCalphas[i]: ", i, refChainCalphas[i]
				contacts_counterpartChainIndices = self.utils.getMatchingStructureSelections(neighborStructureCalpha, contactsOfI, neighborStructureCalpha)
				assert len(contactsOfI) == len(contacts_counterpartChainIndices)
				# access each element contact to create the deltaTerm
				for elementcontact, contacts_counterpartChainIndex in zip(contactsOfI, contacts_counterpartChainIndices):
					overallTerm = np.zeros((3,3))
					#self.utils.assertTwoAtomsAreEqual(refChainCalphas[i], mobChainCalphas[i], useCoords=False, useResname=True)
					#self.utils.assertTwoAtomsAreEqual(elementcontact, boundCounterpartChainCalphas[contacts_counterpartChainIndex], useCoords=False, useResname=True)
					r_ij = calcDistance(element, elementcontact)
					if equilibriumAt == "r_ij":
						r_ij_b = r_ij
					#if element is not in matched reference or contact is not in matched counterpart: r_ij_b = r_ij
					elif not(element in refchain.select('calpha')) or not(elementcontact in unboundCounterpartChain.select('calpha')):
						r_ij_b = r_ij
					else:
						elementPositionInChain = encounter.accessANMs().getCalphaPosition(element, refchain.select('calpha'))
						contactPositionInChain = encounter.accessANMs().getCalphaPosition(elementcontact, unboundCounterpartChain.select('calpha'))
						r_ij_b = calcDistance(mobChain.select('calpha')[elementPositionInChain], boundCounterpartChain.select('calpha')[contactPositionInChain])                        
						self.utils.assertTwoAtomsAreEqual(mobChain.select('calpha')[elementPositionInChain], element, useCoords=False, useResname=True)  
						self.utils.assertTwoAtomsAreEqual(refchain.select('calpha')[elementPositionInChain], element, useCoords=False, useResname=True)  
						self.utils.assertTwoAtomsAreEqual(boundCounterpartChain.select('calpha')[contactPositionInChain], elementcontact, useCoords=False, useResname=True)                       
						self.utils.assertTwoAtomsAreEqual(unboundCounterpartChain.select('calpha')[contactPositionInChain], elementcontact, useCoords=False, useResname=True)                       
#                     
					# make the 3*3 hessian term for this contact (excluding gamma, gamma is multiplied at the end to the sum)  
					deltaTerm = self.make3By3OffDiagonalHessianTermIJ(element, elementcontact, r_ij, r_ij_b)
					overallTerm += deltaTerm
					#print "r_ij, r_ij_b: ", r_ij, r_ij_b
					# multiply the overallTerm with the spring constant gamma
					overallTerm = overallTerm * self.utils.config.customForceConstant
#                   print overallTerm
					offDiagonalHessianMatrix = self.add3By3MatrixtoOffDiagonalHessianMatrixIJ(overallTerm, offDiagonalHessianMatrix, idx*3, contacts_counterpartChainIndex*3)
				#print contactsOfI.numAtoms(), "neighbors, modifying at hessian (loopcounter*3)+1: ", str((loopCounter*3)+1)
				#print str(i)+"'th refchain calpha, hessian line number ", (loopCounter*3)+1, "contacts with ", unboundCounterpartChainCalphas[contacts_counterpartChainIndex], " unboundcounterpartchainindex: ", contacts_counterpartChainIndices
				#print ""
				# add the overallterm to the hessian matrix
				###elementPosition = encounter.accessANMs().getCalphaPosition(element, encounter.getReference().select('calpha'))
		print "added custom terms to offDiagonalHessianMatrix"
		print "total intermolecular contacts: ", contactsCounter
		return offDiagonalHessianMatrix    
	
	def calcCustomH_ANew_U1(self, HR, referenceStructure, neighborStructure, encounter, areStructuresChains, equilibriumAt, workOnReceptor=True, selstr='calpha'):
		""" Modifies the hessian HR or HL by adding additonal terms for intramolecular contacts.
		
			Args:
				HR: The original HR as calculated by prody
				referenceStructure: structure to take calphas from, the hessian HR belongs to it or to its superset if I is a chain
				neighborStructure: structure to apply the neighborhood calculations on
				encounter: object with all encounter information
				areStructuresChains: boolean to describe if the structures are chains (subsets)
				equilibriumAt: is the equilibrium set to r_ij or r_ij_b
				workonReceptor: is the Hessian and the referenceStructure receptor or ligand
				selstr: atomType of the course grained ANM (by default calpha)
		"""
		assert equilibriumAt == "r_ij" or equilibriumAt == "r_ij_b"
		
		if workOnReceptor:
			refchain = encounter.getRefChain()
			mobile = encounter.getMobile()
			mobChain = encounter.getMobChain()
			boundCounterpart = encounter.getBoundCounterpart()
			boundCounterpartChain = encounter.getBoundCounterpartChain()
			unboundCounterpartChain = encounter.getUnboundCounterpartChain()
		else:  
			refchain = encounter.getUnboundCounterpartChain()
			mobile = encounter.getBoundCounterpart()
			mobChain = encounter.getBoundCounterpartChain()
			boundCounterpart = encounter.getMobile()
			boundCounterpartChain = encounter.getMobChain()
			unboundCounterpartChain = encounter.getRefChain()
				
		neighborStructureCalpha = neighborStructure.select('calpha')
				
		for idx, element in enumerate(referenceStructure.select('calpha')):
			contactsOfI = encounter.getIntermolecularNeighborsOfAtom(element, neighborStructure, selstr, str(self.utils.config.customHRdistance))
			# if element has contacts in the neighborStructure, the hessian needs an update in the 3*3 matrix on the diagonal of this element atom
			if contactsOfI:
#               print "contact at i, refChainCalphas[i]: ", i, refChainCalphas[i]
				contacts_counterpartChainIndices = self.utils.getMatchingStructureSelections(neighborStructureCalpha, contactsOfI, neighborStructureCalpha)
				assert len(contactsOfI) == len(contacts_counterpartChainIndices)
				# access each element contact to create the deltaTerm
				overallTerm = np.zeros((3,3))
				for elementcontact, contacts_counterpartChainIndex in zip(contactsOfI, contacts_counterpartChainIndices):
					#self.utils.assertTwoAtomsAreEqual(refChainCalphas[i], mobChainCalphas[i], useCoords=False, useResname=True)
					#self.utils.assertTwoAtomsAreEqual(elementcontact, boundCounterpartChainCalphas[contacts_counterpartChainIndex], useCoords=False, useResname=True)
					if equilibriumAt == "r_ij_b":
						r_ij_b = calcDistance(element, elementcontact)
						elementPositionInChain = encounter.accessANMs().getCalphaPosition(element, mobChain.select('calpha'))
						contactPositionInChain = encounter.accessANMs().getCalphaPosition(elementcontact, boundCounterpartChain.select('calpha'))
						r_ij = calcDistance(refchain.select('calpha')[elementPositionInChain], unboundCounterpartChain.select('calpha')[contactPositionInChain])                        
						self.utils.assertTwoAtomsAreEqual(mobChain.select('calpha')[elementPositionInChain], element, useCoords=False, useResname=True)
						self.utils.assertTwoAtomsAreEqual(refchain.select('calpha')[elementPositionInChain], element, useCoords=False, useResname=True)
						self.utils.assertTwoAtomsAreEqual(boundCounterpartChain.select('calpha')[contactPositionInChain], elementcontact, useCoords=False, useResname=True)                       
						self.utils.assertTwoAtomsAreEqual(unboundCounterpartChain.select('calpha')[contactPositionInChain], elementcontact, useCoords=False, useResname=True)
					else:
						elementPositionInChain = encounter.accessANMs().getCalphaPosition(element, mobChain.select('calpha'))
						contactPositionInChain = encounter.accessANMs().getCalphaPosition(elementcontact, boundCounterpartChain.select('calpha'))
						r_ij = calcDistance(refchain.select('calpha')[elementPositionInChain], unboundCounterpartChain.select('calpha')[contactPositionInChain])                        
						r_ij_b = r_ij
						self.utils.assertTwoAtomsAreEqual(mobChain.select('calpha')[elementPositionInChain], element, useCoords=False, useResname=True)
						self.utils.assertTwoAtomsAreEqual(refchain.select('calpha')[elementPositionInChain], element, useCoords=False, useResname=True)
						self.utils.assertTwoAtomsAreEqual(boundCounterpartChain.select('calpha')[contactPositionInChain], elementcontact, useCoords=False, useResname=True)                       
						self.utils.assertTwoAtomsAreEqual(unboundCounterpartChain.select('calpha')[contactPositionInChain], elementcontact, useCoords=False, useResname=True)
						#r_ij_b = calcDistance(zip(mobile.select('calpha'))[idx][0], zip(boundCounterpart.select('calpha'))[contacts_counterpartChainIndex][0])
					# make the 3*3 hessian term for this contact (excluding gamma, gamma is multiplied at the end to the sum)  
					deltaTerm = self.make3By3HessianTerm(refchain.select('calpha')[elementPositionInChain], unboundCounterpartChain.select('calpha')[contactPositionInChain], r_ij, r_ij_b)
					print refchain.select('calpha')[elementPositionInChain], unboundCounterpartChain.select('calpha')[contactPositionInChain], " r_ij, rij_b: ", r_ij, r_ij_b
					overallTerm += deltaTerm
					#print "r_ij, r_ij_b: ", r_ij, r_ij_b
				# multiply the overallTerm with the spring constant gamma
				overallTerm = overallTerm * self.utils.config.customForceConstant
#                     print overallTerm
				#print contactsOfI.numAtoms(), "neighbors, modifying at hessian (loopcounter*3)+1: ", str((loopCounter*3)+1)
				#print str(i)+"'th refchain calpha, hessian line number ", (loopCounter*3)+1, "contacts with ", unboundCounterpartChainCalphas[contacts_counterpartChainIndex], " unboundcounterpartchainindex: ", contacts_counterpartChainIndices
				#print ""
				# add the overallterm to the hessian matrix
				elementPosition = encounter.accessANMs().getCalphaPosition(refchain.select('calpha')[elementPositionInChain], encounter.getReference().select('calpha'))
				HR = self.add3By3MatrixtoHessian(overallTerm, HR, elementPosition*3)
				print "adding to hessian at: ", (elementPosition*3+1)
		print "added custom terms to hessian"
		return HR    
	
	def calcCustomH_A(self, HR, encounter, workOnReceptor=True, selstr='calpha'):
		""" Modifies the hessian of anm_reference according to calcCustomH_A and returns it. """
		if workOnReceptor:
			refChainCalphas = encounter.getRefChain().select('calpha')
			mobChainCalphas = encounter.getMobChain().select('calpha')
			mobChain = encounter.getMobChain()
			refChain = encounter.getRefChain()
			boundCounterpartChainCalphas = encounter.getBoundCounterpartChain().select('calpha')
			boundCounterpartChain = encounter.getBoundCounterpartChain()
			unboundCounterpartChain = encounter.getUnboundCounterpartChain()
			unboundCounterpartChainCalphas = encounter.getUnboundCounterpartChain().select('calpha')    
			referenceCalphas = encounter.getReference().select('calpha')
		else:
			refChainCalphas = encounter.getUnboundCounterpartChain().select('calpha')
			mobChainCalphas = encounter.getBoundCounterpartChain().select('calpha')
			mobChain = encounter.getBoundCounterpartChain()
			refChain = encounter.getUnboundCounterpartChain()
			boundCounterpartChainCalphas = encounter.getMobChain().select('calpha')
			boundCounterpartChain = encounter.getMobChain()
			unboundCounterpartChain = encounter.getRefChain()
			unboundCounterpartChainCalphas = encounter.getRefChain().select('calpha')    
			referenceCalphas = encounter.getUnboundCounterpart().select('calpha')  
		
		#encounter.printIntermolecularNeighbors(encounter.getReference(), encounter.getUnboundCounterpart(), selstr, str(self.utils.config.customHRdistance))

		# Loop over all calphas in the reference structure (using matched chains) 
		counterUnmatchedCalphas = 0
		loopCounter = 0
		for element in referenceCalphas:
			i = loopCounter - counterUnmatchedCalphas
			if self.utils.doesAtomExistInY(element, refChainCalphas) is None:
				counterUnmatchedCalphas += 1
				loopCounter += 1
				continue
			else:
				contactsOfI = encounter.getIntermolecularNeighbors(refChain, unboundCounterpartChain, i, selstr, str(self.utils.config.customHRdistance))
				# if there are contacts in the unbound counterpart, the hessian needs an update in the 3*3 matrix of the diagonal of this atom
				if contactsOfI:
#                     print "contact at i, refChainCalphas[i]: ", i, refChainCalphas[i]
					contacts_counterpartChainIndices = self.utils.getMatchingStructure(unboundCounterpartChainCalphas, contactsOfI, boundCounterpartChainCalphas)
					assert len(contactsOfI) == len(contacts_counterpartChainIndices)
					# access each element contact to create the deltaTerm
					overallTerm = np.zeros((3,3))
					for elementcontact, contacts_counterpartChainIndex in zip(contactsOfI, contacts_counterpartChainIndices):
						self.utils.assertTwoAtomsAreEqual(refChainCalphas[i], mobChainCalphas[i], useCoords=False, useResname=True)
						self.utils.assertTwoAtomsAreEqual(elementcontact, boundCounterpartChainCalphas[contacts_counterpartChainIndex], useCoords=False, useResname=True)
						r_ij = calcDistance(refChainCalphas[i], elementcontact)
						r_ij_b = calcDistance(mobChainCalphas[i], boundCounterpartChainCalphas[contacts_counterpartChainIndex])
						# make the 3*3 hessian term for this contact (excluding gamma, gamma is multiplied at the end to the sum)  
						deltaTerm = self.make3By3HessianTerm(refChainCalphas[i], elementcontact, r_ij, r_ij_b)
						overallTerm += deltaTerm
						#print "r_ij, r_ij_b: ", r_ij, r_ij_b
					# multiply the overallTerm with the spring constant gamma
					overallTerm = overallTerm * self.utils.config.customForceConstant
#                     print overallTerm
					#print contactsOfI.numAtoms(), "neighbors, modifying at hessian (loopcounter*3)+1: ", str((loopCounter*3)+1)
					print str(i)+"'th refchain calpha, hessian line number ", (loopCounter*3)+1, "contacts with ", unboundCounterpartChainCalphas[contacts_counterpartChainIndex], " unboundcounterpartchainindex: ", contacts_counterpartChainIndices
					print ""
					# add the overallterm to the hessian matrix
					HR = self.add3By3MatrixtoHessian(overallTerm, HR, loopCounter*3)
				loopCounter += 1
		assert(loopCounter-counterUnmatchedCalphas) == refChainCalphas.numAtoms()
		print "added custom terms to hessian"
		return HR
	
	def calcCustomH_A_IJ(self, encounter, workOnReceptor=True, selstr='calpha'):
		""" Modifies the hessian of anm_reference according to calcCustomH_A and returns it. """
		if workOnReceptor:
			refChainCalphas = encounter.getRefChain().select('calpha')
			mobChainCalphas = encounter.getMobChain().select('calpha')
			mobChain = encounter.getMobChain()
			refChain = encounter.getRefChain()
			boundCounterpartChainCalphas = encounter.getBoundCounterpartChain().select('calpha')
			boundCounterpartChain = encounter.getBoundCounterpartChain()
			unboundCounterpartChain = encounter.getUnboundCounterpartChain()
			unboundCounterpartChainCalphas = encounter.getUnboundCounterpartChain().select('calpha')    
			referenceCalphas = encounter.getReference().select('calpha')
			mobileCalphas = encounter.getMobile().select('calpha')
			unboundCounterpart = encounter.getUnboundCounterpart()
			unboundCounterpartCalphas = encounter.getUnboundCounterpart().select('calpha')
		else:
			refChainCalphas = encounter.getUnboundCounterpartChain().select('calpha')
			mobChainCalphas = encounter.getBoundCounterpartChain().select('calpha')
			mobChain = encounter.getBoundCounterpartChain()
			refChain = encounter.getUnboundCounterpartChain()
			boundCounterpartChainCalphas = encounter.getMobChain().select('calpha')
			boundCounterpartChain = encounter.getMobChain()
			unboundCounterpartChain = encounter.getRefChain()
			unboundCounterpartChainCalphas = encounter.getRefChain().select('calpha')    
			referenceCalphas = encounter.getUnboundCounterpart().select('calpha')  
			mobileCalphas = encounter.getBoundCounterpart().select('calpha')
			unboundCounterpart = encounter.getReference()
			unboundCounterpartCalphas = encounter.getReference().select('calpha')
		
		offDiagonalHessianMatrix = np.zeros(((referenceCalphas.numAtoms()*3), (unboundCounterpartCalphas.numAtoms()*3) ))
		
		#encounter.printIntermolecularNeighbors(encounter.getReference(), encounter.getUnboundCounterpart(), selstr, str(self.utils.config.customHRdistance))

		# Loop over all calphas in the reference structure (using matched chains) 
		counterUnmatchedCalphas = 0
		loopCounter = 0
		for element in referenceCalphas:
			i = loopCounter - counterUnmatchedCalphas
			if self.utils.doesAtomExistInY(element, refChainCalphas) is None:
				counterUnmatchedCalphas += 1
				loopCounter += 1
				continue
			else:
				contactsOfI = encounter.getIntermolecularNeighbors(refChain, unboundCounterpartChain, i, selstr, str(self.utils.config.customHRdistance))
				# if there are contacts in the unbound counterpart, the hessian needs an update in the 3*3 matrix of the diagonal of this atom
				if contactsOfI:
#                     print "contact at i, refChainCalphas[i]: ", i, refChainCalphas[i]
					contacts_counterpartChainIndices = self.utils.getMatchingStructure(unboundCounterpartChainCalphas, contactsOfI, boundCounterpartChainCalphas)
					assert len(contactsOfI) == len(contacts_counterpartChainIndices)
					# access each element contact to create the deltaTerm
					for elementcontact, contacts_counterpartChainIndex in zip(contactsOfI, contacts_counterpartChainIndices):
						overallTerm = np.zeros((3,3))
						self.utils.assertTwoAtomsAreEqual(refChainCalphas[i], mobChainCalphas[i], useCoords=False, useResname=True)
						self.utils.assertTwoAtomsAreEqual(elementcontact, boundCounterpartChainCalphas[contacts_counterpartChainIndex], useCoords=False, useResname=True)
						r_ij = calcDistance(refChainCalphas[i], elementcontact)
						r_ij_b = calcDistance(mobChainCalphas[i], boundCounterpartChainCalphas[contacts_counterpartChainIndex])
						# make the 3*3 hessian term for this contact (excluding gamma, gamma is multiplied at the end to the sum)  
						deltaTerm = self.make3By3OffDiagonalHessianTermIJ(refChainCalphas[i], unboundCounterpartChainCalphas[contacts_counterpartChainIndex], r_ij, r_ij)
						overallTerm += deltaTerm
						#print "r_ij, r_ij_b: ", r_ij, r_ij_b
						# multiply the overallTerm with the spring constant gamma
						overallTerm = overallTerm * self.utils.config.customForceConstant
						counterPartCalphaPosition = encounter.accessANMs().getCalphaPosition(unboundCounterpartChainCalphas[contacts_counterpartChainIndex], unboundCounterpart)
						print "off diagonal i,j "+str(loopCounter*3)+" "+str(counterPartCalphaPosition*3)+ " term: ", overallTerm                        
						offDiagonalHessianMatrix = self.add3By3MatrixtoOffDiagonalHessianMatrixIJ(overallTerm, offDiagonalHessianMatrix, loopCounter*3, counterPartCalphaPosition*3)
					#print contactsOfI.numAtoms(), "neighbors, modifying at hessian (loopcounter*3)+1: ", str((loopCounter*3)+1)
					print str(i)+"'th refchain calpha, hessian line number ", (loopCounter*3)+1, "contacts with ", unboundCounterpartChainCalphas[contacts_counterpartChainIndex], " unboundcounterpartchainindex: ", contacts_counterpartChainIndices
					print ""
				loopCounter += 1
		assert(loopCounter-counterUnmatchedCalphas) == refChainCalphas.numAtoms()
		print "added custom terms to hessian"
		return offDiagonalHessianMatrix
	
	def calcCustomH_A_NeighborsBound(self, HR, encounter, selstr='calpha'):
		""" Modifies the hessian of anm_reference according to calcCustomH_A and returns it. """
		refChainCalphas = encounter.getRefChain().select('calpha')
		mobChainCalphas = encounter.getMobChain().select('calpha')
		boundCounterpartChainCalphas = encounter.getBoundCounterpartChain().select('calpha')
		unboundCounterpartChainCalphas = encounter.getUnboundCounterpartChain().select('calpha')    
		
		referenceCalphas = encounter.getReference().select('calpha')
		mobileCalphas = encounter.getMobile().select('calpha')
		
		#encounter.printIntermolecularNeighbors(encounter.getMobile(), encounter.getBoundCounterpart(), selstr, str(self.utils.config.customHRdistance))

		# Loop over all calphas in the reference structure (using matched chains) 
		counterUnmatchedCalphas = 0
		loopCounter = 0
		for element in referenceCalphas:
			i = loopCounter - counterUnmatchedCalphas
			if self.utils.doesAtomExistInY(element, refChainCalphas) is None:
				counterUnmatchedCalphas += 1
				loopCounter += 1
				continue
			else:
				contactsOfI = encounter.getIntermolecularNeighbors(encounter.getMobChain(), encounter.getBoundCounterpartChain(), i, selstr, str(self.utils.config.customHRdistance))
				# if there are contacts in the unbound counterpart, the hessian needs an update in the 3*3 matrix of the diagonal of this atom
				if contactsOfI:
#                     print "contact at i, refChainCalphas[i]: ", i, refChainCalphas[i]
					contacts_counterpartChainIndices = self.utils.getMatchingStructure(boundCounterpartChainCalphas, contactsOfI, unboundCounterpartChainCalphas)
					assert len(contactsOfI) == len(contacts_counterpartChainIndices)
					# access each element contact to create the deltaTerm
					overallTerm = np.zeros((3,3))
					for elementcontact, contacts_counterpartChainIndex in zip(contactsOfI, contacts_counterpartChainIndices):
						self.utils.assertTwoAtomsAreEqual(refChainCalphas[i], mobChainCalphas[i], useCoords=False, useResname=False)
						self.utils.assertTwoAtomsAreEqual(elementcontact, unboundCounterpartChainCalphas[contacts_counterpartChainIndex], useCoords=False, useResname=False)
						self.utils.assertTwoAtomsAreEqual(boundCounterpartChainCalphas[contacts_counterpartChainIndex], elementcontact, useCoords=False, useResname=False)
						self.utils.assertTwoAtomsAreEqual(boundCounterpartChainCalphas[contacts_counterpartChainIndex], unboundCounterpartChainCalphas[contacts_counterpartChainIndex], useCoords=False, useResname=False)
						r_ij = calcDistance(refChainCalphas[i], unboundCounterpartChainCalphas[contacts_counterpartChainIndex])
						r_ij_b = calcDistance(mobChainCalphas[i], elementcontact)
						# make the 3*3 hessian term for this contact (excluding gamma, gamma is multiplied at the end to the sum)  
						# if customHR_B, just use the distance d_0, else use the true distance in the bound pairs for the second derivatives
						if self.utils.config.customHR_B:
							if r_ij >= self.utils.config.customHRdistance:
								deltaTerm = self.make3By3HessianTerm(refChainCalphas[i], unboundCounterpartChainCalphas[contacts_counterpartChainIndex], r_ij, self.utils.config.customHRdistance)
								overallTerm += deltaTerm
						else:
							deltaTerm = self.make3By3HessianTerm(refChainCalphas[i], unboundCounterpartChainCalphas[contacts_counterpartChainIndex], r_ij, r_ij_b)
							overallTerm += deltaTerm
						#print "r_ij, r_ij_b: ", r_ij, r_ij_b
					# multiply the overallTerm with the spring constant gamma
					overallTerm = overallTerm * self.utils.config.customForceConstant
					#print overallTerm
					print contactsOfI.numAtoms(), "neighbors, modifying at hessian (loopcounter*3)+1: ", str((loopCounter*3)+1)
					#print contactsOfI.getSelstr()
					print str(i)+"'th refchain calpha, hessian line number ", (loopCounter*3)+1, "contacts with ", unboundCounterpartChainCalphas[contacts_counterpartChainIndex], " unboundcounterpartchainindex: ", contacts_counterpartChainIndices
					print ""
					# add the overallterm to the hessian matrix
					HR = self.add3By3MatrixtoHessian(overallTerm, HR, loopCounter*3)
				loopCounter += 1
		assert(loopCounter-counterUnmatchedCalphas) == refChainCalphas.numAtoms()
		print "added custom terms to hessian"
		return HR
	
	def calcCustomH_A_NeighborsBoundGeneral(self, HR, encounter, workOnReceptor=True, selstr='calpha'):
		""" Modifies the hessian of anm_reference according to calcCustomH_A and returns it. """
		if workOnReceptor:
			refChainCalphas = encounter.getRefChain().select('calpha')
			mobChainCalphas = encounter.getMobChain().select('calpha')
			mobChain = encounter.getMobChain()
			boundCounterpartChainCalphas = encounter.getBoundCounterpartChain().select('calpha')
			boundCounterpartChain = encounter.getBoundCounterpartChain()
			unboundCounterpartChainCalphas = encounter.getUnboundCounterpartChain().select('calpha')    
			referenceCalphas = encounter.getReference().select('calpha')
		else:
			refChainCalphas = encounter.getUnboundCounterpartChain().select('calpha')
			mobChainCalphas = encounter.getBoundCounterpartChain().select('calpha')
			mobChain = encounter.getBoundCounterpartChain()
			boundCounterpartChainCalphas = encounter.getMobChain().select('calpha')
			boundCounterpartChain = encounter.getMobChain()
			unboundCounterpartChainCalphas = encounter.getRefChain().select('calpha')    
			referenceCalphas = encounter.getUnboundCounterpart().select('calpha')            
		
		#encounter.printIntermolecularNeighbors(encounter.getMobile(), encounter.getBoundCounterpart(), selstr, str(self.utils.config.customHRdistance))

		# Loop over all calphas in the reference structure (using matched chains) 
		counterUnmatchedCalphas = 0
		loopCounter = 0
		for element in referenceCalphas:
			i = loopCounter - counterUnmatchedCalphas
			if self.utils.doesAtomExistInY(element, refChainCalphas) is None:
				counterUnmatchedCalphas += 1
				loopCounter += 1
				continue
			else:
				contactsOfI = encounter.getIntermolecularNeighbors(mobChain, boundCounterpartChain, i, selstr, str(self.utils.config.customHRdistance))
				# if there are contacts in the unbound counterpart, the hessian needs an update in the 3*3 matrix of the diagonal of this atom
				if contactsOfI:
#                     print "contact at i, refChainCalphas[i]: ", i, refChainCalphas[i]
					contacts_counterpartChainIndices = self.utils.getMatchingStructure(boundCounterpartChainCalphas, contactsOfI, unboundCounterpartChainCalphas)
					assert len(contactsOfI) == len(contacts_counterpartChainIndices)
					# access each element contact to create the deltaTerm
					overallTerm = np.zeros((3,3))
					for elementcontact, contacts_counterpartChainIndex in zip(contactsOfI, contacts_counterpartChainIndices):
						self.utils.assertTwoAtomsAreEqual(refChainCalphas[i], mobChainCalphas[i], useCoords=False, useResname=False)
						self.utils.assertTwoAtomsAreEqual(elementcontact, unboundCounterpartChainCalphas[contacts_counterpartChainIndex], useCoords=False, useResname=False)
						self.utils.assertTwoAtomsAreEqual(boundCounterpartChainCalphas[contacts_counterpartChainIndex], elementcontact, useCoords=False, useResname=False)
						self.utils.assertTwoAtomsAreEqual(boundCounterpartChainCalphas[contacts_counterpartChainIndex], unboundCounterpartChainCalphas[contacts_counterpartChainIndex], useCoords=False, useResname=False)
						r_ij = calcDistance(refChainCalphas[i], unboundCounterpartChainCalphas[contacts_counterpartChainIndex])
						r_ij_b = calcDistance(mobChainCalphas[i], elementcontact)
						# make the 3*3 hessian term for this contact (excluding gamma, gamma is multiplied at the end to the sum)  
						# if customHR_B, just use the distance d_0, else use the true distance in the bound pairs for the second derivatives
						if self.utils.config.customHR_B:
							if r_ij >= self.utils.config.customHRdistance:
								deltaTerm = self.make3By3HessianTerm(refChainCalphas[i], unboundCounterpartChainCalphas[contacts_counterpartChainIndex], r_ij, self.utils.config.customHRdistance)
								overallTerm += deltaTerm
						else:
							deltaTerm = self.make3By3HessianTerm(refChainCalphas[i], unboundCounterpartChainCalphas[contacts_counterpartChainIndex], r_ij, r_ij_b)
							overallTerm += deltaTerm
						#print "r_ij, r_ij_b: ", r_ij, r_ij_b
					# multiply the overallTerm with the spring constant gamma
					overallTerm = overallTerm * self.utils.config.customForceConstant
					#print overallTerm
					print contactsOfI.numAtoms(), "neighbors, modifying at hessian (loopcounter*3)+1: ", str((loopCounter*3)+1)
					#print contactsOfI.getSelstr()
					print str(i)+"'th refchain calpha, hessian line number ", (loopCounter*3)+1, "contacts with ", unboundCounterpartChainCalphas[contacts_counterpartChainIndex], " unboundcounterpartchainindex: ", contacts_counterpartChainIndices
					print ""
					# add the overallterm to the hessian matrix
					HR = self.add3By3MatrixtoHessian(overallTerm, HR, loopCounter*3)
				loopCounter += 1
		assert(loopCounter-counterUnmatchedCalphas) == refChainCalphas.numAtoms()
		print "added custom terms to hessian"
		return HR
   
	def calcOffDiagonalHessianBlockMatrixGeneral_IJ(self, encounter, workOnReceptor=True, selstr='calpha'):
		""" Creates the off diagonal hessian block matrix and returns it. """
		if workOnReceptor:
			refChainCalphas = encounter.getRefChain().select('calpha')
			mobChainCalphas = encounter.getMobChain().select('calpha')
			mobChain = encounter.getMobChain()
			boundCounterpartChainCalphas = encounter.getBoundCounterpartChain().select('calpha')
			boundCounterpartChain = encounter.getBoundCounterpartChain()
			unboundCounterpartChainCalphas = encounter.getUnboundCounterpartChain().select('calpha')    
			referenceCalphas = encounter.getReference().select('calpha')
			mobileCalphas = encounter.getMobile().select('calpha')
			unboundCounterpart = encounter.getUnboundCounterpart()   
			unboundCounterpartCalphas = encounter.getUnboundCounterpart().select('calpha')
		else:
			refChainCalphas = encounter.getUnboundCounterpartChain().select('calpha')
			mobChainCalphas = encounter.getBoundCounterpartChain().select('calpha')
			mobChain = encounter.getBoundCounterpartChain()
			boundCounterpartChainCalphas = encounter.getMobChain().select('calpha')
			boundCounterpartChain = encounter.getMobChain()
			unboundCounterpartChainCalphas = encounter.getRefChain().select('calpha')    
			referenceCalphas = encounter.getUnboundCounterpart().select('calpha')  
			mobileCalphas = encounter.getBoundCounterpart().select('calpha')          
			unboundCounterpart = encounter.getReference() 
			unboundCounterpartCalphas = encounter.getReference().select('calpha')

		offDiagonalHessianMatrix = np.zeros(((referenceCalphas.numAtoms()*3), (unboundCounterpartCalphas.numAtoms()*3) ))

		# Loop over all calphas in the reference structure (using matched chains) 
		counterUnmatchedCalphas = 0
		loopCounter = 0
		for element in referenceCalphas:
			i = loopCounter - counterUnmatchedCalphas
			if self.utils.doesAtomExistInY(element, refChainCalphas) is None:
				counterUnmatchedCalphas += 1
				loopCounter += 1
				continue
			else:
				contactsOfI = encounter.getIntermolecularNeighbors(mobChain, boundCounterpartChain, i, selstr, str(self.utils.config.customHRdistance))
				# if there are contacts in the bound counterpart, the off diagonal part of the hessian needs an update in the 3*3 matrix of this atom and its neighbor
				if contactsOfI:
#                     print "contact at i, refChainCalphas[i]: ", i, refChainCalphas[i]
					contacts_counterpartChainIndices = self.utils.getMatchingStructure(boundCounterpartChainCalphas, contactsOfI, unboundCounterpartChainCalphas)
					assert len(contactsOfI) == len(contacts_counterpartChainIndices)
					# access each element contact to create the deltaTerm
					for elementcontact, contacts_counterpartChainIndex in zip(contactsOfI, contacts_counterpartChainIndices):
						overallTerm = np.zeros((3,3))
						self.utils.assertTwoAtomsAreEqual(refChainCalphas[i], mobChainCalphas[i], useCoords=False, useResname=False)
						self.utils.assertTwoAtomsAreEqual(elementcontact, unboundCounterpartChainCalphas[contacts_counterpartChainIndex], useCoords=False, useResname=False)
						self.utils.assertTwoAtomsAreEqual(boundCounterpartChainCalphas[contacts_counterpartChainIndex], elementcontact, useCoords=False, useResname=False)
						self.utils.assertTwoAtomsAreEqual(boundCounterpartChainCalphas[contacts_counterpartChainIndex], unboundCounterpartChainCalphas[contacts_counterpartChainIndex], useCoords=False, useResname=False)
						r_ij = calcDistance(refChainCalphas[i], unboundCounterpartChainCalphas[contacts_counterpartChainIndex])
						r_ij_b = calcDistance(mobChainCalphas[i], elementcontact)
						# make the 3*3 hessian term for this contact
						# if customHR_B, just use the distance d_0, else use the true distance in the bound pairs for the second derivatives
						if self.utils.config.customHR_B:
							if r_ij >= self.utils.config.customHRdistance:
								deltaTerm = self.make3By3OffDiagonalHessianTermIJ(refChainCalphas[i], unboundCounterpartChainCalphas[contacts_counterpartChainIndex], r_ij, self.utils.config.customHRdistance)
								overallTerm += deltaTerm
						else:
							deltaTerm = self.make3By3OffDiagonalHessianTermIJ(refChainCalphas[i], unboundCounterpartChainCalphas[contacts_counterpartChainIndex], r_ij, r_ij_b)
							overallTerm += deltaTerm
						# multiply the overallTerm with the spring constant gamma
						overallTerm = overallTerm * self.utils.config.customForceConstant
						# add the overall Term to the correct off diagonal super element in the hessian
						counterPartCalphaPosition = encounter.accessANMs().getCalphaPosition(unboundCounterpartChainCalphas[contacts_counterpartChainIndex], unboundCounterpart)
						offDiagonalHessianMatrix = self.add3By3MatrixtoOffDiagonalHessianMatrixIJ(overallTerm, offDiagonalHessianMatrix, loopCounter*3, counterPartCalphaPosition*3)
						
						#print "r_ij, r_ij_b: ", r_ij, r_ij_b
						#print overallTerm
						print contactsOfI.numAtoms(), "neighbors, modifying at hessian (loopcounter*3)+1: ", str((loopCounter*3)+1)
					#print contactsOfI.getSelstr()
					#print str(i)+"'th refchain calpha, hessian line number ", (loopCounter*3)+1, "contacts with ", unboundCounterpartChainCalphas[contacts_counterpartChainIndex], " unboundcounterpartchainindex: ", contacts_counterpartChainIndices
					print ""
				loopCounter += 1
		assert(loopCounter-counterUnmatchedCalphas) == refChainCalphas.numAtoms()
		print "added custom terms to hessian"
		return offDiagonalHessianMatrix
   
# origs     
	def secondDerivativeTermOnDiagonal(self, x_i, x_j, r_ij, r_ij_b):
		""" @V / @x_i@x_i (excluding gamma)"""
		result = 1 + (r_ij_b * np.power(x_j - x_i, 2) ) / np.power(r_ij, 3) - r_ij_b/r_ij
		return result
		 
	def secondDerivateTermOffDiagonal(self, x_i, x_j, y_i, y_j, r_ij, r_ij_b):
		"""  @V / @x_i@y_j (excluding gamma) """
		result = r_ij_b * (x_j - x_i) * ((y_j - y_i)/np.power(r_ij, 3))
		return result
		 
	def secondDerivateTermOffDiagonalAtomsIJ(self, x_i, x_j, y_i, y_j, r_ij, r_ij_b):
		""" Equation 21 before reducing, Atilgan paper,  @V / @x_i@y_j (excluding gamma) """
		result = -1.0 * r_ij_b * (x_j - x_i) * ((y_j - y_i)/np.power(r_ij, 3))
		return result
#
# using r_ij_b   
#     def secondDerivativeTermOnDiagonal(self, x_i, x_j, r_ij, r_ij_b):
#         """ @V / @x_i@x_i (excluding gamma) from paper, assume r_ij is at equilibrium r_ij_b.  """
#         result = np.power(x_j - x_i, 2) / np.power(r_ij_b, 2)
#         return result
#            
#     def secondDerivateTermOffDiagonal(self, x_i, x_j, y_i, y_j, r_ij, r_ij_b):
#         """  @V / @x_i@y_j (excluding gamma) from paper, assume r_ij is at equilibrium r_ij_b. """
#         result = ((x_j - x_i)*(y_j - y_i))/ np.power(r_ij_b, 2)
#         return result

# using r_ij   
#     def secondDerivativeTermOnDiagonal(self, x_i, x_j, r_ij, r_ij_b):
#         """ @V / @x_i@x_i (excluding gamma) from paper, assume r_ij is at equilibrium r_ij_b.  """
#         result = np.power(x_j - x_i, 2) / np.power(r_ij, 2)
#         return result
#            
#     def secondDerivateTermOffDiagonal(self, x_i, x_j, y_i, y_j, r_ij, r_ij_b):
#         """  @V / @x_i@y_j (excluding gamma) from paper, assume r_ij is at equilibrium r_ij_b. """
#         result = ((x_j - x_i)*(y_j - y_i))/ np.power(r_ij, 2)
#         return result      
		
	def make3By3HessianTerm(self, refChainCalpha, elementcontact, r_ij, r_ij_b):
		""" Create a 3 by 3 matrix with the added terms for the hessian diagnonal (excluding multiplication with gamma)"""
		x_i = refChainCalpha.getCoords()[0]
		y_i = refChainCalpha.getCoords()[1]
		z_i = refChainCalpha.getCoords()[2]
		x_j = elementcontact.getCoords()[0]
		y_j = elementcontact.getCoords()[1]
		z_j = elementcontact.getCoords()[2]
		deltaTerm = np.zeros((3,3))
		deltaTerm[0][0] = self.secondDerivativeTermOnDiagonal(x_i, x_j, r_ij, r_ij_b)
		deltaTerm[0][1] = self.secondDerivateTermOffDiagonal(x_i, x_j, y_i, y_j, r_ij, r_ij_b)
		deltaTerm[0][2] = self.secondDerivateTermOffDiagonal(x_i, x_j, z_i, z_j, r_ij, r_ij_b)
		deltaTerm[1][0] = deltaTerm[0][1]
		deltaTerm[1][1] = self.secondDerivativeTermOnDiagonal(y_i, y_j, r_ij, r_ij_b)
		deltaTerm[1][2] = self.secondDerivateTermOffDiagonal(y_i, y_j, z_i, z_j, r_ij, r_ij_b)
		deltaTerm[2][0] = deltaTerm[0][2]
		deltaTerm[2][1] = deltaTerm[1][2]
		deltaTerm[2][2] = self.secondDerivativeTermOnDiagonal(z_i, z_j, r_ij, r_ij_b)
		return deltaTerm
	
	def add3By3MatrixtoHessian(self, delta3by3, HR, topleftIndex):
		""" Add the delta3by3 matrix to its corresponding position of HR, located by 
		the topleftIndex. """
		HR[topleftIndex][topleftIndex] += delta3by3[0][0]
		HR[topleftIndex][topleftIndex+1] += delta3by3[0][1]
		HR[topleftIndex][topleftIndex+2] += delta3by3[0][2]
		HR[topleftIndex+1][topleftIndex] += delta3by3[1][0]
		HR[topleftIndex+1][topleftIndex+1] += delta3by3[1][1]
		HR[topleftIndex+1][topleftIndex+2] += delta3by3[1][2]
		HR[topleftIndex+2][topleftIndex] += delta3by3[2][0]
		HR[topleftIndex+2][topleftIndex+1] += delta3by3[2][1]
		HR[topleftIndex+2][topleftIndex+2] += delta3by3[2][2]
		return HR
	
	def add3By3MatrixtoOffDiagonalHessianMatrixIJ(self, delta3by3, offDiagonalHessianMatrix, topleftIndex, counterpartTopleftIndex):
		""" Add the delta3by3 matrix to its corresponding position of HR, located by 
		the topleftIndex. """
		offDiagonalHessianMatrix[topleftIndex][counterpartTopleftIndex] += delta3by3[0][0]
		offDiagonalHessianMatrix[topleftIndex][counterpartTopleftIndex+1] += delta3by3[0][1]
		offDiagonalHessianMatrix[topleftIndex][counterpartTopleftIndex+2] += delta3by3[0][2]
		offDiagonalHessianMatrix[topleftIndex+1][counterpartTopleftIndex] += delta3by3[1][0]
		offDiagonalHessianMatrix[topleftIndex+1][counterpartTopleftIndex+1] += delta3by3[1][1]
		offDiagonalHessianMatrix[topleftIndex+1][counterpartTopleftIndex+2] += delta3by3[1][2]
		offDiagonalHessianMatrix[topleftIndex+2][counterpartTopleftIndex] += delta3by3[2][0]
		offDiagonalHessianMatrix[topleftIndex+2][counterpartTopleftIndex+1] += delta3by3[2][1]
		offDiagonalHessianMatrix[topleftIndex+2][counterpartTopleftIndex+2] += delta3by3[2][2]
		return offDiagonalHessianMatrix
	
	def make3By3OffDiagonalHessianTermIJ(self, refChainCalpha, elementcontact, r_ij, r_ij_b):
		""" Create a 3 by 3 matrix with the added terms for the hessian super element off the diagnonal (excluding multiplication with gamma). """
		x_i = refChainCalpha.getCoords()[0]
		y_i = refChainCalpha.getCoords()[1]
		z_i = refChainCalpha.getCoords()[2]
		x_j = elementcontact.getCoords()[0]
		y_j = elementcontact.getCoords()[1]
		z_j = elementcontact.getCoords()[2]
		deltaTerm = np.zeros((3,3))
		deltaTerm[0][0] = self.secondDerivateTermOffDiagonalAtomsIJ(x_i, x_j, x_i, x_j, r_ij, r_ij_b)
		deltaTerm[0][1] = self.secondDerivateTermOffDiagonalAtomsIJ(x_i, x_j, y_i, y_j, r_ij, r_ij_b)
		deltaTerm[0][2] = self.secondDerivateTermOffDiagonalAtomsIJ(x_i, x_j, z_i, z_j, r_ij, r_ij_b)
		deltaTerm[1][0] = deltaTerm[0][1]
		deltaTerm[1][1] = self.secondDerivateTermOffDiagonalAtomsIJ(y_i, y_j, y_i, y_j, r_ij, r_ij_b)
		deltaTerm[1][2] = self.secondDerivateTermOffDiagonalAtomsIJ(y_i, y_j, z_i, z_j, r_ij, r_ij_b)
		deltaTerm[2][0] = deltaTerm[0][2]
		deltaTerm[2][1] = deltaTerm[1][2]
		deltaTerm[2][2] = self.secondDerivateTermOffDiagonalAtomsIJ(z_i, z_j, z_i, z_j, r_ij, r_ij_b)
		return deltaTerm
	
	def getCalphaPosition(self, atom1, reference):
		""" Returns the position of atom1 among the calphas of reference. Useful if one 
		desires to know the index of an calpha atom in the ANM hessian made from reference calphas.
		
			Args:
				atom1: the calpha atom that the position is desired to know
				reference: the reference structure where the calpha position is obtained from
				
			Returns: Positive integer denoting the calpha position 
		"""
		
		assert atom1.getName() == 'CA'
		referenceCalphas = reference.select('calpha')
#         try:
#             idx = zip(referenceCalphas).index((atom1, ))
#             return idx
#         except ValueError:
#             print "Exception in getCalphaPosition. This calpha cannot be located in the structure provided. "
		for idx, referenceCalpha in enumerate(referenceCalphas):
			if atom1 == referenceCalpha:
				return idx
		raise StopIteration("Exception in getCalphaPosition. This calpha cannot be located in the structure provided. ")
			
	def normalizeM(self, M):
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
				if Mnormed is None:
					Mnormed = modeNormalized.getArray()
				else:
					Mnormed = np.column_stack((Mnormed, modeNormalized.getArray()))
		return Mnormed
	
	def getNoOfZeroEigvals(self, anm):
		""" Return the number of zero eigenvalues, the treshold is defined in the constant ZERO.
		
		Args:
			anm: the anm
			
		Returns: number of zero eigenvalues
		
		"""
		ZERO = 1e-10
		return sum(anm.getEigvals() < ZERO)
	
	def removeInterAtoms(self, arr, interCalphaIndices):
		""" Set x,y,z coordinations of atoms indicated by calphasInterIndices to 0,0,0 in arr. 
		
		Args:
			arr: the array with x,y,z coordinates
			interCalphaIndices: calphas with intermolecular contacts
			
		Returns: arr with x,y,z positions of atoms from interCalphaIndices set to 0,0,0
		"""
		for calphaIndex in interCalphaIndices:
			arr[(calphaIndex*3)] = 0.0
			arr[(calphaIndex*3+1)] = 0.0
			arr[(calphaIndex*3+2)] = 0.0
		return arr
		