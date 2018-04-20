'''
Created on Jan 19, 2014

@author: not
'''

import traceback
import numpy as np
from subprocess import call
import glob
import os
import sys
import gzip
from prody.dynamics.nmdfile import writeNMD
from prody.proteins.pdbfile import writePDB, writePDBStream, parsePDBStream
from prody.dynamics.functions import saveModel
from shutil import copy2
import subprocess
from helperScripts.scriptutils import runBashCommand
import StringIO

class ResultsPrinter(object):
	'''
	This class holds the direct experimental results of a 
	complex bound to unbound protein investigation with regards to NMA.
	'''


	def __init__(self, pdbName):
		'''
		Constructor, set the pdbName of the experimental results.
		'''
		self.pdbName = pdbName

	def setRMSDtostart(self, RMSDtostart):
		""" The sampling ensemble based on predicted RMSD """
		self.RMSDtostart = RMSDtostart

	def setRMSDtobound(self, RMSDtobound):
		""" The sampling ensemble based on predicted RMSD """
		self.RMSDtobound = RMSDtobound

	def setEnsemble(self, Ensemble):
		""" The sampling ensemble based on predicted RMSD """
		self.Ensemble = Ensemble

	def setRMSDPrediction(self, RMSDPrediction):
		""" The predicted RMSD basd on the reference eigenvalues. """
		self.RMSDPrediction = RMSDPrediction 
 
	def setRMSDReductionsWhole(self, RMSDReductionsWhole):
		""" RMSD reductions via the Swarmdock Tapprox calculation. """
		self.RMSDReductionsWhole = RMSDReductionsWhole
	
	def setOverlapTApproxWhole(self, overlapTApproxWhole):
		""" Overlap of the Tapprox as number of modes increases. """
		self.overlapTApproxWhole = overlapTApproxWhole
	
	def setStepPointsReductionWhole(self, stepPointsReductionWhole):
		""" Step points of the modes, how many modes have been used at position i 
		in self.RMSDReductionsWhole """
		self.stepPointsReductionWhole = stepPointsReductionWhole
		
	def setRMSDReductionsInterface(self, RMSDReductionsInterface):
		""" RMSD reductions via the Swarmdock Tapprox calculation. """
		self.RMSDReductionsInterface = RMSDReductionsInterface
	
	def setOverlapTApproxInterface(self, overlapTApproxInterface):
		""" Overlap of the Tapprox as number of modes increases. """
		self.overlapTApproxInterface = overlapTApproxInterface
	
	def setStepPointsReductionInterface(self, stepPointsReductionInterface):
		""" Step points of the modes, how many modes have been used at position i 
		in self.RMSDReductionsInterface """
		self.stepPointsReductionInterface = stepPointsReductionInterface      
 
	def setRMSDReductionsComplex2kWhole(self, RMSDReductionsComplex2kWhole):
		""" RMSD reductions via the Swarmdock Tapprox calculation. """
		self.RMSDReductionsComplex2kWhole = RMSDReductionsComplex2kWhole
	
	def setOverlapTApproxComplex2kWhole(self, overlapTApproxComplex2kWhole):
		""" Overlap of the Tapprox as number of modes increases. """
		self.overlapTApproxComplex2kWhole = overlapTApproxComplex2kWhole
	
	def setStepPointsReductionComplex2kWhole(self, stepPointsReductionComplex2kWhole):
		""" Step points of the modes, how many modes have been used at position i 
		in self.RMSDReductionsComplex2kWhole """
		self.stepPointsReductionComplex2kWhole = stepPointsReductionComplex2kWhole
		
	def setRMSDReductionsComplex1k1kWhole(self, RMSDReductionsComplex1k1kWhole):
		""" RMSD reductions via the Swarmdock Tapprox calculation  . """
		self.RMSDReductionsComplex1k1kWhole = RMSDReductionsComplex1k1kWhole
	
	def setOverlapTApproxComplex1k1kWhole(self, overlapTApproxComplex1k1kWhole):
		""" Overlap of the TapproxInterface as number of modes increases. """
		self.overlapTApproxComplex1k1kWhole = overlapTApproxComplex1k1kWhole
	
	def setStepPointsReductionComplex1k1kWhole(self, stepPointsReductionComplex1k1kWhole):
		""" Step points of the modes, how many modes have been used at position i 
		in self.RMSDReductionsComplex1k1kWhole """
		self.stepPointsReductionComplex1k1kWhole = stepPointsReductionComplex1k1kWhole
		
	def setRMSDReductionsComplex2kInterface(self, RMSDReductionsComplex2kInterface):
		""" RMSD reductions via the Swarmdock Tapprox calculation. """
		self.RMSDReductionsComplex2kInterface = RMSDReductionsComplex2kInterface
	
	def setOverlapTApproxComplex2kInterface(self, overlapTApproxComplex2kInterface):
		""" Overlap of the Tapprox as number of modes increases. """
		self.overlapTApproxComplex2kInterface = overlapTApproxComplex2kInterface
	
	def setStepPointsReductionComplex2kInterface(self, stepPointsReductionComplex2kInterface):
		""" Step points of the modes, how many modes have been used at position i 
		in self.RMSDReductionsComplex2kInterface """
		self.stepPointsReductionComplex2kInterface = stepPointsReductionComplex2kInterface
		
	def setRMSDReductionsComplex1k1kInterface(self, RMSDReductionsComplex1k1kInterface):
		""" RMSD reductions via the Swarmdock Tapprox calculation  . """
		self.RMSDReductionsComplex1k1kInterface = RMSDReductionsComplex1k1kInterface
	
	def setOverlapTApproxComplex1k1kInterface(self, overlapTApproxComplex1k1kInterface):
		""" Overlap of the TapproxInterface as number of modes increases. """
		self.overlapTApproxComplex1k1kInterface = overlapTApproxComplex1k1kInterface
	
	def setStepPointsReductionComplex1k1kInterface(self, stepPointsReductionComplex1k1kInterface):
		""" Step points of the modes, how many modes have been used at position i 
		in self.RMSDReductionsComplex1k1kInterface """
		self.stepPointsReductionComplex1k1kInterface = stepPointsReductionComplex1k1kInterface
		
	def setRMSD_unbound_to_superposed_bound(self, RMSD_unbound_to_superposed_bound):
		""" Initial RMSD from the unbound to bound conformation of the protein. """
		self.RMSD_unbound_to_superposed_bound = RMSD_unbound_to_superposed_bound

	def setC_RMSD_unbound_to_superposed_bound(self, C_RMSD_unbound_to_superposed_bound):
		""" Initial RMSD from the unbound to bound conformation of the protein. """
		self.C_RMSD_unbound_to_superposed_bound = C_RMSD_unbound_to_superposed_bound

	def setRMSD_interface(self, RMSD_interface):
		""" Initial RMSD from the unbound to bound conformation of the interface. """
		self.RMSD_interface = RMSD_interface
		
	def setOverlapArrayWhole(self, overlapArrayWhole):
		""" The overlap of modes with the deformation vector on the protein. """
		self.overlapArrayWhole = overlapArrayWhole
		
	def setOverlapArrayInterface(self, overlapArrayInterface):
		""" The overlap of modes with the deformation vector on the interface. """
		self.overlapArrayInterface = overlapArrayInterface
		
	def setCollectivityArrayWhole(self, collectivityArrayWhole):
		""" The collectivity of modes considering atoms of the protein. """
		self.collectivityArrayWhole = collectivityArrayWhole
		
	def setCollectivityArrayInterface(self, collectivityArrayInterface):
		""" The collectivity of modes considering interface atoms. """
		self.collectivityArrayInterface = collectivityArrayInterface
		
	def setCorrelationArrayWhole(self, correlationArrayWhole):
		""" The correlations of magnitudes of modes on the protein. """
		self.correlationArrayWhole = correlationArrayWhole
		
	def setCorrelationArrayInterface(self, correlationArrayInterface):
		""" The correlations of magnitudes of modes on the interface. """
		self.correlationArrayInterface = correlationArrayInterface
		
	def setCumulOverlapWholePrody(self, cumulOverlapWholePrody):
		""" Cumulative overlap of all modes and the defvec via prodys cumulative overlap method."""
		self.cumulOverlapWholePrody = cumulOverlapWholePrody
		
	def setCumulOverlapInterfacePrody(self, cumulOverlapInterfacePrody):
		""" Cumulative overlap of all modes and the defvecInterface via prodys cumulative overlap method."""
		self.cumulOverlapInterfacePrody = cumulOverlapInterfacePrody
	
	def setEigenvaluesComplex(self, eigenvaluesComplex):
		""" The corresponding eigenvalues of the complex ANM. """
		self.eigenvaluesComplex = eigenvaluesComplex

	def setEigenvaluesReference(self, eigenvaluesReference):
		""" The corresponding eigenvalues of the reference ANM. """
		self.eigenvaluesReference = eigenvaluesReference

	def setEigenvectorsReference(self, eigenvectorsReference):
		""" The corresponding eigenvectors of the reference ANM. """
		self.eigenvectorsReference = eigenvectorsReference

	def setEigenvectorsComplex(self, eigenvectorsComplex):
		""" The corresponding eigenvectors of the reference ANM. """
		self.eigenvectorsComplex = eigenvectorsComplex
	
	def setEigenvaluesComplex2k(self, eigenvaluesComplex2k):
		""" The corresponding eigenvaluesComplex2k of the ANM. """
		self.eigenvaluesComplex2k = eigenvaluesComplex2k
		
	def setEigenvaluesReceptor1k(self, eigenvaluesReceptor1k):
		""" The corresponding eigenvaluesReceptor1k of the ANM. """
		self.eigenvaluesReceptor1k = eigenvaluesReceptor1k
		
	def setEigenvaluesLigand1k(self, eigenvaluesLigand1k):
		""" The corresponding eigenvaluesLigand1k of the ANM. """
		self.eigenvaluesLigand1k = eigenvaluesLigand1k
		
	def setNumberOfModes(self, numberOfModes):
		""" The number of modes, if calculated on calphas, it is the 
		number of calpha atoms * 3 - 6. """
		self.numberOfModes = numberOfModes
		
	def setNumberOfModesComplex(self, numberOfModesComplex):
		""" The number of modes, if calculated on calphas, it is the 
		number of calpha atoms * 3 - 6. """
		self.numberOfModesComplex = numberOfModesComplex
		
	def setNumberOfModesInterface(self, numberOfModesInterface):
		""" If modes were to be calculate on the interface, what is the
		number of calpha atoms * 3 - 6 of the interface. """
		self.numberOfModesInterface = numberOfModesInterface
		
	def setRMSDReductionsWholeByCollectivity(self, RMSDReductionsWholeByCollectivity):
		self.RMSDReductionsWholeByCollectivity = RMSDReductionsWholeByCollectivity
		
	def setRMSDReductionsInterfaceByCollectivity(self, RMSDReductionsInterfaceByCollectivity):
		self.RMSDReductionsInterfaceByCollectivity = RMSDReductionsInterfaceByCollectivity
		
	def setOverlapTApproxWholeByCollectivity(self, overlapTApproxWholeByCollectivity):
		self.overlapTApproxWholeByCollectivity = overlapTApproxWholeByCollectivity
		
	def setOverlapTApproxInterfaceByCollectivity(self, overlapTApproxInterfaceByCollectivity):
		self.overlapTApproxInterfaceByCollectivity = overlapTApproxInterfaceByCollectivity
		
	def setStepPointsReductionWholeByCollectivity(self, stepPointsReductionWholeByCollectivity):
		self.stepPointsReductionWholeByCollectivity = stepPointsReductionWholeByCollectivity
		
	def setStepPointsReductionInterfaceByCollectivity(self, stepPointsReductionInterfaceByCollectivity):
		self.stepPointsReductionInterfaceByCollectivity = stepPointsReductionInterfaceByCollectivity
		
	def setBoundComplex(self, boundComplex):
		self.boundComplex = boundComplex
		
	def setUnboundComplexAligned(self, unboundComplexAligned):
		self.unboundComplexAligned = unboundComplexAligned
		
	def setLRMS(self, L_rms):
		self.L_rms = L_rms
		
	def setCounterpartRMS(self, counterpart_rms):
		self.counterpart_rms = counterpart_rms
		
	def setIRMSbeforeAlign(self, I_rms_before_align):
		self.I_rms_before_align = I_rms_before_align
		
	def setIRMSafterAlign(self, I_rms_after_align):
		self.I_rms_after_align = I_rms_after_align
		
	def setcCase(self, cCase):
		self.cCase = cCase
		
	def setPathOfConfigFile(self, pathOfConfigFile):
		self.pathOfConfigFile = pathOfConfigFile
			
	def setOverlapTable(self, overlapTable):
		self.overlapTable = overlapTable
		
	def setSubSpaceOverlaps(self, subspaceOverlaps, subspaceOverlapsRanges):
		self.subspaceOverlaps = subspaceOverlaps
		self.subspaceOverlapsRanges = subspaceOverlapsRanges
		
	def setCovarianceOverlap(self, covarianceOverlaps, covarianceOverlapsRanges):
		self.covarianceOverlaps = covarianceOverlaps
		self.covarianceOverlapsRanges = covarianceOverlapsRanges        
		
	def setZeroEigvecsProtein1(self, zeroEigvecs):
		self.zeroEigvecsProtein1 = zeroEigvecs
		
	def setZeroEigvecsProtein2(self, zeroEigvecs):
		self.zeroEigvecsProtein2 = zeroEigvecs
		
	def setZeroEigvecsComplex(self, zeroEigvecs):
		self.zeroEigvecsComplex = zeroEigvecs
		
	def setReference(self, reference):
		self.reference = reference
		
	def setMobile(self, mobile):
		self.mobile = mobile
		
	def setUnboundCounterpart(self, unboundCounterpart):
		self.unboundCounterpart = unboundCounterpart
		
	def setBoundCounterpart(self, boundCountertpart):
		self.boundCountertpart = boundCountertpart

	def setRefChain(self, refChain):
		self.refChain = refChain
		
	def setRefChainInterface(self, refChainInterface):
		self.refChainInterface = refChainInterface
		
	def setMobChain(self, mobChain):
		self.mobChain = mobChain
		
	def setMobChainInterface(self, mobChainInterface):
		self.mobChainInterface = mobChainInterface
		
	def setUnboundCounterpartChain(self, unboundCounterpartChain):
		self.unboundCounterpartChain = unboundCounterpartChain
		
	def setUnboundCounterpartChainInterface(self, unboundCounterpartChainInterface):
		self.unboundCounterpartChainInterface = unboundCounterpartChainInterface
		
	def setBoundCounterpartChain(self, boundCounterpartChain):
		self.boundCounterpartChain = boundCounterpartChain
		
	def setBoundCounterpartChainInterface(self, boundCounterpartChainInterface):
		self.boundCounterpartChainInterface = boundCounterpartChainInterface
		
	def setUnboundComplexAlignedChain(self, unboundComplexAlignedChain):
		self.unboundComplexAlignedChain = unboundComplexAlignedChain
		
	def setUnboundComplexChainInterface(self, unboundComplexChainInterface):
		self.unboundComplexChainInterface = unboundComplexChainInterface
	
	def setBoundComplexChain(self, boundComplexChain):
		self.boundComplexChain = boundComplexChain
		
	def setBoundComplexChainInterface(self, boundComplexChainInterface):
		self.boundComplexChainInterface = boundComplexChainInterface
		
	def setL_RMSReductions(self, L_RMSReductions):
		self.L_RMSReductions = L_RMSReductions
		
	def setL_RMSD_unbound_to_superposed_bound(self, L_RMSD_unbound_to_superposed_bound):
		self.L_RMSD_unbound_to_superposed_bound = L_RMSD_unbound_to_superposed_bound
		
	def setSingleModeOverlapsFromSuperset(self, singleModeOverlapsFromSuperset):
		self.singleModeOverlapsFromSuperset = singleModeOverlapsFromSuperset
		
	def setDeformationSnapshots(self, deformationSnapshots):
		self.deformationSnapshots = deformationSnapshots
		
	def setIndicesOfLambdaRSorting(self, indicesOfLambdaRSorting):
		self.indicesOfLambdaRSorting = indicesOfLambdaRSorting


	def writeSampleResults(self, basePath, experimentName, utils, cCase=None):

		np.set_printoptions(threshold=np.nan)
		utils.mkdir_p(basePath)
		utils.mkdir_p(basePath+experimentName+"/")

		if not cCase: 
			path = basePath+experimentName+"/"
		else:
			path = basePath+experimentName+"/"+str(cCase)+"/"
		utils.mkdir_p(path)
		
		try: 
			writePDB(path+"Ensemble.pdb", self.Ensemble)
		except AttributeError, err:
			print "Exception AttributeError occurred: ", err
			print traceback.format_exc()

	def writeRMSDResults(self, basePath, experimentName, utils, cCase=None):

		np.set_printoptions(threshold=np.nan)
		utils.mkdir_p(basePath)
		utils.mkdir_p(basePath+experimentName+"/")

		if not cCase: 
			path = basePath+experimentName+"/"
		else:
			path = basePath+experimentName+"/"+str(cCase)+"/"
		utils.mkdir_p(path)
		
		try: 
			np.savetxt(path+"RMSDPrediction.txt", self.RMSDPrediction, fmt='%15.15f', header='RMSDPrediction')
		except AttributeError, err:
			print "Exception AttributeError occurred: ", err
			print traceback.format_exc()

	def writeComplexRMSD(self, basePath, experimentName, utils, cCase=None):

		np.set_printoptions(threshold=np.nan)
		utils.mkdir_p(basePath)
		utils.mkdir_p(basePath+experimentName+"/")

		if not cCase: 
			path = basePath+experimentName+"/"
		else:
			path = basePath+experimentName+"/"+str(cCase)+"/"
		utils.mkdir_p(path)
		
		try: 
			np.savetxt(path+"C_RMSD_unbound_to_superposed_bound.txt", [self.C_RMSD_unbound_to_superposed_bound], fmt='%15.15f', header='C_RMSD_unbound_to_superposed_bound')
		except AttributeError, err:
			print "Exception AttributeError occurred: ", err
			print traceback.format_exc()
	 
	def writeDirectResults(self, basePath, experimentName, utils, bound_provided, investigation, cCase=None):
		""" Write the direct results to individual text files. 
		
		Args:
			basePath: the base path to the general output folder
			experimentName: the overall experiment name, one folderlevel below basePath
			utils: utilities object
			cCase: deprecated attribute
		
		Returns: path to results
		"""
		np.set_printoptions(threshold=np.nan)
		utils.mkdir_p(basePath)
		utils.mkdir_p(basePath+experimentName+"/")

		if not cCase: 
			path = basePath+experimentName+"/"
		else:
			path = basePath+experimentName+"/"+str(cCase)+"/"
		utils.mkdir_p(path)

		# try: 
		# 	np.savetxt(path+"SamplingRMSDtostart.txt", self.RMSDtostart, fmt='%15.15f', header='SamplingRMSDtostart')
		# except AttributeError, err:
		# 	print "Exception AttributeError occurred: ", err
		# 	print traceback.format_exc()

		# try: 
		# 	np.savetxt(path+"SamplingRMSDtobound.txt", self.RMSDtobound, fmt='%15.15f', header='SamplingRMSDtobound')
		# except AttributeError, err:
		# 	print "Exception AttributeError occurred: ", err
		# 	print traceback.format_exc()

		# try: 
		# 	np.savetxt(path+"RMSDPrediction.txt", self.RMSDPrediction, fmt='%15.15f', header='RMSDPrediction')
		# except AttributeError, err:
		# 	print "Exception AttributeError occurred: ", err
		# 	print traceback.format_exc()

		# try:
		# 	with file(path+"Ensemble.txt", "w") as outfile:
		# 		outfile.write('# Ensemble : {0}\n'.format(self.Ensemble.shape))
		# 		for sample_index in range(0, 10):
		# 			outfile.write('# ensemble#{0}\n'.format(sample_index + 1))
		# 			np.savetxt(outfile, self.Ensemble[sample_index])
		# except AttributeError, err:
		# 	print "Exception AttributeError occurred: ", err
		# 	print traceback.format_exc()



		try: 
			np.savetxt(path+"numberOfModes.txt", [self.numberOfModes], fmt='%d', header='numberOfModes')
		except AttributeError, err:
			print "Exception AttributeError occurred: ", err
			print traceback.format_exc()  

		try: 
			np.savetxt(path+"pdbName.txt", [self.pdbName], fmt='%s', header='pdbName')
		except AttributeError, err:
			print "Exception AttributeError occurred: ", err
			print traceback.format_exc()
#		if investigation == "Complex":

		try: 
			np.savetxt(path+"eigenvaluesComplex.txt", self.eigenvaluesComplex, fmt='%15.15f', header='eigenvaluesComplex')
		except AttributeError, err:
			print "Exception AttributeError occurred: ", err
			print traceback.format_exc()    

		try: 
			np.savetxt(path+"eigenvectorsComplex.txt", self.eigenvectorsComplex, fmt='%15.15f', header='eigenvectorsComplex')
		except AttributeError, err:
			print "Exception AttributeError occurred: ", err
			print traceback.format_exc()
#		else:
		try: 
			np.savetxt(path+"eigenvaluesReference.txt", self.eigenvaluesReference, fmt='%15.15f', header='eigenvaluesReference')
		except AttributeError, err:
			print "Exception AttributeError occurred: ", err
			print traceback.format_exc() 
		try: 
			np.savetxt(path+"eigenvectorsReference.txt", self.eigenvectorsReference, fmt='%15.15f', header='eigenvectorsReference')
		except AttributeError, err:
			print "Exception AttributeError occurred: ", err
			print traceback.format_exc()  
		# try: 
		# 	np.savetxt(path+"zeroEigvecsProtein1.dat", self.zeroEigvecsProtein1.round(3), fmt='%.3f')
		# except AttributeError, err:
		# 	print "Exception AttributeError occurred: ", err
		# 	print traceback.format_exc()
			
		# try: 
		# 	np.savetxt(path+"zeroEigvecsProtein2.dat", self.zeroEigvecsProtein2.round(3), fmt='%.3f')
		# except AttributeError, err:
		# 	print "Exception AttributeError occurred: ", err
		# 	print traceback.format_exc()      
			
		# try:   
		# 	np.savetxt(path+"zeroEigvecsComplex.dat", self.zeroEigvecsComplex.round(3), fmt='%.3f')
		# except AttributeError, err:
		# 	print "Exception AttributeError occurred: ", err
		# 	print traceback.format_exc()           
			 
		if bound_provided == True:
			
			try: 
				np.savetxt(path+"RMSDReductionsWhole.txt", self.RMSDReductionsWhole, fmt='%15.15f', header='RMSDReductionsWhole')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()
			except:
				pass


			try: 
				np.savetxt(path+"overlapTApproxWhole.txt", self.overlapTApproxWhole, fmt='%15.15f', header='overlapTApproxWhole')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()
			except:
				pass

			try: 
				np.savetxt(path+"stepPointsReductionWhole.txt", self.stepPointsReductionWhole, fmt='%d', header='stepPointsReductionWhole')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()
			except:
				pass

			try: 
				np.savetxt(path+"RMSDReductionsInterface.txt", self.RMSDReductionsInterface, fmt='%15.15f', header='RMSDReductionsInterface')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()
			except:
				pass
				
			try: 
				np.savetxt(path+"overlapTApproxInterface.txt", self.overlapTApproxInterface, fmt='%15.15f', header='overlapTApproxInterface')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()
			except:
				pass				
			try: 
				np.savetxt(path+"stepPointsReductionInterface.txt", self.stepPointsReductionInterface, fmt='%d', header='stepPointsReductionInterface')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()
			except:
				pass	 
			# try: 
			# 	np.savetxt(path+"RMSDReductionsComplex2kWhole.txt", self.RMSDReductionsComplex2kWhole, fmt='%15.15f', header='RMSDReductionsComplex2kWhole')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()

			# try: 
			# 	np.savetxt(path+"overlapTApproxComplex2kWhole.txt", self.overlapTApproxComplex2kWhole, fmt='%15.15f', header='overlapTApproxComplex2kWhole')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()
				
			# try: 
			# 	np.savetxt(path+"stepPointsReductionComplex2kWhole.txt", self.stepPointsReductionComplex2kWhole, fmt='%d', header='stepPointsReductionComplex2kWhole')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()

			# try: 
			# 	np.savetxt(path+"RMSDReductionsComplex1k1kWhole.txt", self.RMSDReductionsComplex1k1kWhole, fmt='%15.15f', header='RMSDReductionsComplex1k1kWhole')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()

			# try: 
			# 	np.savetxt(path+"overlapTApproxComplex1k1kWhole.txt", self.overlapTApproxComplex1k1kWhole, fmt='%15.15f', header='overlapTApproxComplex1k1kWhole')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()

			# try: 
			# 	np.savetxt(path+"stepPointsReductionComplex1k1kWhole.txt", self.stepPointsReductionComplex1k1kWhole, fmt='%d', header='stepPointsReductionComplex1k1kWhole')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()
				
			# try: 
			# 	np.savetxt(path+"RMSDReductionsComplex2kInterface.txt", self.RMSDReductionsComplex2kInterface, fmt='%15.15f', header='RMSDReductionsComplex2kInterface')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()

			# try: 
			# 	np.savetxt(path+"overlapTApproxComplex2kInterface.txt", self.overlapTApproxComplex2kInterface, fmt='%15.15f', header='overlapTApproxComplex2kInterface')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()
				
			# try: 
			# 	np.savetxt(path+"stepPointsReductionComplex2kInterface.txt", self.stepPointsReductionComplex2kInterface, fmt='%d', header='stepPointsReductionComplex2kInterface')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()

			# try: 
			# 	np.savetxt(path+"RMSDReductionsComplex1k1kInterface.txt", self.RMSDReductionsComplex1k1kInterface, fmt='%15.15f', header='RMSDReductionsComplex1k1kInterface')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()

			# try: 
			# 	np.savetxt(path+"overlapTApproxComplex1k1kInterface.txt", self.overlapTApproxComplex1k1kInterface, fmt='%15.15f', header='overlapTApproxComplex1k1kInterface')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()

			try: 
				np.savetxt(path+"stepPointsReductionComplex1k1kInterface.txt", self.stepPointsReductionComplex1k1kInterface, fmt='%d', header='stepPointsReductionComplex1k1kInterface')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()
			except:
				pass
			try: 
				np.savetxt(path+"RMSD_unbound_to_superposed_bound.txt", [self.RMSD_unbound_to_superposed_bound], fmt='%15.15f', header='RMSD_unbound_to_superposed_bound')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()
			except:
				pass
			
			try: 
				np.savetxt(path+"RMSD_interface.txt", [self.RMSD_interface], fmt='%15.15f', header='RMSD_interface')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()        
			except:
				pass			
			try: 
				np.savetxt(path+"overlapArrayWhole.txt", self.overlapArrayWhole, fmt='%15.15f', header='overlapArrayWhole')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()        
			except:
				pass			
			try: 
				np.savetxt(path+"overlapArrayInterface.txt", self.overlapArrayInterface, fmt='%15.15f', header='overlapArrayInterface')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()            
			except:
				pass			
			try: 
				np.savetxt(path+"collectivityArrayWhole.txt", self.collectivityArrayWhole, fmt='%15.15f', header='collectivityArrayWhole')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()  
			except:
				pass			
			try: 
				np.savetxt(path+"collectivityArrayInterface.txt", self.collectivityArrayInterface, fmt='%15.15f', header='collectivityArrayInterface')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()
			except:
				pass			
			try: 
				np.savetxt(path+"correlationArrayWhole.txt", self.correlationArrayWhole, fmt='%15.15f', header='correlationArrayWhole')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()         
			except:
				pass			
			try: 
				np.savetxt(path+"correlationArrayInterface.txt", self.correlationArrayInterface, fmt='%15.15f', header='correlationArrayInterface')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()             
			except:
				pass			
			try: 
				np.savetxt(path+"cumulOverlapWholePrody.txt", self.cumulOverlapWholePrody, fmt='%15.15f', header='cumulOverlapWholePrody')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()
			except:
				pass			
			try: 
				np.savetxt(path+"cumulOverlapInterfacePrody.txt", self.cumulOverlapInterfacePrody, fmt='%15.15f', header='cumulOverlapInterfacePrody')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()
			except:
				pass				
			# try: 
			# 	np.savetxt(path+"eigenvaluesComplex2k.txt", self.eigenvaluesComplex2k, fmt='%15.15f', header='eigenvaluesComplex2k')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()   
				
			# try: 
			# 	np.savetxt(path+"eigenvaluesReceptor1k.txt", self.eigenvaluesReceptor1k, fmt='%15.15f', header='eigenvaluesReceptor1k')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()  
				
			try: 
				np.savetxt(path+"eigenvaluesComplex.txt", self.eigenvaluesComplex, fmt='%15.15f', header='eigenvaluesComplex')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()    
			except:
				pass				
			# try: 
			# 	np.savetxt(path+"eigenvaluesLigand1k.txt", self.eigenvaluesLigand1k, fmt='%15.15f', header='eigenvaluesLigand1k')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()  
			   
				
			try: 
				np.savetxt(path+"numberOfModesComplex.txt", [self.numberOfModesComplex], fmt='%d', header='numberOfModesComplex')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()         
			except:
				pass				
			try: 
				np.savetxt(path+"numberOfModesInterface.txt", [self.numberOfModesInterface], fmt='%d', header='numberOfModesInterface')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()               
			except:
				pass				
			# RMSD aided by collectivity
			try: 
				np.savetxt(path+"RMSDReductionsWholeByCollectivity.txt", self.RMSDReductionsWholeByCollectivity, fmt='%15.15f', header='RMSDReductionsWholeByCollectivity')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()
			except:
				pass
			try: 
				np.savetxt(path+"overlapTApproxWholeByCollectivity.txt", self.overlapTApproxWholeByCollectivity, fmt='%15.15f', header='overlapTApproxWholeByCollectivity')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()
			except:
				pass				
			try: 
				np.savetxt(path+"stepPointsReductionWholeByCollectivity.txt", self.stepPointsReductionWholeByCollectivity, fmt='%d', header='stepPointsReductionWholeByCollectivity')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()
			except:
				pass
			try: 
				np.savetxt(path+"RMSDReductionsInterfaceByCollectivity.txt", self.RMSDReductionsInterfaceByCollectivity, fmt='%15.15f', header='RMSDReductionsInterfaceByCollectivity')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()
			except:
				pass
			try: 
				np.savetxt(path+"overlapTApproxInterfaceByCollectivity.txt", self.overlapTApproxInterfaceByCollectivity, fmt='%15.15f', header='overlapTApproxInterfaceByCollectivity')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()
			except:
				pass
			try: 
				np.savetxt(path+"stepPointsReductionInterfaceByCollectivity.txt", self.stepPointsReductionInterfaceByCollectivity, fmt='%d', header='stepPointsReductionInterfaceByCollectivity')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()
			except:
				pass				
			# RMSD complex assessments
			try: 
				np.savetxt(path+"L_rms.txt", [self.L_rms], fmt='%15.15f', header='L_rms')
			except AttributeError, err:
				print "Exception AttributeError occurred: ", err
				print traceback.format_exc()
				
			# try: 
			# 	np.savetxt(path+"counterpart_rms.txt", [self.counterpart_rms], fmt='%15.15f', header='counterpart_rms')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()
				
			try: 
				np.savetxt(path+"I_rms_after_align.txt", [self.I_rms_after_align], fmt='%15.15f', header='I_rms_after_align')
			except AttributeError, err:
				print "Exception AttributeError occurred: ", err
				print traceback.format_exc()
				
			# # measures of comparing hessians
			# try: 
			# 	np.savetxt(path+"subspaceOverlaps.txt", self.subspaceOverlaps, fmt='%15.15f', header='subspaceOverlaps')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()    
			# try: 
			# 	np.savetxt(path+"subspaceOverlapsRanges.txt", self.subspaceOverlapsRanges, fmt='%d', header='subspaceOverlapsRanges')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()      

			# try: 
			# 	np.savetxt(path+"covarianceOverlaps.txt", self.covarianceOverlaps, fmt='%15.15f', header='covarianceOverlaps')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()    
			# try: 
			# 	np.savetxt(path+"covarianceOverlapsRanges.txt", self.covarianceOverlapsRanges, fmt='%d', header='covarianceOverlapsRanges')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()  
				

			# try: 
			# 	with open(path+"overlapTable.prody", "w") as text_file:
			# 		text_file.write(self.overlapTable)
			# except Exception, err:
			# 	print "Exception occurred: ", err
			# 	print traceback.format_exc()       
			# 	pass
			
			try: 
				np.savetxt(path+"L_RMSReductions.txt", self.L_RMSReductions, fmt='%15.15f', header='L_RMSReductions')
			except AttributeError, err:
				print "Exception AttributeError occurred: ", err
				print traceback.format_exc()
				
			try: 
				np.savetxt(path+"L_RMSD_unbound_to_superposed_bound.txt", [self.L_RMSD_unbound_to_superposed_bound], fmt='%15.15f', header='L_RMSD_unbound_to_superposed_bound')
			except AttributeError, err:
				print "Exception AttributeError occurred: ", err
				print traceback.format_exc()

			try: 
				np.savetxt(path+"singleModeOverlapsFromSuperset.txt", self.singleModeOverlapsFromSuperset, fmt='%15.15f', header='singleModeOverlapsFromSuperset')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()
			except:
				pass				
			try: 
				np.savetxt(path+"indicesOfLambdaRSorting.txt", self.indicesOfLambdaRSorting, fmt='%d', header='self.indicesOfLambdaRSorting')
			# except AttributeError, err:
			# 	print "Exception AttributeError occurred: ", err
			# 	print traceback.format_exc()        
			except:
				pass					
		# # save PDBs
		# try: 
		# 	pdbspath = path+"pdbs/"
		# 	utils.mkdir_p(pdbspath)
		# 	writePDB(pdbspath+self.reference.getTitle()+".pdb", self.reference)
		# 	writePDB(pdbspath+self.reference.getTitle()+"_chains.pdb", self.refChain)
		# 	writePDB(pdbspath+self.reference.getTitle()+"_chains_interface.pdb", self.refChainInterface)
		# 	writePDB(pdbspath+self.mobile.getTitle()+".pdb", self.mobile)
		# 	writePDB(pdbspath+self.mobile.getTitle()+"_chains.pdb", self.mobChain)
		# 	writePDB(pdbspath+self.mobile.getTitle()+"_chains_interface.pdb", self.mobChainInterface)
		# 	writePDB(pdbspath+self.unboundCounterpart.getTitle()+".pdb", self.unboundCounterpart)
		# 	writePDB(pdbspath+self.unboundCounterpart.getTitle()+"_chains.pdb", self.unboundCounterpartChain)
		# 	writePDB(pdbspath+self.unboundCounterpart.getTitle()+"_chains_interface.pdb", self.unboundCounterpartChainInterface)
		# 	writePDB(pdbspath+self.boundCountertpart.getTitle()+".pdb", self.boundCountertpart)
		# 	writePDB(pdbspath+self.boundCountertpart.getTitle()+"_chains.pdb", self.boundCounterpartChain)
		# 	writePDB(pdbspath+self.boundCountertpart.getTitle()+"_chains_interface.pdb", self.boundCounterpartChainInterface)
		# 	writePDB(pdbspath+self.unboundComplexAligned.getTitle()+"_ucomplex.pdb", self.unboundComplexAligned)
		# 	writePDB(pdbspath+self.unboundComplexAligned.getTitle()+"_ucomplex_chains.pdb", self.unboundComplexAlignedChain)
		# 	writePDB(pdbspath+self.unboundComplexAligned.getTitle()+"_ucomplex_chains_interface.pdb", self.unboundComplexChainInterface)
		# 	writePDB(pdbspath+self.boundComplex.getTitle()+"_bcomplex.pdb", self.boundComplex)
		# 	writePDB(pdbspath+self.boundComplex.getTitle()+"_bcomplex_chains.pdb", self.boundComplexChain)
		# 	writePDB(pdbspath+self.boundComplex.getTitle()+"_bcomplex_chains_interface.pdb", self.boundComplexChainInterface)
		# except Exception, err:
		# 	print "Exception occurred when writing pdbs: ", err
		# 	print traceback.format_exc()       
		# 	pass
		
		# save deformation snapshots
		# try: 
		# 	deformpath = path+"pdbs/deformationsnapshots/"
		# 	utils.mkdir_p(deformpath)
		# 	for k in self.deformationSnapshots.keys():
		# 		writePDB(deformpath+"protein_deformed_"+str(k)+".pdb", self.deformationSnapshots[k])
		# except Exception, err:
		# 	print "Exception occurred when writing deformation snapshot pdbs: ", err
		# 	print traceback.format_exc()       
		# 	pass
		
		# # make frame PDB
		# try:
		# 	os.chdir(deformpath)
		# 	shellCommand = "../../../../helperScripts/makeFramePDB.sh "+deformpath
		# 	os.system(shellCommand)
		# 	# change the working path back to the directory of the program
		# 	os.chdir(os.path.dirname(os.path.abspath(__file__)))
		# except Exception, err:
		# 	print "Exception occurred when writing the frame PDB files: ", err
		# 	print traceback.format_exc()
			
		# compress the pdb folder and delete the folder
		# try:
		# 	os.chdir(path)
		# 	shellCommand = "tar -jcvf pdbs.tar.bz2 pdbs/"
		# 	os.system(shellCommand)
		# 	shellCommand = "rm -rf pdbs/"
		# 	os.system(shellCommand)            
		# 	# change the working path back to the directory of the program
		# 	os.chdir(os.path.dirname(os.path.abspath(__file__)))
		# except Exception, err:
		# 	print "Exception occurred when bz2ing and deleting the pdbs folder: ", err
		# 	print traceback.format_exc()        
						 
		# combine the files
		# try:
		# 	os.chdir(path)
		# 	filesToCombine = glob.glob("*.txt")
		# 	filesToCombine.sort()
		# 	combineArguments = " ".join(filesToCombine)
		# 	shellCommand = "paste "+combineArguments+" | column -s $'\\t' -t > combined.txt"
		# 	os.system(shellCommand)
		# 	# change the working path back to the directory of the program
		# 	os.chdir(os.path.dirname(os.path.abspath(__file__)))
		# except Exception, err:
		# 	print "Exception occurred when combining direct result files: ", err
		# 	print traceback.format_exc()
			
		# copy the config file to the result folder as a reference
		try:
			os.chdir(path)
			copy2(self.pathOfConfigFile, path)
			# change the working path back to the directory of the program
			os.chdir(os.path.dirname(os.path.abspath(__file__)))
		except Exception, err:
			print "Exception occurred when trying to copy the Configurations file: ", err
			print traceback.format_exc()
	
		return path            
			
			
	def writeNMDResults(self, basePath, experimentName, utils, ANMs, single, cCase=None, storeANMs=True):
		np.set_printoptions(threshold=np.nan)
		if not cCase:
			path = basePath+experimentName+"/"+self.pdbName+"/"
		else:
			path = basePath+experimentName+"/"+self.pdbName+"_"+str(cCase)+"/"
		utils.mkdir_p(path)

		try: 
			writeNMD(path+self.pdbName+"_reference_ANM.nmd", ANMs._anm_reference[0], ANMs._anm_reference[1])
		except Exception, err:
			print "Exception occurred when writing "+path+self.pdbName+"receptor_slc_ANM.nmd file: ", err
			print traceback.format_exc()

       
		if single == False:
			try: 
				fileName = "_complex"
				writeNMD(path+self.pdbName+fileName+"_slc_ANM.nmd", ANMs.getANMComplexSlc()[0], ANMs.getANMComplexSlc()[1])
			except Exception, err:
				print "Exception occurred when writing "+path+self.pdbName+fileName+"_slc_ANM.nmd file: ", err
				print traceback.format_exc()
			try: 
				fileName = "_complex"
				writeNMD(path+self.pdbName+fileName+"_ANM.nmd", ANMs._anm_complex[0], ANMs._anm_complex[1])
			except Exception, err:
				print "Exception occurred when writing "+path+self.pdbName+fileName+"_ANM.nmd file: ", err
				print traceback.format_exc() 
				
				
			# try: 
			# 	writeNMD(path+self.pdbName+"_counterpart_slc_ANM.nmd", ANMs._anm_counterpart_slc[0], ANMs._anm_counterpart_slc[1])
			# except Exception, err:
			# 	print "Exception occurred when writing "+path+self.pdbName+"ligand_slc_ANM.nmd file: ", err
			# 	print traceback.format_exc()
				
			try: 
				writeNMD(path+self.pdbName+"_counterpart_ANM.nmd", ANMs._anm_counterpart[0], ANMs._anm_counterpart[1])
			except Exception, err:
				print "Exception occurred when writing "+path+self.pdbName+"ligand_ANM.nmd file: ", err
				print traceback.format_exc()
				
			try:
				writeNMD(path+self.pdbName+"_complex_slc_interface_ANM.nmd", ANMs._anm_boundcomplex_slc_interface[0], ANMs._anm_boundcomplex_slc_interface[1])
			except Exception, err:
				pass
			#print "Exception occurred when writing "+path+self.pdbName+"_complex_slc_interface.nmd file: ", err
			#print traceback.format_exc()            
			
		# # add the zero eigenvalue modes to the nmd file
		# try:
		# 	os.chdir(path)
		# 	nmdFile = path+self.pdbName+"_reference_ANM.nmd"
		# 	shellCommand = "grep -v \"^mode\" "+nmdFile+" > temp1; grep \"^mode\" "+nmdFile+" > temp2"
		# 	os.system(shellCommand)
		# 	#shellCommand = '''awk \'{printf(\"mode %1d %s\\n\", NR, $0)}\' zeroEigvecs.dat > zerotemp'''
		# 	shellCommand = "awk \'{printf(\"mode %1d %.2f %s\\n\", NR, sqrt(1/(0.0001*NR)), $0)}\'  path + zeroEigvecsProtein1.dat' > zerotemp"
		# 	os.system(shellCommand)
		# 	shellCommand = "cat temp1 zerotemp temp2 > "+path+self.pdbName+"_reference_zeros_ANM.nmd ; rm temp1 zerotemp temp2"
		# 	os.system(shellCommand)
		# 	# change the working path back to the directory of the program
		# 	os.chdir(os.path.dirname(os.path.abspath(__file__)))
		# except Exception, err:
		# 	print "Exception occurred when adding zero eigenvalue eigenvectors to the reference nmd file: ", err
		# 	print traceback.format_exc()
			
		# try:
		# 	os.chdir(path)
		# 	nmdFile = path+self.pdbName+"_counterpart_ANM.nmd"
		# 	shellCommand = "grep -v \"^mode\" "+nmdFile+" > temp1; grep \"^mode\" "+nmdFile+" > temp2"
		# 	os.system(shellCommand)
		# 	#shellCommand = '''awk \'{printf(\"mode %1d %s\\n\", NR, $0)}\' zeroEigvecs.dat > zerotemp'''
		# 	shellCommand = '''awk \'{printf(\"mode %1d %.2f %s\\n\", NR, sqrt(1/(0.0001*NR)), $0)}\' zeroEigvecsProtein2.dat > zerotemp'''
		# 	os.system(shellCommand)
		# 	shellCommand = "cat temp1 zerotemp temp2 > "+path+self.pdbName+"_counterpart_zeros_ANM.nmd ; rm temp1 zerotemp temp2"
		# 	os.system(shellCommand)
		# 	# change the working path back to the directory of the program
		# 	os.chdir(os.path.dirname(os.path.abspath(__file__)))
		# except Exception, err:
		# 	print "Exception occurred when adding zero eigenvalue eigenvectors to the reference nmd file: ", err
		# 	print traceback.format_exc()
			
		# try:
		# 	os.chdir(path)
		# 	nmdFile = path+self.pdbName+"_complex_ANM.nmd"
		# 	shellCommand = "grep -v \"^mode\" "+nmdFile+" > temp1; grep \"^mode\" "+nmdFile+" > temp2"
		# 	os.system(shellCommand)
		# 	#shellCommand = '''awk \'{printf(\"mode %1d %s\\n\", NR, $0)}\' zeroEigvecs.dat > zerotemp'''
		# 	shellCommand = '''awk \'{printf(\"mode %1d %.2f %s\\n\", NR, sqrt(1/(0.0001*NR)), $0)}\' zeroEigvecsComplex.dat > zerotemp'''
		# 	os.system(shellCommand)
		# 	shellCommand = "cat temp1 zerotemp temp2 > "+path+self.pdbName+"_complex_zeros_ANM.nmd ; rm temp1 zerotemp temp2"
		# 	os.system(shellCommand)
		# 	# change the working path back to the directory of the program
		# 	os.chdir(os.path.dirname(os.path.abspath(__file__)))
		# except Exception, err:
		# 	print "Exception occurred when adding zero eigenvalue eigenvectors to the reference nmd file: ", err
		# 	print traceback.format_exc()                       
		
		# compress the nmd files
		try:
			os.chdir(path)
			filesToCombine = glob.glob("*.nmd")
			combineArguments = " ".join(filesToCombine)
			#shellCommand = "paste "+combineArguments+" | column -s $'\\t' -t > combined.txt"
			shellCommand = "tar -jcvf "+self.pdbName+"anms.nmd.tar.bz2 "+combineArguments
			os.system(shellCommand)
			# change the working path back to the directory of the program
			os.chdir(os.path.dirname(os.path.abspath(__file__)))
		except Exception, err:
			print "Exception occurred when gzipping nmd files: ", err
			print traceback.format_exc()
			
		# delete the nmd files
		try:
			os.chdir(path)
			shellCommand = "rm *.nmd"
			os.system(shellCommand)
			# change the working path back to the directory of the program
			os.chdir(os.path.dirname(os.path.abspath(__file__)))
		except Exception, err:
			print "Exception occurred when deleting nmd files: ", err
			print traceback.format_exc()
			
		if storeANMs:
			try: 
				saveModel(ANMs._anm_reference[0], filename=path+self.pdbName+"_reference_ANM", matrices=True)                
			except Exception, err:
				print "Exception occurred saving the reference ANM to a file: ", err
				print traceback.format_exc()
				
			try: 
				saveModel(ANMs._anm_counterpart[0], filename=path+self.pdbName+"_anm_counterpart", matrices=True)
			except Exception, err:
				print "Exception occurred when saving the counterpart ANM to a file: ", err
				print traceback.format_exc()
				
			try: 
				saveModel(ANMs._anm_complex[0], filename=path+self.pdbName+"_anm_complex", matrices=True)
			except Exception, err:
				print "Exception occurred when saving the complex ANM to a file: ", err
				print traceback.format_exc()
			# compress the npz files
			
			try:
				os.chdir(path)
				filesToCombine = glob.glob("*.npz")
				combineArguments = " ".join(filesToCombine)
				#shellCommand = "paste "+combineArguments+" | column -s $'\\t' -t > combined.txt"
				shellCommand = "tar -jcvf "+self.pdbName+"anms.npz.tar.bz2 "+combineArguments
				os.system(shellCommand)
				# change the working path back to the directory of the program
				os.chdir(os.path.dirname(os.path.abspath(__file__)))
			except Exception, err:
				print "Exception occurred when gzipping npz anm files: ", err
				print traceback.format_exc()
				
			# delete the npz files
			try:
				os.chdir(path)
				shellCommand = "rm *.npz"
				os.system(shellCommand)
				# change the working path back to the directory of the program
				os.chdir(os.path.dirname(os.path.abspath(__file__)))
			except Exception, err:
				print "Exception occurred when deleting npz files: ", err
				print traceback.format_exc()            
