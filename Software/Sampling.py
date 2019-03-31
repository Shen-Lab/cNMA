'''
Created on Mar 24, 2019

@author: Yue Cao

This program is for generating sampled complex conformations from complex normal modes. The new conformations contain both rigid-body motion of
the ligand and the flexbility of the complex.

Reference:
Yue Cao and Yang Shen, Bayesian active learning for optimization and uncertainty quantification in protein docking,
arXiv preprint arXiv:1902.00067, 2019
'''

import numpy as np
from prody import LOGGER
from prody.atomic import Atomic, AtomGroup
from prody.ensemble import Ensemble

class complex_sampling():

	def __init__(self, allresidue, allatom, nresir, atoms, rrmsd_pred):
		self.allresidue = allresidue
		self.allatom = allatom
		self.nresir = nresir
		self.atoms = atoms
		self.rrmsd_pred = rrmsd_pred[0][0]


	def Bset(self):
		'''
		Bset is for generating the set of normal modes used for sampling. Here B=I1+I2, where I1 is the top 6 from complex nornal modes and I2 is
		top 6 from receptor normal modes.

		eigvec is n*m, where n is the number of nodes, and m=3*nresic

		'''

		# top k comple normal modes:

		ind=np.arange(0, self.k1).tolist()

		#eigenvalues rescaled:
		rescale_eigval = np.divide(self.eigval , np.linalg.norm(self.eigvec.T[0: self.nresir*3].T, axis=1)**2 )

		#resort them:
		indrec = np.argsort(rescale_eigval)
	#	print len(indrec)

		for i in range(len(indrec)):
			if(indrec[i] not in ind):
				ind.append(indrec[i])
				if(len(ind) == self.k1 + self.k2):
					break

		'''
		ind is the final set of index of the normal modes used for sampling.
		'''
		print "normal modes used: ", ind

		return ind

	'''
	Calculating ||\sum_j \in B \frac{r_j}{\sqrt(\lambda_j)} ||^2
	'''
	def scalar(self, ind, coef, summa):

		s=0.
		for i in range(len(coef)):
			for j in range(len(coef)):
				s += coef[i]*coef[j]/np.sqrt(self.eigval[ind[i]]*self.eigval[ind[j]]) * summa[i][j]

		return s


	def rejectsampling(self, ind):

		sum_r_r=np.zeros((len(ind), len(ind)))
		sum_l_l=np.zeros((len(ind), len(ind)))



		for i in range(len(ind)):
			for j in range(len(ind)):
				sum_r_r[i][j]= np.dot(self.eigvec.T[0: self.nresir*3].T[ind[i]],  self.eigvec.T[0: self.nresir*3].T[ind[j]])
				sum_l_l[i][j]= np.dot(self.eigvec.T[self.nresir*3 : ].T[ind[i]],  self.eigvec.T[self.nresir*3 : ].T[ind[j]])
		

		step=0
		while(1):

			randn = np.random.normal(0, 1, len(ind))
			coef = randn / np.linalg.norm(randn)
			

			if(self.rrmsd == None):
				rrmsd = np.random.normal(0.99, 0.31) * self.rrmsd_pred
			else:
				rrmsd = self.rrmsd

			s=0.
			for i in range(len(ind)):
				s+= coef[i]**2/self.eigval[ind[i]]

			scalar = self.scalar(ind, coef, sum_r_r)
			
			# reject sampling finds a feasible solution and now lets break.
			if( self.nresir * rrmsd**2 * (s - scalar) <= self.nresil * self.th_lrmsd**2 * scalar ):
				break


			step+=1

			if(step>100000):
				self.rrmsd_pred/=2.
	
		#print coef* np.sqrt( rrmsd * self.nresir / scalar), rrmsd, self.rrmsd_pred, step, self.nresil, self.nresir
		#print coef* np.sqrt( rrmsd**2 * self.nresir / scalar)

		return coef* np.sqrt( rrmsd**2 * self.nresir / scalar), rrmsd, np.sqrt((s-scalar)* rrmsd**2 * self.nresir / scalar / self.nresil)





	def generate(self, n_confs=100, rrmsd=None, k1=9, k2=3, th_lrmsd=6):
		'''
		This program is for generating new conformations.
		'''
		self.n_confs = n_confs
		self.rrmsd=rrmsd
		self.k1 = k1
		self.k2 = k2
		self.th_lrmsd = th_lrmsd
		self.eigval = self.allresidue.getEigvals()[6:]
		self.eigvec = self.allresidue.getArray().T[6:]
		self.nresil = len(self.eigvec[0])/3 - self.nresir


		n_confs = int(n_confs)
		n_atoms = self.allatom.numAtoms()
		initial = self.atoms.getCoords()

		print "number of total atoms for complex: ", n_atoms

		ind = self.Bset()
		array = self.allatom._getArray().T[6:]
		confs = []

		'''
		#check
		x, Rrmsd, Lrmsd = self.rejectsampling(ind)

		coords =  np.zeros(3*(self.nresil+self.nresir))

		for j in range(2*k):
			coords +=  x[j]/np.sqrt(self.eigval[ind[j]]) * self.eigvec[ind[j]]

		print ("#%d Rrmsd=%.3f Lrmsd=%.10f" %(1+1, Rrmsd, Lrmsd))
		print np.sqrt( np.sum(coords[0:3*self.nresir]**2)/self.nresir ), np.sqrt( np.sum(coords[3*self.nresir:]**2)/self.nresil )


		print array[0][0:50], np.linalg.norm(array[0]), self.eigvec[0][0:10]
		#checkend
'''

		for i in range(n_confs):
			
			x, Rrmsd, Lrmsd = self.rejectsampling(ind)
			coords = np.zeros(3*n_atoms)			
			
			for j in range(k1+k2):
				
				coords +=  x[j]/np.sqrt(self.eigval[ind[j]]) * array[ind[j]]

			confs.append(coords.reshape((n_atoms, 3)))

			print ("#%d Rrmsd=%.3f Lrmsd=%.3f" %(i+1, Rrmsd, Lrmsd))


		ensemble = Ensemble('Conformations along {0}'.format(self.allatom))

		ensemble.setCoords(initial)

		ensemble.addCoordset(np.array(confs) + initial)

		return ensemble

	
		
    	#print array1.ndim, len(array1), len(array1[0])

