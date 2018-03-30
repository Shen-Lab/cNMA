import numpy as np
import pandas as pd
import os
from sklearn.kernel_ridge import KernelRidge
from sklearn.externals import joblib




class RMSDPrediction(object):


	def __init__(self,protein1_A, protein2_A=None):
		'''
		Constructor
		'''
		self.path = os.path.dirname(__file__) + "/"
		os.environ['paramsPath'] = self.path+'params.txt'
		#self.index = os.popen("if [ -e $paramsPath ]; then d=true; else d=false; fi; echo -n $d").read()
		self.index = os.path.exists("$paramsPath")
		self.protein1_A = protein1_A
		# print paramsPath
		# print self.index


	def modelinput(self):
		''' 
		Input model trained by 2c ZDOCK dataset through KRR for further prediction 
		'''
		if self.index == True:
			self.params = np.loadtxt(self.path+"params.txt", dtype = str)       
			filename = self.params[3] + '.joblib.pkl'
		elif self.index == False:
			filename = 'KRR.joblib.pkl'
		ModelParameters = joblib.load(self.path + filename)
		return ModelParameters

	def eigeninput(self, path):
		'''
		Input eigenvalues file
		'''
		filename = pd.read_csv(path + 'eigenvaluesReference.txt')
		eigenvalues = filename.values
		return eigenvalues

	def calculation(self, eigenvalues, title):
		'''
		calculate top k modes and build feature matrix
		'''
		mult = 1
		summ = 0
		

		phi = np.zeros([1, 8], float)
		eig_len = len(eigenvalues)
		n = np.zeros([eig_len, 2], float)

		#read size

		os.environ['pdbPath'] = self.protein1_A
		os.system("grep ' CA ' $pdbPath | wc -l > Nres.txt")
		size = np.loadtxt(self.path + "Nres.txt")

		#size cutoff
		size_cutoff_index = np.zeros([100], int)

		#rigidity cutoff
		threshold = 250  #threshold range
		rig_cutoff_index = np.zeros([threshold], int)
		count = 0
		cut_index = np.zeros([threshold])
		l_1 = eigenvalues[0]

		#size cutoff
		for i in range(0,100):
			size_cutoff_index[i] = round((3 * size - 6) * i * 0.01)

		for count in range(0, eig_len):
			#factorial
			mult = mult * (1 / float(np.sqrt(eigenvalues[count])))
			n[count, 0] = np.power(mult,(1 / float(count + 1)))


			#summation
			summ = summ + 1 / float(eigenvalues[count])
			n[count, 1] = np.sqrt(summ / (count + 1))

			#rigidity cutoff
			for i in range(0, threshold):
				 if (eigenvalues[count] / l_1 > (i + 1) * 10) and (cut_index[i] == 0):
					rig_cutoff_index[i] = count
					cut_index[i] = 1
		
		pass


		#absolute cutoff feature
		ab_cutoff = 100
		phi[0, 0] = n[ab_cutoff, 1]
		phi[0, 3] = n[ab_cutoff, 0]

		#size_related cutoff feature
		percentage = 100
		X_summ = []
		X_mult = []

		for i in range(0, percentage):
			if size_cutoff_index[i] < 400:
				X_mult.append(n[size_cutoff_index[i], 0])
				X_summ.append(n[size_cutoff_index[i], 1])
			else:
				X_mult.append(n[399, 0])
				X_summ.append(n[399, 1])
			pass

			# size-related cutoff feature (percentage varies)
			if i == 60: #various
				phi[0, 1] = X_summ[i]
				phi[0, 4] = X_mult[i]


		#rigidity_related cutoff feature
		cutoff = 100 #various
		X_summ = 0
		X_mult = 0

		if rig_cutoff_index[cutoff] == 0:
			X_mult = n[399, 0]
			X_summ = n[399, 1]
		else:
			X_mult = n[rig_cutoff_index[cutoff], 0]
			X_summ = n[rig_cutoff_index[cutoff], 1]
		pass

		phi[0, 2] = X_summ
		phi[0, 5] = X_mult

		#size of Nres feature
		phi[0, 6] = size

		#group the features into three sets(phi_mult, phi_summ)
		phi_mult = phi[:, 3:7]
		# phi_summ = phi[0,1,2,6]
		phi_mean = np.loadtxt(self.path + '10AA_phi_mult_mean.txt')
		phi_std = np.loadtxt(self.path + '10AA_phi_mult_std.txt')
		#standardize mult feature
		for i in range(0, phi_mult.shape[1]):
			phi_mult[0, i] = (phi_mult[0, i] - phi_mean[i]) / phi_std[i]
		return phi_mult

	def prediction(self, path, params, feature_matrix):
		RMSD_mean = np.loadtxt(self.path + '10AA_RMSD_mean.txt')
		
		RMSD_predited = params.predict(feature_matrix) + RMSD_mean

		return RMSD_predited

