import random
import matplotlib.pyplot as pyp
import numpy as np
import math
import model
import copy
import pickle
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

class stoch_model:
	def __init__(self, runs=10, sensitivity=False, change=0.1):
		self.vera_model = model.model()
		self.rates = copy.deepcopy(self.vera_model.rates)
		self.probabilities = copy.deepcopy(self.vera_model.probability)
		self.probabilities_G1 = copy.deepcopy(self.vera_model.probability_G1)
		self.probabilities_S_M = copy.deepcopy(self.vera_model.probability_S_M)
		self.comp_weights = copy.deepcopy(self.vera_model.compartment_weights)
		self.fa_weights = copy.deepcopy(self.vera_model.weights_fa)

	
		#calculation fatty acids
		self.fatty_acid_distribution = {}
		self.fatty_acid_std = {}
		self.fatty_acids = ['C16:0', 'C16:1', 'C18:0', 'C18:1']
		self.fa = ['saturated', 'unsaturated']

		#membrane compositions
		self.membrane_comp_mean = {}
		self.membrane_comp_std = {}
		self.lipids = ['CL', 'ES', 'PA', 'PC', 'PE', 'PI', 'PS', 'SE', 'SL', 'TAG']
		self.membranes_length = {}
		self.membranes_length_std = {}

		#membrane lengths
		self.plasma_membrane = []
		self.secretory_vesicles = []
		self.vacuoles = []
		self.nucleus = []
		self.peroxisomes = []
		self.light_microsomes = []
		self.inner_mit_membrane = []
		self.outer_mit_membrane = []
		self.lipid_droplets = []
		self.mean_time_lists = {}
		self.std_time_lists = {}
		self.membrane_lists = {'plasma_membrane': self.plasma_membrane, 'secretory_vesicles': self.secretory_vesicles, 'vacuoles': self.vacuoles, \
								'nucleus': self.nucleus, 'peroxisomes': self.peroxisomes, 'light_microsomes': self.light_microsomes,\
								'inner_mit_membrane': self.inner_mit_membrane, 'outer_mit_membrane': self.outer_mit_membrane, 'lipid_droplets': self.lipid_droplets}
		
		self.membranes_list = ['plasma_membrane', 'secretory_vesicles', 'vacuoles', 'nucleus', 'peroxisomes', 'light_microsomes', 'inner_mit_membrane', 'outer_mit_membrane', 'lipid_droplets']
		
		self.time = []

	def run(self, runs, sensitivity, change):
		
		self.time = [i for i in range(7200)]

		self.fatty_acid_distribution = pickle.load(open("./sens_analysis_start_11-12-15/fa_distr_mean.pkl", "rb"))
		self.fatty_acid_std = pickle.load(open("./sens_analysis_start_11-12-15/fa_distr_std.pkl", "rb"))
		self.mean_time_lists = pickle.load(open("./sens_analysis_start_11-12-15/membranes_length_mean.pkl", "rb"))
		self.std_time_lists = pickle.load(open("./sens_analysis_start_11-12-15/membranes_length_std.pkl", "rb"))
		self.membrane_comp_mean = pickle.load(open("./sens_analysis_start_11-12-15/membrane_comp_mean.pkl", "rb"))
		self.membrane_comp_std = pickle.load(open("./sens_analysis_start_11-12-15/membrane_comp_std.pkl", "rb"))
		#plot heatmaps of sesitivity analysis	
		self.heatmap_fatty_acids()
		#self.heatmap_membrane_comp()
		self.heatmap_membrane_length()


	def heatmap_fatty_acids(self):
		#heatmap of fatty acid distributions: x = fatty acids, y = reaction (4 heatmaps: rates +/- 10%, probabilities +/- 10%)
		heatmaps_change = ['rate+10', 'rate-10', 'probability+10', 'probability-10', 'compartment_weights+10', 'compartment_weights-10', 'fatty_acid_weights+10', 'fatty_acid_weights-10']
		for change in heatmaps_change:
			if 'rate' in change: 
				column_labels = ['Wildtype']
				for par in self.rates:
					column_labels.append(par)
					v_max = 0.015
					v_min = -0.015

			elif 'probability' in change: 
				column_labels = ['Wildtype']
				for par in self.probabilities:
					column_labels.append(par)
					v_max = 0.01
					v_min = -0.01

			elif 'compartment' in change:
				column_labels = self.membranes_list[:-1]
				v_max = 0.0015
				v_min = -0.0015

			elif 'fatty' in change:
				column_labels = self.fa
				v_max = 0.1
				v_min = -0.1

			row_labels = self.fatty_acids
			self.heatmap_fa = []

			i = 0
			for reaction in column_labels:
				if reaction == 'Wildtype':
					self.heatmap_fa.append(i)
					difference = []
					for fa in row_labels:
						diff = math.log(1, 2)
						difference.append(diff)
					self.heatmap_fa[i] = difference
					i += 1

				else:
					self.heatmap_fa.append(i)
					difference = []
					for fa in row_labels:
						if self.fatty_acid_distribution[reaction + change][fa] > 0:
							diff = math.log((self.fatty_acid_distribution[reaction + change][fa] / self.fatty_acid_distribution['Wildtype'][fa]), 2)
						else:
							diff = math.log(0.00001, 2)
						difference.append(diff)
					self.heatmap_fa[i] = difference
					i += 1

			self.heatmap_fa = np.asarray(self.heatmap_fa)
			fig, ax = pyp.subplots()
			general_heatmap = ax.pcolor(self.heatmap_fa, cmap = pyp.cm.seismic, vmax = v_max, vmin = v_min)

			ax.set_xticks(np.arange(self.heatmap_fa.shape[1])+0.5, minor = False)
			ax.set_yticks(np.arange(self.heatmap_fa.shape[0])+0.5, minor = False)

			ax.invert_yaxis()
			ax.xaxis.tick_top()

			ax.set_xticklabels(row_labels, minor = False)
			ax.set_yticklabels(column_labels, minor = False)			

			cbar = pyp.colorbar(general_heatmap)
			cbar.ax.set_ylabel('log2 change', rotation = 270, labelpad = 22)

			ax.set_title('Fatty acids distribution: ' + change, y = 1.07)
			ax.axis('tight')
			
			pyp.show()

	def heatmap_membrane_comp(self):
		#heatmap of membrane compositions: x = lipids, y = reaction (4 heatmaps for each membrane, 4 x 9 heatmaps)
		heatmaps_change = ['rate+10', 'rate-10', 'probability+10', 'probability-10', 'compartment_weights+10', 'compartment_weights-10', 'fatty_acid_weights+10', 'fatty_acid_weights-10']

		for membrane in self.membranes_list:
			for change in heatmaps_change:

				column_labels = ['Wildtype']
				if 'rate' in change:
					for par in self.rates:
						column_labels.append(par)
				elif 'probability' in change:
					for par in self.probabilities:
						column_labels.append(par)
				elif 'compartment' in change:
					column_labels = self.membranes_list[:-1]
				elif 'fatty' in change:
					column_labels = self.fa

				self.heatmap_membrane = []
				i = 0
				for reaction in column_labels:				
					self.heatmap_membrane.append(i)
					difference = []
					for lipid in self.lipids:
						if self.membrane_comp_mean['Wildtype'][membrane][lipid] == 0:
							diff = math.log(1, 2)
						else:
							if reaction == 'Wildtype':
								diff = math.log(1, 2)
							elif self.membrane_comp_mean[reaction + change][membrane][lipid] > 0:
								diff = math.log((self.membrane_comp_mean[reaction + change][membrane][lipid] / self.membrane_comp_mean['Wildtype'][membrane][lipid]), 2)
							else:
								diff = math.log(0.00001, 2)
						difference.append(diff)
					self.heatmap_membrane[i] = difference
					i += 1

				self.heatmap_membrane = np.asarray(self.heatmap_membrane)
				fig, ax = pyp.subplots()
				general_heatmap = ax.pcolor(self.heatmap_membrane, cmap = pyp.cm.Blues)#, vmax = 1.6, vmin = -1.6)

				ax.set_xticks(np.arange(self.heatmap_membrane.shape[1])+0.5, minor = False)
				ax.set_yticks(np.arange(self.heatmap_membrane.shape[0])+0.5, minor = False)

				ax.invert_yaxis()
				ax.xaxis.tick_top()
				ax.set_xticklabels(self.lipids, minor = False)
				column_labels_short = [label.split(' ', 1)[0] for label in column_labels]
				ax.set_yticklabels(column_labels_short, minor = False)

				cbar = pyp.colorbar(general_heatmap)
				cbar.ax.set_ylabel('log2 change', rotation = 270, labelpad = 22)

				ax.set_title(membrane.replace('_', ' ').title() + ' composition: ' + change, y = 1.1)
				ax.axis('tight')
				
				pyp.show()


	def heatmap_membrane_length(self):
		#heatmap of membrane lengths: x = membrane, y = reaction (4 heatmaps: rates +/- 10%, probabilities +/- 10%)
		heatmaps_change = ['rate+10', 'rate-10', 'probability+10', 'probability-10', 'compartment_weights+10', 'compartment_weights-10', 'fatty_acid_weights+10', 'fatty_acid_weights-10']
		row_labels = self.membranes_list[:-1]

		for change in heatmaps_change:
			if 'rate' in change:
				column_labels = ['Wildtype']
				for par in self.rates:
					column_labels.append(par)
			elif 'probability' in change:
				column_labels = ['Wildtype']
				for par in self.probabilities:
					column_labels.append(par)
			elif 'compartment' in change:
				column_labels = self.membranes_list[:-1]
			elif 'fatty' in change:
				column_labels = self.fa

			self.heatmap_mem_length = []

			i = 0
			for reaction in column_labels:
				self.heatmap_mem_length.append(i)
				difference = []
				for membrane in row_labels:
					if reaction == 'Wildtype':
						diff = math.log(1, 2)
					elif self.mean_time_lists[reaction + change][membrane][-1] > 0:
						diff = math.log((self.mean_time_lists[reaction + change][membrane][-1] / self.mean_time_lists['Wildtype'][membrane][-1]), 2)
					else:
						diff = math.log(0.00001, 2)
					difference.append(diff)
				self.heatmap_mem_length[i] = difference
				i += 1

			self.heatmap_mem_length = np.asarray(self.heatmap_mem_length)
			fig, ax = pyp.subplots()
			general_heatmap = ax.pcolor(self.heatmap_mem_length, cmap = pyp.cm.Blues)#, vmax = 0.1, vmin = -0.1)

			ax.set_xticks(np.arange(self.heatmap_mem_length.shape[1])+0.5, minor = False)
			ax.set_yticks(np.arange(self.heatmap_mem_length.shape[0])+0.5, minor = False)

			ax.invert_yaxis()
			ax.xaxis.tick_top()

			ax.set_xticklabels(row_labels, minor = False)
			ax.set_yticklabels(column_labels, minor = False)			

			cbar = pyp.colorbar(general_heatmap)
			cbar.ax.set_ylabel('log2 change', rotation=270, labelpad = 22)

			ax.set_title('Membrane lengths: ' + change, y = 1.4)

			pyp.xticks(rotation = 70)
			ax.axis('tight')
			
			pyp.show()