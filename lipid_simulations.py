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
		# load model and parameters
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
		# simulate
		# sensitivity analysis
		if sensitivity:
			self.results = {}
			self.sensitivity_analysis(runs, change)
			#self.sensitivity_statistics()
		#plot heatmaps of sesitivity analysis	
			#self.heatmap_fatty_acids()
			#self.heatmap_membrane_comp()
			#self.heatmap_membrane_length()

		#no sensitivity analysis
		else:
			self.results = []
			self.run_model(runs)
			self.calc_statistics()		

		# plot results
			self.plot_results_fatty_acids()
			self.plot_results_membranes()
			self.plot_membrane_comp()

	def sensitivity_analysis(self, runs, change):
		#Wildtype data first (run without changes)
		self.results['Wildtype'] = []
		for run in range(runs):
			self.results['Wildtype'].append(self.vera_model.run())
			#if run > 0 and (run+1) % 10 == 0:
		pickle.dump(self.results['Wildtype'], open("./sens_analysis_start_10-12-15/Wildtype.pkl", "wb"))
		
		#sensitivity of rates (+10%)
		for par in self.rates:
			self.vera_model.rates[par] = int(round(self.vera_model.rates[par] + self.vera_model.rates[par] * change))
			self.results[par + 'rate+10'] = []
			for run in range(runs):
				self.results[par + 'rate+10'].append(self.vera_model.run())
				#if run > 0 and (run+1) % 10 == 0:
			pickle.dump(self.results[par + 'rate+10'], open("./sens_analysis_start_10-12-15/"+par+"rate+10.pkl", "wb"))
			self.vera_model.rates[par] = self.rates[par]
		'''
		#sensitivits of rates (-10%)
		for par in self.rates:
			self.vera_model.rates[par] = int(round(self.vera_model.rates[par] - self.vera_model.rates[par] * change))
			self.results[par + 'rate-10'] = []
			for run in range(runs):
				self.results[par + 'rate-10'].append(self.vera_model.run())
				#if run > 0 and (run+1) % 10 == 0:
			pickle.dump(self.results[par + 'rate-10'], open("./sens_analysis_start_10-12-15/"+par+"rate-10.pkl", "wb"))
			self.vera_model.rates[par] = self.rates[par]
		
		#sensitivity of probabilities (+10%)
		for prob in self.probabilities:
			self.vera_model.probability[prob] = self.vera_model.probability[prob] + self.vera_model.probability[prob] * change
			if prob in self.vera_model.probability_G1.keys():
				self.vera_model.probability_G1[prob] = self.vera_model.probability_G1[prob] + self.vera_model.probability_G1[prob] * change
			if prob in self.vera_model.probability_S_M.keys():
				self.vera_model.probability_S_M[prob] = self.vera_model.probability_S_M[prob] + self.vera_model.probability_S_M[prob] * change
			self.results[prob + 'probability+10'] = []
			for run in range(runs):
				self.results[prob + 'probability+10'].append(self.vera_model.run())
				#if run > 0 and (run+1) % 10 == 0:
			pickle.dump(self.results[prob + 'probability+10'], open("./sens_analysis_start_10-12-15/"+prob+"probability+10.pkl", "wb"))
			self.vera_model.probability[prob] = self.probabilities[prob]
			if prob in self.vera_model.probability_G1.keys():
				self.vera_model.probability_G1[prob] = self.probabilities_G1[prob]
			if prob in self.vera_model.probability_S_M.keys():
				self.vera_model.probability_S_M[prob] = self.probabilities_S_M[prob]
		
		#sensitivity of probabilities (-10%)
		for prob in self.probabilities:
			self.vera_model.probability[prob] = self.vera_model.probability[prob] - self.vera_model.probability[prob] * change
			if prob in self.vera_model.probability_G1.keys():
				self.vera_model.probability_G1[prob] = self.vera_model.probability_G1[prob] - self.vera_model.probability_G1[prob] * change
			if prob in self.vera_model.probability_S_M.keys():
				self.vera_model.probability_S_M[prob] = self.vera_model.probability_S_M[prob] - self.vera_model.probability_S_M[prob] * change
			self.results[prob + 'probability-10'] = []
			for run in range(runs):
				self.results[prob + 'probability-10'].append(self.vera_model.run())
				#if run > 0 and (run+1) % 10 == 0:
			pickle.dump(self.results[prob + 'probability-10'], open("./sens_analysis_start_10-12-15/"+prob+"probability-10.pkl", "wb"))
			self.vera_model.probability[prob] = self.probabilities[prob]	
			if prob in self.vera_model.probability_G1.keys():
				self.vera_model.probability_G1[prob] = self.probabilities_G1[prob]
			if prob in self.vera_model.probability_S_M.keys():
				self.vera_model.probability_S_M[prob] = self.probabilities_S_M[prob]
		
		#sensitivity of comp weights (+10%)
		for weight in range(len(self.comp_weights)):
			self.vera_model.compartment_weights[weight] = self.vera_model.compartment_weights[weight] + self.vera_model.compartment_weights[weight] * change
			self.results[self.membranes_list[weight] + 'compartment_weights+10'] = []
			for run in range(runs):
				self.results[self.membranes_list[weight] + 'compartment_weights+10'].append(self.vera_model.run())
				#if run > 0 and (run+1) % 10 == 0:
			pickle.dump(self.results[self.membranes_list[weight] + 'compartment_weights+10'], open("./sens_analysis_start_10-12-15/"+self.membranes_list[weight]+"compartment_weights+10.pkl", "wb"))
			self.vera_model.compartment_weights[weight] = self.comp_weights[weight]	
		
		#sensitivity of comp weights (-10%)
		for weight in range(len(self.comp_weights)):
			self.vera_model.compartment_weights[weight] = self.vera_model.compartment_weights[weight] - self.vera_model.compartment_weights[weight] * change
			self.results[self.membranes_list[weight] + 'compartment_weights-10'] = []
			for run in range(runs):
				self.results[self.membranes_list[weight] + 'compartment_weights-10'].append(self.vera_model.run())
				#if run > 0 and (run+1) % 10 == 0:
				pickle.dump(self.results[self.membranes_list[weight] + 'compartment_weights-10'], open("./sens_analysis_start_10-12-15/"+self.membranes_list[weight]+"compartment_weights-10.pkl", "wb"))
			self.vera_model.compartment_weights[weight] = self.comp_weights[weight]	
		
		#sensitivity of fa weights (+10%)
		for weight in range(len(self.fa_weights)):
			if weight == 0:
				weight2 = 1
			else:
				weight2 = 0
			self.vera_model.weights_fa[weight2] = self.vera_model.weights_fa[weight2] - self.vera_model.weights_fa[weight] * change
			self.vera_model.weights_fa[weight] = self.vera_model.weights_fa[weight] + self.vera_model.weights_fa[weight] * change
			self.results[self.fa[weight] + 'fatty_acid_weights+10'] = []
			for run in range(runs):
				self.results[self.fa[weight] + 'fatty_acid_weights+10'].append(self.vera_model.run())
				#if run > 0 and (run+1) % 10 == 0:
			pickle.dump(self.results[self.fa[weight] + 'fatty_acid_weights+10'], open("./sens_analysis_start_10-12-15/"+self.fa[weight]+"fatty_acid_weights+10.pkl", "wb"))
			self.vera_model.weights_fa[weight] = self.fa_weights[weight]
			self.vera_model.weights_fa[weight2]	= self.fa_weights[weight2]
		
		#sensitivity of fa weights (-10%)
		for weight in range(len(self.fa_weights)):
			if weight == 0:
				weight2 = 1
			else:
				weight2 = 0
			self.vera_model.weights_fa[weight2] = self.vera_model.weights_fa[weight2] + self.vera_model.weights_fa[weight] * change
			self.vera_model.weights_fa[weight] = self.vera_model.weights_fa[weight] - self.vera_model.weights_fa[weight] * change
			self.results[self.fa[weight] + 'fatty_acid_weights-10'] = []
			for run in range(runs):
				self.results[self.fa[weight] + 'fatty_acid_weights-10'].append(self.vera_model.run())
				#if run > 0 and (run+1) % 10 == 0:
			pickle.dump(self.results[self.fa[weight] + 'fatty_acid_weights-10'], open("./sens_analysis_start_10-12-15/"+self.fa[weight]+"fatty_acid_weights-10.pkl", "wb"))
			self.vera_model.weights_fa[weight] = self.fa_weights[weight]
			self.vera_model.weights_fa[weight2]	= self.fa_weights[weight2]
		'''
		#time list imported from model
		self.time = self.vera_model.t
		pickle.dump(self.results.keys(), open("./sens_analysis_start_10-12-15/results_sens_keys.pkl", "wb"))
			

	def sensitivity_statistics(self):
		#fatty acids distribution (mean and std for each change)
		for reaction in self.results:
			self.fatty_acid_distribution[reaction] = {}
			self.fatty_acid_std[reaction] = {}
			for fa in self.fatty_acids:
				fa_relatives = []
				for run in range(len(self.results[reaction])):
					fa_relatives.append(self.results[reaction][run][0][fa])
					self.fatty_acid_distribution[reaction][fa] = np.mean(fa_relatives)
					self.fatty_acid_std[reaction][fa] = np.std(fa_relatives)


		#mean and std membrane lengths after each time step for each change
		for reaction in self.results:
			self.mean_time_lists[reaction] = {}
			self.std_time_lists[reaction] = {}
			for run in range(len(self.results[reaction])):
				self.plasma_membrane.append(self.results[reaction][run][1][0])
				self.secretory_vesicles.append(self.results[reaction][run][1][1])
				self.vacuoles.append(self.results[reaction][run][1][2])
				self.nucleus.append(self.results[reaction][run][1][3])
				self.peroxisomes.append(self.results[reaction][run][1][4])
				self.light_microsomes.append(self.results[reaction][run][1][5])
				self.inner_mit_membrane.append(self.results[reaction][run][1][6])
				self.outer_mit_membrane.append(self.results[reaction][run][1][7])
				self.lipid_droplets.append(self.results[reaction][run][1][8])

			z = 0
			for membrane in self.membrane_lists.values():
				self.mean_time_lists[reaction][self.membrane_lists.keys()[z]] = []
				self.std_time_lists[reaction][self.membrane_lists.keys()[z]] = []
				for i in range(len(self.time)):
					membrane_time = []
					for j in range(len(membrane)):
						membrane_time.append(membrane[j][i])
					mean_time = np.mean(membrane_time)
					std_time = np.std(membrane_time)
					self.mean_time_lists[reaction][self.membrane_lists.keys()[z]].append(mean_time)
					self.std_time_lists[reaction][self.membrane_lists.keys()[z]].append(std_time)
				z += 1


		#mean and std membrane compositions for each change
		for reaction in self.results:
			self.membrane_comp_mean[reaction] = {}
			self.membrane_comp_std[reaction] = {}
			for membrane in self.membrane_lists.keys():
				self.membrane_comp_mean[reaction][membrane] = {}
				self.membrane_comp_std[reaction][membrane] = {}
				for lipid in self.lipids:
					lipid_relatives = []
					for runs in range(len(self.results[reaction])):
						lipid_relatives.append(self.results[reaction][runs][2][membrane][lipid])
					mean_lipid_relatives = np.mean(lipid_relatives)
					std_lipid_relatives = np.std(lipid_relatives)
					self.membrane_comp_mean[reaction][membrane][lipid] = mean_lipid_relatives
					self.membrane_comp_std[reaction][membrane][lipid] = std_lipid_relatives

	def heatmap_fatty_acids(self):
		#heatmap of fatty acid distributions: x = fatty acids, y = reaction (4 heatmaps: rates +/- 10%, probabilities +/- 10%)
		heatmaps_change = [' rate + 10 %', ' rate - 10 %', ' probability + 10 %', ' probability - 10 %', ' compartment weights + 10 %', ' compartment weights - 10 %', ' fatty acid weights + 10 %', ' fatty acid weights - 10 %']
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

			elif 'fatty acid' in change:
				column_labels = self.fa

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
			general_heatmap = ax.pcolor(self.heatmap_fa, cmap = pyp.cm.Blues)#, vmax = 0.2, vmin = -0.2)

			ax.set_xticks(np.arange(self.heatmap_fa.shape[1])+0.5, minor = False)
			ax.set_yticks(np.arange(self.heatmap_fa.shape[0])+0.5, minor = False)

			ax.invert_yaxis()
			ax.xaxis.tick_top()

			ax.set_xticklabels(row_labels, minor = False)
			ax.set_yticklabels(column_labels, minor = False)			

			cbar = pyp.colorbar(general_heatmap)
			cbar.ax.set_ylabel('log2 change', rotation = 270, labelpad = 22)

			ax.set_title('Fatty acids distribution:' + change, y = 1.07)
			ax.axis('tight')
			
			pyp.show()

	def heatmap_membrane_comp(self):
		#heatmap of membrane compositions: x = lipids, y = reaction (4 heatmaps for each membrane, 4 x 9 heatmaps)
		heatmaps_change = [' rate + 10 %', ' rate - 10 %', ' probability + 10 %', ' probability - 10 %', ' compartment weights + 10 %', ' compartment weights - 10 %', ' fatty acid weights + 10 %', ' fatty acid weights - 10 %']

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
					column_labels = self.membranes_list
				elif 'fatty acid' in change:
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

				ax.set_title(membrane.replace('_', ' ').title() + ' composition:' + change, y = 1.1)
				ax.axis('tight')
				
				pyp.show()


	def heatmap_membrane_length(self):
		#heatmap of membrane lengths: x = membrane, y = reaction (4 heatmaps: rates +/- 10%, probabilities +/- 10%)
		heatmaps_change = [' rate + 10 %', ' rate - 10 %', ' probability + 10 %', ' probability - 10 %', ' compartment weights + 10 %', ' compartment weights - 10 %', ' fatty acid weights + 10 %', ' fatty acid weights - 10 %']
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
				column_labels = self.membranes_list
			elif 'fatty acid' in change:
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

			ax.set_title('Membrane lengths:' + change, y = 1.4)

			pyp.xticks(rotation = 70)
			ax.axis('tight')
			
			pyp.show()

	def run_model(self, runs):
		self.results = []
		for run in range(runs):
			self.results.append(self.vera_model.run())
		pickle.dump(self.results, open("new_final_1000runs.pkl", "wb"))
		self.time = self.vera_model.t

	def calc_statistics(self):		

		#mean and std fatty acid distribution
		for fa in self.fatty_acids:
			fa_relatives = []
			for run in range(len(self.results)):
				fa_relatives.append(self.results[run][0][fa])
				self.fatty_acid_distribution[fa] = np.mean(fa_relatives)
				self.fatty_acid_std[fa] = np.std(fa_relatives)

		#mean and std membrane lengths
		for run in range(len(self.results)):
			self.plasma_membrane.append(self.results[run][1][0])
			self.secretory_vesicles.append(self.results[run][1][1])
			self.vacuoles.append(self.results[run][1][2])
			self.nucleus.append(self.results[run][1][3])
			self.peroxisomes.append(self.results[run][1][4])
			self.light_microsomes.append(self.results[run][1][5])
			self.inner_mit_membrane.append(self.results[run][1][6])
			self.outer_mit_membrane.append(self.results[run][1][7])
			self.lipid_droplets.append(self.results[run][1][8])

		z = 0
		for membrane in self.membrane_lists.values():
			self.mean_time_lists[self.membrane_lists.keys()[z]] = []
			self.std_time_lists[self.membrane_lists.keys()[z]] = []
			for i in range(len(self.time)):
				membrane_time = []
				for j in range(len(membrane)):
					membrane_time.append(membrane[j][i])
				mean_time = np.mean(membrane_time)
				std_time = np.std(membrane_time)
				self.mean_time_lists[self.membrane_lists.keys()[z]].append(mean_time)
				self.std_time_lists[self.membrane_lists.keys()[z]].append(std_time)
			z += 1

		for membrane in self.membrane_lists.keys():
			self.membranes_length[membrane] = self.mean_time_lists[membrane][-1]
			self.membranes_length_std[membrane] = self.std_time_lists[membrane][-1]


		#mean and std membrane compositions		
		for membrane in self.membrane_lists.keys():
			self.membrane_comp_mean[membrane] = {}
			self.membrane_comp_std[membrane] = {}
			for lipid in self.lipids:
				lipid_relatives = []
				for runs in range(len(self.results)):
					lipid_relatives.append(self.results[runs][2][membrane][lipid])
				mean_lipid_relatives = np.mean(lipid_relatives)
				std_lipid_relatives = np.std(lipid_relatives)
				self.membrane_comp_mean[membrane][lipid] = mean_lipid_relatives
				self.membrane_comp_std[membrane][lipid] = std_lipid_relatives

			
	def plot_results_fatty_acids(self):
		#histogram of fatty acid distributions: x = fatty acids, y = mean value of all runs
		y = []
		y2 = [50, 10, 30, 10]
		x_names = self.fatty_acid_distribution.keys()
		ind = np.arange(len(x_names))
		width = 0.25
		yerrs = []
		for i in x_names:
			y.append((self.fatty_acid_distribution[i])*100)
			yerrs.append((self.fatty_acid_std[i])*100)
		
		fig, ax = pyp.subplots()
		fig_fa_distr = ax.bar(ind + 0.1, y, width, color = 'b', yerr = yerrs)
		fig_fa_distr2 = ax.bar(ind + width + 0.1, y2, width, color = 'y')
		ax.set_ylabel('%')
		ax.set_title('Fatty acids distribution')
		ax.set_xticks(ind + width + 0.1)
		ax.set_xticklabels(x_names)

		ax.legend((fig_fa_distr[0], fig_fa_distr2[0]), ('Model', 'Klug (2001)'))
		

		pyp.show()


	def plot_results_membranes(self):
		#graph of membrane growth: x = time, y = membrane length
		fig, ax = pyp.subplots()
		for membrane in self.membrane_lists.keys():
			ax.plot(self.time, self.mean_time_lists[membrane], label = membrane, linewidth = 2.0)
		ax.legend([membrane.replace('_', ' ') for membrane in self.membrane_lists.keys()], loc = 'center right', prop = {'size': 11})
		ax.set_title("Membrane growth")
		ax.set_ylabel('# of lipids')
		ax.set_xlabel('time in s')
		pyp.show()

		#histogram of the membrane lengths after the run: x = membrane, y = mean value membrane lengths
		x_names_long = self.membrane_lists.keys()
		x_names = [x_names_long[i].replace('_', ' ') for i in range(len(x_names_long))]
		ind = np.arange(len(x_names))
		width = 0.35
		y = []
		std = []
		for membrane in x_names_long:
			y.append(self.membranes_length[membrane])
			std.append(self.membranes_length_std[membrane])
		fig, ax = pyp.subplots()
		fig_membrane_bars = pyp.bar(ind, y, width, color = 'b', yerr = std)
		ax.set_xticks(ind + (width/2))
		ax.set_xticklabels(x_names)
		ax.set_title("Membrane length after 7200 s", y = 1.03)
		pyp.xticks(rotation = 60)
		pyp.show()

	
	def plot_membrane_comp(self):
		#histogram of the membrane compositions: x = lipids, y = percentages (a histogram for each membrane)
		zinser = {'plasma_membrane': {'PS': 0.17320, 'PI': 0.09124, 'PC': 0.08660, 'PE': 0.10464, 'CL': 0.00103, 'PA': 0.02010, 'ES': 0.48454, 'SE': 0.0, 'TAG': 0.0, 'SL': 0.03557},\
				'secretory_vesicles': {'PS': 0.08205, 'PI': 0.11745, 'PC': 0.20824, 'PE': 0.13573, 'CL': 0.01239, 'PA': 0.01525, 'ES': 0.42900, 'SE': 0.0, 'TAG': 0.0, 'SL': 0.05029},\
				'vacuoles': {'PS': 0.04817, 'PI': 0.16604, 'PC': 0.40517, 'PE': 0.17537, 'CL': 0.02442, 'PA': 0.02866, 'ES': 0.15200, 'SE': 0.0, 'TAG': 0.0, 'SL': 0.06525},\
				'nucleus': {'PS': 0.04038, 'PI': 0.09650, 'PC': 0.27645, 'PE': 0.16848, 'CL': 0.01049, 'PA': 0.01781, 'ES': 0.390, 'SE': 0.0, 'TAG': 0.0, 'SL': 0.02622},\
				'peroxisomes': {'PS': 0.03235, 'PI': 0.11360, 'PC': 0.34656, 'PE': 0.16465, 'CL': 0.05033, 'PA': 0.01150, 'ES': 0.281, 'SE': 0.0, 'TAG': 0.0, 'SL': 0.0},\
				'light_microsomes': {'PS': 0.05304, 'PI': 0.06019, 'PC': 0.40796, 'PE': 0.26583, 'CL': 0.00381, 'PA': 0.00222, 'ES': 0.206, 'SE': 0.0, 'TAG': 0.0, 'SL': 0.00397},\
				'inner_mit_membrane': {'PS': 0.02880, 'PI': 0.12273, 'PC': 0.29107, 'PE': 0.18192, 'CL': 0.12204, 'PA': 0.01137, 'ES': 0.242, 'SE': 0.0, 'TAG': 0.0, 'SL': 0.0},\
				'outer_mit_membrane': {'PS': 0.01189, 'PI': 0.10108, 'PC': 0.45190, 'PE': 0.32307, 'CL': 0.05847, 'PA': 0.04360, 'ES': 0.009, 'SE': 0.0, 'TAG': 0.0, 'SL': 0.0},\
				'lipid_droplets': {'PS': 0.0, 'PI': 0.0, 'PC': 0.0, 'PE': 0.0, 'CL': 0.0, 'PA': 0.0, 'ES': 0.0, 'SE': 0.5, 'TAG': 0.5, 'SL': 0.0}}

		membranes = self.membrane_lists.keys()
		for membrane in membranes:
			y = []
			y2 = []
			std = []
			x_names = self.lipids
			ind = np.arange(len(x_names))
			width = 0.3

			for lipid in x_names:
				y.append((self.membrane_comp_mean[membrane][lipid])*100)
				std.append((self.membrane_comp_std[membrane][lipid])*100)
				y2.append((zinser[membrane][lipid])*100)

			fig, ax = pyp.subplots()
			fig_membrane_comp = ax.bar(ind, y, width, color = 'b', yerr = std)
			fig_membrane_comp2 = ax.bar(ind + width, y2, width, color = 'y')

			ax.set_ylabel('%')
			ax.set_title(membrane.replace('_', ' ').title() + ' lipid composition', y = 1.03)
			ax.set_xticks(ind + width)
			ax.set_xticklabels(x_names)

			ax.legend((fig_membrane_comp[0], fig_membrane_comp2[0]), ('Model', 'Experimental'))
		
			pyp.show()

