import random
import matplotlib.pyplot as pyp
import numpy as np
import model
import copy


class stoch_model:
	def __init__(self, runs=10, sensitivity=False, change=0.1):
		# load model and parameters
		self.vera_model = model.model()
		self.rates = copy.deepcopy(self.vera_model.rates)
		self.probabilities = copy.deepcopy(self.vera_model.probability)

		#calculation fatty acids
		self.fatty_acid_distribution = {}
		self.fatty_acid_std = {}
		self.fatty_acids = ['C16:0', 'C16:1', 'C18:0', 'C18:1']

		#membrane compositions
		self.membrane_comp_mean = {}
		self.membrane_comp_std = {}
		self.lipids = ['CL', 'ES', 'PA', 'PC', 'PE', 'PI', 'PS', 'SE', 'SL', 'TAG']

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
		# results
		
		self.time = []

	def run(self, runs, sensitivity, change):
		# simulate
		# sensitivity analysis
		if sensitivity:
			self.results = {}
			self.sensitivity_analysis(runs, change)
			self.sensitivity_statistics()
			self.heatmap()

		# no sensitivity analysis
		else:
			self.results = []
			self.run_model(runs)
			self.calc_statistics()		

		# plot results
			self.plot_results_fatty_acids()
			self.plot_results_membranes()

	def sensitivity_analysis(self, runs, change):
		self.results['Wildtype'] = []
		for run in range(runs):
			self.results['Wildtype'].append(self.vera_model.run())

		for par in self.rates:
			#old_rate = self.rates[par]
			self.vera_model.rates[par] = int(round(self.vera_model.rates[par] + self.vera_model.rates[par] * change))
			self.results[par + ' rate + 10 %'] = []
			for run in range(runs):
				self.results[par + ' rate + 10 %'].append(self.vera_model.run())
			self.vera_model.rates[par] = self.rates[par]
			
		for par in self.rates:
			#old_rate = self.rates[par]
			self.vera_model.rates[par] = int(round(self.vera_model.rates[par] - self.vera_model.rates[par] * change))
			self.results[par + ' rate - 10 %'] = []
			for run in range(runs):
				self.results[par + ' rate - 10 %'].append(self.vera_model.run())
			self.vera_model.rates[par] = self.rates[par]
		'''
		for prob in self.probabilities:
			self.vera_model.probabilities[prob] = int(round(self.vera_model.probabilities[prob] + self.vera_model.probabilities[prob] * change))
			self.results[prob + ' probability + 10 %'] = []
			for run in range(runs):
				self.results[prob + ' probability + 10 %'].append(self.vera_model.run())
			self.vera_model.probabilities[prob] = self.probabilities[prob]

		for prob in self.probabilities:
			self.vera_model.probabilities[prob] = int(round(self.vera_model.probabilities[prob] - self.vera_model.probabilities[prob] * change))
			self.results[prob + ' probability - 10 %'] = []
			for run in range(runs):
				self.results[prob + ' probability - 10 %'].append(self.vera_model.run())
			self.vera_model.probabilities[prob] = self.probabilities[prob]	
		'''
		self.time = self.vera_model.t	

	def sensitivity_statistics(self):
		#fatty acids
		for reaction in self.results:
			self.fatty_acid_distribution[reaction] = {}
			self.fatty_acid_std[reaction] = {}
			for fa in self.fatty_acids:
				fa_relatives = []
				for run in range(len(self.results[reaction])):
					fa_relatives.append(self.results[reaction][run][0][fa])
					self.fatty_acid_distribution[reaction][fa] = np.mean(fa_relatives)
					self.fatty_acid_std[reaction][fa] = np.std(fa_relatives)


		#mean and std membrane lengths
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


		#mean and std membrane compositions
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

	def heatmap(self):
		heatmaps_change = [' rate + 10 %', ' rate - 10 %']#, ' probability + 10 %', ' probability - 10 %']
		for change in heatmaps_change:

			column_labels = ['Wildtype']
			for par in self.rates:
				column_labels.append(par)
			row_labels = self.fatty_acids
			self.heatmap_fa = []

			i = 0
			for reaction in column_labels:
				if reaction == 'Wildtype':
					self.heatmap_fa.append(i)
					difference = []
					for fa in row_labels:
						diff = self.fatty_acid_distribution[reaction][fa] / self.fatty_acid_distribution['Wildtype'][fa]
						difference.append(diff)
					self.heatmap_fa[i] = difference
					i += 1

				else:
					self.heatmap_fa.append(i)
					difference = []
					for fa in row_labels:
						diff = self.fatty_acid_distribution[reaction + change][fa] / self.fatty_acid_distribution['Wildtype'][fa]
						difference.append(diff)
					self.heatmap_fa[i] = difference
					i += 1

			self.heatmap_fa = np.asarray(self.heatmap_fa)
			fig, ax = pyp.subplots()
			general_heatmap = ax.pcolor(self.heatmap_fa, cmap = pyp.cm.Blues)

			ax.set_xticks(np.arange(self.heatmap_fa.shape[1])+0.5, minor = False)
			ax.set_yticks(np.arange(self.heatmap_fa.shape[0])+0.5, minor = False)

			ax.invert_yaxis()
			ax.xaxis.tick_top()

			ax.set_xticklabels(row_labels, minor = False)
			ax.set_yticklabels(column_labels, minor = False)
			ax.colorbar()

			pyp.show()


	def run_model(self, runs):
		self.results = []
		for run in range(runs):
			self.results.append(self.vera_model.run())
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
		
		self.y = []
		self.x_names = self.fatty_acid_distribution.keys()
		self.x = [j-0.5 for j in range(len(self.x_names))]
		self.yerrs = []
		for i in self.x_names:
			self.y.append(self.fatty_acid_distribution[i])
			self.yerrs.append(self.fatty_acid_std[i])
		
		fig = pyp.bar(self.x, self.y, width = 0.5, yerr = self.yerrs)
		pyp.xticks(self.x, self.x_names)
		pyp.setp(fig, color='red', edgecolor = 'k')
		pyp.show()


	def plot_results_membranes(self):
		
		fig = pyp.figure()
		ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
		for membrane in self.membrane_lists.keys():
			ax.plot(self.time, self.mean_time_lists[membrane], label = membrane)
		ax.legend(bbox_to_anchor = (1.05, 1), loc = 2, borderaxespad = 0.)
		pyp.show()
