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
		# results
		
		self.time = []

	def run(self, runs, sensitivity, change):
		# simulate
		# sensitivity analysis
		if sensitivity:
			self.results = {}
			self.sensitivity_analysis(runs, change)
			self.sensitivity_statistics()
		#plot heatmaps of sesitivity analysis	
			self.heatmap_fatty_acids()
			self.heatmap_membrane_comp()
			self.heatmap_membrane_length()

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

		#sensitivity of rates (+10%)
		for par in self.rates:
			self.vera_model.rates[par] = int(round(self.vera_model.rates[par] + self.vera_model.rates[par] * change))
			self.results[par + ' rate + 10 %'] = []
			for run in range(runs):
				self.results[par + ' rate + 10 %'].append(self.vera_model.run())
			self.vera_model.rates[par] = self.rates[par]
		
		#sensitivits of rates (-10%)
		for par in self.rates:
			self.vera_model.rates[par] = int(round(self.vera_model.rates[par] - self.vera_model.rates[par] * change))
			self.results[par + ' rate - 10 %'] = []
			for run in range(runs):
				self.results[par + ' rate - 10 %'].append(self.vera_model.run())
			self.vera_model.rates[par] = self.rates[par]

		#sensitivity of probabilities (+10%)
		for prob in self.probabilities:
			self.vera_model.probabilities[prob] = int(round(self.vera_model.probabilities[prob] + self.vera_model.probabilities[prob] * change))
			self.results[prob + ' probability + 10 %'] = []
			for run in range(runs):
				self.results[prob + ' probability + 10 %'].append(self.vera_model.run())
			self.vera_model.probabilities[prob] = self.probabilities[prob]

		#sensitivity of probabilities (-10%)
		for prob in self.probabilities:
			self.vera_model.probabilities[prob] = int(round(self.vera_model.probabilities[prob] - self.vera_model.probabilities[prob] * change))
			self.results[prob + ' probability - 10 %'] = []
			for run in range(runs):
				self.results[prob + ' probability - 10 %'].append(self.vera_model.run())
			self.vera_model.probabilities[prob] = self.probabilities[prob]	
		
		#time list imported from model
		self.time = self.vera_model.t	

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
		heatmaps_change = [' rate + 10 %', ' rate - 10 %', ' probability + 10 %', ' probability - 10 %']
		for change in heatmaps_change:
			if 'rate' in change: 
				column_labels = ['Wildtype']
				for par in self.rates:
					column_labels.append(par)

			else: 
				column_labels = ['Wildtype']
				for par in self.probabilities:
					column_labels.append(par)

			row_labels = self.fatty_acids
			self.heatmap_fa = []

			i = 0
			for reaction in column_labels:
				if reaction == 'Wildtype':
					self.heatmap_fa.append(i)
					difference = []
					for fa in row_labels:
						diff = 1
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
			fig, ax = pyp.plot()
			general_heatmap = ax.pcolor(self.heatmap_fa, cmap = pyp.cm.Blues)

			ax.set_xticks(np.arange(self.heatmap_fa.shape[1])+0.5, minor = False)
			ax.set_yticks(np.arange(self.heatmap_fa.shape[0])+0.5, minor = False)

			ax.invert_yaxis()
			ax.xaxis.tick_top()

			ax.set_xticklabels(row_labels, minor = False)
			ax.set_yticklabels(column_labels, minor = False)			

			cbar = pyp.colorbar(general_heatmap)
			cbar.ax.set_ylabel('change', rotation=270)

			ax.set_title('fatty acid distribution' + change, y = 1.05)
			
			pyp.show()

	def heatmap_membrane_comp(self):
		#heatmap of membrane compositions: x = lipids, y = reaction (4 heatmaps for each membrane, 4 x 9 heatmaps)
		column_labels = self.membrane_comp_mean.keys()

		for membrane in self.membrane_lists.keys():
			self.heatmap_membrane = []
			i = 0
			for reaction in column_labels:				
				self.heatmap_membrane.append(i)
				difference = []
				for lipid in self.lipids:
					if self.membrane_comp_mean['Wildtype'][membrane][lipid] == 0:
						diff = 1
					else:
						diff = self.membrane_comp_mean[reaction][membrane][lipid] / self.membrane_comp_mean['Wildtype'][membrane][lipid]
					difference.append(diff)
				self.heatmap_membrane[i] = difference
				i += 1

			self.heatmap_membrane = np.asarray(self.heatmap_membrane)
			fig, ax = pyp.plot()
			general_heatmap = ax.pcolor(self.heatmap_membrane, cmap = pyp.cm.Blues)

			ax.set_xticks(np.arange(self.heatmap_membrane.shape[1])+0.5, minor = False)
			ax.set_yticks(np.arange(self.heatmap_membrane.shape[0])+0.5, minor = False)

			ax.invert_yaxis()
			ax.xaxis.tick_top()
			ax.set_xticklabels(self.lipids, minor = False)
			ax.set_yticklabels(column_labels, minor = False)

			cbar = pyp.colorbar(general_heatmap)
			cbar.ax.set_ylabel('change', rotation=270)

			ax.set_title(membrane + 'composition', y = 1.05)
			
			pyp.show()


	def heatmap_membrane_length(self):
		#heatmap of membrane lengths: x = membrane, y = reaction (4 heatmaps: rates +/- 10%, probabilities +/- 10%)
		heatmaps_change = [' rate + 10 %', ' rate - 10 %', ' probability + 10 %', ' probability - 10 %']
		row_labels = self.membrane_lists.keys()

		for change in heatmaps_change:
			if 'rate' in change:
				column_labels = ['Wildtype']
				for par in self.rates:
					column_labels.append(par)
			else:
				column_labels = ['Wildtype']
				for par in self.probabilities:
					column_labels.append(par)

			self.heatmap_mem_length = []

			i = 0
			for reaction in column_labels:
				if reaction == 'Wildtype':
					self.heatmap_mem_length.append(i)
					difference = []
					for membrane in row_labels:
						diff = 1
						difference.append(diff)
					self.heatmap_mem_length[i] = difference
					i += 1

				else:
					self.heatmap_mem_length.append(i)
					difference = []
					for membrane in row_labels:
						diff = self.mean_time_lists[reaction + change][membrane][-1] / self.mean_time_lists['Wildtype'][membrane][-1]
						difference.append(diff)
					self.heatmap_mem_length[i] = difference
					i += 1

			self.heatmap_mem_length = np.asarray(self.heatmap_mem_length)
			fig, ax = pyp.plot()
			general_heatmap = ax.pcolor(self.heatmap_mem_length, cmap = pyp.cm.Blues)

			ax.set_xticks(np.arange(self.heatmap_mem_length.shape[1])+0.5, minor = False)
			ax.set_yticks(np.arange(self.heatmap_mem_length.shape[0])+0.5, minor = False)

			ax.invert_yaxis()
			ax.xaxis.tick_top()

			ax.set_xticklabels(row_labels, minor = False)
			ax.set_yticklabels(column_labels, minor = False)			

			cbar = pyp.colorbar(general_heatmap)
			cbar.ax.set_ylabel('change', rotation=270)

			ax.set_title('membrane lengths' + change, y = 1.05)
			
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
		x_names = self.fatty_acid_distribution.keys()
		x = [j-0.5 for j in range(len(self.x_names))]
		yerrs = []
		for i in x_names:
			y.append(self.fatty_acid_distribution[i])
			yerrs.append(self.fatty_acid_std[i])
		
		fig = pyp.bar(x, y, width = 0.5, yerr = yerrs)
		pyp.xticks(x, x_names)
		pyp.setp(fig, color='red', edgecolor = 'k')
		pyp.show()


	def plot_results_membranes(self):
		#graph of membrane growth: x = time, y = membrane length
		fig = pyp.figure()
		ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
		for membrane in self.membrane_lists.keys():
			ax.plot(self.time, self.mean_time_lists[membrane], label = membrane)
		ax.legend(bbox_to_anchor = (1.05, 1), loc = 2, borderaxespad = 0.)
		pyp.show()

		#histogram of the membrane lengths after the run: x = membrane, y = mean value membrane lengths
		x_names = self.membrane_lists.keys()
		x = [j-0.5 for j in range(len(x_names))]
		y = []
		std = []
		for membrane in x_names:
			y.append(self.membranes_length[membrane])
			std.append(self.membranes_length_std[membrane])
		fig_membrane_bars = pyp.bar(x, y, width = 0.5, yerr = std)
		pyp.xticks(x, x_names)
		pyp.setp(fig_membrane_bars, color='blue', edgecolor='k')
		pyp.show()

	
	def plot_membrane_comp(self):
		#histogram of the membrane compositions: x = lipids, y = percentages (a histogram for each membrane)
		zinser = {'plasma_membrane': {'PS': 0.17320, 'PI': 0.09124, 'PC': 0.08660, 'PE': 0.10464, 'CL': 0.00103, 'PA': 0.02010, 'ES': 0.48454, 'SE': 0.0, 'TAG': 0.0, 'SL': 0.03557},\
		'secretory_vesicles': {'PS': 0.08205, 'PI': 0.11745, 'PC': 0.20824, 'PE': 0.13573, 'CL': 0.01239, 'PA': 0.01525, 'ES': 0.42900, 'SE': 0.0, 'TAG': 0.0, 'SL': 0.05029},\
		'vacuoles': {'PS': 0.04817, 'PI': 0.16604, 'PC': 0.40517, 'PE': 0.17537, 'CL': 0.02442, 'PA': 0.02866, 'ES': 0.15200, 'SE': 0.0, 'TAG': 0.0, 'SL': 0.06525},\
		'nucleus': {'PS': 0.04038, 'PI': 0.09650, 'PC': 0.27645, 'PE': 0.16848, 'CL': 0.01049, 'PA': 0.01781, 'ES': 0.390, 'SE': 0.0, 'TAG': 0.0, 'SL': 0.02622},\
		'peroxisomes': {'PS': 0.03235, 'PI': 0.11360, 'PC': 0.34656, 'PE': 0.16465, 'CL': 0.05033, 'PA': 0.01150, 'ES': 0.281, 'SE': 0.0, 'TAG': 0.0, 'SL': 0.0},\
		'light_microsomes': {'PS': 0.05304, 'PI': 0.06019, 'PC': 0.40796, 'PE': 0.26583, 'CL': 0.00381, 'PA': 0.00222, 'ES': 0.206, 'SE': 0.0, 'TAG': 0.0, 'SL': 0.00397},\
		'inner_mit_membrane': {'PS': 0.02880, 'PI': 0.06019, 'PC': 0.29107, 'PE': 0.18192, 'CL': 0.12204, 'PA': 0.01137, 'ES': 0.242, 'SE': 0.0, 'TAG': 0.0, 'SL': 0.0},\
		'outer_mit_membrane': {'PS': 0.01189, 'PI': 0.10108, 'PC': 0.45190, 'PE': 0.32307, 'CL': 0.05847, 'PA': 0.04360, 'ES': 0.009, 'SE': 0.0, 'TAG': 0.0, 'SL': 0.0},\
		'lipid_droplets': {'PS': 0.0, 'PI': 0.0, 'PC': 0.0, 'PE': 0.0, 'CL': 0.0, 'PA': 0.0, 'ES': 0.0, 'SE': 0.5, 'TAG': 0.5, 'SL': 0.0}}

		membranes = self.membrane_lists.keys()
		for membrane in membranes:
			y = []
			y2 = []
			std = []
			x_names = self.lipids
			x = [j-0.5 for j in range(len(x_names))]

			for lipid in x_names:
				y.append(self.membrane_comp_mean[membrane][lipid])
				std.append(self.membrane_comp_std[membrane][lipid])
				y2.append(zinser[membrane][lipid])

			fig_membrane_comp = pyp.bar(x, y, width = 0.5, yerr = std)
			fig_membrane_comp2 = pyp.bar(x, y2, width = 0.5, yerr = std)
			pyp.xticks(x, x_names)
			pyp.setp(fig_membrane_comp, color='green', edgecolor = 'k')
			pyp.title(membrane + ' lipid composition')
			pyp.show()

