import random
import matplotlib.pyplot as pyp
import numpy as np
import model

class stoch_model:
	def __init__(self, runs=10, sensitivity=False, change=0.1):
		# load model and parameters
		self.vera_model = model.model()
		self.rates = self.vera_model.rates
		self.probabilities = self.vera_model.probability

		# results
		self.results = []
		self.time = []

	def run(self, runs, sensitivity, change):
		# simulate
		# sensitivity analysis
		if sensitivity:
			self.sensitivity_analysis(runs, change)
		# no sensitivity analysis
		else:
			self.run_model(runs)

		# statistics
		self.calc_statistics()

		# plot results
		self.plot_results_fatty_acids()
		self.plot_results_membranes()

	def sensitivity_analysis(self, runs, change):
		'''for par in self.parameters:
			self.parameters[par] = self.parameters[par] + self.parameters[par] * change
			self.results[par + ' + 10 %'] = []
			for run in range(runs):
				self.results[par + ' + 10 %'].append(self.vera_model.run())
		for par in self.parameters:
			self.parameters[par] = self.parameters[par] - self.parameters[par] * change
			self.results[par + ' - 10 %'].append(self.vera_model.run())
			for run in range(runs):
				self.results[par + ' - 10 %'].append(self.vera_model.run())
		self.time = self.vera_model.t'''

		for par in self.rates:
			self.rates[par] = self.rates[par] + self.rates[par] * change
			self.results[par + ' + 10 %'] = []
			for run in range(runs):
				self.results[par + ' + 10 %'].append(self.vera_model.run())
		for par in self.rates:
			self.rates[par] = self.rates[par] - self.rates[par] * change
			self.results[par + ' - 10 %'].append(self.vera_model.run())
			for run in range(runs):
				self.results[par + ' - 10 %'].append(self.vera_model.run())
		self.time = self.vera_model.t		

	def run_model(self, runs):
		self.results = []
		for run in range(runs):
			self.results.append(self.vera_model.run())
		self.time = self.vera_model.t

	def calc_statistics(self):
		#fatty acids
		self.fatty_acid_distribution = {}
		self.fatty_acid_std = {}
		self.fatty_acids = ['C16:0', 'C16:1', 'C18:0', 'C18:1']

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
		#membrane compositions
		self.membrane_comp_mean = {}
		self.membrane_comp_std = {}
		self.lipids = ['CL', 'ES', 'PA', 'PC', 'PE', 'PI', 'PS', 'SE', 'SL', 'TAG']


		#mean and std fatty acid distribution
		for fa in self.fatty_acids:
			fa_relatives = []
			for runs in range(len(self.results)):
				fa_relatives.append(self.results[runs][0][fa])
				self.fatty_acid_distribution[fa] = np.mean(fa_relatives)
				self.fatty_acid_std[fa] = np.std(fa_relatives)

		#mean and std membrane lengths
		for runs in range(len(self.results)):
			self.plasma_membrane.append(self.results[runs][1][0])
			self.secretory_vesicles.append(self.results[runs][1][1])
			self.vacuoles.append(self.results[runs][1][2])
			self.nucleus.append(self.results[runs][1][3])
			self.peroxisomes.append(self.results[runs][1][4])
			self.light_microsomes.append(self.results[runs][1][5])
			self.inner_mit_membrane.append(self.results[runs][1][6])
			self.outer_mit_membrane.append(self.results[runs][1][7])
			self.lipid_droplets.append(self.results[runs][1][8])

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
		self.x_names = self.plots_fatty_acids.keys()
		self.x = [j-0.5 for j in range(len(self.x_names))]
		self.yerrs = []
		for i in range(len(self.plots_fatty_acids.values())):
			self.y.append(self.plots_fatty_acids.values()[i]['mean'])
			self.yerrs.append(self.plots_fatty_acids.values()[i]['std'])
		
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
