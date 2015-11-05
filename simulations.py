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
		self.plots_total_number = {}
		self.plots_fatty_acids = {'C16:0': {}, 'C16:1': {}, 'C18:0': {}, 'C18:1': {}}
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
		#for plot in self.plots:
		#	self.plot_results(plot, std=None)

	def sensitivity_analysis(self, runs, change):
		for par in self.parameters:
			self.parameters[par] = self.parameters[par] + self.parameters[par] * change
			self.results[par + ' + 10 %'] = []
			for run in range(runs):
				self.results[par + ' + 10 %'].append(self.vera_model.run())
		for par in self.parameters:
			self.parameters[par] = self.parameters[par] - self.parameters[par] * change
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
		self.total_numbers = []
		self.C16_0_relative = []
		self.C16_1_relative = []
		self.C18_0_relative = []
		self.C18_1_relative = []

		for number in range(len(self.results)):	
			self.total_numbers.append(self.results[number][0])
			self.C16_0_relative.append(self.results[number][1]['C16:0'])
			self.C16_1_relative.append(self.results[number][1]['C16:1'])
			self.C18_0_relative.append(self.results[number][1]['C18:0'])
			self.C18_1_relative.append(self.results[number][1]['C18:1'])

		self.plots_total_number['mean'] = np.mean(self.total_numbers)
		self.plots_total_number['std'] = np.std(self.total_numbers)	

		self.plots_fatty_acids['C16:0']['mean'] = np.mean(self.C16_0_relative)
		self.plots_fatty_acids['C16:0']['std'] = np.std(self.C16_0_relative)
		self.plots_fatty_acids['C16:1']['mean'] = np.mean(self.C16_1_relative)
		self.plots_fatty_acids['C16:1']['std'] = np.std(self.C16_1_relative)
		self.plots_fatty_acids['C18:0']['mean'] = np.mean(self.C18_0_relative)
		self.plots_fatty_acids['C18:0']['std'] = np.std(self.C18_0_relative)
		self.plots_fatty_acids['C18:1']['mean'] = np.mean(self.C18_1_relative)
		self.plots_fatty_acids['C18:1']['std'] = np.std(self.C18_1_relative)


		'''
		for name in range(len(self.results)):
			self.plots[name] = {}
			#for value in self.results[name]:
			self.plots[name]['mean'] = np.mean(self.results[name][0])
			self.plots[name]['std'] = np.std(self.results[name][0])
		'''
			
	def plot_results(plot):
		if std:
			pyp.plot()
