import random
import matplotlib.pyplot as pyp
import numpy as np
import model

class stoch_model:
	def __init__(self, runs=1000, sensitivity=False, change=0.1):
		# load model and parameters
		self.vera_model = model.model()
		self.parameters = vera_model.parameters

		# results
		self.results = {}
		self.plots = {}
		self.time = []

	def run(runs, sensitivity, change):
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
		for plot in self.plots:
			self.plot_results(plot, std=None)

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
		self.time = self.vera_model.time

	def run_model(self, runs):
		self.results['result'] = []
		for run in range(runs):
			self.results['result'].append(self.vera_model.run())
		self.time = self.vera_model.time

	def calc_statistics(self):
		for name in self.results:
			self.plots[name] = {}
			self.plots[name]['mean'] = np.mean(self.results[name])
			self.plots[name]['std'] = np.std(self.results[name])

	def plot_results(plot):
		if std:
		pyp.plot()
