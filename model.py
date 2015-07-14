# -*- coding: utf-8 -*-
"""
Created/Started on Wed June 03 2015

@author: Vera
"""
import matplotlib.pyplot as mat
import random
from numpy.random import choice

class lipids(object):
	"""
	general class for all kinds of lipids
	with the head groups 'p', 'inositol', 'serine', 'ethanolamine', 'choline', 'neutral', 
	'cdp'(for cdp-dg) and 'None'(for tag)
	possible ffa for sn2: C14:0, C16:0, C16:1, C18:0, C18:1
	possible ffa for sn1: C16:1 C18:1
	"""
	def __init__(self, head, sn2, sn1, comp):

		self.head_groups = ['p', 'inositol', 'serine', 'ethanolamine', 'choline', 'neutral', 'cdp', 'sterol', None]
		self.sn2_options = ['C16:1', 'C18:1', None]	
		self.sn1_options = ['C14:0', 'C16:0', 'C18:0', None]
		self.compartment_options = ['plasma_membrane', 'secretory_vesicles', 'vacuoles', 'nucleus', 'peroxisomes', 'light_microsomes',\
							'inner_mit_membrane', 'outer_mit_membrane', 'lipid_droplets', None]	
		self.compartment = ['plasma_membrane', 'secretory_vesicles', 'vacuoles', 'nucleus', 'peroxisomes', 'light_microsomes',\
							'inner_mit_membrane', 'outer_mit_membrane']	
		self.compartment_weights = [0.3, 0.1, 0.15, 0.15, 0.05, 0.05, 0.1, 0.1]		
		self.plasma_membrane_comp = {'PS': 0.08062, 'PI': 0.04373, 'PC': 0.04164, 'PE': 0.04976, 'CL': 0.00313, 'PA': 0.01172, 'ES': 0.76800, 'TAG': 0.0}
		self.secretory_vesicles_comp = {'PS': 0.08205, 'PI': 0.11745, 'PC': 0.20824, 'PE': 0.13573, 'CL': 0.01239, 'PA': 0.01525, 'ES': 0.42900, 'TAG': 0.0}
		self.vacuoles_comp = {'PS': 0.04817, 'PI': 0.16604, 'PC': 0.40517, 'PE': 0.17537, 'CL': 0.02442, 'PA': 0.02866, 'ES': 0.15200, 'TAG': 0.0}
		self.nucleus_comp = {'PS': 0.04038, 'PI': 0.09650, 'PC': 0.27645, 'PE': 0.16848, 'CL': 0.01049, 'PA': 0.01781, 'ES': 0.390, 'TAG': 0.0}
		self.peroxisomes_comp = {'PS': 0.03235, 'PI': 0.11360, 'PC': 0.34656, 'PE': 0.16465, 'CL': 0.05033, 'PA': 0.01150, 'ES': 0.281, 'TAG': 0.0}
		self.light_microsomes_comp = {'PS': 0.05304, 'PI': 0.06019, 'PC': 0.40796, 'PE': 0.26583, 'CL': 0.00381, 'PA': 0.00222, 'ES': 0.206, 'TAG': 0.0}
		self.inner_mit_membrane_comp = {'PS': 0.02880, 'PI': 0.06019, 'PC': 0.29107, 'PE': 0.18192, 'CL': 0.12204, 'PA': 0.01137, 'ES': 0.242, 'TAG': 0.0}
		self.outer_mit_membrane_comp = {'PS': 0.01189, 'PI': 0.10108, 'PC': 0.45190, 'PE': 0.32307, 'CL': 0.05847, 'PA': 0.04360, 'ES': 0.009, 'TAG': 0.0}
		self.lipid_droplets_comp = {'PS': 0.0, 'PI': 0.0, 'PC': 0.0, 'PE': 0.0, 'CL': 0.0, 'PA':0.0, 'ES': 0.0, 'TAG': 1.0}
		self.membranes_comp = [self.plasma_membrane_comp, self.secretory_vesicles_comp, self.vacuoles_comp, self.nucleus_comp,\
								self.peroxisomes_comp, self.light_microsomes_comp, self.inner_mit_membrane_comp,\
								self.outer_mit_membrane_comp, self.lipid_droplets_comp]		

		self.head = head
		self.sn2 = sn2
		self.sn1 = sn1
		self.comp = comp

	@property
	def head(self):
		return self.__head
	@head.setter
	def head(self, group):
		if group not in self.head_groups:
			raise TypeError('This is no head group.')
		self.__head = group

	@property
	def sn2(self):
		return self.__sn2
	@sn2.setter
	def sn2(self, chain):
		if chain not in self.sn2_options:
			raise TypeError('This is no possible sn2 chain.')
		self.__sn2 = chain

	@property
	def sn1(self):
		return self.__sn1
	@sn1.setter
	def sn1(self, chain):
		if chain not in self.sn1_options:
			raise TypeError('This is no possible sn1 chain.')
		self.__sn1 = chain

	@property
	def comp(self):
		return self.__comp
	@comp.setter
	def comp(self, local):
		if local not in self.compartment_options:
			raise TypeError('This is no compartment.')
		self.__comp = local


	def comp_choice(self):
		if self.head == 'serine':
			weights = [self.compartment_weights[i] * self.membranes_comp[i]['PS'] for i in range(len(self.compartment_weights))]
			weights_normal = [weight / sum(weights) for weight in weights]
		elif self.head == 'inositol':
			weights = [self.compartment_weights[i] * self.membranes_comp[i]['PI'] for i in range(len(self.compartment_weights))]
			weights_normal = [weight / sum(weights) for weight in weights]
		elif self.head == 'choline':
			weights = [self.compartment_weights[i] * self.membranes_comp[i]['PC'] for i in range(len(self.compartment_weights))]
			weights_normal = [weight / sum(weights) for weight in weights]
		elif self.head == 'ethanolamine':
			weights = [self.compartment_weights[i] * self.membranes_comp[i]['PE'] for i in range(len(self.compartment_weights))]
			weights_normal = [weight / sum(weights) for weight in weights]
		elif self.head == 'p':
			weights = [self.compartment_weights[i] * self.membranes_comp[i]['PA'] for i in range(len(self.compartment_weights))]
			weights_normal = [weight / sum(weights) for weight in weights]
		self.comp = choice(self.compartment, p = weights_normal)


class TAG(lipids):

	def __init__(self, sn3, sn2, sn1):
		super(TAG, self).__init__(sn2, sn1)
		sn3 = None

	def comp_choice(self):
		self.comp = 'lipid_droplets'

class CL(lipids):

	def __init__(self, head, sn2, sn1, sn4, sn3, comp):
		super(CL, self).__init__(head, sn2, sn1, comp)
		sn4 = None
		sn3 = None	
		self.compartment = ['plasma_membrane', 'secretory_vesicles', 'vacuoles', 'nucleus', 'peroxisomes', 'light_microsomes',\
							'inner_mit_membrane', 'outer_mit_membrane']

	def comp_choice(self):
		weights = [self.compartment_weights[i] * self.membranes_comp[i]['CL'] for i in range(len(self.compartment_weights))]
		weights_normal = [weight / sum(weights) for weight in weights]
		self.comp = choice(self.compartment, p = weights_normal)


class fatty_acids(object):	#name als attribut statt der einzelnen unterklassen
	"""
	class for the fatty acids
	attribute C: number of C-Atoms
	attribute saturation: 0 = saturated, 1 = unsaturated
	"""
	def __init__(self, C, saturation):
		self.C = C
		self.saturation = saturation

	@property
	def C(self):
		return self.__C
	@C.setter
	def C(self, number):
		if not isinstance(number, int):
			raise TypeError('Number of C-Atoms must be an int.')
		self.__C = number

	@property
	def saturation(self):
		return self.__saturation
	@saturation.setter
	def saturation(self, number):
		if number not in [0,1]:
			raise TypeError('Saturation must be 0 or 1')
		self.__saturation = number


class enzyme(object):
	'''
	class for the enzymes that are part of the lipid metabolism
	attribute name: enzyme name
	attribute number: amount of the enzyme in one yeast cell
	attribute localisation: localisation of the enzyme and therefore of the reaction
	'''
	def __init__(self, name, number, localisation):
		self.localisation_options = ['cytoplasm', 'ER', 'lipid_particle', 'mitochondrion', 'inner_mit_membrane', 'outer_mit_membrane']
		self.name = name
		self.number = number
		self.localisation = localisation

	@property
	def number(self):
		return self.__number
	@number.setter
	def number(self, number):
		if not isinstance(number, int):
			raise TypeError('Number of enzyme must be an int.')
		self.__number = number


	@property
	def localisation(self):
		return self.__localisation
	@localisation.setter
	def localisation(self, comp):
		if comp not in self.localisation_options:
			raise TypeError('This is no compartment.')
		self.__localisation = comp


class sterol(object):

	def __init__(self, head, comp):
		self.head_groups = ['p', 'inositol', 'serine', 'ethanolamine', 'choline', 'neutral', 'cdp', 'sterol', None]
		self.compartment_options = ['plasma_membrane', 'secretory_vesicles', 'vacuoles', 'nucleus', 'peroxisomes', 'light_microsomes',\
							'inner_mit_membrane', 'outer_mit_membrane', 'lipid_droplets', None]
		self.compartment = ['plasma_membrane', 'secretory_vesicles', 'vacuoles', 'nucleus', 'peroxisomes', 'light_microsomes',\
							'inner_mit_membrane', 'outer_mit_membrane']
		self.compartment_weights = [0.3, 0.1, 0.15, 0.15, 0.05, 0.05, 0.1, 0.1]		
		self.plasma_membrane_comp = {'PS': 0.08062, 'PI': 0.04373, 'PC': 0.04164, 'PE': 0.04976, 'CL': 0.00313, 'PA': 0.01172, 'ES': 0.76800, 'TAG': 0.0}
		self.secretory_vesicles_comp = {'PS': 0.08205, 'PI': 0.11745, 'PC': 0.20824, 'PE': 0.13573, 'CL': 0.01239, 'PA': 0.01525, 'ES': 0.42900, 'TAG': 0.0}
		self.vacuoles_comp = {'PS': 0.04817, 'PI': 0.16604, 'PC': 0.40517, 'PE': 0.17537, 'CL': 0.02442, 'PA': 0.02866, 'ES': 0.15200, 'TAG': 0.0}
		self.nucleus_comp = {'PS': 0.04038, 'PI': 0.09650, 'PC': 0.27645, 'PE': 0.16848, 'CL': 0.01049, 'PA': 0.01781, 'ES': 0.390, 'TAG': 0.0}
		self.peroxisomes_comp = {'PS': 0.03235, 'PI': 0.11360, 'PC': 0.34656, 'PE': 0.16465, 'CL': 0.05033, 'PA': 0.01150, 'ES': 0.281, 'TAG': 0.0}
		self.light_microsomes_comp = {'PS': 0.05304, 'PI': 0.06019, 'PC': 0.40796, 'PE': 0.26583, 'CL': 0.00381, 'PA': 0.00222, 'ES': 0.206, 'TAG': 0.0}
		self.inner_mit_membrane_comp = {'PS': 0.02880, 'PI': 0.06019, 'PC': 0.29107, 'PE': 0.18192, 'CL': 0.12204, 'PA': 0.01137, 'ES': 0.242, 'TAG': 0.0}
		self.outer_mit_membrane_comp = {'PS': 0.01189, 'PI': 0.10108, 'PC': 0.45190, 'PE': 0.32307, 'CL': 0.05847, 'PA': 0.04360, 'ES': 0.009, 'TAG': 0.0}
		self.lipid_droplets_comp = {'PS': 0.0, 'PI': 0.0, 'PC': 0.0, 'PE': 0.0, 'CL': 0.0, 'PA':0.0, 'ES': 0.0, 'TAG': 1.0}
		self.membranes_comp = [self.plasma_membrane_comp, self.secretory_vesicles_comp, self.vacuoles_comp, self.nucleus_comp,\
								self.peroxisomes_comp, self.light_microsomes_comp, self.inner_mit_membrane_comp,\
								self.outer_mit_membrane_comp, self.lipid_droplets_comp]		
		self.head = head
		self.comp = comp

	def comp_choice(self):	
		weights = [self.compartment_weights[i] * self.membranes_comp[i]['ES'] for i in range(len(self.compartment_weights))]
		if sum(weights) != 0:
			weights_normal = [weight / sum(weights) for weight in weights]
			self.comp = choice(self.compartment, p = weights_normal)


	@property
	def head(self):
		return self.__head
	@head.setter
	def head(self, group):
		if group not in self.head_groups:
			raise TypeError('This is no head group.')
		self.__head = group

	@property
	def comp(self):
		return self.__comp
	@comp.setter
	def comp(self, local):
		if local not in self.compartment_options:
			raise TypeError('This is no compartment.')
		self.__comp = local



class model():
	"""
	The model. 
	At the beginning there are several lists defined which will contain the produced lipids. 
	The model starts with a given amount of precursors (pyruvate, dhap, ctp, serine, glucose-6-p, SAM, glycerol-3-p).
	After several reactions and the tranport function in the end there are membranes of different compartments with several different lipids.
	"""
	def __init__(self):
		#determining the timesteps
		self.timesteps = 7200			
		self.time = 0
		self.t = [i for i in range(self.timesteps)]

		# number of available precursors
		self.precursors_dict = {'pyruvate_number' : 5000000, 'dhap_number': 1000000, 'ctp_number': 1000000, 'serine_number': 1000000,\
									'glucose_6_p_number': 1000000, 'SAM_number': 1000000, 'SAH_number': 0, 'glycerol_3_p_mito_number': 1000000}
		
		self.inositol_number = 0
		self.acetyl_coa_number = 0

		self.co2_counter = 0
		self.p_counter = 0
		self.counter = 0

		#list of all enzymes that are part of the reactions of the lipid metabolism
		self.enzymes = [enzyme('Gat1_Ayr1', 1000, 'ER'), enzyme('Slc1', 1000, 'lipid_particle'), enzyme('Pis1', 1000, 'outer_mit_membrane'),\
						enzyme('Cho1', 1000, 'outer_mit_membrane'), enzyme('Psd1', 1000, 'inner_mit_membrane'), \
						enzyme('Cho2', 1000, 'ER'), enzyme('Opi3', 1000, 'mitochondrion'), enzyme('Pah1', 1000, 'cytoplasm'), \
						enzyme('Dga1', 1000, 'lipid_particle')]

		#list of the 4 cell cycle phases
		self.cell_cycle_phases = ['G1', 'S', 'G2', 'M']

		#empty lists for the produced lipids
		self.acyl_coa_list = []
		self.acyl_coa_list_saturated = []
		self.acyl_coa_list_unsaturated = []
		self.lyso_pa_list = []
		self.PA_list = []
		self.CDP_DG_list = []
		self.DAG_list = []
		self.TAG_list = []
		self.PS_list = []
		self.PI_list = []
		self.PE_list = []
		self.PC_list = []
		self.CL_list = []
		self.Ergosterol_list = []

		self.precursor_list = [self.acyl_coa_list, self.PA_list, self.CDP_DG_list, self.TAG_list, self.PS_list, self.PI_list,\
								self.PE_list, self.PC_list, self.CL_list, self.Ergosterol_list]

		self.lipid_lists = [self.PS_list, self.PI_list, self.PC_list, self.PE_list, self.CL_list, self.PA_list, self.Ergosterol_list, self.TAG_list]

		#lists to collect the transported lipids
		self.compartment = ['plasma_membrane', 'secretory_vesicles', 'vacuoles', 'nucleus', 'peroxisomes', 'light_microsomes',\
							'inner_mit_membrane', 'outer_mit_membrane', 'lipid_droplets']	

		self.plasma_membrane = []
		self.secretory_vesicles = []
		self.vacuoles = []
		self.nucleus = []
		self.peroxisomes = []
		self.light_microsomes = []
		self.inner_mit_membrane = []
		self.outer_mit_membrane = []
		self.lipid_droplets = []



		self.compartment_lists = [self.plasma_membrane, self.secretory_vesicles, self.vacuoles, self.nucleus, \
									self.peroxisomes, self.light_microsomes, self.inner_mit_membrane, \
									self.outer_mit_membrane, self.lipid_droplets]

		#collecting the products of every timestep
		self.number_acetyl_coa = [0]
		self.number_acyl_coa = [0]
		self.number_pa = [0]
		self.number_cdp_dg = [0]
		self.number_tag = [0]
		self.number_PS = [0]
		self.number_PI = [0]
		self.number_PE = [0]
		self.number_PC = [0]
		self.number_CL = [0]
		self.number_Ergosterol = [0]

		self.number_lipids_list = [self.number_acyl_coa, self.number_pa, self.number_cdp_dg, self.number_tag,\
									self.number_PS, self.number_PI, self.number_PE, self.number_PC,	self.number_CL,\
									self.number_Ergosterol]

		#counting the lipids in each membrane after every timestep
		self.number_plasma_membrane = [0]
		self.number_secretory_vesicles = [0]
		self.number_vacuoles = [0]
		self.number_nucleus = [0]
		self.number_peroxisomes = [0]
		self.number_light_microsomes = [0]
		self.number_inner_mit_membrane = [0]
		self.number_outer_mit_membrane = [0]
		self.number_lipid_droplets = [0]

		self.number_membranes_list = [self.number_plasma_membrane, self.number_secretory_vesicles, self.number_vacuoles, \
										self.number_nucleus, self.number_peroxisomes, self.number_light_microsomes,\
										self.number_inner_mit_membrane,	self.number_outer_mit_membrane, self.number_lipid_droplets]

		self.chainlength_saturated = {14: 'C14:0', 16: 'C16:0', 18: 'C18:0'}		
		self.chainlength_unsaturated = {16: 'C16:1', 18: 'C18:1'}

		
		self.membrane_lipids = ['PS', 'PI', 'PC', 'PE', 'CL', 'PA', 'ES', 'TAG']

		self.compartment_relatives_dict = {comp: dict(zip(self.membrane_lipids, [0.0 for z in range(8)])) for comp in self.compartment}

		self.rates = {'inositol_synthesis': 250, 'acetyl_coa_synthase': 350, 'acyl_synthase': 250, 'PA_synthese': 200, \
						'CDP_DG_synthase': 180, 'TAG_synthese': 80, 'TAG_lipase': 15, 'PS_synthase': 100, 'PI_synthase': 50,\
						'PE_synthase': 75, 'PC_synthase': 40, 'CL_synthase': 30, 'Ergosterol_synthase': 15}

		#functions to run the model
		for t in range(self.timesteps):
			self.time += 1
			self.cell_cycle()
			self.function_list = [self.inositol_synthesis,
								self.acetyl_coa_synthase,
								self.acyl_synthase,
								self.PA_synthese,
								self.CDP_DG_synthase,
								self.TAG_synthese,
								self.PS_synthase,
								self.PI_synthase,
								self.PE_synthase,
								self.PC_synthase,
								self.CL_synthase,
								self.TAG_lipase,
								self.Ergosterol_synthase,
								self.transport]

			for i in self.function_list:
				func = random.choice(self.function_list)
				func()
				self.function_list.remove(func)
			self.numbers()
			self.membrane_compositions()

		print 'CL: ' + str(self.number_CL[-1]), 'PS: ' + str(self.number_PS[-1]), 'PI: ' + str(self.number_PI[-1]), 'PE: ' + str(self.number_PE[-1]), \
				'PC: ' + str(self.number_PC[-1]), 'PA: ' + str(self.number_pa[-1]), 'TAG: ' + str(self.number_tag[-1]), 'CDP-DG: ' + str(self.number_cdp_dg[-1]),\
				'ES:' + str(self.number_Ergosterol[-1])
		print self.number_CL[-1] + self.number_PS[-1] + self.number_PI[-1] + self.number_PE[-1] + self.number_PC[-1] +\
				self.number_pa[-1] + self.number_tag[-1] + self.number_cdp_dg[-1] + self.number_Ergosterol[-1]
		print len(self.plasma_membrane) + len(self.secretory_vesicles) + len(self.vacuoles) + len(self.nucleus)+\
				len(self.peroxisomes) + len(self.light_microsomes) + len(self.inner_mit_membrane) + \
				len(self.outer_mit_membrane) + len(self.lipid_droplets)
		
#random die Funktionen hintereinander oder Pool vorher aufteilen und Anteile verteilen oder alle an Gesamtpool, aber Ausführen am Ende



	def plot_lipids(self):
		'''
		Plotting the produced lipids before they are transported into the membranes
		'''
		fig = mat.figure()
		ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
		ax.plot(self.t, self.number_tag[:-1], label = 'tag')
		ax.plot(self.t, self.number_PS[:-1], label = 'ps')
		ax.plot(self.t, self.number_PI[:-1], label = 'pi')
		ax.plot(self.t, self.number_PE[:-1], label = 'pe')
		ax.plot(self.t, self.number_PC[:-1], label = 'pc')
		ax.plot(self.t, self.number_CL[:-1], label = 'cl')
		ax.plot(self.t, self.number_Ergosterol[:-1], label = 'es')
		ax.legend(bbox_to_anchor = (1.05, 1), loc = 2, borderaxespad = 0.)
		mat.show()


	def plot_precursors(self):
		fig = mat.figure()
		ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
		ax.plot(self.t, self.number_acetyl_coa[:-1], label = 'acetyl_coa')
		ax.plot(self.t, self.number_acyl_coa[:-1], label = 'acyl_coa')
		ax.plot(self.t, self.number_pa[:-1], label = 'pa')
		ax.plot(self.t, self.number_cdp_dg[:-1], label = 'cdp-dg')
		ax.legend(bbox_to_anchor = (1.05, 1), loc = 2, borderaxespad = 0.)
		mat.show()


	def plot_membranes(self):
		'''
		Plotting the number of lipids in the membranes of different compartments
		'''
		fig = mat.figure()
		ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
		ax.plot(self.t, self.number_plasma_membrane[:-1], label = 'plasma membrane')
		ax.plot(self.t, self.number_secretory_vesicles[:-1], label = 'secretory vesicles')
		ax.plot(self.t, self.number_vacuoles[:-1], label = 'vacuoles')
		ax.plot(self.t, self.number_nucleus[:-1], label = 'nucleus')
		ax.plot(self.t, self.number_peroxisomes[:-1], label = 'peroxisomes')
		ax.plot(self.t, self.number_light_microsomes[:-1], label = 'light_microsomes')
		ax.plot(self.t, self.number_inner_mit_membrane[:-1], label = 'inner_mit_membrane')
		ax.plot(self.t, self.number_outer_mit_membrane[:-1], label = 'outer_mit_membrane')
		ax.plot(self.t, self.number_lipid_droplets[:-1], label = 'lipid_droplets')
		ax.legend(bbox_to_anchor = (1.05, 1), loc = 2, borderaxespad = 0.)
		mat.show()


	def cell_cycle(self):
		'''
		Function to determine the cell cycle phases depending on the elapsed time.
		'''
		if self.time <= 1800:
			self.phase = self.cell_cycle_phases[0]
		elif self.time <= 4500:
			self.phase = self.cell_cycle_phases[1]
		elif self.time <= 6300:
			self.phase = self.cell_cycle_phases[2]
		else:
			self.phase = self.cell_cycle_phases[3]	


	def inositol_synthesis(self):
		'''
		Synthesis of inositol from glucose-6-p by the Myo-inositol-1-p synthase and the Myo-inositol 1
		phosphatase.
		'''
		if self.precursors_dict['glucose_6_p_number'] > self.rates['inositol_synthesis']:
			self.inositol_number += self.rates['inositol_synthesis']
			self.precursors_dict['glucose_6_p_number'] -= self.rates['inositol_synthesis']
			self.p_counter -= self.rates['inositol_synthesis']


	def acetyl_coa_synthase(self):
		'''
		Synthesis of Acetyl-CoA: pyruvate dehydrogenase drives the reaction pyruvate to Acetyl-CoA, CO2 is released
		'''
		if self.precursors_dict['pyruvate_number'] >= self.rates['acetyl_coa_synthase']:			# transformation from pyruvate to acetyl_coa
			self.acetyl_coa_number += self.rates['acetyl_coa_synthase']
			self.precursors_dict['pyruvate_number'] -= self.rates['acetyl_coa_synthase']				
			self.co2_counter += self.rates['acetyl_coa_synthase']


	def acyl_synthase(self):
		'''
		Simplified synthesis of acyl_coa: several Acetyl-CoA are building a fatty acid (C14:0, C16:0, C16:1, C18:0 or C18:1)
		The intermediate Malonyl-CoA is leaved out.
		'''
		choice_list = [0, 1]
		for i in range(self.rates['acetyl_coa_synthase']):
			x = random.random()						#5 reactions in 1 timestep but only with a probability of 90%
			if self.acetyl_coa_number >= 2:		#control if at least 2 Acetyl-CoA are available
				if len(self.acyl_coa_list) == 0:		#starting the first reaction
					new_acyl = fatty_acids(2, choice(choice_list))
					self.acyl_coa_list.append(new_acyl)
					self.acyl_coa_list[-1].C += 2
					self.acetyl_coa_number -= 2

				elif self.acyl_coa_list[-1].C >=14 and x <= 0.1:
					self.acyl_coa_list[-1].saturation = 0
					new_acyl = fatty_acids(2, choice(choice_list))
					self.acyl_coa_list.append(new_acyl)
					self.acyl_coa_list[-1].C += 2
					self.acetyl_coa_number -= 2

				elif self.acyl_coa_list[-1].C >= 16 and x <= 0.45:	#stop the reaction cycle and starting a new one
					new_acyl = fatty_acids(2, choice(choice_list))
					self.acyl_coa_list.append(new_acyl)
					self.acyl_coa_list[-1].C += 2
					self.acetyl_coa_number -= 2

				elif self.acyl_coa_list[-1].C >= 18:				#stop the reaction cycle and starting a new one
					new_acyl = fatty_acids(2, choice(choice_list))
					self.acyl_coa_list.append(new_acyl)
					self.acyl_coa_list[-1].C += 2
					self.acetyl_coa_number -= 2
		
				else:									#adding an Acetyl_CoA to the growing ffa
					self.acyl_coa_list[-1].C += 2
					self.acetyl_coa_number -= 1

		if len(self.acyl_coa_list) > 1:
			for j in range(len(self.acyl_coa_list)-1):
				if self.acyl_coa_list[j].saturation == 0:
					self.acyl_coa_list_saturated.append(self.acyl_coa_list[j])
				elif self.acyl_coa_list[j].saturation == 1:
					self.acyl_coa_list_unsaturated.append(self.acyl_coa_list[j])
			del self.acyl_coa_list[:-1]



	def PA_synthese(self):
		'''
		Synthesis of PA in two reaction steps.
		'''
		for i in range(self.rates['PA_synthese']):
			self.lyso_PA_synthase()
			self.PA_synthase()


	def lyso_PA_synthase(self):
		'''
		Production of Lyso-PA by adding one acyl-coa to DHAP (sn1: always unsaturated) --> DHAP acyltransferase/acyl-DHAP reductase
		'''
		x = random.random()
		if x < 0.95 and len(self.acyl_coa_list_saturated) > 0 and self.precursors_dict['dhap_number']: 	#at least 1 ffa has to be unsaturated 
			sn1_chain = random.randint(0, (len(self.acyl_coa_list_saturated)-1))
			chainlength_sn1 = self.acyl_coa_list_saturated[sn1_chain].C
			lyso_pa = lipids('p', None, self.chainlength_saturated[chainlength_sn1], None)
			self.lyso_pa_list.append(lyso_pa)
			self.precursors_dict['dhap_number'] -= 1
			del self.acyl_coa_list_saturated[sn1_chain]


	def PA_synthase(self):	
		'''
		Synthesis of PA by adding the second fatty acid to DHAP (sn2: saturated or unsaturated) --> 1-acyl-sn-glycerol-3-phosphate acyltransferase
		'''
		for i in range (20):
			x = random.random()
			if x < 0.95 and len(self.acyl_coa_list_unsaturated) > 0 and len(self.lyso_pa_list) > 0:
				sn2_chain = random.randint(0, (len(self.acyl_coa_list_unsaturated)-1))		
				chainlength_sn2 = self.acyl_coa_list_unsaturated[sn2_chain].C
				self.lyso_pa_list[-1].sn2 = self.chainlength_unsaturated[chainlength_sn2]
				self.PA_list.append(self.lyso_pa_list[-1])			
				del self.acyl_coa_list_unsaturated[sn2_chain]		# deletion of the consumed ffa
				del self.lyso_pa_list[-1]


	def CDP_DG_synthase(self):
		'''
		PA is processed to CDP-DG (CDP-diacylglycerol synthase), that further reacts to the phospholipids
		'''
		for i in range(self.rates['CDP_DG_synthase']):
			x = random.random()
			if x <= 0.7 and self.precursors_dict['ctp_number'] > 0:
				if len(self.PA_list) > 0:
					self.PA_list[0].head = 'cdp'
					self.CDP_DG_list.append(self.PA_list[0])		#CDP-DG production from PA
					del self.PA_list[0]
					self.precursors_dict['ctp_number'] -= 1
					self.p_counter -= 2


	def TAG_synthese(self):
		'''
		Function for TAG synthesis divided in production of DAG and TAG afterwards
		'''
		for i in range(self.rates['TAG_synthese']):
			self.DAG_synthase()
			self.TAG_synthase()


	def DAG_synthase(self):
		'''
		DAG synthesis: Removing the head of the lipid and adding the lipid to the DAG list.
		'''
		if self.phase != 'G1':
			weights = [0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.3, 0.3]
		else:
			weights = [ 0.1 for p in range(10)]
		x = choice([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], p = weights)
		if x <= 0.8:
			if len(self.PA_list) > 0:
				self.PA_list[0].head = None
				self.DAG_list.append(self.PA_list[0])
				del self.PA_list[0]


	def TAG_synthase(self):
		'''
		DAG is processed to TAG by adding a third acyl-chain at position sn3.
		'''
		x = random.random()
		if x <= 0.85:
			if len(self.DAG_list) > 0 and len(self.acyl_coa_list_saturated) > 0 and len(self.acyl_coa_list_unsaturated) > 0:
				self.TAG_list.append(self.DAG_list[0])		
				self.TAG_list[-1].__class__ = TAG
				if x <= 0.4:# and len(self.acyl_coa_list_saturated) > 0:
					chainlength_sn3 = self.acyl_coa_list_saturated[0].C
					self.TAG_list[-1].sn3 = self.chainlength_saturated[chainlength_sn3]
					del self.acyl_coa_list_saturated[0]
				else:# len(self.acyl_coa_list_unsaturated) > 0:
					chainlength_sn3 = self.acyl_coa_list_unsaturated[0].C
					self.TAG_list[-1].sn3 = self.chainlength_unsaturated[chainlength_sn3]
					del self.acyl_coa_list_unsaturated[0]
				del self.DAG_list[0]
				self.p_counter -= 1


	def TAG_lipase(self):
		'''
		Cdk1/Cdc28-dependent activation of the major triacylglycerol lipase
		''' 
		if self.phase != 'G1' and len(self.lipid_droplets) > self.rates['TAG_lipase']:
			for i in range(self.rates['TAG_lipase']):
				x = random.random()
				if x <= 0.8:
					self.PA_list.append(self.lipid_droplets[0])
					self.PA_list[-1].__class__ = lipids
					self.PA_list[-1].head = 'p'
					if ':0' in self.lipid_droplets[0].sn3:
						for key, value in self.chainlength_unsaturated.items():
							if value == self.lipid_droplets[0].sn3:
								self.acyl_coa_list_saturated.append(fatty_acids(key, 0))
					elif ':1' in self.lipid_droplets[0].sn3:
						for key, value in self.chainlength_saturated.items():
							if value == self.lipid_droplets[0].sn3:
								self.acyl_coa_list_unsaturated.append(fatty_acids(key, 1))
					del self.lipid_droplets[0]
					self.ctp_number -= 1


	def PS_synthase(self):
		'''
		CDP-DG is processed to PS (PS synthase).
		'''
		for i in range(self.rates['PS_synthase']):
			x = random.random()
			if x <= 0.85 and len(self.CDP_DG_list) >= 1 and self.precursors_dict['serine_number'] > 0:
				self.CDP_DG_list[0].head = 'serine'				#PS synthesis from CDP-DG
				self.PS_list.append(self.CDP_DG_list[0])
				del self.CDP_DG_list[0]
				self.precursors_dict['serine_number'] -= 1


	def PI_synthase(self):
		'''
 		CDP-DG is processed to PI (PI synthase)
		'''
		for i in range(self.rates['PI_synthase']):
			x = random.random()
			if x <= 0.4 and len(self.CDP_DG_list) >= 1 and self.inositol_number > 0:
				self.CDP_DG_list[0].head = 'inositol'			#PI synthesis from CDP-DG
				self.PI_list.append(self.CDP_DG_list[0])
				del self.CDP_DG_list[0]
				self.inositol_number -= 1


	def PE_synthase(self):
		'''
		PE is derived from PS by releasing 1 CO2 --> PS decarboxylase.
		'''
		for i in range(self.rates['PE_synthase']):
			x = random.random()
			if x <= 0.9 and len(self.PS_list) >= 10:
				self.PS_list[0].head = 'ethanolamine'				#PE synthesis from PS
				self.PE_list.append(self.PS_list[0])
				del self.PS_list[0]
				self.co2_counter += 1
		

	def PC_synthase(self):
		'''
		PC is derived from PE. As enzymes serve 3 methyltransferases which need SAM and produce SAH as a side product.
		'''
		for i in range(self.rates['PC_synthase']):
			x = random.random()
			if x <= 0.4 and len(self.PE_list) >= 5 and self.precursors_dict['SAM_number'] >= 3:
				self.PE_list[0].head = 'choline'								#PC synthesis from PE
				self.PC_list.append(self.PE_list[0])
				del self.PE_list[0]
				self.precursors_dict['SAM_number'] -= 3
				self.precursors_dict['SAH_number'] += 3


	def CL_synthase(self):
		'''
		Synthesis of cardiolipin, for which 2 CDP-DG are needed. Different enzymes are needed.
		'''
		for i in range(self.rates['CL_synthase']):
			x = random.random()
			if x <= 0.4 and self.precursors_dict['glycerol_3_p_mito_number'] > 0 and len(self.CDP_DG_list) >= 2:
				self.CDP_DG_list[0].head = 'neutral'
				self.CL_list.append(self.CDP_DG_list[0])
				self.CL_list[-1].__class__ = CL
				self.CL_list[-1].sn4, self.CL_list[-1].sn3 = self.CDP_DG_list[1].sn2, self.CDP_DG_list[1].sn1
				del self.CDP_DG_list[0:2]
				self.precursors_dict['glycerol_3_p_mito_number'] -= 1


	def Ergosterol_synthase(self):
		'''
		Synthesis of the most existing sterol in yeast: ergosterol
		'''
		for i in range(self.rates['Ergosterol_synthase']):
			x = random.random()
			if x <= 0.7 and self.acetyl_coa_number >= 18:
				self.Ergosterol_list.append(sterol('sterol', None))
				self.acetyl_coa_number -= 18
				self.p_counter += 2



	def transport(self):
		'''
		General transport function for all produced lipids.
		'''
		i = 0
		for lipid in self.lipid_lists:
			for j in range(len(lipid)/10):
				lipid[0].comp_choice()
				if lipid[0].comp == 'plasma_membrane':
					self.plasma_membrane.append(lipid[0])
				elif lipid[0].comp == 'secretory_vesicles':
					self.secretory_vesicles.append(lipid[0])
				elif lipid[0].comp == 'vacuoles':
					self.vacuoles.append(lipid[0])
				elif lipid[0].comp == 'nucleus':
					self.nucleus.append(lipid[0])
				elif lipid[0].comp == 'peroxisomes':
					self.peroxisomes.append(lipid[0])
				elif lipid[0].comp == 'light_microsomes':
					self.light_microsomes.append(lipid[0])
				elif lipid[0].comp == 'inner_mit_membrane':
					self.inner_mit_membrane.append(lipid[0])
				elif lipid[0].comp == 'outer_mit_membrane':
					self.outer_mit_membrane.append(lipid[0])
				elif lipid[0].comp == 'lipid_droplets':
					self.lipid_droplets.append(lipid[0])
				del lipid[0]
			i += 1


	def membrane_compositions(self):
		'''
		Function to calculate the lipid composition of all membranes.
		'''
		x = 0
		for comp in self.compartment_lists:
			if len(comp) > 0:
				self.relatives_list = [(float(sum(j.head == 'serine' for j in comp)) / len(comp)),\
										(float(sum(j.head == 'inositol' for j in comp)) / len(comp)),\
										(float(sum(j.head == 'choline' for j in comp)) / len(comp)),\
										(float(sum(j.head == 'ethanolamine' for j in comp)) / len(comp)),\
										(float(sum(j.head == 'neutral' for j in comp)) / len(comp)),\
										(float(sum(j.head == 'p' for j in comp)) / len(comp)),\
										(float(sum(j.head == 'sterol' for j in comp)) / len(comp)),\
										(float(sum(j.comp == 'lipid_droplets' for j in comp)) / len(comp))]
				for i in range(len(self.relatives_list)):
					self.compartment_relatives_dict[self.compartment[x]][self.membrane_lipids[i]] = self.relatives_list[i]
			x += 1


	def numbers(self):		
		#for plotting the production of the lipids
		self.number_acetyl_coa.append(self.acetyl_coa_number)
		for current_precursor_number, number_of_precursor in zip(self.number_lipids_list, self.precursor_list):
			current_precursor_number.append(len(number_of_precursor))

		#for plotting the number of lipids in a certain membrane
		for current_membrane_number, number_of_membrane in zip(self.number_membranes_list, self.compartment_lists):
			current_membrane_number.append(len(number_of_membrane))


'''
if __name__ == '__main__':
	print 'meep'
	model()
	print 'fertsch'
'''