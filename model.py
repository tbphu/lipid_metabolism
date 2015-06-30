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

		self.head_groups = ['p', 'inositol', 'serine', 'ethanolamine', 'choline', 'neutral', 'cdp', None]
		self.sn2_options = ['C14:0', 'C16:0', 'C16:1', 'C18:0', 'C18:1', None]	
		self.sn1_options = ['C16:1', 'C18:1', None]
		self.compartment_options = ['plasma_membrane', 'secretory_vesicles', 'vacuoles', 'nucleus', 'peroxisomes', 'light_microsomes',\
							'inner_mit_membrane', 'outer_mit_membrane', 'lipid_droplets', None]					

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

class TAG(lipids):

	def __init__(self, sn3, sn2, sn1):
		super(TAG, self).__init__(sn2, sn1)
		sn3 = None

class CL(lipids):

	def __init__(self, head, sn2, sn1, sn4, sn3, comp):
		super(CL, self).__init__(head, sn2, sn1, comp)
		sn4 = None
		sn3 = None	


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


class model():
	"""
	The model. 
	At the beginning there are several lists defined which will contain the produced lipids. 
	The model starts with a given amount of precursors (pyruvate, dhap, ctp, serine, glucose-6-p, SAM, glycerol-3-p).
	After several reactions and the tranport function in the end there are membranes of different compartments with several different lipids.
	"""
	def __init__(self):
		#determining the timesteps
		self.timesteps = 720			
		self.time = 0
		self.t = [i for i in range(self.timesteps)]

		# number of available precursors
		self.precursors_dict = {'pyruvate_number' : 500000, 'dhap_number': 10000, 'ctp_number': 10000, 'serine_number': 10000,\
									'glucose_6_p_number': 10000, 'SAM_number': 10000, 'SAH_number': 0, 'glycerol_3_p_mito_number': 10000}
		
		self.inositol_number = 0
		self.acetyl_coa_number = 0

		self.co2_counter = 0
		self.p_counter = 0
		self.counter = 0

		self.compartment = ['plasma_membrane', 'secretory_vesicles', 'vacuoles', 'nucleus', 'peroxisomes', 'light_microsomes',\
							'inner_mit_membrane', 'outer_mit_membrane']	

		#list of all enzymes that are part of the reactions of the lipid metabolism
		self.enzymes = [enzyme('Gat1_Ayr1', 1000, 'ER'), enzyme('Slc1', 1000, 'lipid_particle'), enzyme('Pis1', 1000, 'outer_mit_membrane'),\
						enzyme('Cho1', 1000, 'outer_mit_membrane'), enzyme('Psd1', 1000, 'inner_mit_membrane'), \
						enzyme('Cho2', 1000, 'ER'), enzyme('Opi3', 1000, 'mitochondrion'), enzyme('Pah1', 1000, 'cytoplasm'), \
						enzyme('Dga1', 1000, 'lipid_particle')]

		#list of the 4 cell cycle phases
		self.cell_cycle_phases = ['G1', 'S', 'G2', 'M']

		#empty lists for the produced lipids
		self.acyl_coa_list = []
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

		#lists to collect the transported lipids
		self.plasma_membrane = []
		self.secretory_vesicles = []
		self.vacuoles = []
		self.nucleus = []
		self.peroxisomes = []
		self.light_microsomes = []
		self.inner_mit_membrane = []
		self.outer_mit_membrane = []
		self.lipid_droplets = []

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

		self.chainlength_saturated = {14: 'C14:0', 16: 'C16:0', 18: 'C18:0'}		
		self.chainlength_unsaturated = {16: 'C16:1', 18: 'C18:1'}

		self.membrane_lipids = ['PI', 'PS', 'PC', 'PE', 'CL']

		
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
								self.TAG_lipase]
								#self.transport]

			for i in self.function_list:
				func = random.choice(self.function_list)
				func()
				self.function_list.remove(func)
			self.numbers()

		print 'CL: ' + str(self.number_CL[-1]), 'PS: ' + str(self.number_PS[-1]), 'PI: ' + str(self.number_PI[-1]), 'PE: ' + str(self.number_PE[-1]), \
				'PC: ' + str(self.number_PC[-1]), 'PA: ' + str(self.number_pa[-1]), 'TAG: ' + str(self.number_tag[-1]), 'CDP-DG: ' + str(self.number_cdp_dg[-1])
		print self.number_CL[-1] + self.number_PS[-1] + self.number_PI[-1] + self.number_PE[-1] + self.number_PC[-1] +\
				self.number_pa[-1] + self.number_tag[-1] + self.number_cdp_dg[-1]
		#self.membranes_composition()
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
		if self.time <= 180:
			self.phase = self.cell_cycle_phases[0]
		elif self.time <= 450:
			self.phase = self.cell_cycle_phases[1]
		elif self.time <= 630:
			self.phase = self.cell_cycle_phases[2]
		else:
			self.phase = self.cell_cycle_phases[3]	


	def inositol_synthesis(self):
		'''
		Synthesis of inositol from glucose-6-p by the Myo-inositol-1-p synthase and the Myo-inositol 1
		phosphatase.
		'''
		for i in range(100):
			if self.precursors_dict['glucose_6_p_number'] > 0:
				self.inositol_number += 1
				self.precursors_dict['glucose_6_p_number'] -= 1
				self.p_counter -= 1


	def acetyl_coa_synthase(self):
		'''
		Synthesis of Acetyl-CoA: pyruvate dehydrogenase drives the reaction pyruvate to Acetyl-CoA, CO2 is released
		'''
		for i in range(200):			
			if self.precursors_dict['pyruvate_number'] >= 1:			# transformation from pyruvate to acetyl_coa
				self.acetyl_coa_number += 1
				self.precursors_dict['pyruvate_number'] -= 1				
				self.co2_counter += 1


	def acyl_synthase(self):
		'''
		Simplified synthesis of acyl_coa: several Acetyl-CoA are building a fatty acid (C14, C16:0, C16:1, C18:0 or C18:1)
		The intermediate Malonyl-CoA is leaved out.
		'''
		choice_list = [0, 1]
		weights = [0.05, 0.95]
		for i in range(150):
			x = random.random()						#5 reactions in 1 timestep but only with a probability of 90%
			if self.acetyl_coa_number >= 2:		#control if at least 2 Acetyl-CoA are available
				if len(self.acyl_coa_list) == 0:		#starting the first reaction
					new_acyl = fatty_acids(2, choice(choice_list, p = weights))
					self.acyl_coa_list.append(new_acyl)
					self.acyl_coa_list[-1].C += 2
					self.acetyl_coa_number -= 2

				elif self.acyl_coa_list[-1].C >=14 and x <= 0.1:
					self.acyl_coa_list[-1].saturation = 0
					new_acyl = fatty_acids(2, choice(choice_list, p = weights))
					self.acyl_coa_list.append(new_acyl)
					self.acyl_coa_list[-1].C += 2
					self.acetyl_coa_number -= 2

				elif self.acyl_coa_list[-1].C >= 16 and x <= 0.45:	#stop the reaction cycle and starting a new one
					new_acyl = fatty_acids(2, choice(choice_list, p = weights))
					self.acyl_coa_list.append(new_acyl)
					self.acyl_coa_list[-1].C += 2
					self.acetyl_coa_number -= 2

				elif self.acyl_coa_list[-1].C >= 18:				#stop the reaction cycle and starting a new one
					new_acyl = fatty_acids(2, choice(choice_list, p = weights))
					self.acyl_coa_list.append(new_acyl)
					self.acyl_coa_list[-1].C += 2
					self.acetyl_coa_number -= 2
		
				else:									#adding an Acetyl_CoA to the growing ffa
					self.acyl_coa_list[-1].C += 2
					self.acetyl_coa_number -= 1


	def PA_synthese(self):
		'''
		Synthesis of PA in two reaction steps.
		'''
		for i in range(80):
			self.lyso_PA_synthase()
			self.PA_synthase()


	def lyso_PA_synthase(self):
		'''
		Production of Lyso-PA by adding one acyl-coa to DHAP (sn1: always unsaturated) --> DHAP acyltransferase/acyl-DHAP reductase
		'''
		x = random.random()
		if x < 0.95 and 1 in [self.acyl_coa_list[z].saturation for z in range(len(self.acyl_coa_list))]: 	#at least 1 ffa has to be unsaturated 
			if self.precursors_dict['dhap_number'] > 0 and len(self.acyl_coa_list) > 1:
				sn1_chain = random.randint(0, (len(self.acyl_coa_list)-2))
				if self.acyl_coa_list[sn1_chain].saturation == 1:
					chainlength_sn1 = self.acyl_coa_list[sn1_chain].C
					lyso_pa = lipids('p', None, self.chainlength_unsaturated[chainlength_sn1], None)
					self.lyso_pa_list.append(lyso_pa)
					self.precursors_dict['dhap_number'] -= 1
					del self.acyl_coa_list[sn1_chain]
				else:
					self.lyso_PA_synthase()
		else:
			self.counter +=1


	def PA_synthase(self):	
		'''
		Synthesis of PA by adding the second fatty acid to DHAP (sn2: saturated or unsaturated) --> 1-acyl-sn-glycerol-3-phosphate acyltransferase
		'''
		x = random.random()
		if x < 0.95: 		
			if len(self.lyso_pa_list) > 0 and len(self.acyl_coa_list) > 1:		# available ffa
				sn2_chain = random.randint(0, (len(self.acyl_coa_list)-2))		
				chainlength_sn2 = self.acyl_coa_list[sn2_chain].C
				if self.acyl_coa_list[sn2_chain].saturation == 0:
					self.lyso_pa_list[-1].sn2 = self.chainlength_saturated[chainlength_sn2]	# creating a new lipid: PA
				elif self.acyl_coa_list[sn2_chain].saturation == 1:
					self.lyso_pa_list[-1].sn2 = self.chainlength_unsaturated[chainlength_sn2]
				self.PA_list.append(self.lyso_pa_list[-1])			
				del self.acyl_coa_list[sn2_chain]		# deletion of the consumed ffa
				del self.lyso_pa_list[-1]


	def CDP_DG_synthase(self):
		'''
		PA is processed to CDP-DG (CDP-diacylglycerol synthase), that further reacts to the phospholipids
		'''
		for i in range(30):
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
		self.DAG_synthase()
		self.TAG_synthase()


	def DAG_synthase(self):
		'''
		DAG synthesis: Removing the head of the lipid and adding the lipid to the DAG list.
		'''
		for i in range(10):
			x = random.random()
			if x <= 0.8:
				if len(self.PA_list) > 0:
					self.PA_list[0].head = None
					self.DAG_list.append(self.PA_list[0])
					del self.PA_list[0]


	def TAG_synthase(self):
		'''
		DAG is processed to TAG by adding a third acyl-chain at position sn3.
		'''
		for i in range(10):
			x = random.random()
			if x <= 0.85:
				if len(self.acyl_coa_list) > 1 and len(self.DAG_list) > 0:
					self.TAG_list.append(self.DAG_list[0])		
					self.TAG_list[-1].__class__ = TAG
					chainlength_sn3 = self.acyl_coa_list[0].C
					if self.acyl_coa_list[0].saturation == 0:
						self.TAG_list[-1].sn3 = self.chainlength_saturated[chainlength_sn3]
					elif self.acyl_coa_list[0].saturation == 1:
						self.TAG_list[-1].sn3 = self.chainlength_unsaturated[chainlength_sn3]
					del self.DAG_list[0]
					del self.acyl_coa_list[0]
					self.p_counter -= 1


	def TAG_lipase(self):
		'''
		Cdk1/Cdc28-dependent activation of the major triacylglycerol lipase
		''' 
		if self.phase != 'G1' and len(self.TAG_list) > 2:
			for i in range(3):
				self.PA_list.append(self.TAG_list[0])
				self.PA_list[-1].__class__ = lipids
				self.PA_list[-1].head = 'p'
				if ':0' in self.TAG_list[0].sn3:
					for key, value in self.chainlength_unsaturated.items():
						if value == self.TAG_list[0].sn3:
							self.acyl_coa_list.append(fatty_acids(key, 0))
				elif ':1' in self.TAG_list[0].sn3:
					for key, value in self.chainlength_saturated.items():
						if value == self.TAG_list[0].sn3:
							self.acyl_coa_list.append(fatty_acids(key, 1))
				del self.TAG_list[0]

	def PS_synthase(self):
		'''
		CDP-DG is processed to PS (PS synthase).
		'''
		for i in range(6):
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
		for i in range(2):
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
		for i in range(5):
			x = random.random()
			if x <= 0.9 and len(self.PS_list) >= 1:
				self.PS_list[0].head = 'ethanolamine'				#PE synthesis from PS
				self.PE_list.append(self.PS_list[0])
				del self.PS_list[0]
				self.co2_counter += 1
		

	def PC_synthase(self):
		'''
		PC is derived from PE. As enzymes serve 3 methyltransferases which need SAM and produce SAH as a side product.
		'''
		for i in range(4):
			x = random.random()
			if x <= 0.7 and len(self.PE_list) >= 1 and self.precursors_dict['SAM_number'] >= 3:
				self.PE_list[0].head = 'choline'								#PC synthesis from PE
				self.PC_list.append(self.PE_list[0])
				del self.PE_list[0]
				self.precursors_dict['SAM_number'] -= 3
				self.precursors_dict['SAH_number'] += 3


	def CL_synthase(self):
		'''
		Synthesis of cardiolipin, for which 2 CDP-DG are needed. Different enzymes are needed.
		'''
		x = random.random()
		if x <= 0.3 and self.precursors_dict['glycerol_3_p_mito_number'] > 0 and len(self.CDP_DG_list) >= 2:
			self.CDP_DG_list[0].head = 'neutral'
			self.CL_list.append(self.CDP_DG_list[0])
			self.CL_list[-1].__class__ = CL
			self.CL_list[-1].sn4, self.CL_list[-1].sn3 = self.CDP_DG_list[1].sn2, self.CDP_DG_list[1].sn1
			del self.CDP_DG_list[0:2]
			self.precursors_dict['glycerol_3_p_mito_number'] -= 1


	def transport(self):
		'''
		General transport function for all produced lipids.
		'''
		self.transport_PS()
		self.transport_PI()
		self.transport_PE()
		self.transport_PC()
		self.transport_CL()
		self.transport_TAG()


	def transport_PS(self):
		'''
		Transport of PS to the different compartment membranes. The distribution to the different membranes is random.
		'''
		weights = [0.7, 0.15, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025]
		if len(self.PS_list) > 3:
			z = len(self.PS_list)/10
			for i in range(z):
				self.PS_list[0].comp = choice(self.compartment, p = weights)
				if self.PS_list[0].comp == 'plasma_membrane':
					self.plasma_membrane.append(self.PS_list[0])
				elif self.PS_list[0].comp == 'secretory_vesicles':
					self.secretory_vesicles.append(self.PS_list[0])
				elif self.PS_list[0].comp == 'vacuoles':
					self.vacuoles.append(self.PS_list[0])
				elif self.PS_list[0].comp == 'nucleus':
					self.nucleus.append(self.PS_list[0])
				elif self.PS_list[0].comp == 'peroxisomes':
					self.peroxisomes.append(self.PS_list[0])
				elif self.PS_list[0].comp == 'light_microsomes':
					self.light_microsomes.append(self.PS_list[0])
				elif self.PS_list[0].comp == 'inner_mit_membrane':
					self.inner_mit_membrane.append(self.PS_list[0])
				elif self.PS_list[0].comp == 'outer_mit_membrane':
					self.outer_mit_membrane.append(self.PS_list[0])
				del self.PS_list[0]



	def transport_PI(self):
		'''
		Transport of PI to the different compartment membranes. The distribution to the different membranes is random.
		'''
		weights = [0.18, 0.15, 0.12, 0.11, 0.11, 0.11, 0.11, 0.11]
		if len(self.PI_list) > 3:
			z = len(self.PI_list)/10
			for i in range(z):
				self.PI_list[0].comp = choice(self.compartment, p = weights)
				if self.PI_list[0].comp == 'plasma_membrane':
					self.plasma_membrane.append(self.PI_list[0])
				elif self.PI_list[0].comp == 'secretory_vesicles':
					self.secretory_vesicles.append(self.PI_list[0])
				elif self.PI_list[0].comp == 'vacuoles':
					self.vacuoles.append(self.PI_list[0])
				elif self.PI_list[0].comp == 'nucleus':
					self.nucleus.append(self.PI_list[0])
				elif self.PI_list[0].comp == 'peroxisomes':
					self.peroxisomes.append(self.PI_list[0])
				elif self.PI_list[0].comp == 'light_microsomes':
					self.light_microsomes.append(self.PI_list[0])
				elif self.PI_list[0].comp == 'inner_mit_membrane':
					self.inner_mit_membrane.append(self.PI_list[0])
				elif self.PI_list[0].comp == 'outer_mit_membrane':
					self.outer_mit_membrane.append(self.PI_list[0])
				del self.PI_list[0]



	def transport_PE(self):
		'''
		Transport of PE to the different compartment membranes. The distribution to the different membranes is random.
		'''
		weights = [0.15, 0.18, 0.22, 0.09, 0.09, 0.09, 0.09, 0.09]
		if len(self.PE_list) > 3:
			z = len(self.PE_list)/10
			for i in range(z):
				self.PE_list[0].comp = choice(self.compartment, p = weights)
				if self.PE_list[0].comp == 'plasma_membrane':
					self.plasma_membrane.append(self.PE_list[0])
				elif self.PE_list[0].comp == 'secretory_vesicles':
					self.secretory_vesicles.append(self.PE_list[0])
				elif self.PE_list[0].comp == 'vacuoles':
					self.vacuoles.append(self.PE_list[0])
				elif self.PE_list[0].comp == 'nucleus':
					self.nucleus.append(self.PE_list[0])
				elif self.PE_list[0].comp == 'peroxisomes':
					self.peroxisomes.append(self.PE_list[0])
				elif self.PE_list[0].comp == 'light_microsomes':
					self.light_microsomes.append(self.PE_list[0])
				elif self.PE_list[0].comp == 'inner_mit_membrane':
					self.inner_mit_membrane.append(self.PE_list[0])
				elif self.PE_list[0].comp == 'outer_mit_membrane':
					self.outer_mit_membrane.append(self.PE_list[0])
				del self.PE_list[0]



	def transport_PC(self):
		'''
		Transport of PC to the different compartment membranes. The distribution to the different membranes is random.
		'''
		weights = [0.1, 0.25, 0.11, 0.06, 0.12, 0.12, 0.12, 0.12]
		if len(self.PC_list) > 3:
			z = len(self.PC_list)/10
			for i in range(z):
				self.PC_list[0].comp = choice(self.compartment, p = weights)
				if self.PC_list[0].comp == 'plasma_membrane':
					self.plasma_membrane.append(self.PC_list[0])
				elif self.PC_list[0].comp == 'secretory_vesicles':
					self.secretory_vesicles.append(self.PC_list[0])
				elif self.PC_list[0].comp == 'vacuoles':
					self.vacuoles.append(self.PC_list[0])
				elif self.PC_list[0].comp == 'nucleus':
					self.nucleus.append(self.PC_list[0])
				elif self.PC_list[0].comp == 'peroxisomes':
					self.peroxisomes.append(self.PC_list[0])
				elif self.PC_list[0].comp == 'light_microsomes':
					self.light_microsomes.append(self.PC_list[0])
				elif self.PC_list[0].comp == 'inner_mit_membrane':
					self.inner_mit_membrane.append(self.PC_list[0])
				elif self.PC_list[0].comp == 'outer_mit_membrane':
					self.outer_mit_membrane.append(self.PC_list[0])
				del self.PC_list[0]



	def transport_CL(self):
		'''
		Transport of CL to the different compartment membranes. The distribution to the different membranes is random.
		'''
		weights = [0.05, 0.05, 0.20, 0.14, 0.14, 0.14, 0.14, 0.14]
		if len(self.CL_list) > 1:
			z = len(self.CL_list)/10
			for i in range(z):
				self.CL_list[0].comp = choice(self.compartment, p = weights)
				if self.CL_list[0].comp == 'plasma_membrane':
					self.plasma_membrane.append(self.CL_list[0])
				elif self.CL_list[0].comp == 'secretory_vesicles':
					self.secretory_vesicles.append(self.CL_list[0])
				elif self.CL_list[0].comp == 'vacuoles':
					self.vacuoles.append(self.CL_list[0])
				elif self.CL_list[0].comp == 'nucleus':
					self.nucleus.append(self.CL_list[0])
				elif self.CL_list[0].comp == 'peroxisomes':
					self.peroxisomes.append(self.CL_list[0])
				elif self.CL_list[0].comp == 'light_microsomes':
					self.light_microsomes.append(self.CL_list[0])
				elif self.CL_list[0].comp == 'inner_mit_membrane':
					self.inner_mit_membrane.append(self.CL_list[0])
				elif self.CL_list[0].comp == 'outer_mit_membrane':
					self.outer_mit_membrane.append(self.CL_list[0])
				del self.CL_list[0]

	def transport_TAG(self):
		'''
		Transport of TAG to lipid droplets
		'''
		if len(self.TAG_list) > 3:
			z = len(self.TAG_list)/10
			for i in range(z):
				self.TAG_list[0].comp = 'lipid_droplets'
				self.lipid_droplets.append(self.TAG_list[0])
				del self.TAG_list[0]


	def membranes_composition(self):
		'''
		Function to calculate the relative composition of all membranes. Each is calculated by counting the lipids of every kind
		and dividing this number by the the number of all lipids in this membrane.
		'''
		self.plasma_membrane_composition()
		self.secretory_vesicles_composition()
		self.vacuoles_composition()
		self.nucleus_composition()
		self.peroxisomes_composition()
		self.light_microsomes_composition()
		self.inner_mit_membrane_composition()
		self.outer_mit_membrane_composition()


	def plasma_membrane_composition(self):
		if len(self.plasma_membrane) > 0:
			self.PI_plasma_membrane_absolut = float(sum(i.head == 'inositol' for i in self.plasma_membrane))
			self.PS_plasma_membrane_absolut = float(sum(i.head == 'serine' for i in self.plasma_membrane))
			self.PC_plasma_membrane_absolut = float(sum(i.head == 'choline' for i in self.plasma_membrane))
			self.PE_plasma_membrane_absolut = float(sum(i.head == 'ethanolamine' for i in self.plasma_membrane))
			self.CL_plasma_membrane_absolut = float(sum(i.head == 'neutral' for i in self.plasma_membrane))

			self.PI_plasma_membrane_relative = self.PI_plasma_membrane_absolut / len(self.plasma_membrane)
			self.PS_plasma_membrane_relative = self.PS_plasma_membrane_absolut / len(self.plasma_membrane)
			self.PC_plasma_membrane_relative = self.PC_plasma_membrane_absolut / len(self.plasma_membrane)
			self.PE_plasma_membrane_relative = self.PE_plasma_membrane_absolut / len(self.plasma_membrane)
			self.CL_plasma_membrane_relative = self.CL_plasma_membrane_absolut / len(self.plasma_membrane)

			self.plasma_membrane_relatives = [self.PI_plasma_membrane_relative, self.PS_plasma_membrane_relative,\
										 	  self.PC_plasma_membrane_relative, self.PE_plasma_membrane_relative,\
										 	  self.CL_plasma_membrane_relative]

			self.plasma_membrane_comp = dict(zip(self.membrane_lipids, self.plasma_membrane_relatives))


	def secretory_vesicles_composition(self):
		if len(self.secretory_vesicles) >0:
			self.PI_secretory_vesicles_absolut = float(sum(i.head == 'inositol' for i in self.secretory_vesicles))
			self.PS_secretory_vesicles_absolut = float(sum(i.head == 'serine' for i in self.secretory_vesicles))
			self.PC_secretory_vesicles_absolut = float(sum(i.head == 'choline' for i in self.secretory_vesicles))
			self.PE_secretory_vesicles_absolut = float(sum(i.head == 'ethanolamine' for i in self.secretory_vesicles))
			self.CL_secretory_vesicles_absolut = float(sum(i.head == 'neutral' for i in self.secretory_vesicles))

			self.PI_secretory_vesicles_relative = self.PI_secretory_vesicles_absolut / len(self.secretory_vesicles)
			self.PS_secretory_vesicles_relative = self.PS_secretory_vesicles_absolut / len(self.secretory_vesicles)
			self.PC_secretory_vesicles_relative = self.PC_secretory_vesicles_absolut / len(self.secretory_vesicles)
			self.PE_secretory_vesicles_relative = self.PE_secretory_vesicles_absolut / len(self.secretory_vesicles)
			self.CL_secretory_vesicles_relative = self.CL_secretory_vesicles_absolut / len(self.secretory_vesicles)

			self.secretory_vesicles_relatives = [self.PI_secretory_vesicles_relative, self.PS_secretory_vesicles_relative,\
												  self.PC_secretory_vesicles_relative, self.PE_secretory_vesicles_relative,\
												  self.CL_secretory_vesicles_relative]

			self.secretory_vesicles_comp = dict(zip(self.membrane_lipids, self.secretory_vesicles_relatives))


	def vacuoles_composition(self):
		if len(self.vacuoles) >0:
			self.PI_vacuoles_absolut = float(sum(i.head == 'inositol' for i in self.vacuoles))
			self.PS_vacuoles_absolut = float(sum(i.head == 'serine' for i in self.vacuoles))
			self.PC_vacuoles_absolut = float(sum(i.head == 'choline' for i in self.vacuoles))
			self.PE_vacuoles_absolut = float(sum(i.head == 'ethanolamine' for i in self.vacuoles))
			self.CL_vacuoles_absolut = float(sum(i.head == 'neutral' for i in self.vacuoles))

			self.PI_vacuoles_relative = self.PI_vacuoles_absolut / len(self.vacuoles)
			self.PS_vacuoles_relative = self.PS_vacuoles_absolut / len(self.vacuoles)
			self.PC_vacuoles_relative = self.PC_vacuoles_absolut / len(self.vacuoles)
			self.PE_vacuoles_relative = self.PE_vacuoles_absolut / len(self.vacuoles)
			self.CL_vacuoles_relative = self.CL_vacuoles_absolut / len(self.vacuoles)

			self.vacuoles_relatives = [self.PI_vacuoles_relative, self.PS_vacuoles_relative, self.PC_vacuoles_relative,\
										self.PE_vacuoles_relative, self.CL_vacuoles_relative]

			self.vacuoles_comp = dict(zip(self.membrane_lipids, self.vacuoles_relatives))


	def nucleus_composition(self):
		if len(self.nucleus) >0:
			self.PI_nucleus_absolut = float(sum(i.head == 'inositol' for i in self.nucleus))
			self.PS_nucleus_absolut = float(sum(i.head == 'serine' for i in self.nucleus))
			self.PC_nucleus_absolut = float(sum(i.head == 'choline' for i in self.nucleus))
			self.PE_nucleus_absolut = float(sum(i.head == 'ethanolamine' for i in self.nucleus))
			self.CL_nucleus_absolut = float(sum(i.head == 'neutral' for i in self.nucleus))

			self.PI_nucleus_relative = self.PI_nucleus_absolut / len(self.nucleus)
			self.PS_nucleus_relative = self.PS_nucleus_absolut / len(self.nucleus)
			self.PC_nucleus_relative = self.PC_nucleus_absolut / len(self.nucleus)
			self.PE_nucleus_relative = self.PE_nucleus_absolut / len(self.nucleus)
			self.CL_nucleus_relative = self.CL_nucleus_absolut / len(self.nucleus)

			self.nucleus_relatives = [self.PI_nucleus_relative, self.PS_nucleus_relative, self.PC_nucleus_relative,\
										self.PE_nucleus_relative, self.CL_nucleus_relative]

			self.nucleus_comp = dict(zip(self.membrane_lipids, self.nucleus_relatives))


	def peroxisomes_composition(self):
		if len(self.peroxisomes) >0:
			self.PI_peroxisomes_absolut = float(sum(i.head == 'inositol' for i in self.peroxisomes))
			self.PS_peroxisomes_absolut = float(sum(i.head == 'serine' for i in self.peroxisomes))
			self.PC_peroxisomes_absolut = float(sum(i.head == 'choline' for i in self.peroxisomes))
			self.PE_peroxisomes_absolut = float(sum(i.head == 'ethanolamine' for i in self.peroxisomes))
			self.CL_peroxisomes_absolut = float(sum(i.head == 'neutral' for i in self.peroxisomes))

			self.PI_peroxisomes_relative = self.PI_peroxisomes_absolut / len(self.peroxisomes)
			self.PS_peroxisomes_relative = self.PS_peroxisomes_absolut / len(self.peroxisomes)
			self.PC_peroxisomes_relative = self.PC_peroxisomes_absolut / len(self.peroxisomes)
			self.PE_peroxisomes_relative = self.PE_peroxisomes_absolut / len(self.peroxisomes)
			self.CL_peroxisomes_relative = self.CL_peroxisomes_absolut / len(self.peroxisomes)

			self.peroxisomes_relatives = [self.PI_peroxisomes_relative, self.PS_peroxisomes_relative,\
											self.PC_peroxisomes_relative, self.PE_peroxisomes_relative,\
											self.CL_peroxisomes_relative]

			self.peroxisomes_comp = dict(zip(self.membrane_lipids, self.peroxisomes_relatives))


	def light_microsomes_composition(self):
		if len(self.light_microsomes) >0:
			self.PI_light_microsomes_absolut = float(sum(i.head == 'inositol' for i in self.light_microsomes))
			self.PS_light_microsomes_absolut = float(sum(i.head == 'serine' for i in self.light_microsomes))
			self.PC_light_microsomes_absolut = float(sum(i.head == 'choline' for i in self.light_microsomes))
			self.PE_light_microsomes_absolut = float(sum(i.head == 'ethanolamine' for i in self.light_microsomes))
			self.CL_light_microsomes_absolut = float(sum(i.head == 'neutral' for i in self.light_microsomes))

			self.PI_light_microsomes_relative = self.PI_light_microsomes_absolut / len(self.light_microsomes)
			self.PS_light_microsomes_relative = self.PS_light_microsomes_absolut / len(self.light_microsomes)
			self.PC_light_microsomes_relative = self.PC_light_microsomes_absolut / len(self.light_microsomes)
			self.PE_light_microsomes_relative = self.PE_light_microsomes_absolut / len(self.light_microsomes)
			self.CL_light_microsomes_relative = self.CL_light_microsomes_absolut / len(self.light_microsomes)

			self.light_microsomes_relatives = [self.PI_light_microsomes_relative, self.PS_light_microsomes_relative,\
											self.PC_light_microsomes_relative, self.PE_light_microsomes_relative,\
											self.CL_light_microsomes_relative]

			self.light_microsomes_comp = dict(zip(self.membrane_lipids, self.light_microsomes_relatives))


	def inner_mit_membrane_composition(self):
		if len(self.inner_mit_membrane) >0:
			self.PI_inner_mit_membrane_absolut = float(sum(i.head == 'inositol' for i in self.inner_mit_membrane))
			self.PS_inner_mit_membrane_absolut = float(sum(i.head == 'serine' for i in self.inner_mit_membrane))
			self.PC_inner_mit_membrane_absolut = float(sum(i.head == 'choline' for i in self.inner_mit_membrane))
			self.PE_inner_mit_membrane_absolut = float(sum(i.head == 'ethanolamine' for i in self.inner_mit_membrane))
			self.CL_inner_mit_membrane_absolut = float(sum(i.head == 'neutral' for i in self.inner_mit_membrane))

			self.PI_inner_mit_membrane_relative = self.PI_inner_mit_membrane_absolut / len(self.inner_mit_membrane)
			self.PS_inner_mit_membrane_relative = self.PS_inner_mit_membrane_absolut / len(self.inner_mit_membrane)
			self.PC_inner_mit_membrane_relative = self.PC_inner_mit_membrane_absolut / len(self.inner_mit_membrane)
			self.PE_inner_mit_membrane_relative = self.PE_inner_mit_membrane_absolut / len(self.inner_mit_membrane)
			self.CL_inner_mit_membrane_relative = self.CL_inner_mit_membrane_absolut / len(self.inner_mit_membrane)

			self.inner_mit_membrane_relatives = [self.PI_inner_mit_membrane_relative, self.PS_inner_mit_membrane_relative,\
											self.PC_inner_mit_membrane_relative, self.PE_inner_mit_membrane_relative,\
											self.CL_inner_mit_membrane_relative]

			self.inner_mit_membrane_comp = dict(zip(self.membrane_lipids, self.inner_mit_membrane_relatives))


	def outer_mit_membrane_composition(self):
		if len(self.outer_mit_membrane) >0:
			self.PI_outer_mit_membrane_absolut = float(sum(i.head == 'inositol' for i in self.outer_mit_membrane))
			self.PS_outer_mit_membrane_absolut = float(sum(i.head == 'serine' for i in self.outer_mit_membrane))
			self.PC_outer_mit_membrane_absolut = float(sum(i.head == 'choline' for i in self.outer_mit_membrane))
			self.PE_outer_mit_membrane_absolut = float(sum(i.head == 'ethanolamine' for i in self.outer_mit_membrane))
			self.CL_outer_mit_membrane_absolut = float(sum(i.head == 'neutral' for i in self.outer_mit_membrane))

			self.PI_outer_mit_membrane_relative = self.PI_outer_mit_membrane_absolut / len(self.outer_mit_membrane)
			self.PS_outer_mit_membrane_relative = self.PS_outer_mit_membrane_absolut / len(self.outer_mit_membrane)
			self.PC_outer_mit_membrane_relative = self.PC_outer_mit_membrane_absolut / len(self.outer_mit_membrane)
			self.PE_outer_mit_membrane_relative = self.PE_outer_mit_membrane_absolut / len(self.outer_mit_membrane)
			self.CL_outer_mit_membrane_relative = self.CL_outer_mit_membrane_absolut / len(self.outer_mit_membrane)

			self.outer_mit_membrane_relatives = [self.PI_outer_mit_membrane_relative, self.PS_outer_mit_membrane_relative,\
											self.PC_outer_mit_membrane_relative, self.PE_outer_mit_membrane_relative,\
											self.CL_outer_mit_membrane_relative]

			self.outer_mit_membrane_comp = dict(zip(self.membrane_lipids, self.outer_mit_membrane_relatives))


	def numbers(self):							
		#for plotting the production of the lipids
		self.number_acetyl_coa.append(self.acetyl_coa_number)
		self.number_acyl_coa.append(len(self.acyl_coa_list))
		self.number_pa.append(len(self.PA_list))
		self.number_cdp_dg.append(len(self.CDP_DG_list))
		self.number_tag.append(len(self.TAG_list))
		self.number_PS.append(len(self.PS_list))
		self.number_PI.append(len(self.PI_list))
		self.number_PE.append(len(self.PE_list))
		self.number_PC.append(len(self.PC_list))
		self.number_CL.append(len(self.CL_list))

		#for plotting the number of lipids in a certain membrane
		self.number_plasma_membrane.append(len(self.plasma_membrane))
		self.number_secretory_vesicles.append(len(self.secretory_vesicles))
		self.number_vacuoles.append(len(self.vacuoles))
		self.number_nucleus.append(len(self.nucleus))
		self.number_peroxisomes.append(len(self.peroxisomes))
		self.number_light_microsomes.append(len(self.light_microsomes))
		self.number_inner_mit_membrane.append(len(self.inner_mit_membrane))
		self.number_outer_mit_membrane.append(len(self.outer_mit_membrane))
		self.number_lipid_droplets.append(len(self.lipid_droplets))

