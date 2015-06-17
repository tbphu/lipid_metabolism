# -*- coding: utf-8 -*-
"""
Created/Started on Wed June 03 2015

@author: Vera
"""
import matplotlib.pyplot as mat
import random

class lipids(object):
	"""
	general class for all kinds of lipids
	with the head groups 'p', 'inositol', 'serine', 'ethanolamine', 'choline', 'neutral', 
	'cdp'(for cdp-dg) and 'None'(for tag)
	possible ffa for sn2: C14:0, C16:0, C16:1, C18:0, C18:1
	possible ffa for sn1: C16:1 C18:1
	"""
	def __init__(self, head, sn2, sn1, comp):

		self.head_groups = ['p', 'inositol', 'serine', 'ethanolamine', 'choline', 'neutral', 'cdp', 'None']
		self.sn2_options = ['C14:0', 'C16:0', 'C16:1', 'C18:0', 'C18:1', None]	
		self.sn1_options = ['C16:1', 'C18:1', None]
		self.compartment_options = ['plasma_membrane', 'secretory_vesicles', 'vacuoles', 'nucleus', 'peroxisomes', 'light_microsomes',\
							'inner_mit_membrane', 'outer_mit_membrane', None]					

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


class CL(lipids):

	def __init__(self, head, sn2, sn1, sn4, sn3):
		super(CL, self).__init__(head, sn2, sn1)
		sn4 = None
		sn3 = None	


class fatty_acids(object):	#name als attribut statt der einzelnen unterklassen
	"""
	class for all the small molecules that are needed to create the lipids
	subclasses will be pyruvate, acetyl_coa, acyl_coa, dhap
	attribute C: number of C-Atoms
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
	

class model():
	"""
	The model. 
	At the beginning there are several lists defined which will contain the produced lipids. 
	The model starts with a given amount of precursors (pyruvate, dhap, ctp, serine, glucose-6-p, SAM, glycerol-3-p).
	After several reactions and the tranport function in the end there are membranes of different compartments with several different lipids.
	"""
	def __init__(self):
		#determining the timesteps
		self.timesteps = 500			
		self.t = [i for i in range(self.timesteps)]

		# number of available precursors
		self.pyruvate_number = 10000	
		self.dhap_number = 1000	
		self.ctp_number = 1000
		self.serine_number = 1000
		self.glucose_6_p_number = 1000
		self.inositol_number = 0
		self.acetyl_coa_number = 0
		self.SAM_number = 1000
		self.SAH_number = 0
		self.glycerol_3_p_mito_number = 1000
		self.co2_counter = 0
		self.p_counter = 0

		self.compartment = ['plasma_membrane', 'secretory_vesicles', 'vacuoles', 'nucleus', 'peroxisomes', 'light_microsomes',\
							'inner_mit_membrane', 'outer_mit_membrane']	

		#empty lists for the produced lipids
		self.acyl_coa_list = []
		self.lyso_pa_list = []
		self.PA_list = []
		self.CDP_DG_list = []
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

		self.chainlength_saturated = {14: 'C14:0', 16: 'C16:0', 18: 'C18:0'}		
		self.chainlength_unsaturated = {16: 'C16:1', 18: 'C18:1'}

		#functions to run the model
		for t in range(self.timesteps):
			self.inositol_synthesis()
			self.acetyl_coa_synthase()
			self.acyl_synthase()
			self.PA_synthese()
			self.CDP_DG_synthase()
			self.phospholipid_synthase()
			self.PE_synthase()
			self.PC_synthase()
			self.CL_synthase()
			self.transport()
			self.numbers()
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
		ax.legend(bbox_to_anchor = (1.05, 1), loc = 2, borderaxespad = 0.)
		mat.show()


	def inositol_synthesis(self):
		'''
		Synthesis of inositol from glucose-6-p by the Myo-inositol-1-p synthase and the Myo-inositol 1
		phosphatase.
		'''
		for i in range(random.randint(0,3)):
			if self.glucose_6_p_number > 0:
				self.inositol_number += 1
				self.glucose_6_p_number -= 1
				self.p_counter -= 1


	def acetyl_coa_synthase(self):
		'''
		Synthesis of Acetyl-CoA: pyruvate dehydrogenase drives the reaction pyruvate to Acetyl-CoA, CO2 is released
		'''
		for i in range(random.randint(0,10)):			
			if self.pyruvate_number >= 1:			# transformation from pyruvate to acetyl_coa
				self.acetyl_coa_number += 1
				self.pyruvate_number -= 1				
				self.co2_counter += 1


	def acyl_synthase(self):
		'''
		Simplified synthesis of acyl_coa: several Acetyl-CoA are building a fatty acid (C14, C16:0, C16:1, C18:0 or C18:1)
		The intermediate Malonyl-CoA is leaved out.
		'''
		for i in range(random.randint(0,10)):
			x = random.random()						#5 reactions in 1 timestep but only with a probability of 90%
			if self.acetyl_coa_number >= 2:		#control if at least 2 Acetyl-CoA are available
				if len(self.acyl_coa_list) == 0:		#starting the first reaction
					new_acyl = fatty_acids(2, random.randint(0,1))
					self.acyl_coa_list.append(new_acyl)
					self.acyl_coa_list[-1].C += 2
					self.acetyl_coa_number -= 2

				elif self.acyl_coa_list[-1].C >=14 and x <= 0.1:
					self.acyl_coa_list[-1].saturation = 0
					new_acyl = fatty_acids(2, random.randint(0,1))
					self.acyl_coa_list.append(new_acyl)
					self.acyl_coa_list[-1].C += 2
					self.acetyl_coa_number -= 2

				elif self.acyl_coa_list[-1].C >= 16 and x <= 0.45:	#stop the reaction cycle and starting a new one
					new_acyl = fatty_acids(2, random.randint(0,1))
					self.acyl_coa_list.append(new_acyl)
					self.acyl_coa_list[-1].C += 2
					self.acetyl_coa_number -= 2

				elif self.acyl_coa_list[-1].C >= 18:				#stop the reaction cycle and starting a new one
					new_acyl = fatty_acids(2, random.randint(0,1))
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
		for i in range(0,30):
			self.lyso_PA_synthase()
			self.PA_synthase()


	def lyso_PA_synthase(self):
		'''
		Production of Lyso-PA by adding one acyl-coa to DHAP (sn2: always unsaturated) --> DHAP acyltransferase/acyl-DHAP reductase
		'''
		x = random.random()
		if x < 0.95 and 1 in [self.acyl_coa_list[z].saturation for z in range(len(self.acyl_coa_list))]: 	#at least 1 ffa has to be unsaturated 
			if self.dhap_number > 0 and len(self.acyl_coa_list) > 1:
				sn1_chain = random.randint(0, (len(self.acyl_coa_list)-2))
				if self.acyl_coa_list[sn1_chain].saturation == 1:
					chainlength_sn1 = self.acyl_coa_list[sn1_chain].C
					lyso_pa = lipids('p', None, self.chainlength_unsaturated[chainlength_sn1], None)
					self.lyso_pa_list.append(lyso_pa)
					self.dhap_number -= 1
					del self.acyl_coa_list[sn1_chain]
				else:
					self.lyso_PA_synthase()


	def PA_synthase(self):	
		'''
		Synthesis of PA by adding the second fatty acid to DHAP (sn1: saturated or unsaturated) --> 1-acyl-sn-glycerol-3-phosphate acyltransferase
		'''
		x = random.random()
		if x < 0.95: 		
			if len(self.lyso_pa_list) > 0 and len(self.acyl_coa_list) > 1:		# available ffa
				sn2_chain = random.randint(0, (len(self.acyl_coa_list)-2))		
				chainlength_sn2 = self.acyl_coa_list[sn2_chain].C
				if self.acyl_coa_list[sn2_chain].saturation == 0:
					self.lyso_pa_list[0].sn2 = self.chainlength_saturated[chainlength_sn2]	# creating a new lipid: PA
				elif self.acyl_coa_list[sn2_chain].saturation == 1:
					self.lyso_pa_list[0].sn2 = self.chainlength_unsaturated[chainlength_sn2]
				self.PA_list.append(self.lyso_pa_list[0])			
				del self.acyl_coa_list[sn2_chain]		# deletion of the consumed ffa
				del self.lyso_pa_list[0]


	def CDP_DG_synthase(self):
		'''
		PA is processed to CDP-DG (CDP-diacylglycerol synthase), that further reacts to the phospholipids, or to TAG (PA phosphatase/DAG acyltransferase) 
		for fatty acid storage
		'''
		x = random.random()
		if x <= 0.7 and self.ctp_number > 0:
			if len(self.PA_list) > 0:
				self.PA_list[0].head = 'cdp'
				self.CDP_DG_list.append(self.PA_list[0])		#CDP-DG production from PA
				del self.PA_list[0]
				self.ctp_number -= 1
				self.p_counter -= 2

		elif x <= 0.8:
			if len(self.acyl_coa_list) > 0 and len(self.PA_list) > 0:
				self.PA_list[0].head = 'None'
				self.TAG_list.append(self.PA_list[0])			#TAG production from PA
				del self.PA_list[0]
				del self.acyl_coa_list[0]
				self.p_counter -= 1


	def phospholipid_synthase(self):
		'''
		CDP-DG is processed to PS (PS synthase) or PI (PI synthase) as phospholipids.
		'''
		x = random.random()
		if len(self.CDP_DG_list) >= 1 and self.serine_number > 0:
			if x <= 0.4:
				self.CDP_DG_list[0].head = 'serine'				#PS synthesis from CDP-DG
				self.PS_list.append(self.CDP_DG_list[0])
				del self.CDP_DG_list[0]
				self.serine_number -= 1

			elif x <= 0.6 and self.inositol_number > 0:
				self.CDP_DG_list[0].head = 'inositol'			#PI synthesis from CDP-DG
				self.PI_list.append(self.CDP_DG_list[0])
				del self.CDP_DG_list[0]
				self.inositol_number -= 1


	def PE_synthase(self):
		'''
		PE is derived from PS by releasing 1 CO2 --> PS decarboxylase.
		'''
		x = random.random()
		if x <= 0.15 and len(self.PS_list) >= 1:
			self.PS_list[0].head = 'ethanolamine'				#PE synthesis from PS
			self.PE_list.append(self.PS_list[0])
			del self.PS_list[0]
			self.co2_counter += 1

		

	def PC_synthase(self):
		'''
		PC is derived from PE. As enzymes serve 3 methyltransferases which need SAM and produce SAH as a side product.
		'''
		x = random.random()
		if x <= 0.1 and len(self.PE_list) >= 1 and self.SAM_number >= 3:
			self.PE_list[0].head = 'choline'								#PC synthesis from PE
			self.PC_list.append(self.PE_list[0])
			del self.PE_list[0]
			self.SAM_number -= 3
			self.SAH_number += 3


	def CL_synthase(self):
		'''
		Synthesis of cardiolipin, for which 2 CDP-DG are needed. Different enzymes are needed.
		'''
		x = random.random()
		if x <= 0.1 and self.glycerol_3_p_mito_number > 0 and len(self.CDP_DG_list) >= 2:
			self.CDP_DG_list[0].head = 'neutral'
			self.CL_list.append(self.CDP_DG_list[0])
			self.CL_list[-1].sn4, self.CL_list[-1].sn3 = self.CDP_DG_list[1].sn2, self.CDP_DG_list[1].sn1
			del self.CDP_DG_list[0:2]
			self.glycerol_3_p_mito_number -= 1


	def transport(self):
		'''
		General transport function for all produced lipids.
		'''
		self.transport_PS()
		self.transport_PI()
		self.transport_PE()
		self.transport_PC()
		self.transport_CL()


	def transport_PS(self):
		'''
		Transport of PS to the different compartiment membranes. The distribution to the different membranes is random.
		'''
		if len(self.PS_list) > 3:
			self.PS_list[0].comp = self.compartment[random.randint(0, 7)]
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
		Transport of PI to the different compartiment membranes. The distribution to the different membranes is random.
		'''
		if len(self.PI_list) > 3:
			self.PI_list[0].comp = self.compartment[random.randint(0, 7)]
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
		Transport of PE to the different compartiment membranes. The distribution to the different membranes is random.
		'''
		if len(self.PE_list) > 3:
			self.PE_list[0].comp = self.compartment[random.randint(0, 7)]
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
		Transport of PC to the different compartiment membranes. The distribution to the different membranes is random.
		'''
		if len(self.PC_list) > 3:
			self.PC_list[0].comp = self.compartment[random.randint(0, 7)]
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
		Transport of CL to the different compartiment membranes. The distribution to the different membranes is random.
		'''
		i = 0
		if len(self.CL_list) > 3:
			self.CL_list[i].comp = self.compartment[random.randint(0, 7)]
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

