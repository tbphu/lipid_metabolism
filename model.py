# -*- coding: utf-8 -*-
"""
Created on Wed June 03 2015

@author: Vera
"""
import matplotlib.pyplot as mat
import random

class lipids(object):
	"""
	general class for all kinds of lipids
	with the head groups 'p', 'inositol', 'serine', 'ethanolamine', 'choline', 'neutral', 
	'cdp'(for cdp-dg) and 'None'(for tag)
	possible ffa for sn2: C16:0, C16:1, C18:0, C18:1
	possible ffa for sn1: C16:1 C18:1
	"""
	def __init__(self, head, sn2, sn1):

		self.head_groups = ['p', 'inositol', 'serine', 'ethanolamine', 'choline', 'neutral', 'cdp', 'None']
		self.sn2_options = ['C14:0', 'C16:0', 'C16:1', 'C18:0', 'C18:1', None]	
		self.sn1_options = ['C16:1', 'C18:1', None]						

		self.head = head
		self.sn2 = sn2
		self.sn1 = sn1

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
	At the beginning we set the available number of pyruvate and DHAP in the start function
	acetyl_synthase produces Acetyl-CoA out of Pyruvate by releasing 1 CO2
	In the acyl_synthase reaction we use Acetyl-CoA to synthesise Acyl-CoA --> acyl_coa_list
	These ffa and the available DHAP are used for building the PAs --> PA_list
	PA is further transformed to CDP-DG or TAG --> CDP_DG_synthase
	"""
	def __init__(self):
		self.timesteps = 500			#average value machen mit z.B. 100 durchläufen (standard deviation)
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

		self.chainlength_saturated = {14: 'C14:0', 16: 'C16:0', 18: 'C18:0'}		
		self.chainlength_unsaturated = {16: 'C16:1', 18: 'C18:1'}

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
			self.numbers()
#random die Funktionen hintereinander oder Pool vorher aufteilen und Anteile verteilen oder alle an Gesamtpool, aber Ausführen am Ende
		if len(self.TAG_list) >= 10:
			for i in range(10):
				print 'sn1: ' + str(self.TAG_list[i].sn1)
				print 'sn2: ' + str(self.TAG_list[i].sn2)


	def plot(self):
		fig = mat.figure()
		ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
		ax.plot(self.t, self.number_acetyl_coa[:-1], label = 'acetyl_coa')
		ax.plot(self.t, self.number_acyl_coa[:-1], label = 'acyl_coa')
		ax.plot(self.t, self.number_pa[:-1], label = 'pa')
		ax.plot(self.t, self.number_cdp_dg[:-1], label = 'cdp-dg')
		ax.plot(self.t, self.number_tag[:-1], label = 'tag')
		ax.plot(self.t, self.number_PS[:-1], label = 'ps')
		ax.plot(self.t, self.number_PI[:-1], label = 'pi')
		ax.plot(self.t, self.number_PE[:-1], label = 'pe')
		ax.plot(self.t, self.number_PC[:-1], label = 'pc')
		ax.plot(self.t, self.number_CL[:-1], label = 'cl')
		ax.legend(bbox_to_anchor = (1.05, 1), loc = 2, borderaxespad = 0.)
		mat.show()


	def inositol_synthesis(self):
		for i in range(random.randint(0,3)):
			if self.glucose_6_p_number > 0:
				self.inositol_number += 1
				self.glucose_6_p_number -= 1
				self.p_counter -= 1


	def acetyl_coa_synthase(self):
		for i in range(random.randint(0,10)):			
			if self.pyruvate_number >= 1:			# transformation from pyruvate to acetyl_coa
				self.acetyl_coa_number += 1
				self.pyruvate_number -= 1				
				self.co2_counter += 1


	def acyl_synthase(self):
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
		for i in range(0,20):
			self.lyso_PA_synthase()
			self.PA_synthase()


	def lyso_PA_synthase(self):
		x = random.random()
		if x < 0.9 and 1 in [self.acyl_coa_list[z].saturation for z in range(len(self.acyl_coa_list))]: 	#at least 1 ffa has to be unsaturated 
			if self.dhap_number > 0 and len(self.acyl_coa_list) > 1:
				sn1_chain = random.randint(0, (len(self.acyl_coa_list)-2))
				if self.acyl_coa_list[sn1_chain].saturation == 1:
					chainlength_sn1 = self.acyl_coa_list[sn1_chain].C
					lyso_pa = lipids('p', None, self.chainlength_unsaturated[chainlength_sn1])
					self.lyso_pa_list.append(lyso_pa)
					self.dhap_number -= 1
					del self.acyl_coa_list[sn1_chain]
				else:
					self.lyso_PA_synthase()


	def PA_synthase(self):	
		x = random.random()
		if x < 0.9: 		
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
		x = random.random()
		if x <= 0.15 and len(self.PS_list) >= 1:
			self.PS_list[0].head = 'ethanolamine'				#PE synthesis from PS
			self.PE_list.append(self.PS_list[0])
			del self.PS_list[0]
			self.co2_counter += 1

		

	def PC_synthase(self):
		x = random.random()
		if x <= 0.1 and len(self.PE_list) >= 1 and self.SAM_number >= 3:
			self.PE_list[0].head = 'choline'					#PC synthesis from PE
			self.PC_list.append(self.PE_list[0])
			del self.PE_list[0]
			self.SAM_number -= 3
			self.SAH_number += 3


	def CL_synthase(self):
		x = random.random()
		if x <= 0.1 and self.glycerol_3_p_mito_number > 0 and len(self.CDP_DG_list) >= 2:
			self.CDP_DG_list[0].head = 'neutral'
			self.CL_list.append(self.CDP_DG_list[0])
			self.CL_list[-1].sn4, self.CL_list[-1].sn3 = self.CDP_DG_list[1].sn2, self.CDP_DG_list[1].sn1
			del self.CDP_DG_list[0:2]
			self.glycerol_3_p_mito_number -= 1


	def numbers(self):											#collecting the number of the particular products
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

