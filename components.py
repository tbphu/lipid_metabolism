"""
Component classes for the lipid metabolism model. 

Created/Started on Wed June 03 2015
@author: Vera Schuetzhold - vera.schue@gmail.com
"""
from numpy.random import choice

class lipids(object):
	"""
	Lipid class
	===========
	Main class for all lipid components in the model, phopholipids, TAG & CL.

	Head groups: p, inositol, serine, ethanolamine, choline, neutral, cdp (only cdp-dg) and None (only tag)
	Ffa for sn2: C16:1, C18:1
	Ffa for sn1: C16:0, C16:1, C18:0, C18:1
	"""
	def __init__(self, head, sn2, sn1, comp, comp_weights):

		self.head_groups = ['p', 'inositol', 'serine', 'ethanolamine', 'choline', 'neutral', 'cdp', 'sterol', 'ceramide', None]
		self.sn2_options = ['C16:1', 'C18:1', None]	
		self.sn1_options = ['C16:0', 'C16:1', 'C18:0', 'C18:1', None]
		self.compartment_options = ['plasma_membrane', 'secretory_vesicles', 'vacuoles', 'nucleus', 'peroxisomes', 'light_microsomes',\
							'inner_mit_membrane', 'outer_mit_membrane', 'lipid_droplets', None]	
		self.compartment = ['plasma_membrane', 'secretory_vesicles', 'vacuoles', 'nucleus', 'peroxisomes', 'light_microsomes',\
							'inner_mit_membrane', 'outer_mit_membrane']	
		

		self.plasma_membrane_comp = {'PS': 0.17320, 'PI': 0.09124, 'PC': 0.08660, 'PE': 0.10464, 'CL': 0.00103, 'PA': 0.02010}
		self.secretory_vesicles_comp = {'PS': 0.08205, 'PI': 0.11745, 'PC': 0.20824, 'PE': 0.13573, 'CL': 0.01239, 'PA': 0.01525}
		self.vacuoles_comp = {'PS': 0.04817, 'PI': 0.16604, 'PC': 0.40517, 'PE': 0.17537, 'CL': 0.02442, 'PA': 0.02866}
		self.nucleus_comp = {'PS': 0.04038, 'PI': 0.09650, 'PC': 0.27645, 'PE': 0.16848, 'CL': 0.01049, 'PA': 0.01781}
		self.peroxisomes_comp = {'PS': 0.03235, 'PI': 0.11360, 'PC': 0.34656, 'PE': 0.16465, 'CL': 0.05033, 'PA': 0.01150}
		self.light_microsomes_comp = {'PS': 0.05304, 'PI': 0.06019, 'PC': 0.40796, 'PE': 0.26583, 'CL': 0.00381, 'PA': 0.00222}
		self.inner_mit_membrane_comp = {'PS': 0.02880, 'PI': 0.12273, 'PC': 0.29107, 'PE': 0.18192, 'CL': 0.12204, 'PA': 0.01137}
		self.outer_mit_membrane_comp = {'PS': 0.01189, 'PI': 0.10108, 'PC': 0.45190, 'PE': 0.32307, 'CL': 0.05847, 'PA': 0.04360}
		self.lipid_droplets_comp = {'PS': 0.0, 'PI': 0.0, 'PC': 0.0, 'PE': 0.0, 'CL': 0.0, 'PA':0.0}
		self.membranes_comp = [self.plasma_membrane_comp, self.secretory_vesicles_comp, self.vacuoles_comp, self.nucleus_comp,\
								self.peroxisomes_comp, self.light_microsomes_comp, self.inner_mit_membrane_comp,\
								self.outer_mit_membrane_comp, self.lipid_droplets_comp]		

		self.head = head
		self.sn2 = sn2
		self.sn1 = sn1
		self.comp = comp
		self.comp_weights = comp_weights

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
			weights = [self.comp_weights[i] * self.membranes_comp[i]['PS'] for i in range(len(self.comp_weights))]
			weights_normal = [weight / sum(weights) for weight in weights]
		elif self.head == 'inositol':
			weights = [self.comp_weights[i] * self.membranes_comp[i]['PI'] for i in range(len(self.comp_weights))]
			weights_normal = [weight / sum(weights) for weight in weights]
		elif self.head == 'choline':
			weights = [self.comp_weights[i] * self.membranes_comp[i]['PC'] for i in range(len(self.comp_weights))]
			weights_normal = [weight / sum(weights) for weight in weights]
		elif self.head == 'ethanolamine':
			weights = [self.comp_weights[i] * self.membranes_comp[i]['PE'] for i in range(len(self.comp_weights))]
			weights_normal = [weight / sum(weights) for weight in weights]
		elif self.head == 'p':
			weights = [self.comp_weights[i] * self.membranes_comp[i]['PA'] for i in range(len(self.comp_weights))]
			weights_normal = [weight / sum(weights) for weight in weights]
		self.comp = choice(self.compartment, p = weights_normal)

class TAG(lipids):
	def __init__(self, head, sn3, sn2, sn1, comp, comp_weights):
		super(TAG, self).__init__(head, sn2, sn1, comp, comp_weights)
		sn3 = None

	def comp_choice(self):
		self.comp = 'lipid_droplets'


class CL(lipids):
	def __init__(self, head, sn2, sn1, sn4, sn3, comp, comp_weights):
		super(CL, self).__init__(head, sn2, sn1, comp, comp_weights)
		sn4 = None
		sn3 = None	

	def comp_choice(self):
		weights = [self.comp_weights[i] * self.membranes_comp[i]['CL'] for i in range(len(self.comp_weights))]
		weights_normal = [weight / sum(weights) for weight in weights]
		self.comp = choice(self.compartment, p = weights_normal)


class fatty_acids(object):	#name als attribut statt der einzelnen unterklassen
	"""
	Fatty acid class
	================

	Attribute C: number of C-Atoms
	Attribute saturation: 0 = saturated, 1 = unsaturated
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


class sterol(object):
	"""
	Sterol class
	============
	Head options: sterol, sterylester
	"""
	def __init__(self, head, comp, comp_weights):
		self.head_options = ['sterol', 'sterylester']
		self.compartment_options = ['plasma_membrane', 'secretory_vesicles', 'vacuoles', 'nucleus', 'peroxisomes', 'light_microsomes',\
							'inner_mit_membrane', 'outer_mit_membrane', 'lipid_droplets', None]	
		self.compartment = ['plasma_membrane', 'secretory_vesicles', 'vacuoles', 'nucleus', 'peroxisomes', 'light_microsomes',\
							'inner_mit_membrane', 'outer_mit_membrane']	
		self.compartment_weights = [0.665, 0.01, 0.155, 0.10, 0.005, 0.01, 0.04, 0.015]	
		self.plasma_membrane_comp = {'ES': 0.48454}
		self.secretory_vesicles_comp = {'ES': 0.42900}
		self.vacuoles_comp = {'ES': 0.15200}
		self.nucleus_comp = {'ES': 0.390}
		self.peroxisomes_comp = {'ES': 0.281}
		self.light_microsomes_comp = {'ES': 0.206}
		self.inner_mit_membrane_comp = {'ES': 0.242}
		self.outer_mit_membrane_comp = {'ES': 0.009}
		self.lipid_droplets_comp = {'ES': 0.0}
		self.membranes_comp = [self.plasma_membrane_comp, self.secretory_vesicles_comp, self.vacuoles_comp, self.nucleus_comp,\
								self.peroxisomes_comp, self.light_microsomes_comp, self.inner_mit_membrane_comp,\
								self.outer_mit_membrane_comp, self.lipid_droplets_comp]

		self.head = head
		self.comp = comp
		self.comp_weights = comp_weights

	def comp_choice(self):	
		weights = [self.comp_weights[i] * self.membranes_comp[i]['ES'] for i in range(len(self.comp_weights))]
		if sum(weights) != 0:
			weights_normal = [weight / sum(weights) for weight in weights]
			self.comp = choice(self.compartment, p = weights_normal)

	@property
	def head(self):
		return self.__head
	@head.setter
	def head(self, group):
		if group not in self.head_options:
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


class sterylester(object):
	"""
	Sterylester class
	=================
	Head options: sterylester
	"""
	def __init__(self, head, FA, comp, comp_weights):
		self.head_options = ['sterylester']
		self.compartment_options = ['lipid_droplets', None]
		self.head = head
		self.FA = FA
		self.comp = None
		self.comp_weights = comp_weights
		
	def comp_choice(self):
		self.comp = 'lipid_droplets'

	@property
	def head(self):
		return self.__head
	@head.setter
	def head(self, group):
		if group not in self.head_options:
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


class sphingolipid(object):
	"""
	Sphingolipid class
	==================
	Head options: ceramide
	"""
	def __init__(self, head, comp, comp_weights):
		self.head_options = ['ceramide']
		self.compartment_options = ['plasma_membrane', 'secretory_vesicles', 'vacuoles', 'nucleus', 'peroxisomes', 'light_microsomes',\
							'inner_mit_membrane', 'outer_mit_membrane', None]	
		self.compartment = ['plasma_membrane', 'secretory_vesicles', 'vacuoles', 'nucleus', 'peroxisomes', 'light_microsomes',\
							'inner_mit_membrane', 'outer_mit_membrane']	
		self.compartment_weights = [0.665, 0.01, 0.155, 0.10, 0.005, 0.01, 0.04, 0.015]		
		self.plasma_membrane_comp = {'SL': 0.03557}
		self.secretory_vesicles_comp = {'SL': 0.05029}
		self.vacuoles_comp = {'SL': 0.06525}
		self.nucleus_comp = {'SL': 0.02622}
		self.peroxisomes_comp = {'SL': 0.0}
		self.light_microsomes_comp = {'SL': 0.00397}
		self.inner_mit_membrane_comp = {'SL': 0.0}
		self.outer_mit_membrane_comp = {'SL': 0.0}
		self.lipid_droplets_comp = {'SL': 0.0}
		self.membranes_comp = [self.plasma_membrane_comp, self.secretory_vesicles_comp, self.vacuoles_comp, self.nucleus_comp,\
								self.peroxisomes_comp, self.light_microsomes_comp, self.inner_mit_membrane_comp,\
								self.outer_mit_membrane_comp, self.lipid_droplets_comp]
	
		self.head = head
		self.comp = comp
		self.comp_weights = comp_weights

	def comp_choice(self):	
		weights = [self.comp_weights[i] * self.membranes_comp[i]['SL'] for i in range(len(self.comp_weights))]
		if sum(weights) != 0:
			weights_normal = [weight / sum(weights) for weight in weights]
			self.comp = choice(self.compartment, p = weights_normal)

	@property
	def head(self):
		return self.__head
	@head.setter
	def head(self, group):
		if group not in self.head_options:
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


