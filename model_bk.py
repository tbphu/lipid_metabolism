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
    possible ffa for sn2: C16:1, C18:1
    possible ffa for sn1: C16:0, C16:1, C18:0, C18:1
    """
    def __init__(self, head, sn2, sn1, comp, comp_weights):

        self.head_groups = ['p', 'inositol', 'serine', 'ethanolamine', 'choline', 'neutral', 'cdp', 'sterol', None, 'ceramide']
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


class sterol(object):

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



class model:
    """
    The model.
    At the beginning there are several lists defined which will contain the produced lipids.
    The model starts with a given amount of precursors (pyruvate, dhap, ctp, serine, glucose-6-p, SAM, glycerol-3-p).
    After several reactions and the tranport function in the end there are membranes of different compartments with several different lipids.
    """
    def __init__(self):
        self.timesteps = 10
        self.volume = 35

        self.rates = {'glycerol_3_p_synthesis': 8, 'inositol_synthesis': 5, 'ceramide_synthesis': 2, 'acetyl_coa_synthase': 650, 'acyl_synthase': 450, 'PA_synthese': 17, \
                        'CDP_DG_synthase': 20, 'TAG_synthese': 30, 'TAG_lipase': 23, 'DAG_kinase': 40, 'PS_synthase': 18, 'PI_synthase': 6,\
                        'PE_synthase': 12, 'PC_synthase': 5, 'CL_synthase': 2, 'Ergosterol_synthase': 25, 'Sterylester_synthase': 25, 'Sphingolipid_synthase': 2}
        self.probability = {'glycerol_3_p_synthesis': 0.5, 'inositol_synthesis': 0.5, 'ceramide_synthesis': 0.5, 'acetyl_coa_synthase': 0.8, 'acyl_synthase': 0.5, \
                            'acyl_synthase_C16': 0.375, 'acyl_synthase_C18': 0.998, 'lyso_PA_synthase': 0.45, 'PA_synthase': 0.2, 'CDP_DG_synthase': 0.8, \
                            'DAG_synthase': 0.01, 'TAG_synthase': 0.2, 'TAG_lipase': 0.8, 'DAG_kinase': 0.1, 'PS_synthase': 0.5, 'PI_synthase': 0.5, \
                            'PE_synthase': 0.5, 'PC_synthase': 0.5, 'CL_synthase': 0.05, 'Ergosterol_synthase': 0.6, 'Sterylester_synthase': 0.4, 'Sphingolipid_synthase': 0.2}

        self.probability_G1 = {'DAG_synthase': 0.3, 'TAG_synthase': 0.2, 'TAG_lipase': 0.05, 'DAG_kinase': 0.03}

        self.probability_S_M = {'DAG_synthase': 0.01, 'TAG_synthase': 0.2, 'TAG_lipase': 0.6, 'DAG_kinase': 0.1, 'Sterylester_synthase': 0.2}

        self.manual_threshold = {reaction: 1-prob for reaction, prob in self.probability.iteritems()}

        self.compartment_weights = [0.67, 0.01, 0.155, 0.101, 0.007, 0.007, 0.03, 0.015]

        self.weights_fa = [0.4, 0.6]

    def run(self):
        #determining the timesteps, self.t for plotting, self.time for cell cycle

        self.time = 0
        self.t = [i for i in range(self.timesteps)]


        #number of small molecules in the cell
        self.precursors_dict = {'pyruvate' : 2200., 'acetyl_coa': 1000, 'glycerol-3-p': 1000., 'DHAP': 1000., 'serine': 250., 'glucose_6_p': 1000., 'SAM': 331., 'SAH': 0.,\
                                'glycerol_3_p_mito': 50., 'ceramide': 100, 'GDP-mannose': 0, 'NAD': 0, 'NADH': 0, 'NADP': 0, 'NADPH': 0, 'O2': 0, 'H2O': 0,\
                                'CO2': 0, 'Pi': 0, 'CTP': 1000, 'CMP': 0, 'inositol': 350, 'ATP': 0, 'ADP': 0}


        #number of small molecules that is produced from anywhere in the cell and will be added every 10 seconds
        self.precursors_production_G1 = {'pyruvate' : 1500., 'acetyl_coa': 0, 'glycerol-3-p': 5., 'DHAP': 30., 'serine': 20., 'glucose_6_p': 8., 'SAM': 45., 'SAH': 0.,\
                                        'glycerol_3_p_mito': 5., 'ceramide': 0, 'GDP-mannose': 10, 'NAD': 0, 'NADH': 0, 'NADP': 0, 'NADPH': 0, 'O2': 0, 'H2O': 0,\
                                        'CO2': 0, 'Pi': 0, 'CTP': 20, 'CMP': 0, 'inositol': 0, 'ATP': 0, 'ADP': 0}

        self.precursors_production_S_M = {'pyruvate' : 1700., 'acetyl_coa': 0, 'glycerol-3-p': 5., 'DHAP': 35., 'serine': 30., 'glucose_6_p': 12., 'SAM': 55., 'SAH': 0.,\
                                        'glycerol_3_p_mito': 5., 'ceramide': 0, 'GDP-mannose': 10, 'NAD': 0, 'NADH': 0, 'NADP': 0, 'NADPH': 0, 'O2': 0, 'H2O': 0,\
                                        'CO2': 0, 'Pi': 0, 'CTP': 45, 'CMP': 0, 'inositol': 0, 'ATP': 0, 'ADP': 0}



        #lists of all the precursors for plotting
        self.precursors = {'pyruvate' : [], 'acetyl_coa': [], 'glycerol-3-p': [], 'DHAP': [], 'serine': [], 'glucose_6_p': [], 'SAM': [], 'SAH': [],\
                            'glycerol_3_p_mito': [], 'ceramide': [], 'GDP-mannose': [], 'NAD': [], 'NADH': [], 'NADP': [], 'NADPH': [], 'O2': [], 'H2O': [], 'CO2': [],\
                            'Pi': [], 'CTP': [], 'CMP': [], 'inositol': [], 'ATP': [], 'ADP': []}

        #names of the small molecules for plotting
        self.precursor_keys = ['pyruvate', 'acetyl_coa', 'glycerol-3-p', 'DHAP', 'serine', 'glucose_6_p', 'SAM', 'SAH', 'glycerol_3_p_mito', 'ceramide',\
                                'GDP-mannose', 'NAD', 'NADH', 'NADP', 'NADPH', 'O2', 'H2O', 'CO2', 'Pi', 'CTP', 'CMP', 'inositol', 'ATP', 'ADP']


        #list of the 4 cell cycle phases
        self.cell_cycle_phases = ['G1', 'S', 'G2', 'M']

        #empty lists for the produced fatty acids and lipids
        self.acyl_coa_list = []
        self.acyl_coa_list_C26 = []
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
        self.Sterylester_list = []
        self.Sphingolipid_list = []

        #lists to collect the lists of all produced species and membrane lipids (not transported yet) for plotting
        self.precursor_list = [self.acyl_coa_list, self.PA_list, self.CDP_DG_list, self.TAG_list, self.PS_list, self.PI_list,\
                                self.PE_list, self.PC_list, self.CL_list, self.Ergosterol_list, self.Sterylester_list, self.DAG_list, self.Sphingolipid_list]

        self.lipid_lists = [self.PS_list, self.PI_list, self.PC_list, self.PE_list, self.CL_list, self.PA_list, self.Ergosterol_list, self.Sterylester_list, self.TAG_list, self.Sphingolipid_list]

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

        #relatives list for calculation of the relative parts of each lipid in a membrane for the compartment_relatives_dict
        self.relatives_list = []

        #collecting the products of every timestep. Lipids that are produced in the start function don't need the 0 for plotting
        self.number_acetyl_coa = []
        self.number_acyl_coa = []
        self.number_pa = []
        self.number_cdp_dg = []
        self.number_tag = []
        self.number_PS = []
        self.number_PI = []
        self.number_PE = []
        self.number_PC = []
        self.number_CL = []
        self.number_Ergosterol = []
        self.number_Sterylester = []
        self.number_DAG = []
        self.number_Sphingolipid = []

        self.number_lipids_list = [self.number_acyl_coa, self.number_pa, self.number_cdp_dg, self.number_tag,\
                                    self.number_PS, self.number_PI, self.number_PE, self.number_PC,	self.number_CL,\
                                    self.number_Ergosterol, self.number_Sterylester, self.number_DAG, self.number_Sphingolipid]

        #counting the lipids in each membrane after every timestep
        self.number_plasma_membrane = []
        self.number_secretory_vesicles = []
        self.number_vacuoles = []
        self.number_nucleus = []
        self.number_peroxisomes = []
        self.number_light_microsomes = []
        self.number_inner_mit_membrane = []
        self.number_outer_mit_membrane = []
        self.number_lipid_droplets = []

        self.number_membranes_list = [self.number_plasma_membrane, self.number_secretory_vesicles, self.number_vacuoles, \
                                        self.number_nucleus, self.number_peroxisomes, self.number_light_microsomes,\
                                        self.number_inner_mit_membrane,	self.number_outer_mit_membrane, self.number_lipid_droplets]

        #the possible fatty acids and weights to reach the biological proportions
        self.chainlength_saturated = {16: 'C16:0', 18: 'C18:0'}
        self.chainlength_unsaturated = {16: 'C16:1', 18: 'C18:1'}
        self.unsaturated_weights = [0.375, 0.625]
        self.chainlength_saturated_unsaturated = ['C16:0', 'C18:0', 'C16:1', 'C18:1']
        self.saturation_weights_total = [0.2, 0.2, 0.167, 0.433]

        #names of the membrane lipids
        self.membrane_lipids = ['PS', 'PI', 'PC', 'PE', 'CL', 'PA', 'ES', 'SE', 'TAG', 'SL']

        self.compartment_relatives_dict = {comp: dict(zip(self.membrane_lipids, [0.0 for z in range(10)])) for comp in self.compartment}

        self.Km = {}

        self.start()	#function that produces the lipids and membranes that are existing at the beginning of the cell cycle
        self.Km_calculation()
        for t in range(self.timesteps):
            self.time += 1 		#counting the seconds for knowing the cell cycle phase

            #function that can change the cell cycle phase
            self.cell_cycle()

            if self.phase == 'G1':
                self.precursors_production = self.precursors_production_G1
            else:
                self.precursors_production = self.precursors_production_S_M

            if self.time % 10 == 0:
                for key in self.precursors_dict:
                    self.precursors_dict[key] += self.precursors_production[key]

            self.thresholds = {'glycerol_3_p_synthesis': 1- (self.precursors_dict['DHAP'] / (self.Km['glycerol_3_p_synthesis']['DHAP'] + self.precursors_dict['DHAP'])), \

                                'inositol_synthesis': 1- (self.precursors_dict['glucose_6_p'] / (self.Km['inositol_synthesis']['glucose_6_p'] + self.precursors_dict['glucose_6_p'])), \

                                'ceramide_synthesis': 1- (self.precursors_dict['serine'] / (self.Km['ceramide_synthesis']['serine'] + self.precursors_dict['serine'])),\

                                'acetyl_coa_synthase': 1- (self.precursors_dict['pyruvate'] / (self.Km['acetyl_coa_synthase']['pyruvate'] + self.precursors_dict['pyruvate'])), \

                                'acyl_synthase': 1- (self.precursors_dict['acetyl_coa'] / (self.Km['acyl_synthase']['acetyl_coa'] + self.precursors_dict['acetyl_coa'])),\

                                'acyl_synthase_C16': 0.625,\

                                'acyl_synthase_C18': 0.002,\

                                'lyso_PA_synthase': 1- (((float(len(self.acyl_coa_list_saturated)) + float(len(self.acyl_coa_list_unsaturated))) / (self.Km['lyso_PA_synthase']['acyl_coa'] \
                                                + (float(len(self.acyl_coa_list_saturated)) + float(len(self.acyl_coa_list_unsaturated))))) \
                                                * (self.precursors_dict['DHAP'] / (self.Km['lyso_PA_synthase']['DHAP'] + self.precursors_dict['DHAP']))),\

                                'PA_synthase': 1- ((float(len(self.lyso_pa_list)) / (self.Km['PA_synthase']['lyso-PA'] + float(len(self.lyso_pa_list)))) \
                                                * ((float(len(self.acyl_coa_list_saturated)) + float(len(self.acyl_coa_list_unsaturated))) / (self.Km['PA_synthase']['acyl_coa'] \
                                                + (float(len(self.acyl_coa_list_saturated)) + float(len(self.acyl_coa_list_unsaturated)))))),\

                                'CDP_DG_synthase': 1- ((float(len(self.PA_list)) / (self.Km['CDP_DG_synthase']['PA'] + float(len(self.PA_list)))) \
                                                    * (self.precursors_dict['CTP'] / (self.Km['CDP_DG_synthase']['CTP'] + self.precursors_dict['CTP']))),\

                                'DAG_synthase': 1- (float(len(self.PA_list)) / (self.Km['DAG_synthase']['PA'] + float(len(self.PA_list)))), \

                                'TAG_synthase': 1- (((float(len(self.acyl_coa_list_saturated)) + float(len(self.acyl_coa_list_unsaturated))) / (self.Km['TAG_synthase']['acyl_coa'] \
                                                + (float(len(self.acyl_coa_list_saturated)) + float(len(self.acyl_coa_list_unsaturated))))) * (float(len(self.DAG_list)) / \
                                                (self.Km['TAG_synthase']['DAG'] + float(len(self.DAG_list))))), \

                                'TAG_lipase': 1- (float(len(self.lipid_droplets)) / (self.Km['TAG_lipase']['lipid_droplets'] + float(len(self.lipid_droplets)))), \

                                'DAG_kinase': 1- (float(len(self.DAG_list)) / (self.Km['DAG_kinase']['DAG'] + float(len(self.DAG_list)))), \

                                'PS_synthase': 1- ((float(len(self.CDP_DG_list)) / (self.Km['PS_synthase']['CDP_DG'] + float(len(self.CDP_DG_list)))) \
                                                * (self.precursors_dict['serine'] / (self.Km['PS_synthase']['serine'] + self.precursors_dict['serine']))),\

                                'PI_synthase': 1- ((float(len(self.CDP_DG_list)) / (self.Km['PI_synthase']['CDP_DG'] + float(len(self.CDP_DG_list)))) \
                                                * (self.precursors_dict['inositol'] / (self.Km['PI_synthase']['inositol'] + self.precursors_dict['inositol']))), \

                                'PE_synthase': 1- (float(len(self.PS_list)) / (self.Km['PE_synthase']['PS'] + float(len(self.PS_list)))), \

                                'PC_synthase': 1- ((float(len(self.PE_list)) / (self.Km['PC_synthase']['PE'] + float(len(self.PE_list)))) \
                                                * (self.precursors_dict['SAM'] / (self.Km['PC_synthase']['SAM'] + self.precursors_dict['SAM']))), \

                                'CL_synthase': 1- ((float(len(self.CDP_DG_list)) / (self.Km['CL_synthase']['CDP_DG'] + float(len(self.CDP_DG_list)))) \
                                                * (self.precursors_dict['glycerol_3_p_mito'] / (self.Km['CL_synthase']['glycerol_3_p_mito'] + self.precursors_dict['glycerol_3_p_mito']))),\

                                'Ergosterol_synthase': 1- (self.precursors_dict['acetyl_coa'] / (self.Km['Ergosterol_synthase']['acetyl_coa'] + self.precursors_dict['acetyl_coa'])), \

                                'Sterylester_synthase': 1- ((float(len(self.Ergosterol_list)) / (self.Km['Sterylester_synthase']['ergosterol'] + float(len(self.Ergosterol_list)))) \
                                                        * ((float(len(self.acyl_coa_list_saturated)) + float(len(self.acyl_coa_list_unsaturated))) / (self.Km['Sterylester_synthase']['acyl_coa'] \
                                                        + (float(len(self.acyl_coa_list_saturated)) + float(len(self.acyl_coa_list_unsaturated)))))), \

                                'Sphingolipid_synthase': 1- ((float(len(self.PI_list)) / (self.Km['Sphingolipid_synthase']['PI'] + float(len(self.PI_list))))\
                                                        * (self.precursors_dict['ceramide'] / (self.Km['Sphingolipid_synthase']['ceramide'] + self.precursors_dict['ceramide'])))}


            self.probabilities = self.thresholds 						#manual_threshold

            if self.phase == 'G1':
                self.probabilities['DAG_synthase'] = 1- self.probability_G1['DAG_synthase']
                self.probabilities['TAG_synthase'] = 1- self.probability_G1['TAG_synthase']
                self.probabilities['TAG_lipase'] = 1- self.probability_G1['TAG_lipase']
                self.probabilities['DAG_kinase'] = 1- self.probability_G1['DAG_kinase']
            else:
                self.probabilities['DAG_synthase'] = 1- self.probability_S_M['DAG_synthase']
                self.probabilities['TAG_synthase'] = 1- self.probability_S_M['TAG_synthase']
                self.probabilities['TAG_lipase'] = 1- self.probability_S_M['TAG_lipase']
                self.probabilities['DAG_kinase'] = 1- self.probability_S_M['DAG_kinase']
                self.probabilities['Sterylester_synthase'] = 1- self.probability_S_M['Sterylester_synthase']


            self.function_list = [self.glycerol_3_p_synthesis,
                                self.inositol_synthesis,
                                self.ceramide_synthesis,
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
                                self.DAG_kinase,
                                self.Ergosterol_synthase,
                                self.Sterylester_synthase,
                                self.Sphingolipid_synthase,
                                self.transport]
            #all reactions that take place during one second in a random order
            for i in self.function_list:
                func = random.choice(self.function_list)
                func()
                self.function_list.remove(func)
            self.numbers()				#calculating the produced lipids after each time step
        self.membrane_compositions()	#calculation of the membrane compositions at the end of the run
        self.saturation_counter()		#calculating the percentages of each fatty acid taht was used


        return self.saturation_composition_total, self.number_membranes_list, self.compartment_relatives_dict

        '''
        #printing the numbers of free lipids
        print 'CL: ' + str(self.number_CL[-1]), 'PS: ' + str(self.number_PS[-1]), 'PI: ' + str(self.number_PI[-1]), 'PE: ' + str(self.number_PE[-1]), \
                'PC: ' + str(self.number_PC[-1]), 'PA: ' + str(self.number_pa[-1]), 'TAG: ' + str(self.number_tag[-1]), 'CDP-DG: ' + str(self.number_cdp_dg[-1]),\
                'ES: ' + str(self.number_Ergosterol[-1]), 'SE: ' + str(self.number_Sterylester[-1]), 'DAG: ' + str(self.number_DAG[-1]), 'SL: ' + str(self.number_Sphingolipid[-1])
        #all free lipids added
        print self.number_CL[-1] + self.number_PS[-1] + self.number_PI[-1] + self.number_PE[-1] + self.number_PC[-1] +\
                self.number_pa[-1] + self.number_tag[-1] + self.number_cdp_dg[-1] + self.number_Ergosterol[-1] + self.number_Sterylester[-1] + self.number_Sphingolipid[-1]
        #all membrane lipids that are existing in the membranes together
        print len(self.plasma_membrane) + len(self.secretory_vesicles) + len(self.vacuoles) + len(self.nucleus)+\
                len(self.peroxisomes) + len(self.light_microsomes) + len(self.inner_mit_membrane) + \
                len(self.outer_mit_membrane) + len(self.lipid_droplets)

        print len(self.plasma_membrane), len(self.secretory_vesicles), len(self.vacuoles), len(self.nucleus), \
                len(self.peroxisomes), len(self.light_microsomes), len(self.inner_mit_membrane), \
                len(self.outer_mit_membrane), len(self.lipid_droplets)
        '''

    def plot_precursors(self):
        '''
        Plotting the precursor molecules from the precursors_dict.
        '''
        fig = mat.figure()
        ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
        i = 0
        while i < 6:#in range(len(self.precursor_keys)):
            ax.plot(self.t, self.precursors[self.precursor_keys[i]], label = self.precursor_keys[i])
            i += 1
        ax.plot(self.t, self.precursors['inositol'], label = 'inositol')
        ax.plot(self.t, self.precursors['serine'], label = 'serine')
        ax.plot(self.t, self.precursors['CTP'], label = 'CTP')
        ax.legend(bbox_to_anchor = (1.05, 1), loc = 2, borderaxespad = 0.)
        mat.show()


    def plot_lipids(self):
        '''
        Plotting the produced lipids before they are transported into the membranes
        '''
        fig = mat.figure()
        ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
        ax.plot(self.t, self.number_tag, label = 'tag')
        ax.plot(self.t, self.number_PS, label = 'ps')
        ax.plot(self.t, self.number_PI, label = 'pi')
        ax.plot(self.t, self.number_PE, label = 'pe')
        ax.plot(self.t, self.number_PC, label = 'pc')
        ax.plot(self.t, self.number_CL, label = 'cl')
        ax.plot(self.t, self.number_Ergosterol, label = 'es')
        ax.plot(self.t, self.number_Sterylester, label = 'se')
        ax.plot(self.t, self.number_Sphingolipid, label = 'sl')
        ax.legend(bbox_to_anchor = (1.05, 1), loc = 2, borderaxespad = 0.)
        mat.show()


    def plot_precursor_lipids(self):
        '''
        Plotting some free precursor molecules.
        '''
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
        Plotting the number of lipids in the membranes of different compartments.
        '''
        fig = mat.figure()
        ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
        ax.plot(self.t, self.number_plasma_membrane, label = 'plasma membrane')
        ax.plot(self.t, self.number_secretory_vesicles, label = 'secretory vesicles')
        ax.plot(self.t, self.number_vacuoles, label = 'vacuoles')
        ax.plot(self.t, self.number_nucleus, label = 'nucleus')
        ax.plot(self.t, self.number_peroxisomes, label = 'peroxisomes')
        ax.plot(self.t, self.number_light_microsomes, label = 'light_microsomes')
        ax.plot(self.t, self.number_inner_mit_membrane, label = 'inner_mit_membrane')
        ax.plot(self.t, self.number_outer_mit_membrane, label = 'outer_mit_membrane')
        ax.plot(self.t, self.number_lipid_droplets, label = 'lipid_droplets')
        ax.legend(bbox_to_anchor = (1.05, 1), loc = 2, borderaxespad = 0.)
        mat.show()


    def start(self):
        '''
        Function that produces the starting lipids. Membrane compositions are given in self.membrane_compositions_start_relatives.
        Number of starting lipids for each membrane are given in self.start_lipids.
        '''

        self.plasma_membrane_comp_start = [0.17320, 0.09124, 0.08660, 0.10464, 0.00103, 0.02010, 0.48454, 0.0, 0.0, 0.03557]
        self.secretory_vesicles_comp_start = [0.08205, 0.11745, 0.20824, 0.13573, 0.01239, 0.01525, 0.42900, 0.0, 0.0, 0.05029]
        self.vacuoles_comp_start = [0.04817, 0.16604, 0.40517, 0.17537, 0.02442, 0.02866, 0.15200, 0.0, 0.0, 0.06525]
        self.nucleus_comp_start = [0.04038, 0.09650, 0.27645, 0.16848, 0.01049, 0.01781, 0.390, 0.0, 0.0, 0.02622]
        self.peroxisomes_comp_start = [0.03235, 0.11360, 0.34656, 0.16465, 0.05033, 0.01150, 0.281, 0.0, 0.0, 0.0]
        self.light_microsomes_comp_start = [0.05304, 0.06019, 0.40796, 0.26583, 0.00381, 0.00222, 0.206, 0.0, 0.0, 0.00397]
        self.inner_mit_membrane_comp_start = [0.02880, 0.12273, 0.29107, 0.18192, 0.12204, 0.01137, 0.242, 0.0, 0.0, 0.0]
        self.outer_mit_membrane_comp_start = [0.01189, 0.10108, 0.45190, 0.32307, 0.05847, 0.04360, 0.009, 0.0, 0.0, 0.0]
        self.lipid_droplets_comp_start = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.0]

        self.membrane_compositions_start = [self.plasma_membrane_comp_start, self.secretory_vesicles_comp_start, self.vacuoles_comp_start, self.nucleus_comp_start,\
                                    self.peroxisomes_comp_start, self.light_microsomes_comp_start, self.inner_mit_membrane_comp_start,\
                                    self.outer_mit_membrane_comp_start, self.lipid_droplets_comp_start]

        self.membrane_compositions_start_relatives = []

        #numbers of the comp_start lists should yield 1 when summed
        for membrane_comp_start in self.membrane_compositions_start:
            membrane_comp_start_relative = [z / sum(membrane_comp_start) for z in membrane_comp_start]
            self.membrane_compositions_start_relatives.append(membrane_comp_start_relative)

        self.compositions_start = dict(zip(self.compartment, self.membrane_compositions_start_relatives))

        x = 0
        self.start_lipids = [32950, 500, 2500, 6000, 500, 500, 5000, 2500, 1000]		#number of lipids that are produced in the start function for every membrane
        self.membrane_start = dict(zip(self.compartment, self.start_lipids))
        for membrane in self.compartment_lists:
            for i in range(self.membrane_start[self.compartment[x]]):		#producing the lipids for a membrane, probability for a certain lipid from the composition in Zinser
                self.head_groups_start = ['serine', 'inositol', 'choline', 'ethanolamine', 'neutral', 'p', 'sterol', 'sterylester', None, 'ceramide']
                weights_start = self.compositions_start[self.compartment[x]]
                head = choice(self.head_groups_start, p = weights_start)
                if head == 'sterol':
                    new_lipid = sterol(head, self.compartment[x], self.compartment_weights)
                elif head == 'sterylester':
                    new_lipid = sterylester(head, choice(self.chainlength_unsaturated.values(), p = [0.67, 0.33]), self.compartment[x], self.compartment_weights)
                elif head == 'ceramide':
                    new_lipid = sphingolipid(head, self.compartment[x], self.compartment_weights)
                elif head == 'neutral':
                    new_lipid = CL(head, choice(self.chainlength_unsaturated.values(), p = self.unsaturated_weights), choice(self.chainlength_saturated_unsaturated, p = self.saturation_weights_total),\
                                    choice(self.chainlength_unsaturated.values(), p = self.unsaturated_weights), choice(self.chainlength_saturated_unsaturated, p = self.saturation_weights_total), self.compartment[x], self.compartment_weights)
                elif head == 'serine' or head == 'inositol' or head == 'choline' or head == 'ethanolamine' or head == 'p':
                    new_lipid = lipids(head, choice(self.chainlength_unsaturated.values(), p = self.unsaturated_weights), choice(self.chainlength_saturated_unsaturated, p = self.saturation_weights_total), self.compartment[x], self.compartment_weights)
                else:
                    new_lipid = lipids(head, choice(self.chainlength_unsaturated.values(), p = self.unsaturated_weights), choice(self.chainlength_saturated_unsaturated, p = self.saturation_weights_total), self.compartment[x], self.compartment_weights)
                    new_lipid.__class__ = TAG
                    new_lipid.sn3 = choice(self.chainlength_saturated_unsaturated, p = self.saturation_weights_total)
                membrane.append(new_lipid)
            x += 1


        self.lipid_lists_start = [self.PS_list, self.PI_list, self.PC_list, self.PE_list, self.CL_list, self.lyso_pa_list, self.PA_list, self.CDP_DG_list, self.Ergosterol_list, self.Sterylester_list, self.DAG_list, self.TAG_list, self.Sphingolipid_list]
        z = 0
        for lipid_list in self.lipid_lists_start:
            for i in range(20):
                self.head_groups_start_lipids = ['serine', 'inositol', 'choline', 'ethanolamine', 'neutral', 'lyso', 'p', 'cdp', 'sterol', 'sterylester', 'dag', None, 'ceramide']
                head = self.head_groups_start_lipids[z]
                if head == 'sterol':
                    new_lipid = sterol(head, None, self.compartment_weights)
                elif head == 'sterylester':
                    new_lipid = sterylester(head, choice(self.chainlength_unsaturated.values(), p = [0.67, 0.33]), None, self.compartment_weights)
                elif head == 'neutral':
                    new_lipid = CL(head, choice(self.chainlength_unsaturated.values(), p = self.unsaturated_weights), choice(self.chainlength_saturated_unsaturated, p = self.saturation_weights_total),\
                                            choice(self.chainlength_unsaturated.values(), p = self.unsaturated_weights), choice(self.chainlength_saturated_unsaturated, p = self.saturation_weights_total), None, self.compartment_weights)
                elif head == None:
                    new_lipid = lipids(head, choice(self.chainlength_unsaturated.values(), p = self.unsaturated_weights), choice(self.chainlength_saturated_unsaturated, p = self.saturation_weights_total), None, self.compartment_weights)
                    new_lipid.__class__ = TAG
                    new_lipid.sn3 = choice(self.chainlength_saturated_unsaturated, p = self.saturation_weights_total)
                elif head == 'ceramide':
                    new_lipid = sphingolipid(head, None, self.compartment_weights)
                elif head == 'cdp':
                    new_lipid = lipids('p', choice(self.chainlength_unsaturated.values(), p = self.unsaturated_weights), choice(self.chainlength_saturated_unsaturated, p = self.saturation_weights_total), None, self.compartment_weights)
                elif head == 'lyso':
                    new_lipid = lipids('p', None, choice(self.chainlength_saturated_unsaturated, p = self.saturation_weights_total), None, self.compartment_weights)
                elif head == 'dag':
                    new_lipid = lipids(None, choice(self.chainlength_unsaturated.values(), p = self.unsaturated_weights), choice(self.chainlength_saturated_unsaturated, p = self.saturation_weights_total), None, self.compartment_weights)
                else:
                    new_lipid = lipids(head, choice(self.chainlength_unsaturated.values(), p = self.unsaturated_weights), choice(self.chainlength_saturated_unsaturated, p = self.saturation_weights_total), None, self.compartment_weights)
                lipid_list.append(new_lipid)
            z += 1

        choice_list_acyl_start = [0, 1]
        choice_weights_acyl_start = [0.13, 0.87]
        choice_C_acyl_start = [16, 18]
        for i in range(60):
            new_acyl = fatty_acids(choice(choice_C_acyl_start), choice(choice_list_acyl_start, p = choice_weights_acyl_start))
            if new_acyl.saturation == 0:
                self.acyl_coa_list_saturated.append(new_acyl)
            else:
                self.acyl_coa_list_unsaturated.append(new_acyl)

    def Km_calculation(self):
        #Calculation of Km, 2 substrates: Km for lipids = 5, for acyl-coa = 30
        self.Km = {'glycerol_3_p_synthesis': {'DHAP': 0.}, \
                        'inositol_synthesis': {'glucose_6_p': 0.}, \
                        'ceramide_synthesis': {'serine': 0.},\
                        'acetyl_coa_synthase': {'pyruvate': 0.}, \
                        'acyl_synthase': {'acetyl_coa': 0.},\
                        'acyl_synthase_C16': 0.625,\
                        'acyl_synthase_C18': 0.002,\
                        'lyso_PA_synthase': {'acyl_coa': 30., 'DHAP': 0.}, \
                        'PA_synthase': {'lyso-PA': 5., 'acyl_coa': 30.},\
                        'CDP_DG_synthase': {'PA': 5., 'CTP': 0.},\
                        'DAG_synthase': {'PA': 5.}, \
                        'TAG_synthase': {'DAG': 5., 'acyl_coa': 30.}, \
                        'TAG_lipase': {'lipid_droplets': 0.}, \
                        'DAG_kinase': {'DAG': 5.}, \
                        'PS_synthase': {'CDP_DG': 5., 'serine': 0.},\
                        'PI_synthase': {'CDP_DG': 5., 'inositol': 0.}, \
                        'PE_synthase': {'PS': 5.}, \
                        'PC_synthase': {'PE': 5., 'SAM': 0.}, \
                        'CL_synthase': {'CDP_DG': 5., 'glycerol_3_p_mito': 0.},\
                        'Ergosterol_synthase': {'acetyl_coa': 0.}, \
                        'Sterylester_synthase': {'ergosterol': 5., 'acyl_coa': 30.}, \
                        'Sphingolipid_synthase': {'PI': 5}}

        self.Km = {'glycerol_3_p_synthesis': {'DHAP': self.precursors_dict['DHAP'] / self.probability['glycerol_3_p_synthesis'] - self.precursors_dict['DHAP']}, \
                        'inositol_synthesis': {'glucose_6_p': self.precursors_dict['glucose_6_p'] / self.probability['inositol_synthesis'] - self.precursors_dict['glucose_6_p']}, \
                        'ceramide_synthesis': {'serine': self.precursors_dict['serine'] / self.probability['ceramide_synthesis'] - self.precursors_dict['serine']},\
                        'acetyl_coa_synthase': {'pyruvate': self.precursors_dict['pyruvate'] / self.probability['acetyl_coa_synthase'] - self.precursors_dict['pyruvate']}, \
                        'acyl_synthase': {'acetyl_coa': self.precursors_dict['acetyl_coa'] / self.probability['acyl_synthase'] - self.precursors_dict['acetyl_coa']},\
                        'acyl_synthase_C16': 0.625,\
                        'acyl_synthase_C18': 0.002,\
                        'lyso_PA_synthase': {'acyl_coa': 30., 'DHAP': (self.precursors_dict['DHAP'] / self.probability['lyso_PA_synthase']) * ((float(len(self.acyl_coa_list_saturated)) + float(len(self.acyl_coa_list_unsaturated)))\
                                             / (self.Km['lyso_PA_synthase']['acyl_coa'] + (float(len(self.acyl_coa_list_saturated)) + float(len(self.acyl_coa_list_unsaturated))))) - self.precursors_dict['DHAP']}, \
                        'PA_synthase': {'lyso-PA': 5., 'acyl_coa': 30.},\
                        'CDP_DG_synthase': {'PA': 5., 'CTP': (self.precursors_dict['CTP'] / self.probability['CDP_DG_synthase']) * (float(len(self.PA_list)) / (self.Km['CDP_DG_synthase']['PA'] + float(len(self.PA_list)))) - self.precursors_dict['CTP']},\
                        'DAG_synthase': {'PA': 5.}, \
                        'TAG_synthase': {'DAG': 5., 'acyl_coa': 30.}, \
                        'TAG_lipase': {'lipid_droplets': float(len(self.lipid_droplets)) / self.probability['TAG_lipase'] - float(len(self.lipid_droplets))}, \
                        'DAG_kinase': {'DAG': 5.}, \
                        'PS_synthase': {'CDP_DG': 5., 'serine': (self.precursors_dict['serine'] / self.probability['PS_synthase']) * (float(len(self.CDP_DG_list)) / (self.Km['PS_synthase']['CDP_DG'] + float(len(self.CDP_DG_list)))) - self.precursors_dict['serine']},\
                        'PI_synthase': {'CDP_DG': 5., 'inositol': (self.precursors_dict['inositol'] / self.probability['PI_synthase']) * (float(len(self.CDP_DG_list)) / (self.Km['PI_synthase']['CDP_DG'] + float(len(self.CDP_DG_list)))) - self.precursors_dict['inositol']}, \
                        'PE_synthase': {'PS': 5.}, \
                        'PC_synthase': {'PE': 5., 'SAM': (self.precursors_dict['SAM'] / self.probability['PC_synthase']) * (float(len(self.PE_list)) / (self.Km['PC_synthase']['PE'] + float(len(self.PE_list)))) - self.precursors_dict['SAM']}, \
                        'CL_synthase': {'CDP_DG': 5., 'glycerol_3_p_mito': (self.precursors_dict['glycerol_3_p_mito'] / self.probability['CL_synthase']) * (float(len(self.CDP_DG_list)) / (self.Km['CL_synthase']['CDP_DG'] + float(len(self.CDP_DG_list)))) - self.precursors_dict['glycerol_3_p_mito']},\
                        'Ergosterol_synthase': {'acetyl_coa': self.precursors_dict['acetyl_coa'] / self.probability['Ergosterol_synthase'] - self.precursors_dict['acetyl_coa']}, \
                        'Sterylester_synthase': {'ergosterol': 5., 'acyl_coa': 30.}, \
                        'Sphingolipid_synthase': {'PI': 5, 'ceramide': (self.precursors_dict['ceramide'] / self.probability['Sphingolipid_synthase']) * (float(len(self.PI_list)) / (self.Km['Sphingolipid_synthase']['PI'] + float(len(self.PI_list)))) - self.precursors_dict['ceramide']}}


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


    def glycerol_3_p_synthesis(self):
        '''
        Synthesis of glycerol-3-p out of DHAP.
        '''
        for i in range(self.rates['glycerol_3_p_synthesis']):
            x = random.random()
            if x >= self.probabilities['glycerol_3_p_synthesis']:
                if self.precursors_dict['DHAP'] > 1:
                    self.precursors_dict['glycerol-3-p'] += 1
                    self.precursors_dict['DHAP'] -= 1
                    self.precursors_dict['NADH'] += 1
                    self.precursors_dict['NAD'] -= 1


    def inositol_synthesis(self):
        '''
        Synthesis of inositol from glucose-6-p by the Myo-inositol-1-p synthase and the Myo-inositol 1
        phosphatase.
        '''
        for i in range(self.rates['inositol_synthesis']):
            x = random.random()
            if x >= self.probabilities['inositol_synthesis']:
                if self.precursors_dict['glucose_6_p'] > 1:
                    self.precursors_dict['inositol'] += 1
                    self.precursors_dict['glucose_6_p'] -= 1
                    self.precursors_dict['H2O'] -= 1
                    self.precursors_dict['Pi'] += 1


    def ceramide_synthesis(self):
        '''
        Synthesis of ceramide out of serine and a C16:0 fatty acid
        '''
        for i in range(self.rates['ceramide_synthesis']):
            x = random.random()
            if x >= self.probabilities['ceramide_synthesis']:
                if len(self.acyl_coa_list_C26) > 1 and self.precursors_dict['serine'] > 1 and len(self.acyl_coa_list_saturated) > 1 and any(fa.C == 16 for fa in self.acyl_coa_list_saturated):
                    self.precursors_dict['ceramide'] += 1
                    self.precursors_dict['serine'] -= 1
                    self.precursors_dict['CO2'] += 1
                    self.precursors_dict['NADPH'] += 1
                    self.precursors_dict['NADP'] -= 1
                    j = 0
                    if self.acyl_coa_list_saturated[j].C == 18 and i <= len(self.acyl_coa_list_saturated):
                        j += 1
                    elif self.acyl_coa_list_saturated[j].C == 16:
                        del self.acyl_coa_list_saturated[j]
                        del self.acyl_coa_list_C26[0]


    def acetyl_coa_synthase(self):
        '''
        Synthesis of Acetyl-CoA: pyruvate dehydrogenase drives the reaction pyruvate to Acetyl-CoA, CO2 is released
        '''
        for i in range(self.rates['acetyl_coa_synthase']):
            x = random.random()
            if x >= self.probabilities['acetyl_coa_synthase']:
                if self.precursors_dict['pyruvate'] > 1:			# transformation from pyruvate to acetyl_coa
                    self.precursors_dict['acetyl_coa'] += 1
                    self.precursors_dict['pyruvate'] -= 1
                    self.precursors_dict['NADH'] += 1
                    self.precursors_dict['NAD'] -= 1
                    self.precursors_dict['CO2'] += 1


    def acyl_synthase(self):
        '''
        Simplified synthesis of acyl_coa: several Acetyl-CoA are building a fatty acid (C16:0, C16:1, C18:0 or C18:1)
        The intermediate Malonyl-CoA is leaved out.
        '''
        choice_list = [0, 1]
        choice_weights = [0.12, 0.88]
        for i in range(self.rates['acyl_synthase']):
            x = random.random()						#5 reactions in 1 timestep but only with a probability of 90%
            if self.precursors_dict['acetyl_coa'] > 2:		#control if at least 2 Acetyl-CoA are available
                if len(self.acyl_coa_list) == 0:		#starting the first reaction
                    new_acyl = fatty_acids(2, choice(choice_list, p = choice_weights))
                    self.acyl_coa_list.append(new_acyl)
                    self.acyl_coa_list[-1].C += 2
                    self.precursors_dict['acetyl_coa'] -= 2

                elif self.acyl_coa_list[-1].C == 16 and x >= self.probabilities['acyl_synthase_C16']:	#stop the reaction cycle and starting a new one
                    new_acyl = fatty_acids(2, choice(choice_list, p = choice_weights))
                    self.acyl_coa_list.append(new_acyl)
                    self.acyl_coa_list[-1].C += 2
                    self.precursors_dict['acetyl_coa'] -= 2
                    self.precursors_dict['NADPH'] -= 14
                    self.precursors_dict['NADP'] += 14
                    self.precursors_dict['H2O'] += 7
                    #CO2 production is not mentioned here as onyl acetyl-CoA is used and not malonyl-CoA, so we need all C-atoms we give in the reaction

                elif self.acyl_coa_list[-1].C == 18 and x >= self.probabilities['acyl_synthase_C18']:	#stop the reaction cycle and starting a new one
                    new_acyl = fatty_acids(2, choice(choice_list, p = choice_weights))
                    self.acyl_coa_list.append(new_acyl)
                    self.acyl_coa_list[-1].C += 2
                    self.precursors_dict['acetyl_coa'] -= 2
                    self.precursors_dict['NADPH'] -= 16
                    self.precursors_dict['NADP'] += 16
                    self.precursors_dict['H2O'] += 8

                elif self.acyl_coa_list[-1].C == 26:
                    self.acyl_coa_list[-1].saturation = 0
                    new_acyl = fatty_acids(2, choice(choice_list, p = choice_weights))
                    self.acyl_coa_list.append(new_acyl)
                    self.acyl_coa_list[-1].C += 2
                    self.precursors_dict['acetyl_coa'] -= 2
                    self.precursors_dict['NADPH'] -= 24
                    self.precursors_dict['NADP'] += 24
                    self.precursors_dict['H2O'] += 12

                else:									#adding an Acetyl_CoA to the growing ffa
                    self.acyl_coa_list[-1].C += 2
                    self.precursors_dict['acetyl_coa'] -= 1

        if len(self.acyl_coa_list) > 1:
            for j in range(len(self.acyl_coa_list)-1):
                if self.acyl_coa_list[j].C == 26:
                    self.acyl_coa_list_C26.append(self.acyl_coa_list[j])
                elif self.acyl_coa_list[j].saturation == 0:
                    self.acyl_coa_list_saturated.append(self.acyl_coa_list[j])
                elif self.acyl_coa_list[j].saturation == 1:
                    self.acyl_coa_list_unsaturated.append(self.acyl_coa_list[j])
                    self.precursors_dict['O2'] -= 1
                    self.precursors_dict['H2O'] += 2
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
        choice_list = [0, 1]

        weights_pa = [self.precursors_dict['DHAP'] / (self.precursors_dict['DHAP'] + self.precursors_dict['glycerol-3-p']),\
                    self.precursors_dict['glycerol-3-p'] / (self.precursors_dict['DHAP'] + self.precursors_dict['glycerol-3-p'])]
        x = random.random()
        if x >= self.probabilities['lyso_PA_synthase'] and len(self.acyl_coa_list_saturated) > 1 and len(self.acyl_coa_list_unsaturated) > 1 and (self.precursors_dict['DHAP'] > 1 and self.precursors_dict['glycerol-3-p'] > 1): 	#at least 1 ffa has to be unsaturated
            if choice(choice_list, p = self.weights_fa) == 0:
                sn1_chain = random.randint(0, (len(self.acyl_coa_list_saturated)-1))
                chainlength_sn1 = self.acyl_coa_list_saturated[sn1_chain].C
                lyso_pa = lipids('p', None, self.chainlength_saturated[chainlength_sn1], None, self.compartment_weights)
                del self.acyl_coa_list_saturated[sn1_chain]
            else:
                sn1_chain = random.randint(0, (len(self.acyl_coa_list_unsaturated)-1))
                chainlength_sn1 = self.acyl_coa_list_unsaturated[sn1_chain].C
                lyso_pa = lipids('p', None, self.chainlength_unsaturated[chainlength_sn1], None, self.compartment_weights)
                del self.acyl_coa_list_unsaturated[sn1_chain]
            self.lyso_pa_list.append(lyso_pa)
            i = choice(choice_list, p = weights_pa)
            if i == 0:
                self.precursors_dict['DHAP'] -= 1
                self.precursors_dict['NADPH'] += 1
                self.precursors_dict['NADP'] -= 1
            else:
                self.precursors_dict['glycerol-3-p'] -= 1


    def PA_synthase(self):
        '''
        Synthesis of PA by adding the second fatty acid to lyso_PA (sn2: saturated or unsaturated) --> 1-acyl-sn-glycerol-3-phosphate acyltransferase
        '''
        x = random.random()
        if x >= self.probabilities['PA_synthase'] and len(self.acyl_coa_list_unsaturated) > 1 and len(self.lyso_pa_list) > 1:
            z = random.randint(0, (len(self.lyso_pa_list)-1))
            sn2_chain = random.randint(0, (len(self.acyl_coa_list_unsaturated)-1))
            chainlength_sn2 = self.acyl_coa_list_unsaturated[sn2_chain].C
            self.lyso_pa_list[z].sn2 = self.chainlength_unsaturated[chainlength_sn2]
            self.PA_list.append(self.lyso_pa_list[z])
            del self.acyl_coa_list_unsaturated[sn2_chain]		# deletion of the consumed ffa
            del self.lyso_pa_list[z]


    def CDP_DG_synthase(self):
        '''
        PA is processed to CDP-DG (CDP-diacylglycerol synthase), that further reacts to the phospholipids
        '''
        for i in range(self.rates['CDP_DG_synthase']):
            x = random.random()
            if x >= self.probabilities['CDP_DG_synthase'] and self.precursors_dict['CTP'] > 1 and len(self.PA_list) > 1:
                z = random.randint(0, len(self.PA_list)-1)
                self.PA_list[z].head = 'cdp'
                self.CDP_DG_list.append(self.PA_list[z])		#CDP-DG production from PA
                del self.PA_list[z]
                self.precursors_dict['CTP'] -= 1
                self.precursors_dict['Pi'] += 2


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
        x = random.random()
        if x >= self.probabilities['DAG_synthase'] and len(self.PA_list) > 1:
            z = random.randint(0, len(self.PA_list)-1)
            self.PA_list[z].head = None
            self.DAG_list.append(self.PA_list[z])
            self.precursors_dict['H2O'] -= 1
            self.precursors_dict['Pi'] += 1
            del self.PA_list[z]


    def TAG_synthase(self):
        '''
        DAG is processed to TAG by adding a third acyl-chain at position sn3.
        '''
        x = random.random()
        if x >= self.probabilities['TAG_synthase'] and len(self.DAG_list) > 1 and len(self.acyl_coa_list_saturated) > 1 and len(self.acyl_coa_list_unsaturated) > 1:
            z = random.randint(0, len(self.DAG_list)-1)
            self.TAG_list.append(self.DAG_list[z])
            self.TAG_list[-1].__class__ = TAG
            if x <= 0.575:
                sn3 = random.randint(0, len(self.acyl_coa_list_saturated)-1)
                chainlength_sn3 = self.acyl_coa_list_saturated[sn3].C
                self.TAG_list[-1].sn3 = self.chainlength_saturated[chainlength_sn3]
                del self.acyl_coa_list_saturated[sn3]
            else:
                sn3 = random.randint(0, len(self.acyl_coa_list_unsaturated)-1)
                chainlength_sn3 = self.acyl_coa_list_unsaturated[sn3].C
                self.TAG_list[-1].sn3 = self.chainlength_unsaturated[chainlength_sn3]
                del self.acyl_coa_list_unsaturated[sn3]
            del self.DAG_list[z]


    def TAG_lipase(self):
        '''
        Cdk1/Cdc28-dependent activation of the major triacylglycerol lipase
        '''
        if len(self.lipid_droplets) > self.rates['TAG_lipase']:
            for i in range(self.rates['TAG_lipase']):
                x = random.random()
                if x >= self.probabilities['TAG_lipase']:
                    z = random.randint(0, len(self.lipid_droplets)-1)
                    if self.lipid_droplets[z].head == None:
                        if ':0' in self.lipid_droplets[z].sn3:
                            for key, value in self.chainlength_unsaturated.items():
                                if value == self.lipid_droplets[z].sn3:
                                    self.acyl_coa_list_saturated.append(fatty_acids(key, 0))
                        elif ':1' in self.lipid_droplets[z].sn3:
                            for key, value in self.chainlength_saturated.items():
                                if value == self.lipid_droplets[z].sn3:
                                    self.acyl_coa_list_unsaturated.append(fatty_acids(key, 1))
                        self.DAG_list.append(self.lipid_droplets[z])
                        self.DAG_list[-1].__class__ = lipids
                        delattr(self.DAG_list[-1], 'sn3')
                        self.precursors_dict['H2O'] -= 1
                    elif self.lipid_droplets[z].head == 'sterylester':
                        self.Ergosterol_list.append(sterol('sterol', None, self.compartment_weights))
                        self.precursors_dict['H2O'] -= 1
                        if ':0' in self.lipid_droplets[z].FA:
                            for key, value in self.chainlength_unsaturated.items():
                                if value == self.lipid_droplets[z].FA:
                                    self.acyl_coa_list_saturated.append(fatty_acids(key, 0))
                        elif ':1' in self.lipid_droplets[z].FA:
                            for key, value in self.chainlength_saturated.items():
                                if value == self.lipid_droplets[z].FA:
                                    self.acyl_coa_list_unsaturated.append(fatty_acids(key, 1))
                    del self.lipid_droplets[z]


    def DAG_kinase(self):
        if len(self.DAG_list) > self.rates['DAG_kinase']:
            for i in range(self.rates['DAG_kinase']):
                x = random.random()
                if x >= self.probabilities['DAG_kinase']:
                    z = random.randint(0, len(self.DAG_list)-1)
                    self.PA_list.append(self.DAG_list[z])
                    self.PA_list[-1].head = 'p'
                    self.PA_list[-1].comp = None
                    del self.DAG_list[z]


    def PS_synthase(self):
        '''
        CDP-DG is processed to PS (PS synthase).
        '''
        for i in range(self.rates['PS_synthase']):
            x = random.random()
            if x >= self.probabilities['PS_synthase'] and len(self.CDP_DG_list) > 1 and self.precursors_dict['serine'] > 1:
                z = random.randint(0, len(self.CDP_DG_list)-1)
                self.CDP_DG_list[z].head = 'serine'				#PS synthesis from CDP-DG
                self.PS_list.append(self.CDP_DG_list[z])
                del self.CDP_DG_list[z]
                self.precursors_dict['serine'] -= 1
                self.precursors_dict['CMP'] += 1


    def PI_synthase(self):
        '''
        CDP-DG is processed to PI (PI synthase)
        '''
        for i in range(self.rates['PI_synthase']):
            x = random.random()
            if x >= self.probabilities['PI_synthase'] and len(self.CDP_DG_list) > 1 and self.precursors_dict['inositol'] > 1:
                z = random.randint(0, len(self.CDP_DG_list)-1)
                self.CDP_DG_list[z].head = 'inositol'
                self.PI_list.append(self.CDP_DG_list[z])
                del self.CDP_DG_list[z]
                self.precursors_dict['inositol'] -= 1
                self.precursors_dict['CMP'] += 1


    def PE_synthase(self):
        '''
        PE is derived from PS by releasing 1 CO2 --> PS decarboxylase.
        '''
        for i in range(self.rates['PE_synthase']):
            x = random.random()
            if x >= self.probabilities['PE_synthase'] and len(self.PS_list) >= 10:
                z = random.randint(0, len(self.PS_list)-1)
                self.PS_list[z].head = 'ethanolamine'
                self.PE_list.append(self.PS_list[z])
                self.precursors_dict['CO2'] += 1
                del self.PS_list[z]


    def PC_synthase(self):
        '''
        PC is derived from PE. As enzymes serve 3 methyltransferases which need SAM and produce SAH as a side product.
        '''
        for i in range(self.rates['PC_synthase']):
            x = random.random()
            if x >= self.probabilities['PC_synthase'] and len(self.PE_list) >= 5 and self.precursors_dict['SAM'] >= 4:
                z = random.randint(0, len(self.PE_list)-1)
                self.PE_list[z].head = 'choline'
                self.PC_list.append(self.PE_list[z])
                del self.PE_list[z]
                self.precursors_dict['SAM'] -= 3
                self.precursors_dict['SAH'] += 3


    def CL_synthase(self):
        '''
        Synthesis of cardiolipin, for which 2 CDP-DG are needed. Different enzymes are needed.
        '''
        for i in range(self.rates['CL_synthase']):
            x = random.random()
            if x >= self.probabilities['CL_synthase'] and self.precursors_dict['glycerol_3_p_mito'] > 1 and len(self.CDP_DG_list) > 2:
                z = random.randint(0, len(self.CDP_DG_list)-2)
                self.CDP_DG_list[z].head = 'neutral'
                self.CL_list.append(self.CDP_DG_list[z])
                self.CL_list[-1].__class__ = CL
                self.CL_list[-1].sn4, self.CL_list[-1].sn3 = self.CDP_DG_list[z+1].sn2, self.CDP_DG_list[z+1].sn1
                del self.CDP_DG_list[z:z+1]
                self.precursors_dict['glycerol_3_p_mito'] -= 1
                self.precursors_dict['H2O'] -= 1
                self.precursors_dict['Pi'] += 1
                self.precursors_dict['CMP'] += 2


    def Ergosterol_synthase(self):
        '''
        Synthesis of the most existing sterol in yeast: ergosterol
        '''
        for i in range(self.rates['Ergosterol_synthase']):
            x = random.random()
            if x >= self.probabilities['Ergosterol_synthase'] and self.precursors_dict['acetyl_coa'] > 18:
                self.Ergosterol_list.append(sterol('sterol', None, self.compartment_weights))
                self.precursors_dict['acetyl_coa'] -= 18
                self.precursors_dict['ATP'] -= 3
                self.precursors_dict['ADP'] += 3
                self.precursors_dict['NADPH'] -= 11
                self.precursors_dict['NADP'] += 11
                self.precursors_dict['SAM'] -= 1
                self.precursors_dict['SAH'] += 1
                self.precursors_dict['O2'] -= 7
                self.precursors_dict['H2O'] += 9
                self.precursors_dict['CO2'] += 2


    def Sterylester_synthase(self):
        '''
        Synthesis of sterylesters that are found in lipid droplets out of ergosterol and an unsaturated fatty acid.
        '''
        for i in range(self.rates['Sterylester_synthase']):
            x = random.random()
            if x >= self.probabilities['Sterylester_synthase'] and any(fa.C == 16 for fa in self.acyl_coa_list_unsaturated) and any(fa.C == 18 for fa in self.acyl_coa_list_unsaturated) and len(self.Ergosterol_list) > 1:
                z = random.randint(0, len(self.Ergosterol_list)-1)
                j = 0
                while j < 5:
                    fa_index = random.randint(0, len(self.acyl_coa_list_unsaturated)-1)
                    if self.acyl_coa_list_unsaturated[fa_index].C == 18 and random.random() < 0.33:
                        self.Sterylester_list.append(sterylester('sterylester', 'C18:1', None, self.compartment_weights))
                        del self.Ergosterol_list[z]
                        del self.acyl_coa_list_unsaturated[fa_index]
                        break
                    elif self.acyl_coa_list_unsaturated[fa_index].C == 16:
                        self.Sterylester_list.append(sterylester('sterylester', 'C16:1', None, self.compartment_weights))
                        del self.Ergosterol_list[z]
                        del self.acyl_coa_list_unsaturated[fa_index]
                        break
                    else:
                        j += 1


    def Sphingolipid_synthase(self):
        '''
        Synthesis of the most abundant Sphingolipid mannose-(inositol-phosphate)2-ceramide
        '''
        for i in range(self.rates['Sphingolipid_synthase']):
            x = random.random()
            if x >= self.probabilities['Sphingolipid_synthase'] and len(self.PI_list) >= 2 and self.precursors_dict['ceramide'] > 1 and self.precursors_dict['GDP-mannose'] > 1:
                self.Sphingolipid_list.append(sphingolipid('ceramide', None, self.compartment_weights))
                z= random.randint(0, len(self.PI_list)-2)
                del self.PI_list[z:z+1]
                self.precursors_dict['ceramide'] -= 1
                self.precursors_dict['GDP-mannose'] -= 1


    def transport(self):
        '''
        General transport function for all produced lipids.
        '''
        for lipid in self.lipid_lists:
            if lipid == self.TAG_list or lipid == self.Sterylester_list:
                if len(lipid) > 10:
                    for j in range(len(lipid)/10):
                        z = random.randint(0, len(lipid)-1)
                        lipid[z].comp_choice()
                        if lipid[z].comp == 'lipid_droplets':
                            self.lipid_droplets.append(lipid[z])
                        del lipid[z]
            else:
                if len(lipid) > 5:
                    for j in range(len(lipid)/10):
                        z = random.randint(0, len(lipid)-1)
                        lipid[z].comp_choice()
                        if lipid[z].comp == 'plasma_membrane':
                            self.plasma_membrane.append(lipid[z])
                        elif lipid[z].comp == 'secretory_vesicles':
                            self.secretory_vesicles.append(lipid[z])
                        elif lipid[z].comp == 'vacuoles':
                            self.vacuoles.append(lipid[z])
                        elif lipid[z].comp == 'nucleus':
                            self.nucleus.append(lipid[z])
                        elif lipid[z].comp == 'peroxisomes':
                            self.peroxisomes.append(lipid[z])
                        elif lipid[z].comp == 'light_microsomes':
                            self.light_microsomes.append(lipid[z])
                        elif lipid[z].comp == 'inner_mit_membrane':
                            self.inner_mit_membrane.append(lipid[z])
                        elif lipid[z].comp == 'outer_mit_membrane':
                            self.outer_mit_membrane.append(lipid[z])
                        elif lipid[z].comp == 'lipid_droplets':
                            self.lipid_droplets.append(lipid[z])
                        del lipid[z]


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
                                        (float(sum(j.head == 'sterylester' for j in comp)) / len(comp)),\
                                        (float(sum(j.head == None for j in comp)) / len(comp)),\
                                        (float(sum(j.head == 'ceramide' for j in comp)) / len(comp))]
                for i in range(len(self.relatives_list)):
                    self.compartment_relatives_dict[self.compartment[x]][self.membrane_lipids[i]] = self.relatives_list[i]
            x += 1


    def numbers(self):
        for current_lipid_number, number_of_lipid in zip(self.number_lipids_list, self.precursor_list):
            current_lipid_number.append(len(number_of_lipid))
        #for plotting the number of lipids in a certain membrane
        for i, sp in enumerate(self.compartment_lists):
            self.number_membranes_list[i].append(len(sp))
        # for current_membrane_number, number_of_membrane in zip(self.number_membranes_list, self.compartment_lists):
        #     current_membrane_number.append(len(number_of_membrane))

        for z in range(len(self.precursor_keys)):
            self.precursors[self.precursor_keys[z]].append(self.precursors_dict[self.precursor_keys[z]])


    def saturation_counter(self):
        '''
        composition of fatty acids that should be reached: C16:0 = 10%, C16:1 = 30%, C18:0 = 10%, C18:1 = 50%
        these numbers are from Klug & Daum 2013 'Yeast lipid metabolism at a glance'
        '''
        self.c16_0_sn1 = 0
        self.c16_1_sn1 = 0
        self.c18_0_sn1 = 0
        self.c18_1_sn1 = 0
        self.wrong_fatty_acid = 0
        for c in self.compartment_lists:
            if c == self.lipid_droplets:
                continue
            else:
                for i in range(len(c)):
                    if hasattr(c[i], 'sn1'):
                        if c[i].sn1 == 'C16:0':
                            self.c16_0_sn1 += 1
                        elif c[i].sn1 == 'C16:1':
                            self.c16_1_sn1 += 1
                        elif c[i].sn1 == 'C18:0':
                            self.c18_0_sn1 += 1
                        elif c[i].sn1 == 'C18:1':
                            self.c18_1_sn1 += 1
                        else:
                            self.wrong_fatty_acid += 1
        self.saturation_composition_sn1 = {'C16:0': self.c16_0_sn1, 'C16:1': self.c16_1_sn1, 'C18:0': self.c18_0_sn1, 'C18_1': self.c18_1_sn1}

        self.c16_0_sn2 = 0
        self.c16_1_sn2 = 0
        self.c18_0_sn2 = 0
        self.c18_1_sn2 = 0
        self.wrong_fatty_acid = 0
        for c in self.compartment_lists:
            if c == self.lipid_droplets:
                continue
            else:
                for i in range(len(c)):
                    if hasattr(c[i], 'sn2'):
                        if c[i].sn2 == 'C16:0':
                            self.c16_0_sn2 += 1
                        elif c[i].sn2 == 'C16:1':
                            self.c16_1_sn2 += 1
                        elif c[i].sn2 == 'C18:0':
                            self.c18_0_sn2 += 1
                        elif c[i].sn2 == 'C18:1':
                            self.c18_1_sn2 += 1
                        else:
                            self.wrong_fatty_acid += 1

        self.saturation_composition_sn2 = {'C16:0': self.c16_0_sn2, 'C16:1': self.c16_1_sn2, 'C18:0': self.c18_0_sn2, 'C18_1': self.c18_1_sn2}
        self.total_fatty_acids = self.c16_0_sn1 + self.c16_1_sn1 + self.c18_0_sn1 + self.c18_1_sn1 + self.c16_0_sn2 + self.c16_1_sn2 + self.c18_0_sn2 + self.c18_1_sn2
        #if self.total_fatty_acids > 0:
        #	self.saturation_composition_total = {'C16:0': float(self.c16_0_sn2 + self.c16_0_sn1) / self.total_fatty_acids, 'C16:1': float(self.c16_1_sn2 + self.c16_1_sn1) / self.total_fatty_acids, \
        #										'C18:0': float(self.c18_0_sn2 + self.c18_0_sn1) / self.total_fatty_acids, 'C18:1': float(self.c18_1_sn2 + self.c18_1_sn1) / self.total_fatty_acids}

        self.sterylester_C16 = 0
        self.sterylester_C18 = 0
        for c in self.lipid_droplets:
            if c.head == 'sterylester':
                if c.FA == 'C16:1':
                    self.sterylester_C16 += 1
                elif c.FA == 'C18:1':
                    self.sterylester_C18 += 1

        if self.sterylester_C16 > 0 or self.sterylester_C18 > 0:
            self.composition_sterylester = {'C16:1: ': float(self.sterylester_C16) / (self.sterylester_C16 + self.sterylester_C18),\
                                            'C18:1: ': float(self.sterylester_C18) / (self.sterylester_C16 + self.sterylester_C18)}

        if self.total_fatty_acids > 0:
            self.saturation_composition_total = {'C16:0': float(self.c16_0_sn2 + self.c16_0_sn1) / self.total_fatty_acids,\
                                                 'C16:1': float(self.c16_1_sn2 + self.c16_1_sn1) / self.total_fatty_acids,\
                                                 'C18:0': float(self.c18_0_sn2 + self.c18_0_sn1) / self.total_fatty_acids,\
                                                 'C18:1': float(self.c18_1_sn2 + self.c18_1_sn1) / self.total_fatty_acids}

if __name__ == '__main__':
    # Test run, with runtime tracker
    import time

    st = time.time()
    m = model()
    # test run: 5 sec
    r, mem, s = m.run()
    et = time.time()
    for lili in m.number_lipids_list:
        mat.plot(m.t, lili)
    print "Runtime: " + str(et - st) + "s"
    mat.show()