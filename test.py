# -*- coding: utf-8 -*-
"""
Created/Started on Wed June 03 2015

@author: Vera
"""
import matplotlib.pyplot as mat
import numpy as np
import components


class Model:
    """
    The model.
    At the beginning there are several lists defined which will contain the produced lipids.
    The model starts with a given amount of precursors (pyruvate, dhap, ctp, serine, glucose-6-p, SAM, glycerol-3-p).
    After several reactions and the transport function in the end there are membranes of different compartments 
    with several different lipids.
    """
    def __init__(self):
        self.__init_parameters()
        self.__init_lists()
        self.__init_precursors()
        self.__init_molecules()
        self.__init_simulation()

    def __init_parameters(self):
        """
        Initialise model parameters, all values were adjusted manually to the data of Uchida et al.
        (2011, PMID:21360734) and Zinser et al. (1991, PMID:2002005).
        """
        # VMAX of reactions
        # adjusted manually
        self.rates = {'glycerol_3_p_synthesis': 8,
                      'inositol_synthesis': 5,
                      'ceramide_synthesis': 2,
                      'acetyl_coa_synthase': 650,
                      'acyl_synthase': 450,
                      'PA_synthesis': 17,
                      'CDP_DG_synthase': 20,
                      'TAG_synthesis': 30,
                      'TAG_lipase': 23,
                      'DAG_kinase': 40,
                      'PS_synthase': 18,
                      'PI_synthase': 6,
                      'PE_synthase': 12,
                      'PC_synthase': 5,
                      'CL_synthase': 2,
                      'ergosterol_synthase': 25,
                      'sterylester_synthase': 25,
                      'sphingolipid_synthase': 2}

        # probabilities of reaction to take place
        # adjusted manually
        self.initial_probabilities = {'glycerol_3_p_synthesis': 0.5,
                                      'inositol_synthesis': 0.5,
                                      'ceramide_synthesis': 0.5,
                                      'acetyl_coa_synthase': 0.8,
                                      'acyl_synthase': 0.5,
                                      'acyl_synthase_C16': 0.375,
                                      'acyl_synthase_C18': 0.998,
                                      'lyso_PA_synthase': 0.45,
                                      'PA_synthase': 0.2,
                                      'CDP_DG_synthase': 0.8,
                                      'DAG_synthase': 0.01,
                                      'TAG_synthase': 0.2,
                                      'TAG_lipase': 0.8,
                                      'DAG_kinase': 0.1,
                                      'PS_synthase': 0.5,
                                      'PI_synthase': 0.5,
                                      'PE_synthase': 0.5,
                                      'PC_synthase': 0.5,
                                      'CL_synthase': 0.05,
                                      'ergosterol_synthase': 0.6,
                                      'sterylester_synthase': 0.4,
                                      'sphingolipid_synthase': 0.2}

        # probabilities in CC phase G1
        # adjusted manually
        self.cc_probabilities = {'G1': {'DAG_synthase': 0.3,
                                        'TAG_synthase': 0.2,
                                        'TAG_lipase': 0.05,
                                        'DAG_kinase': 0.03},
                                 'S_M': {'DAG_synthase': 0.01,
                                         'TAG_synthase': 0.2,
                                         'TAG_lipase': 0.6,
                                         'DAG_kinase': 0.1,
                                         'sterylester_synthase': 0.2}}

        self.weights = {'compartments': [0.67, 0.01, 0.155, 0.101, 0.007, 0.007, 0.03, 0.015],
                        # compartment size ratio (molecules) - transport probabilities: adjusted manually
                        # PM, SecVec, Vac, Nuc, Perox, lightMic, MitoMemIn, MitoMemOut
                        'chain_saturated': {16: 'C16:0', 18: 'C18:0'},
                        # possible fatty acids and weights to reach the biological proportions
                        'chain_unsaturated': {16: 'C16:1', 18: 'C18:1'},
                        'chain_complete': ['C16:0', 'C18:0', 'C16:1', 'C18:1'],
                        'FA': [0.4, 0.6],
                        # saturated, unsaturated: adjusted manually
                        'unsaturated': [0.375, 0.625],
                        'saturation_total': [0.2, 0.2, 0.167, 0.433]
                        # sn1: C16:0, C18:0, C16:1, C18:1
                        }

    def __init_lists(self):
        self.components_state = {'acyl_coa': [], 'acyl_coa_C26': [], 'acyl_coa_saturated': [], 'acyl_coa_unsaturated': [], 'lyso_PA': [],
                                 'PA': [], 'CDP_DG': [], 'DAG': [], 'TAG': [], 'PS': [], 'PI': [], 'PE': [], 'PC': [], 'CL': [],
                                 'ergosterol': [], 'sterylester': [], 'sphingolipid': []}

        # lists to collect the transported lipids
        # compartment lists
        self.membranes_state = {'plasma_membrane': [], 'secretory_vesicles': [], 'vacuoles': [], 'nucleus': [], 'peroxisomes': [],
                                'light_microsomes': [], 'inner_mit_membrane': [], 'outer_mit_membrane': [], 'lipid_droplets': []}

        # names of membrane lipids
        self.membrane_lipid_names = ['PS', 'PI', 'PC', 'PE', 'CL', 'PA', 'ES', 'SE', 'TAG', 'SL']
        # model output: membrane ratios for every time step
        self.comp_ratio_dict = \
            {comp: dict(zip(self.membrane_lipid_names, [0] * 10)) for comp in self.membranes_state.keys()}

    def __init_precursors(self): 
        # number of small molecules that is produced from anywhere in the cell and will be added every 10 seconds
        # G1 phase
        self.cc_precursors_production = {'G1': {'pyruvate': 1500., 'acetyl_coa': 0, 'glycerol-3-p': 5., 'DHAP': 30.,
                                                'serine': 20., 'glucose_6_p': 8., 'SAM': 45., 'SAH': 0.,
                                                'glycerol_3_p_mito': 5., 'ceramide': 0, 'GDP-mannose': 10, 'NAD': 0,
                                                'NADH': 0, 'NADP': 0, 'NADPH': 0, 'O2': 0, 'H2O': 0, 'CO2': 0, 'Pi': 0,
                                                'CTP': 20, 'CMP': 0, 'inositol': 0, 'ATP': 0, 'ADP': 0},
                                         'S_M': {'pyruvate': 1700., 'acetyl_coa': 0, 'glycerol-3-p': 5., 'DHAP': 35.,
                                                 'serine': 30., 'glucose_6_p': 12., 'SAM': 55., 'SAH': 0.,
                                                 'glycerol_3_p_mito': 5., 'ceramide': 0, 'GDP-mannose': 10, 'NAD': 0,
                                                 'NADH': 0, 'NADP': 0, 'NADPH': 0, 'O2': 0, 'H2O': 0, 'CO2': 0, 'Pi': 0,
                                                 'CTP': 45, 'CMP': 0, 'inositol': 0, 'ATP': 0, 'ADP': 0}}

        # initialise number of small molecules in the cell
        self.precursors_state = {'pyruvate': 2200., 'acetyl_coa': 1000, 'glycerol-3-p': 1000., 'DHAP': 1000.,
                                 'serine': 250., 'glucose_6_p': 1000., 'SAM': 331., 'SAH': 0., 'glycerol_3_p_mito': 50.,
                                 'ceramide': 100, 'GDP-mannose': 0, 'NAD': 0, 'NADH': 0, 'NADP': 0, 'NADPH': 0, 'O2': 0,
                                 'H2O': 0, 'CO2': 0, 'Pi': 0, 'CTP': 1000, 'CMP': 0, 'inositol': 350, 'ATP': 0, 'ADP': 0}

    def __init_molecules(self):
        # initial lipid molecules, based on Zinser et al. (1991, PMID:2002005)
        # 'other' were interpreted as sphingolipids
        plasma_membrane_comp_start = [0.17320, 0.09124, 0.08660, 0.10464, 0.00103, 0.02010, 0.48454, 0.0, 0.0, 0.03557]
        secretory_vesicles_comp_start = [0.08205, 0.11745, 0.20824, 0.13573, 0.01239, 0.01525, 0.42900, 0.0, 0.0, 0.05029]
        vacuoles_comp_start = [0.04817, 0.16604, 0.40517, 0.17537, 0.02442, 0.02866, 0.15200, 0.0, 0.0, 0.06525]
        nucleus_comp_start = [0.04038, 0.09650, 0.27645, 0.16848, 0.01049, 0.01781, 0.390, 0.0, 0.0, 0.02622]
        peroxisomes_comp_start = [0.03235, 0.11360, 0.34656, 0.16465, 0.05033, 0.01150, 0.281, 0.0, 0.0, 0.0]
        light_microsomes_comp_start = [0.05304, 0.06019, 0.40796, 0.26583, 0.00381, 0.00222, 0.206, 0.0, 0.0, 0.00397]
        inner_mit_membrane_comp_start = [0.02880, 0.12273, 0.29107, 0.18192, 0.12204, 0.01137, 0.242, 0.0, 0.0, 0.0]
        outer_mit_membrane_comp_start = [0.01189, 0.10108, 0.45190, 0.32307, 0.05847, 0.04360, 0.009, 0.0, 0.0, 0.0]
        lipid_droplets_comp_start = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.0]
        # list of lists
        self.membrane_compositions_start = {'plasma_membrane': plasma_membrane_comp_start,
                                            'secretory_vesicles': secretory_vesicles_comp_start,
                                            'vacuoles': vacuoles_comp_start,
                                            'nucleus': nucleus_comp_start,
                                            'peroxisomes': peroxisomes_comp_start,
                                            'light_microsomes': light_microsomes_comp_start,
                                            'inner_mit_membrane': inner_mit_membrane_comp_start,
                                            'outer_mit_membrane': outer_mit_membrane_comp_start,
                                            'lipid_droplets': lipid_droplets_comp_start}

        # number of lipids that are produced in the start function for every membrane
        self.start_lipids = {'plasma_membrane': 32950,
                             'secretory_vesicles': 500,
                             'vacuoles': 2500,
                             'nucleus': 6000,
                             'peroxisomes': 500,
                             'light_microsomes': 500,
                             'inner_mit_membrane': 5000,
                             'outer_mit_membrane': 2500,
                             'lipid_droplets': 1000}

    def __init_simulation(self):
        # time point list for plotting
        self.t = None
        
        # lists of all precursors, needed for plotting
        self.precursors_tc = {'pyruvate': [], 'acetyl_coa': [], 'glycerol-3-p': [], 'DHAP': [], 'serine': [],
                              'glucose_6_p': [], 'SAM': [], 'SAH': [], 'glycerol_3_p_mito': [], 'ceramide': [],
                              'GDP-mannose': [], 'NAD': [], 'NADH': [], 'NADP': [], 'NADPH': [], 'O2': [], 'H2O': [],
                              'CO2': [], 'Pi': [], 'CTP': [], 'CMP': [], 'inositol': [], 'ATP': [], 'ADP': []}

        # collecting the products of every time step. Lipids that are produced in the start function don't need the 0 for plotting
        self.number_lipids_tc = {'acyl_coa': [], 'PA': [], 'CDP_DG': [], 'TAG': [], 'PS': [], 'PI': [],
                                 'PE': [], 'PC': [], 'CL': [], 'ergosterol': [], 'sterylester': [], 'DAG': [], 'sphingolipid': []}

        # counting the lipids in each membrane after every time step
        self.number_membranes_tc = {'plasma_membrane': [], 'secretory_vesicles': [], 'vacuoles': [], 'nucleus': [],
                                    'peroxisomes': [], 'light_microsomes': [], 'inner_mit_membrane': [], 'outer_mit_membrane': [],
                                    'lipid_droplets': []}

        self.probabilities = self.initial_probabilities

    def run(self, timesteps=7200):
        sim_time = 0
        self.t = [i for i in range(timesteps)]

        self.start()  # function that produces the lipids and membranes that are existing at the beginning of the cell cycle

        Km = self.__Km_calculation()
        for t in range(timesteps):
            sim_time += 1  # counting the seconds for knowing the cell cycle phase

            # precursor production is cell cycle dependent
            precursors_production = self.cc_precursors_production[self.cell_cycle(sim_time)]

            if sim_time % 10 == 0:
                for key in self.precursors_state:
                    self.precursors_state[key] += precursors_production[key]

            self.probabilities = self.__calculate_threshold(Km)  # manual_threshold

            if self.cell_cycle(sim_time) == 'G1':
                self.probabilities['DAG_synthase'] = 1 - self.cc_probabilities['G1']['DAG_synthase']
                self.probabilities['TAG_synthase'] = 1 - self.cc_probabilities['G1']['TAG_synthase']
                self.probabilities['TAG_lipase'] = 1 - self.cc_probabilities['G1']['TAG_lipase']
                self.probabilities['DAG_kinase'] = 1 - self.cc_probabilities['G1']['DAG_kinase']
            else:
                self.probabilities['DAG_synthase'] = 1 - self.cc_probabilities['S_M']['DAG_synthase']
                self.probabilities['TAG_synthase'] = 1 - self.cc_probabilities['S_M']['TAG_synthase']
                self.probabilities['TAG_lipase'] = 1 - self.cc_probabilities['S_M']['TAG_lipase']
                self.probabilities['DAG_kinase'] = 1 - self.cc_probabilities['S_M']['DAG_kinase']
                self.probabilities['sterylester_synthase'] = 1 - self.cc_probabilities['S_M']['sterylester_synthase']

            function_list = [self.glycerol_3_p_synthesis,
                             self.inositol_synthesis,
                             self.ceramide_synthesis,
                             self.acetyl_coa_synthase,
                             self.acyl_synthase,
                             self.PA_synthesis,
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
            # all reactions that take place during one second in a random order
            for i in function_list:
                func = np.random.choice(function_list)
                func()
                function_list.remove(func)
            self.numbers()  # calculating the produced lipids after each time step
        self.membrane_compositions()  # calculation of the membrane compositions at the end of the run
        self.saturation_counter()  # calculating the percentages of each fatty acid that was used

        return self.saturation_composition_total, self.number_membranes_tc, self.comp_ratio_dict

    def __calculate_threshold(self, Km):
        """
        Calculate current thresholds based on amount of available precursors.

        Parameter
        ---------
        Km: dict
            current Km values for all reactions

        Return
        ------
        threshold: dict
            current thresholds of all reactions
        """
        threshold = {'glycerol_3_p_synthesis': 1 - (self.precursors_state['DHAP'] / (Km['glycerol_3_p_synthesis']['DHAP'] +
                                                    self.precursors_state['DHAP'])),
                     'inositol_synthesis': 1 - (self.precursors_state['glucose_6_p'] / (Km['inositol_synthesis']['glucose_6_p'] +
                                                self.precursors_state['glucose_6_p'])),
                     'ceramide_synthesis': 1 - (self.precursors_state['serine'] / (Km['ceramide_synthesis']['serine'] +
                                                self.precursors_state['serine'])),
                     'acetyl_coa_synthase': 1 - (self.precursors_state['pyruvate'] / (Km['acetyl_coa_synthase']['pyruvate'] +
                                                 self.precursors_state['pyruvate'])),
                     'acyl_synthase': 1 - (self.precursors_state['acetyl_coa'] / (Km['acyl_synthase']['acetyl_coa'] +
                                           self.precursors_state['acetyl_coa'])),
                     'acyl_synthase_C16': 0.625,
                     'acyl_synthase_C18': 0.002,
                     'lyso_PA_synthase': 1 - (((float(len(self.components_state['acyl_coa_saturated'])) +
                                                float(len(self.components_state['acyl_coa_unsaturated']))) /
                                               (Km['lyso_PA_synthase']['acyl_coa'] +
                                                (float(len(self.components_state['acyl_coa_saturated'])) +
                                                 float(len(self.components_state['acyl_coa_unsaturated']))))) *
                                              (self.precursors_state['DHAP'] /
                                               (Km['lyso_PA_synthase']['DHAP'] + self.precursors_state['DHAP']))),
                     'PA_synthase': 1 - ((float(len(self.components_state['lyso_PA'])) / (Km['PA_synthase']['lyso-PA'] +
                                                                           float(len(self.components_state['lyso_PA'])))) *
                                         ((float(len(self.components_state['acyl_coa_saturated'])) +
                                           float(len(self.components_state['acyl_coa_unsaturated']))) /
                                          (Km['PA_synthase']['acyl_coa'] +
                                           (float(len(self.components_state['acyl_coa_saturated'])) +
                                            float(len(self.components_state['acyl_coa_unsaturated'])))))),
                     'CDP_DG_synthase': 1 - ((float(len(self.components_state['PA'])) / (Km['CDP_DG_synthase']['PA'] +
                                                                          float(len(self.components_state['PA'])))) *
                                             (self.precursors_state['CTP'] / (Km['CDP_DG_synthase']['CTP'] +
                                                                             self.precursors_state['CTP']))),
                     'DAG_synthase': 1 - (float(len(self.components_state['PA'])) / (Km['DAG_synthase']['PA'] +
                                                                      float(len(self.components_state['PA'])))),
                     'TAG_synthase': 1 - (((float(len(self.components_state['acyl_coa_saturated'])) +
                                            float(len(self.components_state['acyl_coa_unsaturated']))) /
                                           (Km['TAG_synthase']['acyl_coa'] +
                                            (float(len(self.components_state['acyl_coa_saturated'])) +
                                             float(len(self.components_state['acyl_coa_unsaturated']))))) *
                                          (float(len(self.components_state['DAG'])) /
                                           (Km['TAG_synthase']['DAG'] + float(len(self.components_state['DAG']))))),
                     'TAG_lipase': 1 - (float(len(self.membranes_state['lipid_droplets'])) / (Km['TAG_lipase']['lipid_droplets'] +
                                                                           float(len(self.membranes_state['lipid_droplets'])))),
                     'DAG_kinase': 1 - (float(len(self.components_state['DAG'])) / (Km['DAG_kinase']['DAG'] +
                                                                     float(len(self.components_state['DAG'])))),
                     'PS_synthase': 1 - ((float(len(self.components_state['CDP_DG'])) / (Km['PS_synthase']['CDP_DG'] +
                                                                          float(len(self.components_state['CDP_DG'])))) *
                                         (self.precursors_state['serine'] / (Km['PS_synthase']['serine'] +
                                                                            self.precursors_state['serine']))),
                     'PI_synthase': 1 - ((float(len(self.components_state['CDP_DG'])) / (Km['PI_synthase']['CDP_DG'] +
                                                                          float(len(self.components_state['CDP_DG'])))) *
                                         (self.precursors_state['inositol'] / (Km['PI_synthase']['inositol'] +
                                                                              self.precursors_state['inositol']))),
                     'PE_synthase': 1 - (float(len(self.components_state['PS'])) / (Km['PE_synthase']['PS'] +
                                                                     float(len(self.components_state['PS'])))),
                     'PC_synthase': 1 - ((float(len(self.components_state['PE'])) / (Km['PC_synthase']['PE'] +
                                                                      float(len(self.components_state['PE'])))) *
                                         (self.precursors_state['SAM'] / (Km['PC_synthase']['SAM'] +
                                                                         self.precursors_state['SAM']))),
                     'CL_synthase': 1 - ((float(len(self.components_state['CDP_DG'])) / (Km['CL_synthase']['CDP_DG'] +
                                                                          float(len(self.components_state['CDP_DG'])))) *
                                         (self.precursors_state['glycerol_3_p_mito'] /
                                          (Km['CL_synthase']['glycerol_3_p_mito'] +
                                           self.precursors_state['glycerol_3_p_mito']))),
                     'ergosterol_synthase': 1 - (self.precursors_state['acetyl_coa'] /
                                                 (Km['ergosterol_synthase']['acetyl_coa'] +
                                                  self.precursors_state['acetyl_coa'])),
                     'sterylester_synthase': 1 - ((float(len(self.components_state['ergosterol'])) /
                                                   (Km['sterylester_synthase']['ergosterol'] +
                                                    float(len(self.components_state['ergosterol'])))) *
                                                  ((float(len(self.components_state['acyl_coa_saturated'])) +
                                                    float(len(self.components_state['acyl_coa_unsaturated']))) /
                                                   (Km['sterylester_synthase']['acyl_coa'] +
                                                    (float(len(self.components_state['acyl_coa_saturated'])) +
                                                     float(len(self.components_state['acyl_coa_unsaturated'])))))),
                     'sphingolipid_synthase': 1 - ((float(len(self.components_state['PI'])) /
                                                    (Km['sphingolipid_synthase']['PI'] +
                                                     float(len(self.components_state['PI'])))) *
                                                   (self.precursors_state['ceramide'] /
                                                    (Km['sphingolipid_synthase']['ceramide'] +
                                                     self.precursors_state['ceramide'])))}
        return threshold

    def plot_membranes(self):
        '''
        Plotting the number of lipids in the membranes of different compartments.
        '''
        fig = mat.figure()
        ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
        ax.plot(self.t, self.number_membranes_tc['plasma_membrane'], label = 'plasma membrane')
        ax.plot(self.t, self.number_membranes_tc['secretory_vesicles'], label = 'secretory vesicles')
        ax.plot(self.t, self.number_membranes_tc['vacuoles'], label = 'vacuoles')
        ax.plot(self.t, self.number_membranes_tc['nucleus'], label = 'nucleus')
        ax.plot(self.t, self.number_membranes_tc['peroxisomes'], label = 'peroxisomes')
        ax.plot(self.t, self.number_membranes_tc['light_microsomes'], label = 'light_microsomes')
        ax.plot(self.t, self.number_membranes_tc['inner_mit_membrane'], label = 'inner_mit_membrane')
        ax.plot(self.t, self.number_membranes_tc['outer_mit_membrane'], label = 'outer_mit_membrane')
        ax.plot(self.t, self.number_membranes_tc['lipid_droplets'], label = 'lipid_droplets')
        ax.legend(bbox_to_anchor = (1.05, 1), loc = 2, borderaxespad = 0.)
        mat.show()

    def start(self):
        """
        Function that produces the initial lipids at t = 0.
        Membrane compositions are given in compositions_start.
        Number of starting lipids for each membrane are given in self.start_lipids.
        """
        compositions_start = {}

        # numbers of the comp_start lists should yield 1 when summed
        for mem in self.membrane_compositions_start:
            membrane_comp_start_ratio = [component / sum(self.membrane_compositions_start[mem]) for component in self.membrane_compositions_start[mem]]
            compositions_start[mem] = membrane_comp_start_ratio

        for membrane in self.membranes_state:
            # producing the lipids for a membrane, probability for a certain lipid from the composition in Zinser
            for j in range(self.start_lipids[membrane]):
                head_groups_start = ['serine', 'inositol', 'choline', 'ethanolamine', 'neutral', 'p', 'sterol',
                                     'sterylester', None, 'ceramide']
                weights_start = compositions_start[membrane]
                head = np.random.choice(head_groups_start, p=weights_start)
                if head == 'sterol':
                    new_lipid = components.Sterol(head, membrane, self.weights['compartments'])
                elif head == 'sterylester':
                    new_lipid = components.Sterylester(head, np.random.choice(self.weights['chain_unsaturated'].values(),
                                                                              p=[0.67, 0.33]),membrane,
                                                       self.weights['compartments'])
                elif head == 'ceramide':
                    new_lipid = components.Sphingolipid(head, membrane, self.weights['compartments'])
                elif head == 'neutral':
                    new_lipid = components.CL(head, np.random.choice(self.weights['chain_unsaturated'].values(),
                                                                     p=self.weights['unsaturated']),
                                              np.random.choice(self.weights['chain_complete'],
                                                               p=self.weights['saturation_total']),
                                              np.random.choice(self.weights['chain_unsaturated'].values(),
                                                               p=self.weights['unsaturated']),
                                              np.random.choice(self.weights['chain_complete'],
                                                               p=self.weights['saturation_total']), membrane,
                                              self.weights['compartments'])
                elif head == 'serine' or head == 'inositol' or head == 'choline' or head == 'ethanolamine' or head == 'p':
                    new_lipid = components.Lipid(head, np.random.choice(self.weights['chain_unsaturated'].values(),
                                                                        p=self.weights['unsaturated']),
                                                 np.random.choice(self.weights['chain_complete'],
                                                                  p=self.weights['saturation_total']),
                                                 membrane, self.weights['compartments'])
                else:
                    new_lipid = components.Lipid(head, np.random.choice(self.weights['chain_unsaturated'].values(),
                                                                        p=self.weights['unsaturated']),
                                                 np.random.choice(self.weights['chain_complete'],
                                                                  p=self.weights['saturation_total']), membrane,
                                                 self.weights['compartments'])
                    new_lipid.__class__ = components.TAG
                    new_lipid.sn3 = np.random.choice(self.weights['chain_complete'],
                                                     p=self.weights['saturation_total'])
                self.membranes_state[membrane].append(new_lipid)

        lipid_lists_start = [self.components_state['PS'], self.components_state['PI'], self.components_state['PC'], self.components_state['PE'], self.components_state['CL'], self.components_state['lyso_PA'], self.components_state['PA'],
                             self.components_state['CDP_DG'], self.components_state['ergosterol'], self.components_state['sterylester'], self.components_state['DAG'], self.components_state['TAG'],
                             self.components_state['sphingolipid']]

        for i, lipid_list in enumerate(lipid_lists_start):
            for j in range(20):
                head_groups_start_lipids = ['serine', 'inositol', 'choline', 'ethanolamine', 'neutral', 'lyso', 'p', 'cdp',
                                            'sterol', 'sterylester', 'dag', None, 'ceramide']
                head = head_groups_start_lipids[i]
                if head == 'sterol':
                    new_lipid = components.Sterol(head, None, self.weights['compartments'])
                elif head == 'sterylester':
                    new_lipid = components.Sterylester(head, np.random.choice(self.weights['chain_unsaturated'].values(),
                                                                              p=[0.67, 0.33]), None, self.weights['compartments'])
                elif head == 'neutral':
                    new_lipid = components.CL(head, np.random.choice(self.weights['chain_unsaturated'].values(), p=self.weights['unsaturated']),
                                              np.random.choice(self.weights['chain_complete'], p=self.weights['saturation_total']),
                                              np.random.choice(self.weights['chain_unsaturated'].values(), p=self.weights['unsaturated']),
                                              np.random.choice(self.weights['chain_complete'], p=self.weights['saturation_total']),
                                              None, self.weights['compartments'])
                elif head is None:
                    new_lipid = components.Lipid(head, np.random.choice(self.weights['chain_unsaturated'].values(),
                                                                        p=self.weights['unsaturated']),
                                                 np.random.choice(self.weights['chain_complete'],
                                                                  p=self.weights['saturation_total']), None, self.weights['compartments'])
                    new_lipid.__class__ = components.TAG
                    new_lipid.sn3 = np.random.choice(self.weights['chain_complete'], p=self.weights['saturation_total'])
                elif head == 'ceramide':
                    new_lipid = components.Sphingolipid(head, None, self.weights['compartments'])
                elif head == 'cdp':
                    new_lipid = components.Lipid('p', np.random.choice(self.weights['chain_unsaturated'].values(),
                                                                       p=self.weights['unsaturated']),
                                                 np.random.choice(self.weights['chain_complete'],
                                                                  p=self.weights['saturation_total']), None, self.weights['compartments'])
                elif head == 'lyso':
                    new_lipid = components.Lipid('p', None, np.random.choice(self.weights['chain_complete'],
                                                                             p=self.weights['saturation_total']), None,
                                                 self.weights['compartments'])
                elif head == 'dag':
                    new_lipid = components.Lipid(None, np.random.choice(self.weights['chain_unsaturated'].values(),
                                                                        p=self.weights['unsaturated']),
                                                 np.random.choice(self.weights['chain_complete'],
                                                                  p=self.weights['saturation_total']), None, self.weights['compartments'])
                else:
                    new_lipid = components.Lipid(head, np.random.choice(self.weights['chain_unsaturated'].values(),
                                                                        p=self.weights['unsaturated']),
                                                 np.random.choice(self.weights['chain_complete'],
                                                                  p=self.weights['saturation_total']), None, self.weights['compartments'])
                lipid_list.append(new_lipid)

        choice_list_acyl_start = [0, 1]
        choice_weights_acyl_start = [0.13, 0.87]
        choice_c_acyl_start = [16, 18]

        for i in range(60):
            new_acyl = components.FattyAcid(np.random.choice(choice_c_acyl_start), np.random.choice(choice_list_acyl_start,
                                                                                                    p=choice_weights_acyl_start))
            if new_acyl.saturation == 0:
                self.components_state['acyl_coa_saturated'].append(new_acyl)
            else:
                self.components_state['acyl_coa_unsaturated'].append(new_acyl)
    # def start(self):
    #     '''
    #     Function that produces the starting lipids. Membrane compositions are given in self.membrane_compositions_start_relatives.
    #     Number of starting lipids for each membrane are given in self.start_lipids.
    #     '''
    #
    #     self.membrane_compositions_start_relatives = []
    #
    #     #numbers of the comp_start lists should yield 1 when summed
    #     for membrane_comp_start in self.membrane_compositions_start:
    #         membrane_comp_start_relative = [z / sum(membrane_comp_start) for z in membrane_comp_start]
    #         self.membrane_compositions_start_relatives.append(membrane_comp_start_relative)
    #
    #     self.compositions_start = dict(zip(self.compartment_names, self.membrane_compositions_start_relatives))
    #
    #     x = 0
    #     # self.start_lipids = [32950, 500, 2500, 6000, 500, 500, 5000, 2500, 1000]  # number of lipids that are produced in the start function for every membrane
    #     self.membrane_start = dict(zip(self.compartment_names, self.start_lipids))
    #     for membrane in self.compartment_lists:
    #         for i in range(self.membrane_start[self.compartment_names[x]]):  # #producing the lipids for a membrane, probability for a certain lipid from the composition in Zinser
    #             self.head_groups_start = ['serine', 'inositol', 'choline', 'ethanolamine', 'neutral', 'p', 'sterol', 'sterylester', None, 'ceramide']
    #             weights_start = self.compositions_start[self.compartment_names[x]]
    #             head = np.random.choice(self.head_groups_start, p = weights_start)
    #             if head == 'sterol':
    #                 new_lipid = components.Sterol(head, self.compartment_names[x], self.weights['compartments'])
    #             elif head == 'sterylester':
    #                 new_lipid = components.Sterylester(head, np.random.choice(self.weights['chain_unsaturated'].values(), p = [0.67, 0.33]), self.compartment_names[x], self.weights['compartments'])
    #             elif head == 'ceramide':
    #                 new_lipid = components.Sphingolipid(head, self.compartment_names[x], self.weights['compartments'])
    #             elif head == 'neutral':
    #                 new_lipid = components.CL(head, np.random.choice(self.weights['chain_unsaturated'].values(), p = self.weights['unsaturated']), np.random.choice(self.weights['chain_complete'], p = self.weights['saturation_total']),\
    #                                 np.random.choice(self.weights['chain_unsaturated'].values(), p = self.weights['unsaturated']), np.random.choice(self.weights['chain_complete'], p = self.weights['saturation_total']), self.compartment_names[x], self.weights['compartments'])
    #             elif head == 'serine' or head == 'inositol' or head == 'choline' or head == 'ethanolamine' or head == 'p':
    #                 new_lipid = components.Lipid(head, np.random.choice(self.weights['chain_unsaturated'].values(), p = self.weights['unsaturated']), np.random.choice(self.weights['chain_complete'], p = self.weights['saturation_total']), self.compartment_names[x], self.weights['compartments'])
    #             else:
    #                 new_lipid = components.Lipid(head, np.random.choice(self.weights['chain_unsaturated'].values(), p = self.weights['unsaturated']), np.random.choice(self.weights['chain_complete'], p = self.weights['saturation_total']), self.compartment_names[x], self.weights['compartments'])
    #                 new_lipid.__class__ = components.TAG
    #                 new_lipid.sn3 = np.random.choice(self.weights['chain_complete'], p = self.weights['saturation_total'])
    #             membrane.append(new_lipid)
    #         x += 1
    #
    #
    #     self.lipid_lists_start = [self.components_state['PS'], self.components_state['PI'], self.components_state['PC'], self.components_state['PE'], self.components_state['CL'], self.components_state['lyso_PA'], self.components_state['PA'], self.components_state['CDP_DG'], self.components_state['ergosterol'], self.components_state['sterylester'], self.components_state['DAG'], self.components_state['TAG'], self.components_state['sphingolipid']]
    #     z = 0
    #     for lipid_list in self.lipid_lists_start:
    #         for i in range(20):
    #             self.head_groups_start_lipids = ['serine', 'inositol', 'choline', 'ethanolamine', 'neutral', 'lyso', 'p', 'cdp', 'sterol', 'sterylester', 'dag', None, 'ceramide']
    #             head = self.head_groups_start_lipids[z]
    #             if head == 'sterol':
    #                 new_lipid = components.Sterol(head, None, self.weights['compartments'])
    #             elif head == 'sterylester':
    #                 new_lipid = components.Sterylester(head, np.random.choice(self.weights['chain_unsaturated'].values(), p = [0.67, 0.33]), None, self.weights['compartments'])
    #             elif head == 'neutral':
    #                 new_lipid = components.CL(head, np.random.choice(self.weights['chain_unsaturated'].values(), p = self.weights['unsaturated']), np.random.choice(self.weights['chain_complete'], p = self.weights['saturation_total']),\
    #                                         np.random.choice(self.weights['chain_unsaturated'].values(), p = self.weights['unsaturated']), np.random.choice(self.weights['chain_complete'], p = self.weights['saturation_total']), None, self.weights['compartments'])
    #             elif head == None:
    #                 new_lipid = components.Lipid(head, np.random.choice(self.weights['chain_unsaturated'].values(), p = self.weights['unsaturated']), np.random.choice(self.weights['chain_complete'], p = self.weights['saturation_total']), None, self.weights['compartments'])
    #                 new_lipid.__class__ = components.TAG
    #                 new_lipid.sn3 = np.random.choice(self.weights['chain_complete'], p = self.weights['saturation_total'])
    #             elif head == 'ceramide':
    #                 new_lipid = components.Sphingolipid(head, None, self.weights['compartments'])
    #             elif head == 'cdp':
    #                 new_lipid = components.Lipid('p', np.random.choice(self.weights['chain_unsaturated'].values(), p = self.weights['unsaturated']), np.random.choice(self.weights['chain_complete'], p = self.weights['saturation_total']), None, self.weights['compartments'])
    #             elif head == 'lyso':
    #                 new_lipid = components.Lipid('p', None, np.random.choice(self.weights['chain_complete'], p = self.weights['saturation_total']), None, self.weights['compartments'])
    #             elif head == 'dag':
    #                 new_lipid = components.Lipid(None, np.random.choice(self.weights['chain_unsaturated'].values(), p = self.weights['unsaturated']), np.random.choice(self.weights['chain_complete'], p = self.weights['saturation_total']), None, self.weights['compartments'])
    #             else:
    #                 new_lipid = components.Lipid(head, np.random.choice(self.weights['chain_unsaturated'].values(), p = self.weights['unsaturated']), np.random.choice(self.weights['chain_complete'], p = self.weights['saturation_total']), None, self.weights['compartments'])
    #             lipid_list.append(new_lipid)
    #         z += 1
    #
    #     choice_list_acyl_start = [0, 1]
    #     choice_weights_acyl_start = [0.13, 0.87]
    #     choice_C_acyl_start = [16, 18]
    #     for i in range(60):
    #         new_acyl = components.FattyAcid(np.random.choice(choice_C_acyl_start), np.random.choice(choice_list_acyl_start, p = choice_weights_acyl_start))
    #         if new_acyl.saturation == 0:
    #             self.components_state['acyl_coa_saturated'].append(new_acyl)
    #         else:
    #             self.components_state['acyl_coa_unsaturated'].append(new_acyl)

    def __Km_calculation(self):
        #Calculation of Km, 2 substrates: Km for lipids = 5, for acyl-coa = 30
        Km = {'glycerol_3_p_synthesis': {'DHAP': 0.}, \
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
                        'ergosterol_synthase': {'acetyl_coa': 0.}, \
                        'sterylester_synthase': {'ergosterol': 5., 'acyl_coa': 30.}, \
                        'sphingolipid_synthase': {'PI': 5}}

        Km = {'glycerol_3_p_synthesis': {'DHAP': self.precursors_state['DHAP'] / self.initial_probabilities['glycerol_3_p_synthesis'] - self.precursors_state['DHAP']}, \
                        'inositol_synthesis': {'glucose_6_p': self.precursors_state['glucose_6_p'] / self.initial_probabilities['inositol_synthesis'] - self.precursors_state['glucose_6_p']}, \
                        'ceramide_synthesis': {'serine': self.precursors_state['serine'] / self.initial_probabilities['ceramide_synthesis'] - self.precursors_state['serine']},\
                        'acetyl_coa_synthase': {'pyruvate': self.precursors_state['pyruvate'] / self.initial_probabilities['acetyl_coa_synthase'] - self.precursors_state['pyruvate']}, \
                        'acyl_synthase': {'acetyl_coa': self.precursors_state['acetyl_coa'] / self.initial_probabilities['acyl_synthase'] - self.precursors_state['acetyl_coa']},\
                        'acyl_synthase_C16': 0.625,\
                        'acyl_synthase_C18': 0.002,\
                        'lyso_PA_synthase': {'acyl_coa': 30., 'DHAP': (self.precursors_state['DHAP'] / self.initial_probabilities['lyso_PA_synthase']) * ((float(len(self.components_state['acyl_coa_saturated'])) + float(len(self.components_state['acyl_coa_unsaturated'])))\
                                             / (Km['lyso_PA_synthase']['acyl_coa'] + (float(len(self.components_state['acyl_coa_saturated'])) + float(len(self.components_state['acyl_coa_unsaturated']))))) - self.precursors_state['DHAP']}, \
                        'PA_synthase': {'lyso-PA': 5., 'acyl_coa': 30.},\
                        'CDP_DG_synthase': {'PA': 5., 'CTP': (self.precursors_state['CTP'] / self.initial_probabilities['CDP_DG_synthase']) * (float(len(self.components_state['PA'])) / (Km['CDP_DG_synthase']['PA'] + float(len(self.components_state['PA'])))) - self.precursors_state['CTP']},\
                        'DAG_synthase': {'PA': 5.}, \
                        'TAG_synthase': {'DAG': 5., 'acyl_coa': 30.}, \
                        'TAG_lipase': {'lipid_droplets': float(len(self.membranes_state['lipid_droplets'])) / self.initial_probabilities['TAG_lipase'] - float(len(self.membranes_state['lipid_droplets']))}, \
                        'DAG_kinase': {'DAG': 5.}, \
                        'PS_synthase': {'CDP_DG': 5., 'serine': (self.precursors_state['serine'] / self.initial_probabilities['PS_synthase']) * (float(len(self.components_state['CDP_DG'])) / (Km['PS_synthase']['CDP_DG'] + float(len(self.components_state['CDP_DG'])))) - self.precursors_state['serine']},\
                        'PI_synthase': {'CDP_DG': 5., 'inositol': (self.precursors_state['inositol'] / self.initial_probabilities['PI_synthase']) * (float(len(self.components_state['CDP_DG'])) / (Km['PI_synthase']['CDP_DG'] + float(len(self.components_state['CDP_DG'])))) - self.precursors_state['inositol']}, \
                        'PE_synthase': {'PS': 5.}, \
                        'PC_synthase': {'PE': 5., 'SAM': (self.precursors_state['SAM'] / self.initial_probabilities['PC_synthase']) * (float(len(self.components_state['PE'])) / (Km['PC_synthase']['PE'] + float(len(self.components_state['PE'])))) - self.precursors_state['SAM']}, \
                        'CL_synthase': {'CDP_DG': 5., 'glycerol_3_p_mito': (self.precursors_state['glycerol_3_p_mito'] / self.initial_probabilities['CL_synthase']) * (float(len(self.components_state['CDP_DG'])) / (Km['CL_synthase']['CDP_DG'] + float(len(self.components_state['CDP_DG'])))) - self.precursors_state['glycerol_3_p_mito']},\
                        'ergosterol_synthase': {'acetyl_coa': self.precursors_state['acetyl_coa'] / self.initial_probabilities['ergosterol_synthase'] - self.precursors_state['acetyl_coa']}, \
                        'sterylester_synthase': {'ergosterol': 5., 'acyl_coa': 30.}, \
                        'sphingolipid_synthase': {'PI': 5, 'ceramide': (self.precursors_state['ceramide'] / self.initial_probabilities['sphingolipid_synthase']) * (float(len(self.components_state['PI'])) / (Km['sphingolipid_synthase']['PI'] + float(len(self.components_state['PI'])))) - self.precursors_state['ceramide']}}
        return Km

    def cell_cycle(self, time):
        '''
        Function to determine the cell cycle phases depending on the elapsed time.
        '''
        if time <= 1800:
            return "G1"
        else:
            return "S_M"


    def glycerol_3_p_synthesis(self):
        '''
        Synthesis of glycerol-3-p out of DHAP.
        '''
        for i in range(self.rates['glycerol_3_p_synthesis']):
            x = np.random.random()
            if x >= self.probabilities['glycerol_3_p_synthesis']:
                if self.precursors_state['DHAP'] > 1:
                    self.precursors_state['glycerol-3-p'] += 1
                    self.precursors_state['DHAP'] -= 1
                    self.precursors_state['NADH'] += 1
                    self.precursors_state['NAD'] -= 1


    def inositol_synthesis(self):
        '''
        Synthesis of inositol from glucose-6-p by the Myo-inositol-1-p synthase and the Myo-inositol 1
        phosphatase.
        '''
        for i in range(self.rates['inositol_synthesis']):
            x = np.random.random()
            if x >= self.probabilities['inositol_synthesis']:
                if self.precursors_state['glucose_6_p'] > 1:
                    self.precursors_state['inositol'] += 1
                    self.precursors_state['glucose_6_p'] -= 1
                    self.precursors_state['H2O'] -= 1
                    self.precursors_state['Pi'] += 1


    def ceramide_synthesis(self):
        '''
        Synthesis of ceramide out of serine and a C16:0 fatty acid
        '''
        for i in range(self.rates['ceramide_synthesis']):
            x = np.random.random()
            if x >= self.probabilities['ceramide_synthesis']:
                if len(self.components_state['acyl_coa_C26']) > 1 and self.precursors_state['serine'] > 1 and len(self.components_state['acyl_coa_saturated']) > 1 and any(fa.C == 16 for fa in self.components_state['acyl_coa_saturated']):
                    self.precursors_state['ceramide'] += 1
                    self.precursors_state['serine'] -= 1
                    self.precursors_state['CO2'] += 1
                    self.precursors_state['NADPH'] += 1
                    self.precursors_state['NADP'] -= 1
                    j = 0
                    if self.components_state['acyl_coa_saturated'][j].C == 18 and i <= len(self.components_state['acyl_coa_saturated']):
                        j += 1
                    elif self.components_state['acyl_coa_saturated'][j].C == 16:
                        del self.components_state['acyl_coa_saturated'][j]
                        del self.components_state['acyl_coa_C26'][0]


    def acetyl_coa_synthase(self):
        '''
        Synthesis of Acetyl-CoA: pyruvate dehydrogenase drives the reaction pyruvate to Acetyl-CoA, CO2 is released
        '''
        for i in range(self.rates['acetyl_coa_synthase']):
            x = np.random.random()
            if x >= self.probabilities['acetyl_coa_synthase']:
                if self.precursors_state['pyruvate'] > 1:			# transformation from pyruvate to acetyl_coa
                    self.precursors_state['acetyl_coa'] += 1
                    self.precursors_state['pyruvate'] -= 1
                    self.precursors_state['NADH'] += 1
                    self.precursors_state['NAD'] -= 1
                    self.precursors_state['CO2'] += 1


    def acyl_synthase(self):
        '''
        Simplified synthesis of acyl_coa: several Acetyl-CoA are building a fatty acid (C16:0, C16:1, C18:0 or C18:1)
        The intermediate Malonyl-CoA is leaved out.
        '''
        choice_list = [0, 1]
        choice_weights = [0.12, 0.88]
        for i in range(self.rates['acyl_synthase']):
            x = np.random.random()						#5 reactions in 1 timestep but only with a probability of 90%
            if self.precursors_state['acetyl_coa'] > 2:		#control if at least 2 Acetyl-CoA are available
                if len(self.components_state['acyl_coa']) == 0:		#starting the first reaction
                    new_acyl = components.FattyAcid(2, np.random.choice(choice_list, p = choice_weights))
                    self.components_state['acyl_coa'].append(new_acyl)
                    self.components_state['acyl_coa'][-1].C += 2
                    self.precursors_state['acetyl_coa'] -= 2

                elif self.components_state['acyl_coa'][-1].C == 16 and x >= self.probabilities['acyl_synthase_C16']:	#stop the reaction cycle and starting a new one
                    new_acyl = components.FattyAcid(2, np.random.choice(choice_list, p = choice_weights))
                    self.components_state['acyl_coa'].append(new_acyl)
                    self.components_state['acyl_coa'][-1].C += 2
                    self.precursors_state['acetyl_coa'] -= 2
                    self.precursors_state['NADPH'] -= 14
                    self.precursors_state['NADP'] += 14
                    self.precursors_state['H2O'] += 7
                    #CO2 production is not mentioned here as onyl acetyl-CoA is used and not malonyl-CoA, so we need all C-atoms we give in the reaction

                elif self.components_state['acyl_coa'][-1].C == 18 and x >= self.probabilities['acyl_synthase_C18']:	#stop the reaction cycle and starting a new one
                    new_acyl = components.FattyAcid(2, np.random.choice(choice_list, p = choice_weights))
                    self.components_state['acyl_coa'].append(new_acyl)
                    self.components_state['acyl_coa'][-1].C += 2
                    self.precursors_state['acetyl_coa'] -= 2
                    self.precursors_state['NADPH'] -= 16
                    self.precursors_state['NADP'] += 16
                    self.precursors_state['H2O'] += 8

                elif self.components_state['acyl_coa'][-1].C == 26:
                    self.components_state['acyl_coa'][-1].saturation = 0
                    new_acyl = components.FattyAcid(2, np.random.choice(choice_list, p = choice_weights))
                    self.components_state['acyl_coa'].append(new_acyl)
                    self.components_state['acyl_coa'][-1].C += 2
                    self.precursors_state['acetyl_coa'] -= 2
                    self.precursors_state['NADPH'] -= 24
                    self.precursors_state['NADP'] += 24
                    self.precursors_state['H2O'] += 12

                else:									#adding an Acetyl_CoA to the growing ffa
                    self.components_state['acyl_coa'][-1].C += 2
                    self.precursors_state['acetyl_coa'] -= 1

        if len(self.components_state['acyl_coa']) > 1:
            for j in range(len(self.components_state['acyl_coa'])-1):
                if self.components_state['acyl_coa'][j].C == 26:
                    self.components_state['acyl_coa_C26'].append(self.components_state['acyl_coa'][j])
                elif self.components_state['acyl_coa'][j].saturation == 0:
                    self.components_state['acyl_coa_saturated'].append(self.components_state['acyl_coa'][j])
                elif self.components_state['acyl_coa'][j].saturation == 1:
                    self.components_state['acyl_coa_unsaturated'].append(self.components_state['acyl_coa'][j])
                    self.precursors_state['O2'] -= 1
                    self.precursors_state['H2O'] += 2
            del self.components_state['acyl_coa'][:-1]


    def PA_synthesis(self):
        '''
        Synthesis of PA in two reaction steps.
        '''
        for i in range(self.rates['PA_synthesis']):
            self.lyso_PA_synthase()
            self.PA_synthase()


    def lyso_PA_synthase(self):
        '''
        Production of Lyso-PA by adding one acyl-coa to DHAP (sn1: always unsaturated) --> DHAP acyltransferase/acyl-DHAP reductase
        '''
        choice_list = [0, 1]

        weights_pa = [self.precursors_state['DHAP'] / (self.precursors_state['DHAP'] + self.precursors_state['glycerol-3-p']),\
                    self.precursors_state['glycerol-3-p'] / (self.precursors_state['DHAP'] + self.precursors_state['glycerol-3-p'])]
        x = np.random.random()
        if x >= self.probabilities['lyso_PA_synthase'] and len(self.components_state['acyl_coa_saturated']) > 1 and len(self.components_state['acyl_coa_unsaturated']) > 1 and (self.precursors_state['DHAP'] > 1 and self.precursors_state['glycerol-3-p'] > 1):  # at least 1 ffa has to be unsaturated
            if np.random.choice(choice_list, p = self.weights['FA']) == 0:
                sn1_chain = np.random.randint(0, (len(self.components_state['acyl_coa_saturated'])-1))
                chainlength_sn1 = self.components_state['acyl_coa_saturated'][sn1_chain].C
                lyso_pa = components.Lipid('p', None, self.weights['chain_saturated'][chainlength_sn1], None, self.weights['compartments'])
                del self.components_state['acyl_coa_saturated'][sn1_chain]
            else:
                sn1_chain = np.random.randint(0, (len(self.components_state['acyl_coa_unsaturated'])-1))
                chainlength_sn1 = self.components_state['acyl_coa_unsaturated'][sn1_chain].C
                lyso_pa = components.Lipid('p', None, self.weights['chain_unsaturated'][chainlength_sn1], None, self.weights['compartments'])
                del self.components_state['acyl_coa_unsaturated'][sn1_chain]
            self.components_state['lyso_PA'].append(lyso_pa)
            i = np.random.choice(choice_list, p = weights_pa)
            if i == 0:
                self.precursors_state['DHAP'] -= 1
                self.precursors_state['NADPH'] += 1
                self.precursors_state['NADP'] -= 1
            else:
                self.precursors_state['glycerol-3-p'] -= 1


    def PA_synthase(self):
        '''
        Synthesis of PA by adding the second fatty acid to lyso_PA (sn2: saturated or unsaturated) --> 1-acyl-sn-glycerol-3-phosphate acyltransferase
        '''
        x = np.random.random()
        if x >= self.probabilities['PA_synthase'] and len(self.components_state['acyl_coa_unsaturated']) > 1 and len(self.components_state['lyso_PA']) > 1:
            z = np.random.randint(0, (len(self.components_state['lyso_PA'])-1))
            sn2_chain = np.random.randint(0, (len(self.components_state['acyl_coa_unsaturated'])-1))
            chainlength_sn2 = self.components_state['acyl_coa_unsaturated'][sn2_chain].C
            self.components_state['lyso_PA'][z].sn2 = self.weights['chain_unsaturated'][chainlength_sn2]
            self.components_state['PA'].append(self.components_state['lyso_PA'][z])
            del self.components_state['acyl_coa_unsaturated'][sn2_chain]		# deletion of the consumed ffa
            del self.components_state['lyso_PA'][z]


    def CDP_DG_synthase(self):
        '''
        PA is processed to CDP-DG (CDP-diacylglycerol synthase), that further reacts to the phospholipids
        '''
        for i in range(self.rates['CDP_DG_synthase']):
            x = np.random.random()
            if x >= self.probabilities['CDP_DG_synthase'] and self.precursors_state['CTP'] > 1 and len(self.components_state['PA']) > 1:
                z = np.random.randint(0, len(self.components_state['PA'])-1)
                self.components_state['PA'][z].head = 'cdp'
                self.components_state['CDP_DG'].append(self.components_state['PA'][z])		#CDP-DG production from PA
                del self.components_state['PA'][z]
                self.precursors_state['CTP'] -= 1
                self.precursors_state['Pi'] += 2


    def TAG_synthese(self):
        '''
        Function for TAG synthesis divided in production of DAG and TAG afterwards
        '''
        for i in range(self.rates['TAG_synthesis']):
            self.DAG_synthase()
            self.TAG_synthase()


    def DAG_synthase(self):
        '''
        DAG synthesis: Removing the head of the lipid and adding the lipid to the DAG list.
        '''
        x = np.random.random()
        if x >= self.probabilities['DAG_synthase'] and len(self.components_state['PA']) > 1:
            z = np.random.randint(0, len(self.components_state['PA'])-1)
            self.components_state['PA'][z].head = None
            self.components_state['DAG'].append(self.components_state['PA'][z])
            self.precursors_state['H2O'] -= 1
            self.precursors_state['Pi'] += 1
            del self.components_state['PA'][z]


    def TAG_synthase(self):
        '''
        DAG is processed to TAG by adding a third acyl-chain at position sn3.
        '''
        x = np.random.random()
        if x >= self.probabilities['TAG_synthase'] and len(self.components_state['DAG']) > 1 and len(self.components_state['acyl_coa_saturated']) > 1 and len(self.components_state['acyl_coa_unsaturated']) > 1:
            z = np.random.randint(0, len(self.components_state['DAG'])-1)
            self.components_state['TAG'].append(self.components_state['DAG'][z])
            self.components_state['TAG'][-1].__class__ = components.TAG
            if x <= 0.575:
                sn3 = np.random.randint(0, len(self.components_state['acyl_coa_saturated'])-1)
                chainlength_sn3 = self.components_state['acyl_coa_saturated'][sn3].C
                self.components_state['TAG'][-1].sn3 = self.weights['chain_saturated'][chainlength_sn3]
                del self.components_state['acyl_coa_saturated'][sn3]
            else:
                sn3 = np.random.randint(0, len(self.components_state['acyl_coa_unsaturated'])-1)
                chainlength_sn3 = self.components_state['acyl_coa_unsaturated'][sn3].C
                self.components_state['TAG'][-1].sn3 = self.weights['chain_unsaturated'][chainlength_sn3]
                del self.components_state['acyl_coa_unsaturated'][sn3]
            del self.components_state['DAG'][z]


    def TAG_lipase(self):
        '''
        Cdk1/Cdc28-dependent activation of the major triacylglycerol lipase
        '''
        if len(self.membranes_state['lipid_droplets']) > self.rates['TAG_lipase']:
            for i in range(self.rates['TAG_lipase']):
                x = np.random.random()
                if x >= self.probabilities['TAG_lipase']:
                    z = np.random.randint(0, len(self.membranes_state['lipid_droplets'])-1)
                    if self.membranes_state['lipid_droplets'][z].head == None:
                        if ':0' in self.membranes_state['lipid_droplets'][z].sn3:
                            for key, value in self.weights['chain_unsaturated'].items():
                                if value == self.membranes_state['lipid_droplets'][z].sn3:
                                    self.components_state['acyl_coa_saturated'].append(components.FattyAcid(key, 0))
                        elif ':1' in self.membranes_state['lipid_droplets'][z].sn3:
                            for key, value in self.weights['chain_saturated'].items():
                                if value == self.membranes_state['lipid_droplets'][z].sn3:
                                    self.components_state['acyl_coa_unsaturated'].append(components.FattyAcid(key, 1))
                        self.components_state['DAG'].append(self.membranes_state['lipid_droplets'][z])
                        self.components_state['DAG'][-1].__class__ = components.Lipid
                        delattr(self.components_state['DAG'][-1], 'sn3')
                        self.precursors_state['H2O'] -= 1
                    elif self.membranes_state['lipid_droplets'][z].head == 'sterylester':
                        self.components_state['ergosterol'].append(components.Sterol('sterol', None, self.weights['compartments']))
                        self.precursors_state['H2O'] -= 1
                        if ':0' in self.membranes_state['lipid_droplets'][z].FA:
                            for key, value in self.weights['chain_unsaturated'].items():
                                if value == self.membranes_state['lipid_droplets'][z].FA:
                                    self.components_state['acyl_coa_saturated'].append(components.FattyAcid(key, 0))
                        elif ':1' in self.membranes_state['lipid_droplets'][z].FA:
                            for key, value in self.weights['chain_saturated'].items():
                                if value == self.membranes_state['lipid_droplets'][z].FA:
                                    self.components_state['acyl_coa_unsaturated'].append(components.FattyAcid(key, 1))
                    del self.membranes_state['lipid_droplets'][z]


    def DAG_kinase(self):
        if len(self.components_state['DAG']) > self.rates['DAG_kinase']:
            for i in range(self.rates['DAG_kinase']):
                x = np.random.random()
                if x >= self.probabilities['DAG_kinase']:
                    z = np.random.randint(0, len(self.components_state['DAG'])-1)
                    self.components_state['PA'].append(self.components_state['DAG'][z])
                    self.components_state['PA'][-1].head = 'p'
                    self.components_state['PA'][-1].comp = None
                    del self.components_state['DAG'][z]


    def PS_synthase(self):
        '''
        CDP-DG is processed to PS (PS synthase).
        '''
        for i in range(self.rates['PS_synthase']):
            x = np.random.random()
            if x >= self.probabilities['PS_synthase'] and len(self.components_state['CDP_DG']) > 1 and self.precursors_state['serine'] > 1:
                z = np.random.randint(0, len(self.components_state['CDP_DG'])-1)
                self.components_state['CDP_DG'][z].head = 'serine'				#PS synthesis from CDP-DG
                self.components_state['PS'].append(self.components_state['CDP_DG'][z])
                del self.components_state['CDP_DG'][z]
                self.precursors_state['serine'] -= 1
                self.precursors_state['CMP'] += 1


    def PI_synthase(self):
        '''
        CDP-DG is processed to PI (PI synthase)
        '''
        for i in range(self.rates['PI_synthase']):
            x = np.random.random()
            if x >= self.probabilities['PI_synthase'] and len(self.components_state['CDP_DG']) > 1 and self.precursors_state['inositol'] > 1:
                z = np.random.randint(0, len(self.components_state['CDP_DG'])-1)
                self.components_state['CDP_DG'][z].head = 'inositol'
                self.components_state['PI'].append(self.components_state['CDP_DG'][z])
                del self.components_state['CDP_DG'][z]
                self.precursors_state['inositol'] -= 1
                self.precursors_state['CMP'] += 1


    def PE_synthase(self):
        '''
        PE is derived from PS by releasing 1 CO2 --> PS decarboxylase.
        '''
        for i in range(self.rates['PE_synthase']):
            x = np.random.random()
            if x >= self.probabilities['PE_synthase'] and len(self.components_state['PS']) >= 10:
                z = np.random.randint(0, len(self.components_state['PS'])-1)
                self.components_state['PS'][z].head = 'ethanolamine'
                self.components_state['PE'].append(self.components_state['PS'][z])
                self.precursors_state['CO2'] += 1
                del self.components_state['PS'][z]


    def PC_synthase(self):
        '''
        PC is derived from PE. As enzymes serve 3 methyltransferases which need SAM and produce SAH as a side product.
        '''
        for i in range(self.rates['PC_synthase']):
            x = np.random.random()
            if x >= self.probabilities['PC_synthase'] and len(self.components_state['PE']) >= 5 and self.precursors_state['SAM'] >= 4:
                z = np.random.randint(0, len(self.components_state['PE'])-1)
                self.components_state['PE'][z].head = 'choline'
                self.components_state['PC'].append(self.components_state['PE'][z])
                del self.components_state['PE'][z]
                self.precursors_state['SAM'] -= 3
                self.precursors_state['SAH'] += 3


    def CL_synthase(self):
        '''
        Synthesis of cardiolipin, for which 2 CDP-DG are needed. Different enzymes are needed.
        '''
        for i in range(self.rates['CL_synthase']):
            x = np.random.random()
            if x >= self.probabilities['CL_synthase'] and self.precursors_state['glycerol_3_p_mito'] > 1 and len(self.components_state['CDP_DG']) > 2:
                z = np.random.randint(0, len(self.components_state['CDP_DG'])-2)
                self.components_state['CDP_DG'][z].head = 'neutral'
                self.components_state['CL'].append(self.components_state['CDP_DG'][z])
                self.components_state['CL'][-1].__class__ = components.CL
                self.components_state['CL'][-1].sn4, self.components_state['CL'][-1].sn3 = self.components_state['CDP_DG'][z+1].sn2, self.components_state['CDP_DG'][z+1].sn1
                del self.components_state['CDP_DG'][z:z+1]
                self.precursors_state['glycerol_3_p_mito'] -= 1
                self.precursors_state['H2O'] -= 1
                self.precursors_state['Pi'] += 1
                self.precursors_state['CMP'] += 2


    def Ergosterol_synthase(self):
        '''
        Synthesis of the most existing sterol in yeast: ergosterol
        '''
        for i in range(self.rates['ergosterol_synthase']):
            x = np.random.random()
            if x >= self.probabilities['ergosterol_synthase'] and self.precursors_state['acetyl_coa'] > 18:
                self.components_state['ergosterol'].append(components.Sterol('sterol', None, self.weights['compartments']))
                self.precursors_state['acetyl_coa'] -= 18
                self.precursors_state['ATP'] -= 3
                self.precursors_state['ADP'] += 3
                self.precursors_state['NADPH'] -= 11
                self.precursors_state['NADP'] += 11
                self.precursors_state['SAM'] -= 1
                self.precursors_state['SAH'] += 1
                self.precursors_state['O2'] -= 7
                self.precursors_state['H2O'] += 9
                self.precursors_state['CO2'] += 2


    def Sterylester_synthase(self):
        '''
        Synthesis of sterylesters that are found in lipid droplets out of ergosterol and an unsaturated fatty acid.
        '''
        for i in range(self.rates['sterylester_synthase']):
            x = np.random.random()
            if x >= self.probabilities['sterylester_synthase'] and any(fa.C == 16 for fa in self.components_state['acyl_coa_unsaturated']) and any(fa.C == 18 for fa in self.components_state['acyl_coa_unsaturated']) and len(self.components_state['ergosterol']) > 1:
                z = np.random.randint(0, len(self.components_state['ergosterol'])-1)
                j = 0
                while j < 5:
                    fa_index = np.random.randint(0, len(self.components_state['acyl_coa_unsaturated'])-1)
                    if self.components_state['acyl_coa_unsaturated'][fa_index].C == 18 and np.random.random() < 0.33:
                        self.components_state['sterylester'].append(components.Sterylester('sterylester', 'C18:1', None, self.weights['compartments']))
                        del self.components_state['ergosterol'][z]
                        del self.components_state['acyl_coa_unsaturated'][fa_index]
                        break
                    elif self.components_state['acyl_coa_unsaturated'][fa_index].C == 16:
                        self.components_state['sterylester'].append(components.Sterylester('sterylester', 'C16:1', None, self.weights['compartments']))
                        del self.components_state['ergosterol'][z]
                        del self.components_state['acyl_coa_unsaturated'][fa_index]
                        break
                    else:
                        j += 1


    def Sphingolipid_synthase(self):
        '''
        Synthesis of the most abundant Sphingolipid mannose-(inositol-phosphate)2-ceramide
        '''
        for i in range(self.rates['sphingolipid_synthase']):
            x = np.random.random()
            if x >= self.probabilities['sphingolipid_synthase'] and len(self.components_state['PI']) >= 2 and self.precursors_state['ceramide'] > 1 and self.precursors_state['GDP-mannose'] > 1:
                self.components_state['sphingolipid'].append(components.Sphingolipid('ceramide', None, self.weights['compartments']))
                z= np.random.randint(0, len(self.components_state['PI'])-2)
                del self.components_state['PI'][z:z+1]
                self.precursors_state['ceramide'] -= 1
                self.precursors_state['GDP-mannose'] -= 1


    def transport(self):
        '''
        General transport function for all produced lipids.
        '''
        # lipids to transport
        transport_lists = [self.components_state['PS'], self.components_state['PI'], self.components_state['PC'], self.components_state['PE'], self.components_state['CL'], self.components_state['PA'],
                           self.components_state['ergosterol'], self.components_state['sterylester'], self.components_state['TAG'], self.components_state['sphingolipid']]

        for lipid in transport_lists:
            if lipid == self.components_state['TAG'] or lipid == self.components_state['sterylester']:
                if len(lipid) > 10:
                    for j in range(len(lipid)/10):
                        z = np.random.randint(0, len(lipid)-1)
                        lipid[z].comp_choice()
                        if lipid[z].comp == 'lipid_droplets':
                            self.membranes_state['lipid_droplets'].append(lipid[z])
                        del lipid[z]
            else:
                if len(lipid) > 5:
                    for j in range(len(lipid)/10):
                        z = np.random.randint(0, len(lipid)-1)
                        lipid[z].comp_choice()
                        if lipid[z].comp == 'plasma_membrane':
                            self.membranes_state['plasma_membrane'].append(lipid[z])
                        elif lipid[z].comp == 'secretory_vesicles':
                            self.membranes_state['secretory_vesicles'].append(lipid[z])
                        elif lipid[z].comp == 'vacuoles':
                            self.membranes_state['vacuoles'].append(lipid[z])
                        elif lipid[z].comp == 'nucleus':
                            self.membranes_state['nucleus'].append(lipid[z])
                        elif lipid[z].comp == 'peroxisomes':
                            self.membranes_state['peroxisomes'].append(lipid[z])
                        elif lipid[z].comp == 'light_microsomes':
                            self.membranes_state['light_microsomes'].append(lipid[z])
                        elif lipid[z].comp == 'inner_mit_membrane':
                            self.membranes_state['inner_mit_membrane'].append(lipid[z])
                        elif lipid[z].comp == 'outer_mit_membrane':
                            self.membranes_state['outer_mit_membrane'].append(lipid[z])
                        elif lipid[z].comp == 'lipid_droplets':
                            self.membranes_state['lipid_droplets'].append(lipid[z])
                        del lipid[z]


    def membrane_compositions(self):
        '''
        Function to calculate the lipid composition of all membranes.
        '''
        x = 0
        for comp in self.membranes_state:
            if len(self.membranes_state[comp]) > 0:
                comp_ratio_list = [(float(sum(j.head == 'serine' for j in self.membranes_state[comp])) / len(self.membranes_state[comp])),\
                                        (float(sum(j.head == 'inositol' for j in self.membranes_state[comp])) / len(self.membranes_state[comp])),\
                                        (float(sum(j.head == 'choline' for j in self.membranes_state[comp])) / len(self.membranes_state[comp])),\
                                        (float(sum(j.head == 'ethanolamine' for j in self.membranes_state[comp])) / len(self.membranes_state[comp])),\
                                        (float(sum(j.head == 'neutral' for j in self.membranes_state[comp])) / len(self.membranes_state[comp])),\
                                        (float(sum(j.head == 'p' for j in self.membranes_state[comp])) / len(self.membranes_state[comp])),\
                                        (float(sum(j.head == 'sterol' for j in self.membranes_state[comp])) / len(self.membranes_state[comp])),\
                                        (float(sum(j.head == 'sterylester' for j in self.membranes_state[comp])) / len(self.membranes_state[comp])),\
                                        (float(sum(j.head == None for j in self.membranes_state[comp])) / len(self.membranes_state[comp])),\
                                        (float(sum(j.head == 'ceramide' for j in self.membranes_state[comp])) / len(self.membranes_state[comp]))]
                for i in range(len(comp_ratio_list)):
                    self.comp_ratio_dict[comp][self.membrane_lipid_names[i]] = comp_ratio_list[i]
            x += 1

    def numbers(self):
        # component list
        components_list = ['acyl_coa', 'PA', 'CDP_DG', 'TAG', 'PS', 'PI', 'PE', 'PC', 'CL', 'ergosterol',
                           'sterylester', 'DAG', 'sphingolipid']

        for sp in components_list:
            self.number_lipids_tc[sp].append(len(self.components_state[sp]))

        for com in self.membranes_state:
            self.number_membranes_tc[com].append(len(self.membranes_state[com]))

        for precursor in self.precursors_tc.keys():
            self.precursors_tc[precursor].append(self.precursors_state[precursor])

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
        for c in self.membranes_state:
            if c == 'lipid_droplets':
                continue
            else:
                for i in range(len(self.membranes_state[c])):
                    if hasattr(self.membranes_state[c][i], 'sn1'):
                        if self.membranes_state[c][i].sn1 == 'C16:0':
                            self.c16_0_sn1 += 1
                        elif self.membranes_state[c][i].sn1 == 'C16:1':
                            self.c16_1_sn1 += 1
                        elif self.membranes_state[c][i].sn1 == 'C18:0':
                            self.c18_0_sn1 += 1
                        elif self.membranes_state[c][i].sn1 == 'C18:1':
                            self.c18_1_sn1 += 1
                        else:
                            self.wrong_fatty_acid += 1
        self.saturation_composition_sn1 = {'C16:0': self.c16_0_sn1, 'C16:1': self.c16_1_sn1, 'C18:0': self.c18_0_sn1, 'C18_1': self.c18_1_sn1}

        self.c16_0_sn2 = 0
        self.c16_1_sn2 = 0
        self.c18_0_sn2 = 0
        self.c18_1_sn2 = 0
        self.wrong_fatty_acid = 0
        for c in self.membranes_state:
            if c == 'lipid_droplets':
                continue
            else:
                for i in range(len(self.membranes_state[c])):
                    if hasattr(self.membranes_state[c][i], 'sn2'):
                        if self.membranes_state[c][i].sn2 == 'C16:0':
                            self.c16_0_sn2 += 1
                        elif self.membranes_state[c][i].sn2 == 'C16:1':
                            self.c16_1_sn2 += 1
                        elif self.membranes_state[c][i].sn2 == 'C18:0':
                            self.c18_0_sn2 += 1
                        elif self.membranes_state[c][i].sn2 == 'C18:1':
                            self.c18_1_sn2 += 1
                        else:
                            self.wrong_fatty_acid += 1

        self.saturation_composition_sn2 = {'C16:0': self.c16_0_sn2, 'C16:1': self.c16_1_sn2, 'C18:0': self.c18_0_sn2, 'C18_1': self.c18_1_sn2}
        self.total_fatty_acids = self.c16_0_sn1 + self.c16_1_sn1 + self.c18_0_sn1 + self.c18_1_sn1 + self.c16_0_sn2 + self.c16_1_sn2 + self.c18_0_sn2 + self.c18_1_sn2

        self.sterylester_C16 = 0
        self.sterylester_C18 = 0
        for c in self.membranes_state['lipid_droplets']:
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
    m = Model()
    # test run: 5 sec
    r, mem, s = m.run()
    et = time.time()
    print len(m.t)
    for lili in m.number_lipids_tc:
        mat.plot(m.t, m.number_lipids_tc[lili])
    print "Runtime: " + str(et - st) + "s"
    mat.show()

    for lili in mem:
        mat.plot(m.t, mem[lili])
    mat.show()

    print m.comp_ratio_dict