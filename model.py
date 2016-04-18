# -*- coding: utf-8 -*-
"""
Created/Started on Wed June 03 2015

@author: Vera
"""
import matplotlib.pyplot as mat
import numpy as np
import components
import reactions


class Model:
    """
    The model.
    At the beginning there are several lists defined which will contain the produced lipids.
    The model starts with a given amount of precursors (pyruvate, dhap, ctp, serine, glucose-6-p, SAM, glycerol-3-p).
    After several reactions and the transport function in the end there are membranes of different compartments 
    with several different lipids.
    """
    def __init__(self):
        # initialise model
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
        self.RATES = {'glycerol_3_p_synthesis': 8,
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
        self.INITIAL_PROBABILITIES = {'glycerol_3_p_synthesis': 0.5,
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
        self.CC_PROBABILITIES = {'G1': {'DAG_synthase': 0.3,
                                        'TAG_synthase': 0.2,
                                        'TAG_lipase': 0.05,
                                        'DAG_kinase': 0.03},
                                 'S_M': {'DAG_synthase': 0.01,
                                         'TAG_synthase': 0.2,
                                         'TAG_lipase': 0.6,
                                         'DAG_kinase': 0.1,
                                         'sterylester_synthase': 0.2}}

        # ratios and weights of FA properties and membrane compositions
        self.WEIGHTS = {  # compartment size ratio (molecules) - transport probabilities: adjusted manually
                          # PM, SecVec, Vac, Nuc, Perox, lightMic, MitoMemIn, MitoMemOut
                        'compartments': [0.67, 0.01, 0.155, 0.101, 0.007, 0.007, 0.03, 0.015],

                        'chain_saturated': {16: 'C16:0', 18: 'C18:0'},
                        # possible fatty acids and weights to reach the biological proportions
                        'chain_unsaturated': {16: 'C16:1', 18: 'C18:1'},
                        'chain_complete': ['C16:0', 'C18:0', 'C16:1', 'C18:1'],
                        'FA': [0.4, 0.6],
                        # C atom fraction C16, C18
                        'unsaturated': [0.375, 0.625],
                        'saturation_total': [0.2, 0.2, 0.167, 0.433]
                        # sn1: C16:0, C18:0, C16:1, C18:1
                        }

    def __init_lists(self):
        """
        Initialise state vectors for components in the model.
        """
        # component state vector of the model
        self.components_state = {'acyl_coa': [], 'acyl_coa_C26': [], 'acyl_coa_saturated': [], 'acyl_coa_unsaturated': [], 'lyso_PA': [],
                                 'PA': [], 'CDP_DG': [], 'DAG': [], 'TAG': [], 'PS': [], 'PI': [], 'PE': [], 'PC': [], 'CL': [],
                                 'ergosterol': [], 'sterylester': [], 'sphingolipid': []}

        # membrane state vector of the model
        self.membranes_state = {'plasma_membrane': [], 'secretory_vesicles': [], 'vacuoles': [], 'nucleus': [], 'peroxisomes': [],
                                'light_microsomes': [], 'inner_mit_membrane': [], 'outer_mit_membrane': [], 'lipid_droplets': []}

        # names of membrane lipids
        self.MEMBRANE_LIPID_NAMES = ['PS', 'PI', 'PC', 'PE', 'CL', 'PA', 'ES', 'SE', 'TAG', 'SL']
        # model output: membrane ratios for every time step
        self.comp_ratio_dict = \
            {comp: dict(zip(self.MEMBRANE_LIPID_NAMES, [0] * 10)) for comp in self.membranes_state.keys()}

    def __init_precursors(self):
        """
        Initialise precursors and precursor production.
        """
        # number of small molecules that is produced from anywhere in the cell and will be added every 10 time steps
        # cell cycle dependent, manually adjusted
        self.CC_PRECURSORS_PRODUCTION = {'G1': {'pyruvate': 1500., 'acetyl_coa': 0, 'glycerol-3-p': 5., 'DHAP': 30.,
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
        """
        Initialise component numbers of the model. Based on experimental data.
        Info
        ----
        Data taken from Zinser et al. (1991, PMID:2002005).
        'Other' were interpreted as sphingolipids
        """
        # initial lipid molecules
        PLASMA_MEMBRANE_COMP_START = [0.17320, 0.09124, 0.08660, 0.10464, 0.00103, 0.02010, 0.48454, 0.0, 0.0, 0.03557]
        SECRETORY_VESICLES_COMP_START = [0.08205, 0.11745, 0.20824, 0.13573, 0.01239, 0.01525, 0.42900, 0.0, 0.0, 0.05029]
        VACUOLES_COMP_START = [0.04817, 0.16604, 0.40517, 0.17537, 0.02442, 0.02866, 0.15200, 0.0, 0.0, 0.06525]
        NUCLEUS_COMP_START = [0.04038, 0.09650, 0.27645, 0.16848, 0.01049, 0.01781, 0.390, 0.0, 0.0, 0.02622]
        PEROXISOMES_COMP_START = [0.03235, 0.11360, 0.34656, 0.16465, 0.05033, 0.01150, 0.281, 0.0, 0.0, 0.0]
        LIGHT_MICROSOMES_COMP_START = [0.05304, 0.06019, 0.40796, 0.26583, 0.00381, 0.00222, 0.206, 0.0, 0.0, 0.00397]
        INNER_MIT_MEMBRANE_COMP_START = [0.02880, 0.12273, 0.29107, 0.18192, 0.12204, 0.01137, 0.242, 0.0, 0.0, 0.0]
        OUTER_MIT_MEMBRANE_COMP_START = [0.01189, 0.10108, 0.45190, 0.32307, 0.05847, 0.04360, 0.009, 0.0, 0.0, 0.0]
        LIPID_DROPLETS_COMP_START = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.0]
        # data in dict
        self.MEMBRANE_COMPOSITIONS_START = {'plasma_membrane': PLASMA_MEMBRANE_COMP_START,
                                            'secretory_vesicles': SECRETORY_VESICLES_COMP_START,
                                            'vacuoles': VACUOLES_COMP_START,
                                            'nucleus': NUCLEUS_COMP_START,
                                            'peroxisomes': PEROXISOMES_COMP_START,
                                            'light_microsomes': LIGHT_MICROSOMES_COMP_START,
                                            'inner_mit_membrane': INNER_MIT_MEMBRANE_COMP_START,
                                            'outer_mit_membrane': OUTER_MIT_MEMBRANE_COMP_START,
                                            'lipid_droplets': LIPID_DROPLETS_COMP_START}

        # number of lipids that are produced in the start function for every membrane
        self.START_LIPIDS = {'plasma_membrane': 32950,
                             'secretory_vesicles': 500,
                             'vacuoles': 2500,
                             'nucleus': 6000,
                             'peroxisomes': 500,
                             'light_microsomes': 500,
                             'inner_mit_membrane': 5000,
                             'outer_mit_membrane': 2500,
                             'lipid_droplets': 1000}

    def __init_simulation(self):
        """
        Initialise simulation and time course vectors.
        """
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

        # initialise probabilities with initial values
        self.probabilities = self.INITIAL_PROBABILITIES
        # initialise reactions of the model
        self.reactions = reactions.Reactions(self.RATES, self.WEIGHTS)

    def run(self, timesteps=7200):
        """
        Start simulation.

        Parameter
        ---------
        timesteps: int
            number of time steps in seconds
        """
        # current model time
        sim_time = 0
        # time range of the model, simulation time
        self.t = [i for i in range(timesteps)]

        # function that produces the lipids and membranes that are existing at the beginning of the cell cycle
        self.start()
        # Calculate Km values of the reactions
        Km = self.__Km_calculation()

        # actual simulation algorithm
        for t in range(timesteps):
            # simulation time
            sim_time += 1

            # precursor production is cell cycle dependent
            precursors_production = self.CC_PRECURSORS_PRODUCTION[self.cell_cycle(sim_time)]

            # add precursors every 10 seconds
            #if sim_time % 10 == 0:
            for key in self.precursors_state:
                self.precursors_state[key] += precursors_production[key]/10.

            # recalculate probabilities based on substrate numbers
            self.probabilities = self.__calculate_threshold(Km)

            # overwrite probabilities with cell cycle dependent values
            if self.cell_cycle(sim_time) == 'G1':
                self.probabilities['DAG_synthase'] = 1 - self.CC_PROBABILITIES['G1']['DAG_synthase']
                self.probabilities['TAG_synthase'] = 1 - self.CC_PROBABILITIES['G1']['TAG_synthase']
                self.probabilities['TAG_lipase'] = 1 - self.CC_PROBABILITIES['G1']['TAG_lipase']
                self.probabilities['DAG_kinase'] = 1 - self.CC_PROBABILITIES['G1']['DAG_kinase']
            else:
                self.probabilities['DAG_synthase'] = 1 - self.CC_PROBABILITIES['S_M']['DAG_synthase']
                self.probabilities['TAG_synthase'] = 1 - self.CC_PROBABILITIES['S_M']['TAG_synthase']
                self.probabilities['TAG_lipase'] = 1 - self.CC_PROBABILITIES['S_M']['TAG_lipase']
                self.probabilities['DAG_kinase'] = 1 - self.CC_PROBABILITIES['S_M']['DAG_kinase']
                self.probabilities['sterylester_synthase'] = 1 - self.CC_PROBABILITIES['S_M']['sterylester_synthase']

            # call reactions of the model
            self.reactions.run_step(self.probabilities, self.components_state, self.membranes_state, self.precursors_state)
            # update model state
            self.components_state = self.reactions.components_state
            self.membranes_state = self.reactions.membranes_state
            self.precursors_state = self.reactions.precursors_state
            # calculating the produced lipids after each time step
            self.numbers()
        # calculation of the membrane compositions at the end of the run
        self.membrane_compositions()
        # calculating the percentages of each fatty acid that was used
        self.saturation_counter()

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

    def start(self):
        """
        Function that produces the initial lipids at t = 0.
        Membrane compositions are given in MEMBRANE_COMPOSITION_START.
        Initial number of lipids for each membrane are given in self.START_LIPIDS.
        """
        compositions_start = {}

        # numbers of the comp_start lists should yield 1 when summed
        # correction of error to ensure sentence above
        # TODO: use original Zinser data and remove section: compostiion_start = self.MEMBRANE_COMPOSITIONS_START
        for mem in self.MEMBRANE_COMPOSITIONS_START:
            membrane_comp_start_ratio = [component / sum(self.MEMBRANE_COMPOSITIONS_START[mem]) for component in self.MEMBRANE_COMPOSITIONS_START[mem]]
            compositions_start[mem] = membrane_comp_start_ratio

        # generate initial membranes, lipids associated with membranes
        for membrane in self.membranes_state:
            # producing the lipids for a membrane, probability for a certain lipid from the composition in Zinser
            for j in range(self.START_LIPIDS[membrane]):
                head_groups_start = ['serine', 'inositol', 'choline', 'ethanolamine', 'neutral', 'p', 'sterol', 'sterylester', None, 'ceramide']
                weights_start = compositions_start[membrane]
                head = np.random.choice(head_groups_start, p=weights_start)
                if head == 'sterol':
                    new_lipid = components.Sterol(head, membrane, self.WEIGHTS['compartments'])
                elif head == 'sterylester':
                    new_lipid = components.Sterylester(head, np.random.choice(self.WEIGHTS['chain_unsaturated'].values(),
                                                                              p=[0.67, 0.33]),membrane,
                                                       self.WEIGHTS['compartments'])
                elif head == 'ceramide':
                    new_lipid = components.Sphingolipid(head, membrane, self.WEIGHTS['compartments'])
                elif head == 'neutral':
                    new_lipid = components.CL(head, np.random.choice(self.WEIGHTS['chain_unsaturated'].values(),
                                                                     p=self.WEIGHTS['unsaturated']),
                                              np.random.choice(self.WEIGHTS['chain_complete'],
                                                               p=self.WEIGHTS['saturation_total']),
                                              np.random.choice(self.WEIGHTS['chain_unsaturated'].values(),
                                                               p=self.WEIGHTS['unsaturated']),
                                              np.random.choice(self.WEIGHTS['chain_complete'],
                                                               p=self.WEIGHTS['saturation_total']), membrane,
                                              self.WEIGHTS['compartments'])
                elif head == 'serine' or head == 'inositol' or head == 'choline' or head == 'ethanolamine' or head == 'p':
                    new_lipid = components.Lipid(head, np.random.choice(self.WEIGHTS['chain_unsaturated'].values(),
                                                                        p=self.WEIGHTS['unsaturated']),
                                                 np.random.choice(self.WEIGHTS['chain_complete'],
                                                                  p=self.WEIGHTS['saturation_total']),
                                                 membrane, self.WEIGHTS['compartments'])
                else:
                    new_lipid = components.Lipid(head, np.random.choice(self.WEIGHTS['chain_unsaturated'].values(),
                                                                        p=self.WEIGHTS['unsaturated']),
                                                 np.random.choice(self.WEIGHTS['chain_complete'],
                                                                  p=self.WEIGHTS['saturation_total']), membrane,
                                                 self.WEIGHTS['compartments'])
                    new_lipid.__class__ = components.TAG
                    new_lipid.sn3 = np.random.choice(self.WEIGHTS['chain_complete'],
                                                     p=self.WEIGHTS['saturation_total'])
                self.membranes_state[membrane].append(new_lipid)

        # generate 'free' lipids, not yet associated with membranes
        lipid_lists_start = [self.components_state['PS'], self.components_state['PI'], self.components_state['PC'], self.components_state['PE'], self.components_state['CL'], self.components_state['lyso_PA'], self.components_state['PA'],
                             self.components_state['CDP_DG'], self.components_state['ergosterol'], self.components_state['sterylester'], self.components_state['DAG'], self.components_state['TAG'],
                             self.components_state['sphingolipid']]

        for i, lipid_list in enumerate(lipid_lists_start):
            for j in range(20):
                head_groups_start_lipids = ['serine', 'inositol', 'choline', 'ethanolamine', 'neutral', 'lyso', 'p', 'cdp',
                                            'sterol', 'sterylester', 'dag', None, 'ceramide']
                head = head_groups_start_lipids[i]
                if head == 'sterol':
                    new_lipid = components.Sterol(head, None, self.WEIGHTS['compartments'])
                elif head == 'sterylester':
                    new_lipid = components.Sterylester(head, np.random.choice(self.WEIGHTS['chain_unsaturated'].values(),
                                                                              p=[0.67, 0.33]), None, self.WEIGHTS['compartments'])
                elif head == 'neutral':
                    new_lipid = components.CL(head, np.random.choice(self.WEIGHTS['chain_unsaturated'].values(), p=self.WEIGHTS['unsaturated']),
                                              np.random.choice(self.WEIGHTS['chain_complete'], p=self.WEIGHTS['saturation_total']),
                                              np.random.choice(self.WEIGHTS['chain_unsaturated'].values(), p=self.WEIGHTS['unsaturated']),
                                              np.random.choice(self.WEIGHTS['chain_complete'], p=self.WEIGHTS['saturation_total']),
                                              None, self.WEIGHTS['compartments'])
                elif head is None:
                    new_lipid = components.Lipid(head, np.random.choice(self.WEIGHTS['chain_unsaturated'].values(),
                                                                        p=self.WEIGHTS['unsaturated']),
                                                 np.random.choice(self.WEIGHTS['chain_complete'],
                                                                  p=self.WEIGHTS['saturation_total']), None, self.WEIGHTS['compartments'])
                    new_lipid.__class__ = components.TAG
                    new_lipid.sn3 = np.random.choice(self.WEIGHTS['chain_complete'], p=self.WEIGHTS['saturation_total'])
                elif head == 'ceramide':
                    new_lipid = components.Sphingolipid(head, None, self.WEIGHTS['compartments'])
                elif head == 'cdp':
                    new_lipid = components.Lipid('p', np.random.choice(self.WEIGHTS['chain_unsaturated'].values(),
                                                                       p=self.WEIGHTS['unsaturated']),
                                                 np.random.choice(self.WEIGHTS['chain_complete'],
                                                                  p=self.WEIGHTS['saturation_total']), None, self.WEIGHTS['compartments'])
                elif head == 'lyso':
                    new_lipid = components.Lipid('p', None, np.random.choice(self.WEIGHTS['chain_complete'],
                                                                             p=self.WEIGHTS['saturation_total']), None,
                                                 self.WEIGHTS['compartments'])
                elif head == 'dag':
                    new_lipid = components.Lipid(None, np.random.choice(self.WEIGHTS['chain_unsaturated'].values(),
                                                                        p=self.WEIGHTS['unsaturated']),
                                                 np.random.choice(self.WEIGHTS['chain_complete'],
                                                                  p=self.WEIGHTS['saturation_total']), None, self.WEIGHTS['compartments'])
                else:
                    new_lipid = components.Lipid(head, np.random.choice(self.WEIGHTS['chain_unsaturated'].values(),
                                                                        p=self.WEIGHTS['unsaturated']),
                                                 np.random.choice(self.WEIGHTS['chain_complete'],
                                                                  p=self.WEIGHTS['saturation_total']), None, self.WEIGHTS['compartments'])
                lipid_list.append(new_lipid)

        # generate 60 free FAs (acyl-coa)
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

    def __Km_calculation(self):
        """
        Calculate current Km values based on substrate concentrations.

        Return
        ------
        Km: dict
            current Km values of all reactions
        """
        # Calculation of Km, 2 substrates: Km for lipids = 5, for acyl-coa = 30
        Km = {'glycerol_3_p_synthesis': {'DHAP': 0.},
              'inositol_synthesis': {'glucose_6_p': 0.},
              'ceramide_synthesis': {'serine': 0.},
              'acetyl_coa_synthase': {'pyruvate': 0.},
              'acyl_synthase': {'acetyl_coa': 0.},
              'acyl_synthase_C16': 0.625,
              'acyl_synthase_C18': 0.002,
              'lyso_PA_synthase': {'acyl_coa': 30., 'DHAP': 0.},
              'PA_synthase': {'lyso-PA': 5., 'acyl_coa': 30.},
              'CDP_DG_synthase': {'PA': 5., 'CTP': 0.},
              'DAG_synthase': {'PA': 5.},
              'TAG_synthase': {'DAG': 5., 'acyl_coa': 30.},
              'TAG_lipase': {'lipid_droplets': 0.},
              'DAG_kinase': {'DAG': 5.},
              'PS_synthase': {'CDP_DG': 5., 'serine': 0.},
              'PI_synthase': {'CDP_DG': 5., 'inositol': 0.},
              'PE_synthase': {'PS': 5.},
              'PC_synthase': {'PE': 5., 'SAM': 0.},
              'CL_synthase': {'CDP_DG': 5., 'glycerol_3_p_mito': 0.},
              'ergosterol_synthase': {'acetyl_coa': 0.},
              'sterylester_synthase': {'ergosterol': 5., 'acyl_coa': 30.},
              'sphingolipid_synthase': {'PI': 5}}

        Km = {'glycerol_3_p_synthesis': {'DHAP': self.precursors_state['DHAP'] / self.INITIAL_PROBABILITIES['glycerol_3_p_synthesis'] - self.precursors_state['DHAP']}, \
              'inositol_synthesis': {'glucose_6_p': self.precursors_state['glucose_6_p'] / self.INITIAL_PROBABILITIES['inositol_synthesis'] - self.precursors_state['glucose_6_p']}, \
              'ceramide_synthesis': {'serine': self.precursors_state['serine'] / self.INITIAL_PROBABILITIES['ceramide_synthesis'] - self.precursors_state['serine']},\
              'acetyl_coa_synthase': {'pyruvate': self.precursors_state['pyruvate'] / self.INITIAL_PROBABILITIES['acetyl_coa_synthase'] - self.precursors_state['pyruvate']}, \
              'acyl_synthase': {'acetyl_coa': self.precursors_state['acetyl_coa'] / self.INITIAL_PROBABILITIES['acyl_synthase'] - self.precursors_state['acetyl_coa']},\
              'acyl_synthase_C16': 0.625,\
              'acyl_synthase_C18': 0.002,\
              'lyso_PA_synthase': {'acyl_coa': 30., 'DHAP': (self.precursors_state['DHAP'] / self.INITIAL_PROBABILITIES['lyso_PA_synthase']) * ((float(len(self.components_state['acyl_coa_saturated'])) + float(len(self.components_state['acyl_coa_unsaturated'])))\
                                     / (Km['lyso_PA_synthase']['acyl_coa'] + (float(len(self.components_state['acyl_coa_saturated'])) + float(len(self.components_state['acyl_coa_unsaturated']))))) - self.precursors_state['DHAP']}, \
              'PA_synthase': {'lyso-PA': 5., 'acyl_coa': 30.},\
              'CDP_DG_synthase': {'PA': 5., 'CTP': (self.precursors_state['CTP'] / self.INITIAL_PROBABILITIES['CDP_DG_synthase']) * (float(len(self.components_state['PA'])) / (Km['CDP_DG_synthase']['PA'] + float(len(self.components_state['PA'])))) - self.precursors_state['CTP']},\
              'DAG_synthase': {'PA': 5.}, \
              'TAG_synthase': {'DAG': 5., 'acyl_coa': 30.}, \
              'TAG_lipase': {'lipid_droplets': float(len(self.membranes_state['lipid_droplets'])) / self.INITIAL_PROBABILITIES['TAG_lipase'] - float(len(self.membranes_state['lipid_droplets']))}, \
              'DAG_kinase': {'DAG': 5.}, \
              'PS_synthase': {'CDP_DG': 5., 'serine': (self.precursors_state['serine'] / self.INITIAL_PROBABILITIES['PS_synthase']) * (float(len(self.components_state['CDP_DG'])) / (Km['PS_synthase']['CDP_DG'] + float(len(self.components_state['CDP_DG'])))) - self.precursors_state['serine']},\
              'PI_synthase': {'CDP_DG': 5., 'inositol': (self.precursors_state['inositol'] / self.INITIAL_PROBABILITIES['PI_synthase']) * (float(len(self.components_state['CDP_DG'])) / (Km['PI_synthase']['CDP_DG'] + float(len(self.components_state['CDP_DG'])))) - self.precursors_state['inositol']}, \
              'PE_synthase': {'PS': 5.}, \
              'PC_synthase': {'PE': 5., 'SAM': (self.precursors_state['SAM'] / self.INITIAL_PROBABILITIES['PC_synthase']) * (float(len(self.components_state['PE'])) / (Km['PC_synthase']['PE'] + float(len(self.components_state['PE'])))) - self.precursors_state['SAM']}, \
              'CL_synthase': {'CDP_DG': 5., 'glycerol_3_p_mito': (self.precursors_state['glycerol_3_p_mito'] / self.INITIAL_PROBABILITIES['CL_synthase']) * (float(len(self.components_state['CDP_DG'])) / (Km['CL_synthase']['CDP_DG'] + float(len(self.components_state['CDP_DG'])))) - self.precursors_state['glycerol_3_p_mito']},\
              'ergosterol_synthase': {'acetyl_coa': self.precursors_state['acetyl_coa'] / self.INITIAL_PROBABILITIES['ergosterol_synthase'] - self.precursors_state['acetyl_coa']}, \
              'sterylester_synthase': {'ergosterol': 5., 'acyl_coa': 30.}, \
              'sphingolipid_synthase': {'PI': 5, 'ceramide': (self.precursors_state['ceramide'] / self.INITIAL_PROBABILITIES['sphingolipid_synthase']) * (float(len(self.components_state['PI'])) / (Km['sphingolipid_synthase']['PI'] + float(len(self.components_state['PI'])))) - self.precursors_state['ceramide']}}
        return Km

    def cell_cycle(self, time):
        """
        Function to determine the cell cycle phase depending on the elapsed time.
        """
        if time <= 1800:
            return "G1"
        else:
            return "S_M"

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
                    self.comp_ratio_dict[comp][self.MEMBRANE_LIPID_NAMES[i]] = comp_ratio_list[i]
            x += 1

    def numbers(self):
        """
        Calculate numbers for plotting, write time courses.
        """
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