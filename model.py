"""
Simulation class for the lipid metabolism model.

@date: 03/06/2015
@author: Vera Schuetzhold - vera.schue@gmail.com
"""
import matplotlib.pyplot as mat
import numpy as np
import components


class Model:
    """
    Model class
    ===========
    At the beginning there are several lists defined which will contain the produced lipids.
    The model starts with a given amount of precursors (pyruvate, dhap, ctp, serine, glucose-6-p, SAM, glycerol-3-p).
    After several reactions and the transport function in the end there are membranes of different compartments with
    several different lipids.
    """
    def __init__(self):
        # initialise parameters and lists
        # model parameters
        self._init_parameters()
        # output lists
        self._init_lists()
        # initialise precursor production
        self._init_precursor_production()
        # initial state
        self._init_molecules()
        # lists for plotting
        self._init_plot_lists()

    def _init_parameters(self):
        """
        Initialise model parameters, all values were adjusted manually to the data of Uchida et al.
        (2011, PMID:21360734) and Zinser et al. (1991, PMID:2002005).
        """
        # VMAX of reactions
        # adjusted manually
        self._rates = {'glycerol_3_p_synthesis': 8, 
                       'inositol_synthesis':     5, 
                       'ceramide_synthesis':     2, 
                       'acetyl_coa_synthase':    650, 
                       'acyl_synthase':          450, 
                       'PA_synthase':            17,
                       'CDP_DG_synthase':        20, 
                       'TAG_synthase':           30,
                       'TAG_lipase':             23, 
                       'DAG_kinase':             40, 
                       'PS_synthase':            18, 
                       'PI_synthase':            6, 
                       'PE_synthase':            12, 
                       'PC_synthase':            5, 
                       'CL_synthase':            2, 
                       'Ergosterol_synthase':    25, 
                       'Sterylester_synthase':   25, 
                       'Sphingolipid_synthase':  2}

        # probabilities of reaction to take place
        # adjusted manually
        self._probability = {'glycerol_3_p_synthesis': 0.5, 
                             'inositol_synthesis':     0.5, 
                             'ceramide_synthesis':     0.5, 
                             'acetyl_coa_synthase':    0.8, 
                             'acyl_synthase':          0.5, 
                             'acyl_synthase_C16':      0.375, 
                             'acyl_synthase_C18':      0.998, 
                             'lyso_PA_synthase':       0.45, 
                             'PA_synthase':            0.2, 
                             'CDP_DG_synthase':        0.8,
                             'DAG_synthase':           0.01, 
                             'TAG_synthase':           0.2, 
                             'TAG_lipase':             0.8, 
                             'DAG_kinase':             0.1, 
                             'PS_synthase':            0.5, 
                             'PI_synthase':            0.5, 
                             'PE_synthase':            0.5, 
                             'PC_synthase':            0.5, 
                             'CL_synthase':            0.05, 
                             'Ergosterol_synthase':    0.6, 
                             'Sterylester_synthase':   0.4, 
                             'Sphingolipid_synthase':  0.2}

        # probabilities in CC phase G1
        # adjusted manually
        self._probability_G1 = {'DAG_synthase': 0.3, 
                                'TAG_synthase':  0.2,
                                'TAG_lipase':    0.05,
                                'DAG_kinase':    0.03}
        # probabilities in CC phases S-M
        # adjusted manually
        self._probability_S_M = {'DAG_synthase': 0.01, 
                                 'TAG_synthase':  0.2,
                                 'TAG_lipase':    0.6,
                                 'DAG_kinase':    0.1,
                                 'Sterylester_synthase': 0.2}

        # thresholds calculated from probabilities: (1-probability)
        self._manual_threshold = {reaction: 1-prob for reaction, prob in self._probability.iteritems()}
        # compartment weights (probabilities)
        # adjusted manually
        self._compartment_weights = [0.67, 0.01, 0.155, 0.101, 0.007, 0.007, 0.03, 0.015]
        # FA probabilities
        # adjusted manually
        # saturated, unsaturated
        self._weights_fa = [0.4, 0.6]
        # C16, C18
        self._unsaturated_weights = [0.375, 0.625]
        # C16:0, C18:0, C16:1, C18:1
        self._saturation_weights_total = [0.2, 0.2, 0.167, 0.433]

    def _init_plot_lists(self):
        # component list
        self.components_list = [self.acyl_coa_list, self.PA_list, self.CDP_DG_list, self.TAG_list, self.PS_list, 
                                self.PI_list, self.PE_list, self.PC_list, self.CL_list, self.Ergosterol_list, 
                                self.Sterylester_list, self.DAG_list, self.Sphingolipid_list]

    def _init_lists(self):
        # lists of all precursors
        self.precursors = {'pyruvate': [], 'acetyl_coa': [], 'glycerol-3-p': [], 'DHAP': [], 'serine': [],
                           'glucose_6_p': [], 'SAM': [], 'SAH': [], 'glycerol_3_p_mito': [], 'ceramide': [],
                           'GDP-mannose': [], 'NAD': [], 'NADH': [], 'NADP': [], 'NADPH': [], 'O2': [], 'H2O': [],
                           'CO2': [], 'Pi': [], 'CTP': [], 'CMP': [], 'inositol': [], 'ATP': [], 'ADP': []}

        # empty lists for the produced fatty acids and lipids
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
        
        # lipids to transport
        self.transport_list = [self.PS_list, self.PI_list, self.PC_list, self.PE_list, self.CL_list, self.PA_list, 
                               self.Ergosterol_list, self.Sterylester_list, self.TAG_list, self.Sphingolipid_list]

        # lists to collect the transported lipids
        # compartment names
        self.compartment = ['plasma_membrane', 'secretory_vesicles', 'vacuoles', 'nucleus', 'peroxisomes',
                            'light_microsomes', 'inner_mit_membrane', 'outer_mit_membrane', 'lipid_droplets']
        # compartment lists
        self.plasma_membrane = []
        self.secretory_vesicles = []
        self.vacuoles = []
        self.nucleus = []
        self.peroxisomes = []
        self.light_microsomes = []
        self.inner_mit_membrane = []
        self.outer_mit_membrane = []
        self.lipid_droplets = []
        # list of lists
        self.compartment_lists = [self.plasma_membrane, self.secretory_vesicles, self.vacuoles, self.nucleus,
                                  self.peroxisomes, self.light_microsomes, self.inner_mit_membrane,
                                  self.outer_mit_membrane, self.lipid_droplets]

        # relatives list for calculation of the relative parts of each lipid in a membrane for the
        # compartment_relatives_dict
        self.relatives_list = []
        # collecting the products of every time step. Lipids that are produced in the start function
        # don't need the 0 for plotting
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
        # list of lists
        self.number_lipids_list = [self.number_acyl_coa, self.number_pa, self.number_cdp_dg, self.number_tag,
                                   self.number_PS, self.number_PI, self.number_PE, self.number_PC, self.number_CL,
                                   self.number_Ergosterol, self.number_Sterylester, self.number_DAG,
                                   self.number_Sphingolipid]
        # counting the lipids in each membrane after every timestep
        self.number_plasma_membrane = []
        self.number_secretory_vesicles = []
        self.number_vacuoles = []
        self.number_nucleus = []
        self.number_peroxisomes = []
        self.number_light_microsomes = []
        self.number_inner_mit_membrane = []
        self.number_outer_mit_membrane = []
        self.number_lipid_droplets = []
        # list of lists
        self.number_membranes_list = [self.number_plasma_membrane, self.number_secretory_vesicles, self.number_vacuoles,
                                      self.number_nucleus, self.number_peroxisomes, self.number_light_microsomes,
                                      self.number_inner_mit_membrane, self.number_outer_mit_membrane,
                                      self.number_lipid_droplets]

        # possible fatty acids and weights to reach the biological proportions
        self.chainlength_saturated = {16: 'C16:0', 18: 'C18:0'}
        self.chainlength_unsaturated = {16: 'C16:1', 18: 'C18:1'}
        self.chainlength_saturated_unsaturated = ['C16:0', 'C18:0', 'C16:1', 'C18:1']
        # names of membrane lipids
        self.membrane_lipids = ['PS', 'PI', 'PC', 'PE', 'CL', 'PA', 'ES', 'SE', 'TAG', 'SL']

        self.compartment_relatives_dict = \
            {comp: dict(zip(self.membrane_lipids, [0.0 for i in range(10)])) for comp in self.compartment}

        # dict for Km values
        self.Km = {}

        # create list of functions
        self.function_list = [self.glycerol_3_p_synthesis,
                              self.inositol_synthesis,
                              self.ceramide_synthesis,
                              self.acetyl_coa_synthase,
                              self.acyl_synthase,
                              self.PA_synthase,
                              self.CDP_DG_synthase,
                              self.TAG_synthase,
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

    def _init_precursor_production(self):

        # number of small molecules that is produced from anywhere in the cell and will be added every 10 seconds
        # G1 phase
        self.precursors_production_G1 = {'pyruvate': 1500., 'acetyl_coa': 0, 'glycerol-3-p': 5., 'DHAP': 30.,
                                         'serine': 20., 'glucose_6_p': 8., 'SAM': 45., 'SAH': 0.,
                                         'glycerol_3_p_mito': 5., 'ceramide': 0, 'GDP-mannose': 10, 'NAD': 0,
                                         'NADH': 0, 'NADP': 0, 'NADPH': 0, 'O2': 0, 'H2O': 0, 'CO2': 0, 'Pi': 0,
                                         'CTP': 20, 'CMP': 0, 'inositol': 0, 'ATP': 0, 'ADP': 0}
        # S-M phase
        self.precursors_production_S_M = {'pyruvate': 1700., 'acetyl_coa': 0, 'glycerol-3-p': 5., 'DHAP': 35.,
                                          'serine': 30., 'glucose_6_p': 12., 'SAM': 55., 'SAH': 0.,
                                          'glycerol_3_p_mito': 5., 'ceramide': 0, 'GDP-mannose': 10, 'NAD': 0,
                                          'NADH': 0, 'NADP': 0, 'NADPH': 0, 'O2': 0, 'H2O': 0, 'CO2': 0, 'Pi': 0,
                                          'CTP': 45, 'CMP': 0, 'inositol': 0, 'ATP': 0, 'ADP': 0}

    def _init_molecules(self):

        # initialise number of small molecules in the cell
        self.precursors_dict = {'pyruvate': 2200., 'acetyl_coa': 1000, 'glycerol-3-p': 1000., 'DHAP': 1000.,
                                'serine': 250., 'glucose_6_p': 1000., 'SAM': 331., 'SAH': 0., 'glycerol_3_p_mito': 50.,
                                'ceramide': 100, 'GDP-mannose': 0, 'NAD': 0, 'NADH': 0, 'NADP': 0, 'NADPH': 0, 'O2': 0,
                                'H2O': 0, 'CO2': 0, 'Pi': 0, 'CTP': 1000, 'CMP': 0, 'inositol': 350, 'ATP': 0, 'ADP': 0}

        # initial lipid molecules, based on Zinser et al. (1991, PMID:2002005)
        # 'other' were interpreted as sphingolipids
        self.plasma_membrane_comp_start = [0.17320, 0.09124, 0.08660, 0.10464, 0.00103, 0.02010, 0.48454, 0.0, 0.0, 0.03557]
        self.secretory_vesicles_comp_start = [0.08205, 0.11745, 0.20824, 0.13573, 0.01239, 0.01525, 0.42900, 0.0, 0.0, 0.05029]
        self.vacuoles_comp_start = [0.04817, 0.16604, 0.40517, 0.17537, 0.02442, 0.02866, 0.15200, 0.0, 0.0, 0.06525]
        self.nucleus_comp_start = [0.04038, 0.09650, 0.27645, 0.16848, 0.01049, 0.01781, 0.390, 0.0, 0.0, 0.02622]
        self.peroxisomes_comp_start = [0.03235, 0.11360, 0.34656, 0.16465, 0.05033, 0.01150, 0.281, 0.0, 0.0, 0.0]
        self.light_microsomes_comp_start = [0.05304, 0.06019, 0.40796, 0.26583, 0.00381, 0.00222, 0.206, 0.0, 0.0, 0.00397]
        self.inner_mit_membrane_comp_start = [0.02880, 0.12273, 0.29107, 0.18192, 0.12204, 0.01137, 0.242, 0.0, 0.0, 0.0]
        self.outer_mit_membrane_comp_start = [0.01189, 0.10108, 0.45190, 0.32307, 0.05847, 0.04360, 0.009, 0.0, 0.0, 0.0]
        self.lipid_droplets_comp_start = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.0]
        # list of lists
        self.membrane_compositions_start = [self.plasma_membrane_comp_start, self.secretory_vesicles_comp_start,
                                            self.vacuoles_comp_start, self.nucleus_comp_start,
                                            self.peroxisomes_comp_start, self.light_microsomes_comp_start,
                                            self.inner_mit_membrane_comp_start, self.outer_mit_membrane_comp_start, 
                                            self.lipid_droplets_comp_start]

    def run(self, timesteps=7200):
        """
        Start the simulation.

        Parameters
        ----------
        timesteps: int
            number of simulation steps in sec

        Returns
        -------
        self.saturation_composition_total: 
        self.number_membranes_list:
        self.compartment_relatives_dict: 
        """
        # determine the timesteps, self.t for plotting, self.time for cell cycle
        self.time = 0
        self.t = [i for i in range(timesteps)]
        # function that produces the lipids and membranes that are existing at the beginning of the cell cycle
        self.start()
        self.Km_calculation()
        for t in range(timesteps):
            # counting the seconds for knowing the cell cycle phase
            self.time += 1

            # cell cycle phase dependent precursor production
            if self.cell_cycle() == 'G1':
                self.precursors_production = self.precursors_production_G1
            else:
                self.precursors_production = self.precursors_production_S_M

            # produce precursor every 10 time steps 
            if self.time % 10 == 0:
                for key in self.precursors_dict:
                    self.precursors_dict[key] += self.precursors_production[key]

            # calculate new thresholds
            self.thresholds = self._calculate_threshold()
            # manual_threshold
            self.probabilities = self.thresholds  

            # change rates dependent on cell cycle phase
            if self.cell_cycle() == 'G1':
                self.probabilities['DAG_synthase'] = 1 - self._probability_G1['DAG_synthase']
                self.probabilities['TAG_synthase'] = 1 - self._probability_G1['TAG_synthase']
                self.probabilities['TAG_lipase'] = 1 - self._probability_G1['TAG_lipase']
                self.probabilities['DAG_kinase'] = 1 - self._probability_G1['DAG_kinase']
            else:
                self.probabilities['DAG_synthase'] = 1 - self._probability_S_M['DAG_synthase']
                self.probabilities['TAG_synthase'] = 1 - self._probability_S_M['TAG_synthase']
                self.probabilities['TAG_lipase'] = 1 - self._probability_S_M['TAG_lipase']
                self.probabilities['DAG_kinase'] = 1 - self._probability_S_M['DAG_kinase']
                self.probabilities['Sterylester_synthase'] = 1 - self._probability_S_M['Sterylester_synthase']

            # all reactions that take place during one second in a random order
            for i in range(len(self.function_list)):
                func = np.random.choice(self.function_list)
                func()
                self.function_list.remove(func)
            self.numbers()  # calculating the produced lipids after each time step
        self.membrane_compositions()  # calculation of the membrane compositions at the end of the run
        self.saturation_counter()  # calculating the percentages of each fatty acid that was used

        return self.saturation_composition_total, self.number_membranes_list, self.compartment_relatives_dict

    def _calculate_threshold(self):
        """
        Calculate current thresholds based on amount of available precursors.
        """
        threshold = {'glycerol_3_p_synthesis': 1 - (self.precursors_dict['DHAP'] /
                                                    (self.Km['glycerol_3_p_synthesis']['DHAP'] +
                                                     self.precursors_dict['DHAP'])),
                     'inositol_synthesis': 1 - (self.precursors_dict['glucose_6_p'] /
                                                (self.Km['inositol_synthesis']['glucose_6_p'] +
                                                 self.precursors_dict['glucose_6_p'])),
                     'ceramide_synthesis': 1 - (self.precursors_dict['serine'] /
                                                (self.Km['ceramide_synthesis']['serine'] +
                                                 self.precursors_dict['serine'])),
                     'acetyl_coa_synthase': 1 - (self.precursors_dict['pyruvate'] /
                                                 (self.Km['acetyl_coa_synthase']['pyruvate'] +
                                                  self.precursors_dict['pyruvate'])),
                     'acyl_synthase': 1 - (self.precursors_dict['acetyl_coa'] /
                                           (self.Km['acyl_synthase']['acetyl_coa'] +
                                            self.precursors_dict['acetyl_coa'])),
                     'acyl_synthase_C16': 0.625,
                     'acyl_synthase_C18': 0.002,
                     'lyso_PA_synthase': 1 - (((float(len(self.acyl_coa_list_saturated)) +
                                                float(len(self.acyl_coa_list_unsaturated))) /
                                               (self.Km['lyso_PA_synthase']['acyl_coa'] +
                                                (float(len(self.acyl_coa_list_saturated)) +
                                                 float(len(self.acyl_coa_list_unsaturated))))) *
                                              (self.precursors_dict['DHAP'] /
                                               (self.Km['lyso_PA_synthase']['DHAP'] + self.precursors_dict['DHAP']))),
                     'PA_synthase': 1 - ((float(len(self.lyso_pa_list)) / (self.Km['PA_synthase']['lyso-PA'] +
                                                                           float(len(self.lyso_pa_list)))) *
                                         ((float(len(self.acyl_coa_list_saturated)) +
                                           float(len(self.acyl_coa_list_unsaturated))) /
                                          (self.Km['PA_synthase']['acyl_coa'] +
                                           (float(len(self.acyl_coa_list_saturated)) +
                                            float(len(self.acyl_coa_list_unsaturated)))))),
                     'CDP_DG_synthase': 1 - ((float(len(self.PA_list)) / (self.Km['CDP_DG_synthase']['PA'] +
                                                                          float(len(self.PA_list)))) *
                                             (self.precursors_dict['CTP'] / (self.Km['CDP_DG_synthase']['CTP'] +
                                                                             self.precursors_dict['CTP']))),
                     'DAG_synthase': 1 - (float(len(self.PA_list)) / (self.Km['DAG_synthase']['PA'] +
                                                                      float(len(self.PA_list)))),
                     'TAG_synthase': 1 - (((float(len(self.acyl_coa_list_saturated)) +
                                            float(len(self.acyl_coa_list_unsaturated))) /
                                           (self.Km['TAG_synthase']['acyl_coa'] +
                                            (float(len(self.acyl_coa_list_saturated)) +
                                             float(len(self.acyl_coa_list_unsaturated))))) *
                                          (float(len(self.DAG_list)) /
                                           (self.Km['TAG_synthase']['DAG'] + float(len(self.DAG_list))))),
                     'TAG_lipase': 1 - (float(len(self.lipid_droplets)) / (self.Km['TAG_lipase']['lipid_droplets'] +
                                                                           float(len(self.lipid_droplets)))),
                     'DAG_kinase': 1 - (float(len(self.DAG_list)) / (self.Km['DAG_kinase']['DAG'] +
                                                                     float(len(self.DAG_list)))),
                     'PS_synthase': 1 - ((float(len(self.CDP_DG_list)) / (self.Km['PS_synthase']['CDP_DG'] +
                                                                          float(len(self.CDP_DG_list)))) *
                                         (self.precursors_dict['serine'] / (self.Km['PS_synthase']['serine'] +
                                                                            self.precursors_dict['serine']))),
                     'PI_synthase': 1 - ((float(len(self.CDP_DG_list)) / (self.Km['PI_synthase']['CDP_DG'] +
                                                                          float(len(self.CDP_DG_list)))) *
                                         (self.precursors_dict['inositol'] / (self.Km['PI_synthase']['inositol'] +
                                                                              self.precursors_dict['inositol']))),
                     'PE_synthase': 1 - (float(len(self.PS_list)) / (self.Km['PE_synthase']['PS'] +
                                                                     float(len(self.PS_list)))),
                     'PC_synthase': 1 - ((float(len(self.PE_list)) / (self.Km['PC_synthase']['PE'] +
                                                                      float(len(self.PE_list)))) *
                                         (self.precursors_dict['SAM'] / (self.Km['PC_synthase']['SAM'] +
                                                                         self.precursors_dict['SAM']))),
                     'CL_synthase': 1 - ((float(len(self.CDP_DG_list)) / (self.Km['CL_synthase']['CDP_DG'] +
                                                                          float(len(self.CDP_DG_list)))) *
                                         (self.precursors_dict['glycerol_3_p_mito'] /
                                         (self.Km['CL_synthase']['glycerol_3_p_mito'] +
                                          self.precursors_dict['glycerol_3_p_mito']))),
                     'Ergosterol_synthase': 1 - (self.precursors_dict['acetyl_coa'] /
                                                 (self.Km['Ergosterol_synthase']['acetyl_coa'] +
                                                  self.precursors_dict['acetyl_coa'])),
                     'Sterylester_synthase': 1 - ((float(len(self.Ergosterol_list)) /
                                                  (self.Km['Sterylester_synthase']['ergosterol'] +
                                                   float(len(self.Ergosterol_list)))) *
                                                  ((float(len(self.acyl_coa_list_saturated)) +
                                                   float(len(self.acyl_coa_list_unsaturated))) /
                                                  (self.Km['Sterylester_synthase']['acyl_coa'] +
                                                   (float(len(self.acyl_coa_list_saturated)) +
                                                    float(len(self.acyl_coa_list_unsaturated)))))),
                     'Sphingolipid_synthase': 1 - ((float(len(self.PI_list)) /
                                                   (self.Km['Sphingolipid_synthase']['PI'] +
                                                    float(len(self.PI_list)))) *
                                                   (self.precursors_dict['ceramide'] /
                                                   (self.Km['Sphingolipid_synthase']['ceramide'] +
                                                    self.precursors_dict['ceramide'])))}
        return threshold

    def start(self):
        """
        Function that produces the starting lipids. Membrane compositions are given in self.membrane_compositions_start_relatives.
        Number of starting lipids for each membrane are given in self.start_lipids.
        """
        self.membrane_compositions_start_relatives = []

        # numbers of the comp_start lists should yield 1 when summed
        for membrane_comp_start in self.membrane_compositions_start:
            membrane_comp_start_relative = [z / sum(membrane_comp_start) for z in membrane_comp_start]
            self.membrane_compositions_start_relatives.append(membrane_comp_start_relative)

        self.compositions_start = dict(zip(self.compartment, self.membrane_compositions_start_relatives))

        x = 0
        # number of lipids that are produced in the start function for every membrane
        self.start_lipids = [32950, 500, 2500, 6000, 500, 500, 5000, 2500, 1000]
        self.membrane_start = dict(zip(self.compartment, self.start_lipids))
        for membrane in self.compartment_lists:
            # producing the lipids for a membrane, probability for a certain lipid from the composition in Zinser
            for i in range(self.membrane_start[self.compartment[x]]):
                self.head_groups_start = ['serine', 'inositol', 'choline', 'ethanolamine', 'neutral', 'p', 'sterol',
                                          'sterylester', None, 'ceramide']
                weights_start = self.compositions_start[self.compartment[x]]
                head = np.random.choice(self.head_groups_start, p=weights_start)
                if head == 'sterol':
                    new_lipid = components.Sterol(head, self.compartment[x], self._compartment_weights)
                elif head == 'sterylester':
                    new_lipid = components.Sterylester(head, np.random.choice(self.chainlength_unsaturated.values(),
                                                                              p=[0.67, 0.33]), self.compartment[x],
                                                       self._compartment_weights)
                elif head == 'ceramide':
                    new_lipid = components.Sphingolipid(head, self.compartment[x], self._compartment_weights)
                elif head == 'neutral':
                    new_lipid = components.CL(head, np.random.choice(self.chainlength_unsaturated.values(),
                                                                     p=self._unsaturated_weights),
                                              np.random.choice(self.chainlength_saturated_unsaturated,
                                                               p=self._saturation_weights_total),
                                    np.random.choice(self.chainlength_unsaturated.values(),
                                                     p=self._unsaturated_weights),
                                              np.random.choice(self.chainlength_saturated_unsaturated,
                                                               p=self._saturation_weights_total), self.compartment[x],
                                              self._compartment_weights)
                elif head == 'serine' or head == 'inositol' or head == 'choline' or head == 'ethanolamine' or head == 'p':
                    new_lipid = components.Lipid(head, np.random.choice(self.chainlength_unsaturated.values(),
                                                                        p=self._unsaturated_weights),
                                                 np.random.choice(self.chainlength_saturated_unsaturated,
                                                                  p=self._saturation_weights_total),
                                                 self.compartment[x], self._compartment_weights)
                else:
                    new_lipid = components.Lipid(head, np.random.choice(self.chainlength_unsaturated.values(),
                                                                        p=self._unsaturated_weights),
                                                 np.random.choice(self.chainlength_saturated_unsaturated,
                                                                  p=self._saturation_weights_total), self.compartment[x],
                                                 self._compartment_weights)
                    new_lipid.__class__ = components.TAG
                    new_lipid.sn3 = np.random.choice(self.chainlength_saturated_unsaturated,
                                                     p=self._saturation_weights_total)
                membrane.append(new_lipid)
            x += 1

        self.lipid_lists_start = [self.PS_list, self.PI_list, self.PC_list, self.PE_list, self.CL_list, self.lyso_pa_list, self.PA_list,
                                  self.CDP_DG_list, self.Ergosterol_list, self.Sterylester_list, self.DAG_list, self.TAG_list,
                                  self.Sphingolipid_list]
        z = 0
        for lipid_list in self.lipid_lists_start:
            for i in range(20):
                self.head_groups_start_lipids = ['serine', 'inositol', 'choline', 'ethanolamine', 'neutral', 'lyso', 'p', 'cdp',
                                                 'sterol', 'sterylester', 'dag', None, 'ceramide']
                head = self.head_groups_start_lipids[z]
                if head == 'sterol':
                    new_lipid = components.Sterol(head, None, self._compartment_weights)
                elif head == 'sterylester':
                    new_lipid = components.Sterylester(head, np.random.choice(self.chainlength_unsaturated.values(),
                                                                              p=[0.67, 0.33]), None, self._compartment_weights)
                elif head == 'neutral':
                    new_lipid = components.CL(head, np.random.choice(self.chainlength_unsaturated.values(), p=self._unsaturated_weights),
                                              np.random.choice(self.chainlength_saturated_unsaturated, p=self._saturation_weights_total),
                                              np.random.choice(self.chainlength_unsaturated.values(), p=self._unsaturated_weights),
                                              np.random.choice(self.chainlength_saturated_unsaturated, p=self._saturation_weights_total),
                                              None, self._compartment_weights)
                elif head is None:
                    new_lipid = components.Lipid(head, np.random.choice(self.chainlength_unsaturated.values(),
                                                                        p=self._unsaturated_weights),
                                                 np.random.choice(self.chainlength_saturated_unsaturated,
                                                                  p=self._saturation_weights_total), None, self._compartment_weights)
                    new_lipid.__class__ = components.TAG
                    new_lipid.sn3 = np.random.choice(self.chainlength_saturated_unsaturated, p=self._saturation_weights_total)
                elif head == 'ceramide':
                    new_lipid = components.Sphingolipid(head, None, self._compartment_weights)
                elif head == 'cdp':
                    new_lipid = components.Lipid('p', np.random.choice(self.chainlength_unsaturated.values(),
                                                                       p=self._unsaturated_weights),
                                                 np.random.choice(self.chainlength_saturated_unsaturated,
                                                                  p=self._saturation_weights_total), None, self._compartment_weights)
                elif head == 'lyso':
                    new_lipid = components.Lipid('p', None, np.random.choice(self.chainlength_saturated_unsaturated,
                                                                             p=self._saturation_weights_total), None,
                                                 self._compartment_weights)
                elif head == 'dag':
                    new_lipid = components.Lipid(None, np.random.choice(self.chainlength_unsaturated.values(),
                                                                        p=self._unsaturated_weights),
                                                 np.random.choice(self.chainlength_saturated_unsaturated,
                                                                  p=self._saturation_weights_total), None, self._compartment_weights)
                else:
                    new_lipid = components.Lipid(head, np.random.choice(self.chainlength_unsaturated.values(),
                                                                        p=self._unsaturated_weights),
                                                 np.random.choice(self.chainlength_saturated_unsaturated,
                                                                  p=self._saturation_weights_total), None, self._compartment_weights)
                lipid_list.append(new_lipid)
            z += 1

        choice_list_acyl_start = [0, 1]
        choice_weights_acyl_start = [0.13, 0.87]
        choice_c_acyl_start = [16, 18]
        for i in range(60):
            new_acyl = components.FattyAcid(np.random.choice(choice_c_acyl_start), np.random.choice(choice_list_acyl_start,
                                                                                                    p=choice_weights_acyl_start))
            if new_acyl.saturation == 0:
                self.acyl_coa_list_saturated.append(new_acyl)
            else:
                self.acyl_coa_list_unsaturated.append(new_acyl)

    def Km_calculation(self):
        # Calculation of Km, 2 substrates: Km for lipids = 5, for acyl-coa = 30
        self.pre_Km = { 'lyso_PA_synthase': {'acyl_coa': 30.},
                        'CDP_DG_synthase': {'PA': 5.},
                        'PS_synthase': {'CDP_DG': 5.},
                        'PI_synthase': {'CDP_DG': 5.},
                        'PC_synthase': {'PE': 5.}, 
                        'CL_synthase': {'CDP_DG': 5.}, 
                        'Sphingolipid_synthase': {'PI': 5.}}

        self.Km = {'glycerol_3_p_synthesis': {'DHAP': self.precursors_dict['DHAP'] / self._probability['glycerol_3_p_synthesis'] -
                                              self.precursors_dict['DHAP']},
                        'inositol_synthesis': {'glucose_6_p': self.precursors_dict['glucose_6_p'] / self._probability['inositol_synthesis'] -
                                               self.precursors_dict['glucose_6_p']},
                        'ceramide_synthesis': {'serine': self.precursors_dict['serine'] / self._probability['ceramide_synthesis'] -
                                               self.precursors_dict['serine']},
                        'acetyl_coa_synthase': {'pyruvate': self.precursors_dict['pyruvate'] / self._probability['acetyl_coa_synthase'] -
                                                self.precursors_dict['pyruvate']},
                        'acyl_synthase': {'acetyl_coa': self.precursors_dict['acetyl_coa'] / self._probability['acyl_synthase'] -
                                          self.precursors_dict['acetyl_coa']},
                        'acyl_synthase_C16': 0.625,
                        'acyl_synthase_C18': 0.002,
                        'lyso_PA_synthase': {'acyl_coa': 30., 
                                             'DHAP': (self.precursors_dict['DHAP'] / self._probability['lyso_PA_synthase']) *
                                                     ((float(len(self.acyl_coa_list_saturated)) +
                                                       float(len(self.acyl_coa_list_unsaturated))) /
                                                      (self.pre_Km['lyso_PA_synthase']['acyl_coa'] +
                                                       (float(len(self.acyl_coa_list_saturated)) +
                                                        float(len(self.acyl_coa_list_unsaturated))))) - self.precursors_dict['DHAP']},
                        'PA_synthase': {'lyso-PA': 5., 
                                        'acyl_coa': 30.},
                        'CDP_DG_synthase': {'PA': 5., 
                                            'CTP': (self.precursors_dict['CTP'] / self._probability['CDP_DG_synthase']) *
                                                   (float(len(self.PA_list)) / (self.pre_Km['CDP_DG_synthase']['PA'] +
                                                                                float(len(self.PA_list)))) - self.precursors_dict['CTP']},
                        'DAG_synthase': {'PA': 5.}, 
                        'TAG_synthase': {'DAG': 5., 
                                         'acyl_coa': 30.}, 
                        'TAG_lipase': {'lipid_droplets': float(len(self.lipid_droplets)) / self._probability['TAG_lipase'] -
                                       float(len(self.lipid_droplets))},
                        'DAG_kinase': {'DAG': 5.}, 
                        'PS_synthase': {'CDP_DG': 5., 
                                        'serine': (self.precursors_dict['serine'] / self._probability['PS_synthase']) *
                                                  (float(len(self.CDP_DG_list)) / (self.pre_Km['PS_synthase']['CDP_DG'] +
                                                                                   float(len(self.CDP_DG_list)))) -
                                        self.precursors_dict['serine']},
                        'PI_synthase': {'CDP_DG': 5., 
                                        'inositol': (self.precursors_dict['inositol'] / self._probability['PI_synthase']) *
                                        (float(len(self.CDP_DG_list)) / (self.pre_Km['PI_synthase']['CDP_DG'] +
                                                                         float(len(self.CDP_DG_list)))) - self.precursors_dict['inositol']},
                        'PE_synthase': {'PS': 5.}, 
                        'PC_synthase': {'PE': 5., 
                                        'SAM': (self.precursors_dict['SAM'] / self._probability['PC_synthase']) *
                                               (float(len(self.PE_list)) / (self.pre_Km['PC_synthase']['PE'] + float(len(self.PE_list)))) -
                                        self.precursors_dict['SAM']},
                        'CL_synthase': {'CDP_DG': 5., 
                                        'glycerol_3_p_mito': (self.precursors_dict['glycerol_3_p_mito'] /
                                                              self._probability['CL_synthase']) * (float(len(self.CDP_DG_list)) /
                                                                                                   (self.pre_Km['CL_synthase']['CDP_DG'] +
                                                                                                    float(len(self.CDP_DG_list)))) -
                                        self.precursors_dict['glycerol_3_p_mito']},
                        'Ergosterol_synthase': {'acetyl_coa': self.precursors_dict['acetyl_coa'] / self._probability['Ergosterol_synthase'] -
                                                self.precursors_dict['acetyl_coa']},
                        'Sterylester_synthase': {'ergosterol': 5., 
                                                 'acyl_coa': 30.}, 
                        'Sphingolipid_synthase': {'PI': 5, 
                                                  'ceramide': (self.precursors_dict['ceramide'] / self._probability['Sphingolipid_synthase']) *
                                                  (float(len(self.PI_list)) / (self.pre_Km['Sphingolipid_synthase']['PI'] +
                                                                               float(len(self.PI_list)))) - self.precursors_dict['ceramide']}}

    def cell_cycle(self):
        """
        Function to determine the cell cycle phases depending on the elapsed time.
        """
        if self.time <= 1800:
            return 'G1'
        else:
            return 'S/G2/M'

    def glycerol_3_p_synthesis(self):
        """
        Synthesis of glycerol-3-p out of DHAP.
        """
        for i in range(self._rates['glycerol_3_p_synthesis']):
            x = np.random.random()
            if x >= self.probabilities['glycerol_3_p_synthesis']:
                if self.precursors_dict['DHAP'] > 1:
                    self.precursors_dict['glycerol-3-p'] += 1
                    self.precursors_dict['DHAP'] -= 1
                    self.precursors_dict['NADH'] += 1
                    self.precursors_dict['NAD'] -= 1

    def inositol_synthesis(self):
        """
        Synthesis of inositol from glucose-6-p by the Myo-inositol-1-p synthase and the Myo-inositol 1 phosphatase.
        """
        for i in range(self._rates['inositol_synthesis']):
            x = np.random.random()
            if x >= self.probabilities['inositol_synthesis']:
                if self.precursors_dict['glucose_6_p'] > 1:
                    self.precursors_dict['inositol'] += 1
                    self.precursors_dict['glucose_6_p'] -= 1
                    self.precursors_dict['H2O'] -= 1
                    self.precursors_dict['Pi'] += 1

    def ceramide_synthesis(self):
        """
        Synthesis of ceramide out of serine and a C16:0 fatty acid
        """
        for i in range(self._rates['ceramide_synthesis']):
            x = np.random.random()
            if x >= self.probabilities['ceramide_synthesis']:
                if len(self.acyl_coa_list_C26) > 1 and self.precursors_dict['serine'] > 1 and len(self.acyl_coa_list_saturated) > 1 and \
                        any(fa.C == 16 for fa in self.acyl_coa_list_saturated):
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
        """
        Synthesis of Acetyl-CoA: pyruvate dehydrogenase drives the reaction pyruvate to Acetyl-CoA, CO2 is released
        """
        for i in range(self._rates['acetyl_coa_synthase']):
            x = np.random.random()
            if x >= self.probabilities['acetyl_coa_synthase']:
                if self.precursors_dict['pyruvate'] > 1:  # transformation from pyruvate to acetyl_coa
                    self.precursors_dict['acetyl_coa'] += 1
                    self.precursors_dict['pyruvate'] -= 1
                    self.precursors_dict['NADH'] += 1
                    self.precursors_dict['NAD'] -= 1
                    self.precursors_dict['CO2'] += 1

    def acyl_synthase(self):
        """
        Simplified synthesis of acyl_coa: several Acetyl-CoA are building a fatty acid (C16:0, C16:1, C18:0 or C18:1)
        The intermediate Malonyl-CoA is leaved out.
        """
        choice_list = [0, 1]
        choice_weights = [0.12, 0.88]
        for i in range(self._rates['acyl_synthase']):
            # 5 reactions in 1 time step but only with a probability of 90%
            x = np.random.random()
            # control if at least 2 Acetyl-CoA are available
            if self.precursors_dict['acetyl_coa'] > 2:
                # starting the first reaction
                if len(self.acyl_coa_list) == 0:
                    new_acyl = components.FattyAcid(2, np.random.choice(choice_list, p=choice_weights))
                    self.acyl_coa_list.append(new_acyl)
                    self.acyl_coa_list[-1].C += 2
                    self.precursors_dict['acetyl_coa'] -= 2
                # stop the reaction cycle and starting a new one
                elif self.acyl_coa_list[-1].C == 16 and x >= self.probabilities['acyl_synthase_C16']:
                    new_acyl = components.FattyAcid(2, np.random.choice(choice_list, p=choice_weights))
                    self.acyl_coa_list.append(new_acyl)
                    self.acyl_coa_list[-1].C += 2
                    self.precursors_dict['acetyl_coa'] -= 2
                    self.precursors_dict['NADPH'] -= 14
                    self.precursors_dict['NADP'] += 14
                    self.precursors_dict['H2O'] += 7
                    # CO2 production is not mentioned here as only acetyl-CoA is used and not malonyl-CoA, so we need all C-atoms
                    # we give in the reaction
                # stop the reaction cycle and starting a new one
                elif self.acyl_coa_list[-1].C == 18 and x >= self.probabilities['acyl_synthase_C18']:
                    new_acyl = components.FattyAcid(2, np.random.choice(choice_list, p=choice_weights))
                    self.acyl_coa_list.append(new_acyl)
                    self.acyl_coa_list[-1].C += 2
                    self.precursors_dict['acetyl_coa'] -= 2
                    self.precursors_dict['NADPH'] -= 16
                    self.precursors_dict['NADP'] += 16
                    self.precursors_dict['H2O'] += 8

                elif self.acyl_coa_list[-1].C == 26:
                    self.acyl_coa_list[-1].saturation = 0
                    new_acyl = components.FattyAcid(2, np.random.choice(choice_list, p=choice_weights))
                    self.acyl_coa_list.append(new_acyl)
                    self.acyl_coa_list[-1].C += 2
                    self.precursors_dict['acetyl_coa'] -= 2
                    self.precursors_dict['NADPH'] -= 24
                    self.precursors_dict['NADP'] += 24
                    self.precursors_dict['H2O'] += 12
                # adding an Acetyl_CoA to the growing ffa
                else:
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

    def PA_synthase(self):
        """
        Synthesis of PA in two reaction steps.
        """
        for i in range(self._rates['PA_synthase']):
            self.lyso_PA_synthase()
            self.PA_synthase()

    def lyso_PA_synthase(self):
        """
        Production of Lyso-PA by adding one acyl-coa to DHAP (sn1: always unsaturated) --> DHAP acyltransferase/acyl-DHAP reductase
        """
        choice_list = [0, 1]

        weights_pa = [self.precursors_dict['DHAP'] / (self.precursors_dict['DHAP'] + self.precursors_dict['glycerol-3-p']),
                      self.precursors_dict['glycerol-3-p'] / (self.precursors_dict['DHAP'] + self.precursors_dict['glycerol-3-p'])]
        x = np.random.random()
        # at least 1 ffa has to be unsaturated
        if x >= self.probabilities['lyso_PA_synthase'] and len(self.acyl_coa_list_saturated) > 1 and \
                len(self.acyl_coa_list_unsaturated) > 1 and (self.precursors_dict['DHAP'] > 1 and
                                                             self.precursors_dict['glycerol-3-p'] > 1):
            if np.random.choice(choice_list, p=self._weights_fa) == 0:
                sn1_chain = np.random.randint(0, (len(self.acyl_coa_list_saturated)-1))
                chainlength_sn1 = self.acyl_coa_list_saturated[sn1_chain].C
                lyso_pa = components.Lipid('p', None, self.chainlength_saturated[chainlength_sn1], None, self._compartment_weights)
                del self.acyl_coa_list_saturated[sn1_chain]
            else:
                sn1_chain = np.random.randint(0, (len(self.acyl_coa_list_unsaturated)-1))
                chainlength_sn1 = self.acyl_coa_list_unsaturated[sn1_chain].C
                lyso_pa = components.Lipid('p', None, self.chainlength_unsaturated[chainlength_sn1], None, self._compartment_weights)
                del self.acyl_coa_list_unsaturated[sn1_chain]
            self.lyso_pa_list.append(lyso_pa)
            i = np.random.choice(choice_list, p=weights_pa)
            if i == 0:
                self.precursors_dict['DHAP'] -= 1
                self.precursors_dict['NADPH'] += 1
                self.precursors_dict['NADP'] -= 1
            else:
                self.precursors_dict['glycerol-3-p'] -= 1

    def PA_synthase(self):
        """
        Synthesis of PA by adding the second fatty acid to lyso_PA (sn2: saturated or unsaturated) -->
        1-acyl-sn-glycerol-3-phosphate acyltransferase
        """
        x = np.random.random()
        if x >= self.probabilities['PA_synthase'] and len(self.acyl_coa_list_unsaturated) > 1 and len(self.lyso_pa_list) > 1:
            z = np.random.randint(0, (len(self.lyso_pa_list)-1))
            sn2_chain = np.random.randint(0, (len(self.acyl_coa_list_unsaturated)-1))
            chainlength_sn2 = self.acyl_coa_list_unsaturated[sn2_chain].C
            self.lyso_pa_list[z].sn2 = self.chainlength_unsaturated[chainlength_sn2]
            self.PA_list.append(self.lyso_pa_list[z])
            # deletion of the consumed ffa
            del self.acyl_coa_list_unsaturated[sn2_chain]
            del self.lyso_pa_list[z]

    def CDP_DG_synthase(self):
        """
        PA is processed to CDP-DG (CDP-diacylglycerol synthase), that further reacts to the phospholipids
        """
        for i in range(self._rates['CDP_DG_synthase']):
            x = np.random.random()
            if x >= self.probabilities['CDP_DG_synthase'] and self.precursors_dict['CTP'] > 1 and len(self.PA_list) > 1:
                z = np.random.randint(0, len(self.PA_list)-1)
                self.PA_list[z].head = 'cdp'
                # CDP-DG production from PA
                self.CDP_DG_list.append(self.PA_list[z])
                del self.PA_list[z]
                self.precursors_dict['CTP'] -= 1
                self.precursors_dict['Pi'] += 2

    def TAG_synthase(self):
        """
        Function for TAG synthesis divided in production of DAG and TAG afterwards
        """
        for i in range(self._rates['TAG_synthase']):
            self.DAG_synthase()
            self.TAG_synthase()

    def DAG_synthase(self):
        """
        DAG synthesis: Removing the head of the lipid and adding the lipid to the DAG list.
        """
        x = np.random.random()
        if x >= self.probabilities['DAG_synthase'] and len(self.PA_list) > 1:
            z = np.random.randint(0, len(self.PA_list)-1)
            self.PA_list[z].head = None
            self.DAG_list.append(self.PA_list[z])
            self.precursors_dict['H2O'] -= 1
            self.precursors_dict['Pi'] += 1
            del self.PA_list[z]

    def TAG_synthase(self):
        """
        DAG is processed to TAG by adding a third acyl-chain at position sn3.
        """
        x = np.random.random()
        if x >= self.probabilities['TAG_synthase'] and len(self.DAG_list) > 1 and len(self.acyl_coa_list_saturated) > 1 and len(self.acyl_coa_list_unsaturated) > 1:
            z = np.random.randint(0, len(self.DAG_list)-1)
            self.TAG_list.append(self.DAG_list[z])
            self.TAG_list[-1].__class__ = components.TAG
            if x <= 0.575:
                sn3 = np.random.randint(0, len(self.acyl_coa_list_saturated)-1)
                chainlength_sn3 = self.acyl_coa_list_saturated[sn3].C
                self.TAG_list[-1].sn3 = self.chainlength_saturated[chainlength_sn3]
                del self.acyl_coa_list_saturated[sn3]
            else:
                sn3 = np.random.randint(0, len(self.acyl_coa_list_unsaturated)-1)
                chainlength_sn3 = self.acyl_coa_list_unsaturated[sn3].C
                self.TAG_list[-1].sn3 = self.chainlength_unsaturated[chainlength_sn3]
                del self.acyl_coa_list_unsaturated[sn3]
            del self.DAG_list[z]

    def TAG_lipase(self):
        """
        Cdk1/Cdc28-dependent activation of the major triacylglycerol lipase
        """
        if len(self.lipid_droplets) > self._rates['TAG_lipase']:
            for i in range(self._rates['TAG_lipase']):
                x = np.random.random()
                if x >= self.probabilities['TAG_lipase']:
                    z = np.random.randint(0, len(self.lipid_droplets)-1)
                    if self.lipid_droplets[z].head is None:
                        if ':0' in self.lipid_droplets[z].sn3:
                            for key, value in self.chainlength_unsaturated.items():
                                if value == self.lipid_droplets[z].sn3:
                                    self.acyl_coa_list_saturated.append(components.FattyAcid(key, 0))
                        elif ':1' in self.lipid_droplets[z].sn3:
                            for key, value in self.chainlength_saturated.items():
                                if value == self.lipid_droplets[z].sn3:
                                    self.acyl_coa_list_unsaturated.append(components.FattyAcid(key, 1))
                        self.DAG_list.append(self.lipid_droplets[z])
                        self.DAG_list[-1].__class__ = components.Lipid
                        delattr(self.DAG_list[-1], 'sn3')
                        self.precursors_dict['H2O'] -= 1
                    elif self.lipid_droplets[z].head == 'sterylester':
                        self.Ergosterol_list.append(components.Sterol('sterol', None, self._compartment_weights))
                        self.precursors_dict['H2O'] -= 1
                        if ':0' in self.lipid_droplets[z].FA:
                            for key, value in self.chainlength_unsaturated.items():
                                if value == self.lipid_droplets[z].FA:
                                    self.acyl_coa_list_saturated.append(components.FattyAcid(key, 0))
                        elif ':1' in self.lipid_droplets[z].FA:
                            for key, value in self.chainlength_saturated.items():
                                if value == self.lipid_droplets[z].FA:
                                    self.acyl_coa_list_unsaturated.append(components.FattyAcid(key, 1))
                    del self.lipid_droplets[z]

    def DAG_kinase(self):
        """
        DAG kinase
        """
        if len(self.DAG_list) > self._rates['DAG_kinase']:
            for i in range(self._rates['DAG_kinase']):
                x = np.random.random()
                if x >= self.probabilities['DAG_kinase']:
                    z = np.random.randint(0, len(self.DAG_list)-1)
                    self.PA_list.append(self.DAG_list[z])
                    self.PA_list[-1].head = 'p'
                    self.PA_list[-1].comp = None
                    del self.DAG_list[z]

    def PS_synthase(self):
        """
        CDP-DG is processed to PS (PS synthase).
        """
        for i in range(self._rates['PS_synthase']):
            x = np.random.random()
            if x >= self.probabilities['PS_synthase'] and len(self.CDP_DG_list) > 1 and self.precursors_dict['serine'] > 1:
                z = np.random.randint(0, len(self.CDP_DG_list)-1)
                # PS synthesis from CDP-DG
                self.CDP_DG_list[z].head = 'serine'
                self.PS_list.append(self.CDP_DG_list[z])
                del self.CDP_DG_list[z]
                self.precursors_dict['serine'] -= 1
                self.precursors_dict['CMP'] += 1

    def PI_synthase(self):
        """
        CDP-DG is processed to PI (PI synthase)
        """
        for i in range(self._rates['PI_synthase']):
            x = np.random.random()
            if x >= self.probabilities['PI_synthase'] and len(self.CDP_DG_list) > 1 and self.precursors_dict['inositol'] > 1:
                z = np.random.randint(0, len(self.CDP_DG_list)-1)
                self.CDP_DG_list[z].head = 'inositol'
                self.PI_list.append(self.CDP_DG_list[z])
                del self.CDP_DG_list[z]
                self.precursors_dict['inositol'] -= 1
                self.precursors_dict['CMP'] += 1

    def PE_synthase(self):
        """
        PE is derived from PS by releasing 1 CO2 --> PS decarboxylase.
        """
        for i in range(self._rates['PE_synthase']):
            x = np.random.random()
            if x >= self.probabilities['PE_synthase'] and len(self.PS_list) >= 10:
                z = np.random.randint(0, len(self.PS_list)-1)
                self.PS_list[z].head = 'ethanolamine'
                self.PE_list.append(self.PS_list[z])
                self.precursors_dict['CO2'] += 1
                del self.PS_list[z]

    def PC_synthase(self):
        """
        PC is derived from PE. As enzymes serve 3 methyltransferases which need SAM and produce SAH as a side product.
        """
        for i in range(self._rates['PC_synthase']):
            x = np.random.random()
            if x >= self.probabilities['PC_synthase'] and len(self.PE_list) >= 5 and self.precursors_dict['SAM'] >= 4:
                z = np.random.randint(0, len(self.PE_list)-1)
                self.PE_list[z].head = 'choline'
                self.PC_list.append(self.PE_list[z])
                del self.PE_list[z]
                self.precursors_dict['SAM'] -= 3
                self.precursors_dict['SAH'] += 3

    def CL_synthase(self):
        """
        Synthesis of cardiolipin, for which 2 CDP-DG are needed. Different enzymes are needed.
        """
        for i in range(self._rates['CL_synthase']):
            x = np.random.random()
            if x >= self.probabilities['CL_synthase'] and self.precursors_dict['glycerol_3_p_mito'] > 1 and len(self.CDP_DG_list) > 2:
                z = np.random.randint(0, len(self.CDP_DG_list)-2)
                self.CDP_DG_list[z].head = 'neutral'
                self.CL_list.append(self.CDP_DG_list[z])
                self.CL_list[-1].__class__ = components.CL
                self.CL_list[-1].sn4, self.CL_list[-1].sn3 = self.CDP_DG_list[z+1].sn2, self.CDP_DG_list[z+1].sn1
                del self.CDP_DG_list[z:z+1]
                self.precursors_dict['glycerol_3_p_mito'] -= 1
                self.precursors_dict['H2O'] -= 1
                self.precursors_dict['Pi'] += 1
                self.precursors_dict['CMP'] += 2

    def Ergosterol_synthase(self):
        """
        Synthesis of the most existing sterol in yeast: ergosterol
        """
        for i in range(self._rates['Ergosterol_synthase']):
            x = np.random.random()
            if x >= self.probabilities['Ergosterol_synthase'] and self.precursors_dict['acetyl_coa'] > 18:
                self.Ergosterol_list.append(components.Sterol('sterol', None, self._compartment_weights))
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
        """
        Synthesis of sterylesters that are found in lipid droplets out of ergosterol and an unsaturated fatty acid.
        """
        for i in range(self._rates['Sterylester_synthase']):
            x = np.random.random()
            if x >= self.probabilities['Sterylester_synthase'] and any(fa.C == 16 for fa in self.acyl_coa_list_unsaturated) and \
                    any(fa.C == 18 for fa in self.acyl_coa_list_unsaturated) and len(self.Ergosterol_list) > 1:
                z = np.random.randint(0, len(self.Ergosterol_list)-1)
                j = 0
                while j < 5:
                    fa_index = np.random.randint(0, len(self.acyl_coa_list_unsaturated)-1)
                    if self.acyl_coa_list_unsaturated[fa_index].C == 18 and np.random.random() < 0.33:
                        self.Sterylester_list.append(components.Sterylester('sterylester', 'C18:1', None, self._compartment_weights))
                        del self.Ergosterol_list[z]
                        del self.acyl_coa_list_unsaturated[fa_index]
                        break
                    elif self.acyl_coa_list_unsaturated[fa_index].C == 16:
                        self.Sterylester_list.append(components.Sterylester('sterylester', 'C16:1', None, self._compartment_weights))
                        del self.Ergosterol_list[z]
                        del self.acyl_coa_list_unsaturated[fa_index]
                        break
                    else:
                        j += 1

    def Sphingolipid_synthase(self):
        """
        Synthesis of the most abundant Sphingolipid mannose-(inositol-phosphate)2-ceramide
        """
        for i in range(self._rates['Sphingolipid_synthase']):
            x = np.random.random()
            if x >= self.probabilities['Sphingolipid_synthase'] and len(self.PI_list) >= 2 and self.precursors_dict['ceramide'] > 1 and \
                    self.precursors_dict['GDP-mannose'] > 1:
                self.Sphingolipid_list.append(components.Sphingolipid('ceramide', None, self._compartment_weights))
                z= np.random.randint(0, len(self.PI_list)-2)
                del self.PI_list[z:z+1]
                self.precursors_dict['ceramide'] -= 1
                self.precursors_dict['GDP-mannose'] -= 1

    def transport(self):
        """
        General transport function for all produced lipids.
        """
        for lipid in self.transport_list:
            if lipid == self.TAG_list or lipid == self.Sterylester_list:
                if len(lipid) > 10:
                    for j in range(len(lipid)/10):
                        z = np.random.randint(0, len(lipid)-1)
                        lipid[z].comp_choice()
                        if lipid[z].comp == 'lipid_droplets':
                            self.lipid_droplets.append(lipid[z])
                        del lipid[z]
            else:
                if len(lipid) > 5:
                    for j in range(len(lipid)/10):
                        z = np.random.randint(0, len(lipid)-1)
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
        """
        Function to calculate the lipid composition of all membranes.
        """
        x = 0
        for comp in self.compartment_lists:
            if len(comp) > 0:
                self.relatives_list = [(float(sum(j.head == 'serine' for j in comp)) / len(comp)),
                                        (float(sum(j.head == 'inositol' for j in comp)) / len(comp)),
                                        (float(sum(j.head == 'choline' for j in comp)) / len(comp)),
                                        (float(sum(j.head == 'ethanolamine' for j in comp)) / len(comp)),
                                        (float(sum(j.head == 'neutral' for j in comp)) / len(comp)),
                                        (float(sum(j.head == 'p' for j in comp)) / len(comp)),
                                        (float(sum(j.head == 'sterol' for j in comp)) / len(comp)),
                                        (float(sum(j.head == 'sterylester' for j in comp)) / len(comp)),
                                        (float(sum(j.head is None for j in comp)) / len(comp)),
                                        (float(sum(j.head == 'ceramide' for j in comp)) / len(comp))]
                for i in range(len(self.relatives_list)):
                    self.compartment_relatives_dict[self.compartment[x]][self.membrane_lipids[i]] = self.relatives_list[i]
            x += 1

    def numbers(self):
        for current_lipid_number, number_of_lipid in zip(self.number_lipids_list, self.components_list):
            current_lipid_number.append(len(number_of_lipid))
        # for plotting the number of lipids in a certain membrane
        for current_membrane_number, number_of_membrane in zip(self.number_membranes_list, self.compartment_lists):
            current_membrane_number.append(len(number_of_membrane))

        for precursor in self.precursors.keys():
            self.precursors[precursor].append(self.precursors_dict[precursor])

    def saturation_counter(self):
        """
        composition of fatty acids that should be reached: C16:0 = 10%, C16:1 = 30%, C18:0 = 10%, C18:1 = 50%
        these numbers are from Klug & Daum 2013 'Yeast lipid metabolism at a glance'
        """
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
        self.saturation_composition_sn1 = {'C16:0': self.c16_0_sn1, 'C16:1': self.c16_1_sn1, 'C18:0': self.c18_0_sn1,
                                           'C18_1': self.c18_1_sn1}

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

        self.saturation_composition_sn2 = {'C16:0': self.c16_0_sn2, 'C16:1': self.c16_1_sn2, 'C18:0': self.c18_0_sn2,
                                           'C18_1': self.c18_1_sn2}
        self.total_fatty_acids = self.c16_0_sn1 + self.c16_1_sn1 + self.c18_0_sn1 + self.c18_1_sn1 + self.c16_0_sn2 + self.c16_1_sn2 + \
            self.c18_0_sn2 + self.c18_1_sn2

        self.sterylester_C16 = 0
        self.sterylester_C18 = 0
        for c in self.lipid_droplets:
            if c.head == 'sterylester':
                if c.FA == 'C16:1':
                    self.sterylester_C16 += 1
                elif c.FA == 'C18:1':
                    self.sterylester_C18 += 1

        if self.sterylester_C16 > 0 or self.sterylester_C18 > 0:
            self.composition_sterylester = {'C16:1: ': float(self.sterylester_C16) / (self.sterylester_C16 + self.sterylester_C18),
                                            'C18:1: ': float(self.sterylester_C18) / (self.sterylester_C16 + self.sterylester_C18)}

        if self.total_fatty_acids > 0:
            self.saturation_composition_total = {'C16:0': float(self.c16_0_sn2 + self.c16_0_sn1) / self.total_fatty_acids,
                                                 'C16:1': float(self.c16_1_sn2 + self.c16_1_sn1) / self.total_fatty_acids,
                                                 'C18:0': float(self.c18_0_sn2 + self.c18_0_sn1) / self.total_fatty_acids,
                                                 'C18:1': float(self.c18_1_sn2 + self.c18_1_sn1) / self.total_fatty_acids}

    def plot_precursors(self):
        """
        Plotting the precursor molecules from the precursors_dict.
        """
        fig = mat.figure()
        ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
        i = 0
        while i < 6:  # in range(len(self.precursors.keys())):
            ax.plot(self.t, self.precursors[self.precursors.keys()[i]], label=self.precursors.keys()[i])
            i += 1
        ax.plot(self.t, self.precursors['inositol'], label='inositol')
        ax.plot(self.t, self.precursors['serine'], label='serine')
        ax.plot(self.t, self.precursors['CTP'], label='CTP')
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        mat.show()

    def plot_lipids(self):
        """
        Plotting the produced lipids before they are transported into the membranes
        """
        fig = mat.figure()
        ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
        ax.plot(self.t, self.number_tag, label='tag')
        ax.plot(self.t, self.number_PS, label='ps')
        ax.plot(self.t, self.number_PI, label='pi')
        ax.plot(self.t, self.number_PE, label='pe')
        ax.plot(self.t, self.number_PC, label='pc')
        ax.plot(self.t, self.number_CL, label='cl')
        ax.plot(self.t, self.number_Ergosterol, label='es')
        ax.plot(self.t, self.number_Sterylester, label='se')
        ax.plot(self.t, self.number_Sphingolipid, label='sl')
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        mat.show()

    def plot_precursor_lipids(self):
        """
        Plotting some free precursor molecules.
        """
        fig = mat.figure()
        ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
        ax.plot(self.t, self.number_acetyl_coa[:-1], label='acetyl_coa')
        ax.plot(self.t, self.number_acyl_coa[:-1], label='acyl_coa')
        ax.plot(self.t, self.number_pa[:-1], label='pa')
        ax.plot(self.t, self.number_cdp_dg[:-1], label='cdp-dg')
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        mat.show()

    def plot_membranes(self):
        """
        Plotting the number of lipids in the membranes of different compartments.
        """
        fig = mat.figure()
        ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
        ax.plot(self.t, self.number_plasma_membrane, label='plasma membrane')
        ax.plot(self.t, self.number_secretory_vesicles, label='secretory vesicles')
        ax.plot(self.t, self.number_vacuoles, label='vacuoles')
        ax.plot(self.t, self.number_nucleus, label='nucleus')
        ax.plot(self.t, self.number_peroxisomes, label='peroxisomes')
        ax.plot(self.t, self.number_light_microsomes, label='light_microsomes')
        ax.plot(self.t, self.number_inner_mit_membrane, label='inner_mit_membrane')
        ax.plot(self.t, self.number_outer_mit_membrane, label='outer_mit_membrane')
        ax.plot(self.t, self.number_lipid_droplets, label='lipid_droplets')
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        mat.show()


if __name__ == '__main__':
    # Test run, with runtime tracker
    import time
    st = time.time()
    m = Model()
    # test run: 5 sec
    m.run(5)
    et = time.time()
    print "Runtime: " + str(et - st) + "s"
