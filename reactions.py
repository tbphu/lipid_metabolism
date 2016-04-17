import numpy as np
import components


class Reactions:
    def __init__(self, rates, weights):
        """
        Description

        Parameters
        ----------
        rates: dict
            VMAX values of all reactions
        weights: dict
            weights for FA properties and transport reactions
        """
        self.RATES = rates
        self.WEIGHTS = weights
        self.probabilities = None
        self.components_state = None
        self.membranes_state = None
        self.precursors_state = None

    def run_step(self, probabilities, component_states, membrane_states, precursor_states):
        """
        Simulate a time step.

        Parameters
        ----------
        probabilities: dict
            probabilities of all reactions
        component_states: dict
            current cell state - lipids & FAs
        membrane_states: dict
            current cell state - membrane compositions
        precursor_states: dict
            current cell state - precursors
        """
        self.probabilities = probabilities
        self.components_state = component_states
        self.membranes_state = membrane_states
        self.precursors_state = precursor_states

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

    def glycerol_3_p_synthesis(self):
        '''
        Synthesis of glycerol-3-p out of DHAP.
        '''
        for i in range(self.RATES['glycerol_3_p_synthesis']):
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
        for i in range(self.RATES['inositol_synthesis']):
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
        for i in range(self.RATES['ceramide_synthesis']):
            x = np.random.random()
            if x >= self.probabilities['ceramide_synthesis']:
                if len(self.components_state['acyl_coa_C26']) > 1 and self.precursors_state['serine'] > 1 and len(
                        self.components_state['acyl_coa_saturated']) > 1 and any(
                                fa.C == 16 for fa in self.components_state['acyl_coa_saturated']):
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
        for i in range(self.RATES['acetyl_coa_synthase']):
            x = np.random.random()
            if x >= self.probabilities['acetyl_coa_synthase']:
                if self.precursors_state['pyruvate'] > 1:  # transformation from pyruvate to acetyl_coa
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
        for i in range(self.RATES['acyl_synthase']):
            x = np.random.random()  # 5 reactions in 1 timestep but only with a probability of 90%
            if self.precursors_state['acetyl_coa'] > 2:  # control if at least 2 Acetyl-CoA are available
                if len(self.components_state['acyl_coa']) == 0:  # starting the first reaction
                    new_acyl = components.FattyAcid(2, np.random.choice(choice_list, p=choice_weights))
                    self.components_state['acyl_coa'].append(new_acyl)
                    self.components_state['acyl_coa'][-1].C += 2
                    self.precursors_state['acetyl_coa'] -= 2

                elif self.components_state['acyl_coa'][-1].C == 16 and x >= self.probabilities[
                    'acyl_synthase_C16']:  # stop the reaction cycle and starting a new one
                    new_acyl = components.FattyAcid(2, np.random.choice(choice_list, p=choice_weights))
                    self.components_state['acyl_coa'].append(new_acyl)
                    self.components_state['acyl_coa'][-1].C += 2
                    self.precursors_state['acetyl_coa'] -= 2
                    self.precursors_state['NADPH'] -= 14
                    self.precursors_state['NADP'] += 14
                    self.precursors_state['H2O'] += 7
                    # CO2 production is not mentioned here as onyl acetyl-CoA is used and not malonyl-CoA, so we need all C-atoms we give in the reaction

                elif self.components_state['acyl_coa'][-1].C == 18 and x >= self.probabilities[
                    'acyl_synthase_C18']:  # stop the reaction cycle and starting a new one
                    new_acyl = components.FattyAcid(2, np.random.choice(choice_list, p=choice_weights))
                    self.components_state['acyl_coa'].append(new_acyl)
                    self.components_state['acyl_coa'][-1].C += 2
                    self.precursors_state['acetyl_coa'] -= 2
                    self.precursors_state['NADPH'] -= 16
                    self.precursors_state['NADP'] += 16
                    self.precursors_state['H2O'] += 8

                elif self.components_state['acyl_coa'][-1].C == 26:
                    self.components_state['acyl_coa'][-1].saturation = 0
                    new_acyl = components.FattyAcid(2, np.random.choice(choice_list, p=choice_weights))
                    self.components_state['acyl_coa'].append(new_acyl)
                    self.components_state['acyl_coa'][-1].C += 2
                    self.precursors_state['acetyl_coa'] -= 2
                    self.precursors_state['NADPH'] -= 24
                    self.precursors_state['NADP'] += 24
                    self.precursors_state['H2O'] += 12

                else:  # adding an Acetyl_CoA to the growing ffa
                    self.components_state['acyl_coa'][-1].C += 2
                    self.precursors_state['acetyl_coa'] -= 1

        if len(self.components_state['acyl_coa']) > 1:
            for j in range(len(self.components_state['acyl_coa']) - 1):
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
        for i in range(self.RATES['PA_synthesis']):
            self.lyso_PA_synthase()
            self.PA_synthase()

    def lyso_PA_synthase(self):
        '''
        Production of Lyso-PA by adding one acyl-coa to DHAP (sn1: always unsaturated) --> DHAP acyltransferase/acyl-DHAP reductase
        '''
        choice_list = [0, 1]

        weights_pa = [self.precursors_state['DHAP'] / (self.precursors_state['DHAP'] + self.precursors_state['glycerol-3-p']), \
                      self.precursors_state['glycerol-3-p'] / (self.precursors_state['DHAP'] + self.precursors_state['glycerol-3-p'])]
        x = np.random.random()
        if x >= self.probabilities['lyso_PA_synthase'] and len(self.components_state['acyl_coa_saturated']) > 1 and len(
                self.components_state['acyl_coa_unsaturated']) > 1 and (
                self.precursors_state['DHAP'] > 1 and self.precursors_state['glycerol-3-p'] > 1):  # at least 1 ffa has to be unsaturated
            if np.random.choice(choice_list, p=self.WEIGHTS['FA']) == 0:
                sn1_chain = np.random.randint(0, (len(self.components_state['acyl_coa_saturated']) - 1))
                chainlength_sn1 = self.components_state['acyl_coa_saturated'][sn1_chain].C
                lyso_pa = components.Lipid('p', None, self.WEIGHTS['chain_saturated'][chainlength_sn1], None, self.WEIGHTS['compartments'])
                del self.components_state['acyl_coa_saturated'][sn1_chain]
            else:
                sn1_chain = np.random.randint(0, (len(self.components_state['acyl_coa_unsaturated']) - 1))
                chainlength_sn1 = self.components_state['acyl_coa_unsaturated'][sn1_chain].C
                lyso_pa = components.Lipid('p', None, self.WEIGHTS['chain_unsaturated'][chainlength_sn1], None,
                                           self.WEIGHTS['compartments'])
                del self.components_state['acyl_coa_unsaturated'][sn1_chain]
            self.components_state['lyso_PA'].append(lyso_pa)
            i = np.random.choice(choice_list, p=weights_pa)
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
        if x >= self.probabilities['PA_synthase'] and len(self.components_state['acyl_coa_unsaturated']) > 1 and len(
                self.components_state['lyso_PA']) > 1:
            z = np.random.randint(0, (len(self.components_state['lyso_PA']) - 1))
            sn2_chain = np.random.randint(0, (len(self.components_state['acyl_coa_unsaturated']) - 1))
            chainlength_sn2 = self.components_state['acyl_coa_unsaturated'][sn2_chain].C
            self.components_state['lyso_PA'][z].sn2 = self.WEIGHTS['chain_unsaturated'][chainlength_sn2]
            self.components_state['PA'].append(self.components_state['lyso_PA'][z])
            del self.components_state['acyl_coa_unsaturated'][sn2_chain]  # deletion of the consumed ffa
            del self.components_state['lyso_PA'][z]

    def CDP_DG_synthase(self):
        '''
        PA is processed to CDP-DG (CDP-diacylglycerol synthase), that further reacts to the phospholipids
        '''
        for i in range(self.RATES['CDP_DG_synthase']):
            x = np.random.random()
            if x >= self.probabilities['CDP_DG_synthase'] and self.precursors_state['CTP'] > 1 and len(self.components_state['PA']) > 1:
                z = np.random.randint(0, len(self.components_state['PA']) - 1)
                self.components_state['PA'][z].head = 'cdp'
                self.components_state['CDP_DG'].append(self.components_state['PA'][z])  # CDP-DG production from PA
                del self.components_state['PA'][z]
                self.precursors_state['CTP'] -= 1
                self.precursors_state['Pi'] += 2

    def TAG_synthese(self):
        '''
        Function for TAG synthesis divided in production of DAG and TAG afterwards
        '''
        for i in range(self.RATES['TAG_synthesis']):
            self.DAG_synthase()
            self.TAG_synthase()

    def DAG_synthase(self):
        '''
        DAG synthesis: Removing the head of the lipid and adding the lipid to the DAG list.
        '''
        x = np.random.random()
        if x >= self.probabilities['DAG_synthase'] and len(self.components_state['PA']) > 1:
            z = np.random.randint(0, len(self.components_state['PA']) - 1)
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
        if x >= self.probabilities['TAG_synthase'] and len(self.components_state['DAG']) > 1 and len(
                self.components_state['acyl_coa_saturated']) > 1 and len(self.components_state['acyl_coa_unsaturated']) > 1:
            z = np.random.randint(0, len(self.components_state['DAG']) - 1)
            self.components_state['TAG'].append(self.components_state['DAG'][z])
            self.components_state['TAG'][-1].__class__ = components.TAG
            if x <= 0.575:
                sn3 = np.random.randint(0, len(self.components_state['acyl_coa_saturated']) - 1)
                chainlength_sn3 = self.components_state['acyl_coa_saturated'][sn3].C
                self.components_state['TAG'][-1].sn3 = self.WEIGHTS['chain_saturated'][chainlength_sn3]
                del self.components_state['acyl_coa_saturated'][sn3]
            else:
                sn3 = np.random.randint(0, len(self.components_state['acyl_coa_unsaturated']) - 1)
                chainlength_sn3 = self.components_state['acyl_coa_unsaturated'][sn3].C
                self.components_state['TAG'][-1].sn3 = self.WEIGHTS['chain_unsaturated'][chainlength_sn3]
                del self.components_state['acyl_coa_unsaturated'][sn3]
            del self.components_state['DAG'][z]

    def TAG_lipase(self):
        '''
        Cdk1/Cdc28-dependent activation of the major triacylglycerol lipase
        '''
        if len(self.membranes_state['lipid_droplets']) > self.RATES['TAG_lipase']:
            for i in range(self.RATES['TAG_lipase']):
                x = np.random.random()
                if x >= self.probabilities['TAG_lipase']:
                    z = np.random.randint(0, len(self.membranes_state['lipid_droplets']) - 1)
                    if self.membranes_state['lipid_droplets'][z].head == None:
                        if ':0' in self.membranes_state['lipid_droplets'][z].sn3:
                            for key, value in self.WEIGHTS['chain_unsaturated'].items():
                                if value == self.membranes_state['lipid_droplets'][z].sn3:
                                    self.components_state['acyl_coa_saturated'].append(components.FattyAcid(key, 0))
                        elif ':1' in self.membranes_state['lipid_droplets'][z].sn3:
                            for key, value in self.WEIGHTS['chain_saturated'].items():
                                if value == self.membranes_state['lipid_droplets'][z].sn3:
                                    self.components_state['acyl_coa_unsaturated'].append(components.FattyAcid(key, 1))
                        self.components_state['DAG'].append(self.membranes_state['lipid_droplets'][z])
                        self.components_state['DAG'][-1].__class__ = components.Lipid
                        delattr(self.components_state['DAG'][-1], 'sn3')
                        self.precursors_state['H2O'] -= 1
                    elif self.membranes_state['lipid_droplets'][z].head == 'sterylester':
                        self.components_state['ergosterol'].append(components.Sterol('sterol', None, self.WEIGHTS['compartments']))
                        self.precursors_state['H2O'] -= 1
                        if ':0' in self.membranes_state['lipid_droplets'][z].FA:
                            for key, value in self.WEIGHTS['chain_unsaturated'].items():
                                if value == self.membranes_state['lipid_droplets'][z].FA:
                                    self.components_state['acyl_coa_saturated'].append(components.FattyAcid(key, 0))
                        elif ':1' in self.membranes_state['lipid_droplets'][z].FA:
                            for key, value in self.WEIGHTS['chain_saturated'].items():
                                if value == self.membranes_state['lipid_droplets'][z].FA:
                                    self.components_state['acyl_coa_unsaturated'].append(components.FattyAcid(key, 1))
                    del self.membranes_state['lipid_droplets'][z]

    def DAG_kinase(self):
        if len(self.components_state['DAG']) > self.RATES['DAG_kinase']:
            for i in range(self.RATES['DAG_kinase']):
                x = np.random.random()
                if x >= self.probabilities['DAG_kinase']:
                    z = np.random.randint(0, len(self.components_state['DAG']) - 1)
                    self.components_state['PA'].append(self.components_state['DAG'][z])
                    self.components_state['PA'][-1].head = 'p'
                    self.components_state['PA'][-1].comp = None
                    del self.components_state['DAG'][z]

    def PS_synthase(self):
        '''
        CDP-DG is processed to PS (PS synthase).
        '''
        for i in range(self.RATES['PS_synthase']):
            x = np.random.random()
            if x >= self.probabilities['PS_synthase'] and len(self.components_state['CDP_DG']) > 1 and self.precursors_state['serine'] > 1:
                z = np.random.randint(0, len(self.components_state['CDP_DG']) - 1)
                self.components_state['CDP_DG'][z].head = 'serine'  # PS synthesis from CDP-DG
                self.components_state['PS'].append(self.components_state['CDP_DG'][z])
                del self.components_state['CDP_DG'][z]
                self.precursors_state['serine'] -= 1
                self.precursors_state['CMP'] += 1

    def PI_synthase(self):
        '''
        CDP-DG is processed to PI (PI synthase)
        '''
        for i in range(self.RATES['PI_synthase']):
            x = np.random.random()
            if x >= self.probabilities['PI_synthase'] and len(self.components_state['CDP_DG']) > 1 and self.precursors_state[
                'inositol'] > 1:
                z = np.random.randint(0, len(self.components_state['CDP_DG']) - 1)
                self.components_state['CDP_DG'][z].head = 'inositol'
                self.components_state['PI'].append(self.components_state['CDP_DG'][z])
                del self.components_state['CDP_DG'][z]
                self.precursors_state['inositol'] -= 1
                self.precursors_state['CMP'] += 1

    def PE_synthase(self):
        '''
        PE is derived from PS by releasing 1 CO2 --> PS decarboxylase.
        '''
        for i in range(self.RATES['PE_synthase']):
            x = np.random.random()
            if x >= self.probabilities['PE_synthase'] and len(self.components_state['PS']) >= 10:
                z = np.random.randint(0, len(self.components_state['PS']) - 1)
                self.components_state['PS'][z].head = 'ethanolamine'
                self.components_state['PE'].append(self.components_state['PS'][z])
                self.precursors_state['CO2'] += 1
                del self.components_state['PS'][z]

    def PC_synthase(self):
        '''
        PC is derived from PE. As enzymes serve 3 methyltransferases which need SAM and produce SAH as a side product.
        '''
        for i in range(self.RATES['PC_synthase']):
            x = np.random.random()
            if x >= self.probabilities['PC_synthase'] and len(self.components_state['PE']) >= 5 and self.precursors_state['SAM'] >= 4:
                z = np.random.randint(0, len(self.components_state['PE']) - 1)
                self.components_state['PE'][z].head = 'choline'
                self.components_state['PC'].append(self.components_state['PE'][z])
                del self.components_state['PE'][z]
                self.precursors_state['SAM'] -= 3
                self.precursors_state['SAH'] += 3

    def CL_synthase(self):
        '''
        Synthesis of cardiolipin, for which 2 CDP-DG are needed. Different enzymes are needed.
        '''
        for i in range(self.RATES['CL_synthase']):
            x = np.random.random()
            if x >= self.probabilities['CL_synthase'] and self.precursors_state['glycerol_3_p_mito'] > 1 and len(
                    self.components_state['CDP_DG']) > 2:
                z = np.random.randint(0, len(self.components_state['CDP_DG']) - 2)
                self.components_state['CDP_DG'][z].head = 'neutral'
                self.components_state['CL'].append(self.components_state['CDP_DG'][z])
                self.components_state['CL'][-1].__class__ = components.CL
                self.components_state['CL'][-1].sn4, self.components_state['CL'][-1].sn3 = self.components_state['CDP_DG'][z + 1].sn2, \
                                                                                           self.components_state['CDP_DG'][z + 1].sn1
                del self.components_state['CDP_DG'][z:z + 1]
                self.precursors_state['glycerol_3_p_mito'] -= 1
                self.precursors_state['H2O'] -= 1
                self.precursors_state['Pi'] += 1
                self.precursors_state['CMP'] += 2

    def Ergosterol_synthase(self):
        '''
        Synthesis of the most existing sterol in yeast: ergosterol
        '''
        for i in range(self.RATES['ergosterol_synthase']):
            x = np.random.random()
            if x >= self.probabilities['ergosterol_synthase'] and self.precursors_state['acetyl_coa'] > 18:
                self.components_state['ergosterol'].append(components.Sterol('sterol', None, self.WEIGHTS['compartments']))
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
        for i in range(self.RATES['sterylester_synthase']):
            x = np.random.random()
            if x >= self.probabilities['sterylester_synthase'] and any(
                            fa.C == 16 for fa in self.components_state['acyl_coa_unsaturated']) and any(
                            fa.C == 18 for fa in self.components_state['acyl_coa_unsaturated']) and len(
                    self.components_state['ergosterol']) > 1:
                z = np.random.randint(0, len(self.components_state['ergosterol']) - 1)
                j = 0
                while j < 5:
                    fa_index = np.random.randint(0, len(self.components_state['acyl_coa_unsaturated']) - 1)
                    if self.components_state['acyl_coa_unsaturated'][fa_index].C == 18 and np.random.random() < 0.33:
                        self.components_state['sterylester'].append(
                            components.Sterylester('sterylester', 'C18:1', None, self.WEIGHTS['compartments']))
                        del self.components_state['ergosterol'][z]
                        del self.components_state['acyl_coa_unsaturated'][fa_index]
                        break
                    elif self.components_state['acyl_coa_unsaturated'][fa_index].C == 16:
                        self.components_state['sterylester'].append(
                            components.Sterylester('sterylester', 'C16:1', None, self.WEIGHTS['compartments']))
                        del self.components_state['ergosterol'][z]
                        del self.components_state['acyl_coa_unsaturated'][fa_index]
                        break
                    else:
                        j += 1

    def Sphingolipid_synthase(self):
        '''
        Synthesis of the most abundant Sphingolipid mannose-(inositol-phosphate)2-ceramide
        '''
        for i in range(self.RATES['sphingolipid_synthase']):
            x = np.random.random()
            if x >= self.probabilities['sphingolipid_synthase'] and len(self.components_state['PI']) >= 2 and self.precursors_state[
                'ceramide'] > 1 and self.precursors_state['GDP-mannose'] > 1:
                self.components_state['sphingolipid'].append(components.Sphingolipid('ceramide', None, self.WEIGHTS['compartments']))
                z = np.random.randint(0, len(self.components_state['PI']) - 2)
                del self.components_state['PI'][z:z + 1]
                self.precursors_state['ceramide'] -= 1
                self.precursors_state['GDP-mannose'] -= 1

    def transport(self):
        '''
        General transport function for all produced lipids.
        '''
        # lipids to transport
        transport_lists = [self.components_state['PS'], self.components_state['PI'], self.components_state['PC'],
                           self.components_state['PE'], self.components_state['CL'], self.components_state['PA'],
                           self.components_state['ergosterol'], self.components_state['sterylester'], self.components_state['TAG'],
                           self.components_state['sphingolipid']]

        for lipid in transport_lists:
            if lipid == self.components_state['TAG'] or lipid == self.components_state['sterylester']:
                if len(lipid) > 10:
                    for j in range(len(lipid) / 10):
                        z = np.random.randint(0, len(lipid) - 1)
                        lipid[z].comp_choice()
                        if lipid[z].comp == 'lipid_droplets':
                            self.membranes_state['lipid_droplets'].append(lipid[z])
                        del lipid[z]
            else:
                if len(lipid) > 5:
                    for j in range(len(lipid) / 10):
                        z = np.random.randint(0, len(lipid) - 1)
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
