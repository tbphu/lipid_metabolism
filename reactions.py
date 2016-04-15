import numpy as np
import components


class Reactions:
    def __init__(self, state, rates, probabilities, precursors):
        """
        Description

        Parameters
        ----------
        precursors: dict
            number of precursors to add
        state: dict
            current membrane states

        Return
        ------
        state: dict
            actual list of objects - membrane compositions
        """
        self.rates = rates
        self.probabilities = probabilities
        self.precursors_dict = precursors

        self.run_step()

    def run_step(self):
        function_list = [self.glycerol_3_p_synthesis,
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
        for i in function_list:
            func = np.random.choice(function_list)
            func()
            function_list.remove(func)

    def glycerol_3_p_synthesis(self):
        '''
        Synthesis of glycerol-3-p out of DHAP.
        '''
        for i in range(self.rates['glycerol_3_p_synthesis']):
            x = np.random.random()
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
            x = np.random.random()
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
            x = np.random.random()
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
            x = np.random.random()
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
            x = np.random.random()						#5 reactions in 1 timestep but only with a probability of 90%
            if self.precursors_dict['acetyl_coa'] > 2:		#control if at least 2 Acetyl-CoA are available
                if len(self.acyl_coa_list) == 0:		#starting the first reaction
                    new_acyl = components.FattyAcid(2, np.random.choice(choice_list, p = choice_weights))
                    self.acyl_coa_list.append(new_acyl)
                    self.acyl_coa_list[-1].C += 2
                    self.precursors_dict['acetyl_coa'] -= 2

                elif self.acyl_coa_list[-1].C == 16 and x >= self.probabilities['acyl_synthase_C16']:	#stop the reaction cycle and starting a new one
                    new_acyl = components.FattyAcid(2, np.random.choice(choice_list, p = choice_weights))
                    self.acyl_coa_list.append(new_acyl)
                    self.acyl_coa_list[-1].C += 2
                    self.precursors_dict['acetyl_coa'] -= 2
                    self.precursors_dict['NADPH'] -= 14
                    self.precursors_dict['NADP'] += 14
                    self.precursors_dict['H2O'] += 7
                    #CO2 production is not mentioned here as onyl acetyl-CoA is used and not malonyl-CoA, so we need all C-atoms we give in the reaction

                elif self.acyl_coa_list[-1].C == 18 and x >= self.probabilities['acyl_synthase_C18']:	#stop the reaction cycle and starting a new one
                    new_acyl = components.FattyAcid(2, np.random.choice(choice_list, p = choice_weights))
                    self.acyl_coa_list.append(new_acyl)
                    self.acyl_coa_list[-1].C += 2
                    self.precursors_dict['acetyl_coa'] -= 2
                    self.precursors_dict['NADPH'] -= 16
                    self.precursors_dict['NADP'] += 16
                    self.precursors_dict['H2O'] += 8

                elif self.acyl_coa_list[-1].C == 26:
                    self.acyl_coa_list[-1].saturation = 0
                    new_acyl = components.FattyAcid(2, np.random.choice(choice_list, p = choice_weights))
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
        x = np.random.random()
        if x >= self.probabilities['lyso_PA_synthase'] and len(self.acyl_coa_list_saturated) > 1 and len(self.acyl_coa_list_unsaturated) > 1 and (self.precursors_dict['DHAP'] > 1 and self.precursors_dict['glycerol-3-p'] > 1): 	#at least 1 ffa has to be unsaturated
            if np.random.choice(choice_list, p = self.weights_fa) == 0:
                sn1_chain = np.random.randint(0, (len(self.acyl_coa_list_saturated)-1))
                chainlength_sn1 = self.acyl_coa_list_saturated[sn1_chain].C
                lyso_pa = components.Lipid('p', None, self.chainlength_saturated[chainlength_sn1], None, self.compartment_weights)
                del self.acyl_coa_list_saturated[sn1_chain]
            else:
                sn1_chain = np.random.randint(0, (len(self.acyl_coa_list_unsaturated)-1))
                chainlength_sn1 = self.acyl_coa_list_unsaturated[sn1_chain].C
                lyso_pa = components.Lipid('p', None, self.chainlength_unsaturated[chainlength_sn1], None, self.compartment_weights)
                del self.acyl_coa_list_unsaturated[sn1_chain]
            self.lyso_pa_list.append(lyso_pa)
            i = np.random.choice(choice_list, p = weights_pa)
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
        x = np.random.random()
        if x >= self.probabilities['PA_synthase'] and len(self.acyl_coa_list_unsaturated) > 1 and len(self.lyso_pa_list) > 1:
            z = np.random.randint(0, (len(self.lyso_pa_list)-1))
            sn2_chain = np.random.randint(0, (len(self.acyl_coa_list_unsaturated)-1))
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
            x = np.random.random()
            if x >= self.probabilities['CDP_DG_synthase'] and self.precursors_dict['CTP'] > 1 and len(self.PA_list) > 1:
                z = np.random.randint(0, len(self.PA_list)-1)
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
        x = np.random.random()
        if x >= self.probabilities['DAG_synthase'] and len(self.PA_list) > 1:
            z = np.random.randint(0, len(self.PA_list)-1)
            self.PA_list[z].head = None
            self.DAG_list.append(self.PA_list[z])
            self.precursors_dict['H2O'] -= 1
            self.precursors_dict['Pi'] += 1
            del self.PA_list[z]

    def TAG_synthase(self):
        '''
        DAG is processed to TAG by adding a third acyl-chain at position sn3.
        '''
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
        '''
        Cdk1/Cdc28-dependent activation of the major triacylglycerol lipase
        '''
        if len(self.lipid_droplets) > self.rates['TAG_lipase']:
            for i in range(self.rates['TAG_lipase']):
                x = np.random.random()
                if x >= self.probabilities['TAG_lipase']:
                    z = np.random.randint(0, len(self.lipid_droplets)-1)
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
                        self.DAG_list[-1].__class__ = components.Lipid
                        delattr(self.DAG_list[-1], 'sn3')
                        self.precursors_dict['H2O'] -= 1
                    elif self.lipid_droplets[z].head == 'sterylester':
                        self.Ergosterol_list.append(components.Sterol('sterol', None, self.compartment_weights))
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
                x = np.random.random()
                if x >= self.probabilities['DAG_kinase']:
                    z = np.random.randint(0, len(self.DAG_list)-1)
                    self.PA_list.append(self.DAG_list[z])
                    self.PA_list[-1].head = 'p'
                    self.PA_list[-1].comp = None
                    del self.DAG_list[z]


    def PS_synthase(self):
        '''
        CDP-DG is processed to PS (PS synthase).
        '''
        for i in range(self.rates['PS_synthase']):
            x = np.random.random()
            if x >= self.probabilities['PS_synthase'] and len(self.CDP_DG_list) > 1 and self.precursors_dict['serine'] > 1:
                z = np.random.randint(0, len(self.CDP_DG_list)-1)
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
            x = np.random.random()
            if x >= self.probabilities['PI_synthase'] and len(self.CDP_DG_list) > 1 and self.precursors_dict['inositol'] > 1:
                z = np.random.randint(0, len(self.CDP_DG_list)-1)
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
            x = np.random.random()
            if x >= self.probabilities['PE_synthase'] and len(self.PS_list) >= 10:
                z = np.random.randint(0, len(self.PS_list)-1)
                self.PS_list[z].head = 'ethanolamine'
                self.PE_list.append(self.PS_list[z])
                self.precursors_dict['CO2'] += 1
                del self.PS_list[z]


    def PC_synthase(self):
        '''
        PC is derived from PE. As enzymes serve 3 methyltransferases which need SAM and produce SAH as a side product.
        '''
        for i in range(self.rates['PC_synthase']):
            x = np.random.random()
            if x >= self.probabilities['PC_synthase'] and len(self.PE_list) >= 5 and self.precursors_dict['SAM'] >= 4:
                z = np.random.randint(0, len(self.PE_list)-1)
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
            x = np.random.random()
            if x >= self.probabilities['CL_synthase'] and self.precursors_dict['glycerol_3_p_mito'] > 1 and len(self.CDP_DG_list) > 2:
                z = np.random.randint(0, len(self.CDP_DG_list)-2)
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
            x = np.random.random()
            if x >= self.probabilities['Ergosterol_synthase'] and self.precursors_dict['acetyl_coa'] > 18:
                self.Ergosterol_list.append(components.Sterol('sterol', None, self.compartment_weights))
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
            x = np.random.random()
            if x >= self.probabilities['Sterylester_synthase'] and any(fa.C == 16 for fa in self.acyl_coa_list_unsaturated) and any(fa.C == 18 for fa in self.acyl_coa_list_unsaturated) and len(self.Ergosterol_list) > 1:
                z = np.random.randint(0, len(self.Ergosterol_list)-1)
                j = 0
                while j < 5:
                    fa_index = np.random.randint(0, len(self.acyl_coa_list_unsaturated)-1)
                    if self.acyl_coa_list_unsaturated[fa_index].C == 18 and np.random.random() < 0.33:
                        self.Sterylester_list.append(components.Sterylester('sterylester', 'C18:1', None, self.compartment_weights))
                        del self.Ergosterol_list[z]
                        del self.acyl_coa_list_unsaturated[fa_index]
                        break
                    elif self.acyl_coa_list_unsaturated[fa_index].C == 16:
                        self.Sterylester_list.append(components.Sterylester('sterylester', 'C16:1', None, self.compartment_weights))
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
            x = np.random.random()
            if x >= self.probabilities['Sphingolipid_synthase'] and len(self.PI_list) >= 2 and self.precursors_dict['ceramide'] > 1 and self.precursors_dict['GDP-mannose'] > 1:
                self.Sphingolipid_list.append(sphingolipid('ceramide', None, self.compartment_weights))
                z= np.random.randint(0, len(self.PI_list)-2)
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