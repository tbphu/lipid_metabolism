import random
import matplotlib.pyplot as pyp
import numpy as np
import math
import model
import copy
import pickle
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})


class StochModel:
    def __init__(self):
        self.model = model.Model()

        # make copy of parameter sets
        self.rates = copy.deepcopy(self.model.rates)
        self.probabilities = copy.deepcopy(self.model.probability)
        self.probabilities_G1 = copy.deepcopy(self.model.probability_G1)
        self.probabilities_S_M = copy.deepcopy(self.model.probability_S_M)
        self.comp_weights = copy.deepcopy(self.model.compartment_weights)
        self.fa_weights = copy.deepcopy(self.model.weights_fa)

        # calculation fatty acids
        self.fatty_acid_distribution = {}
        self.fatty_acid_std = {}
        self.fatty_acids = ['C16:0', 'C16:1', 'C18:0', 'C18:1']
        self.fa = ['saturated', 'unsaturated']

        # membrane compositions
        self.membrane_comp_mean = {}
        self.membrane_comp_std = {}
        self.lipids = ['CL', 'ES', 'PA', 'PC', 'PE', 'PI', 'PS', 'SE', 'SL', 'TAG']
        self.membranes_length = {}
        self.membranes_length_std = {}

        # membrane lengths
        self.plasma_membrane = []
        self.secretory_vesicles = []
        self.vacuoles = []
        self.nucleus = []
        self.peroxisomes = []
        self.light_microsomes = []
        self.inner_mit_membrane = []
        self.outer_mit_membrane = []
        self.lipid_droplets = []
        self.mean_time_lists = {}
        self.std_time_lists = {}
        self.membrane_lists = {'plasma_membrane': self.plasma_membrane, 'secretory_vesicles': self.secretory_vesicles,
                               'vacuoles': self.vacuoles, 'nucleus': self.nucleus, 'peroxisomes': self.peroxisomes,
                               'light_microsomes': self.light_microsomes, 'inner_mit_membrane': self.inner_mit_membrane,
                               'outer_mit_membrane': self.outer_mit_membrane, 'lipid_droplets': self.lipid_droplets}

        self.membranes_list = ['plasma_membrane', 'secretory_vesicles', 'vacuoles', 'nucleus', 'peroxisomes',
                               'light_microsomes', 'inner_mit_membrane', 'outer_mit_membrane', 'lipid_droplets']

        self.time = []

    def run(self, runs=10, sensitivity=False, change=0.1):
        # simulate
        # sensitivity analysis
        if sensitivity:
            self.results = {}
            self.sensitivity_analysis(runs, change)
            self.sensitivity_statistics()
        # plot heatmaps of sesitivity analysis
            # self.heatmap_fatty_acids()
            # self.heatmap_membrane_comp()
            self.heatmap_membrane_length()

    def sensitivity_analysis(self, runs, change):

        self.results_keys = pickle.load(open("./sens_analysis_start_10-12-15/results_sens_keys.pkl", "rb"))
        for i in self.results_keys:
            self.results[i] = pickle.load(open("./sens_analysis_start_10-12-15/"+i+".pkl", "rb"))
        self.time = [i for i in range(7200)]


    def sensitivity_statistics(self):
        #fatty acids distribution (mean and std for each change)
        for reaction in self.results:
            self.fatty_acid_distribution[reaction] = {}
            self.fatty_acid_std[reaction] = {}
            for fa in self.fatty_acids:
                fa_relatives = []
                for run in range(len(self.results[reaction])):
                    fa_relatives.append(self.results[reaction][run][0][fa])
                    self.fatty_acid_distribution[reaction][fa] = np.mean(fa_relatives)
                    self.fatty_acid_std[reaction][fa] = np.std(fa_relatives)


        # mean and std membrane lengths after each time step for each change
        for reaction in self.results:
            self.mean_time_lists[reaction] = {}
            self.std_time_lists[reaction] = {}
            for run in range(len(self.results[reaction])):
                self.plasma_membrane.append(self.results[reaction][run][1][0])
                self.secretory_vesicles.append(self.results[reaction][run][1][1])
                self.vacuoles.append(self.results[reaction][run][1][2])
                self.nucleus.append(self.results[reaction][run][1][3])
                self.peroxisomes.append(self.results[reaction][run][1][4])
                self.light_microsomes.append(self.results[reaction][run][1][5])
                self.inner_mit_membrane.append(self.results[reaction][run][1][6])
                self.outer_mit_membrane.append(self.results[reaction][run][1][7])
                self.lipid_droplets.append(self.results[reaction][run][1][8])

            z = 0
            for membrane in self.membrane_lists.values():
                self.mean_time_lists[reaction][self.membrane_lists.keys()[z]] = []
                self.std_time_lists[reaction][self.membrane_lists.keys()[z]] = []
                for i in range(len(self.time)):
                    membrane_time = []
                    for j in range(len(membrane)):
                        membrane_time.append(membrane[j][i])
                    mean_time = np.mean(membrane_time)
                    std_time = np.std(membrane_time)
                    self.mean_time_lists[reaction][self.membrane_lists.keys()[z]].append(mean_time)
                    self.std_time_lists[reaction][self.membrane_lists.keys()[z]].append(std_time)
                z += 1


        # mean and std membrane compositions for each change
        for reaction in self.results:
            self.membrane_comp_mean[reaction] = {}
            self.membrane_comp_std[reaction] = {}
            for membrane in self.membrane_lists.keys():
                self.membrane_comp_mean[reaction][membrane] = {}
                self.membrane_comp_std[reaction][membrane] = {}
                for lipid in self.lipids:
                    lipid_relatives = []
                    for runs in range(len(self.results[reaction])):
                        lipid_relatives.append(self.results[reaction][runs][2][membrane][lipid])
                    mean_lipid_relatives = np.mean(lipid_relatives)
                    std_lipid_relatives = np.std(lipid_relatives)
                    self.membrane_comp_mean[reaction][membrane][lipid] = mean_lipid_relatives
                    self.membrane_comp_std[reaction][membrane][lipid] = std_lipid_relatives

        pickle.dump(self.fatty_acid_distribution, open("./sens_analysis_start_10-12-15/fa_distr_mean.pkl", "wb"))
        pickle.dump(self.fatty_acid_std, open("./sens_analysis_start_10-12-15/fa_distr_std.pkl", "wb"))
        pickle.dump(self.mean_time_lists, open("./sens_analysis_start_10-12-15/membranes_length_mean.pkl", "wb"))
        pickle.dump(self.std_time_lists, open("./sens_analysis_start_10-12-15/membranes_length_std.pkl", "wb"))
        pickle.dump(self.membrane_comp_mean, open("./sens_analysis_start_10-12-15/membrane_comp_mean.pkl", "wb"))
        pickle.dump(self.membrane_comp_std, open("./sens_analysis_start_10-12-15/membrane_comp_std.pkl", "wb"))

    def heatmap_fatty_acids(self):
        # heatmap of fatty acid distributions: x = fatty acids, y = reaction (4 heatmaps: rates +/- 10%, probabilities +/- 10%)
        # , 'rate-10', 'probability+10', 'probability-10', 'compartment_weights+10', 'compartment_weights-10',
        # 'fatty_acid_weights+10', 'fatty_acid_weights-10']
        heatmaps_change = ['rate+10']
        for change in heatmaps_change:
            if 'rate' in change:
                column_labels = ['Wildtype']
                for par in self.rates:
                    column_labels.append(par)

            elif 'probability' in change:
                column_labels = ['Wildtype']
                for par in self.probabilities:
                    column_labels.append(par)

            elif 'compartment' in change:
                column_labels = self.membranes_list[:-1]

            elif 'fatty acid' in change:
                column_labels = self.fa

            row_labels = self.fatty_acids
            self.heatmap_fa = []

            i = 0
            for reaction in column_labels:
                if reaction == 'Wildtype':
                    self.heatmap_fa.append(i)
                    difference = []
                    for fa in row_labels:
                        diff = math.log(1, 2)
                        difference.append(diff)
                    self.heatmap_fa[i] = difference
                    i += 1

                else:
                    self.heatmap_fa.append(i)
                    difference = []
                    for fa in row_labels:
                        if self.fatty_acid_distribution[reaction + change][fa] > 0:
                            diff = math.log((self.fatty_acid_distribution[reaction + change][fa] / self.fatty_acid_distribution['Wildtype'][fa]), 2)
                        else:
                            diff = math.log(0.00001, 2)
                        difference.append(diff)
                    self.heatmap_fa[i] = difference
                    i += 1

            self.heatmap_fa = np.asarray(self.heatmap_fa)
            fig, ax = pyp.subplots()
            general_heatmap = ax.pcolor(self.heatmap_fa, cmap = pyp.cm.Blues)#, vmax = 0.2, vmin = -0.2)

            ax.set_xticks(np.arange(self.heatmap_fa.shape[1])+0.5, minor = False)
            ax.set_yticks(np.arange(self.heatmap_fa.shape[0])+0.5, minor = False)

            ax.invert_yaxis()
            ax.xaxis.tick_top()

            ax.set_xticklabels(row_labels, minor = False)
            ax.set_yticklabels(column_labels, minor = False)

            cbar = pyp.colorbar(general_heatmap)
            cbar.ax.set_ylabel('log2 change', rotation = 270, labelpad = 22)

            ax.set_title('Fatty acids distribution:' + change, y = 1.07)
            ax.axis('tight')

            pyp.show()

    def heatmap_membrane_comp(self):
        # heatmap of membrane compositions: x = lipids, y = reaction (4 heatmaps for each membrane, 4 x 9 heatmaps)
        heatmaps_change = ['rate+10']#, ' rate - 10 %', ' probability + 10 %', ' probability - 10 %', ' compartment weights + 10 %', ' compartment weights - 10 %', ' fatty acid weights + 10 %', ' fatty acid weights - 10 %']

        for membrane in self.membranes_list:
            for change in heatmaps_change:

                column_labels = ['Wildtype']
                if 'rate' in change:
                    for par in self.rates:
                        column_labels.append(par)
                elif 'probability' in change:
                    for par in self.probabilities:
                        column_labels.append(par)
                elif 'compartment' in change:
                    column_labels = self.membranes_list
                elif 'fatty acid' in change:
                    column_labels = self.fa

                self.heatmap_membrane = []
                i = 0
                for reaction in column_labels:
                    self.heatmap_membrane.append(i)
                    difference = []
                    for lipid in self.lipids:
                        if self.membrane_comp_mean['Wildtype'][membrane][lipid] == 0:
                            diff = math.log(1, 2)
                        else:
                            if reaction == 'Wildtype':
                                diff = math.log(1, 2)
                            elif self.membrane_comp_mean[reaction + change][membrane][lipid] > 0:
                                diff = math.log((self.membrane_comp_mean[reaction + change][membrane][lipid] / self.membrane_comp_mean['Wildtype'][membrane][lipid]), 2)
                            else:
                                diff = math.log(0.00001, 2)
                        difference.append(diff)
                    self.heatmap_membrane[i] = difference
                    i += 1

                self.heatmap_membrane = np.asarray(self.heatmap_membrane)
                fig, ax = pyp.subplots()
                general_heatmap = ax.pcolor(self.heatmap_membrane, cmap = pyp.cm.Blues)#, vmax = 1.6, vmin = -1.6)

                ax.set_xticks(np.arange(self.heatmap_membrane.shape[1])+0.5, minor = False)
                ax.set_yticks(np.arange(self.heatmap_membrane.shape[0])+0.5, minor = False)

                ax.invert_yaxis()
                ax.xaxis.tick_top()
                ax.set_xticklabels(self.lipids, minor = False)
                column_labels_short = [label.split(' ', 1)[0] for label in column_labels]
                ax.set_yticklabels(column_labels_short, minor = False)

                cbar = pyp.colorbar(general_heatmap)
                cbar.ax.set_ylabel('log2 change', rotation = 270, labelpad = 22)

                ax.set_title(membrane.replace('_', ' ').title() + ' composition:' + change, y = 1.1)
                ax.axis('tight')

                pyp.show()


    def heatmap_membrane_length(self):
        #heatmap of membrane lengths: x = membrane, y = reaction (4 heatmaps: rates +/- 10%, probabilities +/- 10%)
        heatmaps_change = ['rate+10']#, ' rate - 10 %', ' probability + 10 %', ' probability - 10 %', ' compartment weights + 10 %', ' compartment weights - 10 %', ' fatty acid weights + 10 %', ' fatty acid weights - 10 %']
        row_labels = self.membranes_list[:-1]

        for change in heatmaps_change:
            if 'rate' in change:
                column_labels = ['Wildtype']
                for par in self.rates:
                    column_labels.append(par)
            elif 'probability' in change:
                column_labels = ['Wildtype']
                for par in self.probabilities:
                    column_labels.append(par)
            elif 'compartment' in change:
                column_labels = self.membranes_list
            elif 'fatty acid' in change:
                column_labels = self.fa

            self.heatmap_mem_length = []

            i = 0
            for reaction in column_labels:
                self.heatmap_mem_length.append(i)
                difference = []
                for membrane in row_labels:
                    if reaction == 'Wildtype':
                        diff = math.log(1, 2)
                    elif self.mean_time_lists[reaction + change][membrane][-1] > 0:
                        diff = math.log((self.mean_time_lists[reaction + change][membrane][-1] / self.mean_time_lists['Wildtype'][membrane][-1]), 2)
                    elif self.mean_time_lists[reaction + change][membrane][-1] <= 0:
                        diff = math.log(0.00001, 2)
                    difference.append(diff)
                self.heatmap_mem_length[i] = difference
                i += 1

            self.heatmap_mem_length = np.asarray(self.heatmap_mem_length)
            fig, ax = pyp.subplots()
            general_heatmap = ax.pcolor(self.heatmap_mem_length, cmap = pyp.cm.Blues)#, vmax = 0.1, vmin = -0.1)

            ax.set_xticks(np.arange(self.heatmap_mem_length.shape[1])+0.5, minor = False)
            ax.set_yticks(np.arange(self.heatmap_mem_length.shape[0])+0.5, minor = False)

            ax.invert_yaxis()
            ax.xaxis.tick_top()

            ax.set_xticklabels(row_labels, minor = False)
            ax.set_yticklabels(column_labels, minor = False)

            cbar = pyp.colorbar(general_heatmap)
            cbar.ax.set_ylabel('log2 change', rotation=270, labelpad = 22)

            ax.set_title('Membrane lengths:' + change, y = 1.4)

            pyp.xticks(rotation = 70)
            ax.axis('tight')

            pyp.show()