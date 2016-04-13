import matplotlib.pyplot as mat 

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
