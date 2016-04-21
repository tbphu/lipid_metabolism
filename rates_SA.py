import model
import sys

def SA(rate):
    changes_rates = {'glycerol_3_p_synthesis': 1,
                      'inositol_synthesis': 1,
                      'ceramide_synthesis': 1,
                      'acetyl_coa_synthase': 65,
                      'acyl_synthase': 45,
                      'PA_synthesis': 2,
                      'CDP_DG_synthase': 2,
                      'TAG_synthesis': 3,
                      'TAG_lipase': 2,
                      'DAG_kinase': 4,
                      'PS_synthase': 2,
                      'PI_synthase': 1,
                      'PE_synthase': 2,
                      'PC_synthase': 1,
                      'CL_synthase': 1,
                      'ergosterol_synthase': 3,
                      'sterylester_synthase': 3,
                      'sphingolipid_synthase': 1}

    changes_rates_neg = {'glycerol_3_p_synthesis': -1,
                      'inositol_synthesis': -1,
                      'ceramide_synthesis': -1,
                      'acetyl_coa_synthase': -65,
                      'acyl_synthase': -45,
                      'PA_synthesis': -2,
                      'CDP_DG_synthase': -2,
                      'TAG_synthesis': -3,
                      'TAG_lipase': -2,
                      'DAG_kinase': -4,
                      'PS_synthase': -2,
                      'PI_synthase': -1,
                      'PE_synthase': -2,
                      'PC_synthase': -1,
                      'CL_synthase': -1,
                      'ergosterol_synthase': 3,
                      'sterylester_synthase': 3,
                      'sphingolipid_synthase': 1}
    # create files to store time courses
    m = model.Model()
    # membrane sizes
    for membrane in m.membranes_state:
        with open('./' + membrane + "_" + rate + '+' + str(changes_rates_neg[rate]) + '.txt', 'wb') as lipid_file:
            line = ""
            timerange = range(0, 7200)
            for i in timerange:
                line += str(i) + ", "
            line += "\n"
            lipid_file.write(line)

    # membrane compositions
    for membrane in m.comp_ratio_dict:
        with open('./' + membrane + "_" + rate + '+' + str(changes_rates_neg[rate]) + '_comp.txt', 'wb') as comp_file:
            line = ""
            for lip in m.MEMBRANE_LIPID_NAMES:
                line += lip + ", "
            line += "\n"
            comp_file.write(line)

    # FA distribution
    with open('./fa_distribution_' + rate + '+' + str(changes_rates_neg[rate]) + '.txt', 'wb') as fa_file:
        line = ""
        for fa in m.saturation_composition_total:
            line += fa + ", "
        line += "\n"
        fa_file.write(line)
    m = None

    # make simulations and save result
    i = 0
    for i in range(1000):
        m = model.Model()
        m.RATES[rate] += changes_rates_neg[rate]
        i += 1
        print i
        r, mem, s = m.run(7200)
        # save membrane size time courses
        for membrane in mem:
            with open('./' + membrane + "_" + rate + '+' + str(changes_rates_neg[rate]) + '.txt', 'a') as lipid_file:
                line = ""
                for tp in mem[membrane]:
                    line += str(tp) + ", "
                line += "\n"
                lipid_file.write(line)
        # save membrane compositions at t = 7200
        for membrane in m.comp_ratio_dict:
            with open('./' + membrane + "_" + rate + '+' + str(changes_rates_neg[rate]) + '_comp.txt', 'a') as comp_file:
                line = ""
                for lipid in m.MEMBRANE_LIPID_NAMES:
                    line += str(m.comp_ratio_dict[membrane][lipid]) + ", "
                line += '\n'
                comp_file.write(line)
        # save fa composition at t = 7200
        with open('./fa_distribution_' + rate + '+' + str(changes_rates_neg[rate]) + '.txt', 'a') as fa_file:
            line = ""
            for fa in m.saturation_composition_total:
                line += str(m.saturation_composition_total[fa])
            line += '\n'
            fa_file.write(line)
        m = None
        r = None
        mem = None
        s = None


if __name__ == '__main__':
    SA(sys.argv[1])
