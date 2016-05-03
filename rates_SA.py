import model
import sys

# script to perform sensitivity analysis - called with argument
# @author: Jens Hahn

def SA(rate):
    changes_rates_pos = {'glycerol_3_p_synthesis': 1,
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
    mod = model.Model()
    # membrane sizes
    for membrane in mod.membranes_state:
        with open('./' + membrane + "_" + rate + '+' + str(changes_rates_neg[rate]) + '.txt', 'wb') as lipid_file:
            line = ""
            timerange = range(0, 7200)
            for i in timerange:
                line += str(i) + ", "
            line += "\n"
            lipid_file.write(line)

    # membrane compositions
    for membrane in mod.comp_ratio_dict:
        with open('./' + membrane + "_" + rate + '+' + str(changes_rates_neg[rate]) + '_comp.txt', 'wb') as comp_file:
            line = ""
            for lip in mod.MEMBRANE_LIPID_NAMES:
                line += lip + ", "
            line += "\n"
            comp_file.write(line)

    # FA distribution
    with open('./fa_distribution_' + rate + '+' + str(changes_rates_neg[rate]) + '.txt', 'wb') as fa_file:
        line = ""
        for fa in mod.saturation_composition_total:
            line += fa + ", "
        line += "\n"
        fa_file.write(line)
    mod = None

    # make simulations and save result
    i = 0
    for i in range(1000):
        mod = model.Model()
        mod.RATES[rate] += changes_rates_neg[rate]
        i += 1
        print i
        sat, mem, comp = mod.run(7200)
        # save membrane size time courses
        for membrane in mem:
            with open('./' + membrane + "_" + rate + '+' + str(changes_rates_neg[rate]) + '.txt', 'a') as lipid_file:
                line = ""
                for tp in mem[membrane]:
                    line += str(tp) + ", "
                line += "\n"
                lipid_file.write(line)
        # save membrane compositions at t = 7200
        for membrane in mod.comp_ratio_dict:
            with open('./' + membrane + "_" + rate + '+' + str(changes_rates_neg[rate]) + '_comp.txt', 'a') as comp_file:
                line = ""
                for lipid in mod.MEMBRANE_LIPID_NAMES:
                    line += str(mod.comp_ratio_dict[membrane][lipid]) + ", "
                line += '\n'
                comp_file.write(line)
        # save fa composition at t = 7200
        with open('./fa_distribution_' + rate + '+' + str(changes_rates_neg[rate]) + '.txt', 'a') as fa_file:
            line = ""
            for fa in mod.saturation_composition_total:
                line += str(mod.saturation_composition_total[fa]) + ", "
            line += '\n'
            fa_file.write(line)
        # garbage collector
        mod = None
        sat = None
        mem = None
        comp = None


if __name__ == '__main__':
    SA(sys.argv[1])
