import model

# create files to store time courses
mod = model.Model()
# membrane sizes
for membrane in mod.membranes_state:
    with open('./' + membrane + '_erg_test.txt', 'wb') as lipid_file:
        line = ""
        timerange = range(0, 7200)
        for i in timerange:
            line += str(i) + ", "
        line += "\n"
        lipid_file.write(line)

# membrane compositions
for membrane in mod.comp_ratio_dict:
    with open('./' + membrane + '_comp_erg_test.txt', 'wb') as comp_file:
        line = ""
        for lip in mod.MEMBRANE_LIPID_NAMES:
            line += lip + ", "
        line += "\n"
        comp_file.write(line)

# FA distribution
with open('./fa_distribution_erg_test.txt', 'wb') as fa_file:
    line = ""
    for fa in mod.saturation_composition_total:
        line += fa + ", "
    line += "\n"
    fa_file.write(line)
mod = None

# make 1000 simulations and save result
i = 0
for i in range(1000):
    mod = model.Model()
    # ---------------------------------------------------
    # actual test condition
    mod.RATES['ergosterol_synthase'] = 0
    # ---------------------------------------------------
    i += 1
    print i
    sat, mem, comp = mod.run(7200)
    # save membrane size time courses
    for membrane in mem:
        with open('./' + membrane + '_erg_test.txt', 'a') as lipid_file:
            line = ""
            for tp in mem[membrane]:
                line += str(tp) + ", "
            line += "\n"
            lipid_file.write(line)
    # save membrane compositions at t = 7200
    for membrane in mod.comp_ratio_dict:
        with open('./' + membrane + '_comp_erg_test.txt', 'a') as comp_file:
            line = ""
            for lipid in mod.MEMBRANE_LIPID_NAMES:
                line += str(mod.comp_ratio_dict[membrane][lipid]) + ", "
            line += '\n'
            comp_file.write(line)
    # save fa composition at t = 7200
    with open('./fa_distribution_erg_test.txt', 'a') as fa_file:
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