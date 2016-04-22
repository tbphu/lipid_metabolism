import model

# create files to store time courses
m = model.Model()
# membrane sizes
for membrane in m.membranes_state:
    with open('./' + membrane + '.txt', 'wb') as lipid_file:
        line = ""
        timerange = range(0, 7200)
        for i in timerange:
            line += str(i) + ", "
        line += "\n"
        lipid_file.write(line)

# membrane compositions
for membrane in m.comp_ratio_dict:
    with open('./' + membrane + '_comp.txt', 'wb') as comp_file:
        line = ""
        for lip in m.MEMBRANE_LIPID_NAMES:
            line += lip + ", "
        line += "\n"
        comp_file.write(line)

# FA distribution
with open('./fa_distribution.txt', 'wb') as fa_file:
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
    m.RATES['ergosterol_synthase'] = 0
    i += 1
    print i
    r, mem, s = m.run(7200)
    # save membrane size time courses
    for membrane in mem:
        with open('./' + membrane + '.txt', 'a') as lipid_file:
            line = ""
            for tp in mem[membrane]:
                line += str(tp) + ", "
            line += "\n"
            lipid_file.write(line)
    # save membrane compositions at t = 7200
    for membrane in m.comp_ratio_dict:
        with open('./' + membrane + '_comp.txt', 'a') as comp_file:
            line = ""
            for lipid in m.MEMBRANE_LIPID_NAMES:
                line += str(m.comp_ratio_dict[membrane][lipid]) + ", "
            line += '\n'
            comp_file.write(line)
    # save fa composition at t = 7200
    with open('./fa_distribution.txt', 'a') as fa_file:
        line = ""
        for fa in m.saturation_composition_total:
            line += str(m.saturation_composition_total[fa])
        line += '\n'
        fa_file.write(line)
    m = None
    r = None
    mem = None
    s = None