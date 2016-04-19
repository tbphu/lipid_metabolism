import model

# create files
m = model.Model()
for membrane in m.membranes_state:
    with open('./' + membrane + '.txt', 'wb') as lipid_file:
        line = "Time, "
        timerange = range(0, 7200)
        for i in timerange:
            line += str(i) + ", "
        line += "\n"
        lipid_file.write(line)

for membrane in m.comp_ratio_dict:
    with open('./' + membrane + '_comp.txt', 'wb') as comp_file:
        line = "Lipid:, "
        for lip in m.MEMBRANE_LIPID_NAMES:
            line += lip + ", "
        line += "\n"
        comp_file.write(line)
m = None

i = 0
for i in range(3):
    m = model.Model()
    i += 1
    print i
    r, mem, s = m.run(0)
    for membrane in m.membranes_state:
        with open('./' + membrane + '.txt', 'a') as lipid_file:
            line = str(len(m.membranes_state[membrane])) + '\n'
            lipid_file.write(line)

    for membrane in m.comp_ratio_dict:
        with open('./' + membrane + '_comp.txt', 'a') as comp_file:
            line = ""
            for lipid in m.comp_ratio_dict[mem]:
                line += m.comp_ratio_dict[mem][lipid] + ", "
            line += '\n'
            comp_file.write(line)
    m = None