import model

# create files
m = model.Model()
for membrane in m.membranes_state:
    with open('./' + membrane + '.txt', 'wb') as lipid_file:
        line = ""
        timerange = range(0, 7200)
        for i in timerange:
            line += str(i) + ", "
        line += "\n"
        lipid_file.write(line)

for membrane in m.comp_ratio_dict:
    with open('./' + membrane + '_comp.txt', 'wb') as comp_file:
        line = ""
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
    r, mem, s = m.run(7200)
    for membrane in mem:
        with open('./' + membrane + '.txt', 'a') as lipid_file:
            line = ""
            for tp in mem[membrane]:
                line += str(tp) + ", "
            line += "\n"
            lipid_file.write(line)

    for membrane in m.comp_ratio_dict:
        with open('./' + membrane + '_comp.txt', 'a') as comp_file:
            line = ""
            for lipid in m.MEMBRANE_LIPID_NAMES:
                line += str(m.comp_ratio_dict[membrane][lipid]) + ", "
            line += '\n'
            comp_file.write(line)
    m = None
    r = None
    mem = None
    s = None