import model

i = 0
for i in range(10):
    m = model.Model()
    i += 1
    print i
    r, mem, s = m.run(0)
    for membrane in m.membranes_state:
        with open('./' + membrane + '.txt', 'a') as lipid_file:
            line = str(len(m.membranes_state[membrane])) + '\n'
            lipid_file.write(line)
    for component in m.components_state:
        with open('./' + component + '.txt', 'a') as comp_file:
            line = str(len(m.components_state[component])) + '\n'
            comp_file.write(line)