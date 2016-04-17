import model

for i in range(1000):

    m = model.Model()
    # test run: 5 sec
    try:
        r, mem, s = m.run()
    except:
        r, mem, s = m.run()

    for membrane in mem:
        with open('./' + membrane + '.txt', 'a') as lipid_file:
            line = str(mem[membrane]).strip('[').strip(']') + '\n'
            lipid_file.write(line)
