import model


def run_model():
    m = model.Model()
    # test run: 5 sec
    try:
        r, mem, s = m.run()
        return r, mem, s
    except:
        print 'fault'
        run_model()

i = 0
for i in range(1000):
    i += 1
    print i
    r, mem, s = run_model()
    for membrane in mem:
        with open('./' + membrane + '.txt', 'a') as lipid_file:
            line = str(mem[membrane]).strip('[').strip(']') + '\n'
            lipid_file.write(line)





