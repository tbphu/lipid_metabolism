import random
import matplotlib.pyplot as pyp

# x axis molecule numbers
x = range(100000)
x = [x/100. for x in x]

# stochastic model reaction
def rate(h, km, s):
    result = 0
    # only discrete molecule numbers
    sz = round(s)
    # convert vmax to discrete integer
    hint = int(round(h))
    # actual reaction rate
    for i in range(hint):
        rand = random.random()
        if rand < (sz/(km + sz)):
            result += 1
    return result

# parameters
vmax = 100.
km = 30.

# compute y axis molecules / sec
y_determ = [vmax * (xs)/(km + xs) for xs in x]
y_stoch = [rate(vmax, km, xs) for xs in x]

# plot result
pyp.plot(x,y_stoch, 'b', label='stochastic', linewidth=1)
pyp.plot(x,y_determ, 'r', label='deterministic', linewidth=3)
pyp.xlabel('number of molecules [#]')
pyp.ylabel('number of molecules per s [#/s]')
pyp.title('Compare deterministic and stochastic approach')
pyp.show()







