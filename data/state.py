from matplotlib import pyplot as plt
import numpy as np

R = 15
state = np.fromfile('state' + str(R) + '.bin')[::100]
pot = np.fromfile('potential' + str(R) + '.bin')[::100]
domain = np.linspace(R / len(state), R + R / len(state), len(state))

plt.rcParams.update({'font.size': 15})
fig = plt.figure()
ax = fig.subplots()

ax.grid()
sline, = ax.plot(domain, state**2)
pline, = ax.plot(domain, pot)
ax.legend((sline, pline), ('$|\\psi (r)|^2$', 'U(r)'))
ax.set_xlabel('$r$')

plt.show()
