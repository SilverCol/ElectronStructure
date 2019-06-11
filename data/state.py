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
xline, = ax.plot(domain, (3 / (2 * np.pi**2))**(1/3) * np.power(np.absolute(state), 2/3))
ax.legend((sline, pline, xline), ('$|\\psi (r)|^2$', '$U(r)$', '$V_{xc}(r)$'))
ax.set_xlabel('$r$')

plt.show()
