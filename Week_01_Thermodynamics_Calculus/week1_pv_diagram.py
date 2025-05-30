import numpy as np
import matplotlib.pyplot as plt

def isothermal(V, n=1, R=8.314, T=300):
    return n * R * T / V

def adiabatic(V, n=1, R=8.314, T0=300, gamma=1.4):
    return (n * R * T0) * (V[0] / V)**gamma

V = np.linspace(0.001, 0.01, 100)
p1 = isothermal(V)
p2 = adiabatic(V)

plt.plot(V, p1, label='Isothermal')
plt.plot(V, p2, label='Adiabatic')
plt.xlabel('Volume [m^3]')
plt.ylabel('Pressure [Pa]')
plt.title('p-V Diagram')
plt.legend()
plt.grid(True)
plt.show()
