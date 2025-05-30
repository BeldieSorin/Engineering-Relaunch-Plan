import numpy as np
import matplotlib.pyplot as plt

# Parameters
L = 1.0        # Length of the rod
T = 0.5        # Total time
nx = 20        # Number of spatial steps
nt = 100       # Number of time steps
alpha = 0.01   # Thermal diffusivity

dx = L / (nx - 1)
dt = T / nt
r = alpha * dt / dx**2

# Initial and boundary conditions
u = np.zeros(nx)
u[int(nx/2)] = 100  # Initial heat pulse at center

# Time stepping
u_record = [u.copy()]
for _ in range(nt):
    u_new = u.copy()
    for i in range(1, nx - 1):
        u_new[i] = u[i] + r * (u[i+1] - 2*u[i] + u[i-1])
    u = u_new
    u_record.append(u.copy())

# Plot result
plt.figure(figsize=(8, 4))
for i in range(0, nt+1, 20):
    plt.plot(np.linspace(0, L, nx), u_record[i], label=f't={i*dt:.2f}s')
plt.title("1D Heat Conduction (Explicit FDM)")
plt.xlabel("Position along rod")
plt.ylabel("Temperature")
plt.legend()
plt.grid(True)
plt.show()
