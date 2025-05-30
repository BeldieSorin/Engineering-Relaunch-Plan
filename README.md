# 🚀 Engineering Relaunch Study Plan (Month 1)

A guided month-long syllabus to rebuild engineering mastery with integrated math and simulation. Designed for self-paced, practical relearning.

---

## 🗅️ Week 1: Thermodynamics + Calculus Refresher
**Goal:** Reconnect with physical laws through calculus

### 🔧 Topics
- 1st & 2nd Law of Thermodynamics
- Ideal gas law, energy conservation
- Derivatives: pressure-volume relationships
- Integrals: work done in expansion/compression
- Visualizing p–v–T surfaces with Python

### 🖖️ Tasks
- [ ] Review thermodynamic state properties: p, V, T, U, H, S
- [ ] Derive work done for isothermal and adiabatic processes (ideal gas)
- [ ] Code a Python script to plot:
  - Isothermal process in p-V space
  - Adiabatic process in p-V space
- [ ] Practice partial derivatives: ∂U/∂T, ∂p/∂V, ∂H/∂p
- [ ] Solve integrals involving work and internal energy
- [ ] Plot T-s diagram for a simple Rankine-like cycle (qualitative)

### 💊 Python Starter Script Outline
```python
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
```

🔗 [Go to Week 1 folder](./Week_01_Thermodynamics_Calculus)

### 📚 Resources
- Moran & Shapiro – *Fundamentals of Engineering Thermodynamics* (to be replaced if needed)
- MIT OCW Thermodynamics
- Paul's Online Math Notes (Calculus I & II)

---

## 🗅️ Week 2: Probability, Statistics & Data Foundations
**Goal:** Build statistical intuition for engineering simulations

### 🔧 Topics
- Probability axioms and rules
- Discrete distributions: binomial, Poisson
- Continuous distributions: normal, uniform
- Mean, variance, skewness, kurtosis
- Central Limit Theorem and its importance
- Confidence intervals and significance
- Monte Carlo simulation (first use)
- Visualizing uncertainty with Python (histograms, scatter plots)

### 🖖️ Tasks
- [ ] Review probability axioms and basic rules (addition/multiplication)
- [ ] Code and plot:
  - A binomial distribution (n=20, p=0.5)
  - A normal distribution (mean=0, std=1)
- [ ] Use NumPy to simulate 10,000 dice rolls and plot histogram
- [ ] Perform a basic Monte Carlo estimate of π (random sampling)
- [ ] Calculate confidence intervals for sample means
- [ ] Generate and interpret a correlation matrix in Pandas

### 💊 Python Snippet Starter
```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binom, norm

x_binom = np.arange(0, 21)
pmf = binom.pmf(x_binom, n=20, p=0.5)

x_norm = np.linspace(-4, 4, 1000)
pdf = norm.pdf(x_norm, loc=0, scale=1)

plt.figure(figsize=(12,5))
plt.subplot(1, 2, 1)
plt.stem(x_binom, pmf, basefmt=" ")
plt.title("Binomial PMF (n=20, p=0.5)")

plt.subplot(1, 2, 2)
plt.plot(x_norm, pdf)
plt.title("Normal PDF (mean=0, std=1)")
plt.grid(True)
plt.show()
```

🔗 [Go to Week 2 folder](./Week_02_Statistics_Probability)

### 📚 Resources
- OpenIntro Statistics
- Khan Academy – Statistics & Probability
- Python: NumPy, Pandas, Matplotlib, SciPy

---

## 🗅️ Week 3: Fluid Mechanics + Vector Calculus
**Goal:** Solidify flow intuition and conservation law grounding

### 🔧 Topics
- Continuity, momentum, and energy equations (differential forms)
- Velocity fields, streamlines, streaklines
- Bernoulli’s equation and control volume analysis
- Vector calculus: gradient, divergence, curl
- Stream function and potential function (2D)
- Vorticity and circulation
- Navier–Stokes equation (introductory form)
- Plotting streamlines and flow fields in Python

### 🖖️ Tasks
- [ ] Derive and interpret the continuity equation
- [ ] Visualize 2D flow field (vector plot) in Python
- [ ] Compute divergence and curl for a sample vector field
- [ ] Sketch streamlines using analytical or numeric approach
- [ ] Apply Bernoulli equation to real or ideal flows
- [ ] Explore stream function in irrotational flow

### 💊 Python Snippet Starter
```python
import numpy as np
import matplotlib.pyplot as plt

X, Y = np.meshgrid(np.linspace(-2, 2, 20), np.linspace(-2, 2, 20))
U = -Y / (X**2 + Y**2 + 0.1)
V = X / (X**2 + Y**2 + 0.1)

plt.figure(figsize=(6,6))
plt.streamplot(X, Y, U, V, color=np.sqrt(U**2 + V**2), linewidth=1)
plt.title("Streamlines of Sample Flow Field")
plt.xlabel("x")
plt.ylabel("y")
plt.axis("equal")
plt.grid(True)
plt.show()
```

🔗 [Go to Week 3 folder](./Week_03_Fluid_Vector)

### 📚 Resources
- White – *Fluid Mechanics*
- NPTEL Fluid Mechanics (IIT video series)
- 3Blue1Brown – *Vector Calculus* visuals
- *Fluid Dynamics and Heat Transfer of Turbomachinery* (for transition to Week 5)

---

## 🗅️ Week 4: PDEs + Numerical Methods
**Goal:** Connect math theory to simulation techniques

### 🔧 Topics
- Heat, wave, Laplace equations (intro)
- Finite difference method (1D heat conduction)
- Euler and Runge-Kutta methods
- Stability & accuracy (Courant condition)

### 💊 Python Snippet Starter
```python
import numpy as np
import matplotlib.pyplot as plt

L = 1.0
nx = 50
dx = L / (nx - 1)
alpha = 1e-4
dt = 0.1
total_time = 10
nt = int(total_time / dt)

u = np.zeros(nx)
u[int(nx/4):int(nx/2)] = 100

for n in range(nt):
    un = u.copy()
    for i in range(1, nx-1):
        u[i] = un[i] + alpha * dt / dx**2 * (un[i+1] - 2*un[i] + un[i-1])

plt.plot(np.linspace(0, L, nx), u)
plt.title("Heat Conduction Profile at Final Time")
plt.xlabel("Length")
plt.ylabel("Temperature")
plt.grid(True)
plt.show()
```

🔗 [Go to Week 4 folder](./Week_04_PDEs_Numerical)

### 📚 Resources
- Zill – *Differential Equations with Boundary Value Problems*
- Chapra – *Applied Numerical Methods*
- MIT OCW Numerical Methods
- Python: SciPy, Jupyter, Matplotlib

---

## 🗅️ Week 5: Turbomachinery Fundamentals
**Goal:** Understand energy transfer in rotating machines

### 🔧 Topics
- Energy exchange in compressors and turbines
- Euler’s turbomachinery equation
- Velocity triangles (inlet/outlet)
- Axial vs radial machines
- Blade loading and losses
- Stage efficiency, reaction ratio

### 🖖️ Tasks
- [ ] Draw velocity triangles for axial turbine stage
- [ ] Derive work input/output using Euler’s equation
- [ ] Compare axial vs radial turbine efficiency
- [ ] Plot velocity triangles using Python
- [ ] Solve numerical example for pressure ratio and efficiency

### 💊 Python Snippet Starter
```python
import matplotlib.pyplot as plt

U = 200  # Blade speed [m/s]
V1 = 300  # Absolute velocity at inlet
alpha1 = 30  # Inlet angle [degrees]

# Velocity triangle (inlet)
plt.figure()
plt.quiver(0, 0, V1 * np.cos(np.radians(alpha1)), V1 * np.sin(np.radians(alpha1)), angles='xy', scale_units='xy', scale=1, label='V1')
plt.quiver(0, 0, U, 0, angles='xy', scale_units='xy', scale=1, label='U')
plt.xlim(0, 350)
plt.ylim(0, 200)
plt.gca().set_aspect('equal', adjustable='box')
plt.grid(True)
plt.title('Turbomachinery Inlet Triangle')
plt.legend()
plt.show()
```

🔗 [Go to Week 5 folder](./Week_05_Turbomachinery)

### 📚 Resources
- Dixon & Hall – *Fluid Mechanics and Thermodynamics of Turbomachinery*
- Open-access tutorials on velocity triangles
- Khan Academy – Rotational energy basics

---

## 🗅️ Week 6: Material Resistance & Structural Mechanics
**Goal:** Rebuild fundamentals of stress analysis and structure behavior

### 🔧 Topics
- Stress, strain, Young’s modulus
- Axial/torsional/bending loads
- Beam deflection formulas
- Mohr’s Circle basics
- Factor of safety and yield criteria

### 🖖️ Tasks
- [ ] Plot axial stress/strain graphs
- [ ] Calculate moment of inertia for simple sections
- [ ] Draw shear and moment diagrams (by hand or Python)
- [ ] Construct Mohr’s Circle for a given stress state
- [ ] Analyze cantilever beam deflection with boundary conditions

### 💊 Python Snippet Starter
```python
import numpy as np
import matplotlib.pyplot as plt

L = 1.0  # length of beam [m]
E = 200e9  # Young's modulus [Pa]
I = 8e-6  # moment of inertia [m^4]
P = 1000  # Load at free end [N]

x = np.linspace(0, L, 100)
y = -P * x**2 * (3*L - x) / (6 * E * I)

plt.plot(x, y * 1e3)
plt.title('Deflection of Cantilever Beam')
plt.xlabel('Position [m]')
plt.ylabel('Deflection [mm]')
plt.grid(True)
plt.show()
```

🔗 [Go to Week 6 folder](./Week_06_Material_Resistance)

### 📚 Resources
- Hibbeler – *Mechanics of Materials*
- MIT OCW – Solid Mechanics
- Python: Matplotlib for beam diagrams

---

## 🗅️ Week 7: Combustion Theory & High-Speed Flow
**Goal:** Understand chemical energy release and its interaction with compressible flow

### 🔧 Topics
- Combustion reactions and stoichiometry
- Enthalpy of formation and heat release
- Equivalence ratio and flammability
- Flame speed, ignition delay
- Shock waves and expansion fans
- Supersonic inlets, combustion chambers
- SCRAMJET cycle intro

### 🖖️ Tasks
- [ ] Balance combustion equations (methane, hydrogen)
- [ ] Estimate adiabatic flame temperature
- [ ] Derive 1D shock jump conditions (Rankine–Hugoniot)
- [ ] Visualize Mach number vs pressure/temp across shock
- [ ] Analyze basic SCRAMJET component efficiency

### 💊 Python Snippet Starter
```python
import numpy as np
import matplotlib.pyplot as plt

M1 = np.linspace(1.1, 5.0, 100)
P2_P1 = (2.8 * M1**2 - 0.8) / 2.8

plt.plot(M1, P2_P1)
plt.title('Normal Shock Pressure Ratio vs Mach Number')
plt.xlabel('Mach number (M1)')
plt.ylabel('P2 / P1')
plt.grid(True)
plt.show()
```

🔗 [Go to Week 7 folder](./Week_07_Combustion_Supersonic)

### 📚 Resources
- Turns – *An Introduction to Combustion*
- Anderson – *Modern Compressible Flow*
- *Hypersonic and High-Temperature Gas Dynamics* (real copy)

---

## 🗅️ Week 8: Multiphysics Simulation & Coupled Systems
**Goal:** Bridge thermofluid theory to simulation using real-world tools

### 🔧 Topics
- What is multiphysics?
- Coupling heat + fluid flow + stress
- CFD vs FEM vs FVM methods
- Meshing strategies
- Turbulence modeling basics
- Case: heat exchanger or nozzle under stress

### 🖖️ Tasks
- [ ] Install OpenFOAM or SimScale or Ansys Student
- [ ] Simulate 2D heat conduction with boundary conditions
- [ ] Visualize temperature and velocity fields
- [ ] Explore solver convergence and stability
- [ ] Read a mesh file (.msh or .vtk)

🔗 [Go to Week 8 folder](./Week_08_Multiphysics_CFD)

### 📚 Resources
- COMSOL Multiphysics tutorials
- OpenFOAM beginner guide
- *Numerical Heat Transfer and Fluid Flow* – Patankar

---
