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
