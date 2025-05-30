# 🚀 Engineering Relaunch Study Plan (Month 1)

A guided month-long syllabus to rebuild engineering mastery with integrated math and simulation. Designed for self-paced, practical relearning.

---

## 📅 Week 1: Thermodynamics + Calculus Refresher
**Goal:** Reconnect with physical laws through calculus

### 🔧 Topics
- 1st & 2nd Law of Thermodynamics
- Ideal gas law, energy conservation
- Derivatives: pressure-volume relationships
- Integrals: work done in expansion/compression
- Visualizing p–v–T surfaces with Python

### 📆 Tasks
- [ ] Review thermodynamic state properties: p, V, T, U, H, S
- [ ] Derive work done for isothermal and adiabatic processes (ideal gas)
- [ ] Code a Python script to plot:
  - Isothermal process in p-V space
  - Adiabatic process in p-V space
- [ ] Practice partial derivatives: ∂U/∂T, ∂p/∂V, ∂H/∂p
- [ ] Solve integrals involving work and internal energy
- [ ] Plot T-s diagram for a simple Rankine-like cycle (qualitative)

### 🖊️ Python Starter Script Outline
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

### 📚 Resources
- Moran & Shapiro – *Fundamentals of Engineering Thermodynamics* (to be replaced if needed)
- MIT OCW Thermodynamics
- Paul's Online Math Notes (Calculus I & II)

---

## 📅 Week 2: Probability, Statistics & Data Foundations
**Goal:** Build statistical intuition for engineering simulations

### 🔧 Topics
- Distributions: normal, uniform, binomial
- Expected value, variance, standard deviation
- Error propagation and confidence intervals
- Monte Carlo simulations (basic)
- Histograms and probability visualizations

### 📚 Resources
- OpenIntro Statistics
- Khan Academy – Statistics & Probability
- Python: NumPy, Pandas, Matplotlib

---

## 📅 Week 3: Fluid Mechanics + Vector Calculus
**Goal:** Solidify flow intuition and conservation law grounding

### 🔧 Topics
- Streamlines, pathlines, velocity fields
- Bernoulli’s principle
- Gradient, divergence, curl
- Applying vector calculus to continuity
- Streamline plotting in Python

### 📚 Resources
- White – *Fluid Mechanics*
- NPTEL Fluid Mechanics (IIT video series)
- 3Blue1Brown – *Vector Calculus* visuals
- **NEW**: *Fluid Dynamics and Heat Transfer of Turbomachinery* (for foundational support leading into Week 5)

---

## 📅 Week 4: PDEs + Numerical Methods
**Goal:** Connect math theory to simulation techniques

### 🔧 Topics
- Heat, wave, Laplace equations (intro)
- Finite difference method (1D heat conduction)
- Euler and Runge-Kutta methods
- Stability & accuracy (Courant condition)

### 📚 Resources
- Zill – *Differential Equations with Boundary Value Problems*
- Chapra – *Applied Numerical Methods*
- MIT OCW Numerical Methods
- Python: SciPy, Jupyter, Matplotlib

---

## 🤔 After Month 1
- **Week 5**: Turbomachinery — Axial/radial stages, velocity triangles, efficiencies  
  ☑ Book: *Fluid Dynamics and Heat Transfer of Turbomachinery*

- **Week 6**: Material Resistance — Stress, strain, failure modes  
  ☑ Book suggestions to be confirmed

- **Week 7**: Combustion Theory & High-Speed Flow  
  ☑ Book: *Hypersonic and High-Temperature Gas Dynamics* (Anderson)  
  ☑ Begin focus on **supersonic combustion** and **aerobreathing engines** (SCRAMJET intro)

- **Week 8+**: Advanced Thermodynamics + Multiphysics Coupling  
  ☑ Bridging thermo with CFD/simulations and optimization techniques

Stay consistent. Note what you understand, log code + math in GitHub, and revisit difficult topics iteratively.

> "Relearning with experience is not repetition — it's engineering evolution."

---

To be continued with **Week 2 structure and exercises...**
