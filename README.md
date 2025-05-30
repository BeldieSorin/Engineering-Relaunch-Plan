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

### 📚 Resources
- Zill – *Differential Equations with Boundary Value Problems*
- Chapra – *Applied Numerical Methods*
- MIT OCW Numerical Methods
- Python: SciPy, Jupyter, Matplotlib

---

## 🗅️ Week 5: Turbomachinery Fundamentals
**Goal:** Build a foundational grasp of axial and radial turbomachinery for propulsion and energy applications

### 🔧 Topics
- Types of turbomachines: compressors vs. turbines, radial vs. axial
- Velocity triangles (inlet/outlet, absolute/relative frames)
- Blade angles, reaction ratio, degree of reaction
- Euler’s turbomachinery equation
- Efficiency measures (isentropic vs. polytropic)
- Stage-by-stage analysis
- Intro to cascade flow

### 🖖️ Tasks
- [ ] Sketch velocity triangles for an axial stage
- [ ] Derive Euler turbine equation and apply to a turbine/compressor stage
- [ ] Calculate stage efficiency for given inlet/outlet conditions
- [ ] Visualize pressure and velocity changes across blades
- [ ] Interpret common performance maps (compressor, turbine)

### 💊 Python/Math Starter
- Placeholder: Python plotting of velocity triangles and simplified turbine work calculation (to be added)

### 📚 Resources
- *Fluid Dynamics and Heat Transfer of Turbomachinery*
- Dixon & Hall – *Fluid Mechanics and Thermodynamics of Turbomachinery*
- MIT OCW Gas Turbines Lectures

---

## 🗅️ Week 6: Material Resistance & Structural Mechanics
**Goal:** Understand stress, strain, and deformation in engineering components

### 🔧 Topics
- Normal and shear stress/strain
- Hooke’s Law and elastic modulus
- Poisson’s ratio and material properties
- Mohr’s Circle (intro)
- Bending, torsion, axial loading
- Beam deflection and stress concentration
- Failure criteria (yield, fracture)

### 🖖️ Tasks
- [ ] Compute axial stress/strain for simple bars
- [ ] Analyze torsion in circular shafts
- [ ] Draw shear and bending moment diagrams
- [ ] Estimate max bending stress in beams
- [ ] Intro sketch of Mohr’s Circle for 2D stress state
- [ ] Review material safety factors and strength assumptions

### 📚 Resources
- Hibbeler – *Mechanics of Materials*
- Gere & Goodno – *Mechanics of Materials*
- MIT OCW Solid Mechanics / Statics
- YouTube: Michel van Biezen – Mechanics of Materials series

---

## 🗅️ Week 7: Combustion Theory & High-Speed Flow
**Goal:** Explore fundamentals of reacting flows and supersonic combustion relevant to airbreathing engines

### 🔧 Topics
- Chemical kinetics: reaction rates, Arrhenius law
- Energy release and flame temperature
- Premixed and diffusion flames
- Detonation vs. deflagration
- Combustion in nozzles and ducts
- High-speed reacting flow: compressibility effects
- Introduction to **SCRAMJET** propulsion
- Oblique shocks and heat addition in supersonic flow

### 🖖️ Tasks
- [ ] Calculate adiabatic flame temperature for hydrogen-air
- [ ] Sketch Mach number change with heat addition (Rayleigh flow)
- [ ] Analyze a SCRAMJET schematic for flow path characteristics
- [ ] Plot oblique shock angle vs. Mach number for various deflection angles
- [ ] Compute simple laminar flame speed from kinetic model

### 💊 Python Snippet Starter
```python
import numpy as np
import matplotlib.pyplot as plt

def oblique_shock_angle(M, theta_deg, gamma=1.4):
    theta = np.radians(theta_deg)
    beta_range = np.radians(np.linspace(theta_deg, 90, 500))
    f = lambda beta: 2 * np.tan(theta) * (M**2 * np.sin(beta)**2 - 1) / \
                     (M**2 * (gamma + np.cos(2 * beta)) + 2) - np.tan(beta - theta)
    diff = [abs(f(b)) for b in beta_range]
    best_beta = beta_range[np.argmin(diff)]
    return np.degrees(best_beta)

theta_vals = np.linspace(1, 30, 100)
shock_angles = [oblique_shock_angle(3.0, th) for th in theta_vals]

plt.plot(theta_vals, shock_angles)
plt.xlabel("Deflection Angle (°)")
plt.ylabel("Shock Angle (°)")
plt.title("Oblique Shock Angle vs Deflection (M=3)")
plt.grid(True)
plt.show()
```

### 📚 Resources
- Anderson – *Hypersonic and High-Temperature Gas Dynamics*
- Turns – *An Introduction to Combustion*
- MIT OCW: Unified Engineering (Propulsion Lectures)
- NASA Glenn SCRAMJET resources

---

## 🤔 After Month 1
- **Week 8+**: Advanced Thermodynamics + Multiphysics Coupling  
  ☑ Bridging thermo with CFD/simulations and optimization techniques

Stay consistent. Note what you understand, log code + math in GitHub, and revisit difficult topics iteratively.

> "Relearning with experience is not repetition — it's engineering evolution."

---

To be continued...
