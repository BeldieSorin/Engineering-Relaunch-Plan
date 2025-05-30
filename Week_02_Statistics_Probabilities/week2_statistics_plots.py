import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binom, norm

# Binomial distribution
x_binom = np.arange(0, 21)
pmf = binom.pmf(x_binom, n=20, p=0.5)

# Normal distribution
x_norm = np.linspace(-4, 4, 1000)
pdf = norm.pdf(x_norm, loc=0, scale=1)

# Plot both
plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.stem(x_binom, pmf, basefmt=" ")
plt.title("Binomial PMF (n=20, p=0.5)")

plt.subplot(1, 2, 2)
plt.plot(x_norm, pdf)
plt.title("Normal PDF (mean=0, std=1)")
plt.grid(True)

plt.tight_layout()
plt.show()
