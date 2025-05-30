import numpy as np
import matplotlib.pyplot as plt

# Define meshgrid for X and Y axes
X, Y = np.meshgrid(np.linspace(-2, 2, 20), np.linspace(-2, 2, 20))

# Define the vector field
U = -Y / (X**2 + Y**2 + 0.1)
V = X / (X**2 + Y**2 + 0.1)

# Plot the streamlines
plt.figure(figsize=(6, 6))
plt.streamplot(X, Y, U, V, color=np.sqrt(U**2 + V**2), linewidth=1)
plt.title("Streamlines of Sample Flow Field")
plt.xlabel("x")
plt.ylabel("y")
plt.axis("equal")
plt.grid(True)
plt.show()
