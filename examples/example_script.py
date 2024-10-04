"""
Simple Plot
===========

This script demonstrates how to create a simple plot using matplotlib.

"""

import matplotlib.pyplot as plt
import numpy as np

# Generate data
x = np.linspace(0, 10, 100)
y = np.sin(x)

# Create plot
plt.plot(x, y)
plt.title("Sine Wave")
plt.xlabel("x")
plt.ylabel("sin(x)")

# Display plot
plt.show()
