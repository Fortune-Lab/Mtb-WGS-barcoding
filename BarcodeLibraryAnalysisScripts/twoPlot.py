import matplotlib.pyplot as plt
import numpy as np

# Sample data
x = np.linspace(0, 10, 100)
y1 = np.sin(x)
y2 = x**2

# Create the figure and the first axes
fig, ax1 = plt.subplots()

# Plot the first graph on the first axes
ax1.plot(x, y1, color='red', label='sin(x)')
ax1.set_xlabel('x')
ax1.set_ylabel('sin(x)', color='red')
ax1.tick_params(axis='y', labelcolor='red')

# Create a second axes sharing the same x-axis
ax2 = ax1.twinx()

# Plot the second graph on the second axes
ax2.plot(x, y2, color='blue', label='x^2')
ax2.set_ylabel('x^2', color='blue')
ax2.tick_params(axis='y', labelcolor='blue')

# Add a legend
#fig.legend(loc="upper left", bbox_to_anchor=(0.1, 0.95))

# Display the plot
plt.savefig("test.png")
