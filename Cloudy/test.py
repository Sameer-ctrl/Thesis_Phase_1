import numpy as np
import matplotlib.pyplot as plt
import corner

# Generate some random data for demonstration
np.random.seed(0)
data = np.random.randn(100, 3)

# Create the corner plot with scatter plots only
figure = corner.corner(data[:, :2], labels=['Parameter 1', 'Parameter 2'], plot_contours=False)

axes_fig=figure.axes

# Calculate the median values of the posterior distributions for Parameter 1 and Parameter 2
median_parameter1 = np.median(data[:, 0])
median_parameter2 = np.median(data[:, 1])

# Plot vertical line at the median value of Parameter 1
# plt.axvline(median_parameter1, color='red', linestyle='--')
axes_fig[2].vlines(median_parameter1,-2,10, color='red', linestyle='--')
# Plot horizontal line at the median value of Parameter 2
axes_fig[2].hlines(median_parameter2,-4,4, color='red', linestyle='--')
axes_fig[2].scatter(0,0,marker='s',s=100)
plt.show()
