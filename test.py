from H1_timeseries import *
from H2_timeseries import *

x, y = yearly_H1_mean_global()

x1, y1 = yearly_H2_mean_global()

# Fit a linear model y = mx + b
m, b = np.polyfit(x, y, 1)  # 1 means linear fit

# Generate fitted line
x_fit = np.linspace(min(x), max(x), 100)  # Smooth line
y_fit = m * x_fit + b

print(m)

# Fit a linear model y = mx + b
m1, b1 = np.polyfit(x1, y1, 1)  # 1 means linear fit

# Generate fitted line
x1_fit = np.linspace(min(x1), max(x1), 100)  # Smooth line
y1_fit = m1 * x1_fit + b1


print(m1)

plt.plot(x, y)
plt.plot(x1, y1)
plt.plot(x_fit, y_fit)
plt.plot(x1_fit, y1_fit)
plt.grid(True)
plt.show()
