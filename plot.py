import matplotlib.pyplot as plt
import numpy as np

# Read data from the file
data = np.loadtxt("Continuity.txt")

# Extract x-coordinates and electron concentrations
x_coordinates = data[0, :]
electron_concentration = data[1:, :]

# Plot the data
plt.figure(figsize=(10, 6))
for i in range(electron_concentration.shape[0]):
    if i==0:
        plt.plot(x_coordinates, electron_concentration[i, :], label=f'Timestamp {i+1}')
    if i==99:
        plt.plot(x_coordinates, electron_concentration[i, :], label=f'Timestamp {i+1}')

plt.xlabel('x (m)')
plt.ylabel('Electron Concentration')
plt.title('Electron Concentration vs. x for Multiple Timestamps')
plt.legend()
plt.grid(True)
plt.show()
