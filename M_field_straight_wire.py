import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# CONSTANTS
mu0 = 4 * np.pi * 1e-7

# Parameters for the finite-length straight wire
L_wire = 40  # Length of the wire (meters)
I_wire = 1  # Current in the wire (Amperes)

# Grid for field calculation
Nx, Ny, Nz = 10, 6, 6
xmin = -1
xmax = 1
ymin = -1
ymax = 1
zmin = -1
zmax = 1

# Define the points along a line parallel to the wire at a fixed y-distance
# Points along the line parallel to the wire
x_line_straight = np.linspace(xmin, xmax, Nx)
y_line_straight = 0
z_line_straight = 0

# Points along the line parallel to the wire
line_distance = 2  # Distance of the line from the wire
x_points = np.linspace(-L_wire/2, L_wire/2, Nx)  # X-coordinates of the points

# Initialize the magnetic field array
B_magnetic_field_line = np.zeros(Nx)

# Calculation using the Biot-Savart Law for a finite-length wire
for i, x in enumerate(x_points):
    for wire_x in np.linspace(-L_wire/2, L_wire/2, 1000):  # Discretize the wire
        dl = L_wire / 1000  # Length of each wire segment
        r = np.sqrt((x - wire_x) ** 2 + line_distance ** 2)  # Distance from the wire segment to the observation point

        # Magnetic field contribution from each wire segment (only y-component due to symmetry)
        dB = mu0 / (4 * np.pi) * (I_wire * dl) / (r ** 2)
        B_magnetic_field_line[i] += dB * line_distance / r  # Biot-Savart Law for straight wire


# Plotting the magnetic field along the line parallel to the wire
plt.figure(figsize=(12, 6))
plt.plot(x_points, B_magnetic_field_line, label='Calculated Magnetic Field', color='blue')
plt.axhline(mu0 * I_wire / (2 * np.pi * line_distance), color='red', linestyle='--', label='Idealized Field')
plt.title('Magnetic Field Along a Line Parallel to a Finite-Length Wire')
plt.xlabel('Position along X-axis (m)')
plt.ylabel('Magnetic Field (T)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
