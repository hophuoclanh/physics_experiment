import numpy as np
import matplotlib.pyplot as plt


# CONSTANTS
mu0 = 4 * np.pi * 1e-7

# PARAMETERS
I = 10 # Ampere
R = 0.6e-2/5  # Solenoid's radius
L = 10e-2 # Solenoid's length
N_turn = 50# Number of turns
N_phi = 60
N1 = N_phi * N_turn
phi_max = 2 * np.pi * N_turn
phi = np.linspace(0, phi_max, N1)
x_wire = np.linspace(0, L, N1)
y_wire = R * np.cos(phi)
z_wire = R * np.sin(phi)

Nx, Ny, Nz = 10, 10, 10
xmin, xmax, ymin, ymax, zmin, zmax = 0, L, -1e-2, 1e-2, -1e-2, 1e-2

# CALCULATION
N = len(x_wire) - 1
central_axis_x = np.linspace(xmin, xmax, Nx)
central_axis_y = 0
central_axis_z = 0

X_C, Y_C, Z_C = np.meshgrid(central_axis_x, central_axis_y, central_axis_z)


Bx, By, Bz = np.zeros_like(X_C, dtype=np.float64), np.zeros_like(Y_C, dtype=np.float64), np.zeros_like(Z_C, dtype=np.float64)


kB = mu0 / (4 * np.pi) * I
for i in range(N):
    dl_x = (x_wire[i + 1] - x_wire[i])
    dl_y = (y_wire[i + 1] - y_wire[i])
    dl_z = (z_wire[i + 1] - z_wire[i])

    rx = X_C - x_wire[i]
    ry = Y_C - y_wire[i]
    rz = Z_C - z_wire[i]

    r = np.sqrt(rx ** 2 + ry ** 2 + rz ** 2)
    r[r == 0] = 1e-20  # To avoid division by zero

    dBx = mu0 / (4 * np.pi) * (I * dl_y * rz - I * dl_z * ry) / r ** 3
    dBy = mu0 / (4 * np.pi) * (I * dl_z * rx - I * dl_x * rz) / r ** 3
    dBz = mu0 / (4 * np.pi) * (I * dl_x * ry - I * dl_y * rx) / r ** 3

    Bx += dBx
    By += dBy
    Bz += dBz

# print (Bx, By, Bz)
# Printing the magnetic field (Bx, By, Bz) at each position along the central axis
for position, field_strength in zip(X_C.flatten(), Bx.flatten()):
    print(f"Position: {position:.2f} m, Bx: {field_strength:.2e}")

# Calculate the number of turns per unit length
n = N_turn / L

# Calculate the magnetic field inside the solenoid
B_inside = mu0 * n * I


plt.figure(figsize=(12, 3))
plt.plot(X_C.flatten(), Bx.flatten(), color='r')
plt.axhline(B_inside, color='blue', label='B_inside')
plt.title('Bx Component')
plt.xlabel('Position along X-axis')
plt.ylabel('Bx')
plt.tight_layout()
plt.show()
