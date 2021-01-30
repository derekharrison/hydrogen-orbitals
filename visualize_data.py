import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import *

grid_data_file = 'grid_data.txt'
nr, ntheta, nphi, max_pd = np.genfromtxt(grid_data_file, unpack=True)

nr = int(nr)
ntheta = int(ntheta)
nphi = int(nphi)

num_data_file = 'data' + '.txt'
X, Z, Psi_sq = np.genfromtxt(num_data_file, unpack=True)

X_2d = np.zeros((nr, nphi))
Z_2d = np.zeros((nr, nphi))
Psi_sq_2d = np.zeros((nr, nphi))

X_2d_m = np.zeros((nr, nphi))
Z_2d_m = np.zeros((nr, nphi))
Psi_sq_2d_m = np.zeros((nr, nphi))

num = 0
for i in range(0, nr):
    for j in range(0, nphi):
        X_2d[i][j] = X[num]
        Z_2d[i][j] = Z[num]
        Psi_sq_2d[i][j] = Psi_sq[num]
        num = num + 1 

for i in range(0, nr):
    for j in range(0, nphi):
        X_2d_m[i][j] = X[num]
        Z_2d_m[i][j] = Z[num]
        Psi_sq_2d_m[i][j] = Psi_sq[num]
        num = num + 1

X_t = np.zeros((2*nr, nphi))
Z_t = np.zeros((2*nr, nphi))
Psi_p = np.zeros((2*nr, nphi))

for i in range(0, nr):
    for j in range(0, nphi):
        X_t[i][j] = X_2d[i][j]
        Z_t[i][j] = Z_2d[i][j]
        Psi_p[i][j] = Psi_sq_2d[i][j]

for i in range(0, nr):
    for j in range(0, nphi):
        X_t[i+nr][j] = X_2d_m[i][j]
        Z_t[i+nr][j] = Z_2d_m[i][j]
        Psi_p[i+nr][j] = Psi_sq_2d_m[i][j]
        
cont = plt.contourf(X_t, Z_t, Psi_p, linspace(0, max_pd, 300))
plt.colorbar()
plt.xlabel("X")
plt.ylabel("Z")

plt.show()
