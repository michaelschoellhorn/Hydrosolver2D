import numpy as np 
import matplotlib.pyplot as plt


def loadData(fileName):

    # Get grid dimensions
    Nx, Ny, Nghost, Cfl, deltaX, deltaY = np.loadtxt(fileName, skiprows=1, max_rows=1, unpack=True)
    Nx = int(Nx)
    Ny = int(Ny)
    Nghost = int(Nghost)

    # Load Q's
    Q1 = np.loadtxt(fileName, skiprows=3, max_rows=Nx+2*Nghost)
    Q2x = np.loadtxt(fileName, skiprows=4+Nx+2*Nghost, max_rows=Nx+2*Nghost)
    Q2y = np.loadtxt(fileName, skiprows=5+2*Nx+4*Nghost, max_rows=Nx+2*Nghost)
    Q3 = np.loadtxt(fileName, skiprows=6+3*Nx+6*Nghost, max_rows=Nx+2*Nghost)

    return Nx, Ny, Nghost, deltaX, deltaY, Q1, Q2x, Q2y, Q3


Nx, Ny, Nghost, deltaX, deltay, Q1, Q2x, Q2y, Q3 = loadData('Data.txt')
oneDQ1 = Q1[:, 5]
oneDQ2x = Q2x[:, 5]
oneDQ3 = Q3[:, 5]
plt.plot(oneDQ1, scaley=True)
plt.show()
plt.plot(oneDQ2x, scaley=True)
plt.show()
plt.plot(oneDQ3, scaley=True)
plt.show()