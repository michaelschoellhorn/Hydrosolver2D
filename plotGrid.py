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

def plot1D(oneDArray):
    fig, ax = plt.subplots()
    ax.plot(oneDArray)
    plt.show()


def plotGrid(twoDArray, colorbarTitle, title, XMax, YMax):
    fig, ax = plt.subplots()
    im = ax.imshow(twoDArray, cmap=plt.cm.viridis, origin='lower', aspect='auto', extent=(0, XMax, 0, YMax))
    fig.colorbar(im, label=colorbarTitle, ax=ax)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title(title)
    plt.show()


Nx, Ny, Nghost, deltaX, deltay, Q1, Q2x, Q2y, Q3 = loadData('Data.txt')
oneDQ1 = Q1[:, 5]
oneDQ2x = Q2x[:, 5]
oneDQ3 = Q3[:, 5]
#plot1D(oneDQ1)
plotGrid(Q1[Nghost:Nx+Nghost, Nghost:Ny+Nghost], r'$\rho$', 'mass density', 4, 4)
plotGrid(Q2x[Nghost:Nx+Nghost, Nghost:Ny+Nghost]/Q1[Nghost:Nx+Nghost, Nghost:Ny+Nghost], r'$u_{x}$', 'velocity in x direction', 4, 4)
plotGrid(Q2y[Nghost:Nx+Nghost, Nghost:Ny+Nghost]/Q1[Nghost:Nx+Nghost, Nghost:Ny+Nghost], r'$u_{y}$', 'velocity in y direction', 4, 4)
plotGrid(Q3[Nghost:Nx+Nghost, Nghost:Ny+Nghost], r'$\epsilon$', 'energy density', 4, 4)
