import matplotlib.pyplot as plt
import numpy as np

def plot_steady_solution():
    NX = 21
    NY = 21
    LX = 1
    LY = 1
    dx = LX / (NX - 1)
    dy = LY / (NY - 1)
    dt = 0.01
    X = np.linspace(0, LX, NX)
    Y = np.linspace(0, LY, NY)
    #X = np.array([0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0]); NX = X.shape[0]
    #Y = np.array([0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0]); NY = Y.shape[0]
    XMG, YMG = np.meshgrid(X, Y, indexing="xy")
    
    u = np.loadtxt("u.txt", delimiter=",")
    v = np.loadtxt("v.txt", delimiter=",")
    p = np.loadtxt("p.txt", delimiter=",")

    plt.contourf(X, Y, p)
    plt.colorbar()
    plt.streamplot(XMG, YMG, u, v, color="black")
    plt.quiver(XMG, YMG, u, v)
    plt.title("Lid Driven Cavity")
    plt.savefig("lid_driven_cavity.png")
    plt.close()

    breakpoint()


if __name__ == "__main__":
    plot_steady_solution()