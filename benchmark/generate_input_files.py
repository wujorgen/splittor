import numpy as np
import sys

def lid_driven_cavity(resolution:int=65, Re:float=10):
    NX = resolution
    NY = resolution
    LX = 1
    LY = 1
    dx = LX / (NX - 1)
    dy = LY / (NY - 1)
    dt = 0.01
    X = np.linspace(0, LX, NX)
    Y = np.linspace(0, LY, NY)
    #X = np.array([0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0]); NX = X.shape[0]
    #Y = np.array([0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0]); NY = Y.shape[0]
    XMG, YMG = np.meshgrid(X, Y, indexing="ij")

    rho = 1
    mu  = 1 / Re

    lid_velocity = 1.0

    with open("grid.txt", "w") as f:
        f.write(f"NX  {NX}  NY  {NY}  NZ  {1}  D  {2}\n")
        f.write("X  Y  Z  VBCTYPE  PBCTYPE  UVAL  VVAL  WVAL  PVAL\n")
        for jdx in range(NY):
            for idx in range(NX):
                ij = jdx * NX + idx
                VBCTYPE = 0
                PBCTYPE = 0
                UVAL = 0
                VVAL = 0
                PVAL = 0
                if jdx == NY - 1:
                    VBCTYPE = 1
                    PBCTYPE = 0  # do not specify - no grad default at domain edge
                    UVAL = lid_velocity
                elif jdx == 0 or idx == 0 or idx == NX - 1:
                    VBCTYPE = 1
                    PBCTYPE = 0  # do not specify - no grad default at domain edge
                if jdx == 1 and idx == 1:
                    # pressure poisson needs a fixed pressure value for at least one point, so here it is
                    # it should be noted the output of the solver is affected if the pin is at the top or bottom.
                    # top gave better results for re=100,400, but pinning a single point seems to not work for 1000.
                    # anyways i'm moving the pin to the bottom left interior point
                    PBCTYPE = 1
                    PVAL = 0
                f.write(f"{X[idx]:.5f}  {Y[jdx]:.5f}  {0:.5f}  {VBCTYPE:5d}  {PBCTYPE:5d}  {UVAL:.5f}  {VVAL:.5f}  {0:.5f}  {PVAL:.5f}\n")

    with open("probleminformation.txt", "w") as f:
        f.write(f"mu  {mu}\n")
        f.write(f"rho  {rho}\n")
        f.write(f"dt  {dt}\n")
        f.write("mode  steady\n")
        f.write("steady_state_max_steps  5000")
    
    np.save("XMG.npy", XMG)
    np.save("YMG.npy", YMG)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("Usage: python path/to/this/script.py RESOLUTION REYNOLDS_NUMBER")
        sys.exit()
    if len(sys.argv) > 1:
        resolution = int(sys.argv[1])
    else:
        resolution = 65
    if len(sys.argv) > 2:
        Re = float(sys.argv[2])
    else:
        Re = 10
    print(f"Generating input file for lid driven cavity with resolution of {resolution} x {resolution} and reynolds number of {Re}")
    lid_driven_cavity(resolution, Re)
