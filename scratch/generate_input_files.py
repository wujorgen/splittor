import numpy as np

def lid_driven_cavity():
    NX = 51
    NY = 51
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

    mu = 0.1
    rho = 1.0

    lid_velocity = 2.0

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
                if jdx == NY - 1 and idx == NX//2:
                    # pressure poisson needs a fixed pressure value for at least one point, so here it is
                    PBCTYPE = 1
                    PVAL = 0
                f.write(f"{X[idx]:.5f}  {Y[jdx]:.5f}  {0:.5f}  {VBCTYPE:5d}  {PBCTYPE:5d}  {UVAL:.5f}  {VVAL:.5f}  {0:.5f}  {PVAL:.5f}\n")

    with open("probleminformation.txt", "w") as f:
        f.write(f"mu  {mu}\n")
        f.write(f"rho  {rho}\n")
        f.write(f"dt  {dt}\n")
        f.write("mode  steady\n")


if __name__ == "__main__":
    lid_driven_cavity()
