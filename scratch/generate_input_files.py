import numpy as np

def lid_driven_cavity():
    NX = 11
    NY = 11
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
                    PBCTYPE = 1
                    UVAL = lid_velocity
                elif jdx == 0 or idx == 0 or idx == NX - 1:
                    VBCTYPE = 1
                    PBCTYPE = 0  # do not specify - no grad default at domain edge
                f.write(f"{X[idx]:.2f}  {Y[jdx]:.2f}  {0:.2f}  {VBCTYPE:2d}  {PBCTYPE:2d}  {UVAL:.2f}  {VVAL:.2f}  {0:.2f}  {PVAL:.2f}\n")

    with open("probleminformation.txt", "w") as f:
        f.write(f"mu  {mu}\n")
        f.write(f"rho  {rho}\n")
        f.write(f"dt  {dt}\n")
        f.write("mode  steady\n")


if __name__ == "__main__":
    lid_driven_cavity()
