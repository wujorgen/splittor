# Splittor
This code solves the incompressible Navier-Stokes equations using the projection method, also known as Chorin's projection method.

# TODO
- minimum viable product
    - explicit, semi-implicit, and implicit scheme(s) for intermediate velocity calculation
    - simple bcs where entire side is bc
        - cavity flow
        - pressure pipe flow
- investigate iterating on pressure corrections in chorin projection step
- more complex BCs
    - example: cavity where jet enters left side, cavity is entirely open on right side. no slip conditions top and bottom, and around the left on the jet. I think we could just enforce no left-right velocity gradient at the opening on the right.
    - allow arbitrary velocity boundary conditions anywhere, pressure can only be direchlet or neumann (default zero) at sides
- turbulence model
- unstructured mesh
