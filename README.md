# Splittor
This code solves the incompressible Navier-Stokes equations using the projection method, also known as Chorin's projection method.

# TODO
- minimum viable product
    - explicit scheme for intermediate velocity calculation
    - simple bcs where entire side is bc
        - cavity flow
        - pressure pipe flow
- investigate iterating on pressure corrections in chorin projection step
- if implicit scheme, figure out what picard iteration is
- more complex BCs
    - example: cavity where jet enters left side, cavity is entirely open on right side. no slip conditions top and bottom, and around the left on the jet. I think we could just enforce no left-right velocity gradient at the opening on the right.
- turbulence model
- unstructured mesh
