"""
THINGS TO LOOK UP

norm(), cross(), dot() ...
ode solvers
Make structure non-immutable
Rotation Matricies ? made my own for now
plotting some of this stuff (EarthPlot, ground tracks ...)


FUNCTIONS INCLUDED
    * Orbital Element Structure
    * OE2ECI
    * GravitationalForce
    * Rotation Matricies
    * ECI2OE
    * anom2E, E2anom
    * E2M, M2E
"""

# Read in all functions from their respective files
include("anomalytransformations.jl")
include("OEstruct.jl")
include("RotMat.jl")
include("GravitationalForcePropagator.jl")
include("OE_ECItransformations.jl")
