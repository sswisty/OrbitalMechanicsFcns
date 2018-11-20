include("OrbitMechFcns.jl")

oe1 = OrbitalElementVec(7500,0.015,0,pi/6,0,pi)

R,V = OE2ECI(oe1)
