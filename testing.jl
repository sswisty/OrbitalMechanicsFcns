include("OrbitMechFcns.jl")

oe₁ = OrbitalElementVec(7500,0.015,1*pi/180,pi/6,0,pi/4)

R,V = OE2ECI(oe₁)

oe₂,ang = ECI2OE(R,V)
