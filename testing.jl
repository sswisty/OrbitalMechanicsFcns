include("OrbitMechFcns.jl")

oe₁ = OrbitalElementVec(7500,0.015,5*pi/180,pi/6,0,pi/4)

R,V = OE2ECI(oe₁)

oe₂,ang = ECI2OE(R,V)


ICs = [R;V]
tspan = (0.0,9000.0)
params = [398600.4418;6378.137;1.0826267e-3]
prob = ODEProblem(GravitationalForce,ICs,tspan,params)
sol = solve(prob,saveat = 2)

xdata = zeros(length(sol.u),1) # There should be a better way to extract this ..
ydata = zeros(length(sol.u),1)
zdata = zeros(length(sol.u),1)

for i = 1:length(sol.u)
    xdata[i] = sol.u[i][1]
    ydata[i] = sol.u[i][2]
    zdata[i] = sol.u[i][3]
end
plot(xdata,ydata,zdata,xlim=(-8000,8000), ylim=(-8000,8000), zlim=(-8000,8000),
       title = "Orbit")

# surface(xdata,ydata,zdata)
