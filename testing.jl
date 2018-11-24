include("OrbitMechFcns.jl")

oe₁ = OrbitalElementVec(7500,0.015,1*pi/180,pi/6,0,pi/4)

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
plot1 = plot(xdata,ydata,zdata,xlim=(-8000,8000), ylim=(-8000,8000), zlim=(-8000,8000),
       title = "Orbit via ODEsolver")

# surface(xdata,ydata,zdata)

t = collect(0:5:86400) # one day, every 5 seconds
epoch = UT12MJD(11,21,2018,0,0,0)
Reci,Veci = MeanMotionProp(oe₁,epoch,t/86400)

xdata = zeros(length(Reci),1) # There should be a better way to extract this ..
ydata = zeros(length(Reci),1)
zdata = zeros(length(Reci),1)

for i = 1:length(Reci)
    xdata[i] = Reci[i][1]
    ydata[i] = Reci[i][2]
    zdata[i] = Reci[i][3]
end

plot2 = plot(xdata,ydata,zdata,xlim=(-8000,8000), ylim=(-8000,8000), zlim=(-8000,8000),
       title = "Orbit via Mean Motion Prop")

display(plot1)
display(plot2)


# Load Earth topo for ground track
EP = EarthGroundPlot()
GMST = MJD2GMST(epoch.+t/86400)
Recef,Vecef = ECI2ECEF(Reci,Veci,GMST)
ϕ,λ,h = ECEF2GEO(Recef)

plot!(λ,ϕ)
