#=
Tranisit

All angles are in radians except for latitudes and longitudes

=#
include("OrbitMechFcns.jl") # Hopefully this will be a package soon!
using LsqFit

# Define constants
μ = 398600.441      #
rₑ = 6378.137       # [km] Earth radius
J₂ = 1.0826267e-3   # Oblatness
# Epoch, start/stop time and step size
# ⏰ = UT12MJD(11,22,2018,0,0,0)
⏰ = UT12MJD(10,7,2018,0,0,0)
# tstart = .0215*86400
tstart = 0 # This should always be zero! Otherwise you miss the propagation from epoch to tstart!
# tstop = 2*60*60 # two hours
tstop = .03*86400
dt = 2 # [sec] step size

# Define the orbital elements of the satellite
a = rₑ+500          # 500 km altitude orbit
e = 0.01            # Near circular orbit
i = deg2rad(89)     # Near polar orbit
Ω = deg2rad(90)
ω = 0
ν = 0
oeᵢ = OrbitalElementVec(a,e,i,Ω,ω,ν)
R,V = OE2ECI(oeᵢ)

# Simulate the orbit with ODE propagator
ICs = [R;V]
tspan = (tstart,tstop)
params = [μ;rₑ;J₂]
prob = ODEProblem(GravitationalForce,ICs,tspan,params)
ECIsol = solve(prob,saveat = dt)

# Extract the position and velocity data
# Subscript i for inertial frame
Rᵢ, Vᵢ = [],[]
for vect in ECIsol.u
    push!(Rᵢ,vect[1:3])
    push!(Vᵢ,vect[4:6])
end

# # Making sure the data looks right!
# xdata = zeros(length(ECIsol.u),1) # There should be a better way to extract this ..
# ydata = zeros(length(ECIsol.u),1)
# zdata = zeros(length(ECIsol.u),1)
# for i = 1:length(ECIsol.u)
#     xdata[i] = ECIsol.u[i][1]
#     ydata[i] = ECIsol.u[i][2]
#     zdata[i] = ECIsol.u[i][3]
# end
# plot1 = plot(xdata,ydata,zdata,xlim=(-8000,8000), ylim=(-8000,8000), zlim=(-8000,8000),
#        title = "Orbit via ODEsolver")


# Turn the ECI postitions into ECEF positions
# Subscript x for earth-fixed frame
MJD = ⏰ .+ ECIsol.t./86400
θ = MJD2GMST(MJD)
Rₓ,Vₓ = ECI2ECEF(Rᵢ,Vᵢ,θ)

# # Plot ground track for fun ...
# EP = EarthGroundPlot()
ϕ,λ,h = ECEF2GEO(Rₓ)
# plot!(λ,ϕ)


# Define the ground station
ϕground = 37.42662
λground = -122.173355

Rsat = GroundRange(Rₓ,ϕground,λground)

Az,Alt = SatAzAlt(Rsat)
viz = findall(Alt.>0)

# Plot GS and where sat is visible
EP = EarthGroundPlot()
plot!([λground λground],[ϕground ϕground],markershape=:star5,markersize=3)
plot!(λ[viz],ϕ[viz])

"""
Here are the functions converting the relative position into the relative
velocities and doppler curve
"""

# Combine the ECEF position and velocities into one vector
X = []
for i = 1:length(Rₓ)
    push!(X,hcat(Rₓ[i]',Vₓ[i]'))
end
Xviz = X[viz]
XX = vcat(Xviz...)

function GetRelVel(X,P)
    H = []
    for vect in X
        val = (vect[1:3]-P)'*vect[4:6]/sqrt((vect[1:3]-P)'*(vect[1:3]-P))
        push!(H,val)
    end
    return H
end
function LSQfun(R,P)
    x = [r[1] for r in R]
    y = [r[2] for r in R]
    z = [r[3] for r in R]
    u = [r[4] for r in R]
    v = [r[5] for r in R]
    w = [r[6] for r in R]

    top = ((x.-P[1]).*u .+ (y.-P[2]).*v .+ (z.-P[3]).*w)
    bottom = sqrt.((x.-P[1]).^2 .+ (y.-P[2]).^2 .+ (z.-P[3]).^2)
    h = top./bottom
    return h
end
#
# @. M3(x,p) = ((x[1]-p[1])*x[4]+(x[2]-p[2])*x[5]+(x[3]-p[3])*x[6])/sqrt((x[1]-p[1])^2+(x[2]-p[2])^2+(x[3]-p[3])^2)
# @. model(x,p) = dot((x[1:3]-p),x[4:6])/sqrt(dot((x[1:3]-p),(x[1:3]-p)))
# @. M2(X,P) = diag((X[1:3]-P)'*X[4:6])./sqrt(diag((X[1:3]-P)'*(X[1:3]-P)))
Rgs = rₑ*[cosd(ϕground)*cosd(λground);cosd(ϕground)*sind(λground);sind(ϕground)]

testdata = GetRelVel(X,Rgs)
testdata2 = LSQfun(Xviz,Rgs)
p0 = [rₑ,0,0]
fit = curve_fit(LSQfun, Xviz, testdata2, p0)

estPos = fit.param
