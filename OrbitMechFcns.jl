"""
THINGS TO LOOK UP

norm(), cross(), dot() ...
ode solvers
Make structure non-immutable
if you dont input all variables into function
Rotation Matricies ? made my own for now
Use functions from one script in another
int vs float (error in boolean statement
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

"""
Rotation Matricies
"""
function rotz(γ)
    rotmat = [cos(γ) -sin(γ) 0; sin(γ) cos(γ) 0; 0 0 1];
    return rotmat
end
function rotx(α)
    rotmat = [1 0 0;0 cos(α) -sin(α); 0 sin(α) cos(α)];
    return rotmat
end

"""
Structure for easy OE implementation
    define all in Radians always!!!
"""

struct OrbitalElementVec
    # Semi-major axis
    a::Float64
    # Eccentricity
    e::Float64
    # Inclination
    i::Float64
    # Right Ascension of Ascending Node
    Ω::Float64
    # Argument of Periapsis
    ω::Float64
    # True Anomaly
    ν::Float64
end

"""
Gravitational Force propagation with J2
Run in ode solver
Inputs: time (for solver)
        X - [6x1] position and velocity state vector
        μ - Gravitational Parameter
        rₑ - Radius of center body
        J₂ - Zonal perturbation force
Output: Ẋ - Derivative of state vector
"""
function GravitationalForce(t,X,μ=398600.4418,rₑ=6378.137,J₂=1.0826267e-3)

    Ẋ = zeros(6,1);
    r = X[1:3];
    v = X[4:6];
    accel = -μ*r/norm(r)^3;

    # Acceleration due to oblateness
    R = norm(r);
    r̂ = r/normr;
    k̂ = [0 0 1]';
    z = r[3];
    accJ2 = - (μ*J₂*rₑ*rₑ/2)*((6*z/(R^5))*k̂ + ((3/(R^4)) - (15*z*z/(R^6)))*r̂);

    Ẋ[1:3] = v;
    Ẋ[4:6] = accel+accJ2;

    return Ẋ
end



"""
OE to ECI
"""
function OE2ECI(oe::OrbitalElementVec,μ=398600.4418)
    P = oe.a*(1-oe.e^2);                # Semi-Latus Rectum
    r_mag = P/(1+oe.e*cos(oe.ν));       # Distance from Earth to orbiting body
    # R in perifocial coordinates [P Q W]'
    r_peri = [r_mag*cos(oe.ν); r_mag*sin(oe.ν); 0];
    v_peri = sqrt(μ/P)*[-sin(oe.ν); (oe.e+cos(oe.ν)); 0];
    if oe.i == 0 && oe.e != 0         # Equitorial and Elliptical
        R1 = 1;
        R2 = 1;
        R3 = rotz(oe.ω);
    elseif oe.e == 0 && oe.i != 0     # Circular and Inclined
        R1 = rotz(oe.Ω);
        R2 = rotx(oe.i);
        R3 = 1;
    elseif oe.i == 0 && oe.e == 0     # Equitorial and Circular
        R1 = 1;
        R2 = 1;
        R3 = 1;
    else                              # Not Circular or Inclined
        R1 = rotz(oe.Ω);
        R2 = rotx(oe.i);
        R3 = rotz(oe.ω);
    end

    R = R1*R2*R3;               # Full rotation matrix

    r_eci = R*r_peri;
    v_eci = R*v_peri;

    return r_eci, v_eci
end



"""
ECI to OE
"""
function  ECI2OE( R, V )
    # INPUTS
    #   R - I,J,K components of position
    #   V - I,J,K components of velocity
    # OUTPUTS
    #   oe - OrbitalElementVec
    # Function by
    #   Shawn Swist ~ 2018

    r = norm(R);
    v = norm(V);

    H = cross(R,V);
    h = norm(H);

    N = cross([0;0;1],H);
    n = norm(N);
    e_vec = 1/μ*((v^2-μ/r).*R-dot(R,V).*V);
    e = norm(e_vec);

    P = h^2/μ;
    a = -P/(e^2-1);

    # Orbital Inclination (always less than 180 deg)
    i = acos(H[3]/h);

    # Rignt Ascension of Ascending Node
    Ω = acos(N[1]/n);
    if N[2] < 0             # If Nj is greater than 0 Om is less than 180
        Ω = 2π- Ω;
    end

    # Argument of periapsis
    ω = acos(dot(N,e_vec)/(n*e));
    if e_vec[3] < 0         # If e(k) is greater than 0 w is less than 180
        ω = 2π - ω;
    end

    # True anomaly
    ν = acos(dot(e_vec,R)/(e*r));
    if dot(R,V) < 0         # If R dot V is greater than zero nu is less than 180
        ν = 2π - ν;
    end

    oe = OrbitalElementVec(a,e,i,Ω,ω,ν)
    return oe
end


"""
Anomaly swapping: ν -> E -> M -> E -> ν
"""

function anom2E(ν,e)
    E = acos((e + cos(ν))/(1 + e*cos(ν)));
    if ν > π
        E = 2π - E;
    end
    return E
end

function E2M(E,e)
    M = E - e*sin(E);
    return M
end

function M2E(M,e,tol=1e-8)
    if M == 0 || M == pi
        # Return the known solutions (trivial)
        E = M;
    else
        # Set up the problem based on an initial guess
        E0 = M;
        d = -(E0 - e*sin(E0) - M)/(1 - e*cos(E0));
        # Loop until the solution converges
        while abs(d) > tol
            E1 = E0 + d;
            d = -(E1 - e*sin(E1) - M)/(1 - e*cos(E1));
            E0 = E1;
        end
        E = E0;
    end
    return E
end

function E2anom(E,e)
    ν = acos((cos(E) - e)/(1 - e*cos(E)));
    if E > π
        ν = 2π - ν;
    end
end


# module OrbitMechFunctions
#     export OrbitalElementVec, OE2ECI, GravitationalForce
# end

oe1 = OrbitalElementVec(7500,0.015,0,pi/6,0,pi)

R,V = OE2ECI(oe1)
