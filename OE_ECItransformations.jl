
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
