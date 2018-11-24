"""
Find when a satellite is visible from a specific ground location
    do i want lat/long coords in degrees or radians???
"""


function GroundRange(Recef,ϕground,λground,rₑ=6378.137)

    Rgs = rₑ*[cos.(ϕground).*cos.(λground);
              cos.(ϕground).*sin.(λground);
              sin.(ϕground)]

    # Ground Station and Rotation Matrix
    ê = [-sin(GClong);cos(GClong);0]
    n̂ = [-sin(GClat)*cos(GClong);-sin(GClat)*sin(GClong);cos(GClat)]
    û = [cos(GClat)*cos(GClong);cosd(GClat)*sin(GClong);sin(GClat)]

    Rxyzenu = [ê n̂ û]'

    Rsat = Rxyzenu.*(r_ecef.-Rgs)
