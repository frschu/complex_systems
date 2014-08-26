lorentz_oscillator
==================

Analysis of the Lorentz model 

This program intergrates the loretz model of convection numerically. The model is given by the following 
differential equations:

    dx/dt = sigma * (-x + y)
    dy/dt = -x * z + r * x - y
    dz/dt = x * y - b * z

with
    x = convection
    y = temperature difference for ascending and descending particles
    z = deviation to linear temperature profile
    sigma = Prandtl number = fluid viscosity / thermal diffusivity
    r = imposed T difference
    b = geometrical factor
sigma, r, b > 0

Linear analysis yields three fixed points:

    R0 = (0, 0, 0)
    eigenvalues: l_1 = -b, l_2,3 = -0.5 * ((sigma + 1) +- sqrt((sigma + 1) ** 2 - 4 * (1 -r)))
    for r < 1: l < 0, stable conduction
    for r > 1: l_1,2 < 0, l_3 > 0, unstable conduction

    R+- = (+-a, +-a, r - 1), a = sqrt(b * (r - 1))
    eigenvalues: l_1 < 0, l_2,3 in C
    for r < r_c: Re(l_2,3) < 0, stable convection
    for r > r_c: Re(l_2,3) > 0, unstable convection
