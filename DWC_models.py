# -*- coding: utf-8 -*-
"""
some python functions to calculate dropwise condensation

Created on Wed Feb  6 12:53:01 2019

@author: Jakob Sablowski
"""

import math
from CoolProp.CoolProp import PropsSI
import scipy.integrate as integrate
from functools import partial


def KimKim2011(medium="Water", p_steam=120, deltaT_sub=5, Theta=90, CAH=10,
               Theta_a=None, Theta_r=None, k_coat=15, delta_coat=0, h_i=None,
               c=1, N_s=250, print_properties=False, **kwargs):
    """ main function, calculates dropwise condensation heat flux as described in:
    Kim, S., & Kim, K. J. (2011). Dropwise Condensation Modeling Suitable for Superhydrophobic Surfaces. Journal of
    Heat Transfer, 133(8), 081502–081502. https://doi.org/10.1115/1.4003742
    
    Parameters
    ----------
    medium:     str
                defines the condensing fluid to calculate fluid properties using CoolProp, list of viable fluids:
                http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids
    p_steam:    float
                pressure in mbar = 100 Pa
    deltaT_sub: float
                temperature difference to the cooled wall in K
    Theta:      float    
                static contact angle in deg
    CAH:        float
                contact angle hysteresis in deg, only used as a fallback if no values for Theta_a and Theta_r are given
    Theta_a:    float
                advancing contact angle in deg
    Theta_r:    float
                receding contact angle in deg
    k_coat:     float
                thermal conductivity of the coating in W/(mK)
    delta_coat: float
                thickness of the coating in m
    h_i:        float            
                interfacial heat transfer coefficient in MW/m²K, if no value is given, h_i is calculated
    c:          float
                numerical constant, "depends on the shape of the drop and on the steepness of the substrate surface" 
    N_s:        float
                number of Nucleation sites in 10^9 1/m² 
    print_properties: bool
                if set to true, calculated fluid properties are printed
    r_lower:    float, optional
                sets a lower boundary for the heat flux calculation, only droplets with a larger radii are considered
    r_upper:    float, optional
                sets an upper boundary for the heat flux calculation, only droplets with a smaller radii are considered
    Returns
    ----------
    q:          float
                heat flux density in W/m²
    q_n:        float
                heat flux density through small droplets in W/m²
    q_N:        float
                heat flux density through large droplets in W/m²
    r_min:      float      
                minimum droplet radius in m    
    r_e:        float
                effective drop radius in m
    r_max:      float
                effective maximum drop radius in m
    Q_drop:     partial function
                rate of heat flow in W depending on drop radius in m
    n:          partial function 
                drop size distribution for small drops depending on drop radius r in m      
    N:          partial function
                drop size distribuion for large drops depending on drop radius r in m 
    """
    # get kwargs
    r_lower = kwargs.get("r_lower", False)
    r_upper = kwargs.get("r_upper", False)
    # prepare input parameters
    Theta, Theta_a, Theta_r, h_i, N_s, T_sat, sigma, k_c, h_fg, rho, g, R_s, rho_g \
        = init_parameters(Theta_a=Theta_a, Theta_r=Theta_r, Theta=Theta,
                          CAH=CAH, p_steam=p_steam, h_i=h_i, medium=medium, N_s=N_s)
    # calculate interfacial heat transfer coefficient h_i
    h_i_calc = h_i_Schrage(R_s=R_s, T_sat=T_sat, h_fg=h_fg, rho_g=rho_g, sigma_c=1)
    if not h_i:
        h_i = h_i_calc
    # calculate drop radii
    r_min = r_min_KimKim(T_sat=T_sat, sigma=sigma, h_fg=h_fg,
                         rho=rho, deltaT_sub=deltaT_sub)
    r_e = r_e_KimKim(N_s)
    r_max = r_max_KimKim(c=c, Theta_r=Theta_r, Theta_a=Theta_a,
                         Theta=Theta, sigma=sigma, rho=rho, g=g)
    # define functions for rate of heat flow through a single droplet and drop size distribution
    Q_drop = partial(Q_drop_KimKim, deltaT_sub=deltaT_sub, r_min=r_min, delta_coat=delta_coat, k_coat=k_coat, k_c=k_c,
                     Theta=Theta, h_i=h_i)
    Q_drop.__doc__ = "rate of heat flow in W depending on drop radius r in m"
    n = partial(n_KimKim, deltaT_sub=deltaT_sub, r_min=r_min, delta_coat=delta_coat, k_coat=k_coat, k_c=k_c,
                Theta=Theta, h_i=h_i, rho=rho, h_fg=h_fg, r_e=r_e, r_max=r_max)
    n.__doc__ = "drop size distribution for small drops depending on drop radius r in m"
    N = partial(N_LeFevre, r_max=r_max)
    N.__doc__ = "drop size distribution for large drops depending on drop radius r in m"

    # integrate and calculate heat flux density
    def Q_drop_n(r):
        """small drops"""
        Q_drop_n = Q_drop(r) * n(r)
        return Q_drop_n

    def Q_drop_N(r):
        """large drops"""
        Q_drop_N = Q_drop(r) * N(r)
        return Q_drop_N

    # optional boundaries for integration
    if (not r_lower or r_lower < r_min):
        r_lower = r_min
    if (not r_upper or r_upper > r_max):
        r_upper = r_max
    if r_lower < r_e:
        if r_upper > r_e:
            q_n, q_n_interr = integrate.quad(Q_drop_n, r_lower, r_e)
        else:
            q_n, q_n_interr = integrate.quad(Q_drop_n, r_lower, r_upper)
    else:
        q_n = 0
    if r_upper > r_e:
        if r_lower < r_e:
            q_N, q_N_interr = integrate.quad(Q_drop_N, r_e, r_upper)
        else:
            q_N, q_N_interr = integrate.quad(Q_drop_N, r_lower, r_upper)
    else:
        q_N = 0
    q = q_n + q_N

    # calculate additional values
    misc = {}
    misc["Bo"] = bond_number(r_max, sigma, rho, g)

    # optional output of calculated fluid properties
    if print_properties:
        print("\nfluid properties:")
        print("fluid: \t", medium)
        print("T_sat: \t", T_sat-273.15, "°C")
        print("sigma: \t", sigma*1000, "mN/m")
        print("h_fg: \t", h_fg/1000, "kJ/kg")
        print("rho_l:\t", rho, "kg/m³")
        print("rho_g:\t", rho_g, "kg/m³")
        print("R_s:\t", R_s, "J/(kg*K)")
        print("\ninterfacial heat transfer coefficient:")
        print("h_i: \t ", round(h_i, 1), "W/m²K")
        print("h_i_calc:", round(h_i_calc, 1), "W/m²K")
    return q, q_n, q_N, r_min, r_e, r_max, Q_drop, n, N, misc


def init_parameters(Theta, CAH, p_steam, h_i, medium, N_s, **kwargs):
    """ converts all input parameters to SI-units, calculates fluid properties using CoolProp
    """
    Theta_a = kwargs.get("Theta_a")
    Theta_r = kwargs.get("Theta_r")
    # calculates advancing and receding contact angles from CAH
    if not Theta_a:
        Theta_a = Theta + 0.5*CAH
        if Theta_a > 180:
            Theta_a = 180
    if not Theta_r:
        Theta_r = Theta - 0.5*CAH
        if Theta_r < 0:
            Theta_r = 0
    # conversion to SI units
    Theta, Theta_a, Theta_r = [math.radians(x) for x in [Theta, Theta_a, Theta_r]]  # deg --> rad
    p_steam = p_steam * 100     												    # mbar --> Pa
    if h_i:
        h_i = h_i * 1000*1000       											    # MW/m²K --> W/m²K
    N_s = N_s * 1000*1000*1000  												    # 10^9 1/m² -->  1/m²
    # fluid properties
    T_sat = PropsSI("T", "P", p_steam, "Q", 1, medium)                              # boiling temperature in K
    sigma = PropsSI("surface_tension", "P", p_steam, "Q", 0, medium)                # Surface tension in N/m
    k_c = PropsSI("conductivity", "P", p_steam, "Q", 0, medium)                     # thermal conductivity in W/(mK)
    h_fg = PropsSI("Hmass", "P", p_steam, "Q", 1, medium) - \
           PropsSI("Hmass", "P", p_steam, "Q", 0, medium)                           # enthalpy of evaporation in J/kg
    rho_l = PropsSI("Dmass", "P", p_steam, "Q", 0, medium)                          # density of condensate in kg/m³
    rho_g = PropsSI("Dmass", "P", p_steam, "Q", 1, medium)                          # density of gas in kg/m³
    g = 9.81    														            # gravitational acceleration in m/s²
    Mmass = PropsSI("molar_mass", "P", p_steam, "Q", 1, medium)                     # molar mass in kg/mol
    R_s = 8.3144598 / Mmass                                                         # specific gas constant in J/(kg*K)
    return Theta, Theta_a, Theta_r, h_i, N_s, T_sat, sigma, k_c, h_fg, rho_l, g, R_s, rho_g


def Q_drop_KimKim(r, deltaT_sub, r_min, delta_coat, k_coat, k_c, Theta, h_i):
    """ rate of heat flow in W """
    Q_drop = (deltaT_sub * math.pi * r**2 * (1 - r_min/r)) / \
        (delta_coat/(k_coat * (math.sin(Theta))**2) + r*Theta/(4*k_c*math.sin(Theta)) + 1/(2*h_i*(1-math.cos(Theta))))
    return Q_drop


def deltaT_drop_KimKim(r, deltaT_sub, r_min, delta_coat, k_coat, k_c, Theta, h_i):
    """ temperature drop caused by conduction through the drop """
    deltaT_drop = Q_drop_KimKim(r, deltaT_sub, r_min, delta_coat, k_coat, k_c, Theta, h_i) \
        * Theta / (4*math.pi*k_c*math.sin(Theta))
    return deltaT_drop


def R_total(deltaT_sub, Q_drop):
    """ total thermal resistance of a single drop
    
    Parameters
    ----------
    deltaT_sub: float
                temperature difference to the cooled wall in K
    Q_drop:     float
                rate of heat flow through drop in W
    
    Returns
    ----------
    R_total:    float
                total thermal resistance of drop in K/W
    """
    R_total = deltaT_sub / Q_drop
    return R_total


def R_iphase(h_i, radius, Theta):
    """ interafacial thermal resistance of a single drop
    
    Parameters
    ----------
    h_i:        float            
                interfacial heat transfer coefficient in MW/m²K
    radius:     float
                radius of drop in m
    Theta:      float    
                static contact angle in deg                   
    
    Returns
    ----------
    R_iphase:   float
                interafacial thermal resistance of drop in K/W
    """
    R_iphase = 1 / (2*h_i*1000*1000*math.pi*radius**2*(1-math.cos(math.radians(Theta))))
    return R_iphase


def R_cond(k_c, radius, Theta):
    """ thermal resistance due to conduction through a single drop
    
    Parameters
    ----------
    k_c:        float            
                thermal conductivity of condensate in W/mK
    radius:     float
                radius of drop in m
    Theta:      float    
                static contact angle in deg                   
    
    Returns
    ----------
    R_cond:     float
                thermal resistance due to conduction through a single drop in K/W
    """
    R_cond = math.radians(Theta) / (4*math.pi*radius*k_c*math.sin(math.radians(Theta)))
    return R_cond


def R_coat(delta_coat, k_coat, radius, Theta):
    """ thermal resistance due to conduction through coating layer
    
    Parameters
    ----------
    k_coat:     float
                thermal conductivity of the coating in W/(mK)
    delta_coat: float
                thickness of the coating in m
    radius:     float
                radius of drop in m
    Theta:      float    
                static contact angle in deg                   
    
    Returns
    ----------
    R_coat:     float
                thermal resistance due to conduction through coating layer in K/W
    """
    R_coat = delta_coat / (k_coat*math.pi*radius**2*(math.sin(math.radians(Theta)))**2)
    return R_coat


def R_curv(deltaT_sub, r_min, radius, Q_drop):
    """ thermal resistance due drop curvature
    
    Parameters
    ----------
    deltaT_sub: float
                temperature difference to the cooled wall in K
    r_min:      float      
                minimum droplet radius in m    
    radius:     float
                radius of drop in m
    Q_drop:     float
                rate of heat flow through drop in W                
    
    Returns
    ----------
    R_curv:     float
                thermal resistance due drop curvature in K/W
    """
    R_curv = (deltaT_sub*r_min / radius) / Q_drop
    return R_curv


def n_KimKim(r, deltaT_sub, r_min, delta_coat, k_coat, k_c, Theta, h_i, rho, h_fg, r_e, r_max):
    """ drop size distributions small drops """
    A_1 = deltaT_sub / (2*rho*h_fg)
    A_2 = Theta * (1-math.cos(Theta)) / (4*k_c*math.sin(Theta))
    A_3 = 1/(2*h_i) + delta_coat*(1-math.cos(Theta)) / (k_coat*(math.sin(Theta))**2)
    tau = 3*r_e**2 * (A_2*r_e + A_3)**2 / (A_1*(11*A_2*r_e**2 - 14*A_2*r_e*r_min + 8*A_3*r_e - 11*A_3*r_min))
    B_1 = A_2/(tau*A_1) * ((r_e**2-r**2)/2 + r_min*(r_e-r) - r_min**2*math.log((r-r_min)/(r_e-r_min)))
    B_2 = A_3/(tau*A_1) * (r_e-r - r_min*math.log((r-r_min)/(r_e-r_min)))
    n = 1/(3*math.pi*r_e**3*r_max) * (r_e/r_max)**(-2/3) * r*(r_e-r_min)/(r-r_min) * \
        (A_2*r+A_3)/(A_2*r_e+A_3) * math.exp(B_1+B_2)
    return n


def N_LeFevre(r, r_max):
    """ drop size distributions for large drops """
    N = 1/(3*math.pi*r**2*r_max) * (r/r_max)**(-2/3)
    return N


def r_min_KimKim(T_sat, sigma, h_fg, rho, deltaT_sub):
    """ minimum droplet radius """
    r_min = 2*T_sat*sigma / (h_fg * rho * deltaT_sub)
    return r_min


def r_e_KimKim(N_s):
    """ effective drop radius """
    r_e = (4*N_s)**(-0.5)
    return r_e


def r_max_KimKim(c, Theta_r, Theta_a, Theta, sigma, rho, g):
    """ effective maximum drop radius """
    r_max = math.sqrt(6*c*(math.cos(Theta_r)-math.cos(Theta_a))*math.sin(Theta)*sigma / \
            (math.pi*(2-3*math.cos(Theta)+(math.cos(Theta))**3)*rho*g))
    return r_max


def h_i_Schrage(R_s, T_sat, h_fg, rho_g, sigma_c=1):
    """ interfacial heat transfer coefficient as given in
    Glicksman, L. R., & Hunt Jr., A. W. (1972). Numerical simulation of dropwise condensation. International Journal of
    Heat and Mass Transfer, 15(11), 2251–2269. https://doi.org/10.1016/0017-9310(72)90046-4
    """
    v_g = 1/rho_g
    h_i = 2*sigma_c/(2-sigma_c) * math.sqrt(1/(2*math.pi*R_s*T_sat)) * h_fg**2/(T_sat * v_g)
    return h_i


def q_filmwise(medium="Water", p_steam=120, deltaT_sub=5, Theta=90, CAH=10,
               Theta_a=None, Theta_r=None, k_coat=15, delta_coat=0, h_i=None,
               c=1, N_s=250, h_fw=0.10):
    """ calculates heat flux during filmvise condensation according to "Nusseltsche Wasserhauttheorie",
    Equation 4.13 in: H. D. Baehr and K. Stephan, Wärme- und Stoffübertragung,
    8th ed. Berlin, Heidelberg: Springer Berlin Heidelberg, 2013. Equation 4.13
    """
    g = 9.81    														            # gravitational acceleration in m/s²
    # fluid properties
    _, _, _, _, _, theta_s, _, lambda_l, delta_h_v, rho_l, g, _, rho_g \
        = init_parameters(Theta_a=Theta_a, Theta_r=Theta_r, Theta=Theta,
                            CAH=CAH, p_steam=p_steam, h_i=h_i, medium=medium, N_s=N_s)
    eta_l = PropsSI("viscosity", "P", p_steam*100, "Q", 0, medium)
    theta_0 = theta_s - deltaT_sub
    alpha_m = 0.943 * ((rho_l * (rho_l-rho_g) * g * delta_h_v * lambda_l**3) / (eta_l * (theta_s-theta_0) * h_fw))**0.25
    q = alpha_m * deltaT_sub
    return q


def bond_number(r_max, sigma, rho_l, g):
    """ calculates the Bond number for the largest droplet according to
    Cha, H.; Vahabi, H.; Wu, A.; Chavan, S.; Kim, M.-K.; Sett, S.; Bosch, S. A.; Wang, W.; Kota, A. K.; Miljkovic, N.
    Dropwise Condensation on Solid Hydrophilic Surfaces. Science Advances 2020, 6 (2), eaax0746.
    https://doi.org/10.1126/sciadv.aax0746"""
    l_y = math.sqrt(sigma / (rho_l*g))
    bond = r_max**2 / l_y**2
    return bond
