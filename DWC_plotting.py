# -*- coding: utf-8 -*-
"""
some basic plotting functions using PyLab

Created on Thu Feb  7 09:07:39 2019

@author: Jakob Sablowski
"""

import math
import matplotlib.pyplot as plt
import numpy as np
try:
    from . import DWC_models as DWCmod
except:
    import DWC_models as DWCmod

# some input parameters to use as a fallback if no other parameters are given
# most values as given in: Kim, S., & Kim, K. J. (2011). Dropwise Condensation Modeling Suitable for Superhydrophobic
# Surfaces. Journal of Heat Transfer, 133(8), 081502–081502. https://doi.org/10.1115/1.4003742
KimKim2011 = {"medium": "Water",
              "p_steam": 337.8,
              "deltaT_sub": 5,
              "Theta": 90,
              "CAH": 10,                 # estimated
              "k_coat": 0.2,
              "delta_coat": 0.000001,
              "h_i": 15.7,
              "c": 1,
              "N_s": 250}


def choose_model(modelName):
    """ wrapper to choose DWC model """
    if modelName == "KimKim2011":
        DWC_model = DWCmod.KimKim2011
    else:
        raise BaseException("Model \"" + modelName + "\" is unknown. Please choose one of the implemented models.")
    return DWC_model


def print_results(input_params=KimKim2011, model="KimKim2011"):
    """ prints calculated values and results """
    DWC = choose_model(model)
    q, q_n, q_N, r_min, r_e, r_max, Q_drop, n, N = DWC(print_properties=True, **input_params)
    print("\nresults:")
    print("q:\t", q, "W/m²")
    print("q_n:\t", q_n, "W/m²")
    print("q_N:\t", q_N, "W/m²")
    print("q_N/q:\t", 100 * round(q_N/q, 3), "%")
    print("r_min:\t", r_min, "m")
    print("r_e:\t", r_e, "m")
    print("r_max:\t", r_max, "m")
    print("\nmodel used: ", model)
    
    
def plot_qdrop_theta_r(input_params=KimKim2011, model="KimKim2011", radii=[0.000001, 0.000005, 0.000010]):
    """ plot rate of heat flow through a single droplet and heat flux at the base of the droplet vs. the contact angle
    for specific drop radii. The heat flux refers to the contact area between droplet and surface.
    
    Parameters
    ----------
    input_params:   dict
                    input parameters for the DWC model
    model:          str
                    name of the model that should be used
    radii:          list of floats
                    drop radii in m for which a graph should be drawn
    """
    DWC = choose_model(model)
    input_params = input_params.copy()      # avoid changing global input_params    
    Theta = np.linspace(5, 175, 20)
    # rate of heat flow
    axs = []
    fig = plt.figure()
    plt.subplot(2, 1, 1)
    for r in radii:
        Q_d = []
        for x in Theta:
            input_params["Theta"] = x
            Q_drop = DWC(**input_params)[6]
            Q_d.append(Q_drop(r)*1000)  # Q in mW
        axs.append(plt.plot(Theta, Q_d, label=r"$r = $" + str(r*1000000) + r"$\ \mathrm{\mu m}$"))
    plt.ylabel(r"$\.Q_{\mathrm{drop}} \ \mathrm{in \ mW}$")
    plt.legend()
    # heat flux
    axs = []
    plt.subplot(2, 1, 2)
    for r in radii:
        q_d = []
        for x in Theta:
            input_params["Theta"] = x
            Q_drop = DWC(**input_params)[6]
            # A = math.pi * r**2                                  # projected are of droplet in m²
            A = math.pi * r**2 * (math.sin(math.radians(x)))**2  # contact area between droplet and surface in m²
            q_d.append(Q_drop(r)/A/1000)                         # q in kW/m²  
        axs.append(plt.plot(Theta, q_d, label="r = " + str(r*1000000) + r"$\ \mu m$"))
    plt.ylabel(r"$\.q_{\mathrm{drop}} \ \mathrm{in \ kW/m^2}$")
    plt.xlabel(r"$\theta$ in deg")
    plt.show()
    return fig


def plot_Rdrop(input_params=KimKim2011, model="KimKim2011"):
    """ plot heat transfer resistances over droplet radius
    """
    DWC = choose_model(model)
    l_radius = np.linspace(0.00000001, 0.001, 100)
    input_params = input_params.copy()                  # avoid changing global input_params    
    k_c = DWCmod.init_parameters(**input_params)[7]     # calculate thermal conductivity of condensate
    fig = plt.figure()
    R_total = []
    R_iphase = [] 
    R_cond = []
    R_coat = []
    R_curv = []
    r_e = DWC(**input_params)[4]
    r_min = DWC(**input_params)[3]
    for radius in l_radius:
        Q_drop = DWC(**input_params)[6]
        R_total.append(DWCmod.R_total(deltaT_sub=input_params["deltaT_sub"], Q_drop=Q_drop(radius)))
        R_iphase.append(DWCmod.R_iphase(h_i=input_params["h_i"], radius=radius, Theta=input_params["Theta"]))
        R_cond.append(DWCmod.R_cond(k_c=k_c, radius=radius, Theta=input_params["Theta"]))
        R_coat.append(DWCmod.R_coat(delta_coat=input_params["delta_coat"], k_coat=input_params["k_coat"], radius=radius,
                                    Theta=input_params["Theta"]))
        R_curv.append(DWCmod.R_curv(deltaT_sub=input_params["deltaT_sub"], r_min=r_min, radius=radius,
                                    Q_drop=Q_drop(radius)))
    plt.plot(l_radius, R_total, "-.", label=r"$R_{\mathrm{total}}$")
    plt.plot(l_radius, R_iphase, label=r"$R_{\mathrm{interphase}}$")
    plt.plot(l_radius, R_curv, label=r"$R_{\mathrm{curvature}}$")
    plt.plot(l_radius, R_cond, label=r"$R_{\mathrm{conduction}}$")
    plt.plot(l_radius, R_coat, label=r"$R_{\mathrm{coating}}$")
    plt.axvline(x=r_e, linestyle="--", label=r"$r_{\mathrm{e}}$")
    plt.yscale("log")
    plt.xscale("log")
    plt.legend()
    plt.ylabel(r"$R \ \mathrm{in \ K/W}$")
    plt.xlabel(r"$r \ \mathrm{in \ m}$")
    plt.ylim(top=10**11)
    plt.show()
    return fig

# TODO rewrite plotting functions similar to plot_q_deltaTsub()
def plot_q_theta(input_params=KimKim2011, model="KimKim2011", **kwargs):
    """ wrapper function to plot the heat flux vs. the contact angle 

    Parameters
    ----------
    input_params:   dict
                    input parameters for the DWC model
    model:          str
                    name of the model that should be used
    CAH:            list of floats, optional
                    contact angles hystereses in deg for which a graph should be drawn
    """
    CAH = kwargs.get("CAH", [input_params["CAH"]])
    fig = plot_q_theta_CAH(input_params, model, CAH=CAH)    
    return fig


def plot_q_theta_CAH(input_params=KimKim2011, model="KimKim2011", CAH=[3, 5, 10, 20, 40]):
    """ plot the heat flux vs. the contact angle for specific contact angle hystereses.

    Parameters
    ----------
    input_params:   dict
                    input parameters for the DWC model
    model:          str
                    name of the model that should be used
    CAH:            list of floats
                    contact angles hystereses in deg for which a graph should be drawn
    """
    DWC = choose_model(model)
    input_params = input_params.copy()      # avoid changing global input_params
    Theta = np.linspace(20, 160, 40)
    axs = []
    fig = plt.figure()
    for y in CAH:
        input_params["CAH"] = y
        q = []
        for x in Theta:
            input_params["Theta"] = x
            q.append(DWC(**input_params)[0]/1000)
        axs.append(plt.plot(Theta, q, label="CAH = " + str(input_params["CAH"]) + "°"))
    plt.ylabel(r"$\.q \ \mathrm{in \ kW/m^2}$")
    plt.xlabel(r"$\theta \ \mathrm{in \ deg}$")
    plt.legend()
    plt.show()    
    return fig


def plot_Nr_r(input_params=KimKim2011, model="KimKim2011", **kwargs):
    """ wrapper function to plot the drop size distribution
    
    Parameters
    ----------
    input_params:   dict
                    input parameters for the DWC model
    model:          str
                    name of the model that should be used
    theta:          list of floats, optional
                    contact angles in deg for which a graph should be drawn
    CAH:            list of floats, optional
                    contact angles hystereses in deg for which a graph should be drawn    
    """
    theta = kwargs.get("theta", [input_params["Theta"]])
    CAH = kwargs.get("CAH")
    N_s = kwargs.get("N_s")
    if CAH:
        fig = plot_Nr_r_CAH(input_params, model, CAH=CAH)    
    elif N_s:
        fig = plot_Nr_r_Ns(input_params, model, N_s=N_s) 
    else:
        fig = plot_Nr_r_theta(input_params, model, theta=theta)    
    return fig


def plot_Nr_r_theta(input_params=KimKim2011, model="KimKim2011", theta=[90, 120, 150]):
    """ plot the drop size distribution for specific contact angles.

    Parameters
    ----------
    input_params:   dict
                    input parameters for the DWC model
    model:          str
                    name of the model that should be used
    thetas:         list of floats
                    contact angles in deg for which a graph should be drawn                
    """ 
    DWC = choose_model(model)
    input_params = input_params.copy()      # avoid changing global input_params
    axs = []
    fig = plt.figure()
    for y in theta:
        input_params["Theta"] = y
        q, q_n, q_N, r_min, r_e, r_max, Q_drop, n, N = DWC(**input_params)
        r_n = np.linspace(r_min, r_e, 50)
        r_n = r_n[1:]
        n = [n(x) for x in r_n]
        r_N = np.linspace(r_e, r_max, 50)
        r_N = r_N[1:]
        N = [N(x) for x in r_N]
        r_ges = np.append(r_n, r_N)
        n_ges = np.append(n, N)
        axs.append(plt.loglog(r_ges, n_ges, label=r"$\theta =$" + str(input_params["Theta"]) + "°"))
    plt.ylabel(r"$n(r) \ \mathrm{and} \ N(r) \ \mathrm{in \ m^{-3}}$")
    plt.xlabel(r"$r \ \mathrm{in \ m}$")
    plt.xlim(right=r_max)
    plt.ylim(bottom=10**8, top=10**18)
    plt.axvline(x=r_e, linestyle="--", label=r"$r_{\mathrm{e}}$")
    plt.show()
    return fig             


def plot_Nr_r_CAH(input_params=KimKim2011, model="KimKim2011", CAH=[3, 10, 40]):
    """ plot the drop size distribution for specific contact angle hystereses.

    Parameters
    ----------
    input_params:   dict
                    input parameters for the DWC model
    model:          str
                    name of the model that should be used
    CAH:            list of floats
                    contact angles hystereses in deg for which a graph should be drawn
    """ 
    DWC = choose_model(model)
    input_params = input_params.copy()      # avoid changing global input_params
    axs = []
    fig = plt.figure()
    for y in CAH:
        input_params["CAH"] = y
        q, q_n, q_N, r_min, r_e, r_max, Q_drop, n, N = DWC(**input_params)
        r_n = np.linspace(r_min, r_e, 50)
        r_n = r_n[1:]
        n = [n(x) for x in r_n]
        r_N = np.linspace(r_e, r_max, 50)
        r_N = r_N[1:]
        N = [N(x) for x in r_N]
        r_ges = np.append(r_n, r_N)
        n_ges = np.append(n, N)
        axs.append(plt.loglog(r_ges, n_ges, label="CAH = " + str(input_params["CAH"]) + "°"))
    plt.ylabel(r"$n(r) \ \mathrm{and} \ N(r) \ \mathrm{in \ m^{-3}}$")
    plt.xlabel(r"$r \ \mathrm{in \ m}$")
    plt.xlim(right=r_max)
    plt.ylim(bottom=10**8, top=10**18)
    plt.axvline(x=r_e,  label=r"$r_{\mathrm{e}}$")
    plt.legend()
    plt.show()
    return fig        


def plot_Nr_r_Ns(input_params=KimKim2011, model="KimKim2011", N_s=[2.5, 25, 250, 2500]):
    """ plot the drop size distribution for specific contact angle hystereses.

    Parameters
    ----------
    input_params:   dict
                    input parameters for the DWC model
    model:          str
                    name of the model that should be used
    N_s:            list of floats
                    nucleation site densities for which a graph should be drawn
                    in 10^9 1/m² 
    """ 
    DWC = choose_model(model)
    input_params = input_params.copy()      # avoid changing global input_params
    axs = []
    fig = plt.figure()
    for y in N_s:
        input_params["N_s"] = y
        q, q_n, q_N, r_min, r_e, r_max, Q_drop, n, N = DWC(**input_params)
        r_n = np.linspace(r_min, r_e, 50)
        r_n = r_n[1:]
        n = [n(x) for x in r_n]
        r_N = np.linspace(r_e, r_max, 50)
        r_N = r_N[1:]
        N = [N(x) for x in r_N]
        r_ges = np.append(r_n, r_N)
        n_ges = np.append(n, N)
        axs.append(plt.loglog(r_ges, n_ges,
                              label=r"$N_{\mathrm{s}}=$" + str('%.1e' % (y * 10**9)) + r"$\ \rm{m}^{-2}$"))
        # axs.append(plt.axvline(x=r_e,  label=r"$r_e$ (" + str(y) + ")"))
    plt.ylabel(r"$n(r) \ \mathrm{and} \ N(r) \ \mathrm{in \ m^{-3}}$")
    plt.xlabel(r"$r \ \mathrm{in \ m}$")
    plt.xlim(right=r_max)
    plt.ylim(bottom=10**8, top=10**18)
    plt.legend()
    plt.show()
    return fig        


def plot_q_deltaTsub(input_params=KimKim2011, model="KimKim2011", **kwargs):
    """ wrapper function to plot the heat flux vs. the surface subcooling temperature 
    
    Parameters
    ----------
    input_params:   dict
                    input parameters for the DWC model
    model:          str
                    name of the model that should be used
    theta:          list of floats, optional
                    contact angles in deg for which a graph should be drawn    
    c:              list of floats, optional
                    constants c for which a graph should be drawn    
    h_i:            list of floats, optional
                    interfacial heat transfer coefficients in MW/m²K for which a graph should be drawn   
    CAH:            list of floats, optional
                    contact angle hystereses in deg for which a graph should be drawn 
    N_s:            list of floats, optional
                    number of nucleation sites per unit area in 10^9 1/m² for which a graph should be drawn
    medium:         list of strings
                    fluids, must be included in CoolProp library
    """
    theta = kwargs.get("theta", [input_params["Theta"]])
    c = kwargs.get("c")
    h_i = kwargs.get("h_i")
    CAH = kwargs.get("CAH")
    N_s = kwargs.get("N_s")
    medium = kwargs.get("medium")
    filmwise = kwargs.get("filmwise")
    if c:
        fig = plot_q_deltaTsub_var(input_params, model, var=c, varname="c", filmwise=filmwise)
    elif h_i:
        fig = plot_q_deltaTsub_var(input_params, model, var=h_i, varname="h_i", filmwise=filmwise)
    elif CAH:
        fig = plot_q_deltaTsub_var(input_params, model, var=CAH, varname="CAH", filmwise=filmwise)
    elif N_s:
        fig = plot_q_deltaTsub_var(input_params, model, var=N_s, varname="N_s", filmwise=filmwise)
    elif medium:
        fig = plot_q_deltaTsub_var(input_params, model, var=medium, varname="medium", filmwise=filmwise)
    else:
        fig = plot_q_deltaTsub_var(input_params, model, var=theta, varname="Theta", filmwise=filmwise)
    return fig


def labelnames_q_deltaTsub(varname):
    """ defines lambda functions for naming the labels in plot_q_deltaTsub according to kwargs"""
    if varname == "medium":
        return lambda y : str(y)
    if varname == "c":
        return lambda y : "$c = $" + str(y)
    if varname == "h_i":
        return lambda y : "$h_{\mathrm{i}} = $" + str(y) + " MW/m²K"
    if varname == "CAH":
        return lambda y : "CAH = " + str(y) + "°"
    if varname == "N_s":
        return lambda y : "$N_{\mathrm{s}} =$" + str('%.1e' % (y * 10**9)) + "$\ \mathrm{m^{-2}}$"
    if varname == "Theta":
        return lambda y : str(y) + "°"
    else:
        return lambda y: str(y)


def plot_q_deltaTsub_var(input_params=KimKim2011, model="KimKim2011", var=[90], varname="Theta", filmwise=False):
    """ plot the heat flux vs. the surface subcooling temperature, optional: vary one additional parameter

    Parameters
    ----------
    input_params:   dict
                    input parameters for the DWC model
    var:            list
                    list with values for the addtional parameter
    varname:        string
                    name of the additional parameter according to input_params
    filmwise:       bool
                    calculate filmwise condensation for comparison if true
    """
    label = labelnames_q_deltaTsub(varname)
    DWC = choose_model(model)
    input_params = input_params.copy()      # avoid changing global input_params
    deltaT_sub = np.linspace(0.1, 10, 20)
    axs = []
    fig = plt.figure()
    for y in var:
        input_params[varname] = y
        q = []
        q_fw = []
        for x in deltaT_sub:
            input_params["deltaT_sub"] = x
            q.append(DWC(**input_params)[0]/1000)
            q_fw.append(DWCmod.q_filmwise(**input_params, H=0.010)/1000)
        axs.append(plt.plot(deltaT_sub, q, label=label(y)))
        if filmwise:
            color = axs[-1][0].get_color()
            axs.append(plt.plot(deltaT_sub, q_fw, color=color, linestyle="--", label=label(y) + " film"))
    plt.ylabel(r"$\.q \ \mathrm{in \ kW/m^2}$")
    plt.xlabel(r"$\Delta T \ \mathrm{in \ K}$")
    plt.legend()
    return fig


def plot_all(input_params=KimKim2011, model="KimKim2011"):
    """ plot all the plots! """
    plot_qdrop_theta_r(input_params, model)
    plot_q_theta(input_params, model)
    plot_q_theta_CAH(input_params, model)
    plot_Nr_r(input_params, model)
    plot_Nr_r_theta(input_params, model)
    plot_Nr_r_CAH(input_params, model)
    plot_q_deltaTsub(input_params, model)
    print_results(input_params, model)