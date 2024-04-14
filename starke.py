from grumpy import Load as gpLoad
from colorama import Fore
import numpy as np
import matplotlib.pyplot as plt


def load(file_name = '', path = '', shup = True, **kwargs):
    data = gpLoad(file_name = file_name, path = path, shup = shup, **kwargs)
    if data.get("Type", "").startswith("spin_deflector"):
        print(Fore.BLUE + "This is spin MDC data.")
        print(f"Energy range:    {data['x'][0]:.3f} to {data['x'][-1]:.3f} eV, {len(data['x'])} points.")
        print(f"Deflector range: {data['y'][0]:.3f} to {data['y'][-1]:.3f} eV, {len(data['y'])} points.")
        print(f"Pass energy: {data.get('Experiment', {}).get('Ep', np.NaN)} eV" + Fore.RESET)
        return data
    else:
        print(Fore.RED + "This is not spin MDC data." + Fore.RESET)
        return {}


def plotSpinMDC(data = {}, asym_shift = 0):
    if not data.get("Type", "").startswith("spin_deflector"):
        print(Fore.RED + "This is not spin MDC data." + Fore.RESET)
        return
    #
    fig, ax = plt.subplots(ncols = 5, figsize = (12,3))
    energy = data["x"]
    deflection = data["y"]
    extent = [energy[0], energy[-1], deflection[-1], deflection[0]]
    #
    ax[0].imshow(data["int"], extent = extent)
    ax[1].imshow(data["int_on"], extent = extent)
    #
    asymmetry = (data["int"] - data["int_on"]) / (data["int"] + data["int_on"])
    vmin, vmax = asymmetry.min(), asymmetry.max()
    ax[2].imshow(asymmetry + asym_shift, extent = extent, cmap = "bwr", vmin = vmin, vmax = vmax)
    #
    c1 = (data["int"] + data["int_on"])/2 * (1 + asymmetry)
    c2 = (data["int"] + data["int_on"])/2 * (1 - asymmetry)
    ax[3].imshow(c1, extent = extent)
    ax[4].imshow(c2, extent = extent)
    #
    for i, ttl in enumerate(["Polarity Off", "Polarity On", "Asymmetry", "Component Off", "Component On"]):
        ax[i].set_title(ttl, fontsize = 10)
        ax[i].invert_yaxis()
        ax[i].set_xlabel("Kinetic energy, eV", fontsize = 9)
    ax[0].set_ylabel("Deflection", fontsize = 9)


def extractEDCfromMDC(data = {}, number = 0, plot = True):
    if not data.get("Type", "").startswith("spin_deflector"):
        print(Fore.RED + "This is not spin MDC data." + Fore.RESET)
        return
    #
    if number >= len(data["y"]):
        print(Fore.RED + f"There are only {len(data['y'])} edcs in this data set, so argument number must be 0 - {len(data['y'])-1}." + Fore.RESET)
        return
    #
    edc = {"defl": data["y"][number], "energy": data["x"], "intensity_off": data["int"][number], "intensity_on": data["int_on"][number]}
    asym = (data["int"][number] - data["int_on"][number]) / (data["int"][number] + data["int_on"][number])
    edc.update({"asymmetry": asym})
    #
    if plot:
        fig, ax = plt.subplots(ncols = 3, figsize = (9,3))
        ax[0].plot(edc["energy"], edc["intensity_off"], label = "off")
        ax[0].plot(edc["energy"], edc["intensity_on"], label = "on")
        ax[1].plot(edc["energy"], edc["intensity_on"], label = "asymmetry")
        c1 = (edc["intensity_off"] + edc["intensity_on"])/2 * (1 + asym)
        c2 = (edc["intensity_off"] + edc["intensity_on"])/2 * (1 - asym)
        ax[2].plot(edc["energy"], c1, label = "comp. off")
        ax[2].plot(edc["energy"], c2, label = "comp. on")
        for i, ttl in enumerate(["edc", "asymmetry", "components"]):
            ax[i].set_title(ttl, fontsize = 10)
            if i in [0,2]: ax[i].legend(fontsize = 8)
        fig.tight_layout()
    return edc


def plotMDCasEDC(data = {}, style = "intensities", figsize = (5,10), stack = 1000):
    if not data.get("Type", "").startswith("spin_deflector"):
        print(Fore.RED + "This is not spin MDC data." + Fore.RESET)
        return
    #
    style_opts = ["intensities", "components", "asymmetries"]
    if not style in style_opts:
        print(Fore.RED + f"The argument style must be {style_opts}" + Fore.RESET)
        return
    fig, ax = plt.subplots(figsize = figsize)
    deflection = data.get("y")
    energy = data.get("x")
    for i, defl in enumerate(deflection):
        edc = extractEDCfromMDC(data = data, number = i, plot = False)
        if style == "intensities":
            ax.plot(energy, edc["intensity_off"] + i*stack, color = "tab:blue")
            ax.plot(energy, edc["intensity_on"]  + i*stack, color = "tab:red")
        elif style == "components":
            c1 = (edc["intensity_off"] + edc["intensity_on"])/2 * (1 + edc["asymmetry"])
            c2 = (edc["intensity_off"] + edc["intensity_on"])/2 * (1 - edc["asymmetry"])
            ax.plot(energy, c1 + i*stack, color = "tab:blue")
            ax.plot(energy, c2 + i*stack, color = "tab:red")
        elif style == "asymmetries":
            ax.plot(energy, edc["asymmetry"] + i*stack, color = "blue")
    ax.set_title("Kinetic energy, eV")
    fig.tight_layout()


def plotMDCpolarization(datac2p = {}, datac2m = {}, datac1 = {}, stack = 1, figsize = (5,10), sherman = 0.29, ylim = (None, None), xlim = (None,None)):
    #
    if datac2p.get("Type", "").startswith("spin_deflector") and datac2m.get("Type", "").startswith("spin_deflector") and datac1.get("Type", "").startswith("spin_deflector"):
        thecase = "xyz"
        energy = datac2p["x"]
        deflector = datac2p["y"]
    elif datac2p.get("Type", "").startswith("spin_deflector") and datac2m.get("Type", "").startswith("spin_deflector"):
        thecase = "xy"
        energy = datac2p["x"]
        deflector = datac2p["y"]
    elif datac1.get("Type", "").startswith("spin_deflector"):
        thecase = "z"
        energy = datac1["x"]
        deflector = datac1["y"]
    else:
        return
    #
    fig, ax = plt.subplots(figsize = figsize)
    for i, defl in enumerate(deflector):
        d1, d2, d3 = [], [], []
        if thecase in ["xyz", "xy"]:
            d1 = extractEDCfromMDC(data = datac2p, number = i, plot = False)
            d2 = extractEDCfromMDC(data = datac2m, number = i, plot = False)
        if thecase in ["xyz", "z"]:
            d3 = extractEDCfromMDC(data = datac1, number = i, plot = False)
        #
        if thecase == "xyz":
            tot_int = d1["intensity_off"] + d1["intensity_on"] + d2["intensity_off"] + d2["intensity_on"] + d3["intensity_off"] + d3["intensity_on"]
            tot_curves = 6
        elif thecase == "xy":
            tot_int = d1["intensity_off"] + d1["intensity_on"] + d2["intensity_off"] + d2["intensity_on"]
            tot_curves = 4
        elif thecase == "z":
            tot_int = d3["intensity_off"] + d3["intensity_on"]
            tot_curves = 2
        #
        if thecase == "xyz":
            px = 1/np.sqrt(2)/sherman * (d1["asymmetry"] - d2["asymmetry"])
            py = -1/np.sqrt(2)/sherman * (d1["asymmetry"] + d2["asymmetry"])
            pz = 1/sherman * d3["asymmetry"]
            ax.plot(energy, px + i*stack, color = "tab:blue", label = "$P_x$")
            ax.plot(energy, py + i*stack, color = "tab:green", label = "$P_y$")
            ax.plot(energy, pz + i*stack, color = "tab:red", label = "$P_z$")
        elif thecase == "xy":
            px = 1/np.sqrt(2)/sherman * (d1["asymmetry"] - d2["asymmetry"])
            py = -1/np.sqrt(2)/sherman * (d1["asymmetry"] + d2["asymmetry"])
            ax.plot(energy, px + i*stack, color = "tab:blue", label = "$P_x$")
            ax.plot(energy, py + i*stack, color = "tab:green", label = "$P_y$")
        elif thecase == "z":
            ax.plot(energy, pz + i*stack, color = "tab:red", label = "$P_z$")

    ax.set_title("Kinetic energy, eV")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    fig.tight_layout()







    