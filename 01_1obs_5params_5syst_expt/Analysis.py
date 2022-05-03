# based on JETSCAPE/STAT/Example.ipynb

#--------------------------------------------------------------------------------------------------------------------------------
#
# Analysis.py  syst1[,syst2[,...]]  obs1[,obs2[,...]]  npc  nwalkers  nburnsteps  nprodsteps
#
#--------------------------------------------------------------------------------------------------------------------------------
import os
import re
import subprocess
import sys

from sklearn.gaussian_process import GaussianProcessRegressor as GPR
from sklearn.gaussian_process import kernels
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import matplotlib.cm as cm
import matplotlib.pyplot as plt

from scipy.linalg import lapack
from scipy import stats
import emcee  # note: emcee 2.2.1 is required; run `pip install --user emcee==2.2.1`
import numpy as np

import importlib

import os
import pickle
from pathlib import Path

import src.reader as Reader
from src.design import Design


import datetime


if (len(sys.argv[1:]) != 6):
    print(f"Usage: {sys.argv[0]} system[,system[,...]] observable[,observable[,...]] npc walkers nburnsteps nprodsteps", file=sys.stderr)
    exit(1)

arg_systems, arg_observables, arg_npc, arg_nwalkers, arg_nburnsteps, arg_nprodsteps = sys.argv[1:]


#--------------------------------------------------------------------------------------------------------------------------------
design = Design()

dataSet = set()
predSet = set()
setLookup = { "Data": dataSet, "Prediction": predSet }

for file in os.listdir("processed"):
    match = re.fullmatch(f"^({'|'.join(map(str, setLookup))})(-main)?-([A-Za-z0-9]+)-([A-Za-z0-9_]+)[.]dat$", file)
    if (match is None): continue

    groups = match.groups()
    if (len(groups) != 4): continue

    which, _, systag, obsname = groups

    setLookup[which].add( (systag, obsname) )  # (system, observable)

sysobslist, sysobscheck = [sorted(f(*setLookup.values())) for f in (set.intersection, set.union)]

if (sysobslist != sysobscheck):
    raise Exception(f"inconsistency between data and predictions: {set(sysobscheck) - set(sysobslist)}")

systems = sorted(arg_systems.split(','))

observables_dict = { systag: set() for systag in systems }

for systag, obsname in sysobslist:
    if systag in observables_dict:
        observables_dict[systag].add(obsname)

for systag, obsnames in observables_dict.items():
    print(f"{systag:10}: {obsnames}", end="\n\n", flush=True)

# Read data files
RawData = {
    systag: {
        obs: Reader.ReadData(f"processed/Data-{systag}-{obs}.dat")
        for obs in observables_dict[systag]
    }
    for systag in systems
}

# format of dictionary returned by Reader.ReadData():
#
# . FileName: "processed/Data-PbPb5020-obsname.dat"
# . DOI: []
# . Source: []
# . System: "PbPb5020"
# . Centrality: ["0", "70"]
# . "XY": ["Centrality", "obsname"]
# . "Label": ["x", "y", "stat,low", "stat,high", "sys,low", "sys,high"]
# . "Data":
#   . "x": np.array([centrality bin centers])
#   . "y": np.array([observable values per centrality bin])
#   . "yerr"
#     . "stat": np.array([[low, high], ... per centrality bin])
#     . "sys": np.array([[low, high], ... per centrality bin])
# . "SysLabel": ["sys,low", "sys,high"]

# Read design points
RawDesign = Reader.ReadDesign("processed/Design.dat")
RawDesign["Design"] = RawDesign["Design"][[_ not in design.bad_points for _ in range(len(RawDesign["Design"]))]]

# Read model prediction
RawPredictions = {
    systag: {
        obs: Reader.ReadPrediction(f"processed/Prediction-main-{systag}-{obs}.dat")
        for obs in observables_dict[systag]
    }
    for systag in systems
}


#--------------------------------------------------------------------------------------------------------------------------------
# Initialize empty dictionary
AllData = {}

# Basic information
AllData["systems"] = list(systems)
AllData["keys"] = design.keys
AllData["labels"] = design.labels
AllData["ranges"] = design.range
AllData["observables"] = {
    systag: [ (obs, ["C0"]) for obs in obslist ]
    for systag, obslist in observables_dict.items()
}

# Data points
Data = {
    systag: {
        obs: {
            "C0": RawData[systag][obs]["Data"]
        } for obs in observables_dict[systag]
    }
    for systag in systems
}

# Model predictions
Prediction = {
    systag: {
        obs: {
            "C0": { "Y": RawPredictions[systag][obs]["Prediction"], "x": RawData[systag][obs]["Data"]["x"] },
        } for obs in observables_dict[systag]
    }
    for systag in systems
}

# Covariance matrices - the indices are [system][measurement1][measurement2], each one is a block of matrix
Covariance = Reader.InitializeCovariance(Data)
for systag in systems:
    for obs in observables_dict[systag]:
        Covariance[systag][(obs, "C0")][(obs, "C0")] = Reader.EstimateCovariance(RawData[systag][obs], RawData[systag][obs])

# Assign data to the dictionary
AllData["design"] = RawDesign["Design"]
AllData["model"] = Prediction
AllData["data"] = Data
AllData["cov"] = Covariance

# Save to the desired pickle file
with open("input/default.p", "wb") as handle:
    pickle.dump(AllData, handle, protocol = pickle.HIGHEST_PROTOCOL)


#--------------------------------------------------------------------------------------------------------------------------------
# Clean past MCMC samples
if os.path.exists("cache/mcmc_chain.hdf"):
    os.remove("cache/mcmc_chain.hdf")

# Clean past emulator
for system in AllData["systems"]:
    if os.path.exists("cache/emulator/" + system + ".pkl"):
        os.remove("cache/emulator/" + system + ".pkl")


#--------------------------------------------------------------------------------------------------------------------------------
#! python3 -m src.emulator --retrain --npc 10
subprocess.run(("python3", "-m", "src.emulator", "--retrain", "--npc", arg_npc))


#--------------------------------------------------------------------------------------------------------------------------------
from src import lazydict, emulator
Emulators = {
    systag: emulator.Emulator.from_cache(systag)
    for systag in systems
}


#--------------------------------------------------------------------------------------------------------------------------------
if os.path.exists("cache/mcmc_chain.hdf"):
    os.remove("cache/mcmc_chain.hdf")
#! python3 -m src.mcmc --nwalkers 500 --nburnsteps 500 1500
subprocess.run(("python3", "-m", "src.mcmc", "--nwalkers", arg_nwalkers, "--nburnsteps", arg_nburnsteps, arg_nprodsteps))


#--------------------------------------------------------------------------------------------------------------------------------
import src
src.Initialize()
from src import mcmc
chain = mcmc.Chain()
MCMCSamples = chain.load()

TransformedSamples = np.copy(MCMCSamples)


#--------------------------------------------------------------------------------------------------------------------------------
#! python3 -m src.plots posterior gp diag_emu
subprocess.run(("python3", "-m", "src.plots", "posterior", "gp", "diag_emu"))


#--------------------------------------------------------------------------------------------------------------------------------
with chain.dataset() as d:
    W = d.shape[0]     # number of walkers
    S = d.shape[1]     # number of steps
    N = d.shape[2]     # number of parameters
    T = int(S / 200)   # "thinning"
    A = 20 / W
    figure, axes = plt.subplots(figsize = (15, 2 * N), ncols = 1, nrows = N)
    for i, ax in enumerate(axes):
        for j in range(0, W):
            ax.plot(range(0, S, T), d[j, ::T, i], alpha = A)
    #plt.tight_layout(True)
    plt.tight_layout()
    plt.savefig("plots/MCMCSamples.pdf", dpi = 192)


#--------------------------------------------------------------------------------------------------------------------------------
NDimension = len(AllData["labels"])
Ranges = np.array(AllData["ranges"]).T
figure, axes = plt.subplots(figsize = (3 * NDimension, 3 * NDimension), ncols = NDimension, nrows = NDimension, squeeze = False)
Names = AllData["labels"]
for i, row in enumerate(axes):
    for j, ax in enumerate(row):
        if i==j:
            ax.hist(MCMCSamples[:,i], bins=50,
                    range=Ranges[:,i], histtype="step", color="green")
            ax.set_xlabel(Names[i])
            ax.set_xlim(*Ranges[:,j])
        if i>j:
            ax.hist2d(MCMCSamples[:, j], MCMCSamples[:, i],
                      bins=50, range=[Ranges[:,j], Ranges[:,i]],
                      cmap="Greens")
            ax.set_xlabel(Names[j])
            ax.set_ylabel(Names[i])
            ax.set_xlim(*Ranges[:,j])
            ax.set_ylim(*Ranges[:,i])
        if i<j:
            ax.axis("off")
#plt.tight_layout(True)
plt.tight_layout()
plt.savefig("plots/Correlation.pdf", dpi = 192)
# figure


#--------------------------------------------------------------------------------------------------------------------------------
Names = design.labels
NDimension = len(Names)
Ranges = np.array(AllData["ranges"]).T
figure, axes = plt.subplots(figsize = (15, 15), ncols = NDimension, nrows = NDimension, squeeze = False)
for i, row in enumerate(axes):
    for j, ax in enumerate(row):
        if i==j:
            ax.hist(TransformedSamples[:,i], bins=50,
                    range=Ranges[:,i], histtype="step")
            ax.set_xlabel(Names[i])
            ax.set_xlim(*Ranges[:,j])
        if i>j:
            ax.hist2d(TransformedSamples[:, j], TransformedSamples[:, i],
                      bins=50, range=[Ranges[:,j], Ranges[:,i]],
                      cmap="Blues")
            ax.set_xlabel(Names[j])
            ax.set_ylabel(Names[i])
            ax.set_xlim(*Ranges[:,j])
            ax.set_ylim(*Ranges[:,i])
        if i<j:
            ax.axis("off")
#plt.tight_layout(True)
plt.tight_layout()
plt.savefig("plots/TransformedCorrelation.pdf", dpi = 192)
# figure


#--------------------------------------------------------------------------------------------------------------------------------
Examples = MCMCSamples[ np.random.choice(range(len(MCMCSamples)), 2500), :]

TempPrediction = {
    systag: Emulators[systag].predict(Examples)
    for systag in systems
}

SystemCount = len(AllData["systems"])

figure, axes = plt.subplots(figsize = (15, 5 * SystemCount), ncols = 2, nrows = SystemCount, squeeze = False)

for s1 in range(0, SystemCount):
    for s2 in range(0, 1):
        axes[s1][s2].set_xlabel(r"Centrality")
        axes[s1][s2].set_ylabel(r"$\left<E_\mathrm{T}\right>$")

        S1 = AllData["systems"][s1]
        O  = AllData["observables"][S1][0][0]
        S2 = AllData["observables"][S1][0][1][s2]

        DX = AllData["data"][S1][O][S2]["x"]
        DY = AllData["data"][S1][O][S2]["y"]
        DE = np.sqrt(AllData["data"][S1][O][S2]["yerr"]["stat"][:,0]**2 + AllData["data"][S1][O][S2]["yerr"]["sys"][:,0]**2)

        for i, y in enumerate(TempPrediction[S1][O][S2]):
            axes[s1][s2].plot(DX, y, "b-", alpha=0.005, label="Posterior" if i==0 else "")
        axes[s1][s2].errorbar(DX, DY, yerr = DE, fmt="ro", label="Measurements")

#plt.tight_layout(True)
plt.tight_layout()
figure.savefig("plots/ObservablePosterior.pdf", dpi = 192)
# figure


#--------------------------------------------------------------------------------------------------------------------------------
Examples = AllData["design"]

TempPrediction = {
    systag: Emulators[systag].predict(Examples)
    for systag in systems
}

SystemCount = len(AllData["systems"])

figure, axes = plt.subplots(figsize = (15, 5 * SystemCount), ncols = 2, nrows = SystemCount, squeeze = False)

for s1 in range(0, SystemCount):
    for s2 in range(0, 1):
        axes[s1][s2].set_xlabel(r"Centrality")
        axes[s1][s2].set_ylabel(r"$\left<E_\mathrm{T}\right>$")

        S1 = AllData["systems"][s1]
        O  = AllData["observables"][S1][0][0]
        S2 = AllData["observables"][S1][0][1][s2]

        DX = AllData["data"][S1][O][S2]["x"]
        DY = AllData["data"][S1][O][S2]["y"]
        DE = np.sqrt(AllData["data"][S1][O][S2]["yerr"]["stat"][:,0]**2 + AllData["data"][S1][O][S2]["yerr"]["sys"][:,0]**2)

        for i, y in enumerate(TempPrediction[S1][O][S2]):
            axes[s1][s2].plot(DX, y, "b-", alpha=0.1, label="Posterior" if i==0 else "")
        axes[s1][s2].errorbar(DX, DY, yerr = DE, fmt="ro", label="Measurements")

#plt.tight_layout(True)
plt.tight_layout()
figure.savefig("plots/PredictedDesign.pdf", dpi = 192)
# figure


#--------------------------------------------------------------------------------------------------------------------------------
#TempPrediction = AllData["model"]

SystemCount = len(AllData["systems"])

figure, axes = plt.subplots(figsize = (15, 5 * SystemCount), ncols = 2, nrows = SystemCount, squeeze = False)

for s1 in range(0, SystemCount):
    for s2 in range(0, 1):
        axes[s1][s2].set_xlabel(r"Centrality")
        axes[s1][s2].set_ylabel(r"$\left<E_\mathrm{T}\right>$")

        S1 = AllData["systems"][s1]
        O  = AllData["observables"][S1][0][0]
        S2 = AllData["observables"][S1][0][1][s2]

        DX = AllData["data"][S1][O][S2]["x"]
        DY = AllData["data"][S1][O][S2]["y"]
        DE = np.sqrt(AllData["data"][S1][O][S2]["yerr"]["stat"][:,0]**2 + AllData["data"][S1][O][S2]["yerr"]["sys"][:,0]**2)

        for i, y in enumerate(TempPrediction[S1][O][S2]):
            axes[s1][s2].plot(DX, y, "b-", alpha=0.1, label="Posterior" if i==0 else "")
        axes[s1][s2].errorbar(DX, DY, yerr = DE, fmt="ro", label="Measurements")

#plt.tight_layout(True)
plt.tight_layout()
figure.savefig("plots/Design.pdf", dpi = 192)
# figure


#--------------------------------------------------------------------------------------------------------------------------------
# close all plots to save memory
plt.close("all")
