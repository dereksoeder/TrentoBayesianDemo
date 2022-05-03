import os
from os import path
from sys import stderr
import numpy as np


predclass = "main"  # alternatively could be "validation"


centsdict   = {}  # sorted array of (centlo, centhi) pairs, per observable, per system (first level is by system, second is by observable)
exptoutdict = {}  # processed experimental data file name, per observable, per system (first level is by system, second is by observable)


# Read design file and list of bad points

desfile = "Design.dat"
despath = path.join("processed", desfile)

design = np.loadtxt(despath, ndmin=2)
if (design is None) or (len(design.shape) != 2) or (np.min(design.shape) <= 0):
    raise ValueError(f"design file '{despath}' is missing or invalid")

badptsfile = "badpoints.txt"

bad_points = np.loadtxt(badptsfile, dtype=int, ndmin=1) if path.isfile(badptsfile) else None

if (bad_points is not None):
    bad_points = np.atleast_1d(bad_points).ravel()

if (bad_points is None) or (len(bad_points.shape) != 1):
    print(f"WARNING: bad points file '{badptsfile}' is missing or invalid", file=stderr)
    bad_points = []
else:
    bad_points = bad_points.tolist()

if (len(bad_points) > 0) and ((np.min(bad_points) < 0) or (np.max(bad_points) >= len(design))):
    print(f"WARNING: bad points in file '{badptsfile}' exceed {len(design)} point(s) in design file '{despath}'; ensure that the bad points are zero-based indices of points in the full design", file=stderr)

good_points = [_ for _ in range(len(design)) if _ not in bad_points]


# Process experimental data

def ParseExpt_genericETmid(input):
    # columns are:  centbinmidpct centbinlopct centbinhipct ET [staterr+ staterr-] syserr+ syserr-
    if (input.shape[1] < 6) or (((input.shape[1] - 6) % 2) != 0):
        return None

    if np.any(input[:,1] < 0.) or np.any(input[:,0] < input[:,1]) or np.any(input[:,2] < input[:,0]) or np.any(input[:,2] > 100.) or np.any(input[:,3] < 0.):
        return None

    yerrs = np.abs(input[:,4:].reshape((len(input), -1, 2))).max(axis=2)  # take the greater of each pair \of +/- error columns
    yerrs = np.sqrt(np.square(yerrs).sum(axis=1))                         # and add them all in quadrature

    return input[:,1:3] / 100., input[:,3], yerrs


for exptfile in os.listdir("expt"):
    exptpath = path.join("expt", exptfile)

    sysobs, exptext = path.splitext(exptfile)
    if (exptext != ".txt"): continue

    sysobs = sysobs.split("-")
    if (len(sysobs) != 2): continue

    systag, obsname = sysobs

    if (obsname != "ETmid"):
        print(f"WARNING: skipping experimental data file '{exptpath}' with unrecognized observable name '{obsname}'", file=stderr)
        continue

    exptdata = np.loadtxt(exptpath, ndmin=2)

    if (exptdata is None) or (len(exptdata.shape) != 2) or (len(exptdata) == 0) or (exptdata.shape[1] < 4):  # expect minimum 4 columns: centlo, centhi, y, yerr
        print(f"WARNING: skipping invalid experimental data file '{exptpath}'", file=stderr)
        continue

    parsedexpt = ParseExpt_genericETmid(exptdata)

    if (parsedexpt is None):
        print(f"WARNING: skipping unparseable experimental data file '{exptpath}'", file=stderr)
        continue

    cents, ys, yerrs = parsedexpt

    sortorder = np.lexsort(np.flip(cents, axis=1).T)  # numpy is awkward about sorting multidimensional arrays; this does what I wish argsort would do
    cents, ys, yerrs = (_[sortorder] for _ in (cents, ys, yerrs))  # sort centrality bins in increasing order

    centsbysyst = centsdict.setdefault(systag, {})

    if (obsname in centsbysyst):
        raise KeyError(f"experimental data file '{exptpath}' duplicates system '{systag}' observable '{obsname}'")

    centsdict.setdefault(systag, {})[obsname] = cents

    exptout = f"Data-{systag}-{obsname}.dat"
    exptoutdict.setdefault(systag, {})[obsname] = exptout

    with open(path.join("processed", exptout), "wt") as fout:
        print(
f"""# Version 1.0
# DOI
# Source
# System {systag}
# Centrality {100.*cents.min():g} {100.*cents.max():g}
# XY Centrality {obsname}
# Label x y stat,low stat,high sys,low sys,high""",
            file=fout)

        for x, y, ystaterr, ysyserr in zip(cents.mean(axis=1), ys, yerrs/np.sqrt(2.), yerrs/np.sqrt(2.)):  # we didn't maintain stat and syst errors separately, so we'll just split the combined error equally between the two, since the file format requires both
            print(f"    {x:g} {y:g} {ystaterr:g} {ystaterr:g} {ysyserr:g} {ysyserr:g}", file=fout)


# Process model data

def ParseModel_TrentoETmid(input, centranges):
    # columns are (see also http://qcd.phy.duke.edu/trento/usage.html):  [0]=event_num  [1]=impact_parameter  [2]=num_participants  [3]=multiplicity  [4]=e2  [5]=e3  [6]=e4  [7]=e5

    if (input.shape[1] != 8):
        return None

    input = input[(-input[:,3]).argsort()]  # sort descending by multiplicity for centrality binning
    mults = input[:,3]                      # we could've just reverse-sorted `mults`, but sorting `input` might be more convenient when implementing other observables

    predvalues = []

    for centlo, centhi in centranges:
        ilo, ihi = ( int(round(len(mults) * _)) for _ in (centlo, centhi) )

        if (ilo >= ihi):
            raise ValueError(f"centrality bin [{centlo:g} .. {centhi:g}) would contain {ihi-ilo} event(s)")

        predvalues.append(mults[ilo:ihi].mean())

    return predvalues


predsdict = {}  # list of predictions (model data observables) parallel with list of centrality bins, per design point, per system (first level of dict is by system, second level is by observable, third is by design point)

for eventsfile in os.listdir("data"):
    eventspath = path.join("data", eventsfile)

    eventsdesc, eventsext = path.splitext(eventsfile)
    if (eventsext != ".txt"): continue

    eventsdesc = eventsdesc.split("-")
    if (len(eventsdesc) != 3) or (eventsdesc[0] != "trento") or not eventsdesc[2].isnumeric():
        continue

    systag, despt = eventsdesc[1], int(eventsdesc[2])
    if systag not in centsdict:
        print(f"WARNING: skipping model data file '{eventspath}' with unrecognized system '{systag}'", file=stderr)
        continue

    if (despt in bad_points):
        continue

    events = np.loadtxt(eventspath, ndmin=2)

    if (events is None) or (len(events.shape) != 2) or (np.min(events.shape) <= 0):
        print(f"WARNING: skipping invalid model data file '{eventspath}'", file=stderr)
        continue

    centsforsys = centsdict[systag]

    for obsname, centsforobs in centsforsys.items():
        if (obsname != "ETmid"):
            print(f"WARNING: skipping unrecognized observable '{obsname}' for model data file '{eventspath}'", file=stderr)
            continue

        preds = ParseModel_TrentoETmid(events, centsforobs)

        if (preds is None):
            print(f"WARNING: could not compute observable '{obsname}' for model data file '{eventspath}'", file=stderr)
            continue

        predsforobs = predsdict.setdefault(systag, {}).setdefault(obsname, {})

        if (despt in predsforobs):
            raise KeyError(f"model data file '{eventspath}' duplicates system '{systag}' design point #{despt}")

        predsforobs[despt] = preds


predssystems = set(predsdict.keys()) & set(centsdict.keys())  # only systems that appear in both experiment and model data

if (predssystems != (set(predsdict.keys()) | set(centsdict.keys()))):
    print(f"WARNING: model data systems {list(predsdict.keys())} do not match experiment data systems {list(centsdict.keys())}", file=stderr)

if (len(predssystems) == 0):
    raise ValueError("no model data to process!")

allobs = set()

for systag in predssystems:
    allobs |= set(centsdict[systag].keys()) | set(predsdict[systag].keys())

for systag in predssystems:
    if (set(centsdict[systag].keys()) != allobs):
        raise ValueError(f"experiment data for system '{systag}' is missing observable(s) {allobs - set(centsdict[systag].keys())}")

    if (set(predsdict[systag].keys()) != allobs):
        raise ValueError(f"model data for system '{systag}' is missing observable(s) {allobs - set(predsdict[systag].keys())}")

for systag in predssystems:
    predsforsys = predsdict[systag]

    for obsname, predsforobs in predsforsys.items():
        if (set(predsforobs.keys()) != set(good_points)):
            raise ValueError(f"model data for system '{systag}' observable '{obsname}' is missing design point(s) {set(good_points) - set(predsforobs.keys())}")

        with open(f"processed/Prediction-{predclass}-{systag}-{obsname}.dat", "wt") as fout:
            print(
f"""# Version 1.0
# Data {exptoutdict[systag][obsname]}
# Design {desfile}""",
                file=fout)

            predvalues = np.array([ predsforobs[despt] for despt in good_points ])  # each row is a design point, each column is a centrality bins
            for predsbin in predvalues.T:                                           # but file format wants each row to be a centrality bin, each column a design point
                print(" ".join(map(str, predsbin)), file=fout)
