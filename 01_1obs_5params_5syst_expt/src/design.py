"""
Based on the JETSCAPE/STAT package, which is based on jbernhard/hic-param-est
See also http://qcd.phy.duke.edu/hic-param-est

Latin-hypercube parameter design.

"""

import itertools
import logging
from pathlib import Path
import re
import numpy as np


class Design:
    """
    Latin-hypercube model design.

    Public attributes:

    - ``type``: 'main' or 'validation'
    - ``keys``: list of parameter keys
    - ``labels``: list of parameter display labels (for TeX / matplotlib)
    - ``range``: list of parameter (min, max) tuples
    - ``min``, ``max``: numpy arrays of parameter min and max
    - ``ndim``: number of parameters (i.e. dimensions)
    - ``points``: list of design point names (formatted numbers)
    - ``array``: the actual design array

    The class also implicitly converts to a numpy array.

    This is probably the worst class in this project, and certainly the least
    generic.  It will probably need to be heavily edited for use in any other
    project, if not completely rewritten.

    """
    def __init__(self, validation=False):
        self.type = 'validation' if validation else 'main'

        self.keys, labels, self.range = [], [], []

        paramsfile = "Parameters.txt"

        with open(paramsfile, "rt") as f:
            for raw in f:
                line = raw.split("#")[0].strip()
                if (len(line) == 0): continue

                match  = re.fullmatch(r'([A-Za-z_][A-Za-z0-9_]*)\s+([^"\s]+|"[^"]+")\s+([0-9.Ee+-]+)\s+([0-9.Ee+-]+)', line)
                groups = match.groups() if match is not None else None
                if (groups is None) or (len(groups) != 4):
                    raise ValueError(f"unparseable line in '{paramsfile}': '{raw.strip()}'")

                self.keys.append(groups[0])
                labels.append(groups[1].strip('"'))
                self.range.append(tuple(map(float, groups[2:4])))

        # convert labels into TeX:
        #   - wrap normal text with \mathrm{}
        #   - escape spaces
        #   - surround with $$
        self.labels = [
            re.sub(r'({[A-Za-z]+})', r'\\mathrm\1', i)
            .replace(' ', r'\ ')
            .join('$$')
            for i in labels
        ]

        self.ndim = len(self.range)
        self.min, self.max = map(np.array, zip(*self.range))

        lhsmin = self.min.copy()
        #if not validation:
        #    for k, m in [
        #            ('fluct_std', 1e-3),
        #    ]:
        #        lhsmin[self.keys.index(k)] = m

        desfile = "processed/Design.dat"

        try:
            rawdes = np.loadtxt(desfile, ndmin=2)
            self.array = rawdes
            logging.info(f"[design.py: Design.__init__] loaded design array with shape {rawdes.shape} from file '{desfile}'")
        except Exception as exc:
            logging.warning(f"[design.py: Design.__init__] could not load design array from file '{desfile}': {exc}")
            raise

        if (rawdes is None) or (len(rawdes.shape) != 2) or (np.min(rawdes.shape) <= 0):
            raise ValueError(f"design file '{desfile}' is missing or invalid")

        if (rawdes.shape[1] != len(self.keys)):
            raise ValueError(f"design file '{desfile}' shape {rawdes.shape} conflicts with {len(self.keys)} parameter(s) in file '{paramsfile}'")

        npoints = len(self.array)

        # use padded numbers for design point names
        fmt = '{:03d}'
        self.points = [fmt.format(i) for i in range(npoints)]

        badptsfile = "badpoints.txt"

        badpts = np.loadtxt(badptsfile, dtype=int, ndmin=1)
        self.bad_points = []

        if (badpts is not None):
            badpts = np.atleast_1d(badpts).ravel()

            if (len(badpts.shape) == 1) and (len(badpts) > 0):
                self.bad_points = badpts.tolist()

        if (len(self.bad_points) > 0) and ((np.min(self.bad_points) < 0) or (np.max(self.bad_points) >= npoints)):
            logging.warning(f"bad points in file '{badptsfile}' exceed {npoints} point(s) in design file '{desfile}'; ensure that the bad points are zero-based indices of points in the full design")

        # drop bad design points
        keep = [n not in self.bad_points for n in range(npoints)]
        self.array = self.array[keep]

    def __array__(self):
        return self.array
