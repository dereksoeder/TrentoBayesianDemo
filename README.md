# T<sub>R</sub>ENTo Bayesian Analysis Demo

Based on <https://github.com/JETSCAPE/STAT> / <https://github.com/jbernhard/hic-param-est>.  Uses the [T<sub>R</sub>ENTo](https://github.com/Duke-QCD/trento) heavy-ion collision initial conditions model.

## Setup

Note that a Linux environment is assumed.

- Make a local copy:  `git clone --recursive https://github.com/dereksoeder/TrentoBayesianDemo`
- Build T<sub>R</sub>ENTo:
  - `mkdir trento/build`
  - `cd trento/build`
  - `cmake ..`
  - `make`
  - T<sub>R</sub>ENTo requires the [Boost](https://www.boost.org/) C++ libraries.  If CMake isn't able to find the libraries automatically, you might need to:
    - Load a Boost [environment module](https://en.wikipedia.org/wiki/Environment_Modules_%28software%29) that already exists on the system: `module load boost`  (You might need to specify a version; e.g.: `module load boost/1.76.0`)
    - Tell CMake the path to Boost; e.g.: `cmake .. -DBOOST_ROOT=XXX/boost_1_76_0 -DBoost_INCLUDE_DIR=XXX/boost_1_76_0` (where `XXX` is something you'll need to fill in)
    - Or [build](https://www.boost.org/more/getting_started/index.html) your own copy of Boost
- From one of the demonstration directories, run `RunMe.sh`

## Demonstrations

- `01_1obs_5params_5syst_expt`
  - **Observable:**  Transverse energy (\[_dE<sub>T</sub>_/_d&eta;_\]\(_&eta;_=0\)), for which T<sub>R</sub>ENTo multiplicity is used as a proxy, in various centrality bins per system
  - **Parameters:**  Normalization at 62.4 GeV (_N_<sub>62.4&nbsp;GeV</sub>), normalization at 200 GeV (_N_<sub>200&nbsp;GeV</sub>), normalization at 2.76 TeV (_N_<sub>2.76&nbsp;TeV</sub>), nucleon width (_w_), fluctuation (_k_)
  - **Collision systems:**  Au-Au @ 62.4 GeV, Au-Au @ 200 GeV, Cu-Cu @ 62.4 GeV, Cu-Cu @ 200 GeV, Pb-Pb @ 2.76 TeV
  - **Comparison:**  Model-to-data
  - **Notes:**  Takes a couple hours to run.  To change the design, replace `processed/Design.dat` and modify `Parameters.txt`.

## Experimental data

Experimental data gathered by @keweiyao from [HEPData](https://www.hepdata.net/).  References:

- Transverse energy at midrapidity (\[_dE<sub>T</sub>_/_d&eta;_\]\(_&eta;_=0\))
  - Au-Au @ 62.4 GeV (PHENIX): [Phys. Rev. C __93__ 024901](https://doi.org/10.1103/PhysRevC.93.024901) Table VII, [arXiv:1509.06727](https://arxiv.org/abs/1509.06727) Table VII, HEPData [INSPIRE ID 1394433](https://www.hepdata.net/record/ins1394433) [Table 11](https://www.hepdata.net/record/96609)
  - Au-Au @ 200 GeV (PHENIX): _Ibid._, HEPData INSPIRE ID 1394433 [Table 7](https://www.hepdata.net/record/96602)
  - Cu-Cu @ 62.4 GeV (PHENIX): _Ibid._, HEPData INSPIRE ID 1394433 [Table 25](https://www.hepdata.net/record/96622)
  - Cu-Cu @ 200 GeV (PHENIX): _Ibid._, HEPData INSPIRE ID 1394433 [Table 23](https://www.hepdata.net/record/96620)
  - Pb-Pb @ 2.76 TeV (ALICE): [Phys. Rev. C __94__ 034903](https://doi.org/10.1103/PhysRevC.94.034903), [arXiv:1603.04775](https://arxiv.org/abs/1603.04775), HEPData INSPIRE ID 1427723 [Table 1](https://www.hepdata.net/record/73995)
