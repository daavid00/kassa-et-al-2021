## Description
This repository contains runscripts to produce the figures in paper [A]. They
have been successfully run in Linux and Mac OS (not tested in Windows). This
implementation is based on the OPM release 2021.04. New scripts have been added/
modified in opm-material, opm-models, and opm-simulators folders, corresponding
to the 'wettability' branch. All original tests in opm-models and opm-tests from
the release 2021.04 can be run using this implementation. The folder Tests 
include examples of scripts to run single simulations.

## Requirements
* [OPM](https://opm-project.org)
* [Python](https://www.python.org/downloads/)

## Python dependencies
* [numpy](https://numpy.org)
* [os](https://docs.python.org/3/library/os.html)
* [meshio](https://github.com/nschloe/meshio)
* [matplotlib](https://matplotlib.org)
* [pyvista](https://www.pyvista.org)

## Installation
* Clone all OPM modules from https://github.com/daavid00, check out the
"wettability" branch, and build all opm modules, specially 'wa' from opm-models
and 'flow' from opm-simulators (see/run the bash file).

`./buildopm.bash`
* Edit line 23 of the python scripts with the full path to the 'flow' and 'wa' 
executable respectively.

## Running the scripts
* From the terminal, e.g., for fig7:

`python3 fig7.py`

## Paper
* [A] Kassa, A.M., Gasda, S.E., Landa-Marbán, D., Sandve, T.H., Kumar, K., 2021.
The role of time-dependent wettability alteration in geological co2 storage.

## Contact
David Landa-Marbán (dmar@norceresearch.no).
