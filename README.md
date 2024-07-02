# Carrier transport properties in AlGaN/AlGaN/AlN HEMT
This repository contains scripts for simulating the 2D electron and hole gas of an AlGan/AlGaN/AlN High-Electron Mobility Transistor (HEMT). The simulations are performed using [nextnano software](https://www.nextnano.com). The workflow includes:

* __2DEG density and 2DHG density Nextnano++ 1D Simulations__
* __Sweep the Device Parameters for Optimmizing 2DEG density and 2DHG density in 1D Simulations:__ Sweep were perfomed for barrier thickness, channel thickness, barrier AlGaN composition, channel AlGaN composition, Schottky barrier height, device temperature,
bandgap bowing parameter for AlGaN alloy, AlGaN pyroelctric bowing parameter, composition gradient in the interfaces.
* __Post-Processing__

All three steps are seamlessly integrated into a single workflow using the [nextnanopy](https://github.com/nextnanopy) Python package.

The repository is structured as follows:
 - [SCRIPTs](SCRIPTs) folder: Notebooks for the simulations and post processing.
 - [INPUTs](SCRIPTs) folder: The nextnano input files.
 - [DATAs](DATAs) folder: The post-processed data


## Installation
### 1.1 Requirements
#### 1.1.1 Python dependency
  - python >= 3.10, pathlib, os, sys, shutil, importlib, datetime,
  - nextnanopy, numpy, matplotlib
  - pandas, openpyxl


### 1.2 Scripts executation or run command
1. Untar nextnano_sim.tar.bz2:
```
  $ tar -xvf nextnano_sim.tar.bz2
  $ cd nextnano_sim/SCRIPTS
```
3. Install conda/mamba environment (if not already installed):
```
  $ <install anaconda or miniforge>
  $ <create conda environment>
  $ conda create -n nextnano python
  $ conda activate nextnano
  $ pip install nextnanopy yaml numpy pandas openpyxl matplotlib
```
3. Update required properties in notebook and run:
```
  $ Open 'NoteBooks/.../....ipynb' in jupyter lab/notebook. Update Cell 4 as needed.
```
7. Done. View results:
```
  $ <once everything is finished open NoteBooks/.../....ipynb for the outputs>
  $ jupyter lab NoteBooks/.../....ipynb
```

### 1.3 Folder_Name: Scripts_Name: Description
- `BaseFunctions:`
  - `general_functions.py`: Contains general functions definitions such as delete folder, create folder, link folder, parsers etc.

- `ClusterSubmitScripts:`
  - `submit_script_v0.pbs:` PBS bash script to submit job in cluster queuing system (version: 0).

- `NoteBooks:`
  - `first_principles_simulation.ipynb:` Ipython notebook. This notebook will be executed in the submit_script_v0.pbs.

- `PostProcessing:`
  - `BandUnfold:` Contains helper python scripts for band-unfolding tasks.
  - `MobilityCalculation:` Contains helper python scripts for mobility calculations.

- `PreProcessing:`
   - `SQS:` Contains helper python scripts to generate SQS structures.
      - `general_functions.py:` Contains general functions definitions for SQS generation using 'sqsgenerator' python module.
      - `aux_functions.py:` Contains auxiliary function definitions for SQS generation using 'sqsgenerator' python module.

- `ReadMeFiles:` Contains helper text files that will be symlinked to appropriate simulation folders.

- `VaspPotcars:` Contains VASP pseudopotential/PAW potential/vdw files.
  - `potpaw_PBE:` VASP PAW potentials (PBE.54). This is symlinked to original PAW potentials path.

- `VaspScripts:`
  - `general_functions.py:` Contains general functions definitions for VASP files generation using 'ASE-VASP' python module.
  - `aux_functions.py:` Contains auxiliary function definitions for VASP files generation using 'ASE-VASP' python module.

- `helper.txt:` Text file describing the folder descriptions and details of the scripts in this folder.

## References and citations
If you find this repository usefull for your work, we would appreciate if you cite the appropriate project specific references as highlighted in the previous section along with the software/packages used. Here is the complete list of references [(bibliography file)](docs/REFERENCES.md):

>> [1] TBA

>>

Here is the bibliography file for your convenience: [bibliography file](docs/REFERENCES.md)

## Contact
Contact us: [Email developer/maintainer team](mailto:stefan.schulz@tyndall.ie,badal.mondal@tyndall.ie,badalmondal.chembgc@gmail.com)

If you have new ideas or would like to contribute to further development of this repository or request new functionality, please get in touch with [us](mailto:stefan.schulz@tyndall.ie,badal.mondal@tyndall.ie,badalmondal.chembgc@gmail.com) or open a pull request. We are looking forward to support your request and discuss ideas.

## Acknowledgements

Developer: [Badal Mondal](https://github.com/bmondal94)

Maintainer: [Badal Mondal](https://github.com/bmondal94)

Project collaborators: 

* [Dr. Stefan Schulz, Tyndall National Institute](https://www.tyndall.ie/people/stefan-schulz/)

## Contributors

TBA

## Version release history
Latest release: v0.1.0

Chekout out version release history [here](docs/RELEASE.md) for the full list of updates and upgrades.

# License
* [GNU General Public License v3.0](LICENSE)
