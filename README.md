# Carrier transport properties in AlGaN/AlGaN/AlN HEMT
This repository contains scripts for simulating the 2D electron gas (2DEG) and 2D hole gas (2DHG) of AlGan/AlGaN/AlN High-Electron Mobility Transistors (HEMTs). The simulations are performed using [nextnano software](https://www.nextnano.com). The workflow includes:

* __2DEG density and 2DHG density Nextnano++ Simulations__
* __Sweep the Device Parameters for Optimmizing 2DEG density and 2DHG density in Simulations:__ Sweeps were perfomed for barrier thickness, channel thickness, barrier AlGaN composition, channel AlGaN composition, Schottky barrier height, device temperature, bandgap bowing parameter for AlGaN alloy, AlGaN pyroelctric bowing parameter, composition gradient in the interfaces.
* __Post-Processing__

All three steps are seamlessly integrated into a single workflow using the [nextnanopy](https://github.com/nextnanopy) Python package.

## Projects
This repository is built on the following projects:

* __2DEG_DensityMobilityInterplay__ : Mondal et. al.,  TBA

## Folders
 - [SCRIPTs](SCRIPTs) : Jupyter Notebooks and executable python scripts for the running simulations and post processing.
 - [INPUTs](SCRIPTs)  : The nextnano input files.
 - [DATAs](DATAs)    : The post-processed data.
 - [docs](docs) : Release history of this repo.

__Note:__ Each folder contains special `helper.txt` file to guide you through details of the corresponding folder and file it contains.

## Installation
### 1.1 Requirements
#### 1.1.1 Python dependency
  - python >= 3.10, pathlib, os, sys, shutil, importlib, datetime,
  - nextnanopy, numpy, matplotlib
  - pandas, openpyxl

### 1.2 Scripts executation or run command
1. Untar nextnano_sim.tar.bz2:
```
  $ tar -xvf SchrodingerPoissonSimulations.tar.bz2
  $ cd SchrodingerPoissonSimulations/SCRIPTS/AlGaN_HEMT/Notebooks/AlGaN_AlGaN_AlN_HEMT/
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
  $ Open 'NoteBooks/.../....ipynb' in jupyter lab/notebook. Update Cells as needed.
```
7. Done. View results:
```
  $ <once everything is finished open NoteBooks/.../....ipynb for the outputs>
  $ jupyter lab NoteBooks/.../....ipynb
```

## References and citations
If you find this repository usefull for your work, we would appreciate if you cite the appropriate project specific references along with the software/packages used. Here is the complete list of references [(bibliography file)](docs/REFERENCES.md):

>> [1] Mondal et. al., TBA

>>

Here is the bibliography file for your convenience: [bibliography file](docs/REFERENCES.md)

## Contact
Contact us: [Email developer/maintainer team](mailto:stefan.schulz@tyndall.ie,badal.mondal@tyndall.ie,badalmondal.chembgc@gmail.com)

If you have new ideas or would like to contribute to further development of this repository or request new functionality, please get in touch with [us](mailto:stefan.schulz@tyndall.ie,badal.mondal@tyndall.ie,badalmondal.chembgc@gmail.com) or open a pull request. We are looking forward to support your request and discuss ideas.

## Acknowledgements

Developer: [Badal Mondal](https://github.com/bmondal94)

Maintainer: [Badal Mondal](https://github.com/bmondal94)

## Contributors

Complete list of contributors in this repository can be found [here](https://github.com/SemiconductorTransport/SchrodingerPoissonSimulations/graphs/contributors). We sincerely thank each and every contributor for their valuable input and support.

## Version release history
Latest release: v1.0.0

Chekout out version release history [here](docs/RELEASE.md) for the full list of updates and upgrades.

# License
* [GNU General Public License v3.0](LICENSE)
