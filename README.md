# HDock-Multimer

`HDock-Multimer (HDM)` is a computational tool developed for predicting the structure of large protein assemblies. 

With the stoichiometry file, AlphaFold-predicted monomer and subcomponent structures as input, `HDM` automatically selects suitable modeling strategies, including asymmetric docking, symmetric docking and assembly. The final output is a set of ranked predicted structures of the full complex.

## Demo


# Installation
## System Requirements
### Hardware requirements
For running the `HDM` program, a standard computer with a standard CPU is sufficient.

For AF2-based prediction of monomers and subcomponents, an NVIDIA GPU with >= 40GB VRAM is required considering the number of residues in the subcomponents to be predicted. While the multi-core CPU and sufficient RAM/disk are helpful for preprocessing and I/O. CPU-only AF2 inference is strongly discouraged due to unreasonable time consumption.

### Software requirements
`HDM` is supported for *Linux*, the code has been tested on CentOS 6.

In order to run the `HDM` program, several third-party packages/programs are required as follows.

#### Python package requirements
The supporting scripts depend on these Python3 packages, which can be installed using pip:
```
numpy
pandas
scipy
openmm
```
In this repository, we use **OpenMM** (instead of Amber14 reported in the manuscript) to perform structural relaxation, because OpenMM is open-source and straightforward to set up. Please note that minor numerical differences from Amber may occur due to different force fields and integrator implementations.

#### JSON processing requirements
This project depends on `jq` to analyze the stoichiometry files in JSON format. We recommend downloading the stand-alone binary (e.g., `jq-linux64`) from the official releases: https://github.com/jqlang/jq/releases/. Make it globally accessible by running:
```bash
mv jq-linux64 ~/bin/jq
chmod +x ~/bin/jq
export PATH=$HOME/bin:$PATH
```

#### Docking tool requirements
This project depends on two external docking toolkits: 

**HDOCKlite** for pairwise docking: \
Download: http://huanglab.phys.hust.edu.cn/software/hdocklite/. \
Included two pre-complied executables: `hdock` and `createpl`.

**HSYMDOCKlite** for cyclic and dihedral symmetric docking: \
Download: http://huanglab.phys.hust.edu.cn/software/hsymdock/.
Included several pre-complied executables: `chdock`, `compcn`, `dhdock`, `dhdock.sh` and `compdn`.

Add their directories to your `PATH` (replace with your actual install paths):
```bash
export PATH=/path/to/your/HDOCKlite:$PATH
export PATH=/path/to/your/HSYMDOCKlite:$PATH
```

#### Structural analysis requirements
This project uses two external tools to extract information from AF-predicted models. You should prepare a C++ compiler for compiling them, and make sure the compiled executables are in your path.

**MMalign** for extracting pairwise transformations: \
Download: https://zhanggroup.org/MM-align/.
You can compile the latest C++ source code `MMalign.cpp` by running:
```bash
g++ -static -O3 -ffast-math -o MMalign MMalign.cpp
mv MMalign ~/bin/
export PATH=$HOME/bin:$PATH
```

**STRIDE** for analyzing protein secondary structures: \
Download: https://webclu.bio.wzw.tum.de/stride/.
You can compile the stride stand-alone program source code by running:
```bash
tar -xzf stride.tar.gz
make
mv stride ~/bin/
export PATH=$HOME/bin:$PATH
```
