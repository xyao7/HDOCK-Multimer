# HDock-Multimer

`HDock-Multimer (HDM)` is a computational tool developed for predicting the structure of large protein assemblies. 

With the stoichiometry file, AlphaFold-predicted monomer and subcomponent structures as input, `HDM` automatically selects suitable modeling strategies, including asymmetric docking, symmetric docking and assembly. The final output is a set of ranked predicted structures of the full complex. 

The `HDM` pipeline consists of 4 stages:
1. Defining subunits in stoichiometry file.
2. Predicting monomer and subcomponent structures using AF2/AFM.
3. Determining suitable modeling strategies.
4. Performing selected modeling strategies to model the full complex.


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
In this repository, we use **OpenMM** (instead of Amber14 reported in the manuscript) to perform structural relaxation, because OpenMM is open-source and straightforward to set up.

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



# Running HDM
## Stage 1 - Defining subunits in the stoichiometry file
In this stage, the `stoi.json` file is used to specify the stoichiometry information of the target complex and the subunit definitions. Generally, each subunit should be a complete chain in the complex. If a chain is very long (e.g., >1500 residues), you may split it into regions by domain boundaries predicted with tools like [IUPred3](https://iupred3.elte.hu/), or into subunits of practical lengths based on your compute constraints.

Each subunit definition should contain these fields:
- `entry_id`, `entity_id`, `uniprot_id`: identifiers for complex and chains.
- `asym_ids`: chain names corresponding to the subunit as well as the stoichiometry.
- `start_res`, `end_res`, `sequence`, `length`: residue range and sequence.
- `subunits`: information of the subunits obtained by splitting the long chain. **Notice** that the chain names of subunits obtained by splitting the long chain must differ from those of the complete chains.

The subunits corresponding to the complete chain and the split regions of chain are defined in the following format:
```
{
  "entry_id": "example",
  "entity_id": "1",
  "uniprot_id": "uniprot_id",
  "asym_ids": ["A", "B", "C"],
  "start_res": 1,
  "end_res": 100,
  "sequence": "MSTAKLVKSKATNLLYTRNDVSDSEKKATVELLNRQVIQFIDLSLITKQAHWNMRGANFIAVHEMLDGFRTALIDHLDTMAERAVQLGGVALGTTQVINS",
  "length": 100
}
```
```
{
  "entry_id": "example",
  "entity_id": "1", 
  "uniprot_id": "uniprot_id",
  "asym_ids": ["A", "B"],
  "start_res": 1,
  "end_res": 200,
  "sequence": "MNPHDLEWLNRIGERKDIMLAVLLLAVVFMMVLPLPPLVLDILIAVNMTISVVLLMIAIYINSPLQFSAFPAVLLVTTLFRLALSVSTTRMILLQADAGQIVYTFGNFVVGGNFIVGIVIFLIITIVQFLVITKGSERVAEVSARFSLDAMPGKQMSIDGDMRAGVIDVNEARERRATIEKESQMFGSMDGAMKFVKGDA",
  "length": 200, 
  "subunits": [
    {"subunit_id": "1-0", "asym_ids": ["A", "B"], "start_res": 1, "end_res": 100, "sequence": "MNPHDLEWLNRIGERKDIMLAVLLLAVVFMMVLPLPPLVLDILIAVNMTISVVLLMIAIYINSPLQFSAFPAVLLVTTLFRLALSVSTTRMILLQADAGQ", "length": 100},
    {"subunit_id": "1-1", "asym_ids": ["a", "b"], "start_res": 101, "end_res": 200, "sequence": "IVYTFGNFVVGGNFIVGIVIFLIITIVQFLVITKGSERVAEVSARFSLDAMPGKQMSIDGDMRAGVIDVNEARERRATIEKESQMFGSMDGAMKFVKGDA", "length": 100}
  ]
}
```
The `stoi.json` is a JSON list of subunits, which format is as follows:
```
[
    {"entry_id": "example", "entity_id": "1", "uniprot_id": "uniprot_id1", "asym_ids": ["A", "B"], "start_res": 1, "end_res": 122, "sequence": "MSRIITAPHIGIEKLSAISLEELSCGLPDRYALPPDGHPVEPHLERLYPTAQSKRSLWDFASPGYTFHGLHRAQDYRRELDTLQSLLTTSQSSELQAAAALLKCQQDDDRLLQIILNLLHKV", "length": 122},
    {"entry_id": "example", "entity_id": "2", "uniprot_id": "uniprot_id2", "asym_ids": ["C"], "start_res": 1, "end_res": 114, "sequence": "MNITLTKRQQEFLLLNGWLQLQCGHAERACILLDALLTLNPEHLAGRRCRLVALLNNNQGERAEKEAQWLISHDPLQAGNWLCLSRAQQLNGDLDKARHAYQHYLELKDHNESP", "length": 114},
    {"entry_id": "example", "entity_id": "3", "uniprot_id": "uniprot_id3", "asym_ids": ["D", "E"], "start_res": 1, "end_res": 200, "sequence": "MNPHDLEWLNRIGERKDIMLAVLLLAVVFMMVLPLPPLVLDILIAVNMTISVVLLMIAIYINSPLQFSAFPAVLLVTTLFRLALSVSTTRMILLQADAGQIVYTFGNFVVGGNFIVGIVIFLIITIVQFLVITKGSERVAEVSARFSLDAMPGKQMSIDGDMRAGVIDVNEARERRATIEKESQMFGSMDGAMKFVKGDA", "length": 200, 
    "subunits": [
      {"subunit_id": "3-0", "asym_ids": ["D", "E"], "start_res": 1, "end_res": 100, "sequence": "MNPHDLEWLNRIGERKDIMLAVLLLAVVFMMVLPLPPLVLDILIAVNMTISVVLLMIAIYINSPLQFSAFPAVLLVTTLFRLALSVSTTRMILLQADAGQ", "length": 100},
      {"subunit_id": "3-1", "asym_ids": ["d", "e"], "start_res": 101, "end_res": 200, "sequence": "IVYTFGNFVVGGNFIVGIVIFLIITIVQFLVITKGSERVAEVSARFSLDAMPGKQMSIDGDMRAGVIDVNEARERRATIEKESQMFGSMDGAMKFVKGDA", "length": 100}
    ]}
]
```
Example above describes a complex with 5 chains (A, B, C, D, E), in which A & B are identical chains (length 122), as well as D & E (length 200). Each of D & E is split into two subunits. 

## Stage 2 - Predicting monomer and subcomponent structures
In this stage, we run AF2/AFM to predict models for monomers (each monomer is a full chain in the complex) and subcomponents.

Using this script to generate fasta files for both:
```
python scripts/prepare_fasta.py -s <stoi.json> -p <complex_id> -n 3
```
Notice that the `<complex_id>` specifies the identifier of the target complex, and `-n 3` defines subcomponent size (3 for trimer).

This command will create two output directories in the same location as `stoi.json`: `mono_fastas/` for storing FASTA files for individual monomer chains, and `subcomp_fastas/` for storing FASTA files for multi-chain subcomponents. All generated FASTA files can be used as input for AF2/AFM structure prediction.

You can use [ColabFold](https://github.com/sokrypton/ColabFold) to run AF2/AFM using a Google Colab Notebook. You can also download the code from [AlphaFold](https://github.com/jcheongs/alphafold-multimer) and run predictions locally (which requires a GPU).

Disordered terminal regions of predicted monomers may interfere with the docking process. We recommend using the following script to preprocess the monomers (requires `STRIDE`):
```
python scripts/filter_monofiles.py <stoi.json> <mono_dir>
```
Here, `mono_dir` is the directory storing AF2-predicted monomer PDB files. **Notice** each file must be named as `chain_id.pdb` (e.g., `A.pdb`), where `chain_id` corresponds to the chain ID of the monomer. The preprocessed PDB files (using the same naming convention) are saved in the same location as `stoi.json`, as well as a summary file `res_indices_D.txt` listing the start/end residue indices for each monoer after filtering.

## Stage 3 - Determining suitable modeling strategies
In this stage, we will select appropriate modeling strategies based on the complex homology and pairwise docking results. 

**Available Modeling Strategies:**
1.**Strategy 1**: asymmetric docking (i.e., multi-body docking).
2.**Strategy 2**: symmetric docking (a combination of docking methods for different symmetry types, covering cyclic, dihedral, tetrahedral and octahedral symmetries).
3.**Strategy 3**: assembly.

**Strategy Selection Rules:**
-**Homomeric complexes:** 1, 2 & 3.
-**Heteromeric complexes with <5 unique chains:** 1 & 3.
-**Heteromeric complexes with >=5 unique chains:** 3. 

Using this script to select appropriate strategies (requires `HDOCKlite` and `HSYMDOCKlite`):
```
python scripts/check_build_type.py <stoi.json>
```
The command will output a file `build_types.txt` recording all selected modeling strategies, and create different directories for storing the results of each strategy. The directories are as follows:
-`asym/`: asymmetric docking.
-`c<n>/`: cyclic symmetric docking (e.g., `c5/` for C5 symmetry).
-`d<n>/`: dihedral symmetric docking (e.g., `d5/` for D5 symmetry.)
-`cubicT/`: tetrahedral symmetric docking.
-`cubicO/`: octahedral symmetric docking.
-`assemble/`: assembly.

## Stage 4 - Modeling, scoring and ranking models
In this stage, we will perform all the modeling strategies selected in the previous stage. 