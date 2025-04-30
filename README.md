# GRPhIN: Graphlet Characterization of Regulatory and Physical Interaction Networks

GRPhIN (Graphlet Characterization of Regulatory and Physical Interaction Networks) is an algorithm for counting graphlets and the specific node positions within each graphlet (called orbits) in mixed regulatory and physical interaction networks. Graph representions of regulatory or physical interactions in isolation may obscure the complete functional context of a protein. PPI networks and GRNs do not exist separately; proteins are transcription factors, genes encode proteins, and physical and regulatory interactions mix and coexist forming their own distinct patterns. Graphlets are small, connected, induced subnetworks that describe patterns, local topologies, and organization in networks.

GRPhIN will be presented at Great Lakes Boinformatics Conference 2025. A preprint is available on bioRxiv:

[GRPhIN: Graphlet Characterization of Regulatory and Physical Interaction Networks](https://www.biorxiv.org/content/10.1101/2025.02.19.639099v1)<br>
Altaf Barelvi*, Oliver Anderson*, Anna Ritz<br>
\* *Equal contribution*

GRPhIN takes as input (1) an undirected PPI network and (2) a directed regulatory network and counts all mixed graphlets and their respective orbits ([Figure 6](https://github.com/Reed-CompBio/motifs/blob/main/Complete%20Graphlet%20%26%20Orbit%20Definitions.pdf)). GRPhIN provides additional functional context to the roles a protein may play beyond traditional isolated network types.

## Usage
Install and activate the GRPhIN conda environment in the root directory:
```conda env create -f environment.yml```
Once you have created the conda environment you can activate it (e.g. `conda activate grphin`).

To run the GRPhIN orbit and graphlet counting algorithm on the example networks, run the `grphin.py` script with `python grphin.py`.

```
usage: grphin.py [-h] -u UNDIRECTED -d DIRECTED -o OUTPUT_DIR
                 [-g GRAPHLETS_ONLY]

GRPhIN Algorithm

options:
  -h, --help            show this help message and exit
  -u, --undirected UNDIRECTED
                        Path to the undirected/PPI edges input file.
  -d, --directed DIRECTED
                        Path to the directed/Reg edges input file.
  -o, --output_dir OUTPUT_DIR
                        Path to the output directory.
  -g, --graphlets_only GRAPHLETS_ONLY
                        Run GRPhIN in graphlets only mode.
```

## Example Dataset


## Directories
- **`data/`** – Contains raw data files for case studies.
- **`final_output/`** – Contains output data files for case studies.
## File Descriptions
- **`countRandomizedNetworks.sh`** - Script to run GRPhIN in graphlets-only mode on a user-defined number of networks. Used to count graphlets in 1000 randomized networks for oxidative stress case studies.
- **`enrichment.py`** – Script for calculating graphlet enrichment statistics.
- **`environment.yml`** – Set up the conda environment with all dependencies required to run the project.
- **`grphin.py`** – Script for running the GRPhIN algorithm.
- **`generateNetworks.py`** – Script to generate randomized networks for significance testing.
- **`iterations_swaps.R`** - Script to generate plot showing thresholds for swaps for each species based on the percent randomization.
- **`iterations_swaps.txt`** - Dataset to calculate the percent randomization based on different numbers of swaps for all species. 
- **`orbit_proteins.py`** – Script for finding protein identities and overrepresented orbits in GRPhIN results.
- **`pageRank.py`** – Script for running simple Random Walk with Restart algorithm to capture a subnetwork based on oxidative stress pathways.
- **`README.md`** – This file, providing documentation for the repository.
- **`significance.py`** – Script to calculate the significance of the appearance of mixed graphlets in oxidative stress pathways compared to random networks.
