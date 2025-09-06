# GRPhIN: Graphlet Characterization of Regulatory and Physical Interaction Networks

GRPhIN (Graphlet Characterization of Regulatory and Physical Interaction Networks) is an algorithm for counting graphlets and the specific node positions within each graphlet (called orbits) in mixed regulatory and physical interaction networks. Graph representions of regulatory or physical interactions in isolation may obscure the complete functional context of a protein. PPI networks and GRNs do not exist separately; proteins are transcription factors, genes encode proteins, and physical and regulatory interactions mix and coexist forming their own distinct patterns. Graphlets are small, connected, induced subnetworks that describe patterns, local topologies, and organization in networks.

GRPhIN takes as input (1) an undirected PPI network and (2) a directed regulatory network and counts all mixed graphlets and their respective orbits ([Figure 6](https://github.com/Reed-CompBio/motifs/blob/main/Complete%20Graphlet%20%26%20Orbit%20Definitions.pdf)). GRPhIN provides additional functional context to the roles a protein may play beyond traditional isolated network types.

## Usage

1. Install and activate the GRPhIN conda environment in the root directory with `conda env create -f environment.yml` and `conda activate grphin`.
2. To run the GRPhIN orbit and graphlet counting algorithm on the _B. subtilis_ oxidative stress network, run the `grphin.py` script with

```
python3 grphin.py -u data/oxidative_stress/txid224308/stress_ppi.csv -d data/oxidative_stress/txid224308/stress_reg.csv -o out_dir
```

## Directories

- **`data/`** – Contains raw data files for case studies.
- **`final_output/`** – Contains output data files for case studies.

## File Descriptions

- **`countRandomizedNetworks.sh`** - Script to run GRPhIN in graphlets-only mode on a user-defined number of networks. Used to count graphlets in 1000 randomized networks for oxidative stress case studies.
- **`enrichment.py`** – Script for calculating graphlet enrichment statistics.
- **`environment.yml`** – Set up the conda environment with all dependencies required to run the project.
- **`graphlet_config.csv`** – File containing all graphlets and their orbits.
- **`grphin.py`** – Script for running the GRPhIN algorithm.
- **`generateNetworks.py`** – Script to generate randomized networks for significance testing.
- **`iterations_swaps.R`** - Script to generate plot showing thresholds for swaps for each species based on the percent randomization.
- **`iterations_swaps.txt`** - Dataset to calculate the percent randomization based on different numbers of swaps for all species.
- **`orbit_proteins.py`** – Script for finding protein identities and overrepresented orbits in GRPhIN results.
- **`pageRank.py`** – Script for running simple Random Walk with Restart algorithm to capture a subnetwork based on oxidative stress pathways.
- **`README.md`** – This file, providing documentation for the repository.
- **`significance.py`** – Script to calculate the significance of the appearance of mixed graphlets in oxidative stress pathways compared to random networks.
- **`StressGraphletSignificance.R`** – Script for analyzing results from network perturbation analysis.

## Analyses

1. Generate the stress response subnetwork using Random Walk with Restart (all species):
`python3 pageRank.py`
2. Generate 1000 randomized subnetworks with 1000 edge swaps based on stress subnetworks (all species): `python3 generateNetworks.py -s 1000 -i 1000`
3. Count graphlets (all species): `bash countRandomizedNets.sh -a -n 1000`
4. Calculate significance for graphlet occurences based on 1000 randomized networks (all species): `python3 significance.py -i 1000`