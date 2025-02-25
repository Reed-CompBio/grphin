# GRPhIN: Graphlet Characterization of Regulatory and Physical Interaction Networks

GRPhIN (Graphlet Characterization of Regulatory and Physical Interaction Networks) is an algorithm for counting graphlets and the specific node positions within each graphlet (called orbits) in mixed regulatory and physical interaction networks. Graph representions of regulatory or physical interactions in isolation may obscure the complete functional context of a protein. PPI networks and GRNs do not exist separately; proteins are transcription factors, genes encode proteins, and physical and regulatory interactions mix and coexist forming their own distinct patterns. Graphlets are small, connected, induced subnetworks that describe patterns, local topologies, and organization in networks.

GRPhIN takes as input (1) an undirected PPI network and (2) a directed regulatory network and counts all mixed graphlets and their respective orbits ([Figure 6](https://github.com/Reed-CompBio/motifs/blob/main/Complete%20Graphlet%20%26%20Orbit%20Definitions.pdf)). GRPhIN provides additional functional context to the roles a protein may play beyond traditional isolated network types.

## Usage
1. Activate the GRPhIN conda environment in the root directory with `conda activate grphin`.
2. To run the GRPhIN orbit and graphlet counting algorithm on the example networks, run the `grphin.py` script with `python3 grphin.py`.
3. Follow the menu options for your organism of interest.

## File Descriptions

- **`README.md`** – This file, providing documentation for the repository.
- **`environment.yml`** – Set up the conda environment with all dependencies required to run the project.
- **`grphin.py`** – Script for running the GRPhIN algorithm.
- **`enrichment.py`** – Script for calculating graphlet enrichment statistics.
- **`orbit_proteins.py`** – Script for finding protein identities and overrepresented orbits in GRPhIN results.
- **`pageRank.py`** – Script for running simple Random Walk with Restart algorithm to capture a subnetwork based on oxidative stress pathways.
- **`data/`** – Contains raw data files for case studies.
- **`final_output/`** – Contains output data files for case studies.