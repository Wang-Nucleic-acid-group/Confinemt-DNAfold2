#DNAfold2 : A Coarse-Grained Model for Ab Initio Prediction of DNA Structures under Spatial Confinement


Introduction
================================================================
This package employs a coarse-grained (CG) model to predict the 3D structures and stabilities of complex DNAs in ion solutions under varying spatial confinement. The model is an extension of an RNA CG model originally developed by Tan's group at Wuhan University, which was designed for ab initio prediction of RNA 3D structure, stability, and salt effects.


Model Features:
================================================================

- (1) Three-bead coarse-grained model including sequence-dependent base-pairing/stacking potentials and structure-based electrostatic potential
- (2) Predicts the 3D structure of DNAs with multi-way junctions from their sequences
- (3) Reproduces thermal stability of DNA junctions across diverse sequences and lengths
- (4) Supports simulation of DNA folding under user-defined spatial confinement sizes




Requirements:
================================================================
- (1) Operating System: Linux or MacOS
- (2) CPU: ≥ 8 cores with 2 threads
- (3) C++ Compiler: GCC ≥ 7.5 (with C++11 support)
- (4) Python: ≥ 3.11.5 and Biopython ≥ 1.81



Usage:
================================================================
1. Prepare Input Files
- Put the target sequence in seq.dat (e.g., ATCCTAGTTATAGGAT)
- Modify the Configuration File (config.dat), including:
  - (1) Sampling_type:	Specifies the sampling algorithm: 1 for Replica Exchange Monte Carlo (REMC), 0 for Monte Carlo Simulated Annealing (MCSA).
  - (2) Folding_steps:	Number of steps for the REMC or MCSA simulation during structure folding.
  - (3) Optimizing_steps: Number of steps for energy minimization and structural optimization after the sampling phase.
  - (4) C_Na:	Concentration of Na⁺ ions in mM.
  - (5) C_Mg:	Concentration of Mg²⁺ ions in mM.
  - (6) Confine_size:	(Optional) Diameter of the spatial confinement sphere in nm.


2. Run Simulation
- bash run.sh



3. Output Files
- Results are stored in the results/directory. Key subdirectories include:
  - (1) Folding_trajectory/: Trajectories at different temperatures
  - (2) CG_structures/: Predicted coarse-grained 3D structures
  - (3) Secondary_structure/: Top-N secondary structures in dot-bracket notation
  - (4) All_atom_structure/: All-atom models converted from CG predictions
  - (5) Thermal_Stability/: Thermal stability analysis files


4. Thermal_Stability/contents:
- thermal_stability.dat: Population fractions of Folded/Unfolded/Intermediate states across temperatures
- Secondary_sec_x.dat: Dot-bracket structures at different temperatures



Recommended Settings:
================================================================

1. For 3D Structure Prediction:
- Folding steps: ≥ 750,000 (minimum 500,000)
- Optimizing steps: ≥ 500,000 (minimum 100,000)


2. For Thermal Stability Prediction:
- Folding steps: ≥ 4,000,000 (minimum 2,000,000)
- To compute melting temperature (Tm), fit the folded/unfolded fractions to a two-state model as in references [1–2]

3. For Spatial Confinement Simulations:
- Specify the confinement diameter (in nm) in config.dat to simulate folding under restricted volumes



References:
================================================================
[1] Z.-C. Mu, Y.-L. Tan, B.-G. Zhang, J. Liu, Y.-Z. Shi. Ab initio predictions for 3D structure and stability of single- and double-stranded DNAs in ion solutions. PLoS Comput. Biol. 18, e1010501 (2022).

[2] X. Wang, Y.-Z. Shi. 3D structure and stability prediction of DNA with multi-way junctions in ionic solutions. PLoS Comput. Biol. 8, e1013346 (2025).


