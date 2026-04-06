# High-Efficiency Prediction of Surfactant Properties through Automated Head–Tail Partitioning  

This repository contains the code, trained models, and results associated with the paper, ensuring full reproducibility of the proposed methodology.  

Each property prediction model is provided as a pair of files:  
- `.joblib`: serialized trained model  
- `joblib_parameters`: list of required molecular descriptors and metadata for inference  

---

## Overview  

The core concept of this work is the automated decomposition of surfactant molecules into chemically meaningful **hydrophilic head** and **hydrophobic tail** domains. Descriptors are computed at the fragment level and used in machine learning models to predict key physicochemical and performance-related properties, including:  

- Critical micelle concentration (CMC)  
- Surface tension at CMC  
- Adsorption efficiency (pC20)  
- Interfacial adsorption parameters:  
  - Maximum surface excess concentration (τₘₐₓ)  
  - Langmuir adsorption constant  

Surfactant performance is strongly governed by amphiphilic structure. Conventional whole-molecule descriptors often fail to capture this asymmetry, particularly when comparing molecules with similar global composition but different head–tail architectures.  

This repository addresses this limitation by combining automated molecular partitioning with domain-specific descriptors, enabling more physically meaningful and robust modeling.  

---

## Key Features  

- Automated partitioning of surfactants into head and tail domains from SMILES  
- Computation of fragment-aware molecular descriptors  
- Pretrained models with clearly defined descriptor requirements  
- Reproducible workflows for property prediction  

---

## Quick Start: Partition a single molecule  

The main script accepts either direct SMILES strings A) or an input file B).
### A) Partition one molecule from the command line

```bash
python molecular_partitioning_script.py --smiles "CCCCCCCCCCCC(=O)NCCC[N+](C)(C)C(O)C[S](O)(=O)=O"
```

### B) Partition molecules from a spreadsheet

```bash
python molecular_partitioning_script.py \
	--input your_file.xlsx \
	--sheet Sheet1 \
	--smiles-col smiles \
	--output partitioned_surfactants.xlsx
```

### C) Save colored structure images

	--output partitioned_surfactants.xlsx \
	--draw-dir partition_images


### Main outputs

Typical output columns include:

- `Surfactant`, `Head`, and `Tail` SMILE notations along with the 24 features proposed in this work. 


### How the partitioning works

At a high level, the workflow:

1. identifies hydrophobic carbon candidates
2. groups them into connected components
3. selects the dominant carbon-rich component(s) as tail region(s)
4. grows the head region from hetero atoms and charged atoms
5. optionally refines the boundary using local charge information
6. exports head/tail fragment SMILES and summary features

This is a heuristic workflow, designed for high-throughput screening rather than exact mechanistic decomposition.

---

## Using trained models

Model artifacts are generally stored as pairs:

- `.joblib`: trained model object
- `joblib_parameters`: list of descriptors and metadata required at inference time

This pattern is used throughout the modelling workflows and is one of the main reproducibility mechanisms in the repository.

## Datasets and Results  

In addition to the code and trained models, this repository provides all results generated for the datasets studied in this work.  

The calculated molecular features (descriptors) are stored in `.xlsx` files, while the visual representations of the head–tail partitioning for each molecule and clustering are provided as compressed `.zip` files.  All data, descriptors, models, and clustering results are organized into three main groups according to the source dataset and target property: **Qin**, **Seddon**, and **pC20**.  

- **Qin dataset**  
  Used exclusively for **CMC prediction**. It contains experimental CMC values for a diverse set of surfactants and is employed to train and evaluate machine learning models and clustering analyses focused on micellization behavior.  

- **Seddon dataset**  
  Used for multiple interfacial properties, including:  
  - CMC  
  - Surface tension  
  - Langmuir constant (Kₗ)  
  - Maximum surface excess concentration (τₘₐₓ)   

- **pC20 dataset (Li et al.)**  
  Used for **adsorption efficiency (pC20)** prediction.

Therefore, the properties evaluated are:

### Log(KL) – Langmuir Constant  

The dataset reported by Seddon et al. includes 154 hydrocarbon surfactants categorized into:  

- 44 ethoxylates , 31 sulphates , 20 alcohols , 14 carboxylates, 14 amides, 6 betaines, 6 sulphonates, 5 tetra-alkyl ammonium surfactants, 4 glucosides, 4 pyrrolidinones, 3 pyridinium-based surfactants, 3 glyceryl-based surfactants  

A PLS regression model with 13 principal components was used.  

---

### Maximum Surface Excess Concentration (τₘₐₓ)  

Values of τₘₐₓ for the same set of 154 surfactants were also reported by Seddon et al.  

An MLPRegressor model was employed, consisting of two hidden layers with 8 and 3 neurons, respectively. The hyperbolic tangent (tanh) activation function was used, with a learning rate of 0.1.  

---
### Log(CMC)  

Two datasets are used to build models for predicting critical micelle concentration (CMC):  

- Qin et al. compiled experimental CMC data for 202 surfactants, including 122 nonionic, 35 cationic, 34 anionic, and 11 zwitterionic surfactants. An eXtreme Gradient Boosting (XGBoost) model was trained using a maximum tree depth of 5 and a gamma value of 0.1.  

- Seddon et al. reported CMC data for 91 surfactants. A Partial Least Squares (PLS) regression model with 12 principal components (PCs) was applied.  

### Surface Tension  

Surface tension data for the same set of 91 surfactants reported by Seddon et al. are also included. PLS regression with 16 principal components was employed. Notably, Seddon et al. were among the first to explicitly correlate interfacial properties with molecular descriptors, making this dataset particularly valuable.  

---

### Adsorption Efficiency (pC20)  

Li et al. reported a dataset comprising 124 surfactants, including quaternary ammonium Gemini surfactants, anionic Gemini surfactants with amide groups, and long-chain polyethoxylated surfactants.  

Following their approach, an MLPRegressor model was implemented with two hidden layers of 50 neurons each. The relu activation function and the lbfgs solver were used, corresponding to the optimal hyperparameter combination reported in the original study.  

---
### Additional Dataset: HLB  

In addition to the main datasets, molecular features and head–tail partitioning results are also provided for a dataset of **120 nonionic surfactants** used for hydrophile–lipophile balance (HLB) analysis.  

This dataset includes ethoxylated alcohols, ethoxylated fatty acid esters, ethoxylated dialkyl acids, and ethoxylated alkylphenols. HLB values were calculated following the group contribution method reported by Guo et al.  

Although no predictive models are developed for this dataset, it is included to demonstrate the applicability of the domain-specific descriptor framework and head–tail partitioning methodology to classical surfactant design parameters such as HLB.  

**References:**  
Qin, S.; Jin, T.; Van Lehn, R. C.; Zavala, V. M. *Predicting Critical Micelle Concentrations for Surfactants Using Graph Convolutional Neural Networks*. J. Phys. Chem. B 2021, 125 (37), 10610–10620.

Seddon, D.; Müller, E. A.; Cabral, J. T. *Machine Learning Hybrid Approach for the Prediction of Surface Tension Profiles of Hydrocarbon Surfactants in Aqueous Solution*. J. Colloid Interface Sci. 2022, 625, 328–339.

Li, S.; Mao, X.; Cao, X.; Feng, Y.; Zhang, Y.; Yin, H. *Visualizing Molecular Structure—Adsorption Efficiency Relationship through an Interpretable Machine Learning Strategy*. Colloids Surf. A Physicochem. Eng. Asp. 2025, 712, 136400.
Guo, X.; Rong, Z.; Ying, X. *Calculation of Hydrophile–Lipophile Balance for Polyethoxylated Surfactants by Group Contribution Method*. J. Colloid Interface Sci. 2006, 298 (1), 441–450.  
