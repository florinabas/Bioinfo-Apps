**Machine Learning & Statistical Inference**
Project description: Protein Expression Data Analysis in Mouse Models through applying machine learning and statistical techniques to analyze the Mice Protein Expression dataset (72 samples, 8 experimental groups) and identifying discriminant protein expression patterns across genotype, treatment, and behavioral conditions.

Steps:
- Data preprocessing: Cleaned missing data, removed low-quality protein features, and ensured balanced class distributions across biological groups.
- Dimensionality reduction (PCA): Reduced 71 protein features to 9 principal components explaining 95% of total variance, enabling efficient downstream modeling.
- Decision Tree Classifier: Built interpretable rule-based models to classify mouse groups by protein profiles, exploring feature importance and model complexity.
- k-Nearest Neighbors (kNN): Implemented instance-based learning with cross-validation and bootstrap sampling; achieved high classification accuracy (AUC ≈ 99.4%).
- Artificial Neural Networks (ANN): Trained multi-layer perceptron models to capture nonlinear expression patterns; achieved best performance (AUC ≈ 99.6%).
- Clustering (K-means, Biclustering): Applied unsupervised learning to identify homogeneous protein expression clusters and co-regulated protein subsets across experimental conditions.
- ANOVA & statistical inference: Performed multi-factor ANOVA to detect proteins with significant differential expression, visualized results with Manhattan-style plots, and linked findings to behavioral effects.

**NGS Variant Analysis and Annotation Pipelines (Galaxy, Bcbio-nextgen, R/Bioconductor)**
Project description: Implementing and optimizing end-to-end pipelines for detecting and annotating genomic variants from NGS data in Saccharomyces cerevisiae. Gained hands-on experience applying supervised and unsupervised machine learning methods, dimensionality reduction, and statistical testing to biological datasets.

Steps:
1. Galaxy workflow: Built a complete variant calling pipeline integrating tools such as FastQC, BWA, Samtools, Picard, FreeBayes, SnpEff, bcftools, and MultiQC. The workflow performed quality control, read alignment to the sacCer3 reference genome, variant detection, annotation, and comprehensive reporting for multiple biological replicates.
2. Bcbio-nextgen pipeline: Configured a GATK-based germline variant calling workflow using Bowtie2 for alignment and HaplotypeCaller for variant detection, including automated SnpEff annotation and multi-sample VCF comparison with bcftools.
3. Variant annotation and comparative analysis in R: Developed R scripts leveraging Bioconductor packages (VariantAnnotation, GenomicRanges, BSgenome, TxDb, ReactomePA, topGO) to annotate variants, identify coding and intergenic changes, determine nonsynonymous SNPs, and perform functional enrichment analysis. Used bedr for intersection and Venn analyses to compare results across replicates and between pipelines.

**Biostatistics in R**
Project description: Analyzed clinical data from a study of 29 patients receiving Botox treatment for migraines, focusing on treatment efficacy and factors influencing migraine frequency, intensity and presence of aura. Gained hands-on experience with R for data handling, visualization and statistical analysis, applied non-parametric tests and regression models to real clinical data and learned to interpret results in a biomedical context and communicate findings clearly.

Steps:
- Data Cleaning & Preparation: Standardized variable formats, handled missing values and outliers.
- Exploratory Analysis: Visualized distributions, compared groups by age, gender, and aura.
- Statistical Testing: Applied Wilcoxon tests, ANOVA, and Levene’s test to assess differences in migraine outcomes.
- Regression Modeling: Built linear and logistic regression models to quantify the impact of treatment and patient factors on migraine intensity and aura presence.
- Prediction & Interpretation: Predicted outcomes for missing data points and evaluated model performance.

**Stochastic Modeling and Bayesian Simulation in Computational Biology**
Project description: Developed and applied stochastic and probabilistic models to study biological sequences and genetic data. Combined Markov chains, Hidden Markov Models (HMMs), Monte Carlo simulations, and Markov Chain Monte Carlo (MCMC) techniques to simulate processes, predict hidden states, and estimate biological parameters such as nucleotide distributions and ABO blood group allele frequencies.

Steps:
- Markov and HMM Modeling: Defined model properties, built transition and emission matrices, simulated DNA sequence evolution, and applied the Viterbi algorithm to identify coding regions.
- Random Walk Simulations: Explored 1D and 2D stochastic processes to understand sequence dynamics and probabilistic behavior.
- Monte Carlo Simulations: Generated random samples using inverse CDF and direct sampling to simulate target probability distributions.
- Gibbs Sampling & MCMC: Implemented Gibbs sampling for multivariate distributions, analyzed effects of correlations and initial values, and estimated posterior distributions of allele frequencies in ABO blood groups using Bayesian methods.
- Bayesian Analysis & Inference: Defined priors and likelihoods, ran simulations in R and WinBUGS, and interpreted marginal and stationary distributions for reliable biological inference.

**Systems Biology and Multiscale Modeling**
Project description: building and simulating biochemical and cellular models using SBML, COPASI and Morpheus.
Steps:
- Applied mathematical modeling, differential equations, and agent-based simulations to study intra-cellular processes and multicellular systems.
- Translating complex biological processes into computational models for predictive analysis.

**Algorithms for Biological Sequence Analysis**

Project description: Implemented core bioinformatics algorithms for DNA/RNA processing, pattern detection, and basic genetic simulations, using problems from the Rosalind platform. Strengthened practical skills in sequence manipulation and algorithmic problem-solving relevant to genomic data analysis.

Steps:
- Wrote functions for nucleotide counting, DNA→RNA transcription, strand complement, and motif search.
- Calculated Hamming distance, translated RNA to protein, and built consensus/profile matrices.
- Simulated simple biological models such as the rabbit recurrence problem.
- Implemented the Needleman–Wunsch global alignment algorithm to compare DNA/protein sequences using customizable match, mismatch, and gap penalties.

**Structural Bioinformatics**

Project description: Conducted a structural bioinformatics study of Cyclooxygenase-1 (COX-1), the main predicted molecular target of Ibuprofen. The project involved identifying the target using SwissTargetPrediction, retrieving functional and structural data from UniProt and the Protein Data Bank (PDB), and performing detailed 3D structural analysis with UCSF Chimera. Additionally, protein surface pockets and binding cavities were characterized using CASTp to better understand potential interaction sites.

Steps:
- Predicted Ibuprofen molecular targets using SwissTargetPrediction based on SMILES similarity.
- Selected COX-1 as the primary target and retrieved its sequence and functional information from UniProt (P23219).
- Downloaded the crystallographic structure PDB: 6Y3C from the Protein Data Bank.
- Visualized the structure in UCSF Chimera, examining: backbone and full-atom representations, secondary structure elements, ligand positioning, hydrophobic surface and electrostatic potential, protein surface and volume calculations and structural comparison using RMSD.
- Identified and measured binding pockets and surface cavities using CASTp, including the main cavity.

**Chemoinformatics**
Project description: Performed a chemoinformatics analysis of key organic compounds in Coffea robusta. The project included identifying major constituents (e.g., caffeine, chlorogenic acids, fatty acids), characterizing their physicochemical properties, and predicting ADME/Tox profiles using multiple in-silico tools. A molecular docking study was also carried out to evaluate interactions between octadecanoic acid and a human CYP enzyme.

Steps:
- Collected major coffee compounds from FooDB and retrieved physicochemical data from PubChem.
- Predicted ADME and drug-likeness using SwissADME (GI absorption, BBB permeability, CYP interactions).
- Evaluated toxicity and safety using admetSAR, Pred-hERG, Pred-Skin, and CarcinoPred-EL.
- Performed ligand–protein docking of octadecanoic acid with a human cytochrome using SwissDock, AutoDock Vina, PLIP, and PatchDock/FireDock.

