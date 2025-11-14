Project description: 
- Protein Expression Data Analysis in Mouse Models through applying machine learning and statistical techniques to analyze the Mice Protein Expression dataset (72 samples, 8 experimental groups) and identifying discriminant protein expression patterns across genotype, treatment, and behavioral conditions.

Steps:
- Data preprocessing: Cleaned missing data, removed low-quality protein features, and ensured balanced class distributions across biological groups.
- Dimensionality reduction (PCA): Reduced 71 protein features to 9 principal components explaining 95% of total variance, enabling efficient downstream modeling.
- Decision Tree Classifier: Built interpretable rule-based models to classify mouse groups by protein profiles, exploring feature importance and model complexity.
- k-Nearest Neighbors (kNN): Implemented instance-based learning with cross-validation and bootstrap sampling; achieved high classification accuracy (AUC ≈ 99.4%).
- Artificial Neural Networks (ANN): Trained multi-layer perceptron models to capture nonlinear expression patterns; achieved best performance (AUC ≈ 99.6%).
- Clustering (K-means, Biclustering): Applied unsupervised learning to identify homogeneous protein expression clusters and co-regulated protein subsets across experimental conditions.
- ANOVA & statistical inference: Performed multi-factor ANOVA to detect proteins with significant differential expression, visualized results with Manhattan-style plots, and linked findings to behavioral effects.
