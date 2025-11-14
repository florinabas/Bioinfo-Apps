Project description: 
- Implementing and optimizing end-to-end pipelines for detecting and annotating genomic variants from NGS data in Saccharomyces cerevisiae. Gained hands-on experience applying supervised and unsupervised machine learning methods, dimensionality reduction, and statistical testing to biological datasets.

Steps:
- Galaxy workflow: Built a complete variant calling pipeline integrating tools such as FastQC, BWA, Samtools, Picard, FreeBayes, SnpEff, bcftools, and MultiQC. The workflow performed quality control, read alignment to the sacCer3 reference genome, variant detection, annotation, and comprehensive reporting for multiple biological replicates.
- Bcbio-nextgen pipeline: Configured a GATK-based germline variant calling workflow using Bowtie2 for alignment and HaplotypeCaller for variant detection, including automated SnpEff annotation and multi-sample VCF comparison with bcftools.
- Variant annotation and comparative analysis in R: Developed R scripts leveraging Bioconductor packages (VariantAnnotation, GenomicRanges, BSgenome, TxDb, ReactomePA, topGO) to annotate variants, identify coding and intergenic changes, determine nonsynonymous SNPs, and perform functional enrichment analysis. Used bedr for intersection and Venn analyses to compare results across replicates and between pipelines.
