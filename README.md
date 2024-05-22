Repository for code used in Finlay, Ireland et al *Cancer Cell*, 2024.
[![DOI](https://zenodo.org/badge/770956471.svg)](https://zenodo.org/doi/10.5281/zenodo.10829938)
==========================================================

# Olfactory neuroblastoma mimics molecular heterogeneity and lineage trajectories of small cell lung cancer

## Authors
John B. Finlay<sup>1,7</sup>, Abbie S. Ireland<sup>2,7</sup>, Sarah B. Hawgood<sup>2</sup>, Tony Reyes<sup>2,3</sup>, Tiffany Ko<sup>1</sup>, Rachelle R. Olsen<sup>3</sup>, Ralph Abi Hachem<sup>1</sup>, David Jang<sup>1</sup>, 
Diana Bell<sup>4</sup>, Joseph M. Chan<sup>5</sup>, Bradley J. Goldstein<sup>1,6*</sup>, Trudy G. Oliver<sup>2,3,8*</sup>

<sup>1</sup>1Department of Head and Neck Surgery & Communication Sciences, Duke University, Durham, NC, 27708, USA\
<sup>2</sup>Department of Pharmacology and Cancer Biology, Duke University, Durham, NC, 27708, USA\
<sup>3</sup>Department of Oncological Sciences, University of Utah, Salt Lake City, UT, 84112, USA\
<sup>4</sup>Division of Anatomic Pathology, City of Hope Comprehensive Cancer Center, Duarte, CA 91010, USA\
<sup>5</sup>Human Oncology and Pathogenesis Program, Memorial-Sloan Kettering Cancer Center, New York City, NY, 10065, USA\
<sup>6</sup>Department of Neurobiology, Duke University, Durham, NC, 27708, USA\
<sup>7</sup>These authors contributed equally.\
<sup>8</sup>Lead contact\
<sup>*</sup>Correspondence: Bradley J. Goldstein (bradley.goldstein@duke.edu) and Trudy G. Oliver (tgo@duke.edu)

## Abstract
<p align="center">
<img width="500" alt="image" src="https://github.com/Goldstein-Lab/Finlay_and_Ireland_et_al_mouse_ONB/assets/114078688/9b2f5ecf-aa38-47cf-8b10-20a969e3b0e9">
</p>
The olfactory epithelium undergoes neuronal regeneration from basal stem cells and is susceptible to olfactory neuroblastoma (ONB), a rare tumor of unclear origins. Employing alterations in Rb1/Trp53/Myc (RPM), we establish a genetically-engineered mouse model of high-grade metastatic ONB exhibiting a NEUROD1+ immature neuronal phenotype. We demonstrate that globose basal cells (GBCs) are a permissive cell of origin for ONB, and that ONBs exhibit cell fate heterogeneity that mimics normal GBC developmental trajectories. ASCL1 loss in RPM ONB leads to emergence of non-neuronal histopathologies, including a POU2F3+ microvillar-like state. Similar to small cell lung cancer (SCLC), mouse and human ONB exhibit: mutually exclusive NEUROD1 and POU2F3-like states, an immune-cold tumor microenvironment, intratumoral cell fate heterogeneity comprising neuronal and non-neuronal lineages, and cell fate plasticity—evidenced by barcode-based lineage tracing and single-cell transcriptomics. Collectively, our findings highlight conserved similarities between ONB and neuroendocrine tumors with significant implications for ONB classification and treatment.

## Manuscript
Manuscript is available online [here]

# Data
Data are accessible on GEO, under [Superseries GSE244123]

Raw data are provided as CellRanger Count output matrices. 

Raw fastq files are accessible on the SRA.
  
# Code
Code to replicate analyses performed in [Finlay, Ireland et al. *Cancer Cell* 2024]

1. Download data from NCBI GEO
2. See Jupyter Notebooks and R scripts above for manuscript specific dataset integration methods and analyses.
3. Please refer to [Scanpy](https://scanpy.readthedocs.io/en/stable/), [scvi-tools](https://docs.scvi-tools.org/en/stable/tutorials/index.html), [Single-cell best practices](https://www.sc-best-practices.org/), [Otscomics](https://github.com/cantinilab/OT-scOmics), [Circos](https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_wrapper_circos.md), [TCGAbiolinks](https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html), [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), [GeomxTools](https://bioconductor.org/packages/devel/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html), and [CountCorrect](https://github.com/BayraktarLab/CountCorrect) tutorials for package-specific code that was adapted for datasets in this manuscript. 

# Contact
Please consult methods described in our manuscript for more details or contact corresponding author for specific requests.




