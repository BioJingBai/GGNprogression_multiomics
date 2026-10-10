# Genomic Instability and Coordinated Tumor Core-Margin Immune Remodeling Define Divergent Progression Trajectories in Pulmonary Ground-Glass Nodules


This repository contains the analysis pipeline used to define the longitudinal growth patterns of pulmonary ground-glass nodules (GGNs) from serial CT and to characterise the tumour and tumour-margin microenvironment that distinguishes them.

The pipeline includes:

- Longitudinal CT kinetics and growth-pattern classification (`01_clinical_kinetics/`)
- Bulk transcriptome signature, deconvolution and enrichment analysis (`02_bulk_transcriptome/`)
- Visium HD spatial transcriptomics (`03_visium_hd/`)
- CosMx SMI spatial molecular imaging (`04_cosmx_smi/`)
- Imaging mass cytometry (`05_imc/`)
- Multiplex immunofluorescence image analysis (`06_mif/`)
- TCR repertoire clonality (`07_tcr_clonality/`)

```bash
export PROJECT_ROOT=/path/to/GGN_progression_data
```

```bash
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8
```

The patient data — CT scans, sequencing data, and spatial and imaging datasets — are controlled access and are not distributed with this repository; access is governed by the data availability statement of the associated manuscript. 

Released under the MIT License — see `LICENSE`. The license covers the source code only and does not extend to any patient-derived data.

R 4.4.1 · Python 3.10.12
