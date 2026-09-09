# OngoingCIN_scWGS

This repo holds the code and data needed to reproduce the analysis for *Tracking ongoing chromosomal instability using single-cell whole-genome sequencing*

## Data

The data needed for reproducing analyses and figures are in deposited in [figshare](https://figshare.com/s/40550ebf156588728eea). Data should be added in the */data* folder.

## Helper scripts

The functions needed for running the code are in the [/scripts](https://github.com/bhernando/OngoingCIN_scWGS/tree/main/scripts) folder

## Linking code to manuscript figures (main, extended and supplement)

The code is structured in 7 directories which we recommend to run in the following order to analyze the data and reproduce figures in the manuscript:

-   [00_Adapt_FeatureEncoding](https://github.com/bhernando/OngoingCIN_scWGS/tree/main/00_Adapt_FeatureEncoding): This repository contains code to justify the design of the new feature and the modification of the changepoint feature (Supplementary Note 1 and Methods)

-   [01_Derive_Gold_CINSignatures](https://github.com/bhernando/OngoingCIN_scWGS/tree/main/01_Derive_New_CINSignatures): This repository contains code to derive signatures that are able to be mapped to individual events (Supplementary Note 1 and Methods).

-   [02_Derive_New_scCINSignatures](https://github.com/bhernando/OngoingCIN_scWGS/tree/main/02_Derive_New_scCINSignatures): This repository contains code to adapt signatures to the single-cell scenario, where only integer copy number values are present (no subclonal alterations) and sequencing depth allows for a minimal bin resolution of 100kb (Supplementary Note 1 and Methods).

-   [03_InSilico_Validation](https://github.com/bhernando/OngoingCIN_scWGS/tree/main/03_InSilico_Validation): This repository contains code to validate the new 3-feature encoding to identify CIN type active in a sample (Supplementary Note 1) and the event-to-signature mapping procedure (Figure 1e, Supplementary Figure 5).

-   [04_Process_scWGS](https://github.com/bhernando/OngoingCIN_scWGS/tree/main/04_Process_scWGS): This repository contains the code used to process single-cell profiles derived from whole-genome sequencing (WGS) data. For each sample, only profiles meeting our inclusion criteria were retained, while noisy or over-fragmented profiles were excluded. The same processing pipeline was applied to all samples, and example data are provided to enable analysis of external datasets. Raw data from our experimentally-derived models will be made available upon publication. This repository also includes the code used to derive the null distribution for identifying cell-unique events.

-   [05_Induced_CIN](https://github.com/bhernando/OngoingCIN_scWGS/tree/main/05_Induced_CIN): This repository contains the code for validating our framework to identify ongoing CIN processes using in vitro models with a specific CIN type induced (Figure 2).

-   [06_Ongoing_CIN](https://github.com/bhernando/OngoingCIN_scWGS/tree/main/06_Ongoing_CIN): This repository contains the code for using ongoing CIN signature for predicting drug response (Figure 3) and for dissecting subclonal diversification in TNBC tumours (Figure 3-4).

-   [07_Miscellaneous](https://github.com/bhernando/OngoingCIN_scWGS/tree/main/07_Miscellaneous): This repository contains the code for performing additional analyses

    -   The *01_ncells_unique_fp_detection.R,* *01a_minussi_real_unique.R* and *01b_unique_fp_detection.R* scripts are used for computing the minimum number of cells required to reliably distinguish cell-unique copy number alterations from clonal/shared events due to potential undersampling (Figure 1c).

    -   The *02_cnas_drivers_subMinussi.R* script is used for exploring amplifications and deletions in driver genes in TNBC subclones.

    -   The *03_hmf_vs_tcga_CX10.R* script is used for comparing CX10 activities between primary and metastatic TNBC tumours (Figure 4d).

    -   The *04_paired_tnbc_msk_CX10.R* script is used for comparing CX10 activities between primary-metastases pairs from TNBC tumours (Figure 4e).

    -   The *05_unique_stats.R* script is used for assessing our statistical approach to identify true and false positive cell-unique events (Supplementary Figure 3).

    -   The *06a_nrpbpcc_threshold.R* and *06b_nrpbpcc_threshold.R* scripts are used for assessing the impact of sequencing depth on excluding underpowered cells and minimize sequencing noise (Supplementary Figure 5-7).

    -   The *07_qc_threshold.R* script is used for assessing the impact of the quality controls applied to the different steps of our single-cell framework (Supplementary Figure 8).
