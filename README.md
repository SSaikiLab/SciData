# Comprehensive data for studying serum exosome microRNA transcriptome in Parkinson’s disease patients

Parkinson’s disease (PD), the second most prevalent neurodegenerative disorder, was classically attributed to alpha-synuclein aggregation and consequent loss of dopaminergic neurons in the substantia nigra pars compacta. Recently, emerging evidence suggested a broader spectrum of contributing factors, including exosome-mediated intercellular communication, which can potentially serve as biomarkers and therapeutic targets. However, there is a remarkable lack of comprehensive studies that connect the serum exosome microRNA (miRNA) transcriptome with demographic, clinical, and neuroimaging data in PD patients. Here, we present serum exosome miRNA transcriptome data generated from four cohort studies. Two of these studies include 96 PD patients and 80 age- and gender-matched controls, with anonymised demographic, clinical, and neuroimaging data provided for PD patients. The other two studies involve 96 PD patients who were evaluated both before and after one year of treatment with rasagiline, a widely prescribed anti-parkinsonism drug. Together, the datasets provide a valuable source for understanding pathogenesis and discovering biomarkers and therapeutic targets in PD.

# 20210
# 2021
The folders contain the serum exosome miRNA transcriptomics data generated from PD participants and healthy control from our 2020 and 2021 cohort.

# 2022
# 2023
The folders contain the serum exosome miRNA transcriptomics data generated from PD participants before and after 1 year of rasagiline treatment.

The R scripts perform DESeq2 based analysis for identifying differentilly expressed miRNAs, with significantly (FDR < 0.05) increased miRNAs and decreased miRNAs provided as individual excel files.
The excel files named results or GSEA denote the overall difference in miRNA expression, with miRNAs ranked by |log_2(⁡FDR)|*FC. 
