# OncologyTT-Datasets
### Oncology Therapeutic Targets (OncologyTT) datasets
#### a collection of drugs and target genes associated with 13 cancer types.

![alt text](https://github.com/MahaThafar/OncologyTT-Datasets/blob/main/OncologyTT_logo.png)
---

> OncologyTT Datasets have been created for research purposes to contribute to Oncology research, more specifically "Identifying novel therapeutic targets to treat cancer".


#### This repository includes several datasets related to 13 cancer types.
#### Description of all folders and files as follow:

> All data are divided into files based on cancer type.
> Data Curation and preprocessing steps to obtain some data is explained in each folder,
> and the preprocessing implementations are provided (for some datasets)

----
**It consists of Three main files, which are:**
1. all_positive_targets_wih_info.xlsx,  that are considered as "Positive Data Samples" including all genes that are "approved" target interacted with drugs for the 13 cancer types combined together in one file with other information.
2. all_positive_genes_with_info.xlsx, also considered "Positive Data Samples," including all cancer driver over-regulated genes for the 13 cancer types combined together in one file with other information.
3. all_negative_genes_with_info.xlsx that are considered as "Negative Data Samples" including all human genes randomly generated from a pool of Human Genes after excluding all positive genes for the 13 cancer types combined together in one file with other information.

**Each file has this information:**
- *[Entry, Gene names, Protein names, Amino-acid sequence, Class Label]*


**Each one of the previous files has a corresponding folder:**
1. Positive targets (Approved target)
2. Positive genes (Cancer driver genes) 
3. Negative genes (unlabeled genes thaat have been used as negative samples)\ 
that consist of several files matching the cancer type.

**Other Data Folders:**/
4. Omics features: that include gene expression and gene mutation features for all genes per cancer across multiple tumor samples./
5. Approved drugs per cancer: that are collected from National Cancer Institue (NCI).

-----------------------------------------------
#### If you use any data of OncologyTT datasets, please cite the following:
TBD

--------------------------------------------------------------------
### For any qutions please contact the first author:
---------------------------------------------------------------------
Maha A. Thafar \
Ph.D. Candidate and resercher | Computer Science\
Computer, Electrical and Mathematical Sciences and Engineering (CEMSE) Division\
Computational Bioscience Research Center (CBRC), King Abdullah university of science and technology.\
Collage of Computers and Information Technology, Taif University (TU).\
Email: maha.thafar@kaust.edu.sa

----
