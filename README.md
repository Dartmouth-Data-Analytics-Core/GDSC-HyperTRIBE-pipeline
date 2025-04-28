# GDSC-HyperTRIBE-pipeline
Pipeline for running HyperTRIBE to identify RNA editing sites (In development)

# Table of Contents
- [Introduction](#introduction)
- [Pipeline](#pipeline)
- [Prerequisities](#prerequisites)
- [Implementation](#implementation)
- [Software](#software)
- [References](#references)


## Introduction
HyperTRIBE is a technique used for the identification of the targets of RNA binding proteins (RBP) in vivo. HyperTRIBE couples an RBP to the catalytic domain of the Drosophila RNA editing enzyme ADAR and expresses the fusion protein in vivo. As the RBP-ADARcd (catalytic domain) fusion protein lacks the RNA recognition features of ADAR, the specificity of the RBP should determine the editing specificity of the fusion protein. RBP targets are marked with novel RNA editing events and identified by sequencing RNA. HyperTRIBE identifies RNA editing sites (A to I change, where I is read as G) by comparing RNA sequence from transcriptome with wild type RNA (wtRNA) from the same background. 

## Pipeline
The major steps of this pipeline are as follows:

1. Trim and align sequence libraries to genome
2. Load alignments data to MySQL database
3. Find RNA edit sites wtRNA-RNA approaches

## Prerequisites
This pipeline relies on a series of conda environments and a MySQL Instance through [Dartmouth Dashboard](https://dashboard.dartmouth.edu/login?ticket=ST-53352-NHRbCfD1GUAsNsEkc4n9NgHRZYs-localhost). Should this instance fail, a new instance can be created following [these instructions] (https://services.dartmouth.edu/TDClient/1806/Portal/KB/ArticleDet?ID=150826). Once an instance is generated, credentials can be obtained in the CLI section on Dartmouth Dashboard. 

**Anytime a new instance is generated, the MySQL credentials will need to be edited in the perl source code files in the `CODE` directory**

The current credentials for the MySQL instance are as follows:

|Host|Database|User|Password|
|----|--------|----|--------|
|dmseq-f11b-db.c.dartmouth.edu|dmseq|admin|gdscPass|

## Implementation

## Software
|Tool|Version|Conda|
|----|-------|-----|
|Trimmomatic|0.39|align|
|bowtie2|2.5.4|align|
|

## References
For more details please see:

[HyperTRIBE ReadTheDocs](https://hypertribe.readthedocs.io/en/latest/index.html)

**Xu, W., Rahman, R., Rosbash, M.** Mechanistic Implications of Enhanced Editing by a HyperTRIBE RNA-binding protein. RNA 24, 173-182 (2018). doi:10.1261/rna.064691.117

**McMahon, A.C., Rahman, R., Jin, H., Shen, J.L., Fieldsend, A., Luo, W., Rosbash, M.** TRIBE: Hijacking an RNA-Editing Enzyme to Identify Cell-Specific Targets of RNA-Binding Proteins. Cell 165, 742-753 (2016). doi: 10.1016/j.cell.2016.03.007.

