# Fitting Potts model to an MSA

Code for fitting the Potts model as described in the paper. Code is adapted from [Ekeberge et. al](https://doi.org/10.1016/j.jcp.2014.07.024).

Instructions for downloading and running the code as described by the original authors can be found in the linked [GitHub repository](https://github.com/mskwark/plmDCA_asymmetric_v3). 


## Fitting Potts model 

Fitting the Potts model to get single-site (h) and coupling (J) coefficients requires the following steps:

1. Running the `plmDCA_asymmetric_MS_edit.m` script in MATLAB to generate the intermediate h- and J- coefficient matrices and matrices.
2. Running the `plm_gauge_transform_save.py` script to process, transform, and save the h- and J-coefficient matrices to the zero-sum gauge. The resulting h- and J-matrices from this script are used for all downstream analysis and sequence design.

### Instructions for plmDCA_asymmetric_MS_edit.m
Used to fit Potts model and generate intermediate h- and J-coefficient matrices. 

A sequence weighting threshold of 80% and a regularization strength of 0.01 are hard-coded into the script. 

#### Example run command
```
plmDCA_asymmetric_MS_edit(fastafile, outputDIfile, nr_of_cores)
```

#### Input parameters 
 * `fastafile`: FASTA file of MSA for fitting Potts model. Examples can be found in the [MSAs directory](https://github.com/barricklab-at-jhu/Potts_analysis_Sternke_Tripp_Barrick/tree/main/MSAs) of this repository.
 * `outputDIfile`: filename for output file of Iirect Information as described in the Ekeberg et al study. This file is not used for anything in this current study. 
 * `nr_of_cores`: Number of cores to use for parallel processing

#### Outputs
 * `h.mat`: Matrix of intermediate h-coefficients. Will be further processed to transform to zero-sum gauge by the `plm_gauge_transform_save.py` script. 
 * `J.mat`: Intermedaite matrix of J coefficients in zero-sum gague. Will be saved as numpy array by `plm_gauge_transform_save.py` script. 
 * `Jtemp1.mat` and `Jtemp2.mat`: Intermediate matrices of asymmetric J-coefficients. Used by the `plm_gauge_transform_save.py` script for generating zero-sum gauge transformed h-coefficient matrix.
 * `msa_seq_vectors.mat`: MSA of vector encoded sequences. Order of sequences will be altered from original MSA file.
 * `seq_weights.mat`: Vector of calculated sequence weights. Weights will be in same order as seqs in msa_seq_vectors.mat

### Instructions for plm_gauge_transform_save.py
After fitting plmDCA_asymmetric_only_h_MS_edit.m script, run plm_gauge_transform_save.py script to transform h-coefficents into zero sum gauge, and to save the final h- and J-coefficient matrices as numpy arrays. 
#### Example run command
```
python plm_gauge_transform_save.py
```
#### Input parameters
The `h.mat`, `J.mat`, `Jtemp1.mat` and `Jtemp2.mat` are coded in script. 

#### Output 
 * `h_plm_zero_sum.npy`: Final matrix of h-coefficients. Matrix dimensions are ( L x 21 ) where L is the length of sequences in the MSA and 21 is the amino acid alphabest used (with order -ACDEFGHIKLMNPQRSTUVWY).
 * `j_plm_zero_sum.npy`: Final matrix of J-coefficients. Matrix dimensions are (L x L x 21 x 21). 

These resulting h- and J-matrices  are used for all downstream analysis and sequence design.

## Fitting independent model 
The above script was further adapted to fit the site independent model by removing the coupling coefficients from model fitting. 

To fit the independent model run the `plmDCA_asymmetric_only_h_MS_edit` script in MATLAB by: 

```
plmDCA_asymmetric_ind_MS_edit(fastafile, outputDIfile, nr_of_cores)
```

As with the Potts model script a sequence weighting threshold of 80% and a regularization strength of 0.01 are hard-coded into the script. 

### Inputs
Inputs are the same as described above for the `plmDCA_asymmetric_MS_edit.m` scripts.

### Outputs
 * `h_onlyH.mat`: Matrix of independent model single-site coefficents (I-coefficients). Matrix dimensions are ( L x 21 ). 

This matrix is the final coefficients used for all downstream analysis and sequence design.

---------------------------------------------
 Copyright 2014 - by Magnus Ekeberg (magnus.ekeberg@gmail.com)
 All rights reserved
 
 Permission is granted for anyone to copy, use, or modify this
 software for any uncommercial purposes, provided this copyright 
 notice is retained, and note is made of any changes that have 
 been made. This software is distributed without any warranty, 
 express or implied. In no event shall the author or contributors be 
 liable for any damage arising out of the use of this software.
 
 The publication of research using this software, modified or not, must include 
 appropriate citations to:

 	M. Ekeberg, C. Lövkvist, Y. Lan, M. Weigt, E. Aurell, Improved contact
 	prediction in proteins: Using pseudolikelihoods to infer Potts models, Phys. Rev. E 87, 012707 (2013)

	M. Ekeberg, T. Hartonen, E. Aurell, Fast pseudolikelihood
	maximization for direct-coupling analysis of protein structure
	from many homologous amino-acid sequences, J. Comput. Phys. 276, 341-356 (2014)
