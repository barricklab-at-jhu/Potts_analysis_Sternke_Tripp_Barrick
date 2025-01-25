Fitting asymmetric plmDCA to an MSA

To fit Potts model
Run plmDCA_asymmetric_MS_edit file in Matlab
Example call for hd file in directory: plmDCA_asymmetric_MS_edit('hd_full_gapStrip_nonredundant.txt, DI_output.txt, 2)

    I have modified the original script to hard code in some values
        sequence weighting threshold set at 80%
        regularization coefficient set at 0.01

    Inputs for file
        fastafile - FASTA file for MSA
        output file - filename for output file of direct information
            not used for anything if just using h,J coefficients
        nr_of_cores - number of cores to run parallel processing with

    Program will output a number of matrices/vectors
        msa_seq_vectors.mat - MSA with vector encoded sequences
            note: the script will re-order sequences from original MSA file
        seq_weights.mat - vector of calculated sequence weights
            note: weights will be in same order as seqs in msa_seq_vectors
        h.mat - matrix of h coefficients 
            note: not in zero-sum gauge
            do not use these as final h coefficients
            use 'plm_gauge_transform_save.py' script to transform to zero-sum gauge
            use output of 'plm_gauge_transform_save.py' script as final h coefficients
        Jtemp1.mat and Jtemp2.mat
            note: not in zero-sum gauge
            J coefficients determined from Jij and Jji
            do not use as final J coefficients
            used for h coefficient gauge transformation
        J.mat
            matrix of J coefficients in zero-sum gague
            use these as final J coefficients
            use 'plm_gauge_transform_save.py' to save as numpy array for use in Python

After fitting plmDCA_asymmetric_only_h_MS_edit.m script, run plm_gauge_transform_save.py script to transform h coefficents into zero sum gauge, and to save the h and J coefficient matrices as numpy arrays (h_plm_zero_sum.npy and j_plm_zero_sum.npy). 
	If protein sequences are length L and using 21 letter amino acid alphabet
	Order of residues: -ACDEFGHIKLMNPQRSTUVWY

            h coeff array is L x 21 matrix
            j coeff array is an L x L x 21 x 21 matrix



To fit independent model (removing J from model, "consensus" model)
Run plmDCA_asymmetric_only_h_MS_edit file in Matlab
Example call for hd file in directory: plmDCA_asymmetric_only_h_MS_edit('hd_full_gapStrip_nonredundant.txt, DI_output.txt, 2)

    I have modified the original script to remove J from model and hard code in some values
        sequence weighting threshold set at 80%
        regularization coefficient set at 0.01

    Same inputs as full model

    Outputs matrix of coefficients (called h_onlyH.mat)




 The publication of research using this software, modified or not, must include 
 appropriate citations to:

 	M. Ekeberg, C. Lövkvist, Y. Lan, M. Weigt, E. Aurell, Improved contact
 	prediction in proteins: Using pseudolikelihoods to infer Potts models, Phys. Rev. E 87, 012707 (2013)

	M. Ekeberg, T. Hartonen, E. Aurell, Fast pseudolikelihood
	maximization for direct-coupling analysis of protein structure
	from many homologous amino-acid sequences, J. Comput. Phys. 276, 341-356 (2014)
