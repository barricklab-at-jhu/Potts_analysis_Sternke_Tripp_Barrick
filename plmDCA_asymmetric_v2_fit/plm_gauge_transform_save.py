import numpy as np
import scipy.io as sio

def read_h_file(filename):
    '''
        Reads in matlab matrix of h coefficients
        
        Outputs numpy array of h coefficients
        
        Notes:
            If protein sequences are length L and using 21 letter amino acid alphabet
            h coeff array is L x 21
            
            Order of residues: -ACDEFGHIKLMNPQRSTUVWY
    '''
    
    mat = sio.loadmat(filename)
    hs = mat['h']
    
    return(hs)

def read_j_file(file, j_name, n_residues=21):
    '''
        Reads in matlab matrix of j coefficients
        
        Outputs numpy array of j coefficients
        
        Notes:
            If protein sequences are length L and using 21 letter amino acid alphabet
            j coeff array is an L x L x 21 x 21 matrix
            
            j_name is name of matrix in Matlab
            
            Order of residues: -ACDEFGHIKLMNPQRSTUVWY
    '''

    mat = sio.loadmat(file)
    j_in = np.transpose(mat[j_name])
    
    n_positions = int(-0.5*(1 - np.sqrt(8*len(j_in) + 1))) + 1
    
    j_out = np.zeros((n_positions, n_positions, n_residues, n_residues))
    j_count = 0
    
    for i in range(n_positions):
        for j in range(n_positions):
            if i < j:
                j_out[i][j] = np.transpose(j_in[j_count])
                j_out[j][i] = j_in[j_count]
                j_count += 1
                
    return(j_out)

def h_zero_sum_transformation(hs, js):
    '''
        This shifts the hs into the Ising gauge
        Zero-sum gauge shifts as much of the energy into the fields as available
        
        gauge transformation equation taken from:
        Ekeberg et al. "Fast pseudolikelihood maximization for direct-coupling analysis 
        of protein structure from many homologous amino-acid sequences." Journal of Computational Physics. 2014.
        
        Equation 24 in paper
        
        This shifts the hs into the zero-sum gauge. Jtrans matrix from Matlab code contains Js in zero-sum gauge.
        Zero-sum gauge shifts as much of the energy into the fields as available
    '''

    h_term_h_gauge_transform = np.zeros_like(hs)

    for i in range(len(hs)):
        mean_hs_at_i = np.mean(hs[i])
        h_term_h_gauge_transform[i] = hs[i] - mean_hs_at_i

    j_term_h_gauge_transform = np.zeros_like(hs)

    for position_i in range(len(js)):
        for position_j in range(len(js)):
            if position_i != position_j:
                full_jij_matrix_mean = np.mean(js[position_i][position_j])
                for res_a in range(len(j_term_h_gauge_transform[0])):
                    res_a_jij_mean = np.mean(js[position_i][position_j][res_a, :])
                    j_term_h_gauge_transform[position_i][res_a] += (res_a_jij_mean - full_jij_matrix_mean)

    h_hat = h_term_h_gauge_transform + j_term_h_gauge_transform
    
    return(h_hat)

#Sets file names for matrices
h_file = 'h.mat'
j_temp1_file = 'Jtemp1.mat'
j_temp2_file = 'Jtemp2.mat'
j_zero_sum_file = 'J.mat'

#Reads in all matrices
h_matrix_in = read_h_file(h_file)
j_temp1 = read_j_file(j_temp1_file, 'Jtemp1')
j_temp2 = read_j_file(j_temp2_file, 'Jtemp2')
j_zero_sum = read_j_file(j_zero_sum_file, 'J')

#Transforms h coefficients to zero-sum gague
h_hat_temp1 = h_zero_sum_transformation(h_matrix_in, j_temp1)
h_hat_temp2 = h_zero_sum_transformation(h_matrix_in, j_temp2)
h_zero_sum = (h_hat_temp1 + h_hat_temp2) / 2

#Saves coefficient matrices as numpy arrays
np.save('h_plm_zero_sum', h_zero_sum)
np.save('j_plm_zero_sum', j_zero_sum)
