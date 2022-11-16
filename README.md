# Robust_penalized_BD_high_dim_GLM

============== Readme description for Matlab codes =============================

This is the readme description for producing results (tables, figures) related to the overdispersed Poisson responses in simulation studies in the paper "Robust estimation in regression and classification methods for large dimensional data". 
Matlab codes for all other tables and figures can be similarly and manually adjusted.

All Matlab codes are located in the same directory.

To implement the computation:

(1) run the Matlab main code: "demo_quasi_like.m". 

(2) This main code calls for other Matlab functions listed below:

GLM_BD_CD_penalized_parameter_estimate.m                  
GLM_BD_CD_penalized_parameter_estimate_SCAD.m             
GLM_BD_CD_penalized_parameter_estimate_SCAD_LLA.m         
GLM_BD_initial_weights_PCR.m                              
GLM_BD_parameter_estimate.m                               
G_function.m                                              
Poisson_cdf.m                                             
Poisson_initial_value_beta_0_in_robust_GLM_BD.m           
Poisson_pmf.m                                             
demo_quasi_like.m                                         
deri_1_G.m                                                
deri_G_1_function.m                                       
generate_Poisson.m                                        
hat_phi_of_overdispersed_Poisson_data_in_robust_BD.m      
p_star_function.m                                         
p_star_function_indefinite.m                              
q_1_q_2_BD.m                                              
q_1_q_2_BD_for_hat_V_n.m                                  
robust_GLM_BD_CD_penalized_parameter_estimate.m           
robust_GLM_BD_CD_penalized_parameter_estimate_SCAD.m      
robust_GLM_BD_CD_penalized_parameter_estimate_SCAD_LLA.m  
robust_GLM_BD_parameter_estimate.m                        
robust_GLM_BD_tune_lambda.m                               
robust_GLM_BD_tune_lambda_SCAD.m                          
robust_p_1_p_2_BD.m                                       
robust_p_1_p_2_BD_for_hat_V_n.m                           
robust_true_pp.m                                          
soft_thres.m                                              
true_penalty.m                                            
true_qq.m                                                 
weight_function_in_robust_BD_estimation.m                 
weight_response_BD.m               

================================= inputs used in the computation =================

For Table 1: run Matlab main code “demo_quasi_like.m”.
input sample size n (100; 200) = 100
input # (100; 500) of simulations = 500

For “study 1 (raw data without outliers), and non-robust estimation methods”:
– Input
  input study (1: without contamination; 2: with contamination) = 1
  input index_robust (0: without; 1: with) robust for Y with psi_c = 0
  input index_robust (0: without; 1: with) robust for X with weight = 0

For “study 1 (raw data without outliers), and robust estimation methods”:
– Input
  input study (1: without contamination; 2: with contamination) = 1
  input index_robust (0: without; 1: with) robust for Y with psi_c = 1
  input index_robust (0: without; 1: with) robust for X with weight = 1

----------------------------------------------------------------------------------

For Table 2: run Matlab main code “demo_quasi_like.m”.
input sample size n (100; 200) = 100
input # (100; 500) of simulations = 500

For “study 2 (contaminated data with outliers), and non-robust estimation methods”:
– Input
  input study (1: without contamination; 2: with contamination) = 2
  input index_robust (0: without; 1: with) robust for Y with psi_c = 0
  input index_robust (0: without; 1: with) robust for X with weight = 0

For “study 2 (contaminated data with outliers), and robust estimation methods”:
– Input
  input study (1: without contamination; 2: with contamination) = 2
  input index_robust (0: without; 1: with) robust for Y with psi_c = 1
  input index_robust (0: without; 1: with) robust for X with weight = 1
