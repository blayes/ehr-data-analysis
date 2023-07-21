* Organization
  - There are two directories, each correspond to a simulation in the experiments section of the paper: sim_class_1 (linear decision boundary) and sim_class_2 (quadratic decision boundary). 
  - Each directory has three sub directories: code, data, and qsub.
  - Directory 'code' has files of all the R source code that was used in the analysis. 
  - Directory 'data' has (if any) simulated data that was used in the analysis. This directory may be empty or absent.
  - Directory 'qsub' has SGE files (.q) that were used to submit jobs on a SGE cluster. 
  
* Files
  - 'sim_data_corr.R' contains the code to simulate and partition the data for classification in first simulation.
  - 'sim_data_nolin_corr.R' contains the code to simulate and partition the data for classification in second simulation.
  - 'sim_data_nonlin_corr_reg.R' contains the code to simulate and partition the data for regression in second simulation.
  - 'analyze_result_class.R' contains the code for analyzing the results for classification of GP with the SE kernel, logistic regression and ridge regression and lasso regression with and without SE kernel-based covariates, SVM and KRR with spectrum and boundrange kernels, and competing methods and making tables.
  - 'analyze_result_reg.R' contains the code for analyzing the results for regression of GP with the SE kernel, logistic regression and ridge regression and lasso regression with and without SE kernel-based covariates, SVM and KRR with spectrum and boundrange kernels, and competing methods and making tables.
  - 'gp_logistic.R' contains the code for the GP classification with the SE kernel.
  - 'submit.R' contains the  code for the R code for submitting a job on the cluster. The files in 'qsub' directory use this file for running simulations.
  - The files ending in phecode do regression and classification with phecode and log(1 + phecode count) as predictors.
  
* Citation:
  If you use the code, then please cite the following paper:
  - Sanvesh Srivastava, Zongyi Xu, Yunyi Li, Nick Street, and Stephanie Gilbertson-White (2021+). Gaussian Process Regression and Classification using International Disease Classification Codes as Covariates. https://arxiv.org/abs/2108.01813

* Contact:
  Please email Sanvesh Srivastava (<sanvesh-srivastava@uiowa.edu>) if you have any questions related to the code.

* Acknowledgment
  Office of Naval Research (ONR-BAA N000141812741) and the National Science Foundation (DMS-1854667/1854662).
  
* NOTE:
  We are unable to provide the code for real data analysis due to privacy concerns. 

