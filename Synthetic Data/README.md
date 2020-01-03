## Synthetic data ##

The folder "1st attempt" contains all the synthetic data generated in the first attempt described in Xiao et al. "Calibrating models of cancer invasion and metastasis: parameter optimization using Approximate Bayesian Computation."

## Some key files: 

* **mean_var_obs.txt**: Means and variances of the 900 reference time series.

* **Round 1 parameters 10000 ecm.txt**: Initial parameters, all parameter values sampled from corresponding initial distributions using "r" command in R.

* **Round 10 parameters 10000 all 3.txt**: Final parameters, the means of parameter columns were taken as the final estimates of the parameters. 

## How are the synthetic data produced? ##

* Initially, sample 10000 values for each parameter from its initial distribution, (as stated in the manuscript), form them into a 10000x6 table in a random order. This table is written into "**Round 1 parameters 10000 ecm.txt**"

* Use the file "**Bhattacharyya_distance.m**" in the MATLAB code folder to read the parameter table, obtain the Bhattacharyya distance of each parameter vector, 10000 singular values are written in "**B-C distance ecm r1.txt**". (Note that this procedure can be carried out parallelly, e.g. open "**Bhattacharyya_distance.m**" twice and change the index in the for loop, let the first file calculate the Bhattacharyya distances of vector 1-5000, and let the second one calculate the Bhattacharyya distances of vector 5001-10000.)

* Now open "**ABC-BCD scheme.R**" in the R code folder to carry out the ABC-BCD optimization scheme. Read in "**Round 1 parameters ecm.txt**" and "**B-C distance ecm r1.txt**" as arguments "paras" and "ss_mat" in the "Rejcon_bcd" function. The result of the function is written into "**Round 2 parameters 10000 ecm.txt**"

A full round for the evaluation of ECM profile is complete. After 3 rounds are carried out for the evaluation of ECM profile, we extract the sample of parameter eta, (4th column of "**Round 4 parameters 10000 ecm.txt**"), set it to be the 4th column of "**Round 1 parameters 10000 ecm_mde.txt**", sample the other parameter values from their initial distribution. 

After 5 rounds are carried out for the evaluation of MDE profile, we extract the samples of parameter eta, dm, alpha (4th, 5th, 6th column of "**Round 6 parameters 10000 ecm_mde.txt**"), set them to be the 4th, 5th, 6th column of "**Round 1 parameters 10000 all 3.txt**", sample the other parameter values from their initial distribution. 

After 9 rounds are carried out for the evaluation of tumour cells profile. The avarage of each column in "**Round 10 parameters 10000 all 3.txt**" can be calculated and taken as the final estimations of the parameters. 
