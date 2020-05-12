## First attempt ##

This contains all the synthetic data generated in the first attempt described in Xiao et al.

## Some key files: 

* **mean_var_obs.txt**: Means and variances of the 240 reference time series.

* **Round 1 parameters 10000 ecm.txt**: Initial parameters, all parameter values are sampled from corresponding initial distributions using "runif" command in R. (See **Prior.R**)

* **Round 10 parameters 10000 all 3.txt**: Final parameters, the column means are taken as the final estimates of the parameters. (See **Results.R**) 

## How can the synthetic data be reproduced? ##

* Initially, sample 10000 values for each parameter from its initial distribution in **Prior.R**, form them into a 10000x6 table so we have 10000 vectors of parameter values that can be substituted back to the PDE solver. This table is written in "**Round 1 parameters 10000 ecm.txt**"

* Open the file "**Automatic.R**" in the folder to read in the functions, then open **Results.R** and use function **bcd()** to obtain the Bhattacharyya distance of each parameter vector in relation to the reference one, 10000 singular values are written in "**B-C distance ecm r1.txt**". 

* Now use function **abc_bcd()** to carry out the ABC-BCD optimization scheme. Read in "**Round 1 parameters 10000 ecm.txt**" and "**B-C distance ecm r1.txt**" as arguments "paras" and "ss_mat". The result of the function is written in "**Round 2 parameters 10000 ecm.txt**"

A full round for the evaluation of ECM profile is complete. After 3 rounds are carried out for the evaluation of ECM profile, we extract the sample of parameter eta, (4th column in "**Round 4 parameters 10000 ecm.txt**"), set it to be the 4th column of "**Round 1 parameters 10000 ecm_mde.txt**", sample the other parameter values from their initial distribution. Function **bcd()** needs to be reloeaded with 

After 5 rounds are carried out for the evaluation of MDE profile, we extract the samples of parameter eta, dm, alpha (4th, 5th, 6th column of "**Round 6 parameters 10000 ecm_mde.txt**"), set them to be the 4th, 5th, 6th column of "**Round 1 parameters 10000 all 3.txt**", sample the other parameter values from their initial distribution. 

After 9 rounds are carried out for the evaluation of tumour cells profile. The avarage of each column in "**Round 10 parameters 10000 all 3.txt**" can be calculated and taken as the final estimations of the parameters. 

Commands written in **Results.R** can be ran line by line to reproduce the results. 
