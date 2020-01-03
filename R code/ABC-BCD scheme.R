# R code of ABC-BCD optimization scheme described in the manuscript, and average B-C distance at the end of different rounds.
# (An obvious decreasing trend can be observed in the average B-C distance, our scheme is on the right track of minimizing B-C distance 
# and obtain the parameter values that can give the best fit to the synthetic data.)

# Used in: Xiao et al. "Calibrating models of cancer invasion and metastasis: parameter optimization using
# Approximate Bayesian Computation."

# Author: Yunchen XIAO

# A function that carries out all the steps in the ABC-BCD scheme described in the manuscript.
# It takes two input arguments, the array of Bhattacharya distance results and the array of corresponding parameter vectors, 
# at the end, it returns the parameter values for the next round of optimization. All detailed steps were kept in the function 
# for our readers to have a clear understanding on the function. 

Rejcon_bcd <- function(ss_mat,paras) {
  ss_mat <- as.matrix(ss_mat) 
  # Set the Bhattacharya distance array to be a matrix instead of a data frame.
  
  invalid <- vector()
  # Create an empty vector to collect the index of invalid terms.
  
  for (j in 1:length(ss_mat[1,])) {
    invalid_sep <- which(is.na(ss_mat[,j]) == TRUE)
    invalid <- c(invalid,invalid_sep)
  } # The indices of invalid terms. 
  
  invalid_index<-unique(invalid) 
  # Uniqueness checking, make sure each index only appears once.
  
  if(length(invalid_index) == 0) {
    ss_mat_valid <- ss_mat
  } else {
    ss_mat_valid <- ss_mat[-invalid_index,] 
    # Valid Bhattacharya distance results.
  }

  
  wt <- 1/(ss_mat_valid[,length(ss_mat_valid[1,])]^(1/2)) 
  # Weights of the valid B-C distance results, calculated by 1/sqrt(bcd).
  
  ss_mat_valid_wt <- cbind(ss_mat_valid,wt) 
  # Valid B-C results + Weights, one more column in the information table.
  
  resamp_prob <- rep(0,length(ss_mat_valid_wt[,length(ss_mat_valid_wt[1,])])) 
  # Create an empty vector to store the resampling probabilities of the parameter vectors 
  # that can return a valid B-C distance result. 
  
  for (i in 1:length(ss_mat_valid_wt[,length(ss_mat_valid_wt[1,])])) {
    if (ss_mat_valid_wt[i,length(ss_mat_valid_wt[1,])] == min(ss_mat_valid_wt[,length(ss_mat_valid_wt[1,])])) {
      resamp_prob[i] <- 0
    } else if (ss_mat_valid_wt[i,length(ss_mat_valid_wt[1,])] == max(ss_mat_valid_wt[,length(ss_mat_valid_wt[1,])])) {
      resamp_prob[i] <- 1
    } else {
      resamp_prob[i] <- (ss_mat_valid_wt[i,length(ss_mat_valid_wt[1,])]-min(ss_mat_valid_wt[,length(ss_mat_valid_wt[1,])]))/(max(ss_mat_valid_wt[,length(ss_mat_valid_wt[1,])])-min(ss_mat_valid_wt[,length(ss_mat_valid_wt[1,])]))
    }
  } 
  # Calculation of the resampling probabilities, see the manuscript for more details.
  
  ss_mat_valid_wt_prob <- cbind(ss_mat_valid_wt,resamp_prob) 
  # Valid B-C results + Weight + Resampling probabilities, one more column in the information matrix.
  
  resamp_ind <- sample(ss_mat_valid_wt_prob[,1],size = length(paras[,1]), 
                     replace = TRUE, prob = ss_mat_valid_wt_prob[,length(ss_mat_valid_wt_prob[1,])]) 
  # Resample the indices of the parameter vectors based on their resampling probabilities. 
  
  paras_nr_unperturbed<-paras[resamp_ind,] 
  # Resampled parameter vectors, without perturbation.
  
  paras_nr_perturbed <- matrix(0,nrow = nrow(paras),ncol = ncol(paras)) 
  # An empty matrix to store the perturbed parameter values for the next round.
  
  for (i in 1:length(paras[1,])) {
    for (j in 1:length(paras[,1])){
      h <- sqrt(1-0.05^2)
      paras_nr_perturbed[j,i] <- rnorm(1,h*paras_nr_unperturbed[j,i]+(1-h)*mean(paras_nr_unperturbed[,i]),
                                     0.05*sd(paras_nr_unperturbed[,i]))
    }
  }  
  # Perturb every single parameter value in "paras_nr_unperturbed", these parameter values will be the ones 
  # being optimized in the next round.
  
  paras_nr_perturbed<-as.data.frame(paras_nr_perturbed) 
  # Change the format back to a dataframe. 
  
  
  return(paras_nr_perturbed) 
  # Parameter values for the next round as the result of the function. 

}
########################################################################


bcd_ecm_r1<-read.table("B-C distance ecm r1.txt",sep="")
# B-C distances of parameter vectors in the first round of ECM profile evaluation.
paras_ori_ecm<-read.table("Round 1 parameters 10000 ecm.txt",sep="",header = TRUE)
# Parameter vectors used in the first round of ECM profile evaluation.
ind_bcd_ecm_r1_nan<-which(is.na(bcd_ecm_r1[,2]) == TRUE)
bcd_ecm_r1_nonan<-bcd_ecm_r1[-ind_bcd_ecm_r1_nan,]
mean(bcd_ecm_r1_nonan[,2]) 
# Average B-C distance: 13.70742

bcd_ecm_r2_p1<-read.table("B-C distance ecm r2 p1.txt",sep="")
bcd_ecm_r2_p2<-read.table("B-C distance ecm r2 p2.txt",sep="")
bcd_ecm_r2<-rbind(bcd_ecm_r2_p1,bcd_ecm_r2_p2)
# B-C distances of parameter vectors in the second round of ECM profile evaluation.
paras_r2_ecm<-read.table("Round 2 parameters 10000 ecm.txt",sep="",header = TRUE)
# Parameter vectors used in the second round of ECM profile evaluation.
ind_bcd_ecm_r2_nan<-which(bcd_ecm_r2[,2] == "NaN")
bcd_ecm_r2_nonan<-bcd_ecm_r2[-ind_bcd_ecm_r2_nan,]
mean(bcd_ecm_r2_nonan[,2]) 
# Average B-C distance: 8.691757

bcd_ecm_r3_p1<-read.table("B-C distance ecm r3 p1.txt",sep="")
bcd_ecm_r3_p2<-read.table("B-C distance ecm r3 p2.txt",sep="")
bcd_ecm_r3<-rbind(bcd_ecm_r3_p1,bcd_ecm_r3_p2)
# B-C distances of parameter vectors in the third round of ECM profile evaluation.
paras_r3_ecm<-read.table("Round 3 parameters 10000 ecm.txt",sep="")
# Parameter vectors used in the third round of ECM profile evaluation.
ind_bcd_ecm_r3_nan<-which(bcd_ecm_r3[,2] == "NaN")
bcd_ecm_r3_nonan<-bcd_ecm_r3[-ind_bcd_ecm_r3_nan,]
mean(bcd_ecm_r3_nonan[,2]) 
# Average B-C distance: 4.604086

paras_post_ecm<-read.table("Round 4 parameters 10000 ecm.txt",sep="")
# Parameter values at the end of the ECM profile evaluation. 

bcd_ecm_mde_r1_p1<-read.table("B-C distance ecm mde r1 p1.txt",sep="")
bcd_ecm_mde_r1_p2<-read.table("B-C distance ecm mde r1 p2.txt",sep="")
bcd_ecm_mde_r1_p3<-read.table("B-C distance ecm mde r1 p3.txt",sep="")
bcd_ecm_mde_r1_p4<-read.table("B-C distance ecm mde r1 p4.txt",sep="")
bcd_ecm_mde_r1<-rbind(bcd_ecm_mde_r1_p1,bcd_ecm_mde_r1_p2,bcd_ecm_mde_r1_p3,bcd_ecm_mde_r1_p4)
# B-C distances of parameter vectors in the first round of MDE profile evaluation.
paras_ori_ecm_mde<-read.table("Round 1 parameters 10000 ecm_mde.txt",sep="")
# Parameter vectors used in the first round of MDE profile evaluation.
ind_bcd_ecm_mde_r1_nan<-which(bcd_ecm_mde_r1[,2] == "NaN")
bcd_ecm_mde_r1_nonan<-bcd_ecm_mde_r1[-ind_bcd_ecm_mde_r1_nan,]
mean(bcd_ecm_mde_r1_nonan[,2]) 
# Average B-C distance: 31.85293

bcd_ecm_mde_r2_p1<-read.table("B-C distance ecm mde r2 p1.txt",sep="")
bcd_ecm_mde_r2_p2<-read.table("B-C distance ecm mde r2 p2.txt",sep="")
bcd_ecm_mde_r2_p3<-read.table("B-C distance ecm mde r2 p3.txt",sep="")
bcd_ecm_mde_r2_p4<-read.table("B-C distance ecm mde r2 p4.txt",sep="")
bcd_ecm_mde_r2<-rbind(bcd_ecm_mde_r2_p1,bcd_ecm_mde_r2_p2,bcd_ecm_mde_r2_p3,bcd_ecm_mde_r2_p4)
# B-C distances of parameter vectors in the second round of MDE profile evaluation.
paras_r2_ecm_mde<-read.table("Round 2 parameters 10000 ecm_mde.txt",sep="")
# Parameter vectors used in the second round of MDE profile evaluation.
ind_bcd_ecm_mde_r2_nan<-which(bcd_ecm_mde_r2[,2] == "NaN")
bcd_ecm_mde_r2_nonan<-bcd_ecm_mde_r2[-ind_bcd_ecm_mde_r2_nan,]
mean(bcd_ecm_mde_r2_nonan[,2]) 
# Average B-C distance: 22.96

bcd_ecm_mde_r3_p1<-read.table("B-C distance ecm mde r3 p1.txt",sep="")
bcd_ecm_mde_r3_p2<-read.table("B-C distance ecm mde r3 p2.txt",sep="")
bcd_ecm_mde_r3_p3<-read.table("B-C distance ecm mde r3 p3.txt",sep="")
bcd_ecm_mde_r3_p4<-read.table("B-C distance ecm mde r3 p4.txt",sep="")
bcd_ecm_mde_r3<-rbind(bcd_ecm_mde_r3_p1,bcd_ecm_mde_r3_p2,bcd_ecm_mde_r3_p3,bcd_ecm_mde_r3_p4)
# B-C distances of parameter vectors in the third round of MDE profile evaluation.
paras_r3_ecm_mde<-read.table("Round 3 parameters 10000 ecm_mde.txt",sep="")
# Parameter vectors used in the third round of MDE profile evaluation.
ind_bcd_ecm_mde_r3_nan<-which(bcd_ecm_mde_r3[,2] == "NaN")
bcd_ecm_mde_r3_nonan<-bcd_ecm_mde_r3[-ind_bcd_ecm_mde_r3_nan,]
mean(bcd_ecm_mde_r3_nonan[,2]) 
# Average B-C distance: 13.95495

bcd_ecm_mde_r4_p1<-read.table("B-C distance ecm mde r4 p1.txt",sep="")
bcd_ecm_mde_r4_p2<-read.table("B-C distance ecm mde r4 p2.txt",sep="")
bcd_ecm_mde_r4_p3<-read.table("B-C distance ecm mde r4 p3.txt",sep="")
bcd_ecm_mde_r4_p4<-read.table("B-C distance ecm mde r4 p4.txt",sep="")
bcd_ecm_mde_r4<-rbind(bcd_ecm_mde_r4_p1,bcd_ecm_mde_r4_p2,bcd_ecm_mde_r4_p3,bcd_ecm_mde_r4_p4)
# B-C distances of parameter vectors in the fourth round of MDE profile evaluation.
paras_r4_ecm_mde<-read.table("Round 4 parameters 10000 ecm_mde.txt",sep="")
# Parameter vectors used in the fourth round of MDE profile evaluation.
ind_bcd_ecm_mde_r4_nan<-which(bcd_ecm_mde_r4[,2] == "NaN")
bcd_ecm_mde_r4_nonan<-bcd_ecm_mde_r4[-ind_bcd_ecm_mde_r4_nan,]
mean(bcd_ecm_mde_r4_nonan[,2]) 
# Average B-C distance: 6.906569

bcd_ecm_mde_r5_p1<-read.table("B-C distance ecm mde r5 p1.txt",sep="")
bcd_ecm_mde_r5_p2<-read.table("B-C distance ecm mde r5 p2.txt",sep="")
bcd_ecm_mde_r5_p3<-read.table("B-C distance ecm mde r5 p3.txt",sep="")
bcd_ecm_mde_r5_p4<-read.table("B-C distance ecm mde r5 p4.txt",sep="")
bcd_ecm_mde_r5<-rbind(bcd_ecm_mde_r5_p1,bcd_ecm_mde_r5_p2,bcd_ecm_mde_r5_p3,bcd_ecm_mde_r5_p4)
# B-C distances of parameter vectors in the fifth round of MDE profile evaluation.
paras_r5_ecm_mde<-read.table("Round 5 parameters 10000 ecm_mde.txt",sep="")
# Parameter vectors used in the fifth round of MDE profile evaluation.
ind_bcd_ecm_mde_r5_nan<-which(bcd_ecm_mde_r5[,2] == "NaN")
bcd_ecm_mde_r5_nonan<-bcd_ecm_mde_r5[-ind_bcd_ecm_mde_r5_nan,]
mean(bcd_ecm_mde_r5_nonan[,2]) 
# Average B-C distance: 2.648578

bcd_all3_r1_p1<-read.table("B-C distance all 3 r1 p1.txt",sep="")
bcd_all3_r1_p2<-read.table("B-C distance all 3 r1 p2.txt",sep="")
bcd_all3_r1_p3<-read.table("B-C distance all 3 r1 p3.txt",sep="")
bcd_all3_r1_p4<-read.table("B-C distance all 3 r1 p4.txt",sep="")
bcd_all3_r1<-rbind(bcd_all3_r1_p1,bcd_all3_r1_p2,bcd_all3_r1_p3,bcd_all3_r1_p4)
# B-C distances of parameter vectors in the first round of tumour cells profile evaluation.
paras_ori_all3<-read.table("Round 1 parameters 10000 all 3.txt",sep="")
# Parameter vectors used in the first round of tumour cells profile evaluation.
ind_bcd_all3_r1_nan<-which(bcd_all3_r1[,2] == "NaN")
bcd_all3_r1_nonan<-bcd_all3_r1[-ind_bcd_all3_r1_nan,]
mean(bcd_all3_r1_nonan[,2]) 
# Average B-C distance: 17.17443

bcd_all3_r2_p1<-read.table("B-C distance all 3 r2 p1.txt",sep="")
bcd_all3_r2_p2<-read.table("B-C distance all 3 r2 p2.txt",sep="")
bcd_all3_r2_p3<-read.table("B-C distance all 3 r2 p3.txt",sep="")
bcd_all3_r2_p4<-read.table("B-C distance all 3 r2 p4.txt",sep="")
bcd_all3_r2<-rbind(bcd_all3_r2_p1,bcd_all3_r2_p2,bcd_all3_r2_p3,bcd_all3_r2_p4)
# B-C distances of parameter vectors in the second round of tumour cells profile evaluation.
paras_r2_all3<-read.table("Round 2 parameters 10000 all 3.txt",sep="")
# Parameter vectors used in the second round of tumour cells profile evaluation.
ind_bcd_all3_r2_nan<-which(bcd_all3_r2[,2] == "NaN")
bcd_all3_r2_nonan<-bcd_all3_r2[-ind_bcd_all3_r2_nan,]
mean(bcd_all3_r2_nonan[,2]) 
# Average B-C distance: 11.79248

bcd_all3_r3_p1<-read.table("B-C distance all 3 r3 p1.txt",sep="")
bcd_all3_r3_p2<-read.table("B-C distance all 3 r3 p2.txt",sep="")
bcd_all3_r3_p3<-read.table("B-C distance all 3 r3 p3.txt",sep="")
bcd_all3_r3_p4<-read.table("B-C distance all 3 r3 p4.txt",sep="")
bcd_all3_r3<-rbind(bcd_all3_r3_p1,bcd_all3_r3_p2,bcd_all3_r3_p3,bcd_all3_r3_p4)
# B-C distances of parameter vectors in the third round of tumour cells profile evaluation.
paras_r3_all3<-read.table("Round 3 parameters 10000 all 3.txt",sep="")
# Parameter vectors used in the third round of tumour cells profile evaluation.
ind_bcd_all3_r3_nan<-which(bcd_all3_r3[,2] == "NaN")
bcd_all3_r3_nonan<-bcd_all3_r3[-ind_bcd_all3_r3_nan,]
mean(bcd_all3_r3_nonan[,2]) 
# Average B-C distance: 6.756307

bcd_all3_r4_p1<-read.table("B-C distance all 3 r4 p1.txt",sep="")
bcd_all3_r4_p2<-read.table("B-C distance all 3 r4 p2.txt",sep="")
bcd_all3_r4_p3<-read.table("B-C distance all 3 r4 p3.txt",sep="")
bcd_all3_r4_p4<-read.table("B-C distance all 3 r4 p4.txt",sep="")
bcd_all3_r4<-rbind(bcd_all3_r4_p1,bcd_all3_r4_p2,bcd_all3_r4_p3,bcd_all3_r4_p4)
# B-C distances of parameter vectors in the fourth round of tumour cells profile evaluation.
paras_r4_all3<-read.table("Round 4 parameters 10000 all 3.txt",sep="")
# Parameter vectors used in the fourth round of tumour cells profile evaluation.
ind_bcd_all3_r4_nan<-which(bcd_all3_r4[,2] == "NaN")
bcd_all3_r4_nonan<-bcd_all3_r4[-ind_bcd_all3_r4_nan,]
mean(bcd_all3_r4_nonan[,2]) 
# Average B-C distance: 3.834115

bcd_all3_r5_p1<-read.table("B-C distance all 3 r5 p1.txt",sep="")
bcd_all3_r5_p2<-read.table("B-C distance all 3 r5 p2.txt",sep="")
bcd_all3_r5_p3<-read.table("B-C distance all 3 r5 p3.txt",sep="")
bcd_all3_r5_p4<-read.table("B-C distance all 3 r5 p4.txt",sep="")
bcd_all3_r5<-rbind(bcd_all3_r5_p1,bcd_all3_r5_p2,bcd_all3_r5_p3,bcd_all3_r5_p4)
# B-C distances of parameter vectors in the fifth round of tumour cells profile evaluation.
paras_r5_all3<-read.table("Round 5 parameters 10000 all 3.txt",sep="")
# Parameter vectors used in the fifth round of tumour cells profile evaluation.
ind_bcd_all3_r5_nan<-which(bcd_all3_r5[,2] == "NaN")
bcd_all3_r5_nonan<-bcd_all3_r5[-ind_bcd_all3_r5_nan,]
mean(bcd_all3_r5_nonan[,2]) 
# Average B-C distance: 2.253531

bcd_all3_r6_p1<-read.table("B-C distance all 3 r6 p1.txt",sep="")
bcd_all3_r6_p2<-read.table("B-C distance all 3 r6 p2.txt",sep="")
bcd_all3_r6_p3<-read.table("B-C distance all 3 r6 p3.txt",sep="")
bcd_all3_r6_p4<-read.table("B-C distance all 3 r6 p4.txt",sep="")
bcd_all3_r6<-rbind(bcd_all3_r6_p1,bcd_all3_r6_p2,bcd_all3_r6_p3,bcd_all3_r6_p4)
# B-C distances of parameter vectors in the sixth round of tumour cells profile evaluation.
paras_r6_all3<-read.table("Round 6 parameters 10000 all 3.txt",sep="")
# Parameter vectors used in the sixth round of tumour cells profile evaluation.
ind_bcd_all3_r6_nan<-which(bcd_all3_r6[,2] == "NaN")
bcd_all3_r6_nonan<-bcd_all3_r6[-ind_bcd_all3_r6_nan,]
mean(bcd_all3_r6_nonan[,2]) 
# Average B-C distance: 1.347388

bcd_all3_r7_p1<-read.table("B-C distance all 3 r7 p1.txt",sep="")
bcd_all3_r7_p2<-read.table("B-C distance all 3 r7 p2.txt",sep="")
bcd_all3_r7_p3<-read.table("B-C distance all 3 r7 p3.txt",sep="")
bcd_all3_r7_p4<-read.table("B-C distance all 3 r7 p4.txt",sep="")
bcd_all3_r7<-rbind(bcd_all3_r7_p1,bcd_all3_r7_p2,bcd_all3_r7_p3,bcd_all3_r7_p4)
# B-C distances of parameter vectors in the seventh round of tumour cells profile evaluation.
paras_r7_all3<-read.table("Round 7 parameters 10000 all 3.txt",sep="")
# Parameter vectors used in the seventh round of tumour cells profile evaluation.
ind_bcd_all3_r7_nan<-which(bcd_all3_r7[,2] == "NaN")
bcd_all3_r7_nonan<-bcd_all3_r7[-ind_bcd_all3_r7_nan,]
mean(bcd_all3_r7_nonan[,2]) 
# Average B-C distance: 0.8424036

bcd_all3_r8_p1<-read.table("B-C distance all 3 r8 p1.txt",sep="")
bcd_all3_r8_p2<-read.table("B-C distance all 3 r8 p2.txt",sep="")
bcd_all3_r8_p3<-read.table("B-C distance all 3 r8 p3.txt",sep="")
bcd_all3_r8_p4<-read.table("B-C distance all 3 r8 p4.txt",sep="")
bcd_all3_r8<-rbind(bcd_all3_r8_p1,bcd_all3_r8_p2,bcd_all3_r8_p3,bcd_all3_r8_p4)
# B-C distances of parameter vectors in the eighth round of tumour cells profile evaluation.
paras_r8_all3<-read.table("Round 8 parameters 10000 all 3.txt",sep="")
# Parameter vectors used in the eighth round of tumour cells profile evaluation.
ind_bcd_all3_r8_nan<-which(bcd_all3_r8[,2] == "NaN")
bcd_all3_r8_nonan<-bcd_all3_r8[-ind_bcd_all3_r8_nan,]
mean(bcd_all3_r8_nonan[,2]) 
# Average B-C distance: 0.5700982

bcd_all3_r9_p1<-read.table("B-C distance all 3 r9 p1.txt",sep="")
bcd_all3_r9_p2<-read.table("B-C distance all 3 r9 p2.txt",sep="")
bcd_all3_r9_p3<-read.table("B-C distance all 3 r9 p3.txt",sep="")
bcd_all3_r9_p4<-read.table("B-C distance all 3 r9 p4.txt",sep="")
bcd_all3_r9<-rbind(bcd_all3_r9_p1,bcd_all3_r9_p2,bcd_all3_r9_p3,bcd_all3_r9_p4)
# B-C distances of parameter vectors in the ninth round of tumour cells profile evaluation.
paras_r9_all3<-read.table("Round 9 parameters 10000 all 3.txt",sep="")
# Parameter vectors used in the ninth round of tumour cells profile evaluation.
ind_bcd_all3_r9_nan<-which(bcd_all3_r9[,2] == "NaN")
bcd_all3_r9_nonan<-bcd_all3_r9[-ind_bcd_all3_r9_nan,]
mean(bcd_all3_r9_nonan[,2]) 
# Average B-C distance: 0.4080612
