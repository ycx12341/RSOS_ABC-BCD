bcd_ecm_r1<-read.table("B-C distance ecm r1.txt",sep="")
paras_ori_ecm<-read.table("Round 1 parameters 10000 ecm.txt",sep="",header = TRUE)
ind_bcd_ecm_r1_nan<-which(is.na(bcd_ecm_r1[,2]) == TRUE)
bcd_ecm_r1_nonan<-bcd_ecm_r1[-ind_bcd_ecm_r1_nan,]
mean(bcd_ecm_r1_nonan[,2]) # Avg BCD: 13.70742

bcd_ecm_r2_p1<-read.table("B-C distance ecm r2 p1.txt",sep="")
bcd_ecm_r2_p2<-read.table("B-C distance ecm r2 p2.txt",sep="")
bcd_ecm_r2<-rbind(bcd_ecm_r2_p1,bcd_ecm_r2_p2)
paras_r2_ecm<-read.table("Round 2 parameters 10000 ecm.txt",sep="",header = TRUE)
ind_bcd_ecm_r2_nan<-which(bcd_ecm_r2[,2] == "NaN")
bcd_ecm_r2_nonan<-bcd_ecm_r2[-ind_bcd_ecm_r2_nan,]
mean(bcd_ecm_r2_nonan[,2]) # Avg BCD: 8.691757

bcd_ecm_r3_p1<-read.table("B-C distance ecm r3 p1.txt",sep="")
bcd_ecm_r3_p2<-read.table("B-C distance ecm r3 p2.txt",sep="")
bcd_ecm_r3<-rbind(bcd_ecm_r3_p1,bcd_ecm_r3_p2)
paras_r3_ecm<-read.table("Round 3 parameters 10000 ecm.txt",sep="")
ind_bcd_ecm_r3_nan<-which(bcd_ecm_r3[,2] == "NaN")
bcd_ecm_r3_nonan<-bcd_ecm_r3[-ind_bcd_ecm_r3_nan,]
mean(bcd_ecm_r3_nonan[,2]) # Avg BCD: 4.604086

paras_post_ecm<-read.table("Round 4 parameters 10000 ecm.txt",sep="")

bcd_ecm_mde_r1_p1<-read.table("B-C distance ecm mde r1 p1.txt",sep="")
bcd_ecm_mde_r1_p2<-read.table("B-C distance ecm mde r1 p2.txt",sep="")
bcd_ecm_mde_r1_p3<-read.table("B-C distance ecm mde r1 p3.txt",sep="")
bcd_ecm_mde_r1_p4<-read.table("B-C distance ecm mde r1 p4.txt",sep="")
bcd_ecm_mde_r1<-rbind(bcd_ecm_mde_r1_p1,bcd_ecm_mde_r1_p2,bcd_ecm_mde_r1_p3,bcd_ecm_mde_r1_p4)
paras_ori_ecm_mde<-read.table("Round 1 parameters 10000 ecm_mde.txt",sep="")
ind_bcd_ecm_mde_r1_nan<-which(bcd_ecm_mde_r1[,2] == "NaN")
bcd_ecm_mde_r1_nonan<-bcd_ecm_mde_r1[-ind_bcd_ecm_mde_r1_nan,]
mean(bcd_ecm_mde_r1_nonan[,2]) # Avg BCD: 31.85293

bcd_ecm_mde_r2_p1<-read.table("B-C distance ecm mde r2 p1.txt",sep="")
bcd_ecm_mde_r2_p2<-read.table("B-C distance ecm mde r2 p2.txt",sep="")
bcd_ecm_mde_r2_p3<-read.table("B-C distance ecm mde r2 p3.txt",sep="")
bcd_ecm_mde_r2_p4<-read.table("B-C distance ecm mde r2 p4.txt",sep="")
bcd_ecm_mde_r2<-rbind(bcd_ecm_mde_r2_p1,bcd_ecm_mde_r2_p2,bcd_ecm_mde_r2_p3,bcd_ecm_mde_r2_p4)
paras_r2_ecm_mde<-read.table("Round 2 parameters 10000 ecm_mde.txt",sep="")
ind_bcd_ecm_mde_r2_nan<-which(bcd_ecm_mde_r2[,2] == "NaN")
bcd_ecm_mde_r2_nonan<-bcd_ecm_mde_r2[-ind_bcd_ecm_mde_r2_nan,]
mean(bcd_ecm_mde_r2_nonan[,2]) # Avg BCD: 22.96

bcd_ecm_mde_r3_p1<-read.table("B-C distance ecm mde r3 p1.txt",sep="")
bcd_ecm_mde_r3_p2<-read.table("B-C distance ecm mde r3 p2.txt",sep="")
bcd_ecm_mde_r3_p3<-read.table("B-C distance ecm mde r3 p3.txt",sep="")
bcd_ecm_mde_r3_p4<-read.table("B-C distance ecm mde r3 p4.txt",sep="")
bcd_ecm_mde_r3<-rbind(bcd_ecm_mde_r3_p1,bcd_ecm_mde_r3_p2,bcd_ecm_mde_r3_p3,bcd_ecm_mde_r3_p4)
paras_r3_ecm_mde<-read.table("Round 3 parameters 10000 ecm_mde.txt",sep="")
ind_bcd_ecm_mde_r3_nan<-which(bcd_ecm_mde_r3[,2] == "NaN")
bcd_ecm_mde_r3_nonan<-bcd_ecm_mde_r3[-ind_bcd_ecm_mde_r3_nan,]
mean(bcd_ecm_mde_r3_nonan[,2]) # Avg BCD: 13.95495

bcd_ecm_mde_r4_p1<-read.table("B-C distance ecm mde r4 p1.txt",sep="")
bcd_ecm_mde_r4_p2<-read.table("B-C distance ecm mde r4 p2.txt",sep="")
bcd_ecm_mde_r4_p3<-read.table("B-C distance ecm mde r4 p3.txt",sep="")
bcd_ecm_mde_r4_p4<-read.table("B-C distance ecm mde r4 p4.txt",sep="")
bcd_ecm_mde_r4<-rbind(bcd_ecm_mde_r4_p1,bcd_ecm_mde_r4_p2,bcd_ecm_mde_r4_p3,bcd_ecm_mde_r4_p4)
paras_r4_ecm_mde<-read.table("Round 4 parameters 10000 ecm_mde.txt",sep="")
ind_bcd_ecm_mde_r4_nan<-which(bcd_ecm_mde_r4[,2] == "NaN")
bcd_ecm_mde_r4_nonan<-bcd_ecm_mde_r4[-ind_bcd_ecm_mde_r4_nan,]
mean(bcd_ecm_mde_r4_nonan[,2]) # Avg BCD: 6.906569

bcd_ecm_mde_r5_p1<-read.table("B-C distance ecm mde r5 p1.txt",sep="")
bcd_ecm_mde_r5_p2<-read.table("B-C distance ecm mde r5 p2.txt",sep="")
bcd_ecm_mde_r5_p3<-read.table("B-C distance ecm mde r5 p3.txt",sep="")
bcd_ecm_mde_r5_p4<-read.table("B-C distance ecm mde r5 p4.txt",sep="")
bcd_ecm_mde_r5<-rbind(bcd_ecm_mde_r5_p1,bcd_ecm_mde_r5_p2,bcd_ecm_mde_r5_p3,bcd_ecm_mde_r5_p4)
paras_r5_ecm_mde<-read.table("Round 5 parameters 10000 ecm_mde.txt",sep="")
ind_bcd_ecm_mde_r5_nan<-which(bcd_ecm_mde_r5[,2] == "NaN")
bcd_ecm_mde_r5_nonan<-bcd_ecm_mde_r5[-ind_bcd_ecm_mde_r5_nan,]
mean(bcd_ecm_mde_r5_nonan[,2]) # Avg BCD: 2.648578

bcd_all3_r1_p1<-read.table("B-C distance all 3 r1 p1.txt",sep="")
bcd_all3_r1_p2<-read.table("B-C distance all 3 r1 p2.txt",sep="")
bcd_all3_r1_p3<-read.table("B-C distance all 3 r1 p3.txt",sep="")
bcd_all3_r1_p4<-read.table("B-C distance all 3 r1 p4.txt",sep="")
bcd_all3_r1<-rbind(bcd_all3_r1_p1,bcd_all3_r1_p2,bcd_all3_r1_p3,bcd_all3_r1_p4)
paras_ori_all3<-read.table("Round 1 parameters 10000 all 3.txt",sep="")
ind_bcd_all3_r1_nan<-which(bcd_all3_r1[,2] == "NaN")
bcd_all3_r1_nonan<-bcd_all3_r1[-ind_bcd_all3_r1_nan,]
mean(bcd_all3_r1_nonan[,2]) # Avg BCD: 17.17443

bcd_all3_r2_p1<-read.table("B-C distance all 3 r2 p1.txt",sep="")
bcd_all3_r2_p2<-read.table("B-C distance all 3 r2 p2.txt",sep="")
bcd_all3_r2_p3<-read.table("B-C distance all 3 r2 p3.txt",sep="")
bcd_all3_r2_p4<-read.table("B-C distance all 3 r2 p4.txt",sep="")
bcd_all3_r2<-rbind(bcd_all3_r2_p1,bcd_all3_r2_p2,bcd_all3_r2_p3,bcd_all3_r2_p4)
paras_r2_all3<-read.table("Round 2 parameters 10000 all 3.txt",sep="")
ind_bcd_all3_r2_nan<-which(bcd_all3_r2[,2] == "NaN")
bcd_all3_r2_nonan<-bcd_all3_r2[-ind_bcd_all3_r2_nan,]
mean(bcd_all3_r2_nonan[,2]) # Avg BCD: 11.79248

bcd_all3_r3_p1<-read.table("B-C distance all 3 r3 p1.txt",sep="")
bcd_all3_r3_p2<-read.table("B-C distance all 3 r3 p2.txt",sep="")
bcd_all3_r3_p3<-read.table("B-C distance all 3 r3 p3.txt",sep="")
bcd_all3_r3_p4<-read.table("B-C distance all 3 r3 p4.txt",sep="")
bcd_all3_r3<-rbind(bcd_all3_r3_p1,bcd_all3_r3_p2,bcd_all3_r3_p3,bcd_all3_r3_p4)
paras_r3_all3<-read.table("Round 3 parameters 10000 all 3.txt",sep="")
ind_bcd_all3_r3_nan<-which(bcd_all3_r3[,2] == "NaN")
bcd_all3_r3_nonan<-bcd_all3_r3[-ind_bcd_all3_r3_nan,]
mean(bcd_all3_r3_nonan[,2]) # Avg BCD: 6.756307

bcd_all3_r4_p1<-read.table("B-C distance all 3 r4 p1.txt",sep="")
bcd_all3_r4_p2<-read.table("B-C distance all 3 r4 p2.txt",sep="")
bcd_all3_r4_p3<-read.table("B-C distance all 3 r4 p3.txt",sep="")
bcd_all3_r4_p4<-read.table("B-C distance all 3 r4 p4.txt",sep="")
bcd_all3_r4<-rbind(bcd_all3_r4_p1,bcd_all3_r4_p2,bcd_all3_r4_p3,bcd_all3_r4_p4)
paras_r4_all3<-read.table("Round 4 parameters 10000 all 3.txt",sep="")
ind_bcd_all3_r4_nan<-which(bcd_all3_r4[,2] == "NaN")
bcd_all3_r4_nonan<-bcd_all3_r4[-ind_bcd_all3_r4_nan,]
mean(bcd_all3_r4_nonan[,2]) # Avg BCD: 3.834115

bcd_all3_r5_p1<-read.table("B-C distance all 3 r5 p1.txt",sep="")
bcd_all3_r5_p2<-read.table("B-C distance all 3 r5 p2.txt",sep="")
bcd_all3_r5_p3<-read.table("B-C distance all 3 r5 p3.txt",sep="")
bcd_all3_r5_p4<-read.table("B-C distance all 3 r5 p4.txt",sep="")
bcd_all3_r5<-rbind(bcd_all3_r5_p1,bcd_all3_r5_p2,bcd_all3_r5_p3,bcd_all3_r5_p4)
paras_r5_all3<-read.table("Round 5 parameters 10000 all 3.txt",sep="")
ind_bcd_all3_r5_nan<-which(bcd_all3_r5[,2] == "NaN")
bcd_all3_r5_nonan<-bcd_all3_r5[-ind_bcd_all3_r5_nan,]
mean(bcd_all3_r5_nonan[,2]) # Avg BCD: 2.253531

bcd_all3_r6_p1<-read.table("B-C distance all 3 r6 p1.txt",sep="")
bcd_all3_r6_p2<-read.table("B-C distance all 3 r6 p2.txt",sep="")
bcd_all3_r6_p3<-read.table("B-C distance all 3 r6 p3.txt",sep="")
bcd_all3_r6_p4<-read.table("B-C distance all 3 r6 p4.txt",sep="")
bcd_all3_r6<-rbind(bcd_all3_r6_p1,bcd_all3_r6_p2,bcd_all3_r6_p3,bcd_all3_r6_p4)
paras_r6_all3<-read.table("Round 6 parameters 10000 all 3.txt",sep="")
ind_bcd_all3_r6_nan<-which(bcd_all3_r6[,2] == "NaN")
bcd_all3_r6_nonan<-bcd_all3_r6[-ind_bcd_all3_r6_nan,]
mean(bcd_all3_r6_nonan[,2]) # Avg BCD: 1.347388

bcd_all3_r7_p1<-read.table("B-C distance all 3 r7 p1.txt",sep="")
bcd_all3_r7_p2<-read.table("B-C distance all 3 r7 p2.txt",sep="")
bcd_all3_r7_p3<-read.table("B-C distance all 3 r7 p3.txt",sep="")
bcd_all3_r7_p4<-read.table("B-C distance all 3 r7 p4.txt",sep="")
bcd_all3_r7<-rbind(bcd_all3_r7_p1,bcd_all3_r7_p2,bcd_all3_r7_p3,bcd_all3_r7_p4)
paras_r7_all3<-read.table("Round 7 parameters 10000 all 3.txt",sep="")
ind_bcd_all3_r7_nan<-which(bcd_all3_r7[,2] == "NaN")
bcd_all3_r7_nonan<-bcd_all3_r7[-ind_bcd_all3_r7_nan,]
mean(bcd_all3_r7_nonan[,2]) # Avg BCD: 0.8424036

bcd_all3_r8_p1<-read.table("B-C distance all 3 r8 p1.txt",sep="")
bcd_all3_r8_p2<-read.table("B-C distance all 3 r8 p2.txt",sep="")
bcd_all3_r8_p3<-read.table("B-C distance all 3 r8 p3.txt",sep="")
bcd_all3_r8_p4<-read.table("B-C distance all 3 r8 p4.txt",sep="")
bcd_all3_r8<-rbind(bcd_all3_r8_p1,bcd_all3_r8_p2,bcd_all3_r8_p3,bcd_all3_r8_p4)
paras_r8_all3<-read.table("Round 8 parameters 10000 all 3.txt",sep="")
ind_bcd_all3_r8_nan<-which(bcd_all3_r8[,2] == "NaN")
bcd_all3_r8_nonan<-bcd_all3_r8[-ind_bcd_all3_r8_nan,]
mean(bcd_all3_r8_nonan[,2]) # Avg BCD: 0.5700982

bcd_all3_r9_p1<-read.table("B-C distance all 3 r9 p1.txt",sep="")
bcd_all3_r9_p2<-read.table("B-C distance all 3 r9 p2.txt",sep="")
bcd_all3_r9_p3<-read.table("B-C distance all 3 r9 p3.txt",sep="")
bcd_all3_r9_p4<-read.table("B-C distance all 3 r9 p4.txt",sep="")
bcd_all3_r9<-rbind(bcd_all3_r9_p1,bcd_all3_r9_p2,bcd_all3_r9_p3,bcd_all3_r9_p4)
paras_r9_all3<-read.table("Round 9 parameters 10000 all 3.txt",sep="")
ind_bcd_all3_r9_nan<-which(bcd_all3_r9[,2] == "NaN")
bcd_all3_r9_nonan<-bcd_all3_r9[-ind_bcd_all3_r9_nan,]
mean(bcd_all3_r9_nonan[,2]) # Avg BCD: 0.4080612

###### A function that carries out all the steps in the algorithm, it takes ###########
###### two input arguments, the array of Bhattacharya distance results and  ###########
###### the array of corresponding parameter vectors, at the end, it returns ###########
###### the parameters that will be used in the next round.                  ###########

Rejcon_bcd<-function(ss_mat,paras) {
  ss_mat<-as.matrix(ss_mat) # Set the Bhattacharya distance array
                            # to be a matrix instead of a data frame.
  
  invalid<-vector()
  
  for (j in 1:length(ss_mat[1,])) {
    invalid_sep<-which(is.na(ss_mat[,j]) == TRUE)
    invalid<-c(invalid,invalid_sep)
  } # The invalid terms. 
  
  invalid_index<-unique(invalid) # Uniqueness checking, 
                                 # make sure each term only appears once.
  
  if(length(invalid_index) == 0) {
    ss_mat_valid<-ss_mat
  } else {
    ss_mat_valid<-ss_mat[-invalid_index,] ## Valid Bhattacharya distance results.
  }

  
  wt<-1/(ss_mat_valid[,length(ss_mat_valid[1,])]^(1/2)) # Weights of the valid B-C distance results.
  
  ss_mat_valid_wt<-cbind(ss_mat_valid,wt) # Valid B-C results + Weight
  
  resamp_prob<-rep(0,length(ss_mat_valid_wt[,length(ss_mat_valid_wt[1,])])) # Create an empty vector to store the resampling probabilities
  
  for (i in 1:length(ss_mat_valid_wt[,length(ss_mat_valid_wt[1,])])) {
    if (ss_mat_valid_wt[i,length(ss_mat_valid_wt[1,])] == min(ss_mat_valid_wt[,length(ss_mat_valid_wt[1,])])) {
      resamp_prob[i]<-0
    } else if (ss_mat_valid_wt[i,length(ss_mat_valid_wt[1,])] == max(ss_mat_valid_wt[,length(ss_mat_valid_wt[1,])])) {
      resamp_prob[i]<-1
    } else {
      resamp_prob[i]<-(ss_mat_valid_wt[i,length(ss_mat_valid_wt[1,])]-min(ss_mat_valid_wt[,length(ss_mat_valid_wt[1,])]))/(max(ss_mat_valid_wt[,length(ss_mat_valid_wt[1,])])-min(ss_mat_valid_wt[,length(ss_mat_valid_wt[1,])]))
    }
  } # Calculation of the resampling probabilities, see the manuscript for more details.
  
  ss_mat_valid_wt_prob<-cbind(ss_mat_valid_wt,resamp_prob) # Valid B-C results + Weight + Resampling probabilities.
  
  resamp_ind<-sample(ss_mat_valid_wt_prob[,1],size = length(paras[,1]), 
                     replace = TRUE, prob = ss_mat_valid_wt_prob[,length(ss_mat_valid_wt_prob[1,])]) # Resample the indices. 
  
  paras_nr_unperturbed<-paras[resamp_ind,] # Resampled parameter vectors, without perturbation.
  
  paras_nr_perturbed<-matrix(0,nrow = nrow(paras),ncol = ncol(paras)) # An empty matrix used to store the perturbed parameter values.
  
  for (i in 1:length(paras[1,])) {
    for (j in 1:length(paras[,1])){
      h = sqrt(1-0.05^2)
      paras_nr_perturbed[j,i]<-rnorm(1,h*paras_nr_unperturbed[j,i]+(1-h)*mean(paras_nr_unperturbed[,i]),
                                     0.05*sd(paras_nr_unperturbed[,i]))
    }
  }  # Perturbation, see the manuscript for more details.
  
  paras_nr_perturbed<-as.data.frame(paras_nr_perturbed) # Change the format back to a dataframe. 
  
  
  return(paras_nr_perturbed) # Parameter values for the next round. 

}
########################################################################

