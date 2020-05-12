#### ABC-BCD scheme ####################
#### Author: Yunchen Xiao  #############

paras <- read.table("Round 1 parameters 10000 ecm.txt",sep="")

bcd <- function(paras) {
  
  # Mean_var_obsed table
  mean_var_obs <- read.table("mean_var_obs.txt",sep="")
  
  # Space
  l1 <- 0
  l2 <- 1
  x11 <- seq(l1,l2,length = 80)
  n.x11 <- length(x11)
  h <- x11[2] - x11[1]
  
  # Time
  T <- 10
  dt <- 0.001
  t <- seq(0, T, by = 0.001)
  
  # Paras
  dn <- paras[,1]
  gamma <- paras[,2]
  ita <- paras[,3]
  dm <- paras[,4]
  alpha <- paras[,5]
  r <- paras[,6]
  
  beta <- 0
  eps <- 0.01
  
  # Initial condition
  n0 <- rep(0, length(x11))
  
  for (i in 1:length(x11)) {
    if (x11[i]<=0.25) {
      n0[i]<-exp(-x11[i]^2/eps)
    } else {
      n0[i]<-0
    }
  }
  
  # Initial densities for the numerical solver
  
  f0 <- 1-0.5*n0
  
  m0 <- 0.5*n0
  
  # Initialise the empty vector to store the results
  
  bcd_tab <- vector()
  
  for (q in 1:length(paras[,1])) {
    
    n <- n0
    f <- f0
    m <- m0
    
    p <- 1
    
    res <- vector()
    
    while(p * dt <= T) {
      f[2:(length(x11)-1)] <- -ita[q] * dt * m[2:(length(x11)-1)] * f[2:(length(x11)-1)] + 
        f[2:(length(x11)-1)]
      
      m[2:(length(x11)-1)] <- dm[q] * (m[1:(length(x11)-2)] + m[3:length(x11)] - 2 * m[2:(length(x11)-1)]) * dt / (h ^ 2) +
        alpha[q] * n[2:(length(x11)-1)] * dt -  
        beta * m[2:(length(x11)-1)] * dt + m[2:(length(x11)-1)]
      
      n[2:(length(x11)-1)] <- dn[q] * (n[1:(length(x11)-2)] + n[3:length(x11)] - 2 * n[2:(length(x11)-1)]) * dt / (h ^ 2) -  
        gamma[q] * (n[3:length(x11)] - n[2:(length(x11)-1)]) * (f[3:length(x11)]-f[2:(length(x11)-1)]) * dt / (h ^ 2) - 
        gamma[q] * n[2:(length(x11)-1)] * (f[1:(length(x11)-2)] + f[3:length(x11)] - 2 * f[2:(length(x11)-1)]) * dt / (h ^ 2) + 
        r[q] * (1 - f[2:(length(x11)-1)] - n[2:(length(x11)-1)]) * n[2:(length(x11)-1)] * dt + n[2:(length(x11)-1)]
      
      #No flux boundary condition
      n[1] <- n[2]
      n[n.x11] <- n[n.x11 - 1]
      
      f[1] <- f[2]
      f[n.x11] <- f[n.x11 - 1]
      
      m[1] <- m[2]
      m[n.x11] <- m[n.x11 - 1]
      
      #res_full<-rbind(res_full,n,f,m)
      
      # Save the results at each t = 0.1*k (k as positive integers) time steps
      if(p %% 1000 == 0) {
        res <- rbind(res, n, f, m)
      }
      
      p <- p + 1
    }
    
    
    res_arr <- matrix(0,nrow = 30, ncol = 80)
    
    tc_ind <- seq(1,28, by = 3)
    ecm_ind <- seq(2,29, by = 3)
    mde_ind <- seq(3,30, by = 3)
    
    res_arr[1:10,] = res[tc_ind,]
    res_arr[11:20,] = res[ecm_ind,]
    res_arr[21:30,] = res[mde_ind,]
    
    mean_var = matrix(0,nrow= 240, ncol = 2); 
    
    for (i in 1:80) {
      mean_var[i,1] <- mean(res_arr[1:10,i])
      mean_var[i,2] <- var(res_arr[1:10,i])
      
      mean_var[i+80,1] <- mean(res_arr[11:20,i])
      mean_var[i+80,2] <- var(res_arr[11:20,i])
      
      mean_var[i+160,1] <- mean(res_arr[21:30,i])
      mean_var[i+160,2] <- var(res_arr[21:30,i])
    }
    
    bcd_vec <- vector();
    
    for (j in 1:240) {
      bcd <- 0.25*log(0.25*((mean_var[j,2]/mean_var_obs[j,2])+(mean_var_obs[j,2]/mean_var[j,2])+2))+0.25*(((mean_var[j,1]-mean_var_obs[j,1])^2)/(mean_var_obs[j,2]+mean_var[j,2]));
      bcd_vec <- c(bcd_vec,bcd);
    } 
    
    inv_index <- which(bcd_vec == "Inf")
    #inv_term <- length(inv_index)

    bcd_vec_2 <- bcd_vec;
    bcd_vec_2[inv_index] <- 0
    
    bcd_sum <- sum(bcd_vec_2)
    
    bcd_tab <- rbind(bcd_tab, bcd_sum)
  
  
  #bcd_tab <- as.data.frame(cbind(seq(1,length(paras[,1]), by = 1),bcd_tab))
  print(c(q, bcd_sum))
  }
  bcd_tab <- cbind(seq(1,length(bcd_tab), by = 1), bcd_tab)
  row.names(bcd_tab) <- seq(1,length(bcd_tab[,1]), by = 1)
  return(bcd_tab)
}

abc_bcd <- function(ss_mat,paras) {
  set.seed(123)
  ss_mat <- as.matrix(ss_mat) # Set the Bhattacharya distance array
  # to be a matrix instead of a data frame.
  
  invalid<-vector()
  
  for (j in 1:length(ss_mat[1,])) {
    invalid_sep <- which(is.na(ss_mat[,j]) == TRUE)
    invalid <- c(invalid,invalid_sep)
  } # The invalid terms. 
  
  invalid_index <- unique(invalid) # Uniqueness checking, 
  # make sure each term only appears once.
  
  if(length(invalid_index) == 0) {
    ss_mat_valid <- ss_mat
  } else {
    ss_mat_valid <- ss_mat[-invalid_index,] ## Valid Bhattacharya distance results.
  }
  
  
  wt <- 1/(ss_mat_valid[,length(ss_mat_valid[1,])]^(1/2)) # Weights of the valid B-C distance results.
  
  ss_mat_valid_wt <- cbind(ss_mat_valid,wt) # Valid B-C results + Weight
  
  resamp_prob <- rep(0,length(ss_mat_valid_wt[,length(ss_mat_valid_wt[1,])])) # Create an empty vector to store the resampling probabilities
  
  for (i in 1:length(ss_mat_valid_wt[,length(ss_mat_valid_wt[1,])])) {
    if (ss_mat_valid_wt[i,length(ss_mat_valid_wt[1,])] == min(ss_mat_valid_wt[,length(ss_mat_valid_wt[1,])])) {
      resamp_prob[i] <- 0
    } else if (ss_mat_valid_wt[i,length(ss_mat_valid_wt[1,])] == max(ss_mat_valid_wt[,length(ss_mat_valid_wt[1,])])) {
      resamp_prob[i] <- 1
    } else {
      resamp_prob[i] <- (ss_mat_valid_wt[i,length(ss_mat_valid_wt[1,])]-min(ss_mat_valid_wt[,length(ss_mat_valid_wt[1,])]))/(max(ss_mat_valid_wt[,length(ss_mat_valid_wt[1,])])-min(ss_mat_valid_wt[,length(ss_mat_valid_wt[1,])]))
    }
  } # Calculation of the resampling probabilities, see the manuscript for more details.
  
  ss_mat_valid_wt_prob <- cbind(ss_mat_valid_wt,resamp_prob) # Valid B-C results + Weight + Resampling probabilities.
  
  resamp_ind <- sample(ss_mat_valid_wt_prob[,1],size = length(paras[,1]), 
                     replace = TRUE, prob = ss_mat_valid_wt_prob[,length(ss_mat_valid_wt_prob[1,])]) # Resample the indices. 
  
  paras_nr_unperturbed <- paras[resamp_ind,] # Resampled parameter vectors, without perturbation.
  
  paras_nr_perturbed <- matrix(0,nrow = nrow(paras),ncol = ncol(paras)) # An empty matrix used to store the perturbed parameter values.
  
  for (i in 1:length(paras[1,])) {
    for (j in 1:length(paras[,1])){
      h = sqrt(1-0.05^2)
      paras_nr_perturbed[j,i] <- rnorm(1,h*paras_nr_unperturbed[j,i]+(1-h)*mean(paras_nr_unperturbed[,i]),
                                     0.05*sd(paras_nr_unperturbed[,i]))
    }
  }  # Perturbation, see the manuscript for more details.
  
  paras_nr_perturbed <- as.data.frame(paras_nr_perturbed) # Change the format back to a data frame. 
  
  
  return(paras_nr_perturbed) # Parameter values for the next round. 
  
}
