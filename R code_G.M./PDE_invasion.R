library(mgcv)
library(tictoc)

# 1D dimensionless space
n.x11 <- 50
x11 <- seq(0, 1, length = n.x11)
h <- x11[2] - x11[1]

# Time
T <- 5
dt <- 0.001
t <- seq(0, T, by = 0.001)

# Parameters
dn <- 0.01
gamma <- 0.05
ita <- 10
dm <- 0.01
alpha <- 0.1
r <- 5

beta <- 0
eps <- 0.01

# Initial condition
n0 <- rep(0, length(x11))

for (i in 1:50) {
  if (x11[i]<=0.25) {
    n0[i]<-exp(-x11[i]^2/eps)
  } else {
    n0[i]<-0
  }
}

n <- n0

f0 <- 1-0.5*n0
f <- f0

m0 <- 0.5*n0
m <- m0

# Store the initial values into the result matrix
res <- matrix(0, 3, n.x11)
res [1, ] <- n0
res [2, ] <- f0
res [3, ] <- m0

# res_full<-res

# Run the equation
p <- 1
while(p * dt <= T) {
  f[2:49] <- -ita * dt * m[2:49] * f[2:49] + 
                      f[2:49]

  m[2:49] <- dm * (m[1:48] + m[3:50] - 2 * m[2:49]) * dt / (h ^ 2) +
                      alpha * n[2:49] * dt -  
                      beta * m[2:49] * dt + m[2:49]
  
  n[2:49] <- dn * (n[1:48] + n[3:50] - 2 * n[2:49]) * dt / (h ^ 2) -  
                      gamma * (n[3:50] - n[2:49]) * (f[3:50]-f[2:49]) * dt / (h ^ 2) - 
                      gamma * n[2:49] * (f[1:48] + f[3:50] - 2 * f[2:49]) * dt / (h ^ 2) + 
                      r * (1 - f[2:49] - n[2:49]) * n[2:49] * dt + n[2:49]
  
  #No flux boundary condition
  n[1] <- n[2]
  n[n.x11] <- n[n.x11 - 1]
  
  f[1] <- f[2]
  f[n.x11] <- f[n.x11 - 1]
  
  m[1] <- m[2]
  m[n.x11] <- m[n.x11 - 1]
  
  #res_full<-rbind(res_full,n,f,m)
  
  # Save the results at each t = 0.1*k (k as positive integers) time steps
  if(p %% 100 == 0) {
    res <- rbind(res, n, f, m)
  }
  
  p <- p + 1
}

# Objective function. 
obj <- function(paras, grad.lhs_n, grad.rhs_dn, grad.rhs_gamma, grad.rhs_r, 
                grad.lhs_f, grad.rhs_ita, 
                grad.lhs_m, grad.rhs_dm, grad.rhs_alpha,fitted_val, obs_val) {
  
  disp1 <- sum((grad.lhs_f + paras[4] * grad.rhs_ita)^2)
  disp2 <- sum((grad.lhs_m - (paras[5] * grad.rhs_dm + paras[6] * grad.rhs_alpha))^2)
  disp3 <- sum((grad.lhs_n - (paras[1] * grad.rhs_dn - paras[2] * grad.rhs_gamma + paras[3] * grad.rhs_r))^2)
  
  disp4 <- sum((fitted_val-obs_val)^2)
  res <- sum(disp1,disp2,disp3,disp4)
  
  #if (debug) cat(dn, res, "\n")
  return(res)
}

#scale = 0.02

res_n <- res[seq(1,151,by = 3),]
res_f <- res[seq(2,152,by = 3),]
res_m <- res[seq(3,153,by = 3),] # Unperturbed data

est <- matrix(0, nrow = 1, ncol = 7)

library(doParallel)
cl <- makeCluster(detectCores() - 1) 
registerDoParallel(cl)

tic()
for (i in 1:10) {
 
  res_n_perturbed <- res_n + rnorm(51*50, 0, sd = 0.1*sd(res_n))
  res_f_perturbed <- res_f + rnorm(51*50, 0, sd = 0.1*sd(res_f))
  res_m_perturbed <- res_m + rnorm(51*50, 0, sd = 0.1*sd(res_m))
  
  dat_n <- data.frame(n = as.vector(res_n_perturbed), t = rep(seq(0,5, by = 0.1), times = n.x11), x11 = rep(x11, each = 51))

  dat_f <- data.frame(f = as.vector(res_f_perturbed), t = rep(seq(0,5, by = 0.1), times = n.x11), x11 = rep(x11, each = 51))

  dat_m <- data.frame(m = as.vector(res_m_perturbed), t = rep(seq(0,5, by = 0.1), times = n.x11), x11 = rep(x11, each = 51))

  suppressWarnings(spl <- gam(n ~ s(t, x11, bs = "ad"), data = dat_n)) 
  suppressWarnings(spl2 <- gam(f ~ s(t, x11, bs = "ad"), data = dat_f)) # GAM fitting * 3
  suppressWarnings(spl3 <- gam(m ~ s(t, x11, bs = "ad"), data = dat_m))
  
  # Smoothing GAM to n at t = 0. (Some plots to check how good the fitted values are)
  # plot(dat_n$n[dat_n$t==1]) # Perturbed data
  # lines(res[31,],col="green") # True solution
  # lines(predict(spl)[which(dat_n$t==1)],col="blue") # Fitted smoothing GAM to n at t = 0.
  
  # Smoothing GAM to n at t = 0.
  # plot(dat_f$f[dat_f$t==0]) # Perturbed data
  # lines(res[2,],col="green") # True solution
  # lines(predict(spl2)[which(dat_f$t==0)],col="blue") # Fitted smoothing GAM to n at t = 0.
  
  # plot(dat_f$f[dat_f$t==1]) # Perturbed data
  # lines(res[32,],col="green") # True solution
  # lines(predict(spl2)[which(dat_f$t==1)],col="blue")
  
  # Smoothing GAM to n at t = 0.
  # plot(dat_m$m[dat_m$t==0]) # Perturbed data
  # lines(res[3,],col="green") # True solution
  # lines(predict(spl3)[which(dat_m$t==0)],col="blue") # Fitted smoothing GAM to n at t = 0.
  
  # Seems like the GAM fitted the data quite well at this time point... 
  
  # Gradients matching:
  tp <- seq(0,5, by = 0.1)
  
  # f first:
  grad_lhs_f_data <- matrix(0,49,48)
  
  for (i in 1:49) {
    for (j in 1:48) {
      predict_f_t1 <- predict(spl2, newdata = data.frame(t = tp[i + 1] + dt, x11 = x11[j + 1]))
      predict_f_t2 <- predict(spl2, newdata = data.frame(t = tp[i + 1] - dt, x11 = x11[j + 1]))
      grad_lhs_f_data[i,j] <- (predict_f_t1 - predict_f_t2) / (2 * dt)
    }
  }
  
  # LHS from data
  
  # for (i in 1:49) {
  #   for (j in 1:48) {
  #     grad_lhs_f_data[i,j] <- (res_f_perturbed[i+2,j+1] - res_f_perturbed[i,j+1])/ (2* 0.1)
  #   }
  # }
  
  grad_rhs_f_ita <- matrix(0,49,48)
  
  for (i in 1:49) {
    for (j in 1:48) {
      predict_f <- predict(spl2, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1]))
      predict_m <- predict(spl3, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1]))
      grad_rhs_f_ita[i, j] <- predict_f*predict_m
    }
  }
  
  # disp_f <- sum((grad_lhs_f-grad_rhs_f_ita*(-10))^2) # Sum of squares: 0.375
  
  # Then let's go for m. 
  
  grad_lhs_m_data <- matrix(0,49,48)
  
  for (i in 1:49) {
    for (j in 1:48) {
      predict_m_t1 <- predict(spl3, newdata = data.frame(t = tp[i + 1] + dt, x11 = x11[j + 1]))
      predict_m_t2 <- predict(spl3, newdata = data.frame(t = tp[i + 1] - dt, x11 = x11[j + 1]))
      grad_lhs_m_data[i,j] <- (predict_m_t1 - predict_m_t2) / (2 * dt)
    }
  }
  
  # LHS data
  # for (i in 1:49) {
  #   for (j in 1:48) {
  #     grad_lhs_m_data[i,j] <- (res_m_perturbed[i+2,j+1] - res_m_perturbed[i,j+1])/ (2 * 0.1)
  #   }
  # }
  
  grad_rhs_m_dm <- matrix(0,49,48)
  
  for (i in 1:49) {
    for (j in 1:48) {
      predict_x1 <- predict(spl3, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1] - h))
      predict_x2 <- predict(spl3, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1]))
      predict_x3 <- predict(spl3, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1] + h))
      grad_rhs_m_dm[i, j] <- (predict_x3 + predict_x1 - 2 * predict_x2) / (h ^ 2)
    }
  }
  
  grad_rhs_m_alpha <- matrix(0,49,48)
  
  for (i in 1:49) {
    for (j in 1:48) {
      predict_n <- predict(spl, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1]))
      grad_rhs_m_alpha[i, j] <- predict_n
    }
  }
  
  # disp_m <- sum((grad_lhs_m-grad_rhs_m_dm*dm-grad_rhs_m_alpha*alpha)^2) # Sum of squares: 0.108
  
  # At last, for n. 
  
  grad_lhs_n_data <- matrix(0,49,48)
  
  for (i in 1:49) {
    for (j in 1:48) {
      predict_n_t1 <- predict(spl, newdata = data.frame(t = tp[i + 1] + dt, x11 = x11[j + 1]))
      predict_n_t2 <- predict(spl, newdata = data.frame(t = tp[i + 1] - dt, x11 = x11[j + 1]))
      grad_lhs_n_data[i,j] <- (predict_n_t1 - predict_n_t2) / (2 * dt)
    }
  }
  
  # for (i in 1:49) {
  #   for (j in 1:48) {
  #     grad_lhs_n_data[i,j] <- (res_n_perturbed[i+2,j+1] - res_n_perturbed[i,j+1])/ (2 * 0.1)
  #   }
  # }
  
  grad_rhs_n_dn <- matrix(0,49,48)
  
  for (i in 1:49) {
    for (j in 1:48) {
      predict_x1 <- predict(spl, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1] - h))
      predict_x2 <- predict(spl, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1]))
      predict_x3 <- predict(spl, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1] + h))
      grad_rhs_n_dn[i, j] <- (predict_x3 + predict_x1 - 2 * predict_x2) / (h ^ 2)
    }
  }
  
  grad_rhs_n_gamma <- matrix(0,49,48)
  
  for (i in 1:49) {
    for (j in 1:48) {
      predict_x1_n <- predict(spl, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1] - h))
      predict_x2_n <- predict(spl, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1]))
      predict_x3_n <- predict(spl, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1] + h))
      predict_n.dash <- (predict_x3_n - predict_x1_n) / (2 * h)
      
      predict_x1_f <- predict(spl2, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1] - h))
      predict_x2_f <- predict(spl2, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1]))
      predict_x3_f <- predict(spl2, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1] + h))
      predict_f.dash <- (predict_x3_f - predict_x1_f) / (2 * h)
      predict_f.ddash <- (predict_x3_f + predict_x1_f - 2 * predict_x2_f) / (h ^ 2)
      
      grad_rhs_n_gamma[i, j] <- predict_n.dash * predict_f.dash + + predict_x2_n * predict_f.ddash
    }
  }
  
  grad_rhs_n_r <- matrix(0,49,48)
  
  for (i in 1:49) {
    for (j in 1:48) {
      predict_n <- predict(spl, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1]))
      
      predict_f <- predict(spl2, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1]))
      
      grad_rhs_n_r[i, j] <- predict_n*(1-predict_n-predict_f)
    }
  }
  
  # grad_rhs_n <- dn*grad_rhs_n_dn-gamma*(grad_rhs_n_gamma.1+grad_rhs_n_gamma.2)+r*grad_rhs_n_r
  
  # grad_rhs_n_full <- matrix(0,49,48)
  
  # for (i in 1:49) {
  #  for (j in 1:48) {
  #    predict_n_x1 <- predict(spl, newdata = data.frame(t = tp[i + 1], x11 = x11[j+1] - h))
  #    predict_n_x2 <- predict(spl, newdata = data.frame(t = tp[i + 1], x11 = x11[j+1]))
  #    predict_n_x3 <- predict(spl, newdata = data.frame(t = tp[i + 1], x11 = x11[j+1] + h))
  
  #    predict_f_x1 <- predict(spl2, newdata = data.frame(t = tp[i + 1], x11 = x11[j+1] - h))
  #    predict_f_x2 <- predict(spl2, newdata = data.frame(t = tp[i + 1], x11 = x11[j+1]))
  #    predict_f_x3 <- predict(spl2, newdata = data.frame(t = tp[i + 1], x11 = x11[j+1] + h))
  #    
  #    predict_n.ddash = (predict_n_x3 + predict_n_x1 - 2* predict_n_x2)/ (h ^ 2)
  #    predict_n.dash = (predict_n_x3 - predict_n_x1) / (2 * h)
  
  #    predict_f.ddash = (predict_f_x3 + predict_f_x1 - 2* predict_f_x2)/ (h ^ 2)
  #    predict_f.dash = (predict_f_x3 - predict_f_x1) / (2 * h)
  
  #    rhs_n_dn <- dn * predict_n.ddash
  #    rhs_n_gamma <- gamma * (predict_n.dash * predict_f.dash + predict_n_x2 * predict_f.ddash)
  #    rhs_n_r <- r * (1 - predict_n_x2 - predict_f_x2) * predict_n_x2
  #    grad_rhs_n_full[i,j] <- rhs_n_dn - rhs_n_gamma + rhs_n_r
  #  }
  #}
  
  #grad_pde_lhs_n <- matrix(0,49,48)
  
  #for (i in 1:49) {
  #  for (j in 1:48) {
  #    grad_pde_lhs_n[i,j] <- (res_n[i+2, j+1]-res_n[i, j+1]) / (2*0.1)
  #  }
  #}
  
  #grad_pde_rhs_n <- matrix(0,49,48)
  
  #for (i in 1:49) {
  #  for (j in 1:48) {
  #    fst <- dn * (res_n[i+1, j+2] + res_n[i+1, j] - 2 * res_n[i+1, j+1]) / (h ^ 2)
  #    snd <- gamma * (res_f[i+1, j+2] - res_f[i+1, j]) * (res_n[i+1, j+2] - res_f[i+1, j]) / (4*(h^2))
  #    trd <- gamma * res_n[i+1, j+1] * (res_f[i+1, j+2] + res_f[i+1, j] - 2 * res_f[i+1, j+1]) / (h ^ 2)
  #    fth <- r * res_n[i+1, j+1] * (1 - res_n[i+1, j+1] - res_f[i+1, j+1])
  
  #    grad_pde_rhs_n[i,j] <- fst - snd - trd + fth
  #  }
  #}
  
  
  # disp_n <- sum((grad_lhs_n-grad_rhs_n_full)^2)
  
  
  
  start.values = c(0.02,0.1,10,20,0.02,0.2)
  
  opt <- optim(start.values, obj, grad.lhs_n = grad_lhs_n_data, grad.rhs_dn = grad_rhs_n_dn, grad.rhs_gamma = grad_rhs_n_gamma,
               grad.rhs_r = grad_rhs_n_r, grad.lhs_f = grad_lhs_f_data, grad.rhs_ita = grad_rhs_f_ita,
               grad.lhs_m = grad_lhs_m_data, grad.rhs_dm = grad_rhs_m_dm, grad.rhs_alpha = grad_rhs_m_alpha, obs_val = c(dat_n$n, dat_f$f, dat_m$m),
               fitted_val = c(spl$fitted.values, spl2$fitted.values, spl3$fitted.values))
  
  est <- rbind(est, c(opt$par,opt$value))
}
toc()
# Variance estimation ?

# Gradient matching other papers ? 

 start.values = c(0.03,0.15,15,30,0.03,0.3)

 optim(start.values, obj, grad.lhs_n = grad_lhs_n, grad.rhs_dn = grad_rhs_n_dn, 
             grad.rhs_gamma.1 = grad_rhs_n_gamma.1, grad.rhs_gamma.2 = grad_rhs_n_gamma.2,
             grad.rhs_r = grad_rhs_n_r, grad.lhs_f = grad_lhs_f, grad.rhs_ita = grad_rhs_f_ita,
             grad.lhs_m = grad_lhs_m, grad.rhs_dm = grad_rhs_m_dm, grad.rhs_alpha = grad_rhs_m_alpha)

 # Inference on the estimated values
 
 est_vals <- read.table("Estimated values.txt", sep="")
 
 dn_vals <- est_vals[2:101,1] 
 gamma_vals <- est_vals[2:101,2]
 r_vals <- est_vals[2:101,3]
 ita_vals <- est_vals[2:101,4]
 dm_vals <- est_vals[2:101,5]
 al_vals <- est_vals[2:101,6]
 
 est_vals_valid <- est_vals[2:101,]
 est_vals_valid <- as.matrix(est_vals_valid)
 est_vals_re <- data.frame(Res = as.vector(est_vals_valid), Paras = rep(c("dn", "gamma", "r", "ita", "dm", "alpha"), each = 100))
 
 par(mfrow = c(1,1))
 
 hist(dn_vals, prob = TRUE, xlab = "Estimated dn values", ylim = c(0,250), main = "Histogram of dn values, true value = 0.01")
 lines(density(dn_vals, adjust = 6), col = "blue",lwd = 2)
 abline(v=mean(dn_vals),col="red", lwd = 2)
 abline(v=0.01, col = "green", lwd = 2)
 legend(x=0.018,y=250,c("Density",
                       "Mean", "True value"),
        cex=1.2,col=c("blue","red","green"),lty = 1)
 
 hist(gamma_vals, prob = TRUE, xlab = "Estimated gamma values", ylim=c(0,100),main = "Histogram of gamma values, true value = 0.05")
 lines(density(gamma_vals, adjust = 8), col = "blue",lwd = 2)
 abline(v=mean(gamma_vals),col="red", lwd = 2)
 abline(v=0.05, col="green", lwd = 2)
 legend(x=0.075,y=100,c("Density",
                        "Mean", "True value"),
        cex=1.2,col=c("blue","red","green"),lty = 1)
 
 hist(r_vals, prob = TRUE,xlab = "Estimated rn values",ylim = c(0,0.8),main = "Histogram of rn values, true value = 5")
 lines(density(r_vals, adjust = 8), col = "blue", lwd = 2)
 abline(v=mean(r_vals),col="red", lwd = 2)
 abline(v=5, col="green", lwd = 2)
 legend(x=8,y=0.8,c("Density",
                        "Mean", "True value"),
        cex=1.2,col=c("blue","red","green"),lty = 1)
 
 hist(ita_vals, prob = TRUE,xlab = "Estimated eta values",ylim = c(0,0.4),main = "Histogram of eta values, true value = 10")
 lines(density(ita_vals, adjust = 20), col = "blue", lwd = 2)
 abline(v=mean(ita_vals),col="red", lwd = 2)
 abline(v=10, col="green", lwd = 2)
 legend(x=15,y=0.4,c("Density",
                    "Mean", "True value"),
        cex=1.2,col=c("blue","red","green"),lty = 1)
 
 hist(dm_vals, prob = TRUE,xlab = "Estimated dm values",main = "Histogram of dm values, true value = 0.01")
 lines(density(dm_vals, adjust = 17), col = "blue", lwd = 2)
 abline(v=mean(dm_vals),col="red", lwd = 2)
 abline(v=0.01, col="green", lwd = 2)
 legend(x=0.04,y=60,c("Density",
                     "Mean", "True value"),
        cex=1.2,col=c("blue","red","green"),lty = 1)
 
 hist(al_vals, prob = TRUE,xlab = "Estimated alpha values",main = "Histogram of alpha values, true value = 0.1")
 lines(density(al_vals, adjust = 13), col = "blue", lwd = 2)
 abline(v=mean(al_vals),col="red", lwd = 2)
 abline(v=0.1, col="green", lwd = 2)
 legend(x=0.3,y=12,c("Density",
                      "Mean", "True value"),
        cex=1.2,col=c("blue","red","green"),lty = 1)
