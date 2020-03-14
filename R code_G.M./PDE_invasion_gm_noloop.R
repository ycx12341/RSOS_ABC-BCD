library(tictoc)
library(mgcv)

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
                grad.lhs_m, grad.rhs_dm, grad.rhs_alpha) {
  
  disp1 <- sum((grad.lhs_f + paras[4] * grad.rhs_ita)^2)
  disp2 <- sum((grad.lhs_m - (paras[5] * grad.rhs_dm + paras[6] * grad.rhs_alpha))^2)
  disp3 <- sum((grad.lhs_n - (paras[1] * grad.rhs_dn - paras[2] * grad.rhs_gamma + paras[3] * grad.rhs_r))^2)
  res <- sum(disp1,disp2,disp3)
  
  #if (debug) cat(dn, res, "\n")
  return(res)
}

#scale = 0.02

res_n <- res[seq(1,151,by = 3),]
res_f <- res[seq(2,152,by = 3),]
res_m <- res[seq(3,153,by = 3),] # Unperturbed data

res_n_perturbed <- res_n + rnorm(51*50, 0, sd = 0.02*sd(res_n))
res_f_perturbed <- res_f + rnorm(51*50, 0, sd = 0.02*sd(res_f))
res_m_perturbed <- res_m + rnorm(51*50, 0, sd = 0.02*sd(res_m))

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

grad_rhs_f_ita <- matrix(0,49,48)

for (i in 1:49) {
  for (j in 1:48) {
    predict_f <- predict(spl2, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1]))
    predict_m <- predict(spl3, newdata = data.frame(t = tp[i + 1], x11 = x11[j + 1]))
    grad_rhs_f_ita[i, j] <- predict_f*predict_m
  }
}

 disp_f <- sum((grad_lhs_f_data-grad_rhs_f_ita*(-10))^2) 

# Then let's go for m. 

grad_lhs_m_data <- matrix(0,49,48)

for (i in 1:49) {
  for (j in 1:48) {
    predict_m_t1 <- predict(spl3, newdata = data.frame(t = tp[i + 1] + dt, x11 = x11[j + 1]))
    predict_m_t2 <- predict(spl3, newdata = data.frame(t = tp[i + 1] - dt, x11 = x11[j + 1]))
    grad_lhs_m_data[i,j] <- (predict_m_t1 - predict_m_t2) / (2 * dt)
  }
}

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

disp_m <- sum((grad_lhs_m_data-grad_rhs_m_dm*dm-grad_rhs_m_alpha*alpha)^2) # Sum of squares: 0.108

# At last, for n. 

grad_lhs_n_data <- matrix(0,49,48)

for (i in 1:49) {
  for (j in 1:48) {
    predict_n_t1 <- predict(spl, newdata = data.frame(t = tp[i + 1] + dt, x11 = x11[j + 1]))
    predict_n_t2 <- predict(spl, newdata = data.frame(t = tp[i + 1] - dt, x11 = x11[j + 1]))
    grad_lhs_n_data[i,j] <- (predict_n_t1 - predict_n_t2) / (2 * dt)
  }
}

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

grad_rhs_n <- dn*grad_rhs_n_dn-gamma*(grad_rhs_n_gamma)+r*grad_rhs_n_r
disp_n <- sum((grad_lhs_n_data-grad_rhs_n)^2)


start.values = c(0.02,0.1,10,20,0.02,0.2)

opt <- optim(start.values, obj, grad.lhs_n = grad_lhs_n_data, grad.rhs_dn = grad_rhs_n_dn, grad.rhs_gamma = grad_rhs_n_gamma,
      grad.rhs_r = grad_rhs_n_r, grad.lhs_f = grad_lhs_f_data, grad.rhs_ita = grad_rhs_f_ita,
      grad.lhs_m = grad_lhs_m_data, grad.rhs_dm = grad_rhs_m_dm, grad.rhs_alpha = grad_rhs_m_alpha, hessian = TRUE)
