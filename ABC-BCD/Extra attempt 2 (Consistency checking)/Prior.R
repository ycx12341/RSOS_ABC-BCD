####### Initial distributions of ABC-BCD scheme ############
####### Author: Yunchen Xiao ############################### 

####### ECM equation #######################################
set.seed(119)
dn<-runif(10000,0.000069,0.02)
gamma<-runif(10000,0.005,0.26)
eta<-runif(10000,7,18)
dm<-runif(10000,0.0001,0.033)
alpha<-runif(10000,0.07,0.18)
rn<-runif(10000,3.5,9)

paras_ecm_r1<-cbind(dn,gamma,eta,dm,alpha,rn)
write.table(paras_ecm_r1,"Round 1 parameters 10000 ecm.txt")
#############################################################

####### MDE equation ########################################
paras_ecm_r4 <- read.table("Round 4 parameters 10000 ecm.txt",sep="",header = TRUE)
set.seed(118) # A different set.seed is used to avoid obtaining the same initial parameters.
dn <- runif(10000,0.000069,0.02)
gamma <- runif(10000,0.005,0.26)
eta <- paras_ecm_r4[,3]
dm <- runif(10000,0.0001,0.033)
alpha <- runif(10000,0.07,0.18)
rn <- runif(10000,3.5,9)

paras_ecm_mde_r1 <- cbind(dn,gamma,eta,dm,alpha,rn)
write.table(paras_ecm_mde_r1, "Round 1 parameters 10000 ecm_mde.txt")
##############################################################


####### Tumour cells equations ###############################
paras_ecm_mde_r6 <- read.table("Round 6 parameters 10000 ecm_mde.txt",sep="",header = TRUE)
set.seed(117) # Another set.seed to avoid obtaining the same initial parameters. 
dn <- runif(10000,0.000069,0.02)
gamma <- runif(10000,0.005,0.26)
eta <- paras_ecm_mde_r6[,3]
dm <- paras_ecm_mde_r6[,4]
alpha <- paras_ecm_mde_r6[,5]
rn <- runif(10000, 3.5,9)

paras_all3_r1 <- cbind(dn,gamma,eta,dm,alpha,rn)
write.table(paras_all3_r1, "Round 1 parameters 10000 all 3.txt")
################################################################