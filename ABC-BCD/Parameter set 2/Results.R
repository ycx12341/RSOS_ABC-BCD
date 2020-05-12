######## Results ###############

# To reproduce the results here, set.seed(123) should be used within the abc_bcd algorithm. 
# Cluster formation

library(doParallel)
cl <- makeCluster(detectCores() - 1) 
registerDoParallel(cl)

############# Round 1 ECM ################
bcd_ecm_r1 <- bcd(paras_ecm_r1)
write.table(bcd_ecm_r1, "bcd_ecm_r1.txt") 
bcd_ecm_r1 <- read.table("bcd_ecm_r1.txt",sep="",header=TRUE)
ind_bcd_ecm_r1_nan<-which(is.na(bcd_ecm_r1[,2]) == TRUE)
bcd_ecm_r1_nonan<-bcd_ecm_r1[-ind_bcd_ecm_r1_nan,]
mean(bcd_ecm_r1_nonan[,2]) # Avg BCD: 4.338274

paras_ecm_r2 <- abc_bcd(bcd_ecm_r1,paras_ecm_r1)
write.table(paras_ecm_r2, "Round 2 parameters 10000 ecm.txt")
##########################################

############# Round 2 ECM ################
paras_ecm_r2 <- read.table("Round 2 parameters 10000 ecm.txt", sep="", header = TRUE)
bcd_ecm_r2 <- bcd(paras_ecm_r2)
write.table(bcd_ecm_r2, "bcd_ecm_r2.txt")
bcd_ecm_r2 <- read.table("bcd_ecm_r2.txt",sep="",header = TRUE)
ind_bcd_ecm_r2_nan<-which(is.na(bcd_ecm_r2[,2]) == TRUE)
bcd_ecm_r2_nonan<-bcd_ecm_r2[-ind_bcd_ecm_r2_nan,]
mean(bcd_ecm_r2_nonan[,2]) # Avg BCD: 1.825226

paras_ecm_r3 <- abc_bcd(bcd_ecm_r2, paras_ecm_r2)
write.table(paras_ecm_r3, "Round 3 parameters 10000 ecm.txt")

############# Round 3 ECM #################
paras_ecm_r3 <- read.table("Round 3 parameters 10000 ecm.txt", sep="", header = TRUE)
bcd_ecm_r3 <- bcd(paras_ecm_r3)
write.table(bcd_ecm_r3, "bcd_ecm_r3.txt")
bcd_ecm_r3 <- read.table("bcd_ecm_r3.txt",sep = "", header = TRUE)
ind_bcd_ecm_r3_nan<-which(is.na(bcd_ecm_r3[,2]) == TRUE)
bcd_ecm_r3_nonan<-bcd_ecm_r3[-ind_bcd_ecm_r3_nan,]
mean(bcd_ecm_r3_nonan[,2]) # Avg BCD: 0.4703119

paras_ecm_r4 <- abc_bcd(bcd_ecm_r3,paras_ecm_r3)
write.table(paras_ecm_r4, "Round 4 parameters 10000 ecm.txt")
############################################

########### Round 1 ECM MDE ####################
paras_ecm_mde_r1 <- read.table("Round 1 parameters 10000 ecm_mde.txt",sep="",header = TRUE)
bcd_ecm_mde_r1 <- bcd(paras_ecm_mde_r1)
write.table(bcd_ecm_mde_r1, "bcd_ecm_mde_r1.txt")
bcd_ecm_mde_r1 <- read.table("bcd_ecm_mde_r1.txt",sep="",header = TRUE)
ind_bcd_ecm_mde_r1_nan<-which(is.na(bcd_ecm_mde_r1[,2]) == TRUE)
bcd_ecm_mde_r1_nonan<-bcd_ecm_mde_r1[-ind_bcd_ecm_mde_r1_nan,]
mean(bcd_ecm_mde_r1_nonan[,2]) # Avg BCD: 6.411427

paras_ecm_mde_r2 <- abc_bcd(bcd_ecm_mde_r1, paras_ecm_mde_r1)
write.table(paras_ecm_mde_r2, "Round 2 parameters 10000 ecm_mde.txt")

########## Round 2 ECM MDE #####################
paras_ecm_mde_r2 <- read.table("Round 2 parameters 10000 ecm_mde.txt",sep="",header = TRUE)
bcd_ecm_mde_r2 <- bcd(paras_ecm_mde_r2)
write.table(bcd_ecm_mde_r2, "bcd_ecm_mde_r2.txt")
bcd_ecm_mde_r2 <- read.table("bcd_ecm_mde_r2.txt",sep="",header = TRUE)
ind_bcd_ecm_mde_r2_nan<-which(is.na(bcd_ecm_mde_r2[,2]) == TRUE)
bcd_ecm_mde_r2_nonan<-bcd_ecm_mde_r2[-ind_bcd_ecm_mde_r2_nan,]
mean(bcd_ecm_mde_r2_nonan[,2]) # Avg BCD: 3.247181

paras_ecm_mde_r3 <- abc_bcd(bcd_ecm_mde_r2, paras_ecm_mde_r2)
write.table(paras_ecm_mde_r3, "Round 3 parameters 10000 ecm_mde.txt")

######### Round 3 ECM MDE ######################
paras_ecm_mde_r3 <- read.table("Round 3 parameters 10000 ecm_mde.txt",sep="",header = TRUE)
bcd_ecm_mde_r3 <- bcd(paras_ecm_mde_r3)
write.table(bcd_ecm_mde_r3, "bcd_ecm_mde_r3.txt")
bcd_ecm_mde_r3 <- read.table("bcd_ecm_mde_r3.txt",sep="",header = TRUE)
ind_bcd_ecm_mde_r3_nan <- which(is.na(bcd_ecm_mde_r3[,2]) == TRUE)
bcd_ecm_mde_r3_nonan <- bcd_ecm_mde_r3[-ind_bcd_ecm_mde_r3_nan,]
mean(bcd_ecm_mde_r3_nonan[,2]) # Avg BCD: 1.309731

paras_ecm_mde_r4 <- abc_bcd(bcd_ecm_mde_r3, paras_ecm_mde_r3)
write.table(paras_ecm_mde_r4, "Round 4 parameters 10000 ecm_mde.txt")

######### Round 4 ECM MDE ######################
paras_ecm_mde_r4 <- read.table("Round 4 parameters 10000 ecm_mde.txt",sep="",header = TRUE)
bcd_ecm_mde_r4 <- bcd(paras_ecm_mde_r4)
write.table(bcd_ecm_mde_r4, "bcd_ecm_mde_r4.txt")
bcd_ecm_mde_r4 <- read.table("bcd_ecm_mde_r4.txt",sep="",header=TRUE)
ind_bcd_ecm_mde_r4_nan <- which(is.na(bcd_ecm_mde_r4[,2]) == TRUE)
bcd_ecm_mde_r4_nonan <- bcd_ecm_mde_r4[-ind_bcd_ecm_mde_r4_nan,]
mean(bcd_ecm_mde_r4_nonan[,2]) # Avg BCD: 0.4099569

paras_ecm_mde_r5 <- abc_bcd(bcd_ecm_mde_r4, paras_ecm_mde_r4)
write.table(paras_ecm_mde_r5, "Round 5 parameters 10000 ecm_mde.txt")

######### Round 5 ECM MDE #######################
paras_ecm_mde_r5 <- read.table("Round 5 parameters 10000 ecm_mde.txt",sep="",header = TRUE)
bcd_ecm_mde_r5 <- bcd(paras_ecm_mde_r5)
write.table(bcd_ecm_mde_r5, "bcd_ecm_mde_r5.txt")
bcd_ecm_mde_r5 <- read.table("bcd_ecm_mde_r5.txt",sep="",header=TRUE)
ind_bcd_ecm_mde_r5_nan <- which(is.na(bcd_ecm_mde_r5[,2]) == TRUE)
bcd_ecm_mde_r5_nonan <- bcd_ecm_mde_r5[-ind_bcd_ecm_mde_r5_nan,]
mean(bcd_ecm_mde_r5_nonan[,2]) # Avg BCD: 0.152304

paras_ecm_mde_r6 <- abc_bcd(bcd_ecm_mde_r5, paras_ecm_mde_r5)
write.table(paras_ecm_mde_r6, "Round 6 parameters 10000 ecm_mde.txt")

######### Round 1 ALL 3 ##########################
paras_all3_r1 <- read.table("Round 1 parameters 10000 all 3.txt",sep="",header = TRUE)
bcd_all3_r1 <- bcd(paras_all3_r1)
write.table(bcd_all3_r1, "bcd_all3_r1.txt")
bcd_all3_r1 <- read.table("bcd_all3_r1.txt",sep="",header = TRUE)
ind_bcd_all3_r1_nan <- which(is.na(bcd_all3_r1[,2]) == TRUE)
bcd_all3_r1_nonan <- bcd_all3_r1[-ind_bcd_all3_r1_nan,]
mean(bcd_all3_r1_nonan[,2]) # Avg BCD: 2.886287

paras_all3_r2 <- abc_bcd(bcd_all3_r1, paras_all3_r1)
write.table(paras_all3_r2, "Round 2 parameters 10000 all 3.txt")

######## Round 2 ALL 3 ###########################
paras_all3_r2 <- read.table("Round 2 parameters 10000 all 3.txt",sep="",header = TRUE)
bcd_all3_r2 <- bcd(paras_all3_r2)
write.table(bcd_all3_r2, "bcd_all3_r2.txt")
bcd_all3_r2 <- read.table("bcd_all3_r2.txt",sep="",header = TRUE)
ind_bcd_all3_r2_nan <- which(is.na(bcd_all3_r2[,2]) == TRUE)
bcd_all3_r2_nonan <- bcd_all3_r2[-ind_bcd_all3_r2_nan,]
mean(bcd_all3_r2_nonan[,2]) # Avg BCD: 1.232785

paras_all3_r3 <- abc_bcd(bcd_all3_r2,paras_all3_r2)
write.table(paras_all3_r3, "Round 3 parameters 10000 all 3.txt")

####### Round 3 ALL 3 ############################
paras_all3_r3 <- read.table("Round 3 parameters 10000 all 3.txt",sep="",header = TRUE)
bcd_all3_r3 <- bcd(paras_all3_r3)
write.table(bcd_all3_r3, "bcd_all3_r3.txt")
bcd_all3_r3 <- read.table("bcd_all3_r3.txt",sep="",header = TRUE)
ind_bcd_all3_r3_nan <- which(is.na(bcd_all3_r3[,2]) == TRUE)
bcd_all3_r3_nonan <- bcd_all3_r3[-ind_bcd_all3_r3_nan,]
mean(bcd_all3_r3_nonan[,2]) # Avg BCD: 0.357828

paras_all3_r4 <- abc_bcd(bcd_all3_r3, paras_all3_r3)
write.table(paras_all3_r4, "Round 4 parameters 10000 all 3.txt")

###### Round 4 ALL 3 #############################
paras_all3_r4 <- read.table("Round 4 parameters 10000 all 3.txt",sep = "", header = TRUE)
bcd_all3_r4 <- bcd(paras_all3_r4)
write.table(bcd_all3_r4, "bcd_all3_r4.txt")
bcd_all3_r4 <- read.table("bcd_all3_r4.txt",sep="",header = TRUE)
ind_bcd_all3_r4_nan <- which(is.na(bcd_all3_r4[,2]) == TRUE)
bcd_all3_r4_nonan <- bcd_all3_r4[-ind_bcd_all3_r4_nan,]
mean(bcd_all3_r4_nonan[,2]) # Avg BCD: 0.1289569

paras_all3_r5 <- abc_bcd(bcd_all3_r4, paras_all3_r4)
write.table(paras_all3_r5, "Round 5 parameters 10000 all 3.txt")

###### Round 5 ALL 3 #############################
paras_all3_r5 <- read.table("Round 5 parameters 10000 all 3.txt",sep="", header = TRUE)
bcd_all3_r5 <- bcd(paras_all3_r5)
write.table(bcd_all3_r5, "bcd_all3_r5.txt")
bcd_all3_r5 <- read.table("bcd_all3_r5.txt",sep="",header = TRUE)
ind_bcd_all3_r5_nan <- which(is.na(bcd_all3_r5[,2]) == TRUE)
bcd_all3_r5_nonan <- bcd_all3_r5[-ind_bcd_all3_r5_nan,]
mean(bcd_all3_r5_nonan[,2]) # Avg BCD: 0.07516977

paras_all3_r6 <- abc_bcd(bcd_all3_r5,paras_all3_r5)
write.table(paras_all3_r6, "Round 6 parameters 10000 all 3.txt")

###### Round 6 ALL 3 #############################
paras_all3_r6 <- read.table("Round 6 parameters 10000 all 3.txt",sep="", header = TRUE)
bcd_all3_r6 <- bcd(paras_all3_r6)
write.table(bcd_all3_r6, "bcd_all3_r6.txt")
bcd_all3_r6 <- read.table("bcd_all3_r6.txt",sep="",header = TRUE)
ind_bcd_all3_r6_nan <- which(is.na(bcd_all3_r6[,2]) == TRUE)
bcd_all3_r6_nonan <- bcd_all3_r6[-ind_bcd_all3_r6_nan,]
mean(bcd_all3_r6_nonan[,2]) # Avg BCD: 0.04936465

paras_all3_r7 <- abc_bcd(bcd_all3_r6, paras_all3_r6)
write.table(paras_all3_r7, "Round 7 parameters 10000 all 3.txt")

###### Round 7 ALL 3 #############################
paras_all3_r7 <- read.table("Round 7 parameters 10000 all 3.txt",sep="",header = TRUE)
bcd_all3_r7 <- bcd(paras_all3_r7)
write.table(bcd_all3_r7, "bcd_all3_r7.txt")
bcd_all3_r7 <- read.table("bcd_all3_r7.txt",sep="",header = TRUE)
# ind_bcd_all3_r7_nan <- which(is.na(bcd_all3_r7[,2]) == TRUE)
# bcd_all3_r7_nonan <- bcd_all3_r7[-ind_bcd_all3_r7_nan,]
mean(bcd_all3_r7[,2]) # Avg BCD: 0.03819095

paras_all3_r8 <- abc_bcd(bcd_all3_r7, paras_all3_r7)
write.table(paras_all3_r8, "Round 8 parameters 10000 all 3.txt")

###### Round 8 ALL 3 #############################
paras_all3_r8 <- read.table("Round 8 parameters 10000 all 3.txt", sep="", header = TRUE)
bcd_all3_r8 <- bcd(paras_all3_r8)
write.table(bcd_all3_r8, "bcd_all3_r8.txt")
bcd_all3_r8 <- read.table("bcd_all3_r8.txt",sep="",header = TRUE)
ind_bcd_all3_r8_nan <- which(is.na(bcd_all3_r8[,2]) == TRUE)
bcd_all3_r8_nonan <- bcd_all3_r8[-ind_bcd_all3_r8_nan,]
mean(bcd_all3_r8_nonan[,2]) # Avg BCD: 0.03059211

paras_all3_r9 <- abc_bcd(bcd_all3_r8, paras_all3_r8)
write.table(paras_all3_r9, "Round 9 parameters 10000 all 3.txt")

###### Round 9 ALL 3 #############################
paras_all3_r9 <- read.table("Round 9 parameters 10000 all 3.txt",sep="", header = TRUE)
bcd_all3_r9 <- bcd(paras_all3_r9)
write.table(bcd_all3_r9, "bcd_all3_r9.txt")
bcd_all3_r9 <- read.table("bcd_all3_r9.txt",sep="",header = TRUE)
# ind_bcd_all3_r9_nan <- which(is.na(bcd_all3_r9[,2]) == TRUE)
# bcd_all3_r9_nonan <- bcd_all3_r9[-ind_bcd_all3_r9_nan,]
mean(bcd_all3_r9[,2]) # Avg BCD: 0.02583689

paras_all3_r10 <- abc_bcd(bcd_all3_r9, paras_all3_r9)
write.table(paras_all3_r10, "Round 10 parameters 10000 all 3.txt")

##### Final estimations #########################

paras_all3_r10 <- read.table("Round 10 parameters 10000 all 3.txt",sep="",header = TRUE)

mean(paras_all3_r10[,1]) # dn = 0.01277436
mean(paras_all3_r10[,2]) # gamma = 0.1517392
mean(paras_all3_r10[,3]) # eta = 11.0013
mean(paras_all3_r10[,4]) # dm = 0.02542924
mean(paras_all3_r10[,5]) # alpha = 0.1099737
mean(paras_all3_r10[,6]) # rn = 6.147018
