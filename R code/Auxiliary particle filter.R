bcd_agg_p1<-read.table("B-C distance APF r1 p1.txt",sep="")
bcd_agg_p2<-read.table("B-C distance APF r1 p2.txt",sep="")
bcd_agg_p3<-read.table("B-C distance APF r1 p3.txt",sep="")
bcd_agg_p4<-read.table("B-C distance APF r1 p4.txt",sep="")
bcd_agg<-rbind(bcd_agg_p1,bcd_agg_p2,bcd_agg_p3,bcd_agg_p4)
ind_bcd_agg_nan<-which(bcd_agg[,2] == "NaN")
bcd_agg_nonan<-bcd_agg[-ind_bcd_agg_nan,]
paras_r1_agg<-read.table("Round 10 parameters 10000 all 3.txt",sep="") # Posterior parameters

post_joint<-rep(0,10000)

for (i in 1:length(post_joint)) {
  post_joint[i]<-prod(paras_r1_agg[i,])
} # Joint posterior parameters

ind_bcd_agg_nonan_ltmedian<-which(bcd_agg_nonan[,2] <= median(bcd_agg_nonan[,2]))
ind_bcd_agg_nonan_ltmedian_para<-bcd_agg_nonan[ind_bcd_agg_nonan_ltmedian,1]
bcd_agg_ltmedian<-bcd_agg[ind_bcd_agg_nonan_ltmedian_para,]
paras_r1_agg_ltmedian<-paras_r1_agg[ind_bcd_agg_nonan_ltmedian_para,]

info_mat<-cbind(bcd_agg_ltmedian,paras_r1_agg_ltmedian)

resamp_prob<-rep(0,length(info_mat[,1]))

for (i in 1:length(resamp_prob)) {
  if (info_mat[i,2] == max(info_mat[,2])) {
    resamp_prob[i]<-0
  } else if (info_mat[i,2] == min(info_mat[,2])) {
    resamp_prob[i]<-1
  } else {
    resamp_prob[i]<-1-((info_mat[i,2]-min(info_mat[,2]))/(max(info_mat[,2])-min(info_mat[,2])))
  }
}

info_mat_prob<-cbind(info_mat,resamp_prob)

resamp_ind<-sample(info_mat_prob[,1],size = 10000,replace = TRUE)
paras_samped_unperturbed<-paras_r1_agg[resamp_ind,]

paras_samped_perturbed<-matrix(0,nrow = nrow(paras_samped_unperturbed),ncol = ncol(paras_samped_unperturbed)) # with perturbation

for (i in 1:length(paras_samped_unperturbed[1,])) {
  for (j in 1:length(paras_samped_unperturbed[,1])){
    h = sqrt(1-0.05^2)
    paras_samped_perturbed[j,i]<-rnorm(1,h*paras_samped_unperturbed[j,i]+(1-h)*mean(paras_samped_unperturbed[,i]),
                                       0.05*sd(paras_samped_unperturbed[,i]))
  }
}  # Perturbation, see the manuscript for more details. 

# paras_samped_perturbed: Proposal parameters

prop_joint<-rep(0,length(paras_samped_perturbed[,1]))

for (i in 1:length(prop_joint)) {
  prop_joint[i]<-prod(paras_samped_perturbed[i,])
} # Joint proposal parameters.

paras_prior<-read.table("Round 1 parameters 10000 ecm.txt",sep="",header = TRUE)
prior_joint<-rep(0,length(paras_prior[,1]))

for (i in 1:length(post_joint)) {
  prior_joint[i]<-prod(paras_prior[i,])
} # Prior and joint prior parameters.

plot(density(prior_joint,adjust = 2),ylim=c(0,350000),main = "Prior, posterior and proposal joint densities (full scope)")
lines(density(post_joint,adjust = 2),col="red")
lines(density(prop_joint,adjust = 2),col="blue")
legend(x=0.0025,y=350000,c("Prior joint density","Posterior joint density","Proposal joint density"),
       cex=1.1,col=c("black","red","blue"),lty = 1)
# Plots of the joint densities. 

prior_joint_den<-density(prior_joint,adjust = 2,n=16384)
prop_joint_den<-density(prop_joint,adjust = 2)

prior_joint_den_x<-as.numeric(unlist(prior_joint_den["x"]))
prior_joint_den_y<-as.numeric(unlist(prior_joint_den["y"]))

prop_joint_den_x<-as.numeric(unlist(prop_joint_den["x"]))
prop_joint_den_y<-as.numeric(unlist(prop_joint_den["y"]))

prop_den_inprior<-rep(0,length(prior_joint))

for (i in 1:length(prop_den_inprior)) {
  prop_den_x_inprior<-which(abs(prop_joint[i]-prior_joint_den_x) == min(abs(prop_joint[i]-prior_joint_den_x)))
  prop_den_inprior[i]<-prior_joint_den_y[prop_den_x_inprior]
}

prop_den_inprop<-rep(0,length(prior_joint))

for (i in 1:length(prop_den_inprop)) {
  prop_den_x_inprop<-which(abs(prop_joint[i]-prop_joint_den_x) == min(abs(prop_joint[i]-prop_joint_den_x)))
  prop_den_inprop[i]<-prop_joint_den_y[prop_den_x_inprop]
} # Density of the proposal parameter vectors in prior and proposal density.

ind<-seq(1,10000, by = 1)
den_infomat<-cbind(ind,prop_joint,prop_den_inprior,prop_den_inprop)

weight_prior_prop<-prop_den_inprior/prop_den_inprop
weight_resamp_prob<-weight_prior_prop/max(weight_prior_prop)

den_infomat_wt_prob<-cbind(den_infomat,weight_prior_prop,weight_resamp_prob)

den_resamp_ind<-sample(den_infomat_wt_prob[,1],size = 10000, replace = TRUE,
                       den_infomat_wt_prob[,6])

paras_samped_unperturbed_flat<-paras_samped_perturbed[den_resamp_ind,]

paras_samped_perturbed_flat<-matrix(0,nrow = nrow(paras_samped_unperturbed_flat),ncol = ncol(paras_samped_unperturbed_flat)) # with perturbation

for (i in 1:length(paras_samped_unperturbed_flat[1,])) {
  for (j in 1:length(paras_samped_unperturbed_flat[,1])){
    h = sqrt(1-0.05^2)
    paras_samped_perturbed_flat[j,i]<-rnorm(1,h*paras_samped_unperturbed_flat[j,i]+(1-h)*mean(paras_samped_unperturbed_flat[,i]),
                                            0.05*sd(paras_samped_unperturbed_flat[,i]))
  }
} # Flattened posterior parameters, with perturbation.

write.table(paras_samped_perturbed_flat,"Flattened posterior parameters.txt")

paras_samped_perturbed_flat_joint<-rep(0,10000)

for (i in 1:10000) {
  paras_samped_perturbed_flat_joint[i]<-prod(paras_samped_perturbed_flat[i,])
} # Joint density of flattened posterior parameters. 

plot(density(prior_joint,adjust = 2),xlim = c(-0.00015,0.0007),ylim=c(0,500000),main = "Prior, posterior, proposal and flattened proposal joint densities (full scope)")
lines(density(post_joint,adjust = 2),col="red")
lines(density(prop_joint,adjust = 2),col="blue")
lines(density(paras_samped_perturbed_flat_joint,adjust=2),col="green")
legend(x=0.0004,y=500000,c("Prior joint density","Posterior joint density",
                           "Proposal joint density","Flattened proposal joint density"),
       cex=1.2,col=c("black","red","blue","green"),lty = 1)

plot(density(prior_joint,adjust = 2),main = "Prior, posterior, proposal and flattened proposal joint densities")
lines(density(post_joint,adjust = 2),col="red")
lines(density(prop_joint,adjust = 2),col="blue")
lines(density(paras_samped_perturbed_flat_joint,adjust=2),col="green")
legend(x=0.002,y=3500,c("Prior joint density","Posterior joint density",
                           "Proposal joint density","Flattened proposal joint density"),
       cex=1.2,col=c("black","red","blue","green"),lty = 1)
# Density plots. 