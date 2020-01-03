%%%%%%%%%%%%%%%% Environment setting %%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all
 
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 24)
set(0,'defaultaxeslinewidth',1)
set(0,'defaultpatchlinewidth',1)
set(0,'defaultlinelinewidth',4)
set(0,'defaultTextInterpreter','latex')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% ECM: Density plot of eta while inferring ECM profile %%%%%%%%%%%%%
ecm_r1 = readtable("Round 1 parameters 10000 ecm.txt");
ecm_r2 = readtable("Round 2 parameters 10000 ecm.txt");
ecm_r3 = readtable("Round 3 parameters 10000 ecm.txt");
ecm_post = readtable("Round 4 parameters 10000 ecm.txt");

eta_r1 = table2array(ecm_r1(:,4));
eta_r2 = table2array(ecm_r2(:,4));
eta_r3 = table2array(ecm_r3(:,4));
eta_post = table2array(ecm_post(:,4));

[g,xii] = ksdensity(eta_r2,'Bandwidth',0.3014);
[h,xiii] = ksdensity(eta_r3,'Bandwidth',0.3014);
[j,xiv] = ksdensity(eta_post,'Bandwidth',0.3014);

figure
yline(1/11,'b','Linewidth',3.5);
xlim([7 18]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'k-')
plot(mean(table2array(ecm_r1(:,4))),0,'xb','markersize',20)
plot(mean(table2array(ecm_r2(:,4))),0,'xg','markersize',20)
plot(mean(table2array(ecm_r3(:,4))),0,'xr','markersize',20)
plot(mean(table2array(ecm_post(:,4))),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density(ECM)','Post-round 1 density(ECM)','Post-round 2 density(ECM)','Post-round 3 density(ECM)'},'Location','northeast','Orientation','vertical','Fontsize',12);
xlabel('$\eta$ values')
ylabel('Probability density')
% Read the parameter file, kernel smooth the probability densities and plot them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% MDE: Density plots of dm and alpha while inferring MDE profile %%%%%%%%%%%%%%%
ecm_mde_r1 = readtable("Round 1 parameters 10000 ecm_mde.txt");
ecm_mde_r2 = readtable("Round 2 parameters 10000 ecm_mde.txt");
ecm_mde_r3 = readtable("Round 3 parameters 10000 ecm_mde.txt");
ecm_mde_r4 = readtable("Round 4 parameters 10000 ecm_mde.txt");
ecm_mde_r5 = readtable("Round 5 parameters 10000 ecm_mde.txt");
ecm_mde_post = readtable("Round 6 parameters 10000 ecm_mde.txt");

dm_r1 = table2array(ecm_mde_r1(:,5));
dm_r2 = table2array(ecm_mde_r2(:,5));
dm_r3 = table2array(ecm_mde_r3(:,5));
dm_r4 = table2array(ecm_mde_r4(:,5));
dm_r5 = table2array(ecm_mde_r5(:,5));
dm_post = table2array(ecm_mde_post(:,5));

alpha_r1 = table2array(ecm_mde_r1(:,6));
alpha_r2 = table2array(ecm_mde_r2(:,6));
alpha_r3 = table2array(ecm_mde_r3(:,6));
alpha_r4 = table2array(ecm_mde_r4(:,6));
alpha_r5 = table2array(ecm_mde_r5(:,6));
alpha_post = table2array(ecm_mde_post(:,6));

[s,xxii] = ksdensity(table2array(ecm_mde_r2(:,5)),'Bandwidth',0.000902);
[d,xxiii] = ksdensity(table2array(ecm_mde_r3(:,5)),'Bandwidth',0.000902);
[z,xxiv] = ksdensity(table2array(ecm_mde_r4(:,5)),'Bandwidth',0.000902);
[q,xxv] = ksdensity(table2array(ecm_mde_r5(:,5)),'Bandwidth',0.000902);
[w,xxvi] = ksdensity(table2array(ecm_mde_post(:,5)),'Bandwidth',0.000902);

figure
yline(1/0.0329,'b-','Linewidth',3.5);
xlim([0.0001 0.033])
hold on;
plot(xxii,s,'g-',xxiii,d,'r-',xxiv,z,'k-',xxv,q,'y-',xxvi,w,'c-')
plot(0.01,0,'x','Color',[.69 .56 .74],'markersize',20)
plot(mean(table2array(ecm_mde_r1(:,5))),0,'xb','markersize',20)
plot(mean(table2array(ecm_mde_r2(:,5))),0,'xg','markersize',20)
plot(mean(table2array(ecm_mde_r3(:,5))),0,'xr','markersize',20)
plot(mean(table2array(ecm_mde_r4(:,5))),0,'xk','markersize',20)
plot(mean(table2array(ecm_mde_r5(:,5))),0,'xy','markersize',20)
plot(mean(table2array(ecm_mde_post(:,5))),0,'xc','markersize',20)
hold off;
legend({'Initial density(ECM&MDE)','Post-round 1 density(ECM&MDE)','Post-round 2 density(ECM&MDE)','Post-round 3 density(ECM&MDE)','Post-round 4 density(ECM&MDE)','Post-round 5 density(ECM&MDE)'},'Location','northeast','Orientation','vertical','Fontsize',9)
xlabel('$d_m$ values')
ylabel('Probability density')

[s,xxii] = ksdensity(table2array(ecm_mde_r2(:,6)),'Bandwidth',0.003034);
[d,xxiii] = ksdensity(table2array(ecm_mde_r3(:,6)),'Bandwidth',0.003034);
[z,xxiv] = ksdensity(table2array(ecm_mde_r4(:,6)),'Bandwidth',0.003034);
[q,xxv] = ksdensity(table2array(ecm_mde_r5(:,6)),'Bandwidth',0.003034);
[w,xxvi] = ksdensity(table2array(ecm_mde_post(:,6)),'Bandwidth',0.003034);

figure
yline(1/0.11,'b-','Linewidth',3.5);
xlim([0.07 0.18])
hold on;
plot(xxii,s,'g-',xxiii,d,'r-',xxiv,z,'k-',xxv,q,'y-',xxvi,w,'c-')
plot(mean(table2array(ecm_mde_r1(:,6))),0,'xb','markersize',20)
plot(mean(table2array(ecm_mde_r2(:,6))),0,'xg','markersize',20)
plot(mean(table2array(ecm_mde_r3(:,6))),0,'xr','markersize',20)
plot(mean(table2array(ecm_mde_r4(:,6))),0,'xk','markersize',20)
plot(mean(table2array(ecm_mde_r5(:,6))),0,'xy','markersize',20)
plot(mean(table2array(ecm_mde_post(:,6))),0,'xc','markersize',20)
hold off;
legend({'Initial density(ECM&MDE)','Post-round 1 density(ECM&MDE)','Post-round 2 density(ECM&MDE)','Post-round 3 density(ECM&MDE)','Post-round 4 density(ECM&MDE)','Post-round 5 density(ECM&MDE)'},'Location','northeast','Orientation','vertical','Fontsize',10)
xlabel('$\alpha$ values')
ylabel('Probability density')
% Note that due to the different scales of parameter values, different bandwidths were chosen to kernel smooth the probability
% densities for the best visual effect.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% Tumour cells: Density plots of dn,gamma,rn while inferring the tumour cells profile %%%%%%%%%%%%%%%%%%%%
all3_r1 = readtable("Round 1 parameters 10000 all 3.txt");
all3_r2 = readtable("Round 2 parameters 10000 all 3.txt");
all3_r3 = readtable("Round 3 parameters 10000 all 3.txt");
all3_r4 = readtable("Round 4 parameters 10000 all 3.txt");
all3_r5 = readtable("Round 5 parameters 10000 all 3.txt");
all3_r6 = readtable("Round 6 parameters 10000 all 3.txt");
all3_r7 = readtable("Round 7 parameters 10000 all 3.txt");
all3_r8 = readtable("Round 8 parameters 10000 all 3.txt");
all3_r9 = readtable("Round 9 parameters 10000 all 3.txt");
all3_post = readtable("Round 10 parameters 10000 all 3.txt");

[s,xxxii] = ksdensity(table2array(all3_r2(:,2)),'Bandwidth',0.0001641);
[d,xxxiii] = ksdensity(table2array(all3_r3(:,2)),'Bandwidth',0.0001641);
[z,xxxiv] = ksdensity(table2array(all3_r4(:,2)),'Bandwidth',0.0001641);
[q,xxxv] = ksdensity(table2array(all3_r5(:,2)),'Bandwidth',0.0001641);
[w,xxxvi] = ksdensity(table2array(all3_r6(:,2)),'Bandwidth',0.0001641);
[e,xxxvii] = ksdensity(table2array(all3_r7(:,2)),'Bandwidth',0.0001641);
[r,xxxviii] = ksdensity(table2array(all3_r8(:,2)),'Bandwidth',0.0001641);
[t,xxxix] = ksdensity(table2array(all3_r9(:,2)),'Bandwidth',0.0001641);
[y,xxxx] = ksdensity(table2array(all3_post(:,2)),'Bandwidth',0.0001641);

figure
yline(1/0.019931,'b-','Linewidth',3.5);
xlim([0.000069 0.02])
hold on;
plot(xxxii,s,'g-',xxxiii,d,'r-',xxxiv,z,'k-',xxxv,q,'y-',xxxvi,w,'c-')
plot(xxxvii,e,'Color',[.69 .56 .74])
plot(xxxviii,r,'Color',[.98 .69 .12])
plot(xxxix,t,'Color',[.13 .27 .48])
plot(xxxx,y,'Color',[.45 .62 .54])
plot(mean(table2array(all3_r1(:,2))),0,'xb','markersize',20)
plot(mean(table2array(all3_r2(:,2))),0,'xg','markersize',20)
plot(mean(table2array(all3_r3(:,2))),0,'xr','markersize',20)
plot(mean(table2array(all3_r4(:,2))),0,'xk','markersize',20)
plot(mean(table2array(all3_r5(:,2))),0,'xy','markersize',20)
plot(mean(table2array(all3_r6(:,2))),0,'xc','markersize',20)
plot(mean(table2array(all3_r7(:,2))),0,'x','Color',[.69 .56 .74],'markersize',20)
plot(mean(table2array(all3_r8(:,2))),0,'x','Color',[.98 .69 .12],'markersize',20)
plot(mean(table2array(all3_r9(:,2))),0,'x','Color',[.13 .27 .48],'markersize',20)
plot(mean(table2array(all3_post(:,2))),0,'x','Color',[.45 .62 .54],'markersize',20)
legend({'Initial density(All 3)','Post-round 1 density(All 3)','Post-round 2 density(All 3)','Post-round 3 density(All 3)','Post-round 4 density(All 3)','Post-round 5 density(All 3)','Post-round 6 density(All 3)','Post-round 7 density(All 3)','Post-round 8 density(All 3)','Final density'},'Location','northeast','Orientation','vertical','Fontsize',10)
hold off;
xlabel('$d_n$ values')
ylabel('Probability density')

[s,xxxii] = ksdensity(table2array(all3_r2(:,3)),'Bandwidth',0.00697);
[d,xxxiii] = ksdensity(table2array(all3_r3(:,3)),'Bandwidth',0.00697);
[z,xxxiv] = ksdensity(table2array(all3_r4(:,3)),'Bandwidth',0.00697);
[q,xxxv] = ksdensity(table2array(all3_r5(:,3)),'Bandwidth',0.00697);
[w,xxxvi] = ksdensity(table2array(all3_r6(:,3)),'Bandwidth',0.00697);
[e,xxxvii] = ksdensity(table2array(all3_r7(:,3)),'Bandwidth',0.00697);
[r,xxxviii] = ksdensity(table2array(all3_r8(:,3)),'Bandwidth',0.00697);
[t,xxxix] = ksdensity(table2array(all3_r9(:,3)),'Bandwidth',0.00697);
[y,xxxx] = ksdensity(table2array(all3_post(:,3)),'Bandwidth',0.00697);

figure
yline(1/0.255,'b-','Linewidth',3.5);
xlim([0.005 0.26]);
hold on;
plot(xxxii,s,'g-',xxxiii,d,'r-',xxxiv,z,'k-',xxxv,q,'y-',xxxvi,w,'c-')
plot(xxxvii,e,'Color',[.69 .56 .74])
plot(xxxviii,r,'Color',[.98 .69 .12])
plot(xxxix,t,'Color',[.13 .27 .48])
plot(xxxx,y,'Color',[.45 .62 .54])
plot(mean(table2array(all3_r1(:,3))),0,'xb','markersize',20)
plot(mean(table2array(all3_r2(:,3))),0,'xg','markersize',20)
plot(mean(table2array(all3_r3(:,3))),0,'xr','markersize',20)
plot(mean(table2array(all3_r4(:,3))),0,'xk','markersize',20)
plot(mean(table2array(all3_r5(:,3))),0,'xy','markersize',20)
plot(mean(table2array(all3_r6(:,3))),0,'xc','markersize',20)
plot(mean(table2array(all3_r7(:,3))),0,'x','Color',[.69 .56 .74],'markersize',20)
plot(mean(table2array(all3_r8(:,3))),0,'x','Color',[.98 .69 .12],'markersize',20)
plot(mean(table2array(all3_r9(:,3))),0,'x','Color',[.13 .27 .48],'markersize',20)
plot(mean(table2array(all3_post(:,3))),0,'x','Color',[.45 .62 .54],'markersize',20)
legend({'Initial density(All 3)','Post-round 1 density(All 3)','Post-round 2 density(All 3)','Post-round 3 density(All 3)','Post-round 4 density(All 3)','Post-round 5 density(All 3)','Post-round 6 density(All 3)','Post-round 7 density(All 3)','Post-round 8 density(All 3)','Final density'},'Location','northeast','Orientation','vertical','Fontsize',10)
hold off;
xlabel('$\gamma$ values')
ylabel('Probability density')

[s,xxxii] = ksdensity(table2array(all3_r2(:,7)),'Bandwidth',0.282);
[d,xxxiii] = ksdensity(table2array(all3_r3(:,7)),'Bandwidth',0.282);
[z,xxxiv] = ksdensity(table2array(all3_r4(:,7)),'Bandwidth',0.282);
[q,xxxv] = ksdensity(table2array(all3_r5(:,7)),'Bandwidth',0.282);
[w,xxxvi] = ksdensity(table2array(all3_r6(:,7)),'Bandwidth',0.282);
[e,xxxvii] = ksdensity(table2array(all3_r7(:,7)),'Bandwidth',0.282);
[r,xxxviii] = ksdensity(table2array(all3_r8(:,7)),'Bandwidth',0.282);
[t,xxxix] = ksdensity(table2array(all3_r9(:,7)),'Bandwidth',0.282);
[y,xxxx] = ksdensity(table2array(all3_post(:,7)),'Bandwidth',0.282);

figure
yline(1/5.5,'b-','Linewidth',3.5);
xlim([3.5 9]);
hold on;
plot(xxxii,s,'g-',xxxiii,d,'r-',xxxiv,z,'k-',xxxv,q,'y-',xxxvi,w,'c-')
plot(xxxvii,e,'Color',[.69 .56 .74])
plot(xxxviii,r,'Color',[.98 .69 .12])
plot(xxxix,t,'Color',[.13 .27 .48])
plot(xxxx,y,'Color',[.45 .62 .54])
plot(mean(table2array(all3_r1(:,7))),0,'xb','markersize',20)
plot(mean(table2array(all3_r2(:,7))),0,'xg','markersize',20)
plot(mean(table2array(all3_r3(:,7))),0,'xr','markersize',20)
plot(mean(table2array(all3_r4(:,7))),0,'xk','markersize',20)
plot(mean(table2array(all3_r5(:,7))),0,'xy','markersize',20)
plot(mean(table2array(all3_r6(:,7))),0,'xc','markersize',20)
plot(mean(table2array(all3_r7(:,7))),0,'x','Color',[.69 .56 .74],'markersize',20)
plot(mean(table2array(all3_r8(:,7))),0,'x','Color',[.98 .69 .12],'markersize',20)
plot(mean(table2array(all3_r9(:,7))),0,'x','Color',[.13 .27 .48],'markersize',20)
plot(mean(table2array(all3_post(:,7))),0,'x','Color',[.45 .62 .54],'markersize',20)
legend({'Initial density(All 3)','Post-round 1 density(All 3)','Post-round 2 density(All 3)','Post-round 3 density(All 3)','Post-round 4 density(All 3)','Post-round 5 density(All 3)','Post-round 6 density(All 3)','Post-round 7 density(All 3)','Post-round 8 density(All 3)','Final density'},'Location','northeast','Orientation','vertical','Fontsize',10)
hold off;
xlabel('$r_n$ values')
ylabel('Probability density')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% Initial and final densities only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[k,l] = ksdensity(table2array(all3_post(:,4)),'Bandwidth',0.3014);

figure
yline(1/11,'b','Linewidth',3.5);
xlim([7 18]);
hold on;
plot(l,k,'k-')
plot(10,0,'xr','markersize',20)
plot(mean(table2array(ecm_r1(:,4))),0,'xb','markersize',20)
plot(mean(table2array(all3_post(:,4))),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density','Final density','True value'},'Location','northeast','Orientation','vertical','Fontsize',15);
xlabel('$\eta$ values')
ylabel('Probability density')

[k,l] = ksdensity(table2array(all3_post(:,5)),'Bandwidth',0.000902);

figure
yline(1/0.0329,'b','Linewidth',3.5);
xlim([0.0001 0.033]);
hold on;
plot(l,k,'k-')
plot(0.01,0,'xr','markersize',20)
plot(mean(table2array(ecm_mde_r1(:,5))),0,'xb','markersize',20)
plot(mean(table2array(all3_post(:,5))),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density','Final density','True value'},'Location','northeast','Orientation','vertical','Fontsize',15);
xlabel('$d_m$ values')
ylabel('Probability density')

[k,l] = ksdensity(table2array(all3_post(:,6)),'Bandwidth',0.003034);

figure
yline(1/0.11,'b','Linewidth',3.5);
xlim([0.07 0.18]);
hold on;
plot(l,k,'k-')
plot(0.1,0,'xr','markersize',20)
plot(mean(table2array(ecm_mde_r1(:,6))),0,'xb','markersize',20)
plot(mean(table2array(all3_post(:,6))),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density','Final density','True value'},'Location','northeast','Orientation','vertical','Fontsize',15);
xlabel('$\alpha$ values')
ylabel('Probability density')

[k,l] = ksdensity(table2array(all3_post(:,2)),'Bandwidth',0.0001641);

figure
yline(1/0.019931,'b','Linewidth',3.5);
xlim([0.000069 0.02]);
hold on;
plot(l,k,'k-')
plot(0.001,0,'xr','markersize',20)
plot(mean(table2array(all3_r1(:,2))),0,'xb','markersize',20)
plot(mean(table2array(all3_post(:,2))),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density','Final density','True value'},'Location','northeast','Orientation','vertical','Fontsize',15);
xlabel('$d_n$ values')
ylabel('Probability density')

[k,l] = ksdensity(table2array(all3_post(:,3)),'Bandwidth',0.00697);

figure
yline(1/0.255,'b','Linewidth',3.5);
xlim([0.005 0.26]);
hold on;
plot(l,k,'k-')
plot(0.05,0,'xr','markersize',20)
plot(mean(table2array(all3_r1(:,3))),0,'xb','markersize',20)
plot(mean(table2array(all3_post(:,3))),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density','Final density','True value'},'Location','northeast','Orientation','vertical','Fontsize',15);
xlabel('$\gamma$ values')
ylabel('Probability density')

[k,l] = ksdensity(table2array(all3_post(:,7)),'Bandwidth',0.282);

figure
yline(1/5.5,'b','Linewidth',3.5);
xlim([3.5 9]);
hold on;
plot(l,k,'k-')
plot(5,0,'xr','markersize',20)
plot(mean(table2array(all3_r1(:,7))),0,'xb','markersize',20)
plot(mean(table2array(all3_post(:,7))),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density','Final density','True value'},'Location','northeast','Orientation','vertical','Fontsize',15);
xlabel('$r_n$ values')
ylabel('Probability density')
% Initial and final densities of the parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% ECM: Densities of eta while inferring the MDE profile %%%%%%%%%%%%%%%%%%%%%%%%%%%%
[k,xxi] = ksdensity(table2array(ecm_mde_r1(:,4)),'Bandwidth',0.3014);
[s,xxii] = ksdensity(table2array(ecm_mde_r2(:,4)),'Bandwidth',0.3014);
[d,xxiii] = ksdensity(table2array(ecm_mde_r3(:,4)),'Bandwidth',0.3014);
[z,xxiv] = ksdensity(table2array(ecm_mde_r4(:,4)),'Bandwidth',0.3014);
[q,xxv] = ksdensity(table2array(ecm_mde_r5(:,4)),'Bandwidth',0.3014);
[w,xxvi] = ksdensity(table2array(ecm_mde_post(:,4)),'Bandwidth',0.3014);

figure
plot(xxi,k,'b-')
xlim([7 18])
hold on;
plot(xxii,s,'g-',xxiii,d,'r-',xxiv,z,'k-',xxv,q,'y-',xxvi,w,'c-')
plot(mean(table2array(ecm_mde_r1(:,4))),0,'xb','markersize',20)
plot(mean(table2array(ecm_mde_r2(:,4))),0,'xg','markersize',20)
plot(mean(table2array(ecm_mde_r3(:,4))),0,'xr','markersize',20)
plot(mean(table2array(ecm_mde_r4(:,4))),0,'xk','markersize',20)
plot(mean(table2array(ecm_mde_r5(:,4))),0,'xy','markersize',20)
plot(mean(table2array(ecm_mde_post(:,4))),0,'xc','markersize',20)
hold off;
legend({'Round 1 density(ECM&MDE)','Post-round 1 density(ECM&MDE)','Post-round 2 density(ECM&MDE)','Post-round 3 density(ECM&MDE)','Post-round 4 density(ECM&MDE)','Post-round 5 density(ECM&MDE)'},'Location','northeast','Orientation','vertical','Fontsize',10)
xlabel('$\eta$ values')
ylabel('Probability density')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% ECM: Densities of eta while inferring the tumour cells profile. %%%%%%%%%%%%%%%%%%%%%%
[k,xxxi] = ksdensity(table2array(all3_r1(:,4)),'Bandwidth',0.3014);
[s,xxxii] = ksdensity(table2array(all3_r2(:,4)),'Bandwidth',0.3014);
[d,xxxiii] = ksdensity(table2array(all3_r3(:,4)),'Bandwidth',0.3014);
[z,xxxiv] = ksdensity(table2array(all3_r4(:,4)),'Bandwidth',0.3014);
[q,xxxv] = ksdensity(table2array(all3_r5(:,4)),'Bandwidth',0.3014);
[w,xxxvi] = ksdensity(table2array(all3_r6(:,4)),'Bandwidth',0.3014);
[e,xxxvii] = ksdensity(table2array(all3_r7(:,4)),'Bandwidth',0.3014);
[r,xxxviii] = ksdensity(table2array(all3_r8(:,4)),'Bandwidth',0.3014);
[t,xxxix] = ksdensity(table2array(all3_r9(:,4)),'Bandwidth',0.3014);
[y,xxxx] = ksdensity(table2array(all3_post(:,4)),'Bandwidth',0.3014);

figure
plot(xxxi,k,'b-')
xlim([7 18]);
hold on;
plot(xxxii,s,'g-',xxxiii,d,'r-',xxxiv,z,'k-',xxxv,q,'y-',xxxvi,w,'c-')
plot(xxxvii,e,'Color',[.69 .56 .74])
plot(xxxviii,r,'Color',[.98 .69 .12])
plot(xxxix,t,'Color',[.13 .27 .48])
plot(xxxx,y,'Color',[.45 .62 .54])
plot(mean(table2array(all3_r1(:,4))),0,'xb','markersize',20)
plot(mean(table2array(all3_r2(:,4))),0,'xg','markersize',20)
plot(mean(table2array(all3_r3(:,4))),0,'xr','markersize',20)
plot(mean(table2array(all3_r4(:,4))),0,'xk','markersize',20)
plot(mean(table2array(all3_r5(:,4))),0,'xy','markersize',20)
plot(mean(table2array(all3_r6(:,4))),0,'xc','markersize',20)
plot(mean(table2array(all3_r7(:,4))),0,'x','Color',[.69 .56 .74],'markersize',20)
plot(mean(table2array(all3_r8(:,4))),0,'x','Color',[.98 .69 .12],'markersize',20)
plot(mean(table2array(all3_r9(:,4))),0,'x','Color',[.13 .27 .48],'markersize',20)
plot(mean(table2array(all3_post(:,4))),0,'x','Color',[.45 .62 .54],'markersize',20)
legend({'Round 1 density(All 3)','Post-round 1 density(All 3)','Post-round 2 density(All 3)','Post-round 3 density(All 3)','Post-round 4 density(All 3)','Post-round 5 density(All 3)','Post-round 6 density(All 3)','Post-round 7 density(All 3)','Post-round 8 density(All 3)','Final density'},'Location','northeast','Orientation','vertical','Fontsize',10)
hold off;
xlabel('$\eta$ values')
ylabel('Probability density')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% MDE: Densities of dm, alpha while inferring tumour cells profile %%%%%%%%%%%%%%%%%%%%%%%
[k,xxxi] = ksdensity(table2array(all3_r1(:,5)),'Bandwidth',0.000902);
[s,xxxii] = ksdensity(table2array(all3_r2(:,5)),'Bandwidth',0.000902);
[d,xxxiii] = ksdensity(table2array(all3_r3(:,5)),'Bandwidth',0.000902);
[z,xxxiv] = ksdensity(table2array(all3_r4(:,5)),'Bandwidth',0.000902);
[q,xxxv] = ksdensity(table2array(all3_r5(:,5)),'Bandwidth',0.000902);
[w,xxxvi] = ksdensity(table2array(all3_r6(:,5)),'Bandwidth',0.000902);
[e,xxxvii] = ksdensity(table2array(all3_r7(:,5)),'Bandwidth',0.000902);
[r,xxxviii] = ksdensity(table2array(all3_r8(:,5)),'Bandwidth',0.000902);
[t,xxxix] = ksdensity(table2array(all3_r9(:,5)),'Bandwidth',0.000902);
[y,xxxx] = ksdensity(table2array(all3_post(:,5)),'Bandwidth',0.000902);

figure
plot(xxxi,k,'b-')
xlim([0.0001 0.033]);
hold on;
plot(xxxii,s,'g-',xxxiii,d,'r-',xxxiv,z,'k-',xxxv,q,'y-',xxxvi,w,'c-')
plot(xxxvii,e,'Color',[.69 .56 .74])
plot(xxxviii,r,'Color',[.98 .69 .12])
plot(xxxix,t,'Color',[.13 .27 .48])
plot(xxxx,y,'Color',[.45 .62 .54])
plot(mean(table2array(all3_r1(:,5))),0,'xb','markersize',20)
plot(mean(table2array(all3_r2(:,5))),0,'xg','markersize',20)
plot(mean(table2array(all3_r3(:,5))),0,'xr','markersize',20)
plot(mean(table2array(all3_r4(:,5))),0,'xk','markersize',20)
plot(mean(table2array(all3_r5(:,5))),0,'xy','markersize',20)
plot(mean(table2array(all3_r6(:,5))),0,'xc','markersize',20)
plot(mean(table2array(all3_r7(:,5))),0,'x','Color',[.69 .56 .74],'markersize',20)
plot(mean(table2array(all3_r8(:,5))),0,'x','Color',[.98 .69 .12],'markersize',20)
plot(mean(table2array(all3_r9(:,5))),0,'x','Color',[.13 .27 .48],'markersize',20)
plot(mean(table2array(all3_post(:,5))),0,'x','Color',[.45 .62 .54],'markersize',20)
legend({'Round 1 density(All 3)','Post-round 1 density(All 3)','Post-round 2 density(All 3)','Post-round 3 density(All 3)','Post-round 4 density(All 3)','Post-round 5 density(All 3)','Post-round 6 density(All 3)','Post-round 7 density(All 3)','Post-round 8 density(All 3)','Final density'},'Location','northeast','Orientation','vertical','Fontsize',10)
hold off;
xlabel('$d_m$ values')
ylabel('Probability density')

[k,xxxi] = ksdensity(table2array(all3_r1(:,6)),'Bandwidth',0.003034);
[s,xxxii] = ksdensity(table2array(all3_r2(:,6)),'Bandwidth',0.003034);
[d,xxxiii] = ksdensity(table2array(all3_r3(:,6)),'Bandwidth',0.003034);
[z,xxxiv] = ksdensity(table2array(all3_r4(:,6)),'Bandwidth',0.003034);
[q,xxxv] = ksdensity(table2array(all3_r5(:,6)),'Bandwidth',0.003034);
[w,xxxvi] = ksdensity(table2array(all3_r6(:,6)),'Bandwidth',0.003034);
[e,xxxvii] = ksdensity(table2array(all3_r7(:,6)),'Bandwidth',0.003034);
[r,xxxviii] = ksdensity(table2array(all3_r8(:,6)),'Bandwidth',0.003034);
[t,xxxix] = ksdensity(table2array(all3_r9(:,6)),'Bandwidth',0.003034);
[y,xxxx] = ksdensity(table2array(all3_post(:,6)),'Bandwidth',0.003034);

figure
plot(xxxi,k,'b-')
xlim([0.07 0.18]);
hold on;
plot(xxxii,s,'g-',xxxiii,d,'r-',xxxiv,z,'k-',xxxv,q,'y-',xxxvi,w,'c-')
plot(xxxvii,e,'Color',[.69 .56 .74])
plot(xxxviii,r,'Color',[.98 .69 .12])
plot(xxxix,t,'Color',[.13 .27 .48])
plot(xxxx,y,'Color',[.45 .62 .54])
plot(mean(table2array(all3_r1(:,6))),0,'xb','markersize',20)
plot(mean(table2array(all3_r2(:,6))),0,'xg','markersize',20)
plot(mean(table2array(all3_r3(:,6))),0,'xr','markersize',20)
plot(mean(table2array(all3_r4(:,6))),0,'xk','markersize',20)
plot(mean(table2array(all3_r5(:,6))),0,'xy','markersize',20)
plot(mean(table2array(all3_r6(:,6))),0,'xc','markersize',20)
plot(mean(table2array(all3_r7(:,6))),0,'x','Color',[.69 .56 .74],'markersize',20)
plot(mean(table2array(all3_r8(:,6))),0,'x','Color',[.98 .69 .12],'markersize',20)
plot(mean(table2array(all3_r9(:,6))),0,'x','Color',[.13 .27 .48],'markersize',20)
plot(mean(table2array(all3_post(:,6))),0,'x','Color',[.45 .62 .54],'markersize',20)
legend({'Round 1 density(All 3)','Post-round 1 density(All 3)','Post-round 2 density(All 3)','Post-round 3 density(All 3)','Post-round 4 density(All 3)','Post-round 5 density(All 3)','Post-round 6 density(All 3)','Post-round 7 density(All 3)','Post-round 8 density(All 3)','Final density'},'Location','northeast','Orientation','vertical','Fontsize',10)
hold off;
xlabel('$\alpha$ values')
ylabel('Probability density')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
