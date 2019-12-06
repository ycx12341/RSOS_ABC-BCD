clc
clear all
close all
 
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 24)
set(0,'defaultaxeslinewidth',1)
set(0,'defaultpatchlinewidth',1)
set(0,'defaultlinelinewidth',4)
set(0,'defaultTextInterpreter','latex')

paras = readtable('.txt'); 
% Read the table of parameter values.
mean_var_obs = readtable('mean_var_obs.txt'); 
mean_var_obs = table2array(mean_var_obs); 
% Read the table of observed summary statistics, change it into an array.

%%%%%%% Space %%%%%%%%%%%%%
l1=0; 
l2=1;
x11=linspace(l1,l2,300);
h=abs(x11(2)-x11(1)); 
% Discretize the dimensionless space, 0 to 1, into 300 steps.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Time %%%%%%%%%%%%%%%
T=10;
dt = 0.0001;
time = 0:dt:T;
tv=0; 
% 10 dimensionless time points altogether, dt set to be 0.0001 in the
% numerical solver.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Parameters %%%%%%%%%%%%%%%

dn = table2array(paras(:,2));


gamma = table2array(paras(:,3));


ita = table2array(paras(:,4));


dm = table2array(paras(:,5));


alpha = table2array(paras(:,6));

r = table2array(paras(:,7)); 
% For each parameter in the PDE system, read its values from the table and
% change the format to array. 

beta = 0;
eps = 0.01; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Initial condition %%%%%%%%%%%
n0 = repelem(0,length(x11));

for i = 1:length(x11)
    if x11(i)<=0.25
        n0(i)=exp(-(x11(i)^2)/eps);
    else 
        n0(i)=0;
    end 
end

n = n0;

f0 = 1-0.5*n0;
f=f0;

m0=0.5*n0; 
m = m0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Array of initial values %%%
inits = [n;f;m];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Initial plot %%%%%%%%%%%%%%%%%%%
figure
plot(x11,n,'b-',x11,f,'k-',x11,m,'g-')
axis([min(x11) max(x11) 0 max(f)+0.1])
axis square
xlabel('$x_1$')
title(['$System$ at $t=$',num2str(0)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Open file %%%%%%%%%%%%%%%%%
fileID = fopen('.txt','w'); 
% Open a .txt file, the Bhattacharya distance of each parameter vector in 
% relation to the reference parameter values will be written in it. 
fclose(fileID);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bcd_sum=[]; 

for i = 1:10000 
% For each round, the total number of parameter vectors investigated was 
% 10000. 
n = inits(1,:);
f = inits(2,:);
m = inits(3,:);

p=1;
%%% Results array %%%%%
res = [];
%%%%%%%%%% Numerical PDE solver %%%%%%%%%%%%%
while p*dt<=T
    f(2:length(x11)-1) = -ita(i)*dt*m(2:length(x11)-1).*f(2:length(x11)-1)+f(2:length(x11)-1);
    m(2:length(x11)-1) = dm(i)*(m(1:length(x11)-2)+m(3:length(x11))-2*m(2:length(x11)-1))*dt/(h^2)+alpha(i)*n(2:length(x11)-1)*dt-beta*m(2:length(x11)-1)*dt+m(2:length(x11)-1);
    n(2:length(x11)-1) = dn(i)*(n(1:length(x11)-2)+n(3:length(x11))-2*n(2:length(x11)-1))*dt/(h^2)...
                         -gamma(i)*(n(3:length(x11))-n(2:length(x11)-1)).*(f(3:length(x11))-f(2:length(x11)-1))*dt/(h^2)...
                         -gamma(i)*n(2:length(x11)-1).*(f(1:length(x11)-2)+f(3:length(x11))...
                         -2*f(2:length(x11)-1))*dt/(h^2)+n(2:length(x11)-1)+r(i)*(1-f(2:length(x11)-1)-n(2:length(x11)-1)).*n(2:length(x11)-1)*dt;                     
    
    %%%%%%%%%% Zero Neumann boundary condition %%%%%%%%%%
    n(1) = n(2);
    n(length(x11)) = n(length(x11)-1);
    
    f(1) = f(2);
    f(length(x11)) = f(length(x11)-1);
    
    m(1) = m(2);
    m(length(x11)) = m(length(x11)-1); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    if mod(p,10000)==0
        temp = round(round(time(p+1)*100)/100);
        res = [res;n;f;m];
    end 
     
    p=p+1; 
    % Store the results of different variables at different timepoints and  
    % different locations in the space. As described in the manuscript, for
    % each parameter vector, a 30*300 array should be created. 
     
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

res_arr = zeros(30,300);

res_arr(1:10,:) = res(1:3:28,:);
res_arr(11:20,:) = res(2:3:29,:);
res_arr(21:30,:) = res(3:3:30,:); 
% Rearrange the result array, row 1-10, TC, row 11-20, ECM, row 21-30, MDE.

mean_var = zeros(900,2); 
% Rearrange the summary statistics.

mean_var(1:300,1) = mean(res_arr(1:10,:)); 
% mean of the TC time series.
mean_var(1:300,2) = var(res_arr(1:10,:)); 
% variance of the TC time series.

mean_var(301:600,1) = mean(res_arr(11:20,:)); 
% mean of the ECM time series.
mean_var(301:600,2) = var(res_arr(11:20,:)); 
% variance of the ECM time series.

mean_var(601:900,1) = mean(res_arr(21:30,:)); 
% mean of the MDE time series.
mean_var(601:900,2) = var(res_arr(21:30,:)); 
% variance of the MDE time series

bcd_vec=[];
% Create an empty vector, store the Bhattacharya distance of each time
% series in relation to the corresponding observed one. 

for j = 1:900 
    bcd = 0.25*log(0.25*((mean_var(j,2)/mean_var_obs(j,2))+(mean_var_obs(j,2)/mean_var(j,2))+2))+0.25*(((mean_var(j,1)-mean_var_obs(j,1))^2)/(mean_var_obs(j,2)+mean_var(j,2)));
    bcd_vec = [bcd_vec,bcd];
end 

% Calculation of each time series' Bhattacharya distance in relation to the
% corresponding observed one. As described in the main body of the 
% manuscript, while inferring the ECM parameter (eta), only the time series 
% of ECM were considered (301:600), when we moved on to estimate the MDE 
% parameters (dm, alpha), we summed up the B-C distance among MDE time 
% series and ECM time series that have been evaluated before. (301:900), 
% finally, when the TC parameters were estimated (dn, gamma, rn), 
% we summed up the B-C distance among all the time series. (1:900)

inv_index = find(bcd_vec == Inf);
inv_term = length(inv_index);

bcd_vec_2 = bcd_vec;
bcd_vec_2(inv_index) = []; 
% Locate the invalid terms, exclude them from the B-C distance calculation. 

bcd_sum = [bcd_sum,sum(bcd_vec_2)]; 
% Sum up the Bhattacharya distance of the time series.

fileID = fopen('.txt','a');
fprintf(fileID,'%4d %5.4f\r\n',i,sum(bcd_vec_2));
fclose(fileID); 

% The Bhattacharya distance of each parameter vector in relation to the
% reference parameter values will be written in the .txt file, at the end
% there will be 10000 singular values stored in this .txt file. 

disp(i)

end 