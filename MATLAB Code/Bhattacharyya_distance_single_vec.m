%%%%%%%%%%%%%%% Environment setting %%%%%%%%%%%%%%%
clc
clear all
close all
 
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 24)
set(0,'defaultaxeslinewidth',1)
set(0,'defaultpatchlinewidth',1)
set(0,'defaultlinelinewidth',4)
set(0,'defaultTextInterpreter','latex')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mean_var_obs = readtable('mean_var_obs.txt');
% Read the table of observed summary statistics.
mean_var_obs = table2array(mean_var_obs);
% Change its format from table to an array

%%%%%%%% Space %%%%%%%%%%%%%%%%%%%%%%%
l1=0;
l2=1;
x11=linspace(l1,l2,300);
h=abs(x11(2)-x11(1));
% Discretization of dimensionless space, 0 to 1, 
% 300 steps.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Time %%%%%%%%%%%%%%%%%%%%%%%%%
T=10;
dt = 0.0001;
time = 0:dt:T;
% Dimensionless time points, 0 to 10,
% the stepsize of time in the finite difference
% scheme is defined to be 0.0001.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%
dn = 0.001550216;
gamma = 0.04672252;
ita = 9.796103;
dm = 0.01065483;
alpha = 0.0998378;
r = 4.84149;
% A single parameter vector, here we use the final 
% parameter estimation in the 1st attempt, mentioned in
% the manuscript.

beta = 0;
eps = 0.01; 
% These two parameters are fixed all the time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Initial condition %%%%%%%%%%%
n0 = repelem(0,length(x11));

for i = 1:length(x11)
    if x11(i)<=0.25
        n0(i)=exp(-(x11(i)^2)/eps);
    else 
        n0(i)=0;
    end 
end % Initial condition of tumour cells.

n = n0; % Initialize the tumour cells vector.

f0 = 1-0.5*n0; % Initial condition of ECM.
f=f0; % Initialize the ECM vector.

m0=0.5*n0; % Initial condition of MDE.
m = m0; % Initialize the MDE vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Array of initial values %%%%%%%
inits = [n;f;m];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Open file %%%%%%%%%%%%%%%%%
fileID = fopen('.txt','w');
fclose(fileID); 
% Open a .txt file to store the result of B-C distance 
% calculation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

res = []; 
% Empty vector used to store results later on.

p=1;

%%%%%%%% Numerical PDE solver %%%%%%%%%%%
while p*dt<=T
    f(2:length(x11)-1) = -ita*dt*m(2:length(x11)-1).*f(2:length(x11)-1)+f(2:length(x11)-1);
    m(2:length(x11)-1) = dm*(m(1:length(x11)-2)+m(3:length(x11))-2*m(2:length(x11)-1))*dt/(h^2)+alpha*n(2:length(x11)-1)*dt-beta*m(2:length(x11)-1)*dt+m(2:length(x11)-1);
    n(2:length(x11)-1) = dn*(n(1:length(x11)-2)+n(3:length(x11))-2*n(2:length(x11)-1))*dt/(h^2)...
                         -gamma*(n(3:length(x11))-n(2:length(x11)-1)).*(f(3:length(x11))-f(2:length(x11)-1))*dt/(h^2)...
                         -gamma*n(2:length(x11)-1).*(f(1:length(x11)-2)+f(3:length(x11))...
                         -2*f(2:length(x11)-1))*dt/(h^2)+n(2:length(x11)-1)+r*(1-f(2:length(x11)-1)-n(2:length(x11)-1)).*n(2:length(x11)-1)*dt;                     
    
    %%%%%%%%% Zero Neumann boundary condition %%%%%%%%%%%
    n(1) = n(2);
    n(length(x11)) = n(length(x11)-1);
    
    f(1) = f(2);
    f(length(x11)) = f(length(x11)-1);
    
    m(1) = m(2);
    m(length(x11)) = m(length(x11)-1); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    if mod(p,10000)==0
        temp = round(round(time(p+1)*100)/100);
        res = [res;n;f;m];
    end
     
    p=p+1;
     
end
% Store the results of different variables at different timepoints and 
% different locations in space returned from the PDE solver. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% Rearrange the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res_arr = zeros(30,300);
% A zero array of 30*300, used to store the rearranged results.

res_arr(1:10,:) = res(1:3:28,:); % Time series of tumour cells. 
res_arr(11:20,:) = res(2:3:29,:); % Time series of ECM.
res_arr(21:30,:) = res(3:3:30,:); % Time series of MDE.
% Rearrange the original results into time series of the format mentioned in
% the manuscript. 

mean_var = zeros(900,2); 
% An empty 900*2 array used to store summary statistics.

mean_var(1:300,1) = mean(res_arr(1:10,:));
% Means of tumour cells time series.
mean_var(1:300,2) = var(res_arr(1:10,:)); 
% Variances of tumour cells time series.

mean_var(301:600,1) = mean(res_arr(11:20,:)); 
% Means of ECM time series.
mean_var(301:600,2) = var(res_arr(11:20,:)); 
% Variances of ECM time series.

mean_var(601:900,1) = mean(res_arr(21:30,:)); 
% Means of MDE time series.
mean_var(601:900,2) = var(res_arr(21:30,:)); 
% Variances of MDE time series.

%%%%%%%%%%%%%%%%%% Bhattacharyya distance calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bcd_vec=[]; % An empty vector used to store the results produced later on.

for j = 1:900
    bcd = 0.25*log(0.25*((mean_var(j,2)/mean_var_obs(j,2))+(mean_var_obs(j,2)/mean_var(j,2))+2))+0.25*(((mean_var(j,1)-mean_var_obs(j,1))^2)/(mean_var_obs(j,2)+mean_var(j,2)));
    bcd_vec = [bcd_vec,bcd];
end 
% Calculate the Bhattacharya distance of each time series to its corresponding observed one,
% append the results into the empty vector created above.

inv_index = find(bcd_vec == Inf);
inv_term = length(inv_index);
% Locate the invalid terms in the vector. 
bcd_vec_2 = bcd_vec;
bcd_vec_2(inv_index) = [];
% Exclude the invalid terms.

fileID = fopen('.txt','a');
fprintf(fileID,'%4d %5.4f\r\n',i,sum(bcd_vec_2));
% The Bhattacharya distance of the parameter vector in relation to the
% reference parameter values will be written in the .txt file.
fclose(fileID);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
