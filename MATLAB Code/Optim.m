% Standard optimizer aims to find the set of parameter values that minimize the Bhattacharyya distance between 
% its simulated time series and the observed time series. 

% Used in: Xiao et al. "Calibrating models of cancer invasion and metastasis: parameter optimization using Approximate
% Bayesian Computation."

% Author: Yunchen Xiao

x0 = [0.0100345,0.1325,6.25,12.5,0.01655,0.125];
% Starting values of the optimizing function.
options = optimset('Display','iter');
% Display the number of iterations.
[x,fval,exitflag,output] = fminsearch(@objectivefcn1,x0,options);
% Call the optimizing function. 

%%%%%%%%%%%%%%%%%%%%%%%%%% Objective function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bcd_sum = objectivefcn1(x)

mean_var_obs = readtable('mean_var_obs.txt'); 
mean_var_obs = table2array(mean_var_obs); 
% Read the table of observed summary statistics, change its format into an array.

%%%%%%%%%%%%%%%% Space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l1=0; 
l2=1;
x11=linspace(l1,l2,300);
h=abs(x11(2)-x11(1)); 
% Discretized dimensionless space, 0 to 1, 300 steps.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%% Time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=10;
dt = 0.0001;
time = 0:dt:T; 
% 10 time points, dt set to be 0.0001 in the numerical solver.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beta = 0;
eps = 0.01; 
% These two parameters are fixed at all time. 

%%%%%%%%%%%%%%% Initial condition %%%%%%%%%%%%%%%%%%%%%
n0 = repelem(0,length(x11));

for i = 1:length(x11)
    if x11(i)<=0.25
        n0(i)=exp(-(x11(i)^2)/eps);
    else 
        n0(i)=0;
    end 
end
% Initial condition of tumour cells.
n = n0;
% Initialize tumour cells vector. 

f0 = 1-0.5*n0;
% Initial condition of ECM.
f=f0;
% Initialize ECM vector.

m0=0.5*n0; 
% Initial condition of MDE
m = m0;
% Initialize MDE vector. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Numerical PDE solver %%%%%%%%%%%%%%%%%%%%%%%
p=1;

res = [];
% An empty array used to store the results.

while p*dt<=T
    f(2:length(x11)-1) = -x(4)*dt*m(2:length(x11)-1).*f(2:length(x11)-1)+f(2:length(x11)-1);
    m(2:length(x11)-1) = x(5)*(m(1:length(x11)-2)+m(3:length(x11))-2*m(2:length(x11)-1))*dt/(h^2)+x(6)*n(2:length(x11)-1)*dt-beta*m(2:length(x11)-1)*dt+m(2:length(x11)-1);
    n(2:length(x11)-1) = x(1)*(n(1:length(x11)-2)+n(3:length(x11))-2*n(2:length(x11)-1))*dt/(h^2)...
                         -x(2)*(n(3:length(x11))-n(2:length(x11)-1)).*(f(3:length(x11))-f(2:length(x11)-1))*dt/(h^2)...
                         -x(2)*n(2:length(x11)-1).*(f(1:length(x11)-2)+f(3:length(x11))...
                         -2*f(2:length(x11)-1))*dt/(h^2)+n(2:length(x11)-1)+x(3)*(1-f(2:length(x11)-1)-n(2:length(x11)-1)).*n(2:length(x11)-1)*dt;                     
    
    %%%%%%%%%%%%% Homogeneous & Zero Neumann boundary condition %%%%%%%%%%%%%%%%%
    n(1) = n(2);
    n(length(x11)) = n(length(x11)-1);
    
    f(1) = f(2);
    f(length(x11)) = f(length(x11)-1);
    
    m(1) = m(2);
    m(length(x11)) = m(length(x11)-1); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    if mod(p,10000)==0
        temp = round(round(time(p+1)*100)/100);
        res = [res;n;f;m];
        % Append the results of different variables at each timepoint into the array.
    end 
     
    p=p+1;
     
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Results rearrangements %%%%%%%%%%%%%%%%%%%%%%%%%
res_arr = zeros(30,300);

res_arr(1:10,:) = res(1:3:28,:);
res_arr(11:20,:) = res(2:3:29,:);
res_arr(21:30,:) = res(3:3:30,:); 
% Rearrange the result array, row 1-10, tumour cells, 
% row 11-20, ECM, row 21-30, MDE.

mean_var = zeros(900,2); 
% An 900*2 zero array, used to store the rearranged results as time series.
mean_var(1:300,1) = mean(res_arr(1:10,:)); 
% mean of the TC time series
mean_var(1:300,2) = var(res_arr(1:10,:)); 
% variance of the TC time series
mean_var(301:600,1) = mean(res_arr(11:20,:)); 
% mean of the ECM time series
mean_var(301:600,2) = var(res_arr(11:20,:)); 
% variance of the ECM time series
mean_var(601:900,1) = mean(res_arr(21:30,:)); 
% mean of the MDE time series
mean_var(601:900,2) = var(res_arr(21:30,:)); 
% variance of the MDE time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% Bhattacharyya distance calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bcd_vec=[];
% Create an empty vector to store the Bhattacharyya distance of the time series in relation to their
% corresponding observed ones.

for j = 1:900 
    bcd = 0.25*log(0.25*((mean_var(j,2)/mean_var_obs(j,2))+(mean_var_obs(j,2)/mean_var(j,2))+2))+0.25*(((mean_var(j,1)-mean_var_obs(j,1))^2)/(mean_var_obs(j,2)+mean_var(j,2)));
    bcd_vec = [bcd_vec,bcd];
end 
% 900 time series for every set of parameters, so 900 singular values will be appended into the vector
% bcd_vec.

inv_index = find(bcd_vec == Inf);
inv_term = length(inv_index);

bcd_vec_2 = bcd_vec;
bcd_vec_2(inv_index) = []; 
% Locate the invalid terms, exclude them from the vector. 

bcd_sum = sum(bcd_vec_2); 
% Sum up the Bhattacharya distance of the time series, our aim is to find the parameter values that can
% minimize such result. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
