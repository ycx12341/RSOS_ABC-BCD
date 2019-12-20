x0 = [0.0100345,0.1325,6.25,12.5,0.01655,0.125];
options = optimset('Display','iter');
[x,fval,exitflag,output] = fminsearch(@objectivefcn1,x0,options);

function bcd_sum = objectivefcn1(x)

mean_var_obs = readtable('mean_var_obs.txt'); 
mean_var_obs = table2array(mean_var_obs); % Read the table of observed summary statistics, change it into an array.

% Space
l1=0; 
l2=1;
x11=linspace(l1,l2,300);
h=abs(x11(2)-x11(1)); % Discretized dimensionless space, 0 to 1, 300 steps.
 
% Time
T=10;
dt = 0.0001;
time = 0:dt:T; % 10 time points, dt set to be 0.0001 in the numerical solver.

beta = 0;
eps = 0.01; 

% Initial condition
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

p=1;

res = [];

while p*dt<=T
    f(2:length(x11)-1) = -x(4)*dt*m(2:length(x11)-1).*f(2:length(x11)-1)+f(2:length(x11)-1);
    m(2:length(x11)-1) = x(5)*(m(1:length(x11)-2)+m(3:length(x11))-2*m(2:length(x11)-1))*dt/(h^2)+x(6)*n(2:length(x11)-1)*dt-beta*m(2:length(x11)-1)*dt+m(2:length(x11)-1);
    n(2:length(x11)-1) = x(1)*(n(1:length(x11)-2)+n(3:length(x11))-2*n(2:length(x11)-1))*dt/(h^2)...
                         -x(2)*(n(3:length(x11))-n(2:length(x11)-1)).*(f(3:length(x11))-f(2:length(x11)-1))*dt/(h^2)...
                         -x(2)*n(2:length(x11)-1).*(f(1:length(x11)-2)+f(3:length(x11))...
                         -2*f(2:length(x11)-1))*dt/(h^2)+n(2:length(x11)-1)+x(3)*(1-f(2:length(x11)-1)-n(2:length(x11)-1)).*n(2:length(x11)-1)*dt;                     
    
    %Homogeneous & Zero Neumann boundary condition
    n(1) = n(2);
    n(length(x11)) = n(length(x11)-1);
    
    f(1) = f(2);
    f(length(x11)) = f(length(x11)-1);
    
    m(1) = m(2);
    m(length(x11)) = m(length(x11)-1); 

     
    if mod(p,10000)==0
        temp = round(round(time(p+1)*100)/100);
        res = [res;n;f;m];
    end 
     
    p=p+1;
     
end

res_arr = zeros(30,300);

res_arr(1:10,:) = res(1:3:28,:);
res_arr(11:20,:) = res(2:3:29,:);
res_arr(21:30,:) = res(3:3:30,:); % Rearrange the result array, row 1-10, TC, row 11-20, ECM, row 21-30, MDE.

mean_var = zeros(900,2); 

mean_var(1:300,1) = mean(res_arr(1:10,:)); % mean of the TC time series
mean_var(1:300,2) = var(res_arr(1:10,:)); % variance of the TC time series

mean_var(301:600,1) = mean(res_arr(11:20,:)); % mean of the ECM time series
mean_var(301:600,2) = var(res_arr(11:20,:)); % variance of the ECM time series

mean_var(601:900,1) = mean(res_arr(21:30,:)); % mean of the MDE time series
mean_var(601:900,2) = var(res_arr(21:30,:)); % variance of the MDE time series

bcd_vec=[];

for j = 1:900 
    bcd = 0.25*log(0.25*((mean_var(j,2)/mean_var_obs(j,2))+(mean_var_obs(j,2)/mean_var(j,2))+2))+0.25*(((mean_var(j,1)-mean_var_obs(j,1))^2)/(mean_var_obs(j,2)+mean_var(j,2)));
    bcd_vec = [bcd_vec,bcd];
end 

inv_index = find(bcd_vec == Inf);
inv_term = length(inv_index);

bcd_vec_2 = bcd_vec;
bcd_vec_2(inv_index) = []; % Locate the invalid terms, exclude them from the B-C distance calculation. 

bcd_sum = sum(bcd_vec_2); % Sum up the Bhattacharya distance of the time series.

end 
