clc
clear all
close all
 
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 24)
set(0,'defaultaxeslinewidth',1)
set(0,'defaultpatchlinewidth',1)
set(0,'defaultlinelinewidth',4)
set(0,'defaultTextInterpreter','latex')

mean_var_obs = readtable('mean_var_obs.txt');
mean_var_obs = table2array(mean_var_obs);

%%%%%%%% Space %%%%%%%%%%%%%%%%%%%%%%%
l1=0;
l2=1;
x11=linspace(l1,l2,300);
h=abs(x11(2)-x11(1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Time %%%%%%%%%%%%%%%%%%%%%%%%%
T=10;
dt = 0.0001;
time = 0:dt:T;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%

dn = 0.001550216;
gamma = 0.04672252;
ita = 9.796103;
dm = 0.01065483;
alpha = 0.0998378;
r = 4.84149;

beta = 0;
eps = 0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Initial condition %%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Array of initial values %%%%%%%
inits = [n;f;m];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Open file %%%%%%%%%%%%%%%%%
fileID = fopen('.txt','w');
fclose(fileID);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

res = [];

p=1;

n=n0;
m=m0;
f=f0;

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res_arr = zeros(30,300);

res_arr(1:10,:) = res(1:3:28,:);
res_arr(11:20,:) = res(2:3:29,:);
res_arr(21:30,:) = res(3:3:30,:);
mean_var = zeros(900,2); 

mean_var(1:300,1) = mean(res_arr(1:10,:));
mean_var(1:300,2) = var(res_arr(1:10,:)); 

mean_var(301:600,1) = mean(res_arr(11:20,:)); 
mean_var(301:600,2) = var(res_arr(11:20,:)); 

mean_var(601:900,1) = mean(res_arr(21:30,:)); 
mean_var(601:900,2) = var(res_arr(21:30,:)); 

bcd_vec=[];

for j = 1:900
    bcd = 0.25*log(0.25*((mean_var(j,2)/mean_var_obs(j,2))+(mean_var_obs(j,2)/mean_var(j,2))+2))+0.25*(((mean_var(j,1)-mean_var_obs(j,1))^2)/(mean_var_obs(j,2)+mean_var(j,2)));
    bcd_vec = [bcd_vec,bcd];
end 

inv_index = find(bcd_vec == Inf);
inv_term = length(inv_index);

bcd_vec_2 = bcd_vec;
bcd_vec_2(inv_index) = [];

bcd_sum = [bcd_sum,sum(bcd_vec_2)];
inv_term_vec = [inv_term_vec,inv_term];

fileID = fopen('.txt','a');
fprintf(fileID,'%4d %5.4f\r\n',i,sum(bcd_vec_2));
% The Bhattacharya distance of the parameter vector in relation to the
% reference parameter values will be written in the .txt file.
fclose(fileID);
