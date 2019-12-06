clc
clear all
close all
 
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 24)
set(0,'defaultaxeslinewidth',1)
set(0,'defaultpatchlinewidth',1)
set(0,'defaultlinelinewidth',4)
set(0,'defaultTextInterpreter','latex')

paras = readtable('Round 10 parameters 10000 all 3.txt');
mean_var_obs = readtable('mean_var_obs.txt');
mean_var_obs = table2array(mean_var_obs);

% Space
l1=0;
l2=1;
x11=linspace(l1,l2,300);
h=abs(x11(2)-x11(1));
 
% Timec
T=10;
dt = 0.0001;
time = 0:dt:T;
tv=0;
 
% Parameters

%dn=0.0139215711;
%gamma=0.235221518;
%ita=10.317482;
%dm=0.0315314271;
%alpha=0.13782755;

dn = table2array(paras(:,2));


gamma = table2array(paras(:,3));


ita = table2array(paras(:,4));


dm = table2array(paras(:,5));


alpha = table2array(paras(:,6));

r = table2array(paras(:,7));

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

%%%Inits.matrix%%%
inits = [n;f;m];
%%%%%%%%%%%%%%%%%%

%%%%% Initial plot %%%%%%%%%%%%%%%%%%%
figure
plot(x11,n,'b-',x11,f,'k-',x11,m,'g-')
axis([min(x11) max(x11) 0 max(f)+0.1])
axis square
xlabel('$x_1$')
title(['$System$ at $t=$',num2str(0)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Open file %%%%%%%%%%%%%%%%%
fileID = fopen('B-C distance APF r1 p4.txt','w');
fclose(fileID);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%bcd_vec = [];

bcd_sum=[];
inv_term_vec = [];

for i = 7501:10000
n = inits(1,:);
f = inits(2,:);
m = inits(3,:);

p=1;
%%% Results matrix%%%%%
res = [];
%%%%%%%%%%%%%%%%%%%%%%%
while p*dt<=T
    f(2:length(x11)-1) = -ita(i)*dt*m(2:length(x11)-1).*f(2:length(x11)-1)+f(2:length(x11)-1);
    m(2:length(x11)-1) = dm(i)*(m(1:length(x11)-2)+m(3:length(x11))-2*m(2:length(x11)-1))*dt/(h^2)+alpha(i)*n(2:length(x11)-1)*dt-beta*m(2:length(x11)-1)*dt+m(2:length(x11)-1);
    n(2:length(x11)-1) = dn(i)*(n(1:length(x11)-2)+n(3:length(x11))-2*n(2:length(x11)-1))*dt/(h^2)...
                         -gamma(i)*(n(3:length(x11))-n(2:length(x11)-1)).*(f(3:length(x11))-f(2:length(x11)-1))*dt/(h^2)...
                         -gamma(i)*n(2:length(x11)-1).*(f(1:length(x11)-2)+f(3:length(x11))...
                         -2*f(2:length(x11)-1))*dt/(h^2)+n(2:length(x11)-1)+r(i)*(1-f(2:length(x11)-1)-n(2:length(x11)-1)).*n(2:length(x11)-1)*dt;                     
    
    %Homogeneous & Zero Neumann
    n(1) = n(2);
    n(length(x11)) = n(length(x11)-1);
    
    f(1) = f(2);
    f(length(x11)) = f(length(x11)-1);
    
    m(1) = m(2);
    m(length(x11)) = m(length(x11)-1); 

     
     
     
    %% Plots
    if mod(p,10000)==0
        
        temp = round(round(time(p+1)*100)/100);
        
        %disp(temp);
     
        %figure
        %plot(x11,n,'b-',x11,f,'k-',x11,m,'g-')
        %axis([min(x11) max(x11) 0 1])
        %axis square
        %xlabel('$x$')
        %title(['$System$ at $t=$',num2str(round(time(p+1)*100)/100)])
        
        res = [res;n;f;m];
        
        %p;
        %tv=[tv p*dt];
        %clf
        %plot(x11,n,'b-',x11,f,'k-',x11,m,'g-')
        %axis([min(x11) max(x11) 0 1])
        %axis square
        %xlabel('$x_1$')
        %title(['$System$ at $t=$',num2str(round(time(p+1)*100)/100)])
        %drawnow
    end
     
    p=p+1;
     
end

%size_res = size(res);
res_arr = zeros(30,300);

res_arr(1:10,:) = res(1:3:28,:);
res_arr(11:20,:) = res(2:3:29,:);
res_arr(21:30,:) = res(3:3:30,:); % rearrange the result matrix, 1-5, TC, 6-10, ECM, 11-15, MDE.

mean_var = zeros(900,2); 
%size_mean_var = size(mean_var);

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
bcd_vec_2(inv_index) = [];

bcd_sum = [bcd_sum,sum(bcd_vec_2)];
inv_term_vec = [inv_term_vec,inv_term];

%bcd_vec = [bcd_vec,bcd];

fileID = fopen('B-C distance APF r1 p4.txt','a');
fprintf(fileID,'%4d %5.4f\r\n',i,sum(bcd_vec_2));
fclose(fileID);

disp(i)

end 