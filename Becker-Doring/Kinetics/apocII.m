close all
clear
% beckdoring script
N=1000;
ICs=zeros(N,1);
ICs(1)=.4;
ICderiv=zeros(N,1);
phi=0.0;
b=5*10^-6; % phi = 0 10^-6
factor=2; % phi=0 2.2*10^6
k_nuc=b*factor; %phi=0 2.5*10^4
k_dnuc=b;
scale=1; % time scale

nc=4; %critical nucleus
r=2.1; %protein raidus
rc=.8; %crowder radius ~.8 for phi=0.01875
%
rsc=2.4; %spherocylinder radius


%data=load('act_dex_phi075.txt');
data=load('apo_phi_0.csv');
data(:,1)=data(:,1)*60;

%options = odeset('MaxStep', scale*data(length(data),1)/5000); %force 1000 time steps 
%no crowds
%sol =ode15s(@(t,c)beckdoring(t,c,phi,b,factor,k_nuc),[0,data(length(data),1)],ICs,options);

%crowds


%rc = ~ 1.75, 1.85, 2.25 for .0375, .075, .15 respectively
t=data(length(data));
sol = ode15s(@(t,c)beckdoring_rc(t,c,phi,b,factor,k_nuc,k_dnuc,rc,r,rsc,nc),[0,t],ICs);

T1=data(:,1)*scale;
data(:,1)=data(:,1)*scale;
Y1=deval(sol,T1);
Y1=Y1.';

t=length(T1);
P=zeros(t,1);
M=zeros(t,1);
L=zeros(t,1);
mass=zeros(t,1);
for i=2:N
    P(:)=P(:)+Y1(:,i);
    M(:)=M(:)+i*Y1(:,i);
end
for i=1:N
    mass(:)=mass(:)+i*Y1(:,i);
end
L=M./P;

dextran=figure;

%scale for phi=0 1.28*10^-7 N=1000
%k=1.3*10^-7; %phi=0 N=500
%k=M(length(M))/data(length(data),2);
%s=data(length(data),2)-(M(length(M))/k);

%data=data-data(1,2);
%k=(M(length(M)))/data(length(data),2); %this method gives k=1.277*10^-7
%k=1.2772*10^-7;

% k=(3.66*10^-6)/(data(length(data),2)-data(1,2));
% s=k*data(1,2);





k=1.6;
s=0;

data(1,2);
figure(dextran)
data(:,2)=k*data(:,2)-s;
data(length(data),2);
M(length(M));


plot(T1,M)
hold on
scatter(data(:,1),data(:,2),'s','filled')