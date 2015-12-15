close all
clear
% beckdoring script
N=1000;
ICs=zeros(N,1);
ICs(1)=5*10^-6;
ICderiv=zeros(N,1);
phi=0.075;
b=1*10^-2; % phi = 0 10^-2
factor=10^7; % phi=0 2.2*10^6
k_nuc=10^4; %phi=0 2.5*10^4
scale=1; % time scale


data=load('act_dex_phi075.txt');
%data=load('dez_data2');
data(:,1)=data(:,1)*60;

%options = odeset('MaxStep', scale*data(length(data),1)/5000); %force 1000 time steps 
%no crowds
%sol =ode15s(@(t,c)beckdoring(t,c,phi,b,factor,k_nuc),[0,data(length(data),1)],ICs,options);

%crowds


sol = ode15s(@(t,c)beckdoring_rc(t,c,phi,b,factor,k_nuc,1.7),[0,scale*data(length(data),1)],ICs);

T1=data(:,1)*scale;
data(:,1)=data(:,1)*scale;
Y1=deval(sol,T1);
Y1=Y1.';

t=length(T1);
P=zeros(t,1);
M=zeros(t,1);
L=zeros(t,1);
mass=zeros(t,1);
for i=3:N
    P(:)=P(:)+Y1(:,i);
    M(:)=M(:)+i*Y1(:,i);
end
for i=1:N
    mass(:)=mass(:)+i*Y1(:,i);
end
L=M./P;

dextran=figure;

%scale for phi=0 1.28*10^-7 N=1000
%k=1.128841*10^-7; %phi=0 N=500
k=M(length(M))/data(length(data),2);
s=data(length(data),2)-(M(length(M))/k);

figure(dextran)
data(:,2)=k*data(:,2);
data(length(data),2)



plot(T1,M)
hold on
scatter(data(:,1),data(:,2),'o')
hold off