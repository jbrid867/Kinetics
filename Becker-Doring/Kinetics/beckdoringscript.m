close all
clear
% beckdoring script
N=50;
ICs=zeros(N,1);
ICs(1)=5*10^-6;
ICderiv=zeros(N,1);
phi=0.0;
b=2*10^-2; % phi = 0 10^-2
factor=7.2*10^5; % phi=0 2.2*10^6
k_nuc=5.6*10^4; %phi=0 2.5*10^4
k_dnuc=b;
scale=1; % time scale

nc=3; %critical nucleus
r=2.5; %protein raidus
rc=1.85; %crowder radius
rsc=5; %spherocylinder radius


%data=load('act_dex_phi075.txt');
data=load('dez_data2');
data(:,1)=data(:,1)*60;

%options = odeset('MaxStep', scale*data(length(data),1)/5000); %force 1000 time steps 
%no crowds
%sol =ode15s(@(t,c)beckdoring(t,c,phi,b,factor,k_nuc),[0,data(length(data),1)],ICs,options);

%crowds


%rc = ~ 1.75, 1.85, 2.25 for .0375, .075, .15 respectively

sol = ode23s(@(t,c)beckdoring_rc(t,c,phi,b,factor,k_nuc,k_dnuc,rc,r,rsc,nc),[0,data(length(data),1)],ICs);

T1=data(:,1)*scale;
data(:,1)=data(:,1)*scale;
Y1=deval(sol,T1);
Y1=Y1.';

t=length(T1);
full_time=T1;
P=zeros(t,1);
M=zeros(t,1);
L=zeros(t,1);
mass=zeros(t,1);
for i=3:N
    P(:)=P(:)+Y1(:,i);
    M(:)=M(:)+i*Y1(:,i);
end
for i=3:N
    mass(:)=mass(:)+i*Y1(:,i);
end
T11=T1;
L=mass./P;
L(1)=0;

dextran=figure;
contours=figure;
number=figure;

%scale for phi=0 1.28*10^-7 N=1000
%k=1.3*10^-7; %phi=0 N=500
%k=M(length(M))/data(length(data),2);
%s=data(length(data),2)-(M(length(M))/k);

%data=data-data(1,2);
%k=(M(length(M)))/data(length(data),2); %this method gives k=1.277*10^-7
%k=1.2772*10^-7;

% k=(3.66*10^-6)/(data(length(data),2)-data(1,2));
% s=k*data(1,2);


%surf(3:900,T1,-log(Y1(:,3:900)));
%set(h,'LineStyle','none')


k=1.0624*10^-7;
s=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data(1,2);
figure(dextran)
data(:,2)=k*data(:,2)-s;
data(length(data),2);
M(length(M));


plot(T1./3600,M.*10^6)
hold on
scatter(data(:,1)./3600,data(:,2).*10^6,'s','filled')
figure(contours)
subplot(2,2,1)
contourf(T1./3600,3:N,Y1(:,3:N).',15);
hold on
plot(T1./3600,L,'--r','LineWidth',2)
tit=title('\phi = 0')
xl=xlabel('t(h)');
yl=ylabel('r');
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([tit, xl, yl], ...
    'FontName'   , 'AvantGarde');
set([gca]             , ...
    'FontSize'   , 12           );
set([tit, xl, yl]  , ...
    'FontSize'   , 13         );
set(gca, 'CLim', [0, 2.5*10^-10]);


figure(number)
plot(full_time,P,'LineWidth',2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phi=0.075
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

phi=0.075;

scale=1; % time scale
nc=3; %critical nucleus
r=2.5; %protein raidus
rc=1.95; %crowder radius
rsc=5; %spherocylinder radius


data=load('act_dex_phi075.txt');
%data=load('dez_data2');
data(:,1)=data(:,1)*60;

%options = odeset('MaxStep', scale*data(length(data),1)/5000); %force 1000 time steps 
%no crowds
%sol =ode15s(@(t,c)beckdoring(t,c,phi,b,factor,k_nuc),[0,data(length(data),1)],ICs,options);

%crowds


%rc = ~ 1.75, 1.85, 2.25 for .0375, .075, .15 respectively

sol = ode23s(@(t,c)beckdoring_rc(t,c,phi,b,factor,k_nuc,k_dnuc,rc,r,rsc,nc),[0,full_time(length(full_time))],ICs);

T1=data(:,1)*scale;
data(:,1)=data(:,1)*scale;
Y1=deval(sol,T1);
Y1=Y1.';
Y2=deval(sol,full_time);
Y2=Y2.';

t=length(T1);
t2=length(full_time);
P=zeros(t2,1);
M=zeros(t,1);
M2=zeros(t2,1);

L=zeros(t,1);
mass=zeros(t,1);
for i=3:N
    P(:)=P(:)+Y2(:,i);
    M(:)=M(:)+i*Y1(:,i);
    M2(:)=M2(:)+i*Y2(:,i);
end
for i=1:N
    mass(:)=mass(:)+i*Y1(:,i);
end
L=M2./P;
L(1)=0;

%scale for phi=0 1.28*10^-7 N=1000
%k=1.3*10^-7; %phi=0 N=500
%k=M(length(M))/data(length(data),2);
%s=data(length(data),2)-(M(length(M))/k);

%data=data-data(1,2);
%k=(M(length(M)))/data(length(data),2); %this method gives k=1.277*10^-7
%k=1.2772*10^-7;

% k=(3.66*10^-6)/(data(length(data),2)-data(1,2));
% s=k*data(1,2);





k=1.0624*10^-7;
s=0;


data(:,2)=k*data(:,2)-s;
data(length(data),2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(dextran)
plot(full_time./3600,M2.*10^6,'r')
scatter(data(:,1)./3600,data(:,2).*10^6,'+')

figure(contours)
subplot(2,2,3);
contourf(full_time./3600,3:N,Y2(:,3:N).',15);
hold on
plot(full_time./3600,L,'--r','LineWidth',2)
tit=title('\phi = 0.075')
xl=xlabel('t(h)');
yl=ylabel('r');
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([tit, xl, yl], ...
    'FontName'   , 'AvantGarde');
set([gca]             , ...
    'FontSize'   , 12           );
set([tit, xl, yl]  , ...
    'FontSize'   , 13         );
set(gca, 'CLim', [0, 2.5*10^-10]);

figure(number)
hold on
plot(full_time,P,'r','LineWidth',2);
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phi=0.0375
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi=0.0375;



scale=1; % time scale
nc=3; %critical nucleus
r=2.5; %protein raidus
rc=1.85; %crowder radius
rsc=5; %spherocylinder radius


data=load('act_dex_phi0375.txt');
%data=load('dez_data2');
data(:,1)=data(:,1)*60;

%options = odeset('MaxStep', scale*data(length(data),1)/5000); %force 1000 time steps 
%no crowds
%sol =ode15s(@(t,c)beckdoring(t,c,phi,b,factor,k_nuc),[0,data(length(data),1)],ICs,options);

%crowds


%rc = ~ 1.75, 1.85, 2.25 for .0375, .075, .15 respectively

sol = ode23s(@(t,c)beckdoring_rc(t,c,phi,b,factor,k_nuc,k_dnuc,rc,r,rsc,nc),[0,full_time(length(full_time),1)],ICs);

T1=data(:,1)*scale;
data(:,1)=data(:,1)*scale;
Y1=deval(sol,T1);
Y1=Y1.';
Y2=deval(sol,full_time);
Y2=Y2.';

t=length(T1);
t2=length(full_time);
P=zeros(t2,1);
M=zeros(t,1);
M2=zeros(t2,1);

L=zeros(t,1);
mass=zeros(t,1);
for i=3:N
    P(:)=P(:)+Y2(:,i);
    M(:)=M(:)+i*Y1(:,i);
    M2(:)=M2(:)+i*Y2(:,i);
end
for i=1:N
    mass(:)=mass(:)+i*Y1(:,i);
end
L=M2./P;
L(1)=0;


%scale for phi=0 1.28*10^-7 N=1000
%k=1.3*10^-7; %phi=0 N=500
%k=M(length(M))/data(length(data),2);
%s=data(length(data),2)-(M(length(M))/k);

%data=data-data(1,2);
%k=(M(length(M)))/data(length(data),2); %this method gives k=1.277*10^-7
%k=1.2772*10^-7;

% k=(3.66*10^-6)/(data(length(data),2)-data(1,2));
% s=k*data(1,2);





k=1.0624*10^-7;
s=0;


data(:,2)=k*data(:,2)-s;
data(length(data),2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(dextran)
plot(full_time./3600,M2.*10^6,'m')
scatter(data(:,1)./3600,data(:,2).*10^6,'o','filled')

figure(contours)
subplot(2,2,2);
contourf(full_time./3600,3:N,Y2(:,3:N).',15);
hold on
plot(full_time./3600,L,'--r','LineWidth',2)
tit=title('\phi = 0.0375');
xl=xlabel('t(h)');
yl=ylabel('r');
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([tit, xl, yl], ...
    'FontName'   , 'AvantGarde');
set([gca]             , ...
    'FontSize'   , 12           );
set([tit, xl, yl]  , ...
    'FontSize'   , 13         );
set(gca, 'CLim', [0, 2.5*10^-10]);

figure(number)
hold on
plot(full_time,P,'m','LineWidth',2);
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phi=0.15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi=0.15;

scale=1; % time scale
nc=3; %critical nucleus
r=2.5; %protein raidus
rc=2.2; %crowder radius
rsc=5; %spherocylinder radius


data=load('act_dex_phi150.txt');
%data=load('dez_data2');
data(:,1)=data(:,1)*60;

%options = odeset('MaxStep', scale*data(length(data),1)/5000); %force 1000 time steps 
%no crowds
%sol =ode15s(@(t,c)beckdoring(t,c,phi,b,factor,k_nuc),[0,data(length(data),1)],ICs,options);

%crowds


%rc = ~ 1.75, 1.85, 2.25 for .0375, .075, .15 respectively
opts=odeset('NonNegative',1);
sol = ode23s(@(t,c)beckdoring_rc(t,c,phi,b,factor,k_nuc,k_dnuc,rc,r,rsc,nc),[0,full_time(length(full_time),1)],ICs,opts);

T1=data(:,1)*scale;
data(:,1)=data(:,1)*scale;
Y1=deval(sol,T1);
Y1=Y1.';
Y2=deval(sol,full_time);
Y2=Y2.';

t=length(T1);
t2=length(full_time);
P=zeros(t2,1);
M=zeros(t,1);
M2=zeros(t2,1);

L=zeros(t,1);
mass=zeros(t,1);
for i=3:N
    P(:)=P(:)+Y2(:,i);
    M(:)=M(:)+i*Y1(:,i);
    M2(:)=M2(:)+i*Y2(:,i);
end
for i=1:N
    mass(:)=mass(:)+i*Y1(:,i);
end
L=M2./P;
L(1)=0;

%scale for phi=0 1.28*10^-7 N=1000
%k=1.3*10^-7; %phi=0 N=500
%k=M(length(M))/data(length(data),2);
%s=data(length(data),2)-(M(length(M))/k);

%data=data-data(1,2);
%k=(M(length(M)))/data(length(data),2); %this method gives k=1.277*10^-7
%k=1.2772*10^-7;

% k=(3.66*10^-6)/(data(length(data),2)-data(1,2));
% s=k*data(1,2);





k=1.0624*10^-7;
s=0;



data(:,2)=k*data(:,2)-s;
data(length(data),2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(dextran)
plot(full_time./3600,M2.*10^6,'k')
scatter(data(:,1)./3600,data(:,2).*10^6,'x')




%tit=title('\phi = 0, 0.0375, 0.075, 0.15');
xlab=xlabel('t(hours)');
ylab=ylabel('c_{fib}');
ylab=ylabel('c_{fib} (\muM)');

set(gca,'fontsize',14);
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([xlab, ylab], ...
    'FontName'   , 'AvantGarde');
set([gca]             , ...
    'FontSize'   , 12           );
set([xlab, ylab]  , ...
    'FontSize'   , 12         );

figure(contours)
subplot(2,2,4);
contourf(full_time./3600,3:N,Y2(:,3:N).',15);
hold on
plot(full_time./3600,L,'--r','LineWidth',2)
tit=title('\phi = 0.15')
xl=xlabel('t(h)');
yl=ylabel('r');
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([tit, xl, yl], ...
    'FontName'   , 'AvantGarde');
set([gca]             , ...
    'FontSize'   , 12           );
set([tit, xl, yl]  , ...
    'FontSize'   , 13         );
set(gca, 'CLim', [0, 2.5*10^-10]);

axes('Position', [0.05 0.05 .95 .95], 'Visible', 'off');
h=ylabel('\mu M');
set(gca, 'CLim', [0, 2.5*10^-11]);

colorbar

figure(number)
hold on
plot(full_time,P,'k','LineWidth',2);
hold off

figure(contours)
for ii = 1:4
    subplot(2,2,ii)
    text(-.5,0,['(',char(ii+96),')'],...
        'color','k',...
        'fontw','b')
end

hold off



%plot(T1./3600,L,'--r','LineWidth',2);

