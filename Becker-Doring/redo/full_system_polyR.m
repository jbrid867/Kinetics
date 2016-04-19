% full theory script
close all
clear all
r=2.5; %protein raidus
%rc=1.35; %crowder radius
rsc=5; %spherocylinder radius
phimax=0.15; %highest volume fraction
l=101;
phi=linspace(0,phimax,l);
dextran=figure;

c0=5*10^-6; %initial concentration

N=100; % max tracked polymer size

c=zeros(l,N);
cfib=zeros(1,l);

a=2.4*10^-2 * 2.15*10^6; %a=2.15*10^6
b=1; %ass and diss rates

alpha=zeros(1,l);
lnalpha=alpha;
lng=alpha;
gamma=alpha;
z=alpha;
gamma(1)=1; alpha(1)=1;

c(1,1) = 2*c0/(1+2*(a/b)*c0 + (1+4*(a/b)*c0)^0.5); % no crowder monomer mass
c(1,2) = (a/b)*c(1,1)^2;
cfib(1)=c0-c(1,1)-c(1,2);


% fit for winter dextran data. FITS VERIFIED BY LEAST SQUARE
for i=2:l
    rc=23.7*phi(i)+1.72;
    R=r/rc;
    Rsc=rsc/rc;
    z(i)=phi(i)/(1-phi(i));
    lng(i)=log(1-phi(i))+(R^3+3*R^2+R)*z(i)+(3*R^3+4.5*R^2)*z(i)^2+3*R^3*z(i)^3;
    lnalpha(i)=(2/3)*(r/rsc)^3*(1.5*(Rsc^2+2*Rsc+1)*z(i)+4.5*(Rsc^2+Rsc)*z(i)^2+4.5*Rsc^2*z(i)^3);
    gamma(i)=exp(lng(i));
    alpha(i)=exp(lnalpha(i));
    lambda=gamma(i)*a/(alpha(i)*b);
    c(i,1)=2*c0/(1+2*lambda*c0 + (1+4*lambda*c0)^0.5);
    c(i,2)=lambda*c(i,1)^2;
    cfib(i)=c0-c(i,1)-c(i,2);
end

data_phi_1=[0;.0375;.075;.15];
data_F_1=[30;36.3;41.2;45.2];

%scale
k=.038*10^-7;
s=30-(1*10^-6)/k;
figure(dextran)
data_F_1=k*(data_F_1-s);
data=scatter(data_phi_1,data_F_1,'s','filled');
hold on
plot(phi,cfib);
tit=title('');
xlab=xlabel('\phi');
ylab=ylabel('c_{fib}');
ylab=ylabel('c_{fib} (\muM)');

set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([tit, xlab, ylab], ...
    'FontName'   , 'AvantGarde');
set([gca]             , ...
    'FontSize'   , 8           );
set([xlab, ylab]  , ...
    'FontSize'   , 10          );
set( tit                    , ...
    'FontSize'   , 12          , ...
    'FontWeight' , 'bold'      );
cfib(2)

%% fit for winter TMAO data
tmao=figure;
r=2.5; %protein raidus
rc=1.06; %crowder radius
rsc=5; %spherocylinder radius

alpha=zeros(1,l);
lnalpha=alpha;
lng=alpha;
gamma=alpha;
z=alpha;
gamma(1)=1; alpha(1)=1;

c(1,1) = 2*c0/(1+2*(a/b)*c0 + (1+4*(a/b)*c0)^0.5); % no crowder monomer mass
c(1,2) = (a/b)*c(1,1)^2;
cfib(1)=c0-c(1,1)-c(1,2);
R=r/rc;
Rsc=rsc/rc;

for i=2:l
    z(i)=phi(i)/(1-phi(i));
    lng(i)=log(1-phi(i))+(R^3+3*R^2+R)*z(i)+(3*R^3+4.5*R^2)*z(i)^2+3*R^3*z(i)^3;
    lnalpha(i)=(2/3)*(r/rsc)^3*(1.5*(Rsc^2+2*Rsc+1)*z(i)+4.5*(Rsc^2+Rsc)*z(i)^2+4.5*Rsc^2*z(i)^3);
    gamma(i)=exp(lng(i));
    alpha(i)=exp(lnalpha(i));
    lambda=gamma(i)*a/(alpha(i)*b);
    c(i,1)=2*c0/(1+2*lambda*c0 + (1+4*lambda*c0)^0.5);
    c(i,2)=lambda*c(i,1)^2;
    cfib(i)=c0-c(i,1)-c(i,2);
end

data_Tphi=[0;0.5*72.7;1*72.1;2*70.9]./1000; % molarity * specific molar volume (rosgen)
data_TF=[29.3;56.7;62.7;74.6];


%scale
k=.132*10^-7;
s=data_TF(1)-(4.4*10^-6)/k;
data_TF=k*(data_TF-s);

figure(tmao)
data=scatter(data_Tphi,data_TF,'s','filled');
hold on
plot(phi,cfib);
tit=title('');
xlab=xlabel('\phi');
ylab=ylabel('c_{fib}');
ylab=ylabel('c_{fib} (\muM)');

set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([tit, xlab, ylab], ...
    'FontName'   , 'AvantGarde');
set([gca]             , ...
    'FontSize'   , 8           );
set([xlab, ylab]  , ...
    'FontSize'   , 10          );
set( tit                    , ...
    'FontSize'   , 12          , ...
    'FontWeight' , 'bold'      );


%% fit for hatters & minton turbidity data ApoC-II
minton=figure;
c0=4.5*10^-2; %from paper
r=2.1; %protein raidus
rc=.95; %crowder radius
rsc=2.4; %spherocylinder radius
phi=linspace(0,0.111,l);

alpha=zeros(1,l);
lnalpha=alpha;
lng=alpha;
gamma=alpha;
z=alpha;
gamma(1)=1; alpha(1)=1;
a=0.5;
b=1; %process limited by nucleation

c(1,1) = 2*c0/(1+2*(a/b)*c0 + (1+4*(a/b)*c0)^0.5); % no crowder monomer mass
c(1,2) = (a/b)*c(1,1)^2;
cfib(1)=c0-c(1,1)-c(1,2);
R=r/rc;
Rsc=rsc/rc;

for i=2:l
    z(i)=phi(i)/(1-phi(i));
    lng(i)=log(1-phi(i))+(R^3+3*R^2+R)*z(i)+(3*R^3+4.5*R^2)*z(i)^2+3*R^3*z(i)^3;
    lnalpha(i)=(2/3)*(r/rsc)^3*(1.5*(Rsc^2+2*Rsc+1)*z(i)+4.5*(Rsc^2+Rsc)*z(i)^2+4.5*Rsc^2*z(i)^3);
    gamma(i)=exp(lng(i));
    alpha(i)=exp(lnalpha(i));
    lambda=gamma(i)*a/(alpha(i)*b);
    c(i,1)=2*c0/(1+2*lambda*c0 + (1+4*lambda*c0)^0.5);
    c(i,2)=lambda*c(i,1)^2;
    cfib(i)=c0-c(i,1)-c(i,2);
end

data_Mphi=[0;0.018;0.036;0.055;0.073;0.096;0.11]; %converted with specific volume from paper
data_MT=[.073;0.0894;0.1151;0.13;0.166;0.24;0.332];


%scale
k=4.15*10^-2;
s=data_MT(1)-(cfib(1))/k;
data_MT=k*(data_MT-s);

figure(minton)
data=scatter(data_Mphi,data_MT,'s','filled');
hold on
plot(phi,cfib);
tit=title('');
xlab=xlabel('\phi');
ylab=ylabel('c_{fib}');
ylab=ylabel('c_{fib} (g/L)');

set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([tit, xlab, ylab], ...
    'FontName'   , 'AvantGarde');
set([gca]             , ...
    'FontSize'   , 8           );
set([xlab, ylab]  , ...
    'FontSize'   , 10          );
set( tit                    , ...
    'FontSize'   , 12          , ...
    'FontWeight' , 'bold'      );


