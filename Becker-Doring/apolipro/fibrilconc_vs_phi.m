%fribil concentration vs phi with data
close all
data_phi=[0;.036;0.073;0.11];
data_F=[41.7;44.2;47.2;47.9];
c0=0.4;

%data_A4=5.957*data_A4 - 12.46*(data_A4.^2);



N=100;
l=101;
phi=linspace(0,0.15,l);
cfib=zeros(l,1);
lng=zeros(l,1);
gamma=zeros(l,1);
goa=zeros(l,1);
alpha=zeros(l,1);

c1_eq=zeros(l,1);
c_0=.4;

%scale
k=8.75*10^-3;
s=0;
data_F=k*(data_F-s);

r=1;
rc=3;
R=r/rc;
r_sc=.75*r;
Rsc=r_sc/rc;
a=.2;
b=10^-6;
aob=a/b;
mass=figure;
fibrils=figure;
for i=1:l
    z=phi(i)/(1-phi(i));
    lng(i)=log(1-phi(i))+(R^3 + 3*R^2 + 3*R)*z +(3*R^3 + 4.5*R^2)*z^2 + (3*R^2)*z^3;
    gamma(i)=exp(lng(i));
    lnalpha=(2/3)*(r/r_sc)^3 * (1.5*(Rsc^2+2*Rsc+1)*z+4.5*(Rsc^2+Rsc)*z^2+4.5*Rsc^2*z^3);
    alpha(i)=exp(lnalpha);
    
    lambda=a*gamma(i)/(b*alpha(i));
    c1_eq(i)=2*c_0/(1+2*lambda*c_0+(1+4*lambda*c_0)^0.5);
end

for i=1:l
    goa=gamma(i)/alpha(i);
    cfib(i)=c_0-c1_eq(i);
end

hold on
figure(fibrils)
plot(phi,cfib);
scatter(data_phi,data_F);