% calculate equilibrium mass of monomers
% THERE IS SOME KIND OF ERROR HERE AND I HAVE NO IDEA WHAT!!!!!!!!
N=100;
l=101;
phi=linspace(0,0.4,l);
c=zeros(l,1);
lng=zeros(l,1);
gamma=zeros(l,1);
alpha=zeros(l,1);



c1_eq=zeros(l,1);
c_0=.00001;

r=1;
rc=.5;
R=r/rc;
r_sc=r;
a=1;
b=10^-6;
Rsc=r_sc/rc;

for i=1:l
    z=phi(i)/(1-phi(i));
    lng(i)=log(1-phi(i))+(R^3 + 3*R^2 + 3*R)*z +(3*R^3 + 4.5*R^2)*z^2 + (3*R^2)*z^3;
    gamma(i)=exp(lng(i));
    lnalpha=(2/3)*(r/r_sc)^3 * (1.5*(Rsc^2+2*Rsc+1*z)+4.5*(Rsc^2+Rsc)^z^2+4.5*Rsc^2*z^3);
    alpha(i)=exp(lnalpha);
    
    lambda=a*gamma(i)/(b*alpha(i));
    c1_eq(i)=2*c_0/(1+2*lambda*c_0+(1+4*lambda*c_0)^0.5);
end
plot(phi,c1_eq,'b');
hold on

r=1;
rc=1;
R=r/rc;
r_sc=r;
Rsc=r_sc/rc;

for i=1:l
    z=phi(i)/(1-phi(i));
    lng(i)=log(1-phi(i))+(R^3 + 3*R^2 + 3*R)*z +(3*R^3 + 4.5*R^2)*z^2 + (3*R^2)*z^3;
    gamma(i)=exp(lng(i));
    lnalpha=(2/3)*(r/r_sc)^3 * (1.5*(Rsc^2+2*Rsc+1*z)+4.5*(Rsc^2+Rsc)^z^2+4.5*Rsc^2*z^3);
    alpha(i)=exp(lnalpha);
    
    lambda=a*gamma(i)/(b*alpha(i));
    c1_eq(i)=2*c_0/(1+2*lambda*c_0+(1+4*lambda*c_0)^0.5);
end
plot(phi,c1_eq,'r');

r=1;
rc=1.5;
R=r/rc;
r_sc=r;
Rsc=r_sc/rc;

for i=1:l
    z=phi(i)/(1-phi(i));
    lng(i)=log(1-phi(i))+(R^3 + 3*R^2 + 3*R)*z +(3*R^3 + 4.5*R^2)*z^2 + (3*R^2)*z^3;
    gamma(i)=exp(lng(i));
    lnalpha=(2/3)*(r/r_sc)^3 * (1.5*(Rsc^2+2*Rsc+1*z)+4.5*(Rsc^2+Rsc)^z^2+4.5*Rsc^2*z^3);
    alpha(i)=exp(lnalpha);
    
    lambda=a*gamma(i)/(b*alpha(i));
    c1_eq(i)=2*c_0/(1+2*lambda*c_0+(1+4*lambda*c_0)^0.5);
end
plot(phi,c1_eq,'c');

r=1;
rc=2;
R=r/rc;
r_sc=r;
Rsc=r_sc/rc;

for i=1:l
    z=phi(i)/(1-phi(i));
    lng(i)=log(1-phi(i))+(R^3 + 3*R^2 + 3*R)*z +(3*R^3 + 4.5*R^2)*z^2 + (3*R^2)*z^3;
    gamma(i)=exp(lng(i));
    lnalpha=(2/3)*(r/r_sc)^3 * (1.5*(Rsc^2+2*Rsc+1*z)+4.5*(Rsc^2+Rsc)^z^2+4.5*Rsc^2*z^3);
    alpha(i)=exp(lnalpha);
    
    lambda=a*gamma(i)/(b*alpha(i));
end
plot(phi,c1_eq,'m');








