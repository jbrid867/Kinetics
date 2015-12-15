function cfib = cfib_actin(x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


%x(1) is a/b and x(2) is rc

aob=x(1);
rc=x(2);

r=2.5; %protein raidus
rsc=5; %spherocylinder radius

z=zeros(4,1);
cfib=zeros(4,1);
c0=5*10^-6; %initial concentration

phi=[0,37.5*10^-3,75*10^-3,150*10^-3];

alpha=0;
lnalpha=alpha;
lng=alpha;
gamma=alpha;
z=alpha;
gamma=1; 

if phi==0
    
    c1 = 2*c0/(1+2*(aob)*c0 + (1+4*(aob)*c0)^0.5); % no crowder monomer mass
    c2 = (aob)*c1^2;
    cfib=c0-c1-c2;
else
    for i=1:4
        R=r/rc;
        Rsc=rsc/rc;

        z(i)=phi(i)/(1-phi(i));
        lng=log(1-phi(i))+(R^3+3*R^2+R)*z(i)+(3*R^3+4.5*R^2)*z(i)^2+3*(R^3)*(z(i)^3);
        lnalpha=(2/3)*(r/rsc)^3*(1.5*(Rsc^2+2*Rsc+1)*z(i)+4.5*(Rsc^2+Rsc)*z(i)^2+4.5*Rsc^2*z(i)^3);
        gamma=exp(lng);
        alpha=exp(lnalpha);
        lambda=(gamma/alpha)*aob;
        c1=2*c0/(1+2*lambda*c0 + (1+4*lambda*c0)^0.5);
        c2=lambda*c1^2;
        cfib(i)=c0-c1-c2;
    end
end


end

