alpha=zeros(1,101); %10-mer, 100 phi
gamma=zeros(10,101);
lng=zeros(10,101);
phi=linspace(0,0.4,101);
rs=1;
rc=.5;
r_sc=.75;
R=rs/rc;
Rsc=r_sc/rc;

for i=1:101
    z=phi(i)/(1-phi(i));
    alpha(i)=exp((2/3)*(rs/r_sc)^3 * 1.5*(Rsc^2+2*Rsc+1)*z+4.5*(Rsc^2+Rsc)*z^2+4.5*Rsc^2*z^3);
    lng(1,i)=log(1-phi(i))+(R^3 + 3*R^2 + 3*R)*z +(3*R^3 + 4.5*R^2)*z^2 + (3*R^2)*z^3;
    gamma(1,i)=exp(lng(1,i));
    for j=2:10
        gamma(j,i)=gamma(1,i)*alpha(i)^(j-1);
    end
end

plot(phi,gamma(1,:)./alpha)