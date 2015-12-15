%% becker-doring equilibrium model. constant c_1

N=100;
l=101;
phi=linspace(0,0.4,l);
c=zeros(N,l);
lng=zeros(N,l);
gamma=zeros(N,l);

r_c=1; % crowder radius
r=1; % monomer radius
r_sc=1; % spherocylinder radius
R=r/r_c;
Rsc=r_sc/r_c;

R_1=R^3 + 3*R^2 +3*R;
R_2=3*R^3+4.5*R^2;

for j=1:l
    z=phi(j)/(1-phi(j));
    lg=log(1-phi(j));
    lng(1,j)=log(1-phi(j))+(R^3 + 3*R^2 + 3*R)*z +(3*R^3 + 4.5*R^2)*z^2 + (3*R^2)*z^3;
    R=r_sc/r_c;
    for i=2:N
        L=2/3*(i*(r/r_sc)^3 - 1);
        lng(i,j)=lg+(R_1+1.5*L*(R^2+2*R+1))*z + (R_2+4.5*L*R*(R+1))*z^2+(3*R^3+4.5*L*R^2)*z^3;
    end

    for i=1:N
        gamma(i,j)=exp(lng(i,j));
    end
    
    for i=1:N
        c(i,j)=(gamma(1,j)^i)/gamma(i,j);
    end
end
x=1:N;
plot(x,log(c(:,26)),'b',x, log(c(:,51)), '-',x, log(c(:,76)), 'c',x, log(c(:,101)), 'm');

