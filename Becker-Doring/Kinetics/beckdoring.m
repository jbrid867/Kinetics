function dc = beckdoring(t,c,phi,b,factor, k_nuc)
% Becker-Doring Kinetics
%   Takes phi, concentrations and parameters as input. returns the
%   derivatives of concentration

N=length(c); %number of terms
J=zeros(N,1);
J_1=0;
a=b*factor; % factor determined from equilibrium results

dc=zeros(N,1);%derivative vector

nc=3; %critical nucleus

r=2.5; %protein raidus
r_c=1.35; %crowder radius
rsc=5; %spherocylinder radius
R=r/r_c;
Rsc=rsc/r_c;

z=phi/(1-phi);
lng=log(1-phi)+(R^3 + 3*R^2 + 3*R)*z +(3*R^3 + 4.5*R^2)*z^2 + (3*R^2)*z^3;
gamma=exp(lng);
lnalpha=(2/3)*(r/rsc)^3*(1.5*(Rsc^2+2*Rsc+1)*z+4.5*(Rsc^2+Rsc)*z^2+4.5*Rsc^2*z^3);
alpha=exp(lnalpha);
if phi==0
    alpha=1;
    gamma=1;
end


%build J's
J(nc-1)=k_nuc*(gamma*c(1))^nc - b*(alpha^(nc-1))*gamma*c(nc); %flux between monomers and nucleus
for i=nc:N-1
    J(i)=alpha^(i-1) * gamma *(a*gamma*c(1)*c(i)-b*alpha*c(i+1)); %flux between r and r+1 mer
    J_1=J_1+J(i); %sum for c_1
end

dc(1)=(-nc*J(nc-1)-J_1); % special case for c1
for i=nc:N
    
    dc(i)=(J(i-1)-J(i)); % main stuff auto truncated by J(N)=0

end

end

