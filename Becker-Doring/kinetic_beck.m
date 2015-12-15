function dc = kinetic_beck(t,c, phi)
N=1000;
dc=zeros(N,1);
lng=zeros(N,1);
gamma=zeros(N,1);
J=zeros(N-1,1);
J_1=0;
a=1; % this is the association rate, should be a vector but i will set all to one
b=zeros(N,1);

r_c=1; % crowder radius
r=1; % monomer radius
r_sc=1; % spherocylinder radius
R=r/r_c;

R_1=R^3 + 3*R^2 +3*R;
R_2=3*R^3+4.5*R^2;

%build gamma
z=phi/(1-phi);
lng(1)=log(1-phi)+(R^3 + 3*R^2 + 3*R)*z +(3*R^3 + 4.5*R^2)*z^2 + (3*R^2)*z^3;
R=r_sc/r_c;
for i=2:N
    L=2/3*(i*(r/r_sc)^3 - 1);
    lng(i)=log(1-phi)+(R_1+1.5*L*(R^2+2*R+1))*z + (R_2+4.5*L*R*(R+1))*z^2+(3*R^3+4.5*L*R^2)*z^3;
end

for i=1:N
    gamma(i)=exp(lng(i));
    b(i)=(i+1)^2 / i;
end
%gamma built

%build J's
for i=1:N-1
    J(i)=a*gamma(i)*gamma(1)*c(i)*c(1)-b(i+1)*gamma(i+1)*c(i+1);
    J_1=J_1+J(i); %sum for c_1
end

dc(1)=(-J(1)-J_1); % special case for c1

for i=2:N-1
    dc(i)=(J(i-1)-J(i)); % main stuff
end
dc(N)=J(N-1); %truncate


end

